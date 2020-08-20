! This module holds the functions and subroutines to save and recover the time propagation states
module mod_checkpoint
  implicit none
!----------------------------------------   General recipe and Library --------------------------------------------------!
! Library needed:
! To install library, use: wget https://code.nxg.name/nxg/libsigwatch/libsigwatch-1.0.tar.gz 
! Extract the libsigwatch-1.0.tar.gz  : % tar xvzf libsigwatch-1.0.tar.gz
! Navigate to the folder libsigwatch-1.0: % cd libsigwatch-1.0
! To configure, build and install, use:
! % ./configure. (Or, change configuration path to Home:  ./configure --prefix=$HOME)
! % make
! % make install
!
!************** Then note: **********************************************************************************************!
! 
! Libraries have been installed in:
!    /home/el144313/lib
!
! If you ever happen to want to link against installed libraries
! in a given directory, LIBDIR, you must either use libtool, and
! specify the full pathname of the library, or use the '-LLIBDIR'
! flag during linking and do at least one of the following:
!    - add LIBDIR to the 'LD_LIBRARY_PATH' environment variable
!      during execution
!    - add LIBDIR to the 'LD_RUN_PATH' environment variable
!      during linking
!    - use the '-Wl,-rpath -Wl,LIBDIR' linker flag
!    - have your system administrator add LIBDIR to '/etc/ld.so.conf'

! See any operating system documentation about shared libraries for
! more information, such as the ld(1) and ld.so(8) manual pages.
! 
!----------------------------------------------------------------------------------------------------------------------!
! Implementation in the code: 
! 1- Look for a state file in which all information required to restart the state are included.
!    * If it exists, read it and restore the state.
!    * Else, create an initial state.
! 2- Save the state after recieving a signal.
! 
! Slurm features for signal handeling: see this link: https://nxg.me.uk/dist/sigwatch
! 1- Scancel --signal USR1$JOB_ID
! 2- Sbatch --signal= INT@60  ! 60 Here is an example: the signal is sent 60 seconds before the job is killed.
! " Modify the program to look periodically for this signal, if signal is recieved, checkpoint and exit."
! 
! To emplement signals:
! 1- Register a signal handler (a function that will modify a global variable when recieving a signal.)
! 2- Test the value of the global variable periodically (At a moment when the state is consistent and easy to recreate.)
! 3. If the value indicates so, save state to disk (and optionally stop)
!
!----------------------------------------------------------------------------------------------------------------------!
! These lines should be added to the submission file:
! ---------------------------------------------------------------------------------------------------------------------!
! #SBATCH —output= results
! #SBATCH —open-mode=append ! Note that Append is necessary here
! #SBATCH —signal= INT@90   ! Example: Send SIGINT 90 seconds before job is killed.
!----------------------------------------------------------------------------------------------------------------------!
!
! ** To have the job re-queued automatically, put these lines at the end of submission file:
!----------------------------------------------------------------------------------------------------------------------!
! date
! echo “restarted ${SLURM_RESTART_COUNT-0} time(s)”
! ./crsigvalibrccount || control require $SLURM_JOB_ID   
!----------------------------------------------------------------------------------------------------------------------!

contains

  !! This subroutine is to Save current state
  subroutine save_state(rank,dimH,nkpt,time,step,eval_kn,evec_kn)
    use mod_kind, only: dp, int32, int64
    use mod_parameters, only: output
    use mod_tools, only: itos
    use mod_io, only: log_warning
    use mod_time_propagator_io, only: write_header_time_prop
    implicit none

    integer(int32), intent(in) :: rank
    !! Rank ID of process
    integer(int32), intent(in) :: dimH
    !! Dimension of the hamiltonian and eigenvalues/eigenvectors
    integer(int64), intent(in) :: nkpt
    !! Number of k-points
    real(dp),                               intent(in) :: time
    !! Last calculated time instant
    real(dp),                               intent(in) :: step 
    !! Last step-size
    real(dp),    dimension(dimH,nkpt),      intent(in) :: eval_kn
    !! Eigenvalues for all k-points
    complex(dp), dimension(dimH,dimH,nkpt), intent(in) :: evec_kn
    !! Eigenvectors (in columns) for all k-points

    !! Local variables:
    integer(int32) :: file_unit = 6101
    !! File unit
    character(len=500) :: output_file
    !! Filename
    character(len=100) :: formatvar
    !! Format variable
    integer(int32) :: i,j
    integer(int64) :: k

    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/checkpoint',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),rank,trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=file_unit, file=output_file, status='replace', form='formatted') 

    call write_header_time_prop(file_unit, "# dim = " // trim(itos(dimH)) // " , nkpt = " // trim(itos(nkpt)) )

    ! Writing iteration, time and step size
    write(unit=file_unit,fmt="(2(es16.9,2x))") time, step

    ! Writing eigenvalues and eigenvectors
    write(formatvar,fmt="(a,i0,a)") '(',dimH,'(es16.8e3,2x))'
    do k=1,nkpt
      write(unit=file_unit,fmt=formatvar) (eval_kn(j,k),j=1,dimH)
      do i=1,dimH
        write(unit=file_unit,fmt=formatvar) (evec_kn(i,j,k),j=1,dimH)
      end do
    end do

    close(file_unit) 

    call log_warning("save_state", "Checkpoint saved successfully.")

  end subroutine save_state


  !! This subroutine is to recover current state
  function recover_state(rank,dimH,nkpt,time,step,eval_kn,evec_kn) result(success)
    use mod_kind, only: dp, int32, int64
    use mod_io, only: log_warning
    use mod_parameters, only: output
    use mod_time_propagator_io, only: check_header_time_prop
    implicit none
    integer(int32), intent(in) :: rank
    !! Rank ID of process
    integer(int32), intent(in) :: dimH
    !! Dimension of the hamiltonian and eigenvalues/eigenvectors
    integer(int64), intent(in) :: nkpt
    !! Number of k-points
    real(dp),                               intent(out) :: time
    !! Last calculated time instant
    real(dp),                               intent(out) :: step 
    !! Last step-size
    real(dp),    dimension(dimH,nkpt),      intent(out) :: eval_kn
    !! Eigenvalues for all k-points
    complex(dp), dimension(dimH,dimH,nkpt), intent(out) :: evec_kn
    !! Eigenvectors (in columns) for all k-points
    logical :: success
    !! Indication of when a state was recovered or not

    !! Local variables:
    integer(int32) :: file_unit = 6102
    !! File unit
    character(len=500) :: output_file
    !! Filename
    character(len=100) :: formatvar
    !! Format variable
    logical :: success_header
    !! Indication of when a state was recovered or not
    integer(int32) :: i,j,stat_check, stat_temp
    integer(int64) :: k

    success = .false.

    ! Try to recover current state (Read file)
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/checkpoint',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),rank,trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=file_unit, file=output_file, status='old', action='read', iostat=stat_check) 
    if (stat_check /= 0) then
      call log_warning("recover_state", "No checkpoint file found.")
      return
    end if

    ! Read header and check if it's the same
    call check_header_time_prop(file_unit,success_header)
    if (.not.success_header) then
      call log_warning("recover_state", "Checkpoint file found, but header differs. Cannot continue from previous point.")
      close(file_unit)
      return
    end if

    ! Reading iteration, time and step size
    read(unit=file_unit,fmt=*, iostat=stat_temp) time, step
    if (stat_check /= 0) then
      call log_warning("recover_state", "Checkpoint file with correct header found, but there was a problem reading previous time and step. Cannot continue from previous point.")
      close(file_unit)
      return
    end if

    ! Writing eigenvalues and eigenvectors
    write(formatvar,fmt="(a,i0,a)") '(',dimH,'(es16.8e3,2x))'
    stat_check = 0
    do k=1,nkpt
      read(unit=file_unit,fmt=formatvar, iostat=stat_temp) (eval_kn(j,k),j=1,dimH)
      stat_check = stat_check + stat_temp
      do i=1,dimH
        read(unit=file_unit,fmt=formatvar, iostat=stat_temp) (evec_kn(i,j,k),j=1,dimH)
        stat_check = stat_check + stat_temp
      end do
    end do

    if (stat_check /= 0) then
      call log_warning("recover_state", "Checkpoint file with correct header found, but there was a problem reading eigenvalues and eigenvectors. Cannot continue from previous point.")
      close(file_unit)
      return
    end if

    ! Read header and check if it's the same
    success = .true.
    call log_warning("recover_state", "Checkpoint recovered.")

    close(file_unit) 


  end function recover_state


end module mod_checkpoint