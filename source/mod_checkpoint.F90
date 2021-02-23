! This module holds the functions and subroutines to save and recover the time propagation states
module mod_checkpoint
  implicit none

contains

  !! This subroutine is to Save current state
  subroutine save_state(rank,dimH,nkpt,rtime,step,eval_kn,evec_kn)
    use mod_kind,               only: dp, int32, int64
    use mod_parameters,         only: output
    use mod_tools,              only: itos
    use mod_io,                 only: log_warning
    use mod_time_propagator_io, only: write_header_time_prop
    implicit none

    integer(int32), intent(in) :: rank
    !! Rank ID of process
    integer(int32), intent(in) :: dimH
    !! Dimension of the hamiltonian and eigenvalues/eigenvectors
    integer(int64), intent(in) :: nkpt
    !! Number of k-points
    real(dp),                               intent(in) :: rtime
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
    write(unit=file_unit,fmt="(2(es16.9,2x))") rtime, step

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
  function recover_state(rank,dimH,nkpt,rtime,step,eval_kn,evec_kn) result(success)
    use mod_kind,               only: dp,int32,int64
    use mod_io,                 only: log_warning
    use mod_parameters,         only: output
    use mod_time_propagator_io, only: check_header_time_prop
    implicit none
    integer(int32), intent(in) :: rank
    !! Rank ID of process
    integer(int32), intent(in) :: dimH
    !! Dimension of the hamiltonian and eigenvalues/eigenvectors
    integer(int64), intent(in) :: nkpt
    !! Number of k-points
    real(dp),                               intent(out) :: rtime
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
    read(unit=file_unit,fmt=*, iostat=stat_temp) rtime, step
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