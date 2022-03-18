module mod_checkpoint
!! This module holds the functions and subroutines to save and recover the time propagation states
  implicit none

contains

subroutine save_state(rank,s,dimH,nkpt,rtime,step,eval_kn,evec_kn,mx_t,my_t,mz_t)
  !! This subroutine is to Save current time-propagation state
    use mod_kind,               only: dp, int32, int64
    use mod_parameters,         only: output
    use mod_tools,              only: itos
    use mod_io,                 only: log_warning
    use mod_time_propagator_io, only: write_header_time_prop
    use mod_system,             only: System_type
    implicit none
    integer(int32), intent(in) :: rank
    !! Rank ID of process
    type(System_type), intent(in)   :: s
    !! System derived type containing number of orbitals and number of atoms
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
    real(dp),    dimension(s%nOrb,s%nAtoms),    intent(in) :: mx_t,my_t,mz_t
    !! Orbital- and site-dependent magnetization components

    ! Local variables:
    integer(int32) :: file_unit = 6101
    !! File unit
    character(len=500) :: output_file
    !! Filename
    character(len=100) :: formatvar
    !! Format variable
    integer(int32) :: i,j,mu
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

    ! Writing current orbital-dependent magnetization vector to use for total torque calculation (dM/dt) when recovering
    write(formatvar,fmt="(a,i0,a)") '(',3*s%nOrb*nAtoms,'(es16.8e3,2x))'
    write(unit=file_unit,fmt=formatvar) ((mx_t(mu,i),my_t(mu,i),mz_t(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)

    close(file_unit) 

    call log_warning("save_state", "Checkpoint saved successfully.")

  end subroutine save_state


  function recover_state(rank,s,dimH,nkpt,rtime,step,eval_kn,evec_kn,mx_t,my_t,mz_t) result(success)
    !! This subroutine is to recover current time-propagation state
    !! Note that the recover has to use the same MPI setup as used when saving the state, since different files per rank are written
    use mod_kind,               only: dp,int32,int64
    use mod_io,                 only: log_warning
    use mod_parameters,         only: output
    use mod_time_propagator_io, only: check_header_time_prop
    use mod_system,             only: System_type
    implicit none
    integer(int32), intent(in) :: rank
    !! Rank ID of process
    type(System_type), intent(in)   :: s
    !! System derived type containing number of orbitals and number of atoms
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
    real(dp),    dimension(s%nOrb,s%nAtoms),    intent(out) :: mx_t,my_t,mz_t
    !! Orbital- and site-dependent magnetization components
    logical :: success
    !! Indication of when a state was recovered or not

    ! Local variables:
    integer(int32) :: file_unit = 6102
    !! File unit
    character(len=500) :: output_file
    !! Filename
    character(len=100) :: formatvar
    !! Format variable
    logical :: success_header
    !! Indication of when a state was recovered or not
    integer(int32) :: i,j,mu,stat_check,stat_temp
    integer(int64) :: k

    success = .false.

    ! Try to recover current state (Read file)
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/checkpoint',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),rank,trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=file_unit, file=output_file, status='old', action='read', iostat=stat_check) 
    if (stat_check /= 0) then
      call log_warning("recover_state", "No checkpoint file found.")
      return
    end if

    ! Read header and check if it is the same
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

    ! Reading eigenvalues and eigenvectors
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

    ! Reading current orbital-dependent magnetization vector to use for total torque calculation (dM/dt) when recovering
    write(formatvar,fmt="(a,i0,a)") '(',3*nOrb*nAtoms,'(es16.8e3,2x))'
    read(unit=file_unit,fmt=formatvar, iostat=stat_temp) ((mx_t(mu,i),my_t(mu,i),mz_t(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)
    stat_check = stat_check + stat_temp

    if (stat_check /= 0) then
      call log_warning("recover_state", "Checkpoint file with correct header found, but there was a problem reading eigenvalues and eigenvectors. Cannot continue from previous point.")
      close(file_unit)
      return
    end if

    ! Read header and check if it is the same
    success = .true.
    call log_warning("recover_state", "Checkpoint recovered.")

    close(file_unit) 


  end function recover_state


end module mod_checkpoint