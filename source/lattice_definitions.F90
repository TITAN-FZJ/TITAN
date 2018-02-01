subroutine setup_electric_field(s)
  use mod_f90_kind, only: double
  use mod_constants, only: deg2rad
  use mod_parameters, only: outputunit, dirEfield, dirEfieldvec, EFp, EFt
  use mod_System, only: System
  use mod_mpi_pars, only: myrank, ierr, MPI_Finalize
  implicit none
  type(System), intent(inout) :: s
  real(double), dimension(3) :: zdir

  zdir = [0.d0,0.d0,1.d0]

  if(dirEfield == -1) then
    continue
  else if(dirEfield == -2) then
    dirEfieldvec = dirEfieldvec(1) * s%a1 + dirEfieldvec(2) * s%a2 + dirEfieldvec(3) * s%a3
  else if(dirEfield == -3) then
    dirEfieldvec = [cos(EFp*deg2rad)*sin(EFt*deg2rad), sin(EFp*deg2rad)*sin(EFt*deg2rad), cos(EFt*deg2rad)]
  else if(1 <= dirEfield .and. dirEfield <= s%nAtoms) then
    stop "Not Implemented."
  else
    if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Unknown E-Field direction!')")
    call MPI_Finalize(ierr)
    stop
  end if

  dirEfieldvec = dirEfieldvec / sqrt(dot_product(dirEfieldvec,dirEfieldvec))
  EFp = atan(dirEfieldvec(2) / dirEfieldvec(1))
  EFt = acos(dirEfieldvec(3))
  return
end subroutine

subroutine setup_long_and_trans_current_neighbors(s)
  use mod_f90_kind, only: double
  use mod_system, only: System

  implicit none

  type(System), intent(inout) :: s

  stop "setup_long_and_trans_current_neighbors not Implemented."
    !
    ! ! Calculating vector in-plane perpendicular to electric field direction
    ! if(s%lbulk) then
    !   versor_Eperp = 0.d0
    !   do i = l_nn(1,1), l_nn(1,2) - 1
    !     if(.not. is_parallel(r_nn(:,i), dirEfieldvec)) then
    !       versor_Eperp = cross_unit(r_nn(:,i), dirEfieldvec)
    !       exit
    !     end if
    !   end do
    !   if(0.d0 == dot_product(versor_Eperp, versor_Eperp)) then
    !     !TODO: proper error message
    !     if(0 == myrank) write(outputunit,*) "[lattice_defintions] Error dirEfieldvec"
    !     call MPI_Finalize(ierr)
    !   end if
    ! else
    !   versor_Eperp = cross_unit(zdir, dirEfieldvec)
    ! end if
    !
    !
    ! ! Storing neighbors and angles for longitudinal and transverse currents
    ! sha_longitudinal = 0
    ! sha_transverse   = 0
    ! longitudinal_neighbors = 0
    ! transverse_neighbors = 0
    ! do i=l_nn(1,1), l_nn(1,2)-1
    !   long_cos(i) = dot_product(c_nn(:,i), dirEfieldvec)
    !   if(long_cos(i)>1.d-8) then
    !     longitudinal_neighbors = longitudinal_neighbors+1
    !     sha_longitudinal(longitudinal_neighbors) = i
    !   end if
    !   transv_cos(i) = dot_product(c_nn(:,i),versor_Eperp)
    !   if(transv_cos(i)>1.d-8) then
    !     transverse_neighbors = transverse_neighbors+1
    !     sha_transverse(transverse_neighbors) = i
    !   end if
    ! end do
    ! if(0 == myrank) then
    !   write(outputunit,"('[lattice_definitions] Longitudinal neighbors: ',10(i0,2x))") (sha_longitudinal(i),i=1,longitudinal_neighbors)
    !   write(outputunit,"('[lattice_definitions]   Transverse neighbors: ',10(i0,2x))") (sha_transverse(i),i=1,transverse_neighbors)
    ! end if
    !
end subroutine setup_long_and_trans_current_neighbors


!
! subroutine lattice_definitions(s)
!   use mod_f90_kind, only: double
!   use mod_constants, only: pi
!   use mod_mpi_pars
!   use mod_System, only: System
!   use mod_parameters
!   use mod_tools, only: cross_unit, is_parallel
!   implicit none
!   type(System), intent(inout) :: s
!   integer :: i
!   real(double),dimension(3) :: versor_Eperp
!   real(double), dimension(3) :: zdir
!   zdir = [0.d0,0.d0,1.d0]
!
!   if(dirEfield == -1) then
!     dirEfieldvec = dirEfieldvec / sqrt(dot_product(dirEfieldvec, dirEfieldvec))
!     EFp = atan(dirEfieldvec(2) / dirEfieldvec(1))
!     EFt = acos(dirEfieldvec(3))
!   else if(dirEfield == -2) then
!     dirEfieldvec = dirEfieldvec(1) * pln_a1 + dirEfieldvec(2) * pln_a2
!     dirEfieldvec = dirEfieldvec / sqrt(dot_product(dirEfieldvec, dirEfieldvec))
!     EFp = atan(dirEfieldvec(2) / dirEfieldvec(1))
!     EFt = acos(dirEfieldvec(3))
!   else if(dirEfield == -3) then
!     dirEfieldvec = [cos(EFp*deg2rad)*sin(EFt*deg2rad), sin(EFp*deg2rad)*sin(EFt*deg2rad), cos(EFt*deg2rad)]
!   else if(dirEfield >=1 .and. dirEfield <= pln_cnt(1)) then
!     dirEfieldvec = s%Neighbors(dirEfield)%dirCos
!   else
!     if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Unknown E-Field direction!')")
!     call MPI_Finalize(ierr)
!     stop
!   end if
!
!   ! Calculating out-of-plane unit vector
!   return
! end subroutine lattice_definitions
