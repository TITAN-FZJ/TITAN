! Calculate Electric Field from Input Parameters

module ElectricField
use mod_f90_kind, only: double
implicit none

integer :: ElectricFieldMode                          !< Direction of in-plane applied electric field
real(double), dimension(3) :: ElectricFieldVector(3)  !< Direction vector of the electric field
real(double) :: EFp = 0.d0, EFt = 0.d0                !< Phi and Theta angles of the electric field in spherical coordinates
character(len=50) :: strElectricField = ""            ! EFt, EFp

contains

  subroutine initElectricField(a1,a2,a3)
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_mpi_pars, only: abortProgram, myrank
    implicit none
    real(double), dimension(3), intent(in) :: a1, a2, a3

    if(ElectricFieldMode == -1) then
      continue

    else if(ElectricFieldMode == -2) then
      ElectricFieldVector(1:3) = ElectricFieldVector(1) * a1(1:3) &
                               + ElectricFieldVector(2) * a2(1:3) &
                               + ElectricFieldVector(3) * a3(1:3)

    else if(ElectricFieldMode == -3) then
      ElectricFieldVector(1:3) = [cos(EFp*pi/180)*sin(EFt*pi/180), sin(EFp*pi/180)*sin(EFt*pi*180), cos(EFt*pi/180)]

    else if(1 <= ElectricFieldMode) then ! not Implemented anyway!  .and. dirEfield <= s%nAtoms) then
      if(myrank.eq.0) call abortProgram("[initElectricField] Electric Field along Atom Directions not implemented yet!!")

    else
      if(myrank.eq.0) call abortProgram("[initElectricField] Unknown E-Field direction!")

    end if

    ElectricFieldVector(1:3) = ElectricFieldVector(1:3) / sqrt(dot_product(ElectricFieldVector(1:3), ElectricFieldVector(1:3)))
    EFp = atan(ElectricFieldVector(2) / ElectricFieldVector(1))
    EFt = acos(ElectricFieldVector(3))

    write(strElectricField,"('_EFp=',es8.1,'_EFt=',es8.1)") EFp,EFt

    return
  end subroutine initElectricField


  subroutine initLongAndTransCurentNeighbors(s)
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
  end subroutine initLongAndTransCurentNeighbors


end module ElectricField