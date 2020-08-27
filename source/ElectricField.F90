! Calculate Electric Field from Input Parameters

module ElectricField
use mod_kind, only: dp
implicit none

integer :: ElectricFieldMode                          !< Direction of in-plane applied electric field
real(dp), dimension(3) :: ElectricFieldVector(3)  !< Direction vector of the electric field
real(dp) :: EFp = 0._dp, EFt = 0._dp                !< Phi and Theta angles of the electric field in spherical coordinates

contains

  subroutine initElectricField(a1,a2,a3)
    use mod_parameters, only: output
    use mod_kind, only: dp
    use mod_constants,  only: deg2rad,rad2deg
    use mod_mpi_pars,   only: myrank,abortProgram
    use mod_tools,      only: rtos
    implicit none
    real(dp), dimension(3), intent(in) :: a1, a2, a3

    if(ElectricFieldMode == -1) then
      continue

    else if(ElectricFieldMode == -2) then
      ElectricFieldVector(1:3) = ElectricFieldVector(1) * a1(1:3) &
                               + ElectricFieldVector(2) * a2(1:3) &
                               + ElectricFieldVector(3) * a3(1:3)

    else if(ElectricFieldMode == -3) then
      ElectricFieldVector(1:3) = [cos(EFp*deg2rad)*sin(EFt*deg2rad), sin(EFp*deg2rad)*sin(EFt*deg2rad), cos(EFt*deg2rad)]

    else if(1 <= ElectricFieldMode) then ! not Implemented anyway!  .and. dirEfield <= s%nAtoms) then
      if(myrank==0) call abortProgram("[initElectricField] Electric Field along Atom Directions not implemented yet!!")

    else
      if(myrank==0) call abortProgram("[initElectricField] Unknown E-Field direction!")

    end if

    ElectricFieldVector(1:3) = ElectricFieldVector(1:3) / sqrt(dot_product(ElectricFieldVector(1:3), ElectricFieldVector(1:3)))
    EFt = acos(ElectricFieldVector(3))*rad2deg
    if((abs(Eft)>1.e-8_dp).and.(abs(abs(Eft)-180._dp)>1.e-8_dp)) then
      EFp = atan2(ElectricFieldVector(2),ElectricFieldVector(1))*rad2deg
    else
      EFp = 0._dp
    end if

    write(output%EField,"('_EFp=',a,'_EFt=',a)") trim(rtos(EFp,"(f7.2)")),trim(rtos(EFt,"(f7.2)"))

  end subroutine initElectricField


  subroutine initLongAndTransCurentNeighbors()
    ! use mod_kind, only: dp
    ! use mod_system,   only: System_type
    use mod_mpi_pars, only: abortProgram
    implicit none

    ! type(System_type), intent(inout) :: s

    call abortProgram("[initLongAndTransCurentNeighbors] Not Implemented!")
      !
      ! ! Calculating vector in-plane perpendicular to electric field direction
      ! if(s%isysdim == 3) then
      !   versor_Eperp = 0._dp
      !   do i = l_nn(1,1), l_nn(1,2) - 1
      !     if(.not. is_parallel(r_nn(:,i), dirEfieldvec)) then
      !       versor_Eperp = cross_unit(r_nn(:,i), dirEfieldvec)
      !       exit
      !     end if
      !   end do
      !   if(0._dp == dot_product(versor_Eperp, versor_Eperp)) then
      !     !TODO: proper error message
      !     if(0 == myrank) write(outputunit,*) "[lattice_defintions] Error dirEfieldvec"
      !     call MPI_Finalize(ierr)
      !   end if
      ! else if(s%isysdim == 2)
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
      !   if(long_cos(i)>1.e-8_dp) then
      !     longitudinal_neighbors = longitudinal_neighbors+1
      !     sha_longitudinal(longitudinal_neighbors) = i
      !   end if
      !   transv_cos(i) = dot_product(c_nn(:,i),versor_Eperp)
      !   if(transv_cos(i)>1.e-8_dp) then
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
