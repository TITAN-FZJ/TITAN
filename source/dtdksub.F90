! Calculate the derivative of H0 of a slab
subroutine dtdksub(kp,dtdk)
  use mod_f90_kind,   only: double
  use mod_constants,  only: cZero, cI
  use AtomTypes,      only: NeighborIndex
  use mod_system,     only: s => sys
  use mod_parameters, only: nOrb
  use ElectricField,  only: ElectricFieldVector
  use mod_imRK4_parameters, only: polarization_vector_e
  implicit none
  integer :: i, j, k
  real(double), intent(in)  :: kp(3)
  complex(double),dimension(nOrb,nOrb,s%nAtoms,s%nAtoms),intent(out)  :: dtdk
  complex(double) :: kpExp

  dtdk = cZero
  ! Mouting derivative of slab's hamiltonian
  ! d/dk t_ij * exp(i*R*k) = i*R t_ij exp(i*R*k)
  do k = 1, s%nNeighbors
    j = s%Neighbors(k)%BasisIndex
    ! kpExp = cI * dot_product(ElectricFieldVector, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    kpExp = cI * dot_product(polarization_vector_e, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    do i = 1, s%nAtoms
      if(.not. s%Neighbors(k)%isHopping(i)) cycle
      dtdk(1:nOrb,1:nOrb,j,i) = dtdk(1:nOrb,1:nOrb,j,i) + s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
    end do
  end do

end subroutine dtdksub

! Calculate the derivative of H0 of a slab using the electric vector potential inside
subroutine dtdksub_e(kp,dtdk)
  use mod_f90_kind,   only: double
  use mod_constants,  only: cZero, cI
  use AtomTypes,      only: NeighborIndex
  use mod_system,     only: s => sys
  use mod_parameters, only: nOrb
  use ElectricField,  only: ElectricFieldVector
  use mod_imRK4,      only: evec_potent
  implicit none
  integer :: i, j, k
  real(double), intent(in)  :: kp(3)
  complex(double),dimension(nOrb,nOrb,s%nAtoms,s%nAtoms),intent(out)  :: dtdk
  complex(double) :: kpExp
  real(double)    :: A_t(3), t

  dtdk = cZero
  ! Mouting derivative of slab's hamiltonian
  ! d/dk t_ij * exp(i*R*k) = i*R t_ij exp(i*R*k)
  do k = 1, s%nNeighbors
    j = s%Neighbors(k)%BasisIndex
    ! kpExp = cI * dot_product(ElectricFieldVector, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    call evec_potent(t, A_t)
    kpExp = cI * dot_product(A_t, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    do i = 1, s%nAtoms
      if(.not. s%Neighbors(k)%isHopping(i)) cycle
      dtdk(1:nOrb,1:nOrb,j,i) = dtdk(1:nOrb,1:nOrb,j,i) + s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
    end do
  end do

end subroutine dtdksub_e