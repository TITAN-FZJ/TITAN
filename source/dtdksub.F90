! Calculate the derivative of H0 of a slab
subroutine dtdksub(kp,dtdk)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cI
  use AtomTypes, only: NeighborIndex
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use ElectricField, only: ElectricFieldVector
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
    kpExp = cI * dot_product(ElectricFieldVector, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    do i = 1, s%nAtoms
      dtdk(1:nOrb,1:nOrb,j,i) = s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
    end do
  end do

end subroutine dtdksub
