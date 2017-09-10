! Calculate the derivative of H0 of a slab containing
! Npl layers (ES not included)
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!        1     2     3       Npl-2 Npl-1  Npl
!         <-S-> <S-1>           <S-1> <-S->
subroutine dtdksub(kp,dtdk)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, zi
  use AtomTypes, only: NeighborIndex
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use ElectricField, only: ElectricFieldVector
  implicit none
  integer :: i, j, k
  real(double), intent(in)  :: kp(3)
  complex(double),dimension(s%nAtoms,s%nAtoms,nOrb,nOrb),intent(out)  :: dtdk
  complex(double) :: kpExp
  dtdk = zero

  ! Mouting derivative of slab's hamiltonian

  do k = 1, s%nNeighbors
    j = s%Neighbors(k)%BasisIndex
    kpExp = zi * dot_product(ElectricFieldVector, s%Neighbors(k)%Position) * exp(zi * dot_product(kp,s%Neighbors(k)%CellVector))
    do i = 1, s%nAtoms
      dtdk(j,i,1:nOrb,1:nOrb) = s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
    end do
  end do

  return
end subroutine dtdksub
