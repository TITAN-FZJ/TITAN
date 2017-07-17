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
  use mod_parameters, only: offset
  use ElectricField, only: ElectricFieldVector
  implicit none
  integer :: i, j, l
  real(double), intent(in)  :: kp(3)
  complex(double),dimension(s%nAtoms,s%nAtoms,nOrb,nOrb),intent(out)  :: dtdk
  type(NeighborIndex), pointer :: current
  dtdk = zero

  ! Mouting derivative of slab's hamiltonian

  do i=1, s%nAtoms
    do j = 1, s%nAtoms
      do l = 1, s%nStages
        current => s%Basis(i)%NeighborList(l,j)%head
        do while(associated(current))
          dtdk(j,i,:,:) = dtdk(j,i,:,:) + zi * dot_product(ElectricFieldVector, s%Neighbors(current%index)%Position) * s%Neighbors(current%index)%t0i(:,:,i) * exp(zi * dot_product(kp, s%Neighbors(current%index)%CellVector))
          current => current%next
        end do
      end do
    end do
  end do
  return
end subroutine dtdksub
