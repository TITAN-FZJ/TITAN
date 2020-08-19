! Calculate the derivative of H0 of a slab in the direction ue
subroutine dtdksub(ue,kp,dtdk)
  use mod_kind, only: dp
  use mod_constants,        only: cZero, cI
  use AtomTypes,            only: NeighborIndex
  use mod_system,           only: s => sys
  use mod_parameters,       only: nOrb
  implicit none
  integer :: i, j, k
  real(dp), intent(in)  :: kp(3)
  real(dp), intent(in)  :: ue(3)
  complex(dp),dimension(nOrb,nOrb,s%nAtoms,s%nAtoms),intent(out)  :: dtdk
  complex(dp) :: kpExp

  dtdk = cZero
  ! Mouting derivative of slab's hamiltonian
  ! d/dk t_ij * exp(i*R*k) = i*R t_ij exp(i*R*k)
  do k = 1, s%nNeighbors
    j = s%Neighbors(k)%BasisIndex
    ! kpExp = cI * dot_product(ElectricFieldVector, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    kpExp = cI * dot_product(ue, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    do i = 1, s%nAtoms
      if(.not. s%Neighbors(k)%isHopping(i)) cycle
      dtdk(1:nOrb,1:nOrb,j,i) = dtdk(1:nOrb,1:nOrb,j,i) + s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
    end do
  end do

end subroutine dtdksub