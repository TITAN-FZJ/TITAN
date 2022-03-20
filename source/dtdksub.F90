! Calculate the derivative of H0 of a slab in the direction ue
subroutine dtdksub(ue,kp,dtdk)
  use mod_kind,       only: dp
  use mod_constants,  only: cZero, cI
  use AtomTypes,      only: NeighborIndex
  use mod_system,     only: s => sys
  implicit none
  integer :: i, j, k
  real(dp), intent(in)  :: kp(3)
  real(dp), intent(in)  :: ue(3)
  complex(dp),dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms),intent(out)  :: dtdk
  complex(dp) :: kpExp

  dtdk = cZero
  ! Mouting derivative of slabs hamiltonian
  ! d/dk t_ij * exp(i*R*k) = i*R t_ij exp(i*R*k)
  do k = 1, s%nNeighbors
    j = s%Neighbors(k)%BasisIndex
    ! kpExp = cI * dot_product(ElectricFieldVector, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    kpExp = cI * dot_product(ue, s%Neighbors(k)%CellVector) * exp(cI * dot_product(kp,s%Neighbors(k)%CellVector))
    do i = 1, s%nAtoms
      if(s%Neighbors(k)%isHopping(i)) then
        dtdk(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb,j,i) = dtdk(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb,j,i) + s%Neighbors(k)%t0i(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb,i) * kpExp
      end if
    end do
  end do

end subroutine dtdksub