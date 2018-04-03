! Calculate hamiltonian of the unit cell
subroutine hamiltk(sys,kp,hk)
  use mod_f90_kind,  only: double
  use mod_constants, only: cI, cZero
  use AtomTypes,     only: NeighborIndex
  use mod_System,    only: ia, System
  use TightBinding,  only: nOrb,nOrb2
  use mod_magnet,    only: lb, sb
  use mod_SOC,       only: ls
  use mod_Umatrix,   only: hee
  implicit none
  integer :: i, j, k
  real(double), intent(in) :: kp(3)
  type(System), intent(in) :: sys
  complex(double), dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2), intent(out) :: hk
  complex(double) :: tmp(nOrb,nOrb)
  complex(double) :: kpExp
  hk = cZero

  ! Mouting slab hamiltonian
  !dir$ ivdep:loop
  do i=1, sys%nAtoms
    hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
    hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)

    hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                         + lb(1:nOrb2,1:nOrb2,i) + sb(1:nOrb2,1:nOrb2,i) + hee(1:nOrb2,1:nOrb2,i) &
                                         + ls(1:nOrb2,1:nOrb2,i)
  end do

  !dir$ ivdep:loop
  do k = 1, sys%nNeighbors
    j = sys%Neighbors(k)%BasisIndex
    kpExp = exp(cI * dot_product(kp, sys%Neighbors(k)%CellVector))

    !dir$ ivdep:loop
    do i = 1, sys%nAtoms
      if(.not. sys%Neighbors(k)%isHopping(i)) cycle
      tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
      hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:nOrb,1:nOrb)
      hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:nOrb,1:nOrb)
    end do
  end do

  ! do i = ia(1,1), ia(4,sys%nAtoms)
  !   do j = i, ia(4,sys%nAtoms)
  !     if(abs(hk(j,i)-conjg(hk(i,j))) > 1.d-12) then
  !       print *, i,j,abs(hk(j,i)-conjg(hk(i,j)))
  !     end if
  !   end do
  ! end do

end subroutine hamiltk

! Calculate hamiltonian of the unit cell
! and the spin-orbit coupling contribution separately
subroutine hamiltklinearsoc(sys,kp,hk,vsoc)
  use mod_f90_kind,  only: double
  use mod_constants, only: cZero, cI
  use mod_system,    only: ia, System
  use AtomTypes,     only: NeighborIndex
  use TightBinding,  only: nOrb,nOrb2
  use mod_SOC,       only: ls
  use mod_magnet,    only: lb, sb
  use mod_Umatrix,   only: hee
  implicit none
  integer :: i, j, k
  real(double), intent(in) :: kp(3)
  type(System), intent(in) :: sys
  complex(double),dimension(sys%nAtoms*nOrb2,sys%nAtoms*nOrb2),intent(out)  :: hk,vsoc
  complex(double) :: tmp(nOrb, nOrb)
  complex(double) :: kpExp

  hk = cZero
  vsoc = cZero

  ! Mouting slab hamiltonian
  do i=1,sys%nAtoms
    hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = sys%Types(sys%Basis(i)%Material)%onSite
    hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = sys%Types(sys%Basis(i)%Material)%onSite
    hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                         + lb(:,:,i) + sb(:,:,i) + hee(:,:,i)
    vsoc(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = ls(:,:,i)
  end do

  do k = 1, sys%nNeighbors
    j = sys%Neighbors(k)%BasisIndex
    kpExp = exp(cI * dot_product(kp,sys%Neighbors(k)%CellVector))

    do i = 1, sys%nAtoms
      tmp = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
      hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp
      hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp
    end do
  end do
end subroutine hamiltklinearsoc
