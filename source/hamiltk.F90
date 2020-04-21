! Calculate hamiltonian of the unit cell
subroutine hamiltk(sys,kp,hk)
  use mod_f90_kind,   only: double
  use mod_constants,  only: cI, cZero
  use AtomTypes,      only: NeighborIndex
  use mod_System,     only: ia, System
  use mod_parameters, only: nOrb,nOrb2
  use mod_magnet,     only: lb, sb
  use mod_SOC,        only: ls
  use mod_Umatrix,    only: hee
  implicit none
  integer :: i, j, k
  real(double), intent(in) :: kp(3)
  type(System), intent(in) :: sys
  complex(double), dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2), intent(out) :: hk
  complex(double) :: tmp(nOrb,nOrb)
  complex(double) :: kpExp

  hk = cZero

  ! Mouting slab hamiltonian

  ! On-site terms
  !dir$ ivdep:loop
  do i=1, sys%nAtoms
    ! spin-up on-site tight-binding term
    hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
    ! spin-down on-site tight-binding term
    hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
    ! External magnetic field (orbital + spin) + Electron-electron interaction (Hubbard) + Spin-orbit coupling
    hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                         + lb(1:nOrb2,1:nOrb2,i) + sb(1:nOrb2,1:nOrb2,i) + hee(1:nOrb2,1:nOrb2,i) &
                                         + ls(1:nOrb2,1:nOrb2,i)
  end do

  ! Inter-site hopping terms
  !dir$ ivdep:loop
  do k = 1, sys%nNeighbors
    j = sys%Neighbors(k)%BasisIndex
    ! exp(ik.(R_i-R_j))
    kpExp = exp(cI * dot_product(kp, sys%Neighbors(k)%CellVector))

    !dir$ ivdep:loop
    do i = 1, sys%nAtoms
      if(.not. sys%Neighbors(k)%isHopping(i)) cycle
      tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
      ! Spin-up
      hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:nOrb,1:nOrb)
      ! Spin-down
      hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:nOrb,1:nOrb)
    end do
  end do

  ! ! Test if hamiltonian is Hermitian (to be commented out, uncomment to use it)
  ! do i = ia(1,1), ia(4,sys%nAtoms)
  !   do j = i, ia(4,sys%nAtoms)
  !     if(abs(hk(j,i)-conjg(hk(i,j))) > 1.d-15) then
  !       write(*,"('Hamiltonian not hermitian',i0,2x,i0,2x,es11.4)") i,j,abs(hk(j,i)-conjg(hk(i,j)))
  !     end if
  !   end do
  ! end do

end subroutine hamiltk

! Calculate hamiltonian of the unit cell
! and the spin-orbit coupling contribution separately
subroutine hamiltklinearsoc(sys,kp,hk,vsoc)
  use mod_f90_kind,   only: double
  use mod_constants,  only: cZero, cI
  use mod_system,     only: ia, System
  use AtomTypes,      only: NeighborIndex
  use mod_parameters, only: nOrb,nOrb2
  use mod_magnet,     only: lb, sb
  use mod_SOC,        only: ls
  use mod_Umatrix,    only: hee
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

  ! On-site terms
  !dir$ ivdep:loop
  do i=1,sys%nAtoms
    ! spin-up on-site tight-binding term
    hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
    ! spin-down on-site tight-binding term
    hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
    ! External magnetic field (orbital + spin) + Electron-electron interaction (Hubbard)
    hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                         + lb(1:nOrb2,1:nOrb2,i) + sb(1:nOrb2,1:nOrb2,i) + hee(1:nOrb2,1:nOrb2,i)
    ! Spin-orbit coupling term
    vsoc(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = ls(1:nOrb2,1:nOrb2,i)
  end do

  ! Inter-site hopping terms
  !dir$ ivdep:loop
  do k = 1, sys%nNeighbors
    j = sys%Neighbors(k)%BasisIndex
    ! exp(ik.(R_i-R_j))
    kpExp = exp(cI * dot_product(kp,sys%Neighbors(k)%CellVector))

    !dir$ ivdep:loop
    do i = 1, sys%nAtoms
      if(.not. sys%Neighbors(k)%isHopping(i)) cycle
      tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
      ! Spin-up
      hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:nOrb,1:nOrb)
      ! Spin-down
      hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:nOrb,1:nOrb)
    end do
  end do
end subroutine hamiltklinearsoc
