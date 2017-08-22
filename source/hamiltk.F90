! Calculate hamiltonian of a slab containing
! Npl layers + 1 layer of Empty spheres on each side
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!  1     2     3     4      Npl-1   Npl  Npl+1  Npl+2
!         <-S-> <S-1>           <S-1> <-S->
subroutine hamiltk(kp,hk)
  use mod_f90_kind, only: double
  use mod_constants, only: zi, zero
  use AtomTypes, only: NeighborIndex
  use mod_System, only: ia, s => sys
  use TightBinding, only: nOrb
  use mod_magnet, only: lb, sb
  use mod_SOC, only: ls, socscale
  use mod_mpi_pars, only: abortProgram
  use mod_Umatrix, only: hee
  implicit none
  integer :: i, j, k,l,m
  real(double), intent(in) :: kp(3)
  complex(double), dimension(s%nAtoms*2*nOrb, s%nAtoms*2*nOrb), intent(out) :: hk
  complex(double) :: tmp(nOrb,nOrb)
  complex(double) :: kpExp

  hk = zero

  ! Mouting slab hamiltonian
  do i=1, s%nAtoms
    hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = s%Types(s%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
    hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = s%Types(s%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)

    hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                         + lb(:,:,i) + sb(:,:,i) + hee(:,:,i) &
                                         + socscale * s%Types(s%Basis(i)%Material)%Lambda * ls
  end do

  do k = 1, s%nNeighbors
    j = s%Neighbors(k)%BasisIndex
    kpExp = exp(zi * dot_product(kp, s%Neighbors(k)%CellVector))

    do i = 1, s%nAtoms
      !if(.not. s%Neighbors(k)%isHopping(i)) cycle
      tmp(1:nOrb,1:nOrb) = s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i)
      tmp = tmp * kpExp
      hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:nOrb,1:nOrb)
      hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:nOrb,1:nOrb)
    end do
  end do

  do i = ia(1,1), ia(4,s%nAtoms)
    do j = i, ia(4,s%nAtoms)
      if(abs(hk(j,i)-conjg(hk(i,j))) > 1.d-12) then
        print *, i,j,abs(hk(j,i)-conjg(hk(i,j)))
      end if
    end do
  end do

  return
end subroutine hamiltk

! Calculate hamiltonian of a slab containing
! Npl layers + 1 layer of Empty spheres on each side
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!  1     2     3     4      Npl-1   Npl  Npl+1  Npl+2
!         <-S-> <S-1>           <S-1> <-S->
subroutine hamiltklinearsoc(kp,hk,vsoc)
  use mod_f90_kind,      only: double
  use mod_constants,     only: zero, zi
  use mod_system,        only: ia, s => sys
  use AtomTypes, only: NeighborIndex
  use TightBinding, only: nOrb
  use mod_SOC,    only: socscale, ls
  use mod_magnet,        only: lb, sb
  use mod_Umatrix, only: hee
  implicit none
  integer :: i, j, k
  real(double), intent(in)  :: kp(3)
  complex(double),dimension(s%nAtoms*2*nOrb,s%nAtoms*2*nOrb),intent(out)  :: hk,vsoc
  complex(double) :: tmp(nOrb, nOrb)
  complex(double) :: kpExp

  hk = zero
  vsoc = zero

  ! Mouting slab hamiltonian
  do i=1,s%nAtoms
    hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = s%Types(s%Basis(i)%Material)%onSite
    hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = s%Types(s%Basis(i)%Material)%onSite
    hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                         + lb(:,:,i) + sb(:,:,i) + hee(:,:,i)
    vsoc(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = socscale * s%Types(s%Basis(i)%Material)%Lambda * ls
  end do

  do k = 1, s%nNeighbors
    j = s%Neighbors(k)%BasisIndex
    kpExp = exp(zi * dot_product(kp,s%Neighbors(k)%CellVector))

    do i = 1, s%nAtoms
      tmp = s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
      hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp
      hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp
    end do
  end do
  return
end subroutine hamiltklinearsoc
