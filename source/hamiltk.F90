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
  use mod_System, only: s => sys
  use TightBinding, only: nOrb
  use mod_magnet, only: lb, sb
  use mod_SOC, only: ls, socscale
  implicit none
  integer :: i, j, k, l
  integer :: i_u0, i_u1, i_d0, i_d1, j_u0, j_u1, j_d0, j_d1
  integer :: index
  real(double), intent(in) :: kp(3)
  complex(double) :: hee(2*nOrb, 2*nOrb, s%nAtoms)
  complex(double), dimension(s%nAtoms*2*nOrb, s%nAtoms*2*nOrb), intent(out) :: hk
  complex(double), dimension(nOrb, nOrb) :: tmp
  type(NeighborIndex), pointer :: current
  hk = zero

  call U_matrix(hee, s%nAtoms, nOrb)

  ! Mouting slab hamiltonian
  do i=1, s%nAtoms
    i_u0 = (i-1)*2*nOrb + 1
    i_u1 = i_u0 + nOrb - 1
    i_d0 = i_u0 + nOrb
    i_d1 = i_d0 + nOrb - 1

    hk(i_u0:i_u1, i_u0:i_u1) = s%Types(s%Basis(i)%Material)%onSite ! + sb hee +socscale*lambda*ls
    hk(i_d0:i_d1, i_d0:i_d1) = s%Types(s%Basis(i)%Material)%onSite ! + sb hee +socscale*lambda*ls

    hk(i_u0:i_d1, i_u0:i_d1) = hk(i_u0:i_d1, i_u0:i_d1) &
                             + lb(:,:,i) + sb(:,:,i) + hee(:,:,i) &
                             + socscale * s%Types(s%Basis(i)%Material)%Lambda * ls

    do j = 1, s%nAtoms
      j_u0 = (j-1) * 2 * nOrb + 1
      j_u1 = j_u0 + nOrb - 1
      j_d0 = j_u0 + nOrb
      j_d1 = j_d0 + nOrb - 1
      do k = 1, s%nStages
        current => s%Basis(i)%NeighborList(k,j)%head
        do while(associated(current))
          index = current%index
          !tmp = s%Neighbors(index)%t0i(1:nOrb,1:nOrb,i)
          !tmp = tmp * exp(zi*dot_product(kp,s%Neighbors(index)%CellVector))
           hk(j_u0:j_u1, i_u0:i_u1) = hk(j_u0:j_u1, i_u0:i_u1) + s%Neighbors(index)%t0i(1:nOrb,1:nOrb,i)*exp(zi*dot_product(kp,s%Neighbors(index)%CellVector))
           hk(j_d0:j_d1, i_d0:i_d1) = hk(j_d0:j_d1, i_d0:i_d1) + s%Neighbors(index)%t0i(1:nOrb,1:nOrb,i)*exp(zi*dot_product(kp,s%Neighbors(index)%CellVector))
          !hk(j_u0:j_u1, i_u0:i_u1) = hk(j_u0:j_u1, i_u0:i_u1) + tmp(1:nOrb, 1:nOrb)
          !hk(j_d0:j_d1, i_d0:i_d1) = hk(j_d0:j_d1, i_d0:i_d1) + tmp(1:nOrb, 1:nOrb)
          current => current%next
        end do
      end do
    end do
  end do

  ! do i = 1, 18
  !   do j = 1, 18
  !     print *, i,j,hk(j,i)
  !   end do
  ! end do
  ! stop

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
  use mod_system,        only: s => sys
  use AtomTypes, only: NeighborIndex
  use TightBinding, only: nOrb
  use mod_SOC,    only: socscale, ls
  use mod_magnet,        only: lb, sb
  implicit none
  integer :: i, j, k
  integer :: i_u0, i_u1, i_d0, i_d1, j_u0, j_u1, j_d0, j_d1
  real(double), intent(in)  :: kp(3)
  complex(double) :: hee(2*nOrb, 2*nOrb, s%nAtoms)
  complex(double),dimension(s%nAtoms*2*nOrb,s%nAtoms*2*nOrb),intent(out)  :: hk,vsoc
  type(NeighborIndex), pointer :: current
  hk = zero
  vsoc = zero

  call U_matrix(hee, s%nAtoms, nOrb)

  ! Mouting slab hamiltonian
  do i=1,s%nAtoms
    i_u0 = (i-1)*2*nOrb + 1
    i_u1 = i_u0 + nOrb - 1
    i_d0 = i_u0 + nOrb
    i_d1 = i_d0 + nOrb - 1

    hk(i_u0:i_u1, i_u0:i_u1) = s%Types(s%Basis(i)%Material)%onSite ! + sb hee +socscale*lambda*ls
    hk(i_d0:i_d1, i_d0:i_d1) = s%Types(s%Basis(i)%Material)%onSite ! + sb hee +socscale*lambda*ls

    hk(i_u0:i_d1, i_u0:i_d1) = hk(i_u0:i_d1, i_u0:i_d1) &
                             + lb(:,:,i) + sb(:,:,i) + hee(:,:,i)
    vsoc(i_u0:i_d1, i_u0:i_d1) = socscale * s%Types(s%Basis(i)%Material)%Lambda * ls

    do j = 1, s%nAtoms
      j_u0 = (j-1) * 2 * nOrb + 1
      j_u1 = j_u0 + nOrb - 1
      j_d0 = j_u0 + nOrb
      j_d1 = j_d0 + nOrb - 1
      do k = 1, s%nStages
        current => s%Basis(i)%NeighborList(k,j)%head
        do while(associated(current))
          hk(j_u0:j_u1, i_u0:i_u1) = hk(j_u0:j_u1, i_u0:i_u1) + s%Neighbors(current%index)%t0i(:,:,i)*exp(zi*dot_product(kp,s%Neighbors(current%index)%CellVector))
          hk(j_d0:j_d1, i_d0:i_d1) = hk(j_d0:j_d1, i_d0:i_d1) + s%Neighbors(current%index)%t0i(:,:,i)*exp(zi*dot_product(kp,s%Neighbors(current%index)%CellVector))
          current => current%next
        end do
      end do
    end do
  end do
  return
end subroutine hamiltklinearsoc
