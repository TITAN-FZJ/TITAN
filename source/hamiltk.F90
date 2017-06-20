! Calculate hamiltonian of a slab containing
! Npl layers + 1 layer of Empty spheres on each side
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!  1     2     3     4      Npl-1   Npl  Npl+1  Npl+2
!         <-S-> <S-1>           <S-1> <-S->
subroutine hamiltk(kp,hk)
  use mod_tight_binding, only: lambda, ls, t0, t0i
  use mod_parameters,    only: Npl_total, socscale
  use mod_constants,     only: zi, zero
  use mod_f90_kind,      only: double
  use mod_system,        only: npln, r_nn, l_nn
  use mod_magnet,        only: lb, sb
  implicit none
  integer :: i, j, l, loc_pln
  integer :: i0, i1, j0, j1
  real(double), intent(in) :: kp(3)
  complex(double) :: hee(18,18,Npl_total)
  complex(double), dimension((Npl_total)*18,(Npl_total)*18), intent(out) :: hk

  hk = zero

  call U_matrix(hee)

  ! Mouting slab hamiltonian
  do i=1, Npl_total
     i0 = (i-1)*18+1
     i1 = i0+9

     loc_pln = npln
     if(Npl_total < i + npln) loc_pln = Npl_total - i + 1

     hk(i0:i0+8,i0:i0+8) = t0(:,:,i) ! + sb hee +socscale*lambda*ls
     hk(i1:i1+8,i1:i1+8) = t0(:,:,i) ! + sb hee +socscale*lambda*ls

     hk(i0:i0+17, i0:i0+17) = hk(i0:i0+17, i0:i0+17) + lb(:,:,i) + sb(:,:,i) + hee(:,:,i) + (socscale*lambda(i)*ls)

     do j = 1, loc_pln
        j0 = i0 + 18 * (j-1)
        j1 = i0 + 18 * (j-1) + 9
        do l = l_nn(1,j), l_nn(1,j+1)-1
          hk(j0:j0+8, i0:i0+8) = hk(j0:j0+8, i0:i0+8) + t0i(:,:,l,i)*exp(zi*dot_product(kp,r_nn(:,l)))
          hk(j1:j1+8, i1:i1+8) = hk(j1:j1+8, i1:i1+8) + t0i(:,:,l,i)*exp(zi*dot_product(kp,r_nn(:,l)))
        end do
        if(j > 1) hk(i0:i0+17, j0:j0+17) = transpose(conjg(hk(j0:j0+17, i0:i0+17)))
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
  use mod_tight_binding, only: lambda, ls, t0, t0i
  use mod_parameters,    only: socscale, Npl_total
  use mod_constants,     only: zero, zi
  use mod_f90_kind,      only: double
  use mod_system,        only: npln, r_nn, l_nn
  use mod_magnet,        only: lb, sb
  implicit none
  integer :: i,j,l, loc_pln
  integer :: i0,i1,j0,j1
  real(double), intent(in)  :: kp(3)
  complex(double) :: hee(18,18,Npl_total)
  complex(double),dimension((Npl_total)*18,(Npl_total)*18),intent(out)  :: hk,vsoc

  hk = zero
  vsoc = zero

  call U_matrix(hee)

  ! Mouting slab hamiltonian
  do i=1,Npl_total
     i0 = (i-1)*18+1
     i1 = i0+9

     loc_pln = npln
     if( Npl_total < i + npln ) loc_pln = Npl_total - i + 1

     hk(i0:i0+8,i0:i0+8) = t0(:,:,i) ! + sb hee +socscale*lambda*ls
     hk(i1:i1+8,i1:i1+8) = t0(:,:,i) ! + sb hee +socscale*lambda*ls

     hk  (i0:i0+17, i0:i0+17) = hk(i0:i0+17, i0:i0+17) + lb(:,:,i) + sb(:,:,i) + hee(:,:,i)
     vsoc(i0:i0+17, i0:i0+17) = socscale*lambda(i)*ls

     do j = 1, loc_pln
        j0 = i0 + 18 * (j-1)
        j1 = i0 + 18 * (j-1) + 9
        do l = l_nn(1,j), l_nn(1,j+1)-1
           hk(j0:j0+8, i0:i0+8) = hk(j0:j0+8, i0:i0+8) + t0i(:,:,l,i)*exp(zi*dot_product(kp,r_nn(:,l)))
           hk(j1:j1+8, i1:i1+8) = hk(j1:j1+8, i1:i1+8) + t0i(:,:,l,i)*exp(zi*dot_product(kp,r_nn(:,l)))
        end do
        if(j > 1) hk(i0:i0+17, j0:j0+17) = transpose(conjg(hk(j0:j0+17, i0:i0+17)))
     end do
  end do
  return
end subroutine hamiltklinearsoc
