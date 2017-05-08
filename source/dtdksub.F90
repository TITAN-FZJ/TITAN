! Calculate the derivative of H0 of a slab containing
! Npl layers (ES not included)
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!        1     2     3       Npl-2 Npl-1  Npl
!         <-S-> <S-1>           <S-1> <-S->
subroutine dtdksub(kp,dtdk)
  use mod_system,        only: r_nn, l_nn, npln
  use mod_f90_kind,      only: double
  use mod_constants,     only: zero, zi
  use mod_parameters,    only: Npl, dirEfieldvec
  use mod_tight_binding, only: t0i, tbmode
  implicit none
  integer     :: i, j, l, loc_pln, loc_layers, offset
  real(double), intent(in)  :: kp(3)
  complex(double),dimension(Npl,Npl,9,9),intent(out)  :: dtdk

  offset = 0
  if(tbmode == 2) offset = 1

  dtdk = zero

  ! Mouting derivative of slab's hamiltonian

  do i=1, Npl
     loc_pln = npln
     if( Npl < i + npln ) loc_pln = Npl - i

     do j = 1, loc_pln
        do l = l_nn(1,j), l_nn(1,j+1)-1
           dtdk(i,i+j-1,:,:) = dtdk(i,i+j-1,:,:) + zi * dot_product(dirEfieldvec, r_nn(:, l)) * t0i(i+offset,l,:,:) * exp(zi * dot_product(kp, r_nn(:,l)))
        end do
        if(j > 1) dtdk(i+j-1,i,:,:) = transpose(conjg(dtdk(i,i+j-1,:,:)))
     end do
  end do
  return
end subroutine dtdksub
