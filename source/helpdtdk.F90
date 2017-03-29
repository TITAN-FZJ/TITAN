! Mounts the derivative of hoppings between
! orbitals for plane sites, 1st n.n.
! input:
!   kp = in-plane wave vector k
!   ly = layer
!
! in common:
!   pi, sq2, sq3, zi, zero
!   n0    = # of 1st. and 2nd. in-plane n.n. (n01+n02)
!   n1    = # of 1st. inter-plane n.n.
!   n2    = # of 2nd. inter-plane n.n.
!   r0(i,l) = coord. of the in-plane 1st. and 2nd n.n. (i=1,n0; l = 1,3 (x,y,z));
!   c0(i,l) = direction cosines of r0
!   r1(i,l) = coord. of the inter-plane 1st. n.n. (i=1,n1)
!   c1(i,l) = direction cosines of r1
!   r2(i,l) = coord. of the inter-plane 2nd. n.n. (i=1,n2)
!   c2(i,l) = direction cosines of r2
!
! output:
!   h00 = q-parallel in-plane hoppings
!   h01 = q-parallel inter-plane hoppings
subroutine helpdtdk(kp,h00,h01,h10,ly)
  use mod_f90_kind
  use mod_parameters, only: dirEfieldvec,lattice
  use mod_constants
  use mod_lattice
  use mod_tight_binding
  implicit none
  integer         :: i
  integer,      intent(in)  :: ly
  real(double), intent(in)  :: kp(3)
  real(double)    :: kr0,kr1,kr2,dirEr0,dirEr1,dirEr2
  complex(double) :: expikr0,expikr1,expikr2
  complex(double),dimension(9,9), intent(out)   :: h00,h01,h10

! derivative of in-plane of hoppings h00
  h00  = zero
! First nearest neighbors
  do i=1,n0
!   exponential of i.k.r0(i); i=1,n0
    kr0   = dot_product(kp, r0(i,:))
    expikr0 = exp(zi*kr0)
    ! term from derivative related to the direction of applied field i.u.r0(i); i=1,n0
    dirEr0   = dot_product(dirEfieldvec, r0(i,:))
    h00 = h00 + (zi*dirEr0*expikr0*t00(ly,i,:,:))
  end do

! inter-plane hoppings h01 and h10
  h01  = zero
! First nearest neighbors
  do i=1,n1
!   exp(i.k.r1(i)) ; i=1,n1
    kr1   = dot_product(kp, r1(i,:))
    expikr1 = exp(zi*kr1)
    ! term from derivative related to the direction of applied field i.u.r1(i); i=1,n1
    dirEr1  = dot_product(dirEfieldvec, r1(i,:))
    h01   = h01 + (zi*dirEr1*expikr1*t01(ly,i,:,:))
  end do

  ! Second nearest neighbors
  do i=1,n2
  !   exp(i.k.r2(i)) ; i=1,n2
    kr2   = dot_product(kp, r2(i,:))
    expikr2 = exp(zi*kr2)
    ! term from derivative related to the direction of applied field i.u.r2(i); i=1,n2
    dirEr2  = dot_product(dirEfieldvec, r2(i,:))
    h01   = h01 + (zi*dirEr2*expikr2*t01(ly,n1+i,:,:))
  end do
  h10 = transpose(conjg(h01))
  ! NOTE: The derivative of h02 is always zero for fcc because
  ! dirEfield is in-plane and r2 is perpendicular to the plane

  return
end subroutine helpdtdk
