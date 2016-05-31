! Mounts the hoppings between orbitals for plane sites, 1st n.n. (with SOC)
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
subroutine helphbccsoc(kp,h00so,h01so,h10so,ly)
  use mod_f90_kind
  use mod_constants
  use mod_lattice
  use mod_tight_binding, only: t00, t01
  use mod_parameters, only: Npl
  implicit none
  integer     :: i
  integer,      intent(in)  :: ly
  real(double), intent(in)  :: kp(3)
  real(double)  :: kr0,kr1,kr2
  complex(double) :: expikr0,expikr1,expikr2
  complex(double),dimension(9,9)    :: h00,h01
  complex(double),dimension(18,18), intent(out) :: h00so,h01so,h10so

! in-plane hoppings h00
! on-site
  h00  = t00(ly,0,:,:)
! First and second nearest neighbors
  do i=1,n0
!   exponential of i.k.r0(i); i=1,n0
    kr0   = (kp(1)*r0(i,1))+(kp(2)*r0(i,2))+(kp(3)*r0(i,3))
    expikr0 = exp(zi*kr0)
    h00 = h00 + (expikr0*t00(ly,i,:,:))
  end do
  h00so = zero
  h00so(1:9,1:9)      = h00
  h00so(10:18,10:18)  = h00

! inter-plane hoppings h01
  if(ly.gt.Npl+1) then
    h01so = zero
    h10so = zero
    return
  end if
! Initializing h01
  h01  = zero
! First nearest neighbors
  do i=1,n1
!   exp(i.k.r1(i)) ; i=1,n1
    kr1   = (kp(1)*r1(i,1)) + (kp(2)*r1(i,2)) + (kp(3)*r1(i,3))
    expikr1 = exp(zi*kr1)
    h01   = h01 + (expikr1*t01(ly,i,:,:))
  end do
! Second nearest neighbors
  do i=1,n2
!   exp(i.k.r2(i)) ; i=1,n2
    kr2   = (kp(1)*r2(i,1)) + (kp(2)*r2(i,2)) + (kp(3)*r2(i,3))
    expikr2 = exp(zi*kr2)
    h01   = h01 + (expikr2*t01(ly,n1+i,:,:))
  end do
  h01so = zero
  h01so(1:9,1:9)      = h01
  h01so(10:18,10:18)  = h01
  h10so = transpose(conjg(h01so))

  return
end subroutine helphbccsoc



! Mounts the hoppings between orbitals for plane sites, 1st n.n. and 2nd n.n.(with SOC)
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
!   h02 = q-parallel inter-plane hoppings
subroutine helphfccsoc(kp,h00so,h01so,h10so,h02so,h20so,ly)
  use mod_f90_kind
  use mod_constants
  use mod_lattice
  use mod_tight_binding, only: t00, t01, t02
  use mod_parameters, only: Npl
  implicit none
  integer     :: i
  integer,      intent(in)  :: ly
  real(double), intent(in)  :: kp(3)
  real(double)  :: kr0,kr1,kr2
  complex(double) :: expikr0,expikr1,expikr2
  complex(double),dimension(9,9)    :: h00,h01,h02
  complex(double),dimension(18,18), intent(out) :: h00so,h01so,h10so,h02so,h20so

! in-plane hoppings h00
! on-site
  h00  = t00(ly,0,:,:)
! First and second nearest neighbors
  do i=1,n0
!   exponential of i.k.r0(i); i=1,n0
    kr0   = (kp(1)*r0(i,1))+(kp(2)*r0(i,2))+(kp(3)*r0(i,3))
    expikr0 = exp(zi*kr0)
    h00 = h00 + (expikr0*t00(ly,i,:,:))
  end do
  h00so = zero
  h00so(1:9,1:9)      = h00
  h00so(10:18,10:18)  = h00

! inter-plane hoppings h01
  if(ly.gt.Npl+1) then
    h01so = zero
    h10so = zero
    return
  end if
! Initializing h01
  h01  = zero
! First nearest neighbors
  do i=1,n1
!   exp(i.k.r1(i)) ; i=1,n1
    kr1   = (kp(1)*r1(i,1)) + (kp(2)*r1(i,2)) + (kp(3)*r1(i,3))
    expikr1 = exp(zi*kr1)
    h01   = h01 + (expikr1*t01(ly,i,:,:))
  end do
  h01so = zero
  h01so(1:9,1:9)      = h01
  h01so(10:18,10:18)  = h01
  h10so = transpose(conjg(h01so))

! inter-plane hoppings h02
  if(ly.gt.Npl) then
    h02so = zero
    h20so = zero
    return
  end if
! Initializing h02
  h02  = zero
! Second nearest neighbors
  do i=1,n2
!   exp(i.k.r2(i)) ; i=1,n2
    kr2   = (kp(1)*r2(i,1)) + (kp(2)*r2(i,2)) + (kp(3)*r2(i,3))
    expikr2 = exp(zi*kr2)
    h02   = h02 + (expikr2*t02(ly,i,:,:))
  end do
  h02so = zero
  h02so(1:9,1:9)      = h02
  h02so(10:18,10:18)  = h02
  h20so = transpose(conjg(h02so))

  return
end subroutine helphfccsoc