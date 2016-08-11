module mod_lattice
  use mod_f90_kind
  implicit none
! n01 - number of in-plane 1st nearest neighbors
! n02 - number of in-plane 2nd nearest neighbors
! n1 - number of inter-plane 1st nearest neighbors
! n2 - number of inter-plane 2nd nearest neighbors
! r0 - position of in-plane 1st and 2nd nearest neighbors
! c0 - direction cosines of in-plane 1st and 2nd nearest neighbors
! r1 - position of inter-plane 1st nearest neighbors
! c1 - direction cosines of inter-plane 1st nearest neighbors
! r2 - position of inter-plane 1st nearest neighbors
! c2 - direction cosines of inter-plane 1st nearest neighbors
  integer :: n01,n02,n0,n1,n2
! plnn - number of nearest neighbors planes - to mount (tri-)diagonal matrices
  integer :: plnn
  real(double), dimension(:,:), allocatable :: r0,c0,r1,c1,r2,c2

contains

! BCC(110)
! To calculate the sum over the BZ,
! the cubic axis are rotated are chosen such as
!   z=z'|
!       |  /y'
!       | /\ a
!       |/__|________y
!       /\
!      /__\
!    x/ a  \x'
!
!   a = pi/4
! The direction perpendicular to the layers is y'
!
! n01,n02   = # of 1st. and 2nd. in-plane n.n. respectively
! n1,n2     = # of 1st. and 2nd. inter-plane n.n. respectively
! r0( i, l) = coord. of the in-plane (first & second n.n.)
!       (i=1,n0; l = 1,3 (x,y,z);
!       a0 = latt. const.
! c0( i, l) = direction cosines of r0
! r1( i, l) = coord. of the inter-plane 1st. n.n. (i=1,n1)
! c1( i, l) = direction cosines of r1
! r2( i, l) = coord. of the inter-plane 2nd. n.n. (i=1,n2)
! c2( i, l) = direction cosines of r2

  subroutine bcc110()
    use mod_f90_kind
    use mod_constants
    use mod_parameters, only: a0
    real(double)  :: aux

    n01=4
    n02=2
    n0=n01+n02
    plnn=1
    n1=2
    n2=2

    allocate(r0(n0,3),c0(n0,3),r1(n1,3),c1(n1,3),r2(n2,3),c2(n2,3))

  ! BCC (1 -1 0)
  ! ! in-plane 1st. and 2nd. n.n.
    r0(1,1) = 0.5d0
    r0(1,2) = 0.5d0
    r0(1,3) = 0.5d0

    r0(2,1) =-0.5d0
    r0(2,2) =-0.5d0
    r0(2,3) = 0.5d0

    r0(3,1) =-0.5d0
    r0(3,2) =-0.5d0
    r0(3,3) =-0.5d0

    r0(4,1) = 0.5d0
    r0(4,2) = 0.5d0
    r0(4,3) =-0.5d0

    r0(5,1) = 0.d0
    r0(5,2) = 0.d0
    r0(5,3) = 1.d0

    r0(6,1) = 0.d0
    r0(6,2) = 0.d0
    r0(6,3) =-1.d0

    r0 = r0*a0

    aux = 1.d0/sq3

    c0(1,1) = aux
    c0(1,2) = aux
    c0(1,3) = aux

    c0(2,1) =-aux
    c0(2,2) =-aux
    c0(2,3) = aux

    c0(3,1) =-aux
    c0(3,2) =-aux
    c0(3,3) =-aux

    c0(4,1) = aux
    c0(4,2) = aux
    c0(4,3) =-aux

    c0(5,1) = 0.d0
    c0(5,2) = 0.d0
    c0(5,3) = 1.d0

    c0(6,1) = 0.d0
    c0(6,2) = 0.d0
    c0(6,3) =-1.d0

  !  inter-plane 1st. n.n.

    r1(1,1) =-0.5d0
    r1(1,2) = 0.5d0
    r1(1,3) = 0.5d0

    r1(2,1) =-0.5d0
    r1(2,2) = 0.5d0
    r1(2,3) =-0.5d0

    r1 = r1*a0

    c1(1,1) =-aux
    c1(1,2) = aux
    c1(1,3) = aux

    c1(2,1) =-aux
    c1(2,2) = aux
    c1(2,3) =-aux

  !  inter-plane 2nd. n.n.

    r2(1,1) = 0.d0
    r2(1,2) = 1.d0
    r2(1,3) = 0.d0

    r2(2,1) =-1.d0
    r2(2,2) = 0.d0
    r2(2,3) = 0.d0

    r2 = r2*a0

    c2(1,1) = 0.d0
    c2(1,2) = 1.d0
    c2(1,3) = 0.d0

    c2(2,1) =-1.d0
    c2(2,2) = 0.d0
    c2(2,3) = 0.d0

    return
  end subroutine bcc110

! FCC(100)
! The direction perpendicular to the layers is z'
!
! n01,n02   = # of 1st. and 2nd. in-plane n.n. respectively
! n1,n2     = # of 1st. and 2nd. inter-plane n.n. respectively
! r0( i, l) = coord. of the in-plane (first & second n.n.)
!       (i=1,n0; l = 1,3 (x,y,z);
!       a0 = latt. const.
! c0( i, l) = direction cosines of r0
! r1( i, l) = coord. of the inter-plane 1st. n.n. (i=1,n1)
! c1( i, l) = direction cosines of r1
! r2( i, l) = coord. of the inter-plane 2nd. n.n. (i=1,n2)
! c2( i, l) = direction cosines of r2

  subroutine fcc100()
    use mod_f90_kind
    use mod_constants
    use mod_parameters, only: a0
    real(double)  :: aux

    n01=4
    n02=4
    n0=n01+n02
    plnn=2
    n1=4
    n2=1

    allocate(r0(n0,3),c0(n0,3),r1(n1,3),c1(n1,3),r2(n2,3),c2(n2,3))

  ! FCC (1 0 0)
  ! ! in-plane 1st. n.n.
    r0(1,1) = 0.5d0
    r0(1,2) = 0.5d0
    r0(1,3) = 0.0d0

    r0(2,1) =-0.5d0
    r0(2,2) = 0.5d0
    r0(2,3) = 0.0d0

    r0(3,1) =-0.5d0
    r0(3,2) =-0.5d0
    r0(3,3) = 0.0d0

    r0(4,1) = 0.5d0
    r0(4,2) =-0.5d0
    r0(4,3) = 0.0d0

    r0(5,1) = 1.0d0
    r0(5,2) = 0.0d0
    r0(5,3) = 0.0d0

    r0(6,1) = 0.0d0
    r0(6,2) = 1.0d0
    r0(6,3) = 0.0d0

    r0(7,1) =-1.0d0
    r0(7,2) = 0.0d0
    r0(7,3) = 0.0d0

    r0(8,1) = 0.0d0
    r0(8,2) =-1.0d0
    r0(8,3) = 0.0d0

    aux = sq2
    c0(1:4,:) = r0(1:4,:)*aux
    c0(5:8,:) = r0(5:8,:)

    r0 = r0*a0

  !  inter-plane 1st. n.n.
    r1(1,1) = 0.5d0
    r1(1,2) = 0.0d0
    r1(1,3) = 0.5d0

    r1(2,1) = 0.0d0
    r1(2,2) = 0.5d0
    r1(2,3) = 0.5d0

    r1(3,1) =-0.5d0
    r1(3,2) = 0.0d0
    r1(3,3) = 0.5d0

    r1(4,1) = 0.0d0
    r1(4,2) =-0.5d0
    r1(4,3) = 0.5d0

    c1 = r1*aux

    r1 = r1*a0

  !  inter-plane 2nd. n.n.
    r2(1,1) = 0.d0
    r2(1,2) = 0.d0
    r2(1,3) = 1.d0

    c2 = r2

    r2 = r2*a0

    return
  end subroutine fcc100

! FCC(111)
! The direction perpendicular to the layers is z'
!
! n01,n02   = # of 1st. and 2nd. in-plane n.n. respectively
! n1,n2     = # of 1st. and 2nd. inter-plane n.n. respectively
! r0( i, l) = coord. of the in-plane (first & second n.n.)
!       (i=1,n0; l = 1,3 (x,y,z);
!       a0 = latt. const.
! c0( i, l) = direction cosines of r0
! r1( i, l) = coord. of the inter-plane 1st. n.n. (i=1,n1)
! c1( i, l) = direction cosines of r1
! r2( i, l) = coord. of the inter-plane 2nd. n.n. (i=1,n2)
! c2( i, l) = direction cosines of r2

  subroutine fcc111()
    use mod_f90_kind
    use mod_constants
    use mod_parameters, only: a0
    real(double)  :: aux

    n01=6
    n02=0
    n0=n01+n02
    plnn=1
    n1=3
    n2=3

    allocate(r0(n0,3),c0(n0,3),r1(n1,3),c1(n1,3),r2(n2,3),c2(n2,3))

  ! FCC (1 1 1)
  ! ! in-plane 1st. n.n. (there are no 2nd. n.n.)
    r0(1,1) =-0.5d0
    r0(1,2) = 0.5d0
    r0(1,3) = 0.0d0

    r0(2,1) =-0.5d0
    r0(2,2) = 0.0d0
    r0(2,3) = 0.5d0

    r0(3,1) = 0.0d0
    r0(3,2) =-0.5d0
    r0(3,3) = 0.5d0

    r0(4,1) = 0.5d0
    r0(4,2) =-0.5d0
    r0(4,3) = 0.0d0

    r0(5,1) = 0.5d0
    r0(5,2) = 0.0d0
    r0(5,3) =-0.5d0

    r0(6,1) = 0.0d0
    r0(6,2) = 0.5d0
    r0(6,3) =-0.5d0

    aux = sq2
    c0 = aux*r0

    r0 = r0*a0

  !  inter-plane 1st. n.n.
    r1(1,1) = 0.0d0
    r1(1,2) = 0.5d0
    r1(1,3) = 0.5d0

    r1(2,1) = 0.5d0
    r1(2,2) = 0.0d0
    r1(2,3) = 0.5d0

    r1(3,1) = 0.5d0
    r1(3,2) = 0.5d0
    r1(3,3) = 0.0d0

    c1 = aux*r1

    r1 = r1*a0

  !  inter-plane 2nd. n.n.
    r2(1,1) = 0.d0
    r2(1,2) = 0.d0
    r2(1,3) = 1.d0

    r2(2,1) = 1.d0
    r2(2,2) = 0.d0
    r2(2,3) = 0.d0

    r2(3,1) = 0.d0
    r2(3,2) = 1.d0
    r2(3,3) = 0.d0

    c2 = r2

    r2 = r2*a0

    return
  end subroutine fcc111

end module mod_lattice