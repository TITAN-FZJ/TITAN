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

! Out-of-plane unit vector
  real(double), dimension(3) :: versor_oop
! In-plane unit vector perpendicular to current direction
  integer :: transverse_neighbors
  integer :: longitudinal_neighbors
  real(double), dimension(3) :: versor_Eperp

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
    implicit none
    real(double)  :: aux
    real(double) :: a1(3)
    real(double) :: a2(3)
    real(double) :: a3(3)
    real(double) :: p(3)

    p  = [ 1.0d0, -1.0d0,  0.0d0] / sqrt(2.d0)
    a1 = [-0.5d0,  0.5d0,  0.5d0] * a0
    a2 = [ 0.5d0, -0.5d0,  0.5d0] * a0
    a3 = [ 0.5d0,  0.5d0, -0.5d0] * a0

    call next_neighbour_init(a1, a2, a3, p)

    ! n01=4
    ! n02=2
    ! n0=n01+n02
    ! plnn=1
    ! n1=2
    ! n2=2
    !
    ! allocate(r0(n0,3),c0(n0,3),r1(n1,3),c1(n1,3),r2(n2,3),c2(n2,3))
    !
    ! ! BCC (1 -1 0)
    ! ! ! in-plane 1st. and 2nd. n.n.
    ! r0(1,1) = 0.5d0
    ! r0(1,2) = 0.5d0
    ! r0(1,3) = 0.5d0
    !
    ! r0(2,1) =-0.5d0
    ! r0(2,2) =-0.5d0
    ! r0(2,3) = 0.5d0
    !
    ! r0(3,1) =-0.5d0
    ! r0(3,2) =-0.5d0
    ! r0(3,3) =-0.5d0
    !
    ! r0(4,1) = 0.5d0
    ! r0(4,2) = 0.5d0
    ! r0(4,3) =-0.5d0
    !
    ! r0(5,1) = 0.d0
    ! r0(5,2) = 0.d0
    ! r0(5,3) = 1.d0
    !
    ! r0(6,1) = 0.d0
    ! r0(6,2) = 0.d0
    ! r0(6,3) =-1.d0
    !
    ! r0 = r0*a0
    !
    ! aux = 1.d0/sq3
    !
    ! c0(1,1) = aux
    ! c0(1,2) = aux
    ! c0(1,3) = aux
    !
    ! c0(2,1) =-aux
    ! c0(2,2) =-aux
    ! c0(2,3) = aux
    !
    ! c0(3,1) =-aux
    ! c0(3,2) =-aux
    ! c0(3,3) =-aux
    !
    ! c0(4,1) = aux
    ! c0(4,2) = aux
    ! c0(4,3) =-aux
    !
    ! c0(5,1) = 0.d0
    ! c0(5,2) = 0.d0
    ! c0(5,3) = 1.d0
    !
    ! c0(6,1) = 0.d0
    ! c0(6,2) = 0.d0
    ! c0(6,3) =-1.d0
    !
    ! !  inter-plane 1st. n.n.
    !
    ! r1(1,1) =-0.5d0
    ! r1(1,2) = 0.5d0
    ! r1(1,3) = 0.5d0
    !
    ! r1(2,1) =-0.5d0
    ! r1(2,2) = 0.5d0
    ! r1(2,3) =-0.5d0
    !
    ! r1 = r1*a0
    !
    ! c1(1,1) =-aux
    ! c1(1,2) = aux
    ! c1(1,3) = aux
    !
    ! c1(2,1) =-aux
    ! c1(2,2) = aux
    ! c1(2,3) =-aux
    !
    ! !  inter-plane 2nd. n.n.
    !
    ! r2(1,1) = 0.d0
    ! r2(1,2) = 1.d0
    ! r2(1,3) = 0.d0
    !
    ! r2(2,1) =-1.d0
    ! r2(2,2) = 0.d0
    ! r2(2,3) = 0.d0
    !
    ! r2 = r2*a0
    !
    ! c2(1,1) = 0.d0
    ! c2(1,2) = 1.d0
    ! c2(1,3) = 0.d0
    !
    ! c2(2,1) =-1.d0
    ! c2(2,2) = 0.d0
    ! c2(2,3) = 0.d0

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
    implicit none
    real(double)  :: aux
    real(double) :: a1(3)
    real(double) :: a2(3)
    real(double) :: a3(3)
    real(double) :: p(3)

    p  = [ 0.0d0,  0.0d0,  -1.0d0]
    a1 = [ 0.0d0,  0.5d0,  0.5d0] * a0
    a2 = [ 0.5d0,  0.0d0,  0.5d0] * a0
    a3 = [ 0.5d0,  0.5d0,  0.0d0] * a0

    call next_neighbour_init(a1, a2, a3, p)

    ! n01=4
    ! n02=4
    ! n0=n01+n02
    ! plnn=2
    ! n1=4
    ! n2=1
    !
    ! allocate(r0(n0,3),c0(n0,3),r1(n1,3),c1(n1,3),r2(n2,3),c2(n2,3))
    !
    ! ! FCC (1 0 0)
    ! ! ! in-plane 1st. n.n.
    ! r0(1,1) = 0.5d0
    ! r0(1,2) = 0.5d0
    ! r0(1,3) = 0.0d0
    !
    ! r0(2,1) =-0.5d0
    ! r0(2,2) = 0.5d0
    ! r0(2,3) = 0.0d0
    !
    ! r0(3,1) =-0.5d0
    ! r0(3,2) =-0.5d0
    ! r0(3,3) = 0.0d0
    !
    ! r0(4,1) = 0.5d0
    ! r0(4,2) =-0.5d0
    ! r0(4,3) = 0.0d0
    !
    ! r0(5,1) = 1.0d0
    ! r0(5,2) = 0.0d0
    ! r0(5,3) = 0.0d0
    !
    ! r0(6,1) = 0.0d0
    ! r0(6,2) = 1.0d0
    ! r0(6,3) = 0.0d0
    !
    ! r0(7,1) =-1.0d0
    ! r0(7,2) = 0.0d0
    ! r0(7,3) = 0.0d0
    !
    ! r0(8,1) = 0.0d0
    ! r0(8,2) =-1.0d0
    ! r0(8,3) = 0.0d0
    !
    ! aux = sq2
    ! c0(1:4,:) = r0(1:4,:)*aux
    ! c0(5:8,:) = r0(5:8,:)
    !
    ! r0 = r0*a0
    !
    ! !  inter-plane 1st. n.n.
    ! r1(1,1) = 0.5d0
    ! r1(1,2) = 0.0d0
    ! r1(1,3) = 0.5d0
    !
    ! r1(2,1) = 0.0d0
    ! r1(2,2) = 0.5d0
    ! r1(2,3) = 0.5d0
    !
    ! r1(3,1) =-0.5d0
    ! r1(3,2) = 0.0d0
    ! r1(3,3) = 0.5d0
    !
    ! r1(4,1) = 0.0d0
    ! r1(4,2) =-0.5d0
    ! r1(4,3) = 0.5d0
    !
    ! c1 = r1*aux
    !
    ! r1 = r1*a0
    !
    ! !  inter-plane 2nd. n.n.
    ! r2(1,1) = 0.d0
    ! r2(1,2) = 0.d0
    ! r2(1,3) = 1.d0
    !
    ! c2 = r2
    !
    ! r2 = r2*a0

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
    implicit none
    real(double)  :: aux
    real(double) :: a1(3)
    real(double) :: a2(3)
    real(double) :: a3(3)
    real(double) :: p(3)

    p  = -1.d0 * [ 1.0d0,  1.0d0,  1.0d0] / sqrt(3.d0)
    a1 = [ 0.0d0,  0.5d0,  0.5d0] * a0
    a2 = [ 0.5d0,  0.0d0,  0.5d0] * a0
    a3 = [ 0.5d0,  0.5d0,  0.0d0] * a0

    call next_neighbour_init(a1, a2, a3, p)

    ! n01=6
    ! n02=0
    ! n0=n01+n02
    ! plnn=1
    ! n1=3
    ! n2=3
    !
    ! allocate(r0(n0,3),c0(n0,3),r1(n1,3),c1(n1,3),r2(n2,3),c2(n2,3))
    !
    ! ! FCC (1 1 1)
    ! ! ! in-plane 1st. n.n. (there are no 2nd. n.n.)
    ! r0(1,1) =-0.5d0
    ! r0(1,2) = 0.5d0
    ! r0(1,3) = 0.0d0
    !
    ! r0(2,1) =-0.5d0
    ! r0(2,2) = 0.0d0
    ! r0(2,3) = 0.5d0
    !
    ! r0(3,1) = 0.0d0
    ! r0(3,2) =-0.5d0
    ! r0(3,3) = 0.5d0
    !
    ! r0(4,1) = 0.5d0
    ! r0(4,2) =-0.5d0
    ! r0(4,3) = 0.0d0
    !
    ! r0(5,1) = 0.5d0
    ! r0(5,2) = 0.0d0
    ! r0(5,3) =-0.5d0
    !
    ! r0(6,1) = 0.0d0
    ! r0(6,2) = 0.5d0
    ! r0(6,3) =-0.5d0
    !
    ! aux = sq2
    ! c0 = aux*r0
    !
    ! r0 = r0*a0
    !
    ! !  inter-plane 1st. n.n.
    ! r1(1,1) = 0.0d0
    ! r1(1,2) = 0.5d0
    ! r1(1,3) = 0.5d0
    !
    ! r1(2,1) = 0.5d0
    ! r1(2,2) = 0.0d0
    ! r1(2,3) = 0.5d0
    !
    ! r1(3,1) = 0.5d0
    ! r1(3,2) = 0.5d0
    ! r1(3,3) = 0.0d0
    !
    ! c1 = aux*r1
    !
    ! r1 = r1*a0
    !
    ! !  inter-plane 2nd. n.n.
    ! r2(1,1) = 0.d0
    ! r2(1,2) = 0.d0
    ! r2(1,3) = 1.d0
    !
    ! r2(2,1) = 1.d0
    ! r2(2,2) = 0.d0
    ! r2(2,3) = 0.d0
    !
    ! r2(3,1) = 0.d0
    ! r2(3,2) = 1.d0
    ! r2(3,3) = 0.d0
    !
    ! c2 = r2
    !
    ! r2 = r2*a0

    return
  end subroutine fcc111

  subroutine next_neighbour_init(a1,a2,a3,p)
    use mod_f90_kind
    use mod_constants
    implicit none
    real(double), intent(in) :: a1(3)
    real(double), intent(in) :: a2(3)
    real(double), intent(in) :: a3(3)
    real(double), intent(in) :: p(3)
    real(double), allocatable :: tmp(:,:,:)          !> @var Array contains all neighbours in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real(double), allocatable :: dist(:)             !> @var Contains all distances in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real(double), allocatable :: on_plane(:,:,:)     !> @var Array of type (3,2,on_cnt), all coordinates and directional cosines of neighbours in plane described by vector
    real(double), allocatable :: off_plane(:,:,:)    !> @var Array of type (3,2,off_cnt), containing half the intra-plane elements
    real(double), allocatable :: on_plane_dist(:)    !> @var Contains the distance of each point to the Origin
    real(double), allocatable :: off_plane_dist(:,:) !> @var Contains the distance of each point to the Origin and to the plane
    real(double) :: dist_tmp
    real(double) :: pos_tmp(3)
    real(double) :: cos_tmp(3)
    real(double), allocatable :: cnt(:,:)            !> @var Contains distance and stage informations for all elements in tmp
    integer, allocatable :: on_plane_aux(:)          !> @var on_plane_aux(i) contains the first element of the i-1 -th neighbour stage
    integer, allocatable :: off_plane_aux(:)
    integer :: nn_stages, nn_size, nnt, on_cnt, off_cnt
    integer :: i, j, l, m, n
    logical :: pl_flag

    nn_stages = 2
    nnt = 2 * nn_stages

    allocate(tmp(3, 2, (2*nnt+1)**3))
    allocate(dist((2*nnt+1)**3), cnt(2, nn_stages+2))
    allocate(on_plane_aux(nn_stages+2), off_plane_aux(nn_stages+2))

    ! Generate neighbours up to 2 * nn_stages the number of wanted neighbours
    ! Insertion sort approach
    i = 1
    do l = -nnt, nnt
      do m = -nnt, nnt
        do n = -nnt, nnt
          pos_tmp = l*a1 + m*a2 + n*a3
          cos_tmp = pos_tmp / sqrt(pos_tmp(1)**2 + pos_tmp(2)**2 + pos_tmp(3)**2)
          dist_tmp = sqrt(pos_tmp(1)**2 + pos_tmp(2)**2 + pos_tmp(3)**2)
          j = i-1
          do while( j >= 1 .and. dist(j) > dist_tmp)
            tmp(:,:, j+1) = tmp(:,:,j)
            dist(j+1) = dist(j)
            j = j-1
          end do
          tmp(:, 1, j+1) = pos_tmp
          tmp(:, 2, j+1) = cos_tmp
          dist(j+1) = dist_tmp
          i = i+1
        end do
      end do
    end do

    ! Count member of each neighbour stage
    ! cnt(1,i) contains the last member of stage i
    ! cnt(2,i) contains the actual distance of the i-th neighbours
    i = 1
    cnt(1,1) = 1.d0
    do j=1, nn_stages+1
      cnt(2,j) = dist(i)
      do while(cnt(2,j) == dist(i))
        i = i+1
      end do
      cnt(1,j+1) = i
    end do

    nn_size = cnt(1, nn_stages+2) - 1

    ! Sort into on and off plane of interest
    allocate(on_plane(3,2,nn_size), off_plane(3,2,nn_size))
    allocate(on_plane_dist(nn_size), off_plane_dist(2, nn_size))

    on_cnt = 0
    off_cnt = 0
    on_plane_aux = 1
    off_plane_aux = 1

    do i=1, nn_size
      dist_tmp = p(1)*tmp(1,1,i) + p(2)*tmp(2,1,i) + p(3)*tmp(3,1,i)
      if(abs(dist_tmp) < 1d-9) then
        on_cnt = on_cnt + 1
        on_plane(:,:, on_cnt) = tmp(:,:, i)
        on_plane_dist(on_cnt) = dist(i)
        do j=1, nn_stages+1
          if( dist(i) <= cnt(2,j) ) then
            on_plane_aux(j+1) = on_plane_aux(j+1) + 1
          end if
        end do
      else if(dist_tmp < 0) then
        off_cnt = off_cnt + 1
        off_plane(:,:, off_cnt) = tmp(:,:, i)
        off_plane_dist(1,off_cnt) = dist(i)
        off_plane_dist(2,off_cnt) = abs(p(1) * tmp(1,1,i) + p(2) * tmp(2,1,i) + p(3) * tmp(3,1,i))
        do j=1, nn_stages + 1
          if( dist(i) <= cnt(2,j) ) then
            off_plane_aux(j+1) = off_plane_aux(j+1) + 1
          end if
        end do
      end if
    end do

    ! Check number of planes
    plnn = 0
    do i=1, off_cnt
      pl_flag = .true.
      do j = i-1, 1, -1
        if(off_plane_dist(2,i) == off_plane_dist(2,j)) then
          pl_flag = .false.
        end if
      end do
      if(pl_flag) then
        plnn = plnn + 1
      end if
    end do

    n01 = on_plane_aux(3) - on_plane_aux(2)
    n02 = on_plane_aux(4) - on_plane_aux(3)
    n0  = n01 + n02
    n1  = off_plane_aux(3) - off_plane_aux(2)
    n2  = off_plane_aux(4) - off_plane_aux(3)

    allocate(r0(n0,3), c0(n0,3), r1(n1,3), c1(n1,3), r2(n2,3), c2(n2,3))

    do i=1,n0
      r0(i, :) = on_plane(:, 1, on_plane_aux(2) + i - 1)
      c0(i, :) = on_plane(:, 2, on_plane_aux(2) + i - 1)
    end do

    do i=1,n1
      r1(i,:) = off_plane(:, 1, off_plane_aux(2) +i-1)
      c1(i,:) = off_plane(:, 2, off_plane_aux(2) +i-1)
    end do

    do i=1,n2
      r2(i,:) = off_plane(:, 1, off_plane_aux(3) + i - 1)
      c2(i,:) = off_plane(:, 2, off_plane_aux(3) + i - 1)
    end do

    deallocate(tmp, dist, on_plane ,on_plane_dist, on_plane_aux)
    deallocate(off_plane, off_plane_dist, off_plane_aux, cnt)

    return
  end subroutine

end module mod_lattice
