!-------------------------------------------------------------------------------
! juTITAN
!-------------------------------------------------------------------------------
!
! MODULE: mod_lattice
!
!> @author
!> Filipe Guimaraes, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> Calculation of next nearest neighbours in plane and out of plane.
!
! REVISION HISTORY:
! 12 August 2015 - Initial Version
! 14 March 2017 - Last revision
!-------------------------------------------------------------------------------
module mod_lattice
  use mod_f90_kind
  implicit none

  integer :: n01      !< Number of in-plane 1st nearest neighbors
  integer :: n02      !< Number of in-plane 2nd nearest neighbors
  integer :: n0       !< Number of all in-plane nearest neighbours
  integer :: n1       !< Number of inter-plane 1st nearest neighbors
  integer :: n2       !< Number of inter-plane 2nd nearest neighbors
  integer :: plnn     !< Number of nearest neighbors planes - to mount (tri-)diagonal matrices
  real(double), dimension(:,:), allocatable :: r0 !< Position of in-plane 1st and 2nd nearest neighbors
  real(double), dimension(:,:), allocatable :: c0 !< Direction cosines of in-plane 1st and 2nd nearest neighbors
  real(double), dimension(:,:), allocatable :: r1 !< Position of inter-plane 1st nearest neighbors
  real(double), dimension(:,:), allocatable :: c1 !< Direction cosines of inter-plane 1st nearest neighbors
  real(double), dimension(:,:), allocatable :: r2 !< Position of inter-plane 1st nearest neighbors
  real(double), dimension(:,:), allocatable :: c2 !< Direction cosines of inter-plane 1st nearest neighbors

  real(double), dimension(3) :: versor_oop   !< Out-of-plane unit vector
  integer :: transverse_neighbors            !< TBD.
  integer :: longitudinal_neighbors          !< TBD.
  real(double), dimension(3) :: versor_Eperp !< In-plane unit vector perpendicular to current direction

contains
  !-----------------------------------------------------------------------------
  !> @author
  !> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut
  !
  ! DESCRIPTION:
  !> Calculation of next nearest neighbours in plane and out of plane.
  !
  ! REVISION HISTORY:
  ! 14 March 2017 - Current Revision
  !> @param[in] a1 Lattice unit vector
  !> @param[in] a2 Lattice unit vector
  !> @param[in] a3 Lattice unit vector
  !> @param[in] p  Plane direction vector
  !-----------------------------------------------------------------------------
  subroutine next_neighbour_init(a1,a2,a3,p)
    use mod_f90_kind
    use mod_constants
    implicit none
    real(double), intent(in) :: a1(3)
    real(double), intent(in) :: a2(3)
    real(double), intent(in) :: a3(3)
    real(double), intent(in) :: p(3)
    real(double), allocatable :: tmp(:,:,:)          ! Array contains all neighbours in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real(double), allocatable :: dist(:)             ! Contains all distances in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real(double), allocatable :: on_plane(:,:,:)     ! Array of type (3,2,on_cnt), all coordinates and directional cosines of neighbours in plane described by vector
    real(double), allocatable :: off_plane(:,:,:)    ! Array of type (3,2,off_cnt), containing half the intra-plane elements
    real(double), allocatable :: on_plane_dist(:)    ! Contains the distance of each point to the Origin
    real(double), allocatable :: off_plane_dist(:,:) ! Contains the distance of each point to the Origin and to the plane
    real(double) :: dist_tmp
    real(double) :: pos_tmp(3)
    real(double) :: cos_tmp(3)
    real(double), allocatable :: cnt(:,:)            ! Contains distance and stage informations for all elements in tmp
    integer, allocatable :: on_plane_aux(:)          ! on_plane_aux(i) contains the first element of the i-1 -th neighbour stage
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
