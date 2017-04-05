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
  use mod_f90_kind, only: double
  implicit none

  integer :: n01  !< Number of in-plane 1st nearest neighbors
  integer :: n02  !< Number of in-plane 2nd nearest neighbors
  integer :: n0   !< Number of all in-plane nearest neighbours
  integer :: n1   !< Number of inter-plane 1st nearest neighbors
  integer :: n2   !< Number of inter-plane 2nd nearest neighbors
  integer :: plnn !< Number of nearest neighbors planes - to mount (tri-)diagonal matrices
  real(double), dimension(:,:), allocatable :: r0 !< Position of in-plane 1st and 2nd nearest neighbors
  real(double), dimension(:,:), allocatable :: c0 !< Direction cosines of in-plane 1st and 2nd nearest neighbors
  real(double), dimension(:,:), allocatable :: r1 !< Position of inter-plane 1st nearest neighbors
  real(double), dimension(:,:), allocatable :: c1 !< Direction cosines of inter-plane 1st nearest neighbors
  real(double), dimension(:,:), allocatable :: r2 !< Position of inter-plane 1st nearest neighbors
  real(double), dimension(:,:), allocatable :: c2 !< Direction cosines of inter-plane 1st nearest neighbors

  real(double), dimension(3) :: versor_oop   !< Out-of-plane unit vector
  real(double), dimension(3) :: versor_Eperp !< In-plane unit vector perpendicular to current direction
  integer :: transverse_neighbors            !< Description missing.
  integer :: longitudinal_neighbors          !< Description missing.

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
  subroutine next_neighbour_init()
    use mod_f90_kind, only: double
    use mod_mpi_pars
    use mod_tools
    use mod_parameters, only: a1, a2, a3, &
                              a1_pln, a2_pln, pln_dir, &
                              nn_stages, outputunit
    implicit none
    real(double), dimension(:,:,:), allocatable :: nn                      ! Array contains all neighbours in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real(double), dimension(:),     allocatable :: dist                    ! Contains all distances in the ellipsoid 2*nn_stages*(a1,a2,a3)
    real(double), dimension(:,:,:), allocatable :: on_plane                ! Array of type (3,2,on_cnt), all coordinates and directional cosines of neighbours in plane described by vector
    real(double), dimension(:,:,:), allocatable :: off_plane               ! Array of type (3,2,off_cnt), containing half the intra-plane elements
    real(double), dimension(:),     allocatable :: on_plane_dist           ! Contains the distance of each point to the Origin
    real(double), dimension(:,:),   allocatable :: off_plane_dist          ! Contains the distance of each point to the Origin and to the plane
    real(double), dimension(:,:),   allocatable :: cnt                     ! Contains distance and stage informations for all elements in tmp
    integer,      dimension(:),     allocatable :: on_plane_aux            ! on_plane_aux(i) contains the first element of the i-1 -th neighbour stage
    integer,      dimension(:),     allocatable :: off_plane_aux
    real(double), dimension(3)                  :: pos_tmp, cos_tmp
    real(double)                                :: dist_tmp
    integer                                     :: nn_size, nnt
    integer                                     :: on_cnt, off_cnt
    integer                                     :: i, j, l, m, n
    logical                                     :: pl_flag

    nnt = 2 * nn_stages

    allocate(nn(3, 2, (2*nnt+1)**3))
    allocate(dist((2*nnt+1)**3), cnt(2, nn_stages+2))
    allocate(on_plane_aux(nn_stages+2), off_plane_aux(nn_stages+2))

    ! Generate neighbours up to 2*nn_stages the number of wanted neighbours
    ! Insertion sort approach
    i = 1
    do l = -nnt, nnt
       do m = -nnt, nnt
          do n = -nnt, nnt
             pos_tmp  = l*a1 + m*a2 + n*a3
             dist_tmp = sqrt(dot_product(pos_tmp,pos_tmp))
             cos_tmp  = pos_tmp / dist_tmp
             j = i-1
             do while( j >= 1 .and. dist(j) > dist_tmp)
                nn(:,:, j+1) = nn(:,:,j)
                dist(j+1) = dist(j)
                j = j-1
             end do
             nn(:, 1, j+1) = pos_tmp
             nn(:, 2, j+1) = cos_tmp
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
    allocate( on_plane(3,2,nn_size),  on_plane_dist(   nn_size))
    allocate(off_plane(3,2,nn_size), off_plane_dist(2, nn_size))

    pln_dir = pln_dir / sqrt(dot_product(pln_dir, pln_dir))

    on_cnt = 0
    off_cnt = 0
    on_plane_aux = 1
    off_plane_aux = 1

    do i=1, nn_size
       dist_tmp = dot_product(pln_dir, nn(:,1,i))
       if(abs(dist_tmp) < 1d-9) then
          on_cnt = on_cnt + 1
          on_plane(:,:, on_cnt) = nn(:,:, i)
          on_plane_dist(on_cnt) = dist(i)
          do j=1, nn_stages+1
             if( dist(i) <= cnt(2,j) ) then
                on_plane_aux(j+1) = on_plane_aux(j+1) + 1
             end if
          end do
       else if(dist_tmp < 0) then
          off_cnt = off_cnt + 1
          off_plane(:,:, off_cnt)   = nn(:,:, i)
          off_plane_dist(1,off_cnt) = dist(i)
          off_plane_dist(2,off_cnt) = abs(dot_product(pln_dir, nn(:,1,i)))
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

    ! Definition of in-plane basis
    a1_pln = r0(1,:)
    a2_pln = 0.d0

    do i = 2, n0
       if( .not. is_parallel(r0(1,:), r0(i,:))) then
          a2_pln = r0(i,:)
          exit
       end if
    end do
    if(0.d0 == dot_product(a2_pln, a2_pln)) then
       if(myrank.eq.0) write(outputunit,"('[next_neighbour_init] No non-collinear in-plane vectors found')")
       call MPI_Finalize(ierr)
       stop
    end if
    ! end in-plane basis

    deallocate(nn, dist, on_plane ,on_plane_dist, on_plane_aux)
    deallocate(off_plane, off_plane_dist, off_plane_aux, cnt)

    return
  end subroutine next_neighbour_init
end module mod_lattice
