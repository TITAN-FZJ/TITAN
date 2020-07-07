module adaptiveMesh
  use mod_BrillouinZone, only: FractionalBrillouinZone
  use MPI_f08,           only: MPI_Comm
  implicit none
  integer*8,                     dimension(:,:), allocatable :: E_k_imag_mesh
  type(FractionalBrillouinZone), dimension(:),   allocatable :: bzs
  integer*8,                     dimension(:),   allocatable :: all_nkpt,all_nkpt_rep
  integer*8 :: total_points, local_points
  integer*4 :: activeRank, activeSize
  type(MPI_Comm) :: activeComm
  integer   :: minimumBZmesh

  interface get_nkpt
    module procedure get_nkpt_int4, &
               get_nkpt_int8
  end interface get_nkpt

contains

  !! Get the number of points in the Brillouin Zone
  !! for all the energies in the imaginary axis
  subroutine generateAdaptiveMeshes(sys,pn1)
    use mod_f90_kind,      only: double
    use mod_parameters,    only: total_nkpt => kptotal_in
    use EnergyIntegration, only: y
    use mod_System,        only: System
    use mod_BrillouinZone, only: count_3D_BZ, count_2D_BZ, count_1D_BZ
    use mod_mpi_pars
    implicit none
    type(System) :: sys
    integer      :: i,pn1
    integer      :: nx, ny, nz
    integer*8    :: nall

    if(.not.allocated(all_nkpt)) allocate(all_nkpt(pn1),all_nkpt_rep(pn1))
    total_points = 0
    do i = 1, pn1
      nall = get_nkpt(y(i), y(1), total_nkpt, sys%isysdim)
      select case(sys%isysdim)
      case(3)
        nx = ceiling( (dble(nall))**(1.d0/3.d0), kind(nx) )
        ny = ceiling( (dble(nall))**(1.d0/3.d0), kind(ny) )
        nz = ceiling( (dble(nall))**(1.d0/3.d0), kind(nz) )
        nall = int( nx*ny*nz, kind(nall) )
        call count_3D_BZ(nall,sys%a1,sys%a2,sys%a3,all_nkpt(i),all_nkpt_rep(i))
      case(2)
        nx = ceiling( (dble(nall))**(1.d0/2.d0), kind(nx) )
        ny = ceiling( (dble(nall))**(1.d0/2.d0), kind(ny) )
        nz = 0
        nall = int( nx*ny, kind(nall) )
        call count_2D_BZ(nall,sys%a1,sys%a2,all_nkpt(i),all_nkpt_rep(i))
      case default
        nx = ceiling( (dble(nall)), kind(nx) )
        ny = 0
        nz = 0
        nall = int( nx, kind(nall) )
        call count_1D_BZ(nall,sys%a1,all_nkpt(i),all_nkpt_rep(i))
      end select
      total_points = total_points + all_nkpt(i)
    end do
  end subroutine generateAdaptiveMeshes

  !! Generate the distributed combined points {e,kx,ky,kz}
  !! (locally for a given MPI process)
  subroutine genLocalEKMesh(sys,rank,size,comm)
    use mod_parameters,    only: total_nkpt => kptotal_in
    use EnergyIntegration, only: pn1, y
    use mod_System,        only: System
    use mod_mpi_pars,      only: calcWorkload, MPI_Comm
    implicit none
    type(System),   intent(in) :: sys
    integer*4,      intent(in) :: rank
    integer*4,      intent(in) :: size
    type(MPI_Comm), intent(in) :: comm
    integer*8 :: firstPoint, lastPoint
    integer*8 :: j, m, n, p, q, nall
    integer   :: i

    activeComm = comm
    activeRank = rank
    activeSize = size
    call calcWorkload(total_points,activeSize,activeRank,firstPoint,lastPoint)
    local_points = lastPoint - firstPoint + 1
    allocate(E_k_imag_mesh(2,local_points))
    allocate(bzs(pn1))

    m = 0
    n = 0
    do i = 1, pn1
      if(m > lastPoint) exit
      if(m + all_nkpt(i) < firstPoint) then
        m = m + all_nkpt(i)
        cycle
      end if

      if(firstPoint <= m) then
        p = 1
      else
        p = firstPoint - m
      end if

      if(lastPoint < m + all_nkpt(i)) then
        q = lastPoint - m
      else
        q = all_nkpt(i)
      end if

      bzs(i) % rank = rank
      bzs(i) % size = size
      bzs(i) % comm = comm

      bzs(i) % workload = q - p + 1
      bzs(i) % nkpt = all_nkpt(i)
      nall = get_nkpt(y(i), y(1), total_nkpt, sys%isysdim)

      select case(sys%isysdim)
      case(3)
        bzs(i) % nkpt_x = ceiling( (dble(nall))**(1.d0/3.d0), kind(bzs(i) % nkpt_x) )
        bzs(i) % nkpt_y = ceiling( (dble(nall))**(1.d0/3.d0), kind(bzs(i) % nkpt_y) )
        bzs(i) % nkpt_z = ceiling( (dble(nall))**(1.d0/3.d0), kind(bzs(i) % nkpt_z) )
        call bzs(i) % gen3DFraction(sys,p,q)
      case(2)
        bzs(i) % nkpt_x = ceiling( (dble(nall))**(1.d0/2.d0), kind(bzs(i) % nkpt_x) )
        bzs(i) % nkpt_y = ceiling( (dble(nall))**(1.d0/2.d0), kind(bzs(i) % nkpt_y) )
        bzs(i) % nkpt_z = 1
        call bzs(i) % gen2DFraction(sys,p,q)
      case default
        bzs(i) % nkpt_x = ceiling( (dble(nall)), kind(bzs(i) % nkpt_x) )
        bzs(i) % nkpt_y = 1
        bzs(i) % nkpt_z = 1
        call bzs(i) % gen1DFraction(sys,p,q)
      end select

      do j = 1, bzs(i)%workload
        n = n + 1
        E_k_imag_mesh(1,n) = int(i,8)
        E_k_imag_mesh(2,n) = j
      end do
      m = m + all_nkpt(i)

    end do
  end subroutine genLocalEKMesh

  subroutine freeLocalEKMesh()
    implicit none
    integer*8 :: i
    do i = 1, local_points
      call bzs(E_k_imag_mesh(1,i)) % free()
    end do
    deallocate(bzs, E_k_imag_mesh)
  end subroutine freeLocalEKMesh

  !! Calculate the number of k-points for a given energy
  !! by using a power decay of the ratio of the imaginary parts
  !! to the dimension of the system (4-bit integer version)
  integer function get_nkpt_int4(e, e0, nkpt_total, sysdim)
    use mod_f90_kind, only: double
    implicit none
    real(double), intent(in) :: e, e0
    integer,      intent(in) :: sysdim
    integer,      intent(in) :: nkpt_total

    select case(sysdim)
    case(3)
      get_nkpt_int4 = nkpt_total / (e/e0)**sqrt(3.d0) !**log(3.d0)
    case(2)
      get_nkpt_int4 = nkpt_total / (e/e0)**sqrt(2.d0) !**log(2.d0)
    case default
      get_nkpt_int4 = nkpt_total / (e/e0)**sqrt(1.d0) !**log(1.d0)
    end select

    if(get_nkpt_int4 < minimumBZmesh ) get_nkpt_int4 = minimumBZmesh
  end function get_nkpt_int4

  !! Calculate the number of k-points for a given energy
  !! by using a power decay of the ratio of the imaginary parts
  !! to the dimension of the system (8-bit integer version)
  integer*8 function get_nkpt_int8(e, e0, nkpt_total, sysdim)
    use mod_f90_kind, only: double
    implicit none
    real(double), intent(in) :: e, e0
    integer,      intent(in) :: sysdim
    integer*8,    intent(in) :: nkpt_total

    select case(sysdim)
    case(3)
      get_nkpt_int8 = nkpt_total / (e/e0)**sqrt(3.d0) !**log(3.d0)
    case(2)
      get_nkpt_int8 = nkpt_total / (e/e0)**sqrt(2.d0) !**log(2.d0)
    case default
      get_nkpt_int8 = nkpt_total / (e/e0)**sqrt(1.d0) !**log(1.d0)
    end select

    if(get_nkpt_int8 < int(minimumBZmesh,8) ) get_nkpt_int8 = int(minimumBZmesh,8)
  end function get_nkpt_int8

end module adaptiveMesh
