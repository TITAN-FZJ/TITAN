module adaptiveMesh
  use mod_kind,          only: int32, int64
  use mod_BrillouinZone, only: FractionalBrillouinZone
  ! use MPI_f08,           only: MPI_Comm
  implicit none
  integer(int64),                dimension(:,:), allocatable :: E_k_imag_mesh
  type(FractionalBrillouinZone), dimension(:),   allocatable :: bzs
  integer(int64),                dimension(:),   allocatable :: all_nkpt,all_nkpt_norep
  integer(int64) :: total_points, local_points
  integer(int32) :: activeComm, activeRank, activeSize
  ! MPI_f08:
  ! integer(int32) :: activeRank, activeSize
  ! type(MPI_Comm) :: activeComm
  integer   :: minimumBZmesh

  interface get_nkpt
    module procedure get_nkpt_int4, &
               get_nkpt_int8
  end interface get_nkpt

contains

  subroutine generateAdaptiveMeshes(sys,pn1)
  !> Get the number of points in the Brillouin Zone
  !> for all the energies in the imaginary axis
    use mod_kind,          only: dp, int64
    use mod_parameters,    only: total_nkpt => kptotal_in
    use EnergyIntegration, only: y
    use mod_System,        only: System_type
    use mod_BrillouinZone, only: count_3D_BZ, count_2D_BZ, count_1D_BZ
    implicit none
    type(System_type) :: sys
    integer           :: i,pn1
    integer           :: nx, ny, nz
    integer(int64)    :: nall

    if(.not.allocated(all_nkpt)) allocate(all_nkpt(pn1),all_nkpt_norep(pn1))
    total_points = 0
    do i = 1, pn1
      nall = get_nkpt(y(i), y(1), total_nkpt, sys%isysdim)
      select case(sys%isysdim)
      case(3)
        nx = ceiling( (dble(nall))**(0.333333333333333_dp), kind(nx) )
        ny = ceiling( (dble(nall))**(0.333333333333333_dp), kind(ny) )
        nz = ceiling( (dble(nall))**(0.333333333333333_dp), kind(nz) )
        nall = int( nx*ny*nz, kind(nall) )
        call count_3D_BZ(nall,sys%a1,sys%a2,sys%a3,all_nkpt(i),all_nkpt_norep(i))
      case(2)
        nx = ceiling( (dble(nall))**(0.5_dp), kind(nx) )
        ny = ceiling( (dble(nall))**(0.5_dp), kind(ny) )
        nz = 0
        nall = int( nx*ny, kind(nall) )
        call count_2D_BZ(nall,sys%a1,sys%a2,all_nkpt(i),all_nkpt_norep(i))
      case default
        nx = ceiling( (dble(nall)), kind(nx) )
        ny = 0
        nz = 0
        nall = int( nx, kind(nall) )
        call count_1D_BZ(nall,sys%a1,all_nkpt(i),all_nkpt_norep(i))
      end select
      total_points = total_points + all_nkpt(i)
    end do
  end subroutine generateAdaptiveMeshes


  subroutine deallocateAdaptiveMeshes()
  !> Deallocate energy meshes
    implicit none

    if(allocated(all_nkpt))     deallocate( all_nkpt )
    if(allocated(all_nkpt_norep)) deallocate( all_nkpt_norep )
  end subroutine deallocateAdaptiveMeshes

  subroutine genLocalEKMesh(sys,rank,isize,comm)
  !> Generate the distributed combined points {e,kx,ky,kz}
  !> (locally for a given MPI process)
    use mod_kind,          only: dp, int32, int64
    use mod_parameters,    only: total_nkpt => kptotal_in
    use EnergyIntegration, only: pn1, y
    use mod_System,        only: System_type
    use mod_mpi_pars,      only: calcWorkload!, MPI_Comm
    implicit none

    type(System_type), intent(in) :: sys
    integer(int32),    intent(in) :: rank
    integer(int32),    intent(in) :: isize
    integer(int32),    intent(in) :: comm
    ! MPI_f08:
    ! type(System_type),   intent(in) :: sys
    ! integer(int32),      intent(in) :: rank
    ! integer(int32),      intent(in) :: isize
    ! type(MPI_Comm), intent(in) :: comm
    integer(int64) :: firstPoint, lastPoint
    integer(int64) :: j, m, n, p, q, nall
    integer        :: i

    activeComm = comm
    activeRank = rank
    activeSize = isize
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
      bzs(i) % isize = isize
      bzs(i) % comm = comm

      bzs(i) % workload = q - p + 1
      bzs(i) % nkpt = all_nkpt(i)
      nall = get_nkpt(y(i), y(1), total_nkpt, sys%isysdim)

      select case(sys%isysdim)
      case(3)
        bzs(i) % nkpt_x = ceiling( (dble(nall))**(0.333333333333333_dp), kind(bzs(i) % nkpt_x) )
        bzs(i) % nkpt_y = ceiling( (dble(nall))**(0.333333333333333_dp), kind(bzs(i) % nkpt_y) )
        bzs(i) % nkpt_z = ceiling( (dble(nall))**(0.333333333333333_dp), kind(bzs(i) % nkpt_z) )
        call bzs(i) % gen3DFraction(sys,p,q)
      case(2)
        bzs(i) % nkpt_x = ceiling( (dble(nall))**(0.5_dp), kind(bzs(i) % nkpt_x) )
        bzs(i) % nkpt_y = ceiling( (dble(nall))**(0.5_dp), kind(bzs(i) % nkpt_y) )
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
    use mod_kind, only: int64
    implicit none
    integer(int64) :: i

    do i = 1, local_points
      call bzs(E_k_imag_mesh(1,i)) % free()
    end do

    if( allocated(bzs) ) deallocate(bzs)
    if( allocated(E_k_imag_mesh) ) deallocate(E_k_imag_mesh)

  end subroutine freeLocalEKMesh

  integer function get_nkpt_int4(e, e0, nkpt_total, sysdim)
  !> Calculate the number of k-points for a given energy
  !> by using a power decay of the ratio of the imaginary parts
  !> to the dimension of the system (4-bit integer version)
    use mod_kind, only: dp
    implicit none
    real(dp), intent(in) :: e, e0
    integer,  intent(in) :: sysdim
    integer,  intent(in) :: nkpt_total

    select case(sysdim)
    case(3)
      get_nkpt_int4 = ceiling( nkpt_total/((e/e0)**sqrt(3._dp)) ) !**log(3._dp)
    case(2)
      get_nkpt_int4 = ceiling( nkpt_total/((e/e0)**sqrt(2._dp)) ) !**log(2._dp)
    case default
      get_nkpt_int4 = ceiling( nkpt_total/((e/e0)**sqrt(1._dp)) ) !**log(1._dp)
    end select

    if(get_nkpt_int4 < minimumBZmesh ) get_nkpt_int4 = minimumBZmesh
  end function get_nkpt_int4

  integer(int64) function get_nkpt_int8(e, e0, nkpt_total, sysdim)
  !> Calculate the number of k-points for a given energy
  !> by using a power decay of the ratio of the imaginary parts
  !> to the dimension of the system (8-bit integer version)
    use mod_kind, only: dp, int64
    implicit none
    real(dp),       intent(in) :: e, e0
    integer,        intent(in) :: sysdim
    integer(int64), intent(in) :: nkpt_total
    select case(sysdim)
    case(3)
      get_nkpt_int8 = ceiling( nkpt_total/((e/e0)**sqrt(3._dp)) ) !**log(3._dp)
    case(2)
      get_nkpt_int8 = ceiling( nkpt_total/((e/e0)**sqrt(2._dp)) ) !**log(2._dp)
    case default
      get_nkpt_int8 = ceiling( nkpt_total/((e/e0)**sqrt(1._dp)) ) !**log(1._dp)
    end select

    if(get_nkpt_int8 < int(minimumBZmesh,8) ) get_nkpt_int8 = int(minimumBZmesh,8)
  end function get_nkpt_int8

end module adaptiveMesh
