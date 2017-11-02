module adaptiveMesh
   use mod_BrillouinZone
   integer :: total_points, local_points
   integer, dimension(:,:), allocatable :: E_k_imag_mesh
   type(BrillouinZone), dimension(:), allocatable :: bzs
   integer, dimension(:),allocatable :: all_nkpt
   integer :: activeComm, activeRank, activeSize

contains

   subroutine generateAdaptiveMeshes()
      implicit none

      call calcTotalPoints()
      return
   end subroutine generateAdaptiveMeshes

   subroutine calcTotalPoints()
      use mod_parameters, only: total_nkpt
      use EnergyIntegration, only: pn1, y
      use mod_System, only: s => sys
      use mod_BrillouinZone
      use mod_mpi_pars
      implicit none
      integer :: i
      type(BrillouinZone) :: bzone

      allocate(all_nkpt(pn1))
      if(myrank == 0) then
         total_points = 0
         do i = 1, pn1
            all_nkpt(i) = count_3D_BZ(get_nkpt(y(i), y(1), total_nkpt, s%lbulk),s%a1,s%a2,s%a3) !bzone%nkpt
            total_points = total_points + all_nkpt(i)
         end do
      end if
      call MPI_Bcast(total_points, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(all_nkpt, pn1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      return
   end subroutine calcTotalPoints

   subroutine genLocalEKMesh(rank, size, comm)
      use mod_parameters, only: total_nkpt
      use EnergyIntegration, only: pn1, y
      use mod_System, only: s => sys
      use mod_BrillouinZone
      use mod_mpi_pars, only: calcWorkload
      implicit none
      integer, intent(in) :: rank
      integer, intent(in) :: size
      integer, intent(in) :: comm
      integer :: firstPoint, lastPoint
      integer :: i, j, m, n
      integer :: nkpt

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
         do j = 1, all_nkpt(i)
            m = m + 1
            if(m < firstPoint .or. m > lastPoint) cycle
            n = n + 1
            E_k_imag_mesh(1,n) = i
            E_k_imag_mesh(2,n) = j

            if(.not. bzs(i)%isAlloc()) then
               nkpt = get_nkpt(y(i), y(1), total_nkpt, s%lbulk)
               bzs(i) % nkpt = nkpt
               bzs(i) % nkpt_x = 0
               bzs(i) % nkpt_y = 0
               bzs(i) % nkpt_z = 0

               call bzs(i) % setup(s%lbulk, s%a1, s%a2, s%a3)
            end if
         end do
      end do

      return
   end subroutine genLocalEKMesh

   subroutine freeLocalEKMesh()
      implicit none
      integer :: i
      do i = 1, local_points
         call bzs(E_k_imag_mesh(1,i)) % free()
      end do
      deallocate(bzs, E_k_imag_mesh)
      return
   end subroutine freeLocalEKMesh

   ! subroutine generate_real_E_k_mesh()
   !    use EnergyIntegration, only: pn2
   !    implicit none
   !    integer :: i,j,m
   !    allocate(E_k_real_mesh(pn2*bzs(1)%nkpt, 2))
   !    total_points_real = pn2 * bzs(1)%nkpt
   !    m = 0
   !    do i = 1, pn2
   !       do j = 1, bzs(1)%nkpt
   !          m = m + 1
   !          E_k_real_mesh(m,1) = i
   !          E_k_real_mesh(m,2) = j
   !       end do
   !    end do
   !    return
   ! end subroutine generate_real_E_k_mesh

   integer function get_nkpt(e, e0, nkpt_total, bulk)
      use mod_f90_kind, only: double
      implicit none
      real(double), intent(in) :: e, e0
      logical, intent(in) :: bulk
      integer, intent(in) :: nkpt_total

      if(bulk) then
         get_nkpt = nkpt_total / (e/e0)**sqrt(3.d0) !**log(3.d0)
      else
         get_nkpt = nkpt_total / (e/e0)**sqrt(2.d0) !**log(2.d0)
      end if
      if(get_nkpt < 1000) get_nkpt = 1000
      return
   end function get_nkpt

end module adaptiveMesh
