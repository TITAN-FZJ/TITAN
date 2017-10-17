module adaptiveMesh
   use mod_BrillouinZone
   integer :: total_points, local_points
   integer, dimension(:,:), allocatable :: E_k_imag_mesh
   type(BrillouinZone), dimension(:), allocatable :: bzs
   integer, dimension(:),allocatable :: all_nkpt

contains

   subroutine generateAdaptiveMeshes()
      use EnergyIntegration, only: pn1
      implicit none

      call calcTotalPoints()
      call genImagEKMesh()

      return
   end subroutine generateAdaptiveMeshes

   subroutine calcTotalPoints()
      use mod_parameters, only: total_nkpt
      use EnergyIntegration, only: pn1, y
      use mod_System, only: s => sys
      use mod_BrillouinZone
      use mod_mpi_pars
      implicit none
      integer :: i,j,m,n
      type(BrillouinZone) :: bzone

      integer :: work, remainder, start_imag, end_imag
      allocate(all_nkpt(pn1))
      if(myrank == 0) then
         total_points = 0
         do i = 1, pn1
            bzone % nkpt = get_nkpt(y(i), y(1), total_nkpt, s%lbulk)
            bzone % nkpt_x = 0
            bzone % nkpt_y = 0
            bzone % nkpt_z = 0

            call bzone % setup(s%lbulk, s%a1, s%a2, s%a3)
            all_nkpt(i) = bzone%nkpt
            total_points = total_points + all_nkpt(i)
            call bzone % free()
         end do

      else
         call MPI_Bcast(total_points, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(all_nkpt, pn1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
      end if

      ! Calculate workload for each MPI process
      remainder = mod(total_points, numprocs_row)
      if(myrank_row < remainder) then
        work = ceiling(dble(total_points) / dble(numprocs_row))
        start_imag = myrank_row*work + 1
        end_imag = (myrank_row+1) * work
      else
        work = floor(dble(total_points) / dble(numprocs_row))
        start_imag = myrank_row*work + 1 + remainder
        end_imag = (myrank_row+1) * work + remainder
      end if

      local_points = end_imag - start_imag + 1
      allocate(E_k_imag_mesh(2,end_imag-start_imag+1))

      m = 0
      n = 0
      do i = 1, pn1
         if(m > end_imag) exit
         do j = 1, all_nkpt(i)
            m = m + 1
            if(m < start_imag .or. m > end_imag) cycle
            n = n + 1
            E_k_imag_mesh(1,n) = i
            E_k_imag_mesh(2,n) = j
         end do
      end do
      return

   end subroutine calcTotalPoints


   subroutine genImagEKMesh()
      use mod_parameters, only: total_nkpt
      use EnergyIntegration, only: pn1, y
      use mod_System, only: s => sys
      use mod_mpi_pars
      use mod_BrillouinZone
      implicit none
      integer :: nkpt
      integer :: i,j,m
      integer :: work, remainder, start_imag, end_imag

      allocate(bzs(pn1))
      total_points = 0

      do i = 1, local_points
         j = E_k_imag_mesh(1,i)
         m = E_k_imag_mesh(2,i)
         if(bzs(j)%isAlloc()) cycle
         nkpt = get_nkpt(y(j), y(1), total_nkpt, s%lbulk)
         bzs(j) % nkpt = nkpt
         bzs(j) % nkpt_x = 0
         bzs(j) % nkpt_y = 0
         bzs(j) % nkpt_z = 0

         call bzs(j) % setup(s%lbulk, s%a1, s%a2, s%a3)
         if(all_nkpt(j) /= bzs(j)%nkpt) print *, "Something wrong.", i,j,m,bzs(j)%nkpt, all_nkpt(j)

      end do

      return

   end subroutine genImagEKMesh

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
