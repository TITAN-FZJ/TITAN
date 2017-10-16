module adaptiveMesh
   use mod_BrillouinZone
   integer :: total_points, total_points_real, total_points_imag
   integer, dimension(:,:), allocatable :: E_k_imag_mesh
   integer, dimension(:,:), allocatable :: E_k_real_mesh
   type(BrillouinZone), dimension(:), allocatable :: bzs

contains

   subroutine generateAdaptiveMeshes()
      implicit none
      call generate_imag_E_k_mesh()
      call generate_real_E_k_mesh()
      return
   end subroutine generateAdaptiveMeshes

   subroutine generate_imag_E_k_mesh()
      use mod_parameters, only: total_nkpt
      use EnergyIntegration, only: pn1, y
      use mod_System, only: s => sys
      use mod_mpi_pars
      use mod_BrillouinZone
      implicit none
      integer :: nkpt
      integer :: i,j,m

      allocate(bzs(pn1))
      total_points = 0
      do i = 1, pn1
         nkpt = get_nkpt(y(i), y(1), total_nkpt, s%lbulk)
         bzs(i) % nkpt = nkpt
         call bzs(i) % setup(s%lbulk, s%a1, s%a2, s%a3)
         total_points = total_points + bzs(i)%nkpt
      end do

      allocate(E_k_imag_mesh(total_points, 2))
      m = 0
      do i = 1, pn1
         do j = 1, bzs(i)%nkpt
            m = m + 1
            E_k_imag_mesh(m,1) = i
            E_k_imag_mesh(m,2) = j
         end do
      end do
      return

   end subroutine generate_imag_E_k_mesh

   subroutine generate_real_E_k_mesh()
      use EnergyIntegration, only: pn2
      implicit none
      integer :: i,j,m
      allocate(E_k_real_mesh(pn2*bzs(1)%nkpt, 2))
      total_points_real = pn2 * bzs(1)%nkpt
      m = 0
      do i = 1, pn2
         do j = 1, bzs(1)%nkpt
            m = m + 1
            E_k_real_mesh(m,1) = i
            E_k_real_mesh(m,2) = j
         end do
      end do
      return
   end subroutine generate_real_E_k_mesh

   integer function get_nkpt(e, e0, nkpt_total, bulk)
      use mod_f90_kind, only: double
      implicit none
      real(double), intent(in) :: e, e0
      logical, intent(in) :: bulk
      integer, intent(in) :: nkpt_total

      if(bulk) then
         get_nkpt = nkpt_total / (e/e0)**log(3.d0)
      else
         get_nkpt = nkpt_total / (e/e0)**log(2.d0)
      end if
      if(get_nkpt < 10) get_nkpt = 10
      return
   end function get_nkpt

end module adaptiveMesh
