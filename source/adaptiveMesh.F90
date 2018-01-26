module adaptiveMesh
   use mod_BrillouinZone, only: FractionalBrillouinZone
   integer*8 :: total_points, local_points
   integer, dimension(:,:), allocatable :: E_k_imag_mesh
   type(FractionalBrillouinZone), dimension(:), allocatable :: bzs
   integer, dimension(:),allocatable :: all_nkpt
   integer :: activeComm, activeRank, activeSize

contains

   subroutine generateAdaptiveMeshes(pn1)
      use mod_parameters, only: total_nkpt => kptotal_in
      use EnergyIntegration, only: y
      use mod_System, only: s => sys
      use mod_BrillouinZone, only: count_3D_BZ, count_2D_BZ
      use mod_mpi_pars
      implicit none
      integer :: i,pn1

      allocate(all_nkpt(pn1))
      total_points = 0
      do i = 1, pn1
         if(s%lbulk) then
            all_nkpt(i) = count_3D_BZ(get_nkpt(y(i), y(1), total_nkpt, s%lbulk),s%a1,s%a2,s%a3)
         else
            all_nkpt(i) = count_2D_BZ(get_nkpt(y(i), y(1), total_nkpt, s%lbulk),s%a1,s%a2)
         end if
         total_points = total_points + all_nkpt(i)
      end do
      return
   end subroutine generateAdaptiveMeshes

   subroutine genLocalEKMesh(rank, size, comm)
      use mod_parameters, only: total_nkpt => kptotal_in
      use EnergyIntegration, only: pn1, y
      use mod_System, only: s => sys
      use mod_mpi_pars, only: calcWorkload
      implicit none
      integer, intent(in) :: rank
      integer, intent(in) :: size
      integer, intent(in) :: comm
      integer*8 :: firstPoint, lastPoint
      integer*8 :: i, j, m, n, p, q
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
         if(m + all_nkpt(i) < firstPoint) then
            m = m + all_nkpt(i)
            cycle
         end if

         p = firstPoint - m
         if(lastPoint < m + all_nkpt(i)) then
            q = lastPoint - m
         else
            q = all_nkpt(i)
         end if

         bzs(i) % a1 = s % a1
         bzs(i) % a2 = s % a2
         bzs(i) % a3 = s % a3
         bzs(i) % bulk = s%lbulk

         if(s%lbulk) then
            bzs(i) % nkpt_x = ceiling((dble(nkpt))**(1.d0/3.d0))
            bzs(i) % nkpt_y = ceiling((dble(nkpt))**(1.d0/3.d0))
            bzs(i) % nkpt_z = ceiling((dble(nkpt))**(1.d0/3.d0))
            call bzs(i) % generate_3d_fraction(p,q,all_nkpt(i))
         else
            bzs(i) % nkpt_x = ceiling((dble(nkpt))**(1.d0/2.d0))
            bzs(i) % nkpt_y = ceiling((dble(nkpt))**(1.d0/2.d0))
            bzs(i) % nkpt_z = 0
            call bzs(i) % generate_3d_fraction(p,q,all_nkpt(i))
         end if

         do j = 1, bzs(i)%workload
            n = n + 1
            E_k_imag_mesh(1,n) = i
            E_k_imag_mesh(2,n) = j
         end do
         m = m + q

         ! do j = 1, all_nkpt(i)
         !    m = m + 1
         !    if(m < firstPoint .or. m > lastPoint) cycle
         !    n = n + 1
         !    E_k_imag_mesh(1,n) = i
         !    E_k_imag_mesh(2,n) = j
         !
         !    if(.not. bzs(i)%isAlloc()) then
         !       nkpt = get_nkpt(y(i), y(1), total_nkpt, s%lbulk)
         !       bzs(i) % a1 = s % a1
         !       bzs(i) % a2 = s % a2
         !       bzs(i) % a3 = s % a3
         !       bzs(i) % bulk = s%lbulk
         !
         !       if(s%lbulk) then
         !          bzs(i) % nkpt_x = ceiling((dble(nkpt))**(1.d0/3.d0))
         !          bzs(i) % nkpt_y = ceiling((dble(nkpt))**(1.d0/3.d0))
         !          bzs(i) % nkpt_z = ceiling((dble(nkpt))**(1.d0/3.d0))
         !       else
         !          bzs(i) % nkpt_x = ceiling((dble(nkpt))**(1.d0/2.d0))
         !          bzs(i) % nkpt_y = ceiling((dble(nkpt))**(1.d0/2.d0))
         !          bzs(i) % nkpt_z = 0
         !       end if
         !
         !       call bzs(i) % setup()
         !    end if
         ! end do
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
