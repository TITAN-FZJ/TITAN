module adaptiveMesh
   use mod_BrillouinZone, only: FractionalBrillouinZone
   integer,                       dimension(:,:), allocatable :: E_k_imag_mesh
   type(FractionalBrillouinZone), dimension(:),   allocatable :: bzs
   integer,                       dimension(:),   allocatable :: all_nkpt
   integer*8 :: total_points, local_points
   integer*4 :: activeComm, activeRank, activeSize
   integer   :: minimumBZmesh

  interface get_nkpt
     module procedure get_nkpt_int4, &
                      get_nkpt_int8
  end interface get_nkpt

contains

   subroutine generateAdaptiveMeshes(sys,pn1)
      use mod_f90_kind,      only: double
      use mod_parameters,    only: total_nkpt => kptotal_in
      use EnergyIntegration, only: y
      use mod_System,        only: System
      use mod_BrillouinZone, only: count_3D_BZ, count_2D_BZ
      use mod_mpi_pars
      implicit none
      type(System) :: sys
      integer      :: i,pn1
      integer*8    :: nx, ny, nz, nall

      if(.not.allocated(all_nkpt)) allocate(all_nkpt(pn1))
      total_points = 0
      do i = 1, pn1
         nall = get_nkpt(y(i), y(1), total_nkpt, sys%lbulk)
         if(sys%lbulk) then
            nx = ceiling((dble(nall))**(1.d0/3.d0))
            ny = ceiling((dble(nall))**(1.d0/3.d0))
            nz = ceiling((dble(nall))**(1.d0/3.d0))
            all_nkpt(i) = count_3D_BZ(nx*ny*nz,sys%a1,sys%a2,sys%a3)
         else
            nx = ceiling((dble(nall))**(1.d0/2.d0))
            ny = ceiling((dble(nall))**(1.d0/2.d0))
            nz = 0
            all_nkpt(i) = count_2D_BZ(nx*ny,sys%a1,sys%a2)
         end if
         total_points = total_points + all_nkpt(i)
      end do
   end subroutine generateAdaptiveMeshes

   subroutine genLocalEKMesh(sys,rank,size,comm)
      use mod_parameters,    only: total_nkpt => kptotal_in
      use EnergyIntegration, only: pn1, y
      use mod_System,        only: System
      use mod_mpi_pars,      only: calcWorkload
      implicit none
      type(System), intent(in) :: sys
      integer*4,    intent(in) :: rank
      integer*4,    intent(in) :: size
      integer*4,    intent(in) :: comm
      integer*8 :: firstPoint, lastPoint
      integer*8 :: i, j, m, n, p, q
      integer   :: nall

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

         if(firstPoint < m) then
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
         nall = get_nkpt(y(i), y(1), total_nkpt, sys%lbulk)

         if(sys%lbulk) then
            bzs(i) % nkpt_x = ceiling((dble(nall))**(1.d0/3.d0))
            bzs(i) % nkpt_y = ceiling((dble(nall))**(1.d0/3.d0))
            bzs(i) % nkpt_z = ceiling((dble(nall))**(1.d0/3.d0))
            call bzs(i) % generate_3d_fraction(sys,int(p,8),int(q,8))
         else
            bzs(i) % nkpt_x = ceiling((dble(nall))**(1.d0/2.d0))
            bzs(i) % nkpt_y = ceiling((dble(nall))**(1.d0/2.d0))
            bzs(i) % nkpt_z = 0
            call bzs(i) % generate_2d_fraction(sys,int(p,8),int(q,8))
         end if

         do j = 1, bzs(i)%workload
            n = n + 1
            E_k_imag_mesh(1,n) = i
            E_k_imag_mesh(2,n) = j
         end do
         m = m + all_nkpt(i)

         ! do j = 1, all_nkpt(i)
         !    m = m + 1
         !    if(m < firstPoint .or. m > lastPoint) cycle
         !    n = n + 1
         !    E_k_imag_mesh(1,n) = i
         !    E_k_imag_mesh(2,n) = j
         !
         !    if(.not. bzs(i)%isAlloc()) then
         !       nkpt = get_nkpt(y(i), y(1), total_nkpt, sys%lbulk)
         !
         !       if(sys%lbulk) then
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
   end subroutine genLocalEKMesh

   subroutine freeLocalEKMesh()
      implicit none
      integer :: i
      do i = 1, local_points
         call bzs(E_k_imag_mesh(1,i)) % free()
      end do
      deallocate(bzs, E_k_imag_mesh)
   end subroutine freeLocalEKMesh

   integer function get_nkpt_int4(e, e0, nkpt_total, bulk)
      use mod_f90_kind, only: double
      implicit none
      real(double), intent(in) :: e, e0
      logical, intent(in) :: bulk
      integer, intent(in) :: nkpt_total

      if(bulk) then
         get_nkpt_int4 = nkpt_total / (e/e0)**sqrt(3.d0) !**log(3.d0)
      else
         get_nkpt_int4 = nkpt_total / (e/e0)**sqrt(2.d0) !**log(2.d0)
      end if
      if(get_nkpt_int4 < minimumBZmesh ) get_nkpt_int4 = minimumBZmesh
   end function get_nkpt_int4

   integer*8 function get_nkpt_int8(e, e0, nkpt_total, bulk)
      use mod_f90_kind, only: double
      implicit none
      real(double), intent(in) :: e, e0
      logical, intent(in) :: bulk
      integer*8, intent(in) :: nkpt_total

      if(bulk) then
         get_nkpt_int8 = nkpt_total / (e/e0)**sqrt(3.d0) !**log(3.d0)
      else
         get_nkpt_int8 = nkpt_total / (e/e0)**sqrt(2.d0) !**log(2.d0)
      end if
      if(get_nkpt_int8 < minimumBZmesh ) get_nkpt_int8 = minimumBZmesh
   end function get_nkpt_int8

end module adaptiveMesh
