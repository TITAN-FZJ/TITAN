!-------------------------------------------------------------------------------------!
! TITAN - TIme-dependent  Transport and Angular momentum properties of Nanostructures !
!-------------------------------------------------------------------------------------!
!
! MODULE: mod_BrillouinZone
!
!> @author
!> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> Description of the Brillouin Zone in reciprocal space
!> Contains types and routines to setup the Brillouin Zone
!
! REVISION HISTORY:
! 26 September 2017 - Last Changes
!-------------------------------------------------------------------------------------!
module mod_BrillouinZone
   use mod_f90_kind, only: double
   implicit none

   type :: BrillouinZone
      integer*8 :: nkpt   = 0
      integer*8 :: nkpt_x = 0
      integer*8 :: nkpt_y = 0
      integer*8 :: nkpt_z = 0

      real(double), dimension(:,:), allocatable :: kp
      real(double), dimension(:),   allocatable :: w
      real(double), dimension(3) :: b1, b2, b3

   contains
      procedure :: free => deallocate_BrillouinZone
      procedure :: isAlloc => isAlloc_BrillouinZone
      procedure :: count => count_BrillouinZone
      procedure :: print => output_kpoints
   end type BrillouinZone

   type, extends(BrillouinZone) :: FractionalBrillouinZone
     integer*8 :: workload = 0
     integer*8 :: first = 0
     integer*8 :: last = 0
     integer*4 :: size, rank, comm
   contains
     procedure :: setup_fraction => genFraction
     procedure :: generate_2d_fraction => gen2DFraction
     procedure :: generate_3d_fraction => gen3DFraction
   end type FractionalBrillouinZone

   type(FractionalBrillouinZone) :: realBZ

contains

  subroutine genFraction(self, sys, rank, size, comm)
    use mod_mpi_pars, only: calcWorkload
    use mod_System,   only: System
    implicit none
    class(FractionalBrillouinZone) :: self
    type(System),intent(in) :: sys
    integer*4,   intent(in) :: rank, size, comm

    if(allocated(self%kp)) deallocate(self%kp)
    if(allocated(self%w )) deallocate(self%w)

    self % rank = rank
    self % size = size
    self % comm = comm

    call self%count(sys)
    call calcWorkload(self%nkpt, self%size, self%rank, self%first, self%last)
    self%workload = self%last - self%first + 1
    if(sys%lbulk) then
      call self%generate_3d_fraction(sys,self%first,self%last)
    else
      call self%generate_2d_fraction(sys,self%first,self%last)
    end if

  end subroutine genFraction

  subroutine count_BrillouinZone(self,sys)
    use mod_System, only: System
    implicit none
    class(BrillouinZone)     :: self
    type(System), intent(in) :: sys
    integer*8                :: total

    if(sys%lbulk) then
      total = self%nkpt_x * self%nkpt_y * self%nkpt_z
      self%nkpt = count_3D_BZ(total,sys%a1,sys%a2,sys%a3)
    else
      total = self%nkpt_x * self%nkpt_y
      self%nkpt = count_2D_BZ(total,sys%a1,sys%a2)
    end if

  end subroutine count_BrillouinZone

  subroutine deallocate_BrillouinZone(self)
    implicit none
    class(BrillouinZone) :: self
    if(allocated(self%kp)) deallocate(self%kp)
    if(allocated(self%w)) deallocate(self%w)
  end subroutine deallocate_BrillouinZone

  logical function isAlloc_BrillouinZone(self)
     implicit none
     class(BrillouinZone) :: self
     if(allocated(self%kp)) then
        isAlloc_BrillouinZone = .true.
     else
        isAlloc_BrillouinZone = .false.
     end if
  end function isAlloc_BrillouinZone

  subroutine gen3DFraction(self,sys,first,last)
    use mod_f90_kind,  only: double
    use mod_constants, only: tpi
    use mod_tools,     only: cross
    use mod_mpi_pars,  only: abortProgram
    use mod_System,    only: System
    use mod_tools,     only: itos
    implicit none
    class(FractionalBrillouinZone) :: self
    type(System), intent(in)       :: sys
    integer*8,    intent(in)       :: first, last
    real(double), dimension(3,8)   :: bz_vec
    real(double), dimension(3)     :: diff
    real(double), dimension(3)     :: kp, b1, b2, b3
    real(double) :: vol
    real(double) :: smallest_dist, distance, ini_smallest_dist
    integer*8    :: nkpt, l
    integer*8    :: nx, ny, nz
    integer*8    :: count, added, weight, range
    integer      :: j

    allocate( self%w(self%workload), self%kp(3,self%workload) )
    self%w = 1.d0
    nkpt = self%nkpt_x * self%nkpt_y * self%nkpt_z
    vol  = tpi / dot_product(sys%a1, cross(sys%a2,sys%a3))
    b1   = vol * cross(sys%a2,sys%a3)
    b2   = vol * cross(sys%a3,sys%a1)
    b3   = vol * cross(sys%a1,sys%a2)

    bz_vec(1:3,1) = 0.d0
    bz_vec(1:3,2) = b1
    bz_vec(1:3,3) = b2
    bz_vec(1:3,4) = b1 + b2
    bz_vec(1:3,5) = b3
    bz_vec(1:3,6) = b1 + b3
    bz_vec(1:3,7) = b2 + b3
    bz_vec(1:3,8) = b1 + b2 + b3

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(b1 + b2 + b3, b1 + b2 + b3))
    count = 0
    added = 0
    !Run over all the kpoints generated initially.
    do l = 1, nkpt
      if(count > last) exit
      weight = 0
      range = 0
      nx = mod(floor(dble(l-1) / dble(self%nkpt_y * self%nkpt_z)), self%nkpt_x)
      ny = mod(floor(dble(l-1) / dble(self%nkpt_z)), self%nkpt_y)
      nz = mod(l-1, self%nkpt_z)
      kp = dble(nx)*b1 / dble(self%nkpt_x) + dble(ny)*b2 / dble(self%nkpt_y) + dble(nz)*b3 / dble(self%nkpt_z)

      smallest_dist = ini_smallest_dist
      !Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer.
      do j = 1, 8
        diff = kp - bz_vec(:,j)
        distance = sqrt(dot_product(diff, diff))
        if(distance < smallest_dist) smallest_dist = distance
      end do

      !Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 8
        diff = kp - bz_vec(:,j)
        distance = sqrt(dot_product(diff,diff))
        if( abs(distance-smallest_dist) < 1.d-12 ) then
          count = count + 1
          weight = weight + 1
          if(count >= first .and. count <= last ) then
            added = added + 1
            range = range + 1
            self%kp(:,added) = diff
          end if
        end if
      end do
      self%w(added-range+1:added) = 1.0 / dble(weight)
    end do
    self%w = self%w / dble(nkpt)
    if(added > self%workload) call abortProgram("[gen3DFraction] Generated more points than it should have! ")
    if(added < self%workload) call abortProgram("[gen3DFraction] Generated less points than it should have! ")
  end subroutine gen3DFraction

  integer*8 function count_3D_BZ(nkpt_in, a1, a2, a3)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    implicit none

    integer*8, intent(in) :: nkpt_in
    real(double), dimension(3), intent(in) :: a1, a2, a3
    real(double) :: vol
    real(double), dimension(3,8) :: bz_vec
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) ::  diff
    real(double), dimension(3) :: kp, b1, b2, b3
    integer*8 :: l, nkpt_x, nkpt_y, nkpt_z, nkpt, numextrakbz
    integer*8 :: nkpt_perdim !n. of k point per dimension
    integer*8 :: nx, ny, nz
    integer   :: j, smallest_index

    nkpt_perdim = ceiling((dble(nkpt_in))**(1.d0/3.d0))
    nkpt_x = nkpt_perdim
    nkpt_y = nkpt_perdim
    nkpt_z = nkpt_perdim

    nkpt = nkpt_x * nkpt_y * nkpt_z

    vol = tpi / dot_product(a1, cross(a2,a3))
    b1 = vol * cross(a2, a3)
    b2 = vol * cross(a3, a1)
    b3 = vol * cross(a1, a2)

    bz_vec(1:3,1) = 0.d0
    bz_vec(1:3,2) = b1
    bz_vec(1:3,3) = b2
    bz_vec(1:3,4) = b1 + b2
    bz_vec(1:3,5) = b3
    bz_vec(1:3,6) = b1 + b3
    bz_vec(1:3,7) = b2 + b3
    bz_vec(1:3,8) = b1 + b2 + b3

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(b1 + b2 + b3, b1 + b2 + b3))
    numextrakbz = 0
    !Run over all the kpoints generated initially.
    !$omp parallel do default(none) reduction(+:numextrakbz) if(nkpt > 1000000) &
    !$omp& private(l, j, nx, ny, nz, kp, smallest_dist, smallest_index, diff, distance) &
    !$omp& shared(nkpt, ini_smallest_dist, bz_vec, b1, b2, b3, nkpt_x, nkpt_y, nkpt_z)
    do l=1, nkpt
      nx = mod(floor(dble(l-1) / dble(nkpt_y * nkpt_z)), nkpt_x)
      ny = mod(floor(dble(l-1) / dble(nkpt_z)), nkpt_y)
      nz = mod(l-1, nkpt_z)
      kp = dble(nx)*b1 / dble(nkpt_x) + dble(ny)*b2 / dble(nkpt_y) + dble(nz)*b3 / dble(nkpt_z)

      smallest_dist=ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer.
      do j=1, 8
        diff = kp - bz_vec(:,j)
        distance = sqrt(dot_product(diff, diff))
        if(distance < smallest_dist) then
          smallest_dist = distance
          smallest_index = j
        end if
      end do
      ! Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 8
        diff=kp - bz_vec(:,j)
        distance=sqrt(dot_product(diff,diff))
        if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. j/=smallest_index ) numextrakbz=numextrakbz+1
      end do
    end do
    !$omp end parallel do
    count_3D_BZ = nkpt + numextrakbz

  end function count_3D_BZ

  subroutine gen2DFraction(self,sys,first,last)
    use mod_f90_kind,  only: double
    use mod_constants, only: tpi
    use mod_tools,     only: cross
    use mod_mpi_pars,  only: abortProgram
    use mod_System,    only: System
    implicit none
    class(FractionalBrillouinZone) :: self
    type(System), intent(in)       :: sys
    integer*8,    intent(in)       :: first, last
    real(double), dimension(3,4)   :: bz_vec
    real(double), dimension(3)     ::  diff, zdir
    real(double), dimension(3)     :: kp, b1, b2
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double) :: vol
    integer*8    :: l, nkpt
    integer*8    :: nx, ny
    integer*8    :: count, added, weight, range
    integer      :: j

    zdir = [0.0,0.0,1.0]
    allocate( self%w(self%workload), self%kp(3,self%workload) )
    self%w = 1.d0
    nkpt = self%nkpt_x * self%nkpt_y
    vol  = tpi/dot_product(zdir, cross(sys%a1,sys%a2))
    b1   = vol*cross(sys%a1,zdir)
    b2   = vol*cross(zdir,sys%a2)

    bz_vec(:,1) = 0.d0
    bz_vec(:,2) = b1
    bz_vec(:,3) = b2
    bz_vec(:,4) = b1 + b2

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(b1 + b2, b1 + b2))
    count = 0
    added = 0
    !Run over all the kpoints generated initially.
    do l=1, nkpt
      if(count > last) exit
      weight = 0
      range = 0
      nx = mod(floor(dble(l-1) / dble(self%nkpt_y)), self%nkpt_x)
      ny = mod(l-1, self%nkpt_y)
      kp = dble(nx)*b1 / dble(self%nkpt_x) + dble(ny)*b2 / dble(self%nkpt_y)

      smallest_dist = ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer.
      do j = 1, 4
        diff = kp - bz_vec(:,j)
        distance = sqrt(dot_product(diff, diff))
        if(distance < smallest_dist) smallest_dist = distance
      end do

      ! Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 4
        diff = kp - bz_vec(:,j)
        distance=sqrt(dot_product(diff,diff))
        if( abs(distance-smallest_dist) < 1.d-12 ) then
          count = count + 1
          weight = weight + 1
          if(count >= first .and. count <= last ) then
            added = added + 1
            range = range + 1
            self%kp(:,added) = diff
          end if
        end if
      end do
      self%w(added-range+1:added) = 1.0 / dble(weight)
    end do
    self%w = self%w / dble(nkpt)
    if(added > self%workload) call abortProgram("[gen2DFraction] Generated more points than it should have!")
    if(added < self%workload) call abortProgram("[gen2DFraction] Generated less points than it should have!")
  end subroutine gen2DFraction

  integer*8 function count_2D_BZ(nkpt_in, a1, a2)
    use mod_f90_kind,  only: double
    use mod_constants, only: tpi
    use mod_tools,     only: cross
    implicit none

    integer*8, intent(in) :: nkpt_in
    real(double), dimension(3), intent(in) :: a1,a2
    real(double), dimension(3) :: kp, b1, b2
    real(double), dimension(3) :: zdir
    real(double), dimension(3,4) :: bz_vec
    real(double), dimension(3)   :: diff
    real(double) :: smallest_dist, distance, ini_smallest_dist, vol
    integer   :: j, smallest_index
    integer*8 :: l, nkpt_x, nkpt_y, nx, ny, nkpt, nkpt_perdim, numextrakbz

    zdir = [0,0,1]
    nkpt_perdim = ceiling(sqrt(dble(nkpt_in)))
    nkpt_x = nkpt_perdim
    nkpt_y = nkpt_perdim

    nkpt = nkpt_x * nkpt_y

    vol = tpi / dot_product(zdir, cross(a1,a2))
    b1 = vol * cross(a1, zdir)
    b2 = vol * cross(zdir, a2)

    bz_vec(:,1) = 0.d0
    bz_vec(:,2) = b1
    bz_vec(:,3) = b2
    bz_vec(:,4) = b1 + b2

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(b1 + b2, b1 + b2))
    numextrakbz = 0
    !Run over all the kpoints generated initially.
    !$omp parallel do default(none) reduction(+:numextrakbz) if(nkpt > 1000000) &
    !$omp& private(l,j,nx,ny,kp,smallest_dist, smallest_index, diff, distance) &
    !$omp& shared(nkpt, ini_smallest_dist, bz_vec, b1, b2, nkpt_x, nkpt_y)
    do l = 1, nkpt

      nx = mod(floor(dble(l-1) / dble(nkpt_y)), nkpt_x)
      ny = mod(l-1, nkpt_y)
      kp = dble(nx)*b1 / dble(nkpt_x) + dble(ny)*b2 / dble(nkpt_y)

      smallest_dist = ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer.
      do j = 1, 4
        diff = kp - bz_vec(:,j)
        distance = sqrt(dot_product(diff, diff))
        if(distance < smallest_dist) then
          smallest_dist = distance
          smallest_index = j
        end if
      end do
      !Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j = 1, 4
        diff = kp - bz_vec(:, j)
        distance = sqrt(dot_product(diff, diff))
        if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. j /= smallest_index ) numextrakbz = numextrakbz + 1
      end do
    end do
    !$omp end parallel do
    count_2D_BZ = nkpt + numextrakbz
  end function count_2D_BZ

  subroutine output_kpoints(self)
    implicit none
    class(BrillouinZone) :: self
    integer :: i

    open (unit=3333, file='kpoints',status='replace')
    write(unit=3333,fmt="(a)") ' #      kx            ky            kz            wk'
    do i=1,self%nkpt
       write(unit=3333,fmt="(4(f12.9,2x))") self%kp(1:3,i),self%w(i)
    end do
    close(3333)

  end subroutine output_kpoints


! This subroutine generates all the k-points
! and store the different values of the
! component "component" (1,2,3)
  subroutine store_diff(nkpt_in, b1, b2, b3, component, ndiffk, diff_k_unsrt)
    use mod_f90_kind,  only: double
    use mod_constants, only: tpi
    use mod_tools,     only: itos
    use mod_mpi_pars,  only: abortProgram,rField
    implicit none

    integer*8, intent(in) :: nkpt_in
    !! Initial kpoints (as given by realBZ%nkpt_x * realBZ%nkpt_y * realBZ%nkpt_z)
    real(double), dimension(3), intent(in)  :: b1, b2, b3
    !! Reciprocal vectors
    integer,   intent(in)  :: component
    !! Which component to store the different values and the map
    integer*8, intent(out) :: ndiffk
    !! Number of different kp(component)
    real(double), dimension(:), allocatable, intent(out) :: diff_k_unsrt
    !! Different values of kp(component)
    real(double), dimension(:), allocatable              :: diff_k_temp
    !! Temporary array for different values of kp(component)
    real(double), dimension(3,8) :: bz_vec
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) :: diff
    real(double), dimension(3) :: kp
    integer*8  :: l, count, nkpt_x, nkpt_y, nkpt_z, nkpt
    integer*8  :: nkpt_perdim
    !! Number of k points per dimension
    integer*8  :: nx, ny, nz
    integer    :: j, smallest_index

    nkpt_perdim = ceiling((dble(nkpt_in))**(1.d0/3.d0))
    nkpt_x = nkpt_perdim
    nkpt_y = nkpt_perdim
    nkpt_z = nkpt_perdim

    nkpt = nkpt_x * nkpt_y * nkpt_z

    !! Allocating large temporary vector to store different values of kp(component)
    allocate( diff_k_temp(2*nkpt_perdim) )
    diff_k_temp(:) = 999.d0
    ndiffk = 0

    ! Counter of total number of points
    count = 0

    bz_vec(1:3,1) = 0.d0
    bz_vec(1:3,2) = b1
    bz_vec(1:3,3) = b2
    bz_vec(1:3,4) = b1 + b2
    bz_vec(1:3,5) = b3
    bz_vec(1:3,6) = b1 + b3
    bz_vec(1:3,7) = b2 + b3
    bz_vec(1:3,8) = b1 + b2 + b3

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(b1 + b2 + b3, b1 + b2 + b3))
    !Run over all the kpoints generated initially.
    !$omp parallel do default(none) reduction(+:count) if(nkpt > 1000000) &
    !$omp& private(l, j, nx, ny, nz, kp, smallest_dist, smallest_index, diff, distance) &
    !$omp& shared(nkpt, ini_smallest_dist, bz_vec, b1, b2, b3, nkpt_x, nkpt_y, nkpt_z, component, diff_k_temp, ndiffk)
    do l=1, nkpt
      nx = mod(floor(dble(l-1) / dble(nkpt_y * nkpt_z)), nkpt_x)
      ny = mod(floor(dble(l-1) / dble(nkpt_z)), nkpt_y)
      nz = mod(l-1, nkpt_z)
      kp = dble(nx)*b1 / dble(nkpt_x) + dble(ny)*b2 / dble(nkpt_y) + dble(nz)*b3 / dble(nkpt_z)

      smallest_dist=ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer.
      do j=1, 8
        diff = kp - bz_vec(:,j)
        distance = sqrt(dot_product(diff, diff))
        if(distance < smallest_dist) then
          smallest_dist = distance
          smallest_index = j
        end if
      end do
      !Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 8
        diff=kp - bz_vec(:,j)
        distance=sqrt(dot_product(diff,diff))
        if( ( abs(distance-smallest_dist) < 1.d-12 )  ) then

          count = count + 1 ! Counting the total number of points generated
          ! Checking if all numbers are different than current list (if any of them is zero, the point is already stored)
          !$omp critical
          if( all(abs( diff_k_temp(:) - abs(diff(component)) )> 1.d-12,1)  ) then
            ndiffk = ndiffk + 1
            diff_k_temp(ndiffk) = abs(diff(component))
          end if
          !$omp end critical

        end if

      end do
    end do
    !$omp end parallel do

    if(count/=realBZ%nkpt) &
      call abortProgram("[store_diff] Generated different number of points than it should have! count = " // trim(itos(count)) // ", nkpt = " // trim(itos(realBZ%nkpt)))

    allocate(diff_k_unsrt(ndiffk))
    diff_k_unsrt(:) = diff_k_temp(1:ndiffk)

  end subroutine store_diff

end module mod_BrillouinZone
