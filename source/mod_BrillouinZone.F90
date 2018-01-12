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
      integer*8 :: nkpt = 0
      integer :: nkpt_x = 0
      integer :: nkpt_y = 0
      integer :: nkpt_z = 0

      real(double), dimension(:,:), allocatable :: kp
      real(double), dimension(:), allocatable :: w
      real(double), dimension(3) :: a1, a2, a3
      real(double), dimension(3) :: b1, b2, b3
      logical :: bulk
   contains
      procedure :: setup => setup_BrillouinZone
      procedure :: free => deallocate_BrillouinZone
      procedure :: isAlloc => isAlloc_BrillouinZone
      procedure :: count => count_BrillouinZone
      procedure, private :: generate_2d => generate_2D_BZ
      procedure, private :: generate_3d => generate_3D_BZ
      procedure :: print => output_kpoints
   end type BrillouinZone

   type, extends(BrillouinZone) :: FractionalBrillouinZone
     integer*8 :: workload = 0
     integer*8 :: first = 0
     integer*8 :: last = 0
     integer :: size, rank, comm
   contains
     procedure :: setup_fraction => genFraction
     procedure, private :: generate_2d_fraction => gen2DFraction
     procedure, private :: generate_3d_fraction => gen3DFraction
   end type FractionalBrillouinZone

   type(FractionalBrillouinZone) :: realBZ

contains
  subroutine setup_BrillouinZone(self)
    implicit none
    class(BrillouinZone) :: self

    if(allocated(self%kp)) deallocate(self%kp)
    if(allocated(self%w)) deallocate(self%w)

    if(self%bulk) then
      call self % generate_3d()
    else
      call self % generate_2d()
    end if
    return
  end subroutine setup_BrillouinZone

  subroutine genFraction(self, rank, size, comm)
    use mod_mpi_pars, only: calcWorkload
    implicit none
    class(FractionalBrillouinZone) :: self
    integer, intent(in) :: rank, size, comm

    if(allocated(self%kp)) deallocate(self%kp)
    if(allocated(self%w)) deallocate(self%w)

    self % rank = rank
    self % size = size
    self % comm = comm

    call self%count()
    call calcWorkload(self%nkpt, self%size, self%rank, self%first, self%last)
    self%workload = self%last - self%first + 1
    if(self%bulk) then
      call self%generate_3d_fraction(self%first,self%last, self%nkpt)
    else
      call self%generate_2d_fraction(self%first,self%last, self%nkpt)
    end if

    return
  end subroutine genFraction

  subroutine count_BrillouinZone(self)
    implicit none
    class(BrillouinZone) :: self
    integer :: total

    if(self % bulk) then
      total = self%nkpt_x * self%nkpt_y * self % nkpt_z
      self%nkpt = count_3D_BZ(total,self%a1,self%a2,self%a3)
    else
      total = self%nkpt_x * self%nkpt_y
      self%nkpt = count_2D_BZ(total,self%a1,self%a2)
    end if
    return
  end subroutine count_BrillouinZone

  subroutine deallocate_BrillouinZone(self)
     implicit none
     class(BrillouinZone) :: self
     if(allocated(self%kp)) deallocate(self%kp)
     if(allocated(self%w)) deallocate(self%w)
     return
  end subroutine deallocate_BrillouinZone

  logical function isAlloc_BrillouinZone(self)
     implicit none
     class(BrillouinZone) :: self
     if(allocated(self%kp)) then
        isAlloc_BrillouinZone = .true.
     else
        isAlloc_BrillouinZone = .false.
     end if
     return
  end function isAlloc_BrillouinZone

  subroutine generate_3D_BZ(self)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    implicit none

    class(BrillouinZone) :: self
    real(double) :: vol
    real(double), dimension(3,8) :: bz_vec
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) ::  diff
    real(double), allocatable :: inikbz(:,:), iniwkbz(:)
    real(double), allocatable, dimension(:,:) :: extrakbz
    real(double), allocatable, dimension(:) :: extrawkbz
    integer :: l,j,m, k, smallest_index, numextrakbz

    self%nkpt = self%nkpt_x * self%nkpt_y * self%nkpt_z
    if(self%nkpt < 1000) then
      allocate(extrakbz(3,self%nkpt*10), extrawkbz(self%nkpt*10))
    else
      allocate(extrakbz(3,self%nkpt/2), extrawkbz(self%nkpt/2))
    end if
    allocate(iniwkbz(self%nkpt), inikbz(3, self%nkpt) )

    vol = tpi / dot_product(self%a1, cross(self%a2,self%a3))
    self%b1 = vol * cross(self%a2, self%a3)
    self%b2 = vol * cross(self%a3, self%a1)
    self%b3 = vol * cross(self%a1, self%a2)

    bz_vec(1:3,1) = 0.d0
    bz_vec(1:3,2) = self%b1
    bz_vec(1:3,3) = self%b2
    bz_vec(1:3,4) = self%b1 + self%b2
    bz_vec(1:3,5) = self%b3
    bz_vec(1:3,6) = self%b1 + self%b3
    bz_vec(1:3,7) = self%b2 + self%b3
    bz_vec(1:3,8) = self%b1 + self%b2 + self%b3

    !Generate k-points in the paralelogram determined by b1 and b2
    iniwkbz = 1.d0
    m = 0
    do l = 1, self%nkpt_x
       do j = 1, self%nkpt_y
          do k = 1, self%nkpt_z !Run over the entire BZ
             m = m + 1
             inikbz(:, m)= dble(l-1)*self%b1 / dble(self%nkpt_x) + dble(j-1)*self%b2 / dble(self%nkpt_y) + dble(k-1)*self%b3 / dble(self%nkpt_z)
          end do
       end do
    end do

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(self%b1 + self%b2 + self%b3, self%b1 + self%b2 + self%b3))
    numextrakbz=0
    !Run over all the kpoints generated initially.
    do l=1, self%nkpt

       smallest_dist=ini_smallest_dist
       m=0
       !Checks to each of the 4 BZ's the kpoint belongs by checking
       ! to each BZ it's closer.
       do j=1, 8
          diff = inikbz(:,l) - bz_vec(:,j)
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
          diff=inikbz(:,l)-bz_vec(:,j)
          distance=sqrt(dot_product(diff,diff))
          if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. &
               j/=smallest_index ) then
             m=m+1
             numextrakbz=numextrakbz+1
             extrakbz(:,numextrakbz)=inikbz(:,l)-bz_vec(:,j)
          end if
       end do
       if(m/=0) then
          !The weight of the kpoint in the border is shared with
          ! its clones.
          iniwkbz(l)=1.d0/dble(m+1)
          extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/dble(m+1)
       end if
       !Translate the kpoint to the 1st BZ
       inikbz(:,l)=inikbz(:,l)-bz_vec(:,smallest_index)
    end do

    allocate( self%w(self%nkpt+numextrakbz),self%kp(3,self%nkpt+numextrakbz) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    self%kp(1:3, 1:self%nkpt)=inikbz (1:3,:)
    self%w(1:self%nkpt)  =iniwkbz(:)
    if(numextrakbz/=0) then
       self%kp(:,self%nkpt+1:self%nkpt+numextrakbz) = extrakbz(:,1:numextrakbz)
       self%w(self%nkpt+1:self%nkpt+numextrakbz) = extrawkbz(1:numextrakbz)
    end if
    self%nkpt = self%nkpt + numextrakbz
    self%w = self%w / dble(self%nkpt)

    deallocate(extrakbz, extrawkbz)
    deallocate(iniwkbz, inikbz)

    return
  end subroutine  generate_3D_BZ

  subroutine gen3DFraction(self, first, last, total)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    use mod_mpi_pars, only: abortProgram
    implicit none

    class(FractionalBrillouinZone) :: self
    integer*8, intent(in) :: first, last, total
    real(double) :: vol
    real(double), dimension(3,8) :: bz_vec
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) ::  diff
    real(double), dimension(3) :: kp, b1, b2, b3
    integer :: nkpt
    integer :: l, j
    integer :: nx, ny, nz
    integer :: count, added, weight, range

    allocate( self%w(self%workload), self%kp(3,self%workload) )
    self%w = 1.d0
    nkpt = self%nkpt_x * self%nkpt_y * self%nkpt_z
    vol = tpi / dot_product(self%a1, cross(self%a2,self%a3))
    b1 = vol * cross(self%a2, self%a3)
    b2 = vol * cross(self%a3, self%a1)
    b3 = vol * cross(self%a1, self%a2)

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
      if(count > self%last) exit
      weight = 0
      range = 0
      nx = mod(floor(dble(l-1) / dble(self%nkpt_y * self%nkpt_z)), self%nkpt_x)
      ny = mod(floor(dble(l-1) / dble(self%nkpt_z)), self%nkpt_y)
      nz = mod(l-1, self%nkpt_z)
      kp = dble(nx)*b1 / dble(self%nkpt_x) + dble(ny)*b2 / dble(self%nkpt_y) + dble(nz)*b3 / dble(self%nkpt_z)

      smallest_dist = ini_smallest_dist
      !Checks to each of the 4 BZ's the kpoint belongs by checking
      ! to each BZ it's closer.
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
    self%w = self%w / dble(count)
    print *, added, self%workload
    if(added > self%workload) call abortProgram("[gen3DFraction] Generated more points than it should have!")
    if(added < self%workload) call abortProgram("[gen3DFraction] Generated less points than it should have!")
    return
  end subroutine gen3DFraction

  integer function count_3D_BZ(nkpt_in, a1, a2, a3)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    implicit none

    integer, intent(in) :: nkpt_in
    real(double), dimension(3), intent(in) :: a1, a2, a3
    real(double) :: vol
    real(double), dimension(3,8) :: bz_vec
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) ::  diff
    real(double), dimension(3) :: kp, b1, b2, b3
    integer :: nkpt_x, nkpt_y, nkpt_z, nkpt
    integer :: l,j, smallest_index, numextrakbz
    integer :: nkpt_perdim !n. of k point per dimension
    integer :: nx, ny, nz

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
    numextrakbz=0
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
      !Checks to each of the 4 BZ's the kpoint belongs by checking
      ! to each BZ it's closer.
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
        if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. j/=smallest_index ) numextrakbz=numextrakbz+1
      end do
    end do
    !$omp end parallel do
    count_3D_BZ = nkpt + numextrakbz

    return
  end function count_3D_BZ

  subroutine generate_2D_BZ(self)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    implicit none
    class(BrillouinZone) :: self
    real(double), dimension(3) :: zdir
    real(double), dimension(3,4) :: bz_vec
    real(double), dimension(3)   :: diff
    real(double), dimension(:,:), allocatable :: extrakbz
    real(double), dimension(:), allocatable :: extrawkbz
    real(double), dimension(:,:), allocatable :: inikbz
    real(double), dimension(:), allocatable :: iniwkbz
    real(double) :: smallest_dist, distance, ini_smallest_dist, vol
    integer :: j,l,m, smallest_index, numextrakbz

    zdir = [0,0,1]
    allocate(extrakbz(3,self%nkpt*10))
    allocate(extrawkbz(self%nkpt*10))

    self%nkpt = self%nkpt_x * self%nkpt_y

    allocate(iniwkbz(self%nkpt), inikbz(3,self%nkpt))

    vol = tpi / dot_product(zdir, cross(self%a1,self%a2))
    self%b1 = vol * cross(self%a1, zdir)
    self%b2 = vol * cross(zdir, self%a2)

    bz_vec(:,1) = 0.d0
    bz_vec(:,2) = self%b1
    bz_vec(:,3) = self%b2
    bz_vec(:,4) = self%b1 + self%b2

    !Generate k-points in the paralelogram determined by b1 and b2
    iniwkbz = 1.d0
    m = 0
    do l = 1, self%nkpt_x
       do j = 1, self%nkpt_y !Run over the entire BZ
          m = m + 1
          inikbz(:, m) = dble(l-1)*self%b1 / dble(self%nkpt_x) + dble(j-1)*self%b2 / dble(self%nkpt_y)
       end do
    end do

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(self%b1 + self%b2, self%b1 + self%b2))
    numextrakbz = 0
    !Run over all the kpoints generated initially.
    do l = 1, self%nkpt
       smallest_dist = ini_smallest_dist
       m = 0
       !Checks to each of the 4 BZ's the kpoint belongs by checking
       ! to each BZ it's closer.
       do j = 1, 4
          diff = inikbz(:,l) - bz_vec(:,j)
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
          diff = inikbz(:, l) - bz_vec(:, j)
          distance = sqrt(dot_product(diff, diff))
          if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. j /= smallest_index ) then
             m = m + 1
             numextrakbz = numextrakbz + 1
             extrakbz(:, numextrakbz) = inikbz(:, l) - bz_vec(:,j)
          end if
       end do
       if(m /= 0) then
          !The weight of the kpoint in the border is shared with
          ! its clones.
          iniwkbz(l) = 1.d0 / dble(m+1)
          extrawkbz(numextrakbz - m + 1 : numextrakbz) = 1.d0 / dble(m+1)
       end if
       !Translate the kpoint to the 1st BZ
       inikbz(:, l) = inikbz(:, l) - bz_vec(:, smallest_index)
    end do

    allocate( self%w(self%nkpt + numextrakbz), self%kp(3, self%nkpt + numextrakbz) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    self%kp(:, 1:self%nkpt) = inikbz (:,:)
    self%w(1:self%nkpt) = iniwkbz(:)
    if(numextrakbz/=0) then
       self%kp(:, self%nkpt + 1 : self%nkpt + numextrakbz) = extrakbz(:, 1 : numextrakbz)
       self%w(self%nkpt + 1 : self%nkpt + numextrakbz) = extrawkbz(1 : numextrakbz)
    end if
    self%nkpt = self%nkpt + numextrakbz
    self%w = self%w/dble(self%nkpt)

    deallocate(extrakbz, extrawkbz)
    deallocate(iniwkbz, inikbz)

    return
  end subroutine  generate_2D_BZ

  subroutine gen2DFraction(self, first, last, total)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    use mod_mpi_pars, only: abortProgram
    implicit none

    class(FractionalBrillouinZone) :: self
    integer*8, intent(in) :: first, last, total
    real(double) :: vol
    real(double), dimension(3,4) :: bz_vec
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) ::  diff, zdir
    real(double), dimension(3) :: kp, b1, b2
    integer :: nkpt
    integer :: l, j
    integer :: nx, ny
    integer :: count, added, weight

    zdir = [0.0,0.0,1.0]
    allocate( self%w(self%workload), self%kp(3,self%workload) )
    self%w = 1.d0
    nkpt = self%nkpt_x * self%nkpt_y

    vol = tpi / dot_product(zdir, cross(self%a1,self%a2))
    b1 = vol * cross(self%a1, zdir)
    b2 = vol * cross(zdir, self%a2)

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
      if(count > self%last) exit
      weight = 0
      nx = mod(floor(dble(l-1) / dble(self%nkpt_y)), self%nkpt_x)
      ny = mod(l-1, self%nkpt_y)
      kp = dble(nx)*b1 / dble(self%nkpt_x) + dble(ny)*b2 / dble(self%nkpt_y)

      smallest_dist = ini_smallest_dist
      !Checks to each of the 4 BZ's the kpoint belongs by checking
      ! to each BZ it's closer.
      do j = 1, 4
        diff = kp - bz_vec(:,j)
        if(sqrt(dot_product(diff, diff)) < smallest_dist) smallest_dist = distance
      end do

      !Checks if the kpoint is in the border between two or more
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
            self%kp(:,added) = diff
          end if
        end if
      end do
      self%w(added-weight+1:added) = 1.0 / dble(weight)
    end do
    self%w = self%w / dble(count)
    if(added > self%workload) call abortProgram("[gen2DFraction] Generated more points than it should have!")
    if(added < self%workload) call abortProgram("[gen2DFraction] Generated less points than it should have!")
    return
  end subroutine gen2DFraction

  integer function count_2D_BZ(nkpt_in, a1, a2)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    implicit none

    integer, intent(in) :: nkpt_in
    real(double), dimension(3), intent(in) :: a1,a2
    real(double), dimension(3) :: kp, b1, b2
    real(double), dimension(3) :: zdir
    real(double), dimension(3,4) :: bz_vec
    real(double), dimension(3)   :: diff
    real(double) :: smallest_dist, distance, ini_smallest_dist, vol
    integer :: j,l, smallest_index, numextrakbz, nkpt_perdim
    integer :: nkpt_x, nkpt_y, nx, ny, nkpt

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
      !Checks to each of the 4 BZ's the kpoint belongs by checking
      ! to each BZ it's closer.
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
    return
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

    return
  end subroutine output_kpoints

end module mod_BrillouinZone
