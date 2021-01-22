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
  use mod_kind, only: dp, int32, int64
  ! use MPI_f08,      only: MPI_Comm
  implicit none

  type :: BrillouinZone
    integer(int64) :: nkpt_rep = 0 ! K-points with repetitions on the edges
    integer(int64) :: nkpt     = 0 ! K-points without repetitions on the edges (weights take care of that)
    integer(int32) :: nkpt_x   = 0
    integer(int32) :: nkpt_y   = 0
    integer(int32) :: nkpt_z   = 0

    real(dp), dimension(:,:), allocatable :: kp
    real(dp), dimension(:),   allocatable :: w
    real(dp), dimension(3) :: b1, b2, b3

  contains
    procedure :: free => deallocate_BrillouinZone
    procedure :: isAlloc => isAlloc_BrillouinZone
    procedure :: countBZ
    procedure :: print => output_kpoints
  end type BrillouinZone

  type, extends(BrillouinZone) :: FractionalBrillouinZone
    integer(int64) :: workload = 0
    integer(int64) :: first = 0
    integer(int64) :: last = 0
    integer(int32) :: size, rank, comm
    ! MPI_f08:
    ! integer(int32) :: size, rank
    ! type(MPI_Comm) :: comm
  contains
    procedure :: setup_fraction
    procedure :: gen1DFraction
    procedure :: gen2DFraction
    procedure :: gen3DFraction
  end type FractionalBrillouinZone

  type(FractionalBrillouinZone) :: realBZ

contains

  subroutine setup_fraction(self, sys, rank, size, comm, lkpoints)
    use mod_kind,     only: int32
    use mod_mpi_pars, only: calcWorkload
    use mod_System,   only: System_type
    implicit none
    class(FractionalBrillouinZone) :: self
    type(System_type),intent(in) :: sys
    integer(int32),   intent(in) :: rank, size, comm
    logical,          intent(in), optional :: lkpoints
    ! MPI_f08:
    ! type(System_type),   intent(in) :: sys
    ! integer(int32),      intent(in) :: rank, size
    ! type(MPI_Comm), intent(in) :: comm
    ! logical,        intent(in), optional :: lkpoints

    character(len=50) :: filename

    if(allocated(self%kp)) deallocate(self%kp)
    if(allocated(self%w )) deallocate(self%w)

    self % rank = rank
    self % size = size
    self % comm = comm

    ! call self%countBZ(sys)
    call calcWorkload(self%nkpt, self%size, self%rank, self%first, self%last)
    self%workload = self%last - self%first + 1
    select case(sys%isysdim)
    case(3)
      call self%gen3DFraction(sys,self%first,self%last)
    case(2)
      call self%gen2DFraction(sys,self%first,self%last)
    case default
      call self%gen1DFraction(sys,self%first,self%last)
    end select

    if(present(lkpoints).and.lkpoints) then
      filename = "kpoints_" // trim(sys%Name)
      open (unit=3333, file=filename,status='replace')
      call self%print(3333)
      close(3333)
    end if

  end subroutine setup_fraction

  subroutine countBZ(self,sys)
    use mod_kind,   only: int64
    use mod_System, only: System_type
    implicit none
    class(BrillouinZone)     :: self
    type(System_type), intent(in) :: sys
    integer(int64)                :: total

    select case(sys%isysdim)
    case(3)
      total = self%nkpt_x * self%nkpt_y * self%nkpt_z
      call count_3D_BZ(total,sys%a1,sys%a2,sys%a3,self%nkpt,self%nkpt_rep)
    case(2)
      total = self%nkpt_x * self%nkpt_y
      call count_2D_BZ(total,sys%a1,sys%a2,self%nkpt,self%nkpt_rep)
    case default
      total = self%nkpt_x 
      call count_1D_BZ(total,sys%a1,self%nkpt,self%nkpt_rep)
    end select

  end subroutine countBZ

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
    use mod_kind,      only: dp, int32, int64
    use mod_constants, only: tpi
    use mod_tools,     only: cross, itos, vec_norm
    use mod_mpi_pars,  only: abortProgram
    use mod_System,    only: System_type
    implicit none
    class(FractionalBrillouinZone) :: self
    type(System_type), intent(in)  :: sys
    integer(int64),    intent(in)  :: first, last
    real(dp), dimension(3,8)   :: bz_vec
    real(dp), dimension(3,8)   :: diff
    real(dp), dimension(3)     :: kp, b1, b2, b3, largest
    real(dp) :: vol
    real(dp) :: smallest_dist, distance(8), ini_smallest_dist
    integer(int64)    :: nkpt, l
    integer(int32)    :: nx, ny, nz
    integer(int64)    :: kount, added, weight, range
    integer(int32)    :: j

    allocate( self%w(self%workload), self%kp(3,self%workload) )
    self%w = 1._dp
    nkpt = self%nkpt_x * self%nkpt_y * self%nkpt_z
    vol  = tpi / dot_product(sys%a1, cross(sys%a2,sys%a3))
    b1   = vol * cross(sys%a2,sys%a3)
    b2   = vol * cross(sys%a3,sys%a1)
    b3   = vol * cross(sys%a1,sys%a2)

    bz_vec(1:3,1) = 0._dp
    bz_vec(1:3,2) = b1
    bz_vec(1:3,3) = b2
    bz_vec(1:3,4) = b1 + b2
    bz_vec(1:3,5) = b3
    bz_vec(1:3,6) = b1 + b3
    bz_vec(1:3,7) = b2 + b3
    bz_vec(1:3,8) = b1 + b2 + b3

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2+b3|, bigger than the distance of any genarated kpoint
    largest = b1 + b2 + b3
    ini_smallest_dist = 10._dp * vec_norm(largest, 3)
    kount = 0
    added = 0
    !Run over all the kpoints generated initially.
    do l = 1, nkpt
      if(kount > last) exit
      weight = 0
      range = 0
      nx = int(mod(floor(dble(l-1) / dble(self%nkpt_y*self%nkpt_z),kind(nx)), self%nkpt_x),kind(nx))
      ny = int(mod(floor(dble(l-1) / dble(self%nkpt_z),kind(ny)), self%nkpt_y),kind(ny))
      nz = int(mod(l-1, self%nkpt_z),kind(nz))
      kp = dble(nx)*b1 / dble(self%nkpt_x) + dble(ny)*b2 / dble(self%nkpt_y) + dble(nz)*b3 / dble(self%nkpt_z)

      smallest_dist = ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer and stores the smallest_dist
      do j = 1, 8
        diff(:,j) = kp - bz_vec(:,j)
        distance(j) = vec_norm(diff(:,j), 3)
        if(distance(j) < smallest_dist) smallest_dist = distance(j)
      end do

      ! Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 8
        if( abs(distance(j)-smallest_dist) < 1.e-12_dp ) then
          kount = kount + 1
          weight = weight + 1
          if(kount >= first .and. kount <= last ) then
            added = added + 1
            range = range + 1
            self%kp(:,added) = diff(:,j)
          end if
        end if
      end do
      self%w(added-range+1:added) = 1._dp / dble(weight)
    end do
    self%w = self%w / dble(nkpt)

    if(added == self%workload) return

    if(added > self%workload) then
      call abortProgram("[gen3DFraction] Generated more points than it should have! added = " // trim(itos(added)) // ", self%workload = " // trim(itos(self%workload)))
    else
      call abortProgram("[gen3DFraction] Generated less points than it should have! added = " // trim(itos(added)) // ", self%workload = " // trim(itos(self%workload)))
    end if
 end subroutine gen3DFraction

  subroutine count_3D_BZ(nkpt_in, a1, a2, a3, numextrakbz, nkpt)
    use mod_kind,      only: dp, int32, int64
    use mod_constants, only: tpi
    use mod_tools,     only: cross, vec_norm
    implicit none
    integer(int64),         intent(in)  :: nkpt_in
    real(dp), dimension(3), intent(in)  :: a1, a2, a3
    integer(int64),         intent(out) :: numextrakbz, nkpt
    real(dp), dimension(3)   :: kp, b1, b2, b3, largest
    real(dp), dimension(3,8) :: bz_vec
    real(dp), dimension(3,8) :: diff
    real(dp) :: smallest_dist, distance(8), ini_smallest_dist, vol
    integer(int32)    :: j, nkpt_x, nkpt_y, nkpt_z, nx, ny, nz, nkpt_perdim
    integer(int64)    :: l

    nkpt_perdim = ceiling((dble(nkpt_in))**(1._dp/3._dp))
    nkpt_x = nkpt_perdim
    nkpt_y = nkpt_perdim
    nkpt_z = nkpt_perdim

    nkpt = nkpt_x * nkpt_y * nkpt_z

    vol = tpi / dot_product(a1, cross(a2,a3))
    b1  = vol * cross(a2, a3)
    b2  = vol * cross(a3, a1)
    b3  = vol * cross(a1, a2)

    bz_vec(1:3,1) = 0._dp
    bz_vec(1:3,2) = b1
    bz_vec(1:3,3) = b2
    bz_vec(1:3,4) = b1 + b2
    bz_vec(1:3,5) = b3
    bz_vec(1:3,6) = b1 + b3
    bz_vec(1:3,7) = b2 + b3
    bz_vec(1:3,8) = b1 + b2 + b3

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2+b3|, bigger than the distance of any genarated kpoint
    largest = b1 + b2 + b3
    ini_smallest_dist = 10._dp * vec_norm(largest, 3)
    numextrakbz = 0
    !Run over all the kpoints generated initially.
    !$omp parallel do default(none) reduction(+:numextrakbz) &
    !$omp& private(l, j, nx, ny, nz, kp, smallest_dist, diff, distance) &
    !$omp& shared(nkpt, ini_smallest_dist, bz_vec, b1, b2, b3, nkpt_x, nkpt_y, nkpt_z)
    do l=1, nkpt
      nx = int(mod(floor(dble(l-1) / dble(nkpt_y*nkpt_z),kind(nx)), nkpt_x),kind(nx))
      ny = int(mod(floor(dble(l-1) / dble(nkpt_z),kind(ny)), nkpt_y),kind(ny))
      nz = int(mod(l-1, nkpt_z),kind(nz))
      kp = dble(nx)*b1 / dble(nkpt_x) + dble(ny)*b2 / dble(nkpt_y) + dble(nz)*b3 / dble(nkpt_z)

      smallest_dist=ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer and stores the smallest_dist
      do j = 1, 8
        diff(:,j) = kp - bz_vec(:,j)
        distance(j) = vec_norm(diff(:,j), 3)
        if(distance(j) < smallest_dist) smallest_dist = distance(j)
      end do

      ! Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 8
        if( abs(distance(j)-smallest_dist) < 1.e-12_dp ) numextrakbz=numextrakbz+1
      end do
    end do
    !$omp end parallel do
  end subroutine count_3D_BZ

  subroutine gen2DFraction(self,sys,first,last)
    use mod_kind,      only: dp, int32, int64
    use mod_constants, only: tpi
    use mod_tools,     only: cross, itos, vec_norm
    use mod_mpi_pars,  only: abortProgram
    use mod_System,    only: System_type
    implicit none
    class(FractionalBrillouinZone) :: self
    type(System_type), intent(in)  :: sys
    integer(int64),    intent(in)  :: first, last
    real(dp), dimension(3,4)   :: bz_vec
    real(dp), dimension(3,4)   :: diff
    real(dp), dimension(3)     :: kp, zdir, b1, b2, largest
    real(dp) :: vol
    real(dp) :: smallest_dist, distance(4), ini_smallest_dist
    integer(int64)    :: nkpt, l
    integer(int32)    :: nx, ny
    integer(int64)    :: kount, added, weight, range
    integer(int32)    :: j

    zdir = [0._dp,0._dp,1._dp]
    allocate( self%w(self%workload), self%kp(3,self%workload) )
    self%w = 1._dp
    nkpt = self%nkpt_x * self%nkpt_y
    vol  = tpi/dot_product(zdir, cross(sys%a1,sys%a2))
    b1   = vol*cross(sys%a1,zdir)
    b2   = vol*cross(zdir,sys%a2)

    bz_vec(:,1) = 0._dp
    bz_vec(:,2) = b1
    bz_vec(:,3) = b2
    bz_vec(:,4) = b1 + b2

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    largest = b1 + b2
    ini_smallest_dist = 10._dp * vec_norm(largest, 3)
    kount = 0
    added = 0
    !Run over all the kpoints generated initially.
    do l=1, nkpt
      if(kount > last) exit
      weight = 0
      range = 0
      nx = int(mod(floor(dble(l-1) / dble(self%nkpt_y),kind(nx)), self%nkpt_x),kind(nx))
      ny = int(mod(l-1, self%nkpt_y),kind(ny))
      kp = dble(nx)*b1 / dble(self%nkpt_x) + dble(ny)*b2 / dble(self%nkpt_y)

      smallest_dist = ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer and stores the smallest_dist
      do j = 1, 4
        diff(:,j) = kp - bz_vec(:,j)
        distance(j) = vec_norm(diff(:,j), 3)
        if(distance(j) < smallest_dist) smallest_dist = distance(j)
      end do

      ! Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 4
        if( abs(distance(j)-smallest_dist) < 1.e-12_dp ) then
          kount = kount + 1
          weight = weight + 1
          if(kount >= first .and. kount <= last ) then
            added = added + 1
            range = range + 1
            self%kp(:,added) = diff(:,j)
          end if
        end if
      end do
      self%w(added-range+1:added) = 1.0 / dble(weight)
    end do
    self%w = self%w / dble(nkpt)

    if(added == self%workload) return

    if(added > self%workload) then
      call abortProgram("[gen2DFraction] Generated more points than it should have! added = " // trim(itos(added)) // ", self%workload = " // trim(itos(self%workload)))
    else
      call abortProgram("[gen2DFraction] Generated less points than it should have! added = " // trim(itos(added)) // ", self%workload = " // trim(itos(self%workload)))
    end if
  end subroutine gen2DFraction

  subroutine count_2D_BZ(nkpt_in, a1, a2, numextrakbz, nkpt)
    use mod_kind,      only: dp, int32, int64
    use mod_constants, only: tpi
    use mod_tools,     only: cross, vec_norm
    implicit none
    integer(int64),         intent(in)  :: nkpt_in
    real(dp), dimension(3), intent(in)  :: a1,a2
    integer(int64),         intent(out) :: numextrakbz, nkpt
    real(dp), dimension(3)   :: kp, b1, b2, largest
    real(dp), dimension(3)   :: zdir
    real(dp), dimension(3,4) :: bz_vec
    real(dp), dimension(3,4) :: diff
    real(dp) :: smallest_dist, distance(4), ini_smallest_dist, vol
    integer(int32)    :: j, nkpt_x, nkpt_y, nx, ny, nkpt_perdim
    integer(int64)    :: l

    zdir = [0._dp,0._dp,1._dp]
    nkpt_perdim = ceiling(sqrt(dble(nkpt_in)))
    nkpt_x = nkpt_perdim
    nkpt_y = nkpt_perdim

    nkpt = nkpt_x * nkpt_y

    vol = tpi / dot_product(zdir, cross(a1,a2))
    b1  = vol * cross(a1, zdir)
    b2  = vol * cross(zdir, a2)

    bz_vec(:,1) = 0._dp
    bz_vec(:,2) = b1
    bz_vec(:,3) = b2
    bz_vec(:,4) = b1 + b2

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    largest = b1 + b2
    ini_smallest_dist = 10._dp * vec_norm(largest, 3)
    numextrakbz = 0
    !Run over all the kpoints generated initially.
    !$omp parallel do default(none) reduction(+:numextrakbz) &
    !$omp& private(l,j,nx,ny,kp,smallest_dist, diff, distance) &
    !$omp& shared(nkpt, ini_smallest_dist, bz_vec, b1, b2, nkpt_x, nkpt_y)
    do l = 1, nkpt
      nx = int(mod(floor(dble(l-1) / dble(nkpt_y),kind(nx)), nkpt_x),kind(nx))
      ny = int(mod(l-1, nkpt_y),kind(ny))
      kp = dble(nx)*b1 / dble(nkpt_x) + dble(ny)*b2 / dble(nkpt_y)

      smallest_dist = ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer and stores the smallest_dist
      do j = 1, 4
        diff(:,j) = kp - bz_vec(:,j)
        distance(j) = vec_norm(diff(:,j), 3)
        if(distance(j) < smallest_dist) smallest_dist = distance(j)
      end do

      !Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j = 1, 4
        if( abs(distance(j)-smallest_dist) < 1.e-12_dp ) numextrakbz = numextrakbz + 1
      end do
    end do
    !$omp end parallel do
  end subroutine count_2D_BZ

  subroutine gen1DFraction(self,sys,first,last)
    use mod_kind,      only: dp, int32, int64
    use mod_constants, only: tpi
    use mod_tools,     only: cross, itos, vec_norm
    use mod_mpi_pars,  only: abortProgram
    use mod_System,    only: System_type
    implicit none
    class(FractionalBrillouinZone) :: self
    type(System_type), intent(in)  :: sys
    integer(int64),    intent(in)  :: first, last
    real(dp), dimension(3,2)   :: bz_vec
    real(dp), dimension(3,2)   :: diff
    real(dp), dimension(3)     :: kp, b1, zdir, ydir
    real(dp) :: vol
    real(dp) :: smallest_dist, distance(2), ini_smallest_dist
    integer(int64)    :: nkpt, l
    integer(int32)    :: nx
    integer(int64)    :: kount, added, weight, range
    integer(int32)    :: j

    zdir = [0._dp,0._dp,1._dp]
    ydir = [0._dp,1._dp,0._dp]
    allocate( self%w(self%workload), self%kp(3,self%workload) )
    self%w = 1._dp
    nkpt = self%nkpt_x
    vol  = tpi/dot_product(zdir, cross(sys%a1,ydir))
    b1   = vol*cross(zdir,ydir)

    bz_vec(:,1) = 0._dp
    bz_vec(:,2) = b1

    !Translate the k-points to the 1st BZ.
    !10*|b1|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10._dp * vec_norm(b1, 3)
    kount = 0
    added = 0
    !Run over all the kpoints generated initially.
    do l=1, nkpt
      if(kount > last) exit
      weight = 0
      range = 0
      nx = int(mod(l-1, self%nkpt_x),kind(nx))
      kp = dble(nx)*b1 / dble(self%nkpt_x)

      smallest_dist = ini_smallest_dist
      ! Checks to which of the 2 BZ's the kpoint belongs by checking
      ! to which BZ it's closer and stores the smallest_dist
      do j = 1, 2
        diff(:,j) = kp - bz_vec(:,j)
        distance(j) = vec_norm(diff(:,j), 3)
        if(distance(j) < smallest_dist) smallest_dist = distance(j)
      end do

      ! Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 2
        if( abs(distance(j)-smallest_dist) < 1.e-12_dp ) then
          kount = kount + 1
          weight = weight + 1
          if(kount >= first .and. kount <= last ) then
            added = added + 1
            range = range + 1
            self%kp(:,added) = diff(:,j)
          end if
        end if
      end do
      self%w(added-range+1:added) = 1.0 / dble(weight)
    end do
    self%w = self%w / dble(nkpt)

    if(added == self%workload) return

    if(added > self%workload) then
      call abortProgram("[gen1DFraction] Generated more points than it should have! added = " // trim(itos(added)) // ", self%workload = " // trim(itos(self%workload)))
    else
      call abortProgram("[gen1DFraction] Generated less points than it should have! added = " // trim(itos(added)) // ", self%workload = " // trim(itos(self%workload)))
    end if
  end subroutine gen1DFraction

  subroutine count_1D_BZ(nkpt_in, a1, numextrakbz, nkpt)
    use mod_kind,      only: dp, int32, int64
    use mod_constants, only: tpi
    use mod_tools,     only: cross, vec_norm
    implicit none
    integer(int64),                  intent(in)  :: nkpt_in
    real(dp), dimension(3), intent(in)  :: a1
    integer(int64),                  intent(out) :: numextrakbz, nkpt
    real(dp), dimension(3)   :: kp, b1
    real(dp), dimension(3)   :: zdir, ydir
    real(dp), dimension(3,2) :: bz_vec
    real(dp), dimension(3,2) :: diff
    real(dp) :: smallest_dist, distance(2), ini_smallest_dist, vol
    integer(int32)    :: j, nkpt_x, nkpt_perdim, nx
    integer(int64)    :: l

    zdir = [0._dp,0._dp,1._dp]
    ydir = [0._dp,1._dp,0._dp]
    nkpt_perdim = ceiling(dble(nkpt_in)) 
    nkpt_x = nkpt_perdim

    nkpt = nkpt_x

    vol = tpi / dot_product(zdir, cross(a1,ydir))
    b1  = vol * cross(zdir, ydir)

    bz_vec(:,1) = 0._dp
    bz_vec(:,2) = b1

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10._dp * vec_norm(b1, 3)
    numextrakbz = 0
    !Run over all the kpoints generated initially.
    !$omp parallel do default(none) reduction(+:numextrakbz) &
    !$omp& private(l,j,nx,kp,smallest_dist, diff, distance) &
    !$omp& shared(nkpt, ini_smallest_dist, bz_vec, b1, nkpt_x)
    do l = 1, nkpt
      nx = int(mod(l-1, nkpt_x),kind(nx))
      kp = dble(nx)*b1 / dble(nkpt_x)

      smallest_dist = ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer and stores the smallest_dist
      do j = 1, 2
        diff(:,j) = kp - bz_vec(:,j)
        distance(j) = vec_norm(diff(:,j), 3)
        if(distance(j) < smallest_dist) smallest_dist = distance(j)
      end do
      !Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j = 1, 2
        if( abs(distance(j)-smallest_dist) < 1.e-12_dp ) numextrakbz = numextrakbz + 1
      end do
    end do
    !$omp end parallel do
  end subroutine count_1D_BZ
  

  subroutine output_kpoints(self,unit)
    use mod_kind, only: int32, int64
    implicit none
    class(BrillouinZone) :: self
    integer(int64) :: i
    integer(int32) :: unit

    write(unit=unit,fmt="(a)") ' #      kx            ky            kz            wk'
    do i=1,self%nkpt
       write(unit=unit,fmt="(4(f12.9,2x))") self%kp(1:3,i),self%w(i)
    end do

  end subroutine output_kpoints


! This subroutine generates all the k-points
! and store the different values of the
! component "component" (1,2,3)
  subroutine store_diff(nkpt_in, b1, b2, b3, component, ndiffk, diff_k_unsrt)
    use mod_kind,       only: dp, int32, int64
    use mod_tools,      only: itos
    use mod_mpi_pars,   only: abortProgram, myrank
    use mod_parameters, only: kptotal_in
    !$ use omp_lib
    implicit none

    integer(int64),  intent(in) :: nkpt_in
    !! Initial kpoints (as given by realBZ%nkpt_x * realBZ%nkpt_y * realBZ%nkpt_z)
    real(dp), dimension(3), intent(in)  :: b1, b2, b3
    !! Reciprocal vectors
    integer,   intent(in)  :: component
    !! Which component to store the different values and the map
    integer,   intent(out) :: ndiffk
    !! Number of different kp(component)
    real(dp), dimension(:), allocatable, intent(out) :: diff_k_unsrt
    !! Different values of kp(component)
    real(dp), dimension(:), allocatable              :: diff_k_loc,diff_k_temp,diff_k_temp2
    !! Local and temporary array for different values of kp(component)
    real(dp), dimension(3,8) :: bz_vec
    real(dp) :: smallest_dist, distance, ini_smallest_dist
    real(dp), dimension(3) :: diff
    real(dp), dimension(3) :: kp
    integer, dimension(:), allocatable :: ndiffk_loc
    integer(int64)  :: l, nkpt
    integer(int32)  :: nkpt_perdim, nkpt_x, nkpt_y, nkpt_z
    !! Number of k points per dimension
    integer(int32)  :: nx, ny, nz , start, end
    integer(int32)  :: maxdiffk, ndiffk_max, kount
    !! Number of different k points, maximum value, and counter
    integer(int32)  :: j, mythread, nthreads

    nkpt_perdim = ceiling((dble(nkpt_in))**(1._dp/3._dp))
    nkpt_x = nkpt_perdim
    nkpt_y = nkpt_perdim
    nkpt_z = nkpt_perdim

    nkpt = int(nkpt_x * nkpt_y * nkpt_z,kind(nkpt))

    nthreads = 1
    !$ nthreads = omp_get_max_threads()
    ndiffk_max = 2*nkpt_perdim
    allocate( ndiffk_loc(0:nthreads) )
    ndiffk_loc(:) = 0
    ndiffk_loc(0) = 1

    ! Counter of total number of points
    kount = 0

    bz_vec(1:3,1) = 0._dp
    bz_vec(1:3,2) = b1
    bz_vec(1:3,3) = b2
    bz_vec(1:3,4) = b1 + b2
    bz_vec(1:3,5) = b3
    bz_vec(1:3,6) = b1 + b3
    bz_vec(1:3,7) = b2 + b3
    bz_vec(1:3,8) = b1 + b2 + b3

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10._dp * sqrt(dot_product(b1 + b2 + b3, b1 + b2 + b3))
    !Run over all the kpoints generated initially.
    !$omp parallel default(none) &
    !$omp& private(l, j, nx, ny, nz, kp, smallest_dist, diff, distance, diff_k_loc, start, end, mythread) &
    !$omp& shared(kount, ndiffk_max, nkpt, ini_smallest_dist, bz_vec, b1, b2, b3, nkpt_x, nkpt_y, nkpt_z, component, ndiffk_loc, diff_k_temp, maxdiffk)

    mythread = 1
    !$ mythread = omp_get_thread_num()+1

    !! Allocating large temporary vector to store different values of kp(component)
    allocate( diff_k_loc(ndiffk_max) )
    diff_k_loc(:) = 999._dp

    !$omp do schedule(dynamic) reduction(+:kount)
    do l=1, nkpt
      nx = int(mod(floor(dble(l-1) / dble(nkpt_y*nkpt_z),kind(nx)), nkpt_x),kind(nx))
      ny = int(mod(floor(dble(l-1) / dble(nkpt_z),kind(ny)), nkpt_y),kind(nx))
      nz = int(mod(l-1, nkpt_z),kind(nz))
      kp = dble(nx)*b1 / dble(nkpt_x) + dble(ny)*b2 / dble(nkpt_y) + dble(nz)*b3 / dble(nkpt_z)

      smallest_dist=ini_smallest_dist
      ! Checks to which of the 4 BZ's the kpoint belongs by checking
      ! to which BZ it's closer.
      do j=1, 8
        diff = kp - bz_vec(:,j)
        distance = sqrt(dot_product(diff, diff))
        if(distance < smallest_dist) smallest_dist = distance
      end do
      !Checks if the kpoint is in the border between two or more
      ! BZ's. If yes, create a clone of it to translate later into
      ! the 1st BZ.
      do j=1, 8
        diff=kp - bz_vec(:,j)
        distance=sqrt(dot_product(diff,diff))
        if( abs(distance-smallest_dist) < 1.e-12_dp ) then
          kount = kount + 1 ! Counting the total number of points generated
          ! Checking if all numbers are different than current local list (if any of them is zero, the point is already stored)
          if( all(abs( diff_k_loc(:) - abs(diff(component)) )> 1.e-12_dp,1)  ) then
            ndiffk_loc(mythread) = ndiffk_loc(mythread) + 1
            diff_k_loc(ndiffk_loc(mythread)) = abs(diff(component))
          end if
        end if
      end do
    end do
    !$omp end do

    if(mythread==1) then
      maxdiffk = sum(ndiffk_loc(:))
      allocate( diff_k_temp(maxdiffk) )
    end if
    start = sum(ndiffk_loc(0:mythread-1))
    end   = start + ndiffk_loc(mythread)
    !$omp barrier
    diff_k_temp(start:end) = diff_k_loc(1:ndiffk_loc(mythread))

    !$omp end parallel

    if(kount/=realBZ%nkpt) &
      call abortProgram("[store_diff] Generated different number of points than it should have! kount = " // trim(itos(kount)) // ", nkpt = " // trim(itos(realBZ%nkpt)))

    ! Checking for different kp(component) between the locally generated ones
    allocate( diff_k_temp2(maxdiffk) )
    diff_k_temp2(:) = 999._dp
    ndiffk = 0
    do j=1,maxdiffk

      if( all(abs( diff_k_temp2(:) - diff_k_temp(j) )> 1.e-12_dp,1)  ) then
        ndiffk = ndiffk + 1
        diff_k_temp2(ndiffk) = diff_k_temp(j)
      end if

    end do

    allocate(diff_k_unsrt(ndiffk))
    diff_k_unsrt(:) = diff_k_temp2(1:ndiffk)

    ! Writing different k(component) into file
    if(myrank==0) then
      open (unit=9999, file='diffkz_'//trim(itos(kptotal_in)),status='replace')
      write(unit=9999,fmt="(i0)") ndiffk
      do l=1,ndiffk
        write(unit=9999,fmt="(es23.16)") diff_k_unsrt(l)
      end do
      close(9999)
    end if

  end subroutine store_diff

end module mod_BrillouinZone
