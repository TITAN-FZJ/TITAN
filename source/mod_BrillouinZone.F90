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
      integer :: nkpt
      integer :: nkpt_x = 0
      integer :: nkpt_y = 0
      integer :: nkpt_z = 0

      real(double), dimension(:,:), allocatable :: kp
      real(double), dimension(:), allocatable :: w
      real(double), dimension(3) :: b1, b2, b3
      logical :: is_bulk
   contains
      procedure :: setup => setup_BrillouinZone
      procedure, private :: generate_2d => generate_2D_BZ
      procedure, private :: generate_3d => generate_3D_BZ
      procedure :: print => output_kpoints
   end type BrillouinZone

   type(BrillouinZone) :: BZ ! Move to main at some point.


contains
  subroutine setup_BrillouinZone(self, bulk, a1, a2, a3)
    implicit none
    class(BrillouinZone) :: self
    logical, intent(in) :: bulk
    real(double), dimension(3), intent(in) :: a1, a2, a3

    if(bulk) then
      self%is_bulk = .true.
      call self % generate_3d(a1, a2, a3)
    else
      self%is_bulk = .false.
      call self % generate_2d(a1, a2)
    end if
    return
  end subroutine setup_BrillouinZone

  subroutine generate_3D_BZ(self, a1, a2, a3)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    implicit none

    class(BrillouinZone) :: self
    real(double), dimension(3), intent(in) :: a1, a2, a3
    real(double) :: vol
    real(double), dimension(3,8) :: bz_vec
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) ::  diff
    real(double), allocatable :: inikbz(:,:), iniwkbz(:)
    real(double), allocatable, dimension(:,:) :: extrakbz
    real(double), allocatable, dimension(:) :: extrawkbz
    integer :: l,j,m, k, smallest_index, numextrakbz
    integer :: nkpt_perdim !n. of k point per dimension


    if(abs(self%nkpt_x)+abs(self%nkpt_y)+abs(self%nkpt_z) == 0) then
      nkpt_perdim=ceiling((dble(self%nkpt))**(1.d0/3.d0))
      nkpt_perdim=nkpt_perdim
      self%nkpt_x = nkpt_perdim
      self%nkpt_y = nkpt_perdim
      self%nkpt_z = nkpt_perdim
    end if

    self%nkpt = self%nkpt_x * self%nkpt_y * self%nkpt_z
    allocate(extrakbz(3,self%nkpt*10), extrawkbz(self%nkpt*10))
    allocate(iniwkbz(self%nkpt), inikbz(3, self%nkpt) )

    vol = tpi / dot_product(a1, cross(a2,a3))
    self%b1 = vol * cross(a2, a3)
    self%b2 = vol * cross(a3, a1)
    self%b3 = vol * cross(a1, a2)

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
    ini_smallest_dist = 10.d0 * sqrt(dot_product(self%b1 + self%b2, self%b1 + self%b2))
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
       self%kp(:,self%nkpt+1:self%nkpt+numextrakbz)=extrakbz(:,1:numextrakbz)
       self%w(self%nkpt+1:self%nkpt+numextrakbz) =extrawkbz(1:numextrakbz)
    end if
    self%w = self%w / dble(self%nkpt)
    self%nkpt = self%nkpt + numextrakbz

    deallocate(extrakbz, extrawkbz)

    return
  end subroutine  generate_3D_BZ

  subroutine generate_2D_BZ(self, a1, a2)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    implicit none
    class(BrillouinZone) :: self
    real(double), dimension(3), intent(in) :: a1,a2
    real(double), dimension(3) :: zdir
    real(double), dimension(3,4) :: bz_vec
    real(double), dimension(3)   :: diff
    real(double), dimension(:,:), allocatable :: extrakbz
    real(double), dimension(:), allocatable :: extrawkbz
    real(double), dimension(:,:), allocatable :: inikbz
    real(double), dimension(:), allocatable :: iniwkbz
    real(double) :: smallest_dist, distance, ini_smallest_dist, vol
    integer :: j,l,m, smallest_index, numextrakbz, nkpt_perdim

    zdir = [0,0,1]
    allocate(extrakbz(3,self%nkpt*10))
    allocate(extrawkbz(self%nkpt*10))
    nkpt_perdim = ceiling(sqrt(dble(self%nkpt)))
    nkpt_perdim = nkpt_perdim
    self%nkpt = nkpt_perdim**2

    allocate(iniwkbz(self%nkpt), inikbz(3,self%nkpt))

    vol = tpi / dot_product(zdir, cross(a1,a2))
    self%b1 = vol * cross(a1, zdir)
    self%b2 = vol * cross(zdir, a2)

    bz_vec(:,1) = 0.d0
    bz_vec(:,2) = self%b1
    bz_vec(:,3) = self%b2
    bz_vec(:,4) = self%b1 + self%b2

    !Generate k-points in the paralelogram determined by b1 and b2
    iniwkbz = 1.d0
    m = 0
    do l = 1, nkpt_perdim;
       do j = 1, nkpt_perdim !Run over the entire BZ
          m = m + 1
          inikbz(:, m) = ( dble(l-1)*self%b1 + dble(j-1)*self%b2 ) / dble(nkpt_perdim)
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
    self%w = self%w/dble(self%nkpt)
    self%nkpt = self%nkpt + numextrakbz
    return
  end subroutine  generate_2D_BZ


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
