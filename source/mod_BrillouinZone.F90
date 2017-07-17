module mod_BrillouinZone
implicit none

contains
  subroutine setup_BrillouinZone(s)
    use mod_system, only: System
    implicit none
    type(System), intent(inout) :: s

    if(s%lbulk) then
      call generate_3D_BZ(s)
    else
      call generate_2D_BZ(s)
    end if
    return
  end subroutine setup_BrillouinZone

  subroutine generate_3D_BZ(s)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    use mod_System, only: System
    implicit none
    type(System), intent(inout) :: s

    real(double) :: vol
    real(double), dimension(3,8) :: BZ
    real(double) :: smallest_dist, distance, ini_smallest_dist
    real(double), dimension(3) ::  diff
    real(double), allocatable :: inikbz(:,:), iniwkbz(:)
    real(double), allocatable, dimension(:,:) :: extrakbz
    real(double), allocatable, dimension(:) :: extrawkbz
    integer :: l,j,m, k, smallest_index, numextrakbz
    integer :: nkpt_perdim !n. of k point per dimension

    allocate(extrakbz(3,s%nkpt*10), extrawkbz(s%nkpt*10))

    nkpt_perdim=ceiling((dble(s%nkpt))**(1.d0/3.d0))
    nkpt_perdim=nkpt_perdim
    s%nkpt = nkpt_perdim**3
    allocate( iniwkbz(s%nkpt), inikbz(3, s%nkpt) )

    vol = tpi / dot_product(s%a1, cross(s%a2,s%a3))
    s%b1 = vol * cross(s%a2, s%a3)
    s%b2 = vol * cross(s%a3, s%a1)
    s%b3 = vol * cross(s%a1, s%a2)

    BZ(:,1) = 0.d0
    BZ(:,2) = s%b1
    BZ(:,3) = s%b2
    BZ(:,4) = s%b1 + s%b2
    BZ(:,5) = s%b3
    BZ(:,6) = s%b1 + s%b3
    BZ(:,7) = s%b2 + s%b3
    BZ(:,8) = s%b1 + s%b2 + s%b3

    !Generate k-points in the paralelogram determined by b1 and b2
    iniwkbz = 1.d0
    m = 0
    do l = 1, nkpt_perdim
       do j = 1, nkpt_perdim
          do k = 1, nkpt_perdim !Run over the entire BZ
             m = m + 1
             inikbz(:, m)= ( dble(l-1)*s%b1 + dble(j-1)*s%b2 + dble(k-1)*s%b3)/ dble(nkpt_perdim)
          end do
       end do
    end do

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist=10.d0*sqrt(dot_product(s%b1+s%b2,s%b1+s%b2))
    numextrakbz=0
    !Run over all the kpoints generated initially.
    do l=1, s%nkpt
       smallest_dist=ini_smallest_dist
       m=0
       !Checks to each of the 4 BZ's the kpoint belongs by checking
       ! to each BZ it's closer.
       do j=1, 8
          diff=inikbz(:,l)-BZ(:,j)
          distance=sqrt(dot_product(diff,diff))
          if(distance<smallest_dist) then
             smallest_dist=distance
             smallest_index=j
          end if
       end do
       !Checks if the kpoint is in the border between two or more
       ! BZ's. If yes, create a clone of it to translate later into
       ! the 1st BZ.
       do j=1, 8
          diff=inikbz(:,l)-BZ(:,j)
          distance=sqrt(dot_product(diff,diff))
          if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. &
               j/=smallest_index ) then
             m=m+1
             numextrakbz=numextrakbz+1
             extrakbz(:,numextrakbz)=inikbz(:,l)-BZ(:,j)
          end if
       end do
       if(m/=0) then
          !The weight of the kpoint in the border is shared with
          ! its clones.
          iniwkbz(l)=1.d0/dble(m+1)
          extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/dble(m+1)
       end if
       !Translate the kpoint to the 1st BZ
       inikbz(:,l)=inikbz(:,l)-BZ(:,smallest_index)
    end do

    allocate( s%wkbz(s%nkpt+numextrakbz),s%kbz(3,s%nkpt+numextrakbz) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    s%kbz (:,1:s%nkpt)=inikbz (:,:)
    s%wkbz(1:s%nkpt)  =iniwkbz(:)
    if(numextrakbz/=0) then
       s%kbz(:,s%nkpt+1:s%nkpt+numextrakbz)=extrakbz(:,1:numextrakbz)
       s%wkbz(s%nkpt+1:s%nkpt+numextrakbz) =extrawkbz(1:numextrakbz)
    end if
    s%wkbz=s%wkbz/dble(s%nkpt)
    s%nkpt=s%nkpt + numextrakbz

    deallocate(extrakbz, extrawkbz)

    return
  end subroutine  generate_3D_BZ

  subroutine generate_2D_BZ(s)
    use mod_f90_kind, only: double
    use mod_constants, only: tpi
    use mod_tools, only: cross
    use mod_system, only: System
    implicit none
    type(System), intent(inout) :: s

    real(double), dimension(3,4) :: BZ
    real(double), dimension(3)   :: diff
    real(double), dimension(:,:), allocatable :: extrakbz
    real(double), dimension(:), allocatable :: extrawkbz
    real(double), dimension(:,:), allocatable :: inikbz
    real(double), dimension(:), allocatable :: iniwkbz
    real(double) :: smallest_dist, distance, ini_smallest_dist, vol
    integer :: j,l,m, smallest_index, numextrakbz, nkpt_perdim
    real(double), dimension(3) :: zdir

    zdir = [0,0,1]
    allocate(extrakbz(3,s%nkpt*10))
    allocate(extrawkbz(s%nkpt*10))
    nkpt_perdim = ceiling(sqrt(dble(s%nkpt)))
    nkpt_perdim = nkpt_perdim
    s%nkpt = nkpt_perdim**2

    allocate(iniwkbz(s%nkpt), inikbz(3,s%nkpt))

    vol = tpi / dot_product(zdir, cross(s%a1,s%a2))
    s%b1 = vol * cross(s%a1, zdir)
    s%b2 = vol * cross(zdir, s%a2)

    BZ(:,1) = 0.d0
    BZ(:,2) = s%b1
    BZ(:,3) = s%b2
    BZ(:,4) = s%b1 + s%b2

    !Generate k-points in the paralelogram determined by b1 and b2
    iniwkbz = 1.d0
    m = 0
    do l = 1, nkpt_perdim;
       do j = 1, nkpt_perdim !Run over the entire BZ
          m = m + 1
          inikbz(:, m) = ( dble(l-1)*s%b1 + dble(j-1)*s%b2 ) / dble(nkpt_perdim)
       end do
    end do

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist = 10.d0 * sqrt(dot_product(s%b1+s%b2, s%b1+s%b2))
    numextrakbz = 0
    !Run over all the kpoints generated initially.
    do l = 1, s%nkpt
       smallest_dist = ini_smallest_dist
       m = 0
       !Checks to each of the 4 BZ's the kpoint belongs by checking
       ! to each BZ it's closer.
       do j = 1, 4
          diff = inikbz(:,l) - BZ(:,j)
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
          diff = inikbz(:, l) - BZ(:, j)
          distance = sqrt(dot_product(diff, diff))
          if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. j /= smallest_index ) then
             m = m + 1
             numextrakbz = numextrakbz + 1
             extrakbz(:, numextrakbz) = inikbz(:, l) - BZ(:,j)
          end if
       end do
       if(m /= 0) then
          !The weight of the kpoint in the border is shared with
          ! its clones.
          iniwkbz(l) = 1.d0 / dble(m+1)
          extrawkbz(numextrakbz - m + 1 : numextrakbz) = 1.d0 / dble(m+1)
       end if
       !Translate the kpoint to the 1st BZ
       inikbz(:, l) = inikbz(:, l) - BZ(:, smallest_index)
    end do

    allocate( s%wkbz(s%nkpt + numextrakbz), s%kbz(3, s%nkpt + numextrakbz) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    s%kbz(:, 1:s%nkpt) = inikbz (:,:)
    s%wkbz(1:s%nkpt) = iniwkbz(:)
    if(numextrakbz/=0) then
       s%kbz(:, s%nkpt + 1 : s%nkpt + numextrakbz) = extrakbz(:, 1 : numextrakbz)
       s%wkbz(s%nkpt + 1 : s%nkpt + numextrakbz) = extrawkbz(1 : numextrakbz)
    end if
    s%wkbz = s%wkbz/dble(s%nkpt)
    s%nkpt = s%nkpt + numextrakbz
    return
  end subroutine  generate_2D_BZ


  subroutine output_kpoints(s)
    use mod_System, only: System
    implicit none
    type(System), intent(in) :: s
    integer :: i

    open (unit=3333, file='kpoint_new',status='replace')
    write(unit=3333,fmt="(a)") ' #      kx            ky            kz            wk'
    do i=1,s%nkpt
       !write(unit=2222,fmt="(3(f12.9,2x))") kbz2d(i,1),kbz2d(i,2),wkbz(i)
       write(unit=3333,fmt="(4(f12.9,2x))") s%kbz(:,i),s%wkbz(i)
    end do
    !close(2222)
    close(3333)

    return
  end subroutine output_kpoints

end module mod_BrillouinZone
