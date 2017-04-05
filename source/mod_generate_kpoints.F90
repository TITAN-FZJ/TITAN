!-------------------------------------------------------------------------------
! juTITAN
!-------------------------------------------------------------------------------
!
! MODULE: mod_lattice
!
!> @author
!> Filipe Guimaraes, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> Calculation of next nearest neighbours in plane and out of plane.
!
! REVISION HISTORY:
! 12 August 2015 - Initial Version
! 14 March 2017 - Last revision
!-------------------------------------------------------------------------------
module mod_generate_kpoints
  use mod_f90_kind, only: double
  use mod_parameters, only: a0, outputunit
  use mod_constants
  use mod_mpi_pars
  implicit none

  ! wave vector integration variables
  integer :: nkpoints                                     !> Number of generated k-points
  real(double), dimension(:,:), allocatable    :: kbz     !< k-Point vectors in the first BZ
  real(double), dimension(:), allocatable      :: wkbz    !< Weights of k-Points in kbz

contains
  !-- Generating integration points in the Brillouin Zone --

  !-----------------------------------------------------------------------------
  !> @author
  !> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut
  !
  !> @DESCRIPTION
  !> Wrapper function for k-point generation.
  !> This subroutine will call either the subroutine for the full 3D Brillouin Zone
  !> or if a plane direction is given the subroutine for the surface Brillouin
  !> Zone in the plane direction.
  !
  ! REVISION HISTORY:
  ! 14 March 2017 - Current Revision
  !> @param[in] pln Optional Parameter - Plane direction vector
  !-----------------------------------------------------------------------------
  subroutine generate_kpoints(pln)
    use mod_f90_kind, only: double
    use mod_mpi_pars, only: myrank
    implicit none
    real(double), optional, intent(in) :: pln(3)

    if(present(pln)) then
       call generate_kpoints_2d(pln)
    else
       call generate_kpoints_3d()
    end if

    if((myrank==0))  write(outputunit,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints

    return
  end subroutine generate_kpoints

  !-----------------------------------------------------------------------------
  !> @author
  !> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut
  !
  !> @description
  !> This subroutine will generate the surface Brillouin Zone in
  !> the given plane direction.
  !
  ! REVISION HISTORY:
  ! 14 March 2017 - Current Revision
  !> @param[in] pln Optional Parameter - Plane direction vector
  !-----------------------------------------------------------------------------
  subroutine generate_kpoints_2d(pln)
    use mod_tools,      only: cross, is_parallel
    use mod_constants,  only: tpi
    use mod_lattice,    only: r0, n0
    use mod_parameters, only: nkpt
    implicit none
    real(kind=8), optional, intent(in) :: pln(3)
    real(kind=8) :: b1(3), b2(3), b3(3)
    real(kind=8) :: BZ(4,3), smallest_dist, distance, diff(3), ini_smallest_dist
    real(kind=8), allocatable :: inikbz(:,:), iniwkbz(:)
    real(kind=8) :: extrakbz(nkpt*10,3), extrawkbz(nkpt*10)
    integer :: i,l,j,m,k, smallest_indice, numextrakbz
    integer :: nkpt_perdim

    nkpt_perdim = ceiling(sqrt(dble(nkpt)) / 6.d0)
    nkpt_perdim = nkpt_perdim * 6
    nkpoints    = nkpt_perdim**2
    allocate( iniwkbz(nkpoints), inikbz(nkpoints, 3))

    b2 = cross(r0(1,:),pln)
    b1 = 0.d0

    do i = 2, n0
       if( .not. is_parallel(r0(1,:),r0(i,:))) then
          b1 = cross(r0(i,:), pln)
          exit
       end if
    end do

    if(0.d0 == dot_product(b1,b1)) then
       if(myrank==0) write(outputunit,"('[generate_kpoints_2d] No non-collinear in-plane vector found')")
       call MPI_Finalize(ierr)
       stop
    end if

    b1 = b1  * tpi / (dot_product( b1, r0(1,:)))
    b2 = b2  * tpi / (dot_product( b2, r0(i,:)))
    b3 = pln * tpi / (dot_product(pln, pln    ))
    BZ(1,:)=0.d0
    BZ(2,:)=b1
    BZ(3,:)=b2
    BZ(4,:)=b1+b2
    !call abort
    !Generate k-points in the paralelogram determined by b1 and b2
    iniwkbz=1.d0
    m=0
    do l=1, nkpt_perdim; do j=1, nkpt_perdim !Run over the entire BZ
       m=m+1
       inikbz(m,:)= ( dble(l-1)*b1 + dble(j-1)*b2 )/ dble(nkpt_perdim)
    end do; end do

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist=10.d0*sqrt(DOT_product(b1+b2,b1+b2))
    numextrakbz=0
    !Run over all the kpoints generated initially.
    do l=1, nkpoints
       smallest_dist=ini_smallest_dist
       m=0
       !Checks to each of the 4 BZ's the kpoint belongs by checking
       ! to each BZ it's closer.
       do j=1, 4
          diff=inikbz(l,:)-BZ(j,:)
          distance=sqrt(DOT_product(diff,diff))
          if(distance<smallest_dist) then
             smallest_dist=distance
             smallest_indice=j
          end if
       end do
       !Checks if the kpoint is in the border between two or more
       ! BZ's. If yes, create a clone of it to translate later into
       ! the 1st BZ.
       do j=1, 4
          diff=inikbz(l,:)-BZ(j,:)
          distance=sqrt(DOT_product(diff,diff))
          if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. &
               j/=smallest_indice ) then
             m=m+1
             numextrakbz=numextrakbz+1
             extrakbz(numextrakbz,:)=inikbz(l,:)-BZ(j,:)
          end if
       end do
       if(m/=0) then
          !The weight of the kpoint in the border is shared with
          ! its clones.
          iniwkbz(l)=1.d0/dble(m+1)
          extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/dble(m+1)
       end if
       !Translate the kpoint to the 1st BZ
       inikbz(l,:)=inikbz(l,:)-BZ(smallest_indice,:)
    end do

    allocate( wkbz(nkpoints+numextrakbz),kbz(nkpoints+numextrakbz,3) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    kbz (1:nkpoints,:)=inikbz (:,:)
    wkbz(1:nkpoints)  =iniwkbz(:)
    if(numextrakbz/=0) then
       kbz(nkpoints+1:nkpoints+numextrakbz,:)=extrakbz(1:numextrakbz,:)
       wkbz(nkpoints+1:nkpoints+numextrakbz) =extrawkbz(1:numextrakbz)
    end if
    wkbz=wkbz/dble(nkpoints)
    nkpoints=nkpoints+numextrakbz
    nkpoints=nkpoints

    !PRINT *, "firstBZ = ", firstBZ
    !PRINT "(a,i0,a,i0)", "nkpt=",nkpoints,"  nkpt_perdim=",nkpt_perdim
    return
  end subroutine generate_kpoints_2d

  !-----------------------------------------------------------------------------
  !> @author
  !> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut
  !
  !> @DESCRIPTION
  !> This subroutine will generate the full 3D Brillouin Zone.
  !
  ! REVISION HISTORY:
  ! 14 March 2017 - Current Revision
  !> @param[in] pln Optional Parameter - Plane direction vector
  !-----------------------------------------------------------------------------
  subroutine generate_kpoints_3d()
    use mod_tools,      only: cross
    use mod_constants,  only: tpi
    use mod_parameters, only: nkpt, a1, a2, a3
    implicit none
    real(kind=8) :: b1(3), b2(3), b3(3)
    real(kind=8) :: vol
    real(KIND=8) :: BZ(8,3), smallest_dist, distance, diff(3), ini_smallest_dist
    real(KIND=8), allocatable :: inikbz(:,:), iniwkbz(:)
    real(KIND=8) :: extrakbz(nkpt*10,3), extrawkbz(nkpt*10)
    integer :: l,j,m, k, smallest_indice, numextrakbz
    integer :: nkpt_perdim !n. of k point per dimension

    nkpt_perdim=ceiling((dble(nkpt))**(1.d0/3.d0)/6.d0)
    nkpt_perdim=nkpt_perdim*6
    nkpoints=nkpt_perdim**3
    allocate( iniwkbz(nkpoints), inikbz(nkpoints,3) )

    vol = tpi / dot_product(a1, cross(a2,a3))
    b1 = vol * cross(a2,a3)
    b2 = vol * cross(a3,a1)
    b3 = vol * cross(a1,a2)

    BZ(1,:)=0.d0
    BZ(2,:)=b1
    BZ(3,:)=b2
    BZ(4,:)=b1+b2
    BZ(5,:)=b3
    BZ(6,:)=b1+b3
    BZ(7,:)=b2+b3
    BZ(8,:)=b1+b2+b3

    !Generate k-points in the paralelogram determined by b1 and b2
    iniwkbz=1.d0
    m=0
    do l=1, nkpt_perdim; do j=1, nkpt_perdim; do k=1, nkpt_perdim !Run over the entire BZ
       m=m+1
       inikbz(m,:)= ( dble(l-1)*b1 + dble(j-1)*b2 + dble(k-1)*b3)/ dble(nkpt_perdim)
    end do; end do; end do

    !Translate the k-points to the 1st BZ.
    !10*|b1+b2|, bigger than the distance of any genarated kpoint
    ini_smallest_dist=10.d0*sqrt(DOT_product(b1+b2,b1+b2))
    numextrakbz=0
    !Run over all the kpoints generated initially.
    do l=1, nkpoints
       smallest_dist=ini_smallest_dist
       m=0
       !Checks to each of the 4 BZ's the kpoint belongs by checking
       ! to each BZ it's closer.
       do j=1, 8
          diff=inikbz(l,:)-BZ(j,:)
          distance=sqrt(DOT_product(diff,diff))
          if(distance<smallest_dist) then
             smallest_dist=distance
             smallest_indice=j
          end if
       end do
       !Checks if the kpoint is in the border between two or more
       ! BZ's. If yes, create a clone of it to translate later into
       ! the 1st BZ.
       do j=1, 8
          diff=inikbz(l,:)-BZ(j,:)
          distance=sqrt(DOT_product(diff,diff))
          if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. &
               j/=smallest_indice ) then
             m=m+1
             numextrakbz=numextrakbz+1
             extrakbz(numextrakbz,:)=inikbz(l,:)-BZ(j,:)
          end if
       end do
       if(m/=0) then
          !The weight of the kpoint in the border is shared with
          ! its clones.
          iniwkbz(l)=1.d0/dble(m+1)
          extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/dble(m+1)
       end if
       !Translate the kpoint to the 1st BZ
       inikbz(l,:)=inikbz(l,:)-BZ(smallest_indice,:)
    end do

    allocate( wkbz(nkpoints+numextrakbz),kbz(nkpoints+numextrakbz,3) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    kbz (1:nkpoints,:)=inikbz (:,:)
    wkbz(1:nkpoints)  =iniwkbz(:)
    if(numextrakbz/=0) then
       kbz(nkpoints+1:nkpoints+numextrakbz,:)=extrakbz(1:numextrakbz,:)
       wkbz(nkpoints+1:nkpoints+numextrakbz) =extrawkbz(1:numextrakbz)
    end if
    wkbz=wkbz/dble(nkpoints)
    nkpoints=nkpoints+numextrakbz

    !PRINT *, "firstBZ = ", firstBZ
    !PRINT "(a,i0,a,i0)", "nkpt=",nkpoints,"  nkpt_perdim=",nkpt_perdim
    return
  end subroutine generate_kpoints_3d

  !-----------------------------------------------------------------------------
  !> @author
  !> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut
  !
  !> @DESCRIPTION
  !> Writing the k-points to file "kpoints3d".
  !
  ! REVISION HISTORY:
  ! 14 March 2017 - Current Revision
  !-----------------------------------------------------------------------------
  subroutine write_kpoints_to_file()
    implicit none
    integer :: i

    !open (unit=2222, file='kpoints2d',status='unknown')
    !write(unit=2222,fmt="(a)") ' #      kx            ky            wk'
    open (unit=3333, file='kpoints3d',status='unknown')
    write(unit=3333,fmt="(a)") ' #      kx            ky            kz            wk'
    do i=1,nkpoints
       !write(unit=2222,fmt="(3(f12.9,2x))") kbz2d(i,1),kbz2d(i,2),wkbz(i)
       write(unit=3333,fmt="(4(f12.9,2x))") kbz(i,1),kbz(i,2),kbz(i,3),wkbz(i)
    end do
    close(2222)
    close(3333)

    return
  end subroutine write_kpoints_to_file

end module mod_generate_kpoints
