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
  use mod_f90_kind
  use mod_parameters, only: a0, ncp, outputunit
  use mod_constants
  use mod_mpi_pars
  implicit none

! wave vector integration variables
  integer :: nkpoints                           !> Number of generated k-points
  real(double),allocatable      :: kbz(:,:)     !< k-Point vectors in the first BZ
  real(double),allocatable      :: wkbz(:)      !< Weights of k-Points in kbz
  !real(double),allocatable      :: kbz2d(:,:)   !< Legacy variable, 2D BZ

contains
!   !-- Generating integration points in the Brillouin Zone --
!   subroutine generate_kpoints_bcc110()
!     use mod_mpi_pars
!     implicit none
!     integer         :: i,j
!     integer         :: m0,m1,imax,ki,jmax,kj,iz
!     real(double)    :: beta,beta2,fact
!     real(double)    :: akxp(4),akzp(4),akxt,akzt,akx,aky,akz
!     !   Generation of the Cunningham points for a centred rectangular
!     !   lattice. The units are chosen with the lattice parameter a0.
!
!     allocate( kbz(4**(ncp+1),3),wkbz(4**(ncp+1)),kbz2d(4**(ncp+1),2) )
!
! !--- k integration parameters ------
!     beta  = 0.5d0*sq2
!     beta2 = beta*beta
!     fact  = sq2/4.d0
!     m0 = 2**ncp
!     m1 = m0*2
!     imax = m0 - 0.5d0
! !-----------------------------------
!     ! Generating k points
!     nkpoints   = 0
!     do ki = 0,imax
!       i      = 2*ki+1
!       akxp(1) = tpi*fact*dble(i)/(a0*dble(m0))
!       jmax    = nint(beta2*(-2*ki+m0-1))+m0-0.5d0
!       do kj=0,jmax
!         j = 2*kj+1
!         akzp(1) = tpi*fact*dble(j)/(a0*dble(m1)*beta)
!   !     Generating k_// inside the full BZ
!   !     2nd quadrant
!         akxp(2) =-akxp(1)
!         akzp(2) = akzp(1)
!   !     3rd quadrant
!         akxp(3) =-akxp(1)
!         akzp(3) =-akzp(1)
!   !     4th quadrant
!         akxp(4) = akxp(1)
!         akzp(4) =-akzp(1)
!
!         do iz=1,4
!           nkpoints = nkpoints + 1
!
!           akxt = akxp(iz)
!           akzt = akzp(iz)
!
!           ! Storing the kpoints in the 2DBZ
!           kbz2d(nkpoints,1) = akxt
!           kbz2d(nkpoints,2) = akzt
!
!           ! Transformation to cartesian axis
!           akx = akxt*beta
!           aky = akx
!           akz = akzt
!
!           ! Storing the kpoints in the BZ transformed to the cartesian system
!           kbz(nkpoints,1) = akx
!           kbz(nkpoints,2) = aky
!           kbz(nkpoints,3) = akz
!         end do ! iz
!       end do ! kj
!     end do ! ki
!
!     wkbz(:) = 1.d0/dble(nkpoints)
!
!     if((myrank.eq.0))  write(outputunit,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints
!
!     return
!   end subroutine generate_kpoints_bcc110
!
!   subroutine generate_kpoints_fcc100()
!     use mod_mpi_pars
!     implicit none
!     integer         :: iz
!     integer         :: ki,kj,nkmax,icount
!     real(double)    :: auxk,auxw
!     real(double)    :: akp(8,2)
! !   generation of Cunningham points for a square lattice
! !   the units are chosen using the variable for the lattice parameter a0.
!
!     if(ncp.le.0) then
!       if(myrank.eq.0) write(outputunit,"('[mod_generate_kpoints] ncp must be greater than 0')")
!       call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
!     end if
!     nkpoints = 8*(2**(ncp-1))*(1+2**ncp)
!     allocate( kbz(nkpoints,3),wkbz(nkpoints),kbz2d(nkpoints,2) )
!     nkmax = 2**ncp
!     auxk  = dble(2*nkmax)
!     auxw  = 8.d0*dble(2**(ncp+ncp-1))
!
!     icount = 0
!     do ki=1,nkmax
!       akp(1,1)=pi*(dble(2*ki-1))/(auxk*a0)
!       do kj=1,ki
!         akp(1,2)=pi*(dble(2*kj-1))/(auxk*a0)
!
!         ! generating k_// inside the full bz
!         ! 1st quadrant
!         akp(2,1) = akp(1,2)
!         akp(2,2) = akp(1,1)
!         ! 2nd quadrant
!         akp(3,1) =-akp(2,1)
!         akp(3,2) = akp(2,2)
!         akp(4,1) =-akp(3,2)
!         akp(4,2) =-akp(3,1)
!         ! 3rd quadrant
!         akp(5,1) = akp(4,1)
!         akp(5,2) =-akp(4,2)
!         akp(6,1) = akp(5,2)
!         akp(6,2) = akp(5,1)
!         ! 4th quadrant
!         akp(7,1) =-akp(6,1)
!         akp(7,2) = akp(6,2)
!         akp(8,1) =-akp(7,2)
!         akp(8,2) =-akp(7,1)
!
!         do iz=1,8
!           icount = icount + 1
!
!           kbz2d(icount,1) = akp(iz,1)
!           kbz2d(icount,2) = akp(iz,2)
!
!           kbz(icount,1) = akp(iz,1)+akp(iz,2)
!           kbz(icount,2) =-akp(iz,1)+akp(iz,2)
!           kbz(icount,3) = 0.d0
!
!           ! Weight of k point
!           if(ki.eq.kj) then ! Half of the weight if the point is in the diagonal
!             wkbz(icount)=0.5d0/auxw
!           else
!             wkbz(icount)=1.d0/auxw
!           end if
!
!         end do
!       end do
!     end do
!
!     if(icount.ne.nkpoints) then
!        if(myrank.eq.0) write(outputunit,"('[mod_generate_kpoints] Incorrect number of points: icount = ',i0,', nkpoints = ',i0)") icount,nkpoints
!        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
!     end if
!     if((myrank.eq.0))  write(outputunit,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints
!
!     return
!   end subroutine generate_kpoints_fcc100
!
!   subroutine generate_kpoints_fcc111()
!     ! Generation of the k-points points for an hexagonal lattice.
!     ! The generated k-points units are chosen using a0.
!     use mod_mpi_pars
!     use mod_constants, only: pi,sq3
!     implicit none
!     integer         :: iz
!     integer         :: ki,kj
!     real(kind=8)    :: a1(2),a2(2),b1(2),b2(2),c1(2),c2(2)
!     real(kind=8)    :: weight,dnkpoints
!     real(kind=8)    :: kv(2),angle,rm(2,2)
!     real(double), dimension(3,3) :: rm2,rm3
!     real(kind=8), dimension(100000)   :: wk
!     real(kind=8), dimension(100000,2) :: kvm
!     real(double), dimension(3) :: qvl,k2anx
!
! !   Lattice vetors in real space
!     a1 = [0.5d0, sqrt(3.d0)/2.d0]*a0
!     a2 = [0.5d0,-sqrt(3.d0)/2.d0]*a0
!
! !   Reciprocal lattice determination
!     b1 =         [a2(2),-a2(1)]
!     b1 = 2.d0*pi*[a2(2),-a2(1)]/dot_product(a1,b1)
!     b2 =         [a1(2),-a1(1)]
!     b2 = 2.d0*pi*[a1(2),-a1(1)]/dot_product(a2,b2)
!
!     angle = -pi/6.d0
!     rm(1,1) = cos(angle)
!     rm(2,1) = sin(angle)
!     rm(1,2) = -rm(2,1)
!     rm(2,2) = rm(1,1)
!
!     c1 = (b1+b2)/3.d0
!     c2 = b1-c1
!     c1 = matmul(rm,c1)
!     c2 = matmul(rm,c2)
!
!     kvm = 0.d0
!
!     ! Gamma point
!     nkpoints = 1
!     kvm(nkpoints,:) = [ 0.d0,0.d0 ]
!     wk(nkpoints) = 1.d0
!     dnkpoints = wk(nkpoints)
!
!     do ki=1,ncp
!       do kj=1,ncp
!         weight = 1.d0
!         kv = ( ki*c1 + (kj-1.d0)*c2 ) / DBLE(ncp)
!
!         if(kv(1).gt.2.d0*pi/sq3) cycle
!         if(kv(1).eq.2.d0*pi/sq3) weight = 0.5d0
!
!         if((ki.eq.ncp).and.(kj.eq.1)) weight = 1.d0/3.d0
!
!         dnkpoints = dnkpoints + weight
!         nkpoints = nkpoints + 1
!         kvm(nkpoints,:) = kv
!         wk(nkpoints) = weight
!
!         do iz=1,5
!           angle = dble(iz)*pi/3.d0
!           rm(1,1) = cos(angle)
!           rm(1,2) = -sin(angle)
!           rm(2,1) = sin(angle)
!           rm(2,2) = cos(angle)
!
!           dnkpoints = dnkpoints + weight
!           nkpoints = nkpoints + 1
!           kvm(nkpoints,:) = matmul(rm,kv)
!           wk(nkpoints) = weight
!         end do
!       end do
!     end do
!
!     allocate( kbz(nkpoints,3),wkbz(nkpoints),kbz2d(nkpoints,2) )
!
!     ! Rotation matrices to transform 2D BZ to 3D BZ
!     rm2 = 0.d0
!     rm2(1,1) = 1.d0/sq3
!     rm2(1,3) = sq2/sq3
!     rm2(3,1) = -sq2/sq3
!     rm2(2,2) = 1.d0
!     rm2(3,3) = 1.d0/sq3
!
!     rm3 = 0.d0
!     rm3(1,1) = 1.d0/sq2
!     rm3(1,2) = -1.d0/sq2
!     rm3(2,1) = 1.d0/sq2
!     rm3(2,2) = 1.d0/sq2
!     rm3(3,3) = 1.d0
!
!     weight = 0.d0
!     do iz = 1,nkpoints
!       kbz2d(iz,:) = kvm(iz,:)
!       k2anx(1) = kvm(iz,1)
!       k2anx(2) = kvm(iz,2)
!       k2anx(3) = 0.d0
!       qvl = matmul(rm2,k2anx)
!       kbz(iz,:) = matmul(rm3,qvl)
!       wkbz(iz) = wk(iz)/dnkpoints
!       weight = weight + wkbz(iz)
!     end do
!
!     if((myrank.eq.0))  write(outputunit,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints
!
! !     write(outputunit,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints
! !     write(outputunit,"('[mod_generate_kpoints] ',es9.2,' total weight.')") weight
! !     write(outputunit,"('[mod_generate_kpoints] ',f15.2,' dnkpoints.')") dnkpoints
!
!     return
!   end subroutine generate_kpoints_fcc111

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
    use mod_mpi_pars
    implicit none
    real(double), optional, intent(in) :: pln(3)

    if(present(pln)) then
      call generate_kpoints_2d(pln)
    else
      call generate_kpoints_3d()
    end if

    if((myrank.eq.0))  write(outputunit,"('[mod_generate_kpoints] ',i0,' k-points generated.')") nkpoints

    return
  end subroutine

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
    use mod_tools
    use mod_constants
    use mod_lattice
    use mod_parameters
    implicit none
    integer :: eff_nkpt
    logical :: firstBZ = .true.
    real(kind=8), optional, intent(in) :: pln(3)
    real(kind=8) :: b1(3), b2(3), b3(3)
    real(kind=8) :: vol
    real(kind=8) :: BZ(4,3), smallest_dist, distance, diff(3), ini_smallest_dist
    real(kind=8), allocatable :: inikbz(:,:), iniwkbz(:)
    real(kind=8) :: extrakbz(nkpt*10,3), extrawkbz(nkpt*10)
    integer :: i,l,j,m,k, smallest_indice, numextrakbz
    integer :: nkpt_perdim

    nkpt_perdim = ceiling(sqrt(dble(nkpt)) / 6.d0)
    nkpt_perdim = nkpt_perdim * 6
    eff_nkpt    = nkpt_perdim**2
    allocate( iniwkbz(eff_nkpt), inikbz(eff_nkpt, 3))

    b2 = cross(r0(1,:),pln)
    b1 = 0.d0

    do i = 2, n0
      if( .not. is_parallel(r0(1,:),r0(i,:))) then
        b1 = cross(r0(i,:), pln)
        exit
      end if
    end do

    if(0.d0 == dot_product(b1,b1)) then
      if(myrank.eq.0) write(outputunit,"('[generate_kpoints_2d] No non-collinear in-plane vector found')")
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
    DO l=1, nkpt_perdim; DO j=1, nkpt_perdim !Run over the entire BZ
       m=m+1
       inikbz(m,:)= ( DBLE(l-1)*b1 + DBLE(j-1)*b2 )/ DBLE(nkpt_perdim)
    END DO; END DO

    !Translate the k-points to the 1st BZ.
    IF(firstBZ) THEN
       !10*|b1+b2|, bigger than the distance of any genarated kpoint
       ini_smallest_dist=10.d0*SQRT(DOT_PRODUCT(b1+b2,b1+b2))
       numextrakbz=0
       !Run over all the kpoints generated initially.
       DO l=1, eff_nkpt
          smallest_dist=ini_smallest_dist
          m=0
          !Checks to each of the 4 BZ's the kpoint belongs by checking
          ! to each BZ it's closer.
          DO j=1, 4
             diff=inikbz(l,:)-BZ(j,:)
             distance=SQRT(DOT_PRODUCT(diff,diff))
             IF(distance.LT.smallest_dist) THEN
                smallest_dist=distance
                smallest_indice=j
             END IF
          END DO
          !Checks if the kpoint is in the border between two or more
          ! BZ's. If yes, create a clone of it to translate later into
          ! the 1st BZ.
          DO j=1, 4
             diff=inikbz(l,:)-BZ(j,:)
             distance=SQRT(DOT_PRODUCT(diff,diff))
             IF( ( ABS(distance-smallest_dist) .LT. 1.d-12 ) .AND. &
                                           j.NE.smallest_indice ) THEN
                m=m+1
                numextrakbz=numextrakbz+1
                extrakbz(numextrakbz,:)=inikbz(l,:)-BZ(j,:)
             END IF
          END DO
          IF(m.NE.0) THEN
             !The weight of the kpoint in the border is shared with
             ! its clones.
             iniwkbz(l)=1.d0/DBLE(m+1)
             extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/DBLE(m+1)
          END IF
          !Translate the kpoint to the 1st BZ
          inikbz(l,:)=inikbz(l,:)-BZ(smallest_indice,:)
       END DO
    END IF !IF(firstBZ)

    ALLOCATE( wkbz(eff_nkpt+numextrakbz),kbz(eff_nkpt+numextrakbz,3) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    kbz (1:eff_nkpt,:)=inikbz (:,:)
    wkbz(1:eff_nkpt)  =iniwkbz(:)
    IF(numextrakbz.NE.0) THEN
       kbz(eff_nkpt+1:eff_nkpt+numextrakbz,:)=extrakbz(1:numextrakbz,:)
       wkbz(eff_nkpt+1:eff_nkpt+numextrakbz) =extrawkbz(1:numextrakbz)
    END IF
    wkbz=wkbz/DBLE(eff_nkpt)
    eff_nkpt=eff_nkpt+numextrakbz
    nkpoints=eff_nkpt

    !PRINT *, "firstBZ = ", firstBZ
    !PRINT "(a,i0,a,i0)", "nkpt=",eff_nkpt,"  nkpt_perdim=",nkpt_perdim
    return
  end subroutine

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
    use mod_tools
    use mod_constants
    use mod_parameters
    implicit none
    integer :: eff_nkpt
    logical :: firstBZ = .true.
    real(kind=8) :: b1(3), b2(3), b3(3)
    real(kind=8) :: vol
    REAL(KIND=8) :: BZ(8,3), smallest_dist, distance, diff(3), ini_smallest_dist
    REAL(KIND=8), ALLOCATABLE :: inikbz(:,:), iniwkbz(:)
    REAL(KIND=8) :: extrakbz(nkpt*10,3), extrawkbz(nkpt*10)
    INTEGER :: l,j,m, k, smallest_indice, numextrakbz
    INTEGER :: nkpt_perdim !n. of k point per dimension

    nkpt_perdim=ceiling((dble(nkpt))**(1.d0/3.d0)/6.d0)
    nkpt_perdim=nkpt_perdim*6
    eff_nkpt=nkpt_perdim**3
    ALLOCATE( iniwkbz(eff_nkpt), inikbz(eff_nkpt,3) )

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
    DO l=1, nkpt_perdim; DO j=1, nkpt_perdim; DO k=1, nkpt_perdim !Run over the entire BZ
       m=m+1
       inikbz(m,:)= ( DBLE(l-1)*b1 + DBLE(j-1)*b2 + DBLE(k-1)*b3)/ DBLE(nkpt_perdim)
    END DO; END DO; END DO

    !Translate the k-points to the 1st BZ.
    IF(firstBZ) THEN
       !10*|b1+b2|, bigger than the distance of any genarated kpoint
       ini_smallest_dist=10.d0*SQRT(DOT_PRODUCT(b1+b2,b1+b2))
       numextrakbz=0
       !Run over all the kpoints generated initially.
       DO l=1, eff_nkpt
          smallest_dist=ini_smallest_dist
          m=0
          !Checks to each of the 4 BZ's the kpoint belongs by checking
          ! to each BZ it's closer.
          DO j=1, 8
             diff=inikbz(l,:)-BZ(j,:)
             distance=SQRT(DOT_PRODUCT(diff,diff))
             IF(distance.LT.smallest_dist) THEN
                smallest_dist=distance
                smallest_indice=j
             END IF
          END DO
          !Checks if the kpoint is in the border between two or more
          ! BZ's. If yes, create a clone of it to translate later into
          ! the 1st BZ.
          DO j=1, 8
             diff=inikbz(l,:)-BZ(j,:)
             distance=SQRT(DOT_PRODUCT(diff,diff))
             IF( ( ABS(distance-smallest_dist) .LT. 1.d-12 ) .AND. &
                                           j.NE.smallest_indice ) THEN
                m=m+1
                numextrakbz=numextrakbz+1
                extrakbz(numextrakbz,:)=inikbz(l,:)-BZ(j,:)
             END IF
          END DO
          IF(m.NE.0) THEN
             !The weight of the kpoint in the border is shared with
             ! its clones.
             iniwkbz(l)=1.d0/DBLE(m+1)
             extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/DBLE(m+1)
          END IF
          !Translate the kpoint to the 1st BZ
          inikbz(l,:)=inikbz(l,:)-BZ(smallest_indice,:)
       END DO
    END IF !IF(firstBZ)

    ALLOCATE( wkbz(eff_nkpt+numextrakbz),kbz(eff_nkpt+numextrakbz,3) )
    !The final array of kpoints will be the initial one by the clones
    ! of the ones in the border between BZ's.
    kbz (1:eff_nkpt,:)=inikbz (:,:)
    wkbz(1:eff_nkpt)  =iniwkbz(:)
    IF(numextrakbz.NE.0) THEN
       kbz(eff_nkpt+1:eff_nkpt+numextrakbz,:)=extrakbz(1:numextrakbz,:)
       wkbz(eff_nkpt+1:eff_nkpt+numextrakbz) =extrawkbz(1:numextrakbz)
    END IF
    wkbz=wkbz/DBLE(eff_nkpt)
    eff_nkpt=eff_nkpt+numextrakbz
    nkpoints=eff_nkpt

    !PRINT *, "firstBZ = ", firstBZ
    !PRINT "(a,i0,a,i0)", "nkpt=",eff_nkpt,"  nkpt_perdim=",nkpt_perdim
    return
  end subroutine
end module mod_generate_kpoints
