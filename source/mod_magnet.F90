module mod_magnet
  use mod_f90_kind, only: double
  use mod_parameters, only: dmax
  implicit none

  logical :: lfield !< Turn on/off static magnetic field, option to give in magnetic field in tesla

  integer                     :: iter                                   ! self-consistency iteration
  real(double),allocatable    :: mx(:),my(:),mz(:),mvec_cartesian(:,:),hdel(:)    ! Magnetization and exchange split delta/2
  real(double),allocatable    :: lxm(:),lym(:),lzm(:)                   ! Orbital angular momentum in global frame of reference
  real(double),allocatable    :: lxpm(:),lypm(:),lzpm(:)                ! Orbital angular momentum in local frame of reference
  complex(double),allocatable :: mp(:),hdelp(:)
  complex(double),allocatable :: mm(:),hdelm(:)
  real(double),allocatable    :: mabs(:),mtheta(:),mphi(:),mvec_spherical(:,:)
  real(double),allocatable    :: labs(:),ltheta(:),lphi(:)
  real(double),allocatable    :: lpabs(:),lptheta(:),lpphi(:)
  real(double),allocatable    :: eps1(:)                                ! Center of the bands for each l - eps(Npl)
  real(double), dimension(:), allocatable :: hhwx, hhwy, hhwz           ! Static magnetic fields in each direction
  complex(double), dimension(:,:,:), allocatable :: lb, sb
  !! Zeeman matrices
  complex(double), dimension(:,:,:), allocatable :: lxp, lyp, lzp
  !! Site dependent Angular momentum matrices in local frame
  complex(double), dimension(:,:), allocatable :: lx, ly, lz
  !! Angular momentum matrices in global frame

  !========================================================================================!
  ! Values of magnetic field in cartesian or spherical coordinates
  !integer, parameter :: dmax = 20
  character(len=9), dimension(7) :: dcfield = ["hwa      ","hwt      ","hwp      ","hwahwt   ","hwahwp   ","hwthwp   ","hwahwthwp"]
  character(len=60) :: dc_header
  character(len=60), allocatable :: dc_fields(:),dcprefix(:)
  real(double), allocatable :: hw_list(:,:)
  integer :: dcfield_dependence = 0, dc_count = 0

  real(double) :: hwx = 0.d0, hwy = 0.d0, hwz = 0.d0, tesla = 1.d0

  integer :: hwa_npts = 0, hwa_npt1 = 1
  integer :: hwt_npts = 0, hwt_npt1 = 1
  integer :: hwp_npts = 0, hwp_npt1 = 1
  real(double) :: hwa, hwa_i = 0.d0, hwa_f = 0.d0, hwa_s
  real(double) :: hwt, hwt_i = 0.d0, hwt_f = 0.d0, hwt_s
  real(double) :: hwp, hwp_i = 0.d0, hwp_f = 0.d0, hwp_s

  ! Layer-resolved scale of magnetic field (including empty spheres)
  logical :: lhwscale        = .false.
  real(double) :: hwscale(dmax)   = 1.d0
  logical :: lhwrotate       = .false.
  real(double) :: hwtrotate(dmax) = 0.d0, hwprotate(dmax) = 0.d0

  integer :: skip_steps_hw = 0
  !! How many iterations are to be skipped from the beginning
  integer :: hw_count
  !! Current iteration of magnetic field loop.
  integer :: total_hw_npt1
  !! Total number of magnetic field iterations.

contains

  ! This subroutine sets up external magnetic fields and related loop
  !subroutine prepare_field()
  subroutine setMagneticLoopPoints()
    implicit none
    integer       :: i, j, k, l

    if(lfield) then
      if (hwa_npts==0) then
        hwa_f = hwa_i
        hwa_npts = 1
      end if
      if((abs(hwa_i) + abs(hwa_f)) > 1.d-8) then
        if (hwt_npts==0) then
          hwt_f = hwt_i
          hwt_npts = 1
        end if
        if (hwp_npts==0) then
          hwp_f = hwp_i
          hwp_npts = 1
        end if
        ! External field angular loops steps
        hwa_s = (hwa_f - hwa_i)/hwa_npts
        if(abs(hwa_s) <= 1.d-10) hwa_npt1 = 1
        hwt_s = (hwt_f - hwt_i)/hwt_npts
        if(abs(hwt_s) <= 1.d-10) hwt_npt1 = 1
        hwp_s = (hwp_f - hwp_i)/hwp_npts
        if(abs(hwp_s) <= 1.d-10) hwp_npt1 = 1
      else ! hwa_i and hwa_f = 0
        hwa_i   = sqrt(hwx**2+hwy**2+hwz**2)
        hwa_f   = hwa_i
        if(abs(hwa_i)<1.d-8) then
          lfield = .false.
        else
          hwt_i    = acos(hwz/hwa_i)
          hwt_s    = 0.d0
          hwt_npt1 = 1
          hwp_i    = atan2(hwy,hwx)
          hwp_s    = 0.d0
          hwp_npt1 = 1
        end if
      end if
    else ! lfield
      hwa_i    = 0.d0
      hwa_f    = 0.d0
      hwa_s    = 0.d0
      hwa_npt1 = 1
      hwt_i    = 0.d0
      hwt_f    = 0.d0
      hwt_s    = 0.d0
      hwt_npt1 = 1
      hwp_i    = 0.d0
      hwp_f    = 0.d0
      hwp_s    = 0.d0
      hwp_npt1 = 1
    end if

    ! Total number of points in the loops
    total_hw_npt1 = hwa_npt1*hwt_npt1*hwp_npt1
    if(total_hw_npt1==1) skip_steps_hw = 0
    allocate(hw_list(total_hw_npt1,3))

    ! Creating list of magnetic fields (in spherical coordinates)
    i = 1
    do j = 1, hwa_npt1
      hwa = hwa_i + (j-1)*hwa_s
      do k=1, hwt_npt1
        hwt = hwt_i + (k-1)*hwt_s
        do l = 1, hwp_npt1
          hwp = hwp_i + (l-1)*hwp_s
          hw_list(i,:) = [ hwa , hwt , hwp ]
          i = i+1
        end do
      end do
    end do

    return
  end subroutine setMagneticLoopPoints

  subroutine initMagneticField(nAtoms)
    use mod_f90_kind, only: double
    use mod_constants, only: cZero, pi
    use mod_parameters, only: outputfile_loop
    use mod_mpi_pars, only: abortProgram, myrank_row_hw
    use mod_SOC, only: SOC, llinearsoc
    implicit none

    integer, intent(in) :: nAtoms
    integer :: AllocateStatus
    integer :: i

    allocate( hhwx(nAtoms),hhwy(nAtoms),hhwz(nAtoms), STAT = AllocateStatus )
    if (AllocateStatus /= 0) call abortProgram("[main] Not enough memory for: hhwx,hhwy,hhwz,sb,lb")

    !---------------------- Turning off field for hwa=0 --------------------
    if(abs(hw_list(hw_count,1)) < 1.d-8) then
      lfield = .false.
      if((llinearsoc) .or. (.not.SOC) .and. (myrank_row_hw==0)) write(outputfile_loop,"('[main] WARNING: No external magnetic field is applied and SOC is off/linear order: Goldstone mode is present!')")
    else
      lfield = .true.
    end if
    sb   = cZero
    lb   = cZero
    hhwx = 0.d0
    hhwy = 0.d0
    hhwz = 0.d0

    !--------------------- Defining the magnetic fields -------------------- TODO
    if(lfield) then
      ! Variables of the hamiltonian
      ! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
      ! to take into account the fact that we are considering negative
      ! external fields to get the peak at positive energies
      do i = 1, nAtoms
        hhwx(i) = -0.5d0*hwscale(i)*hw_list(hw_count,1)*sin((hw_list(hw_count,2) + hwtrotate(i))*pi)*cos((hw_list(hw_count,3) + hwprotate(i))*pi)*tesla
        hhwy(i) = -0.5d0*hwscale(i)*hw_list(hw_count,1)*sin((hw_list(hw_count,2) + hwtrotate(i))*pi)*sin((hw_list(hw_count,3) + hwprotate(i))*pi)*tesla
        hhwz(i) = -0.5d0*hwscale(i)*hw_list(hw_count,1)*cos((hw_list(hw_count,2) + hwtrotate(i))*pi)*tesla
        if(abs(hhwx(i))<1.d-8) hhwx(i) = 0.d0
        if(abs(hhwy(i))<1.d-8) hhwy(i) = 0.d0
        if(abs(hhwz(i))<1.d-8) hhwz(i) = 0.d0
      end do
      ! Testing if hwscale is used
      lhwscale = any(abs(hwscale(1:nAtoms)-1.d0) > 1.d-8)
      ! Testing if hwrotate is used
      lhwrotate = (any(abs(hwtrotate(1:nAtoms)) > 1.d-8).or.any(abs(hwprotate(1:nAtoms))>1.d-8))
    end if

  end subroutine initMagneticField

  subroutine lb_matrix(nAtoms, nOrbs)
    use mod_f90_kind,   only: double
    use mod_constants,  only: cZero
    use mod_parameters, only: lnolb
    implicit none
    integer, intent(in) :: nOrbs, nAtoms
    integer :: i
    complex(double), dimension(:,:), allocatable :: lbsigma

    if(allocated(lb)) deallocate(lb)
    allocate(lb(2*nOrbs,2*nOrbs,nAtoms))


    ! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
    ! to take into account the fact that we are considering negative
    ! external fields to get the peak at positive energies
    lb = cZero
    if((lfield).and.(.not.lnolb)) then
      allocate(lbsigma(nOrbs, nOrbs))
      do i=1, nAtoms
        lbsigma(1:nOrbs, 1:nOrbs) = 0.5d0*(lx*hhwx(i) + ly*hhwy(i) + lz*hhwz(i))
        lb(1:nOrbs, 1:nOrbs, i) = lbsigma(:,:)
        lb(nOrbs+1:2*nOrbs, nOrbs+1:2*nOrbs, i) = lbsigma(:,:)
      end do
      deallocate(lbsigma)
    end if
    return
  end subroutine lb_matrix

  ! Spin Zeeman hamiltonian
  subroutine sb_matrix(nAtoms, nOrbs)
    use mod_f90_kind, only: double
    use mod_constants, only: cZero, cI
    implicit none
    integer, intent(in) :: nAtoms, nOrbs
    integer :: i,mu,nu

    if(allocated(sb)) deallocate(sb)
    allocate(sb(2*nOrbs,2*nOrbs,nAtoms))

    ! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
    ! to take into account the fact that we are considering negative
    ! external fields to get the peak at positive energies
    sb = cZero

    if(lfield) then
      do i=1, nAtoms
        do mu=1,nOrbs
          nu=mu+nOrbs
          sb(mu,mu,i) = hhwz(i)
          sb(nu,nu,i) =-hhwz(i)
          sb(nu,mu,i) = hhwx(i)-cI*hhwy(i)
          sb(mu,nu,i) = hhwx(i)+cI*hhwy(i)
        end do
      end do
    end if
    return
  end subroutine sb_matrix

  ! This subroutine calculate the orbital angular momentum matrix in the cubic system of coordinates
  subroutine l_matrix()
    use mod_f90_kind
    use mod_constants
    use TightBinding, only: nOrb
    implicit none
    complex(double), dimension(9,9) :: Lp,Lm

    allocate(lx(nOrb, nOrb), ly(nOrb,nOrb), lz(nOrb,nOrb))

    lz = cZero

    lz(2,3) = -cI
    lz(3,2) = cI

    lz(5,8) = 2.d0*cI
    lz(8,5) = -2.d0*cI
    lz(6,7) = cI
    lz(7,6) = -cI

    Lp = cZero
    Lm = cZero

    Lp(2,4) = -cOne
    Lp(3,4) = -cI
    Lp(4,2) = cOne
    Lp(4,3) = cI

    Lp(5,6) = -cOne
    Lp(5,7) = -cI

    Lp(6,5) = cOne
    Lp(6,8) = -cI
    Lp(6,9) = -sq3*cI

    Lp(7,5) = cI
    Lp(7,8) = cOne
    Lp(7,9) = -sq3

    Lp(8,6) = cI
    Lp(8,7) = -cOne

    Lp(9,6) = sq3*cI
    Lp(9,7) = sq3

    Lm = transpose(conjg(Lp))

    lx = 0.5d0*(Lp+Lm)
    ly = -0.5d0*cI*(Lp-Lm)

    return
  end subroutine l_matrix

  ! This subroutine calculate the orbital angular momentum matrix in the spin system of coordinates
  subroutine lp_matrix(theta, phi)
    use mod_f90_kind, only: double
    use TightBinding, only: nOrb
    use mod_System, only: s => sys
    implicit none
    real(double), dimension(s%nAtoms), intent(in) :: theta, phi
    integer :: i

    if(allocated(lxp)) deallocate(lxp)
    if(allocated(lyp)) deallocate(lyp)
    if(allocated(lzp)) deallocate(lzp)
    allocate(lxp(nOrb,nOrb,s%nAtoms), lyp(nOrb,nOrb,s%nAtoms), lzp(nOrb,nOrb,s%nAtoms))

    do i = 1, s%nAtoms
      lxp(:,:,i) = (lx*cos(theta(i))*cos(phi(i)))+(ly*cos(theta(i))*sin(phi(i)))-(lz*sin(theta(i)))
      lyp(:,:,i) =-(lx*sin(phi(i)))+(ly*cos(phi(i)))
      lzp(:,:,i) = (lx*sin(theta(i))*cos(phi(i)))+(ly*sin(theta(i))*sin(phi(i)))+(lz*cos(theta(i)))
    end do

    return
  end subroutine lp_matrix

  subroutine allocate_magnet_variables(nAtoms)
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms
    integer :: AllocateStatus
    allocate(eps1(nAtoms), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: eps1")

    allocate( mx(nAtoms), my(nAtoms), mz(nAtoms), &
              mvec_cartesian(nAtoms,3), &
              mvec_spherical(nAtoms,3), &
              hdel(nAtoms),hdelp(nAtoms),hdelm(nAtoms), &
              mp(nAtoms),mm(nAtoms), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm")

    allocate( mabs(nAtoms), mtheta(nAtoms), mphi(nAtoms), &
              labs(nAtoms), ltheta(nAtoms), lphi(nAtoms), &
              lpabs(nAtoms), lptheta(nAtoms), lpphi(nAtoms), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi")

    return
  end subroutine

  subroutine deallocate_magnet_variables()
    implicit none

    if(allocated(eps1)) deallocate(eps1)
    if(allocated(mx)) deallocate(mx)
    if(allocated(my)) deallocate(my)
    if(allocated(mz)) deallocate(mz)
    if(allocated(mm)) deallocate(mm)
    if(allocated(mp)) deallocate(mp)
    if(allocated(hdel)) deallocate(hdel)
    if(allocated(hdelm)) deallocate(hdelm)
    if(allocated(hdelp)) deallocate(hdelp)
    if(allocated(mvec_spherical)) deallocate(mvec_spherical)
    if(allocated(mvec_cartesian)) deallocate(mvec_cartesian)
    if(allocated(mabs)) deallocate(mabs)
    if(allocated(mtheta)) deallocate(mtheta)
    if(allocated(mphi)) deallocate(mphi)
    if(allocated(labs)) deallocate(labs)
    if(allocated(ltheta)) deallocate(ltheta)
    if(allocated(lphi)) deallocate(lphi)
    if(allocated(lpabs)) deallocate(lpabs)
    if(allocated(lptheta)) deallocate(lptheta)
    if(allocated(lpphi)) deallocate(lpphi)
    if(allocated(lxm)) deallocate(lxm)
    if(allocated(lym)) deallocate(lym)
    if(allocated(lzm)) deallocate(lzm)
    if(allocated(lxpm)) deallocate(lxpm)
    if(allocated(lypm)) deallocate(lypm)
    if(allocated(lzpm)) deallocate(lzpm)
    if(allocated(hhwx)) deallocate(hhwx)
    if(allocated(hhwy)) deallocate(hhwy)
    if(allocated(hhwz)) deallocate(hhwz)
    if(allocated(sb)) deallocate(sb)
    if(allocated(lb)) deallocate(lb)

    return
  end subroutine

end module mod_magnet
