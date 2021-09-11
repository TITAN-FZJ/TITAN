module mod_magnet
  use mod_kind,       only: dp, int32
  use mod_parameters, only: dmax
  implicit none

  logical                 :: lfield
  !! Turn on/off static magnetic field (default: false)
  logical                 :: lconstraining_field
  !! Turn on/off constraining field (default: false)
  integer(int32)          :: constr_type = -1
  !! Type of constraning field to use (1 - transverse, 2 - longitudinal, 3 - full)
  real(dp)                :: cmix
  !! constraining fields mixing parameter (default 0.01)
  integer(int32)          :: iter = 0
  !! self-consistency iteration
  integer(int32)          :: maxiter
  !! Number of maximum iterations (useful for profiling and tracing)
  real(dp)                :: rhot
  !! Total occupation (summed over sites and orbitals)
  real(dp),allocatable    :: rho(:,:)
  !! orbital-dependent and d-orbital charge density per site
  real(dp),allocatable    :: mx(:,:),my(:,:),mz(:,:)
  !! orbital-dependent magnetization per site in cartesian coordinates
  real(dp),allocatable    :: rhos(:),rhop(:)
  !! s and p -orbital charge density
  real(dp),allocatable    :: rhod(:),mxd(:),myd(:),mzd(:),mzd0(:)
  !! d-orbital charge density and magnetization per site and initial values
  complex(dp),allocatable :: mp(:,:),mpd(:),mpd0(:)
  !! circular components (plus) of the total and d-orbital magnetization and intial value
  real(dp), dimension(:,:), allocatable :: rho0
  !! initial occupation per orbital obtained from hopping parameters only (to be used in the hamiltonian)
  real(dp), dimension(:), allocatable :: rhod0
  !! initial occupation of d orbitals obtained from hopping parameters only (to be used in the hamiltonian)
  real(dp),allocatable    :: lxm(:),lym(:),lzm(:)
  !! Orbital angular momentum in global frame of reference
  real(dp),allocatable    :: lxpm(:),lypm(:),lzpm(:)
  !! Orbital angular momentum in local frame of reference
  real(dp),allocatable    :: mabs(:),mtheta(:),mphi(:),mabsd(:),m_fix(:,:),m_fix_abs(:)
  !! Site-dependent spherical components of magnetization
  real(dp), dimension(:,:),allocatable   :: mvec_spherical,mvec_cartesian,mdvec_cartesian
  !! Magnetization vectors in cartesian and spherical components
  real(dp)                :: mtotal_cartesian(3),mtotal_spherical(3)
  !! Total magnetization in cartesian and spherical components
  real(dp),allocatable    :: labs(:),ltheta(:),lphi(:)
  !! Site-dependent spherical components of orbital magnetization in global frame
  real(dp),allocatable    :: lpabs(:),lptheta(:),lpphi(:)
  !! Site-dependent spherical components of orbital magnetization in local frame
  real(dp), dimension(:,:), allocatable :: hhw
  !! Half of Static magnetic fields in each direction
  real(dp),allocatable    :: bc(:,:)
  !! Constraining fields
  logical :: lrot = .false.
  !! Logical variable that indicates the need to rotate susceptibility matrix

  !========================================================================================!
  ! Values of magnetic field in cartesian or spherical coordinates
  !integer, parameter :: dmax = 20
  character(len=9), dimension(7) :: dcfield = ["hwa      ","hwt      ","hwp      ","hwahwt   ","hwahwp   ","hwthwp   ","hwahwthwp"]
  character(len=60) :: dc_header
  character(len=60), allocatable :: dc_fields(:),dcprefix(:)
  real(dp), allocatable :: hw_list(:,:)
  integer(int32) :: dcfield_dependence = 0, dc_count = 0

  real(dp) :: hwx = 0._dp, hwy = 0._dp, hwz = 0._dp, tesla = 1._dp

  integer(int32) :: hwa_npts = 0, hwa_npt1 = 1
  integer(int32) :: hwt_npts = 0, hwt_npt1 = 1
  integer(int32) :: hwp_npts = 0, hwp_npt1 = 1
  real(dp) :: hwa, hwa_i = 0._dp, hwa_f = 0._dp, hwa_s = 0._dp
  real(dp) :: hwt, hwt_i = 0._dp, hwt_f = 0._dp, hwt_s = 0._dp
  real(dp) :: hwp, hwp_i = 0._dp, hwp_f = 0._dp, hwp_s = 0._dp

  ! Layer-resolved scale of magnetic field (including empty spheres)
  logical :: lhwscale        = .false.
  real(dp) :: hwscale(dmax)   = 1._dp
  logical :: lhwrotate       = .false.
  real(dp) :: hwtrotate(dmax) = 0._dp, hwprotate(dmax) = 0._dp

  integer(int32) :: skip_steps_hw = 0
  !! How many iterations are to be skipped from the beginning
  integer(int32) :: hw_count
  !! Current iteration of magnetic field loop.
  integer(int32) :: total_hw_npt1
  !! Total number of magnetic field iterations.

contains

  ! This subroutine sets up external magnetic fields and related loop
  subroutine setMagneticLoopPoints()
    use mod_constants, only: rad2deg
    implicit none
    !! Amount of points to skip
    integer       :: i, j, k, l

    if(lfield) then
      if (hwa_npts==0) then
        hwa_f = hwa_i
        hwa_npts = 1
      end if
      if((abs(hwa_i) + abs(hwa_f)) > 1.e-8_dp) then
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
        if(abs(hwa_s) <= 1.e-15_dp) hwa_npt1 = 1
        hwt_s = (hwt_f - hwt_i)/hwt_npts
        if(abs(hwt_s) <= 1.e-15_dp) hwt_npt1 = 1
        hwp_s = (hwp_f - hwp_i)/hwp_npts
        if(abs(hwp_s) <= 1.e-15_dp) hwp_npt1 = 1
      else ! hwa_i and hwa_f = 0
        hwa_i   = sqrt(hwx**2+hwy**2+hwz**2)
        hwa_f   = hwa_i
        if(abs(hwa_i)<1.e-8_dp) then
          lfield = .false.
        else
          hwt_i    = acos(hwz/hwa_i)*rad2deg
          hwt_s    = 0._dp
          hwt_npt1 = 1
          hwp_i    = atan2(hwy,hwx)*rad2deg
          hwp_s    = 0._dp
          hwp_npt1 = 1
        end if
      end if
    else ! lfield
      hwa_i    = 0._dp
      hwa_f    = 0._dp
      hwa_s    = 0._dp
      hwa_npt1 = 1
      hwt_i    = 0._dp
      hwt_f    = 0._dp
      hwt_s    = 0._dp
      hwt_npt1 = 1
      hwp_i    = 0._dp
      hwp_f    = 0._dp
      hwp_s    = 0._dp
      hwp_npt1 = 1
    end if

    ! Total number of points in the loops
    total_hw_npt1 = hwa_npt1*hwt_npt1*hwp_npt1
    if(total_hw_npt1==1) skip_steps_hw = 0
    total_hw_npt1 = total_hw_npt1 - skip_steps_hw

    allocate(hw_list(total_hw_npt1,3))

    ! Creating list of magnetic fields (in spherical coordinates)
    i = 1
    do j = 1, hwa_npt1
      hwa = hwa_i + (j-1)*hwa_s
      do k=1, hwt_npt1
        hwt = hwt_i + (k-1)*hwt_s
        do l = 1, hwp_npt1
          if(i <= skip_steps_hw) cycle
          hwp = hwp_i + (l-1)*hwp_s
          hw_list(i - skip_steps_hw, :) = [ hwa , hwt , hwp ]
          i = i+1
        end do
      end do
    end do

  end subroutine setMagneticLoopPoints

  subroutine initMagneticField(s)
    use mod_constants, only: cZero, deg2rad
    use mod_mpi_pars,  only: abortProgram
    use mod_System,    only: System_type
    implicit none
    type(System_type), intent(inout) :: s

    integer :: i, AllocateStatus

    if(allocated(hhw)) deallocate(hhw)
    allocate( hhw(3,s%nAtoms), STAT = AllocateStatus )
    if (AllocateStatus /= 0) call abortProgram("[main] Not enough memory for: hhw")

    !---------------------- Turning off field for hwa=0 --------------------
    if(abs(hw_list(hw_count,1)) < 1.e-8_dp) then
      lfield = .false.
    else
      lfield = .true.
    end if

    do i = 1, s%nAtoms
      s%Basis(i)%sb(:,:) = cZero
      s%Basis(i)%lb(:,:) = cZero
    end do
    hhw  = 0._dp

    !--------------------- Defining the magnetic fields -------------------- TODO
    if(lfield) then
      ! Variables of the hamiltonian
      ! There is an extra  minus sign in the definition of hhw
      ! to take into account the fact that we are considering negative
      ! external fields to get the peak at positive energies
      do i = 1, s%nAtoms
        hhw(1,i) = -0.5_dp*hwscale(i)*hw_list(hw_count,1)*sin((hw_list(hw_count,2) + hwtrotate(i))*deg2rad)*cos((hw_list(hw_count,3) + hwprotate(i))*deg2rad)*tesla
        hhw(2,i) = -0.5_dp*hwscale(i)*hw_list(hw_count,1)*sin((hw_list(hw_count,2) + hwtrotate(i))*deg2rad)*sin((hw_list(hw_count,3) + hwprotate(i))*deg2rad)*tesla
        hhw(3,i) = -0.5_dp*hwscale(i)*hw_list(hw_count,1)*cos((hw_list(hw_count,2) + hwtrotate(i))*deg2rad)*tesla
        !          ^ EXTRA MINUS SIGN
        if(abs(hhw(1,i))<1.e-8_dp) hhw(1,i) = 0._dp
        if(abs(hhw(2,i))<1.e-8_dp) hhw(2,i) = 0._dp
        if(abs(hhw(3,i))<1.e-8_dp) hhw(3,i) = 0._dp
      end do
      ! Testing if hwscale is used
      lhwscale = any(abs(hwscale(1:s%nAtoms)-1._dp) > 1.e-8_dp)
      ! Testing if hwrotate is used
      lhwrotate = (any(abs(hwtrotate(1:s%nAtoms)) > 1.e-8_dp).or.any(abs(hwprotate(1:s%nAtoms))>1.e-8_dp))
    end if
  end subroutine initMagneticField

  subroutine lb_matrix(s)
  !! Orbital Zeeman term
    use mod_parameters, only: lnolb
    use mod_System,     only: System_type
    implicit none
    type(System_type), intent(inout)   :: s
    integer :: i,mu,nu,mup,nup,sigma

    ! There is an extra  minus sign in the definition of hhw
    ! to take into account the fact that we are considering negative
    ! external fields to get the peak at positive energies
    if((lfield).and.(.not.lnolb)) then
      do i=1, s%nAtoms

        do nu=1,s%Types(s%Basis(i)%Material)%nOrb
          nup = nu+s%Types(s%Basis(i)%Material)%nOrb
          do mu=1,s%Types(s%Basis(i)%Material)%nOrb
            mup = mu+s%Types(s%Basis(i)%Material)%nOrb

            ! L.B
            do sigma=1,3
              s%Basis(i)%lb(mu ,nu ) = s%Basis(i)%lb(mu ,nu ) + s%Types(s%Basis(i)%Material)%lvec(mu,nu,sigma)*hhw(sigma,i)
            end do
            s%Basis(i)%lb(mup,nup) = s%Basis(i)%lb(mu,nu)
          end do
        end do
      end do
    end if
  end subroutine lb_matrix


  subroutine sb_matrix(s)
  !! Spin Zeeman hamiltonian
    use mod_constants, only: cZero, cI
    use mod_System,    only: System_type
    implicit none
    type(System_type), intent(inout)   :: s
    real(dp), dimension(3,s%nAtoms) :: b
    integer :: i,mu,nu

    ! There is an extra  minus sign in the definition of hhw
    ! to take into account the fact that we are considering negative
    ! external fields to get the peak at positive energies

    if(lfield.or.lconstraining_field) then
      ! write(output%unit_loop,"('[sb_matrix] hhw,bc = ', 10000es16.8)") hhw, bc
      b = hhw - 0.5_dp*bc
      do i=1, s%nAtoms
        s%Basis(i)%sb(:,:) = cZero
        do mu=1,s%Types(s%Basis(i)%Material)%nOrb
          nu=mu+s%Types(s%Basis(i)%Material)%nOrb
          s%Basis(i)%sb(mu,mu) = b(3,i)
          s%Basis(i)%sb(nu,nu) =-b(3,i)
          s%Basis(i)%sb(nu,mu) = b(1,i)+cI*b(2,i)
          s%Basis(i)%sb(mu,nu) = b(1,i)-cI*b(2,i)
        end do
      end do
    end if
  end subroutine sb_matrix

  ! This subroutine calculate the orbital angular momentum matrix 
  ! for all orbitals, in the cubic system of coordinates
  subroutine l_matrix(s)
    use mod_kind,      only: dp
    use mod_constants, only: cI, cZero, cOne, sq3
    use AtomTypes,     only: default_orbitals
    use mod_System,    only: System_type
    implicit none
    type(System_type), intent(inout)   :: s
    complex(dp), dimension(size(default_orbitals),size(default_orbitals),3) :: lvec
    !! Angular momentum vector matrices in global frame
    complex(dp), dimension(size(default_orbitals),size(default_orbitals)) :: Lp,Lm
    !! L+ and L- matrices in global frame
    integer :: i, mu, nu

    lvec(:,:,3) = cZero

    lvec(2,3,3) = -cI
    lvec(3,2,3) = cI

    lvec(5,8,3) = 2._dp*cI
    lvec(8,5,3) = -2._dp*cI
    lvec(6,7,3) = cI
    lvec(7,6,3) = -cI

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

    lvec(:,:,1) = 0.5_dp*   (Lp(:,:)+Lm(:,:))
    lvec(:,:,2) =-0.5_dp*cI*(Lp(:,:)-Lm(:,:))

    do i = 1, s%nTypes
      do mu = 1, s%Types(i)%nOrb
        do nu = 1, s%Types(i)%nOrb
          s%Types(i)%lvec(mu,nu,:) = lvec(s%Types(i)%Orbs(mu),s%Types(i)%Orbs(nu),:)
        end do
      end do
    end do

  end subroutine l_matrix

  ! This subroutine calculate the orbital angular momentum matrix in the local system of coordinates
  subroutine lp_matrix(theta, phi)
    use mod_kind,      only: dp
    use mod_constants, only: deg2rad
    use mod_System,    only: s => sys
    implicit none
    real(dp), dimension(s%nAtoms), intent(in) :: theta, phi
    integer :: i

    do i = 1, s%nAtoms
      s%Basis(i)%lpvec(:,:,1) = ( s%Types(s%Basis(i)%Material)%lvec(:,:,1)*cos(theta(i)*deg2rad)*cos(phi(i)*deg2rad)) &
                              + ( s%Types(s%Basis(i)%Material)%lvec(:,:,2)*cos(theta(i)*deg2rad)*sin(phi(i)*deg2rad)) &
                              - ( s%Types(s%Basis(i)%Material)%lvec(:,:,3)*sin(theta(i)*deg2rad))
      s%Basis(i)%lpvec(:,:,2) = (-s%Types(s%Basis(i)%Material)%lvec(:,:,1)                      *sin(phi(i)*deg2rad)) &
                              + ( s%Types(s%Basis(i)%Material)%lvec(:,:,2)                      *cos(phi(i)*deg2rad))
      s%Basis(i)%lpvec(:,:,3) = ( s%Types(s%Basis(i)%Material)%lvec(:,:,1)*sin(theta(i)*deg2rad)*cos(phi(i)*deg2rad)) &
                              + ( s%Types(s%Basis(i)%Material)%lvec(:,:,2)*sin(theta(i)*deg2rad)*sin(phi(i)*deg2rad)) &
                              + ( s%Types(s%Basis(i)%Material)%lvec(:,:,3)*cos(theta(i)*deg2rad))
    end do
  end subroutine lp_matrix

  subroutine allocate_magnet_variables(nAtoms,nOrbs)
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms
    integer, intent(in) :: nOrbs
    integer :: AllocateStatus

    allocate( rho(nOrbs,nAtoms), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocate_magnet_variables] Not enough memory for: rho")

    allocate( mx(nOrbs,nAtoms), my(nOrbs,nAtoms), mz(nOrbs,nAtoms), mp(nOrbs,nAtoms), &
              mvec_cartesian(3,nAtoms),mdvec_cartesian(3,nAtoms), &
              mvec_spherical(3,nAtoms), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[allocate_magnet_variables] Not enough memory for: mx,my,mz,mp,mvec_cartesian,mdvec_cartesian,mvec_spherical")

    allocate( mxd(nAtoms), myd(nAtoms), mzd(nAtoms), mpd(nAtoms), &
              rhos(nAtoms), rhop(nAtoms), rhod(nAtoms), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[allocate_magnet_variables] Not enough memory for: mxd,myd,mzd,mpd,rhos,rhop,rhod")

    allocate( mabs(nAtoms), mtheta(nAtoms), mphi(nAtoms), &
              mabsd(nAtoms), m_fix(3,nAtoms), m_fix_abs(nAtoms), &
              labs(nAtoms), ltheta(nAtoms), lphi(nAtoms), &
              lpabs(nAtoms), lptheta(nAtoms), lpphi(nAtoms), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[allocate_magnet_variables] Not enough memory for: mabs,mabsd,m_fix,m_fix_abs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi")

    allocate( hhw(3,nAtoms), STAT = AllocateStatus )
    if (AllocateStatus /= 0) call abortProgram("[allocate_magnet_variables] Not enough memory for: hhw")

    allocate( bc(3,nAtoms), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[allocate_magnet_variables] Not enough memory for: bc")
    bc(:,:) = 0.e0_dp

  end subroutine allocate_magnet_variables

  subroutine deallocate_magnet_variables()
    implicit none

    if(allocated(rho)) deallocate(rho)
    if(allocated(mx)) deallocate(mx)
    if(allocated(my)) deallocate(my)
    if(allocated(mz)) deallocate(mz)
    if(allocated(mp)) deallocate(mp)
    if(allocated(mvec_cartesian)) deallocate(mvec_cartesian)
    if(allocated(mdvec_cartesian)) deallocate(mdvec_cartesian)
    if(allocated(mvec_spherical)) deallocate(mvec_spherical)
    if(allocated(mxd)) deallocate(mxd)
    if(allocated(myd)) deallocate(myd)
    if(allocated(mzd)) deallocate(mzd)
    if(allocated(mpd)) deallocate(mpd)
    if(allocated(rhos)) deallocate(rhos)
    if(allocated(rhop)) deallocate(rhop)
    if(allocated(rhod)) deallocate(rhod)
    if(allocated(mabs)) deallocate(mabs)
    if(allocated(mabsd)) deallocate(mabsd)
    if(allocated(m_fix)) deallocate(m_fix)
    if(allocated(m_fix_abs)) deallocate(m_fix_abs)
    if(allocated(mtheta)) deallocate(mtheta)
    if(allocated(mphi)) deallocate(mphi)
    if(allocated(labs)) deallocate(labs)
    if(allocated(ltheta)) deallocate(ltheta)
    if(allocated(lphi)) deallocate(lphi)
    if(allocated(lpabs)) deallocate(lpabs)
    if(allocated(lptheta)) deallocate(lptheta)
    if(allocated(lpphi)) deallocate(lpphi)
    if(allocated(hhw)) deallocate(hhw)
    if(allocated(bc)) deallocate(bc)

    if(allocated(lxm)) deallocate(lxm)
    if(allocated(lym)) deallocate(lym)
    if(allocated(lzm)) deallocate(lzm)
    if(allocated(lxpm)) deallocate(lxpm)
    if(allocated(lypm)) deallocate(lypm)
    if(allocated(lzpm)) deallocate(lzpm)
  end subroutine

  subroutine set_fieldpart(kount)
    use mod_parameters, only: ltesla, lnolb, output
    use mod_tools,      only: rtos
    implicit none
    integer :: kount

    output%BField = ""
    if(lfield) then
      write(output%BField, "('_hwa=',a,'_hwt=',a,'_hwp=',a)") trim(rtos(hw_list(kount,1),"(es9.2)")),trim(rtos(hw_list(kount,2),"(f7.2)")),trim(rtos(hw_list(kount,3),"(f7.2)"))
      if(ltesla)    output%BField = trim(output%BField) // "_tesla"
      if(lnolb)     output%BField = trim(output%BField) // "_nolb"
      if(lhwscale)  output%BField = trim(output%BField) // "_hwscale"
      if(lhwrotate) output%BField = trim(output%BField) // "_hwrotate"
    end if

    output%dcBField = ""
    if(dcfield_dependence/=7) then
      if((dcfield_dependence/=1).and.(dcfield_dependence/=4).and.(dcfield_dependence/=5)) write(output%dcBField,"(a,'_hwa=',a)") trim(output%dcBField),trim(rtos(hwa,"(es9.2)"))
      if((dcfield_dependence/=2).and.(dcfield_dependence/=4).and.(dcfield_dependence/=6)) write(output%dcBField,"(a,'_hwt=',a)") trim(output%dcBField),trim(rtos(hwt,"(f7.2)"))
      if((dcfield_dependence/=3).and.(dcfield_dependence/=5).and.(dcfield_dependence/=6)) write(output%dcBField,"(a,'_hwp=',a)") trim(output%dcBField),trim(rtos(hwp,"(f7.2)"))
    end if
    if(ltesla)    output%dcBField = trim(output%dcBField) // "_tesla"
    if(lnolb)     output%dcBField = trim(output%dcBField) // "_nolb"
    if(lhwscale)  output%dcBField = trim(output%dcBField) // "_hwscale"
    if(lhwrotate) output%dcBField = trim(output%dcBField) // "_hwrotate"
  end subroutine set_fieldpart

end module mod_magnet
