program DHE
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_io
  use mod_magnet
  use mod_tight_binding
  use mod_prefactors
  use mod_generate_epoints
  use mod_generate_kpoints
  use mod_lattice
  use mod_progress
  use mod_mpi_pars
  use mod_dnsqe
  use mod_susceptibilities
  use mod_disturbances
  use mod_currents
  use mod_beff
  use mod_torques
  use MPI
  implicit none
  character(len=400)            :: varm
  character(len=8)              :: date
  character(len=10)             :: time
  character(len=5)              :: zone
  integer                       :: i,j,iw,err,sigma,sigmap,mu,nu,gamma,xi,neighbor,iflag
  integer                       :: AllocateStatus,values(8),lwa,ifail=0
  real(double)                  :: e,Icabs

  ! Spin diffusion length - not working?
!   integer                       :: pos,inn
!   real(double)                  :: sdl2
!   complex(double),allocatable   :: sdl(:)

  ! Exchange interaction
  real(double),dimension(:,:),allocatable     :: trJij
  real(double),dimension(:,:,:,:),allocatable :: Jij,Jijs,Jija
  ! LDOS
  real(double),dimension(:,:),allocatable :: ldosu,ldosd
  ! Self consistency variables
  logical                       :: selfcon
  real(double),allocatable      :: fvec(:),jac(:,:),wa(:),sc_solu(:)
  real(double),allocatable      :: diag(:),w(:,:),qtf(:)
  real(double)                  :: epsfcn,factor
  integer                       :: neq,maxfev,ml,mr,mode,nfev,njev

  ! Identity and U matrix
  complex(double), dimension(:,:),   allocatable :: identt,Umatorb,temp
!---------------------------- begin MPI vars ---------------------------
#ifndef _UFF
!$  integer :: provided
#endif
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^  end MPI vars ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  external selfconsistency,selfconsistencyjac,selfconsistencyjacnag
!------------------------ begin MPI initialization ---------------------
#ifdef _OPENMP
#ifdef _UFF
  call MPI_Init(ierr)
#else
!   call MPI_Init_thread(MPI_THREAD_MULTIPLE,provided,ierr)
  call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
#endif
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
#else
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
#endif
!-------------------- Useful constants and matrices --------------------
  pi  = 4.d0*atan(1.d0)
  tpi = 2.d0*pi
  sq2 = sqrt(2.d0)
  hsq2 = 0.5d0*sq2
  sq3 = sqrt(3.d0)
  identorb9  = zero
  identorb18 = zero
  do i=1,9
    identorb9(i,i) = zum
  end do
  do i=1,18
    identorb18(i,i) = zum
  end do
!------------------------- Reading parameters --------------------------
  if(myrank.eq.0) then
    call date_and_time(date, time, zone, values)
    write(*,"('[main] Started on: ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") values(3),values(2),values(1),values(5),values(6),values(7)
    start_program = MPI_Wtime()
  end if
  call ioread()
!---------------------------- Getting hostname -------------------------
#ifdef _JUQUEEN
  call hostnm_(host)
#else
  call hostnm(host)
#endif
!-------------------- Define the lattice structure ---------------------
!-------------------- Generating k points in 2D BZ ---------------------
  select case (lattice)
  case("bcc110")
    call bcc110()
    call generate_kpoints_bcc110()
  case("fcc100")
    call fcc100()
    call generate_kpoints_fcc100()
  case("fcc111")
    call fcc111()
    call generate_kpoints_fcc111()
  case default
    if(myrank.eq.0) write(*,"('[main] Lattice not defined: ',a,'!')") lattice
    stop
  end select
  ! Writing BZ points and weights into files
  if((lkpoints).and.(myrank.eq.0)) then
    open (unit=2222, file='kpoints2d',status='unknown')
    write(unit=2222,fmt="(a)") ' #      kx            ky            wk'
    open (unit=3333, file='kpoints3d',status='unknown')
    write(unit=3333,fmt="(a)") ' #      kx            ky            kz            wk'
    do i=1,nkpoints
      write(unit=2222,fmt="(3(f12.9,2x))") kbz2d(i,1),kbz2d(i,2),wkbz(i)
      write(unit=3333,fmt="(4(f12.9,2x))") kbz(i,1),kbz(i,2),kbz(i,3),wkbz(i)
    end do
    close(2222)
    close(3333)
  end if
!---- Generating integration points of the complex energy integral -----
  call generate_imag_epoints()
!--------------------------------- Loops -------------------------------
  number_of_planes: do Npl = Npl_i,Npl_f
    hw_intensity: do hwa_count=1,hwa_npt1
      hwa = hwa_i + (hwa_count-1)*hwa_s
      if(abs(hwa).lt.1.d-8) then
        FIELD = .false.
      else
        FIELD = .true.
      end if
      hw_theta: do hwt_count=1,hwt_npt1
        hwt = hwt_i + (hwt_count-1)*hwt_s
        hw_phi: do hwp_count=1,hwp_npt1
          hwp = hwp_i + (hwp_count-1)*hwp_s
!--------------------- Defining the magnetic fields --------------------
    if(FIELD) then
      hwx  = hwa*sin(hwt*pi)*cos(hwp*pi)
      hwy  = hwa*sin(hwt*pi)*sin(hwp*pi)
      hwz  = hwa*cos(hwt*pi)
      if(abs(hwx).lt.1.d-8) hwx = 0.d0
      if(abs(hwy).lt.1.d-8) hwy = 0.d0
      if(abs(hwz).lt.1.d-8) hwz = 0.d0
      ! Variables of the hamiltonian
      hhwx  = 0.5d0*hwx*tesla
      hhwy  = 0.5d0*hwy*tesla
      hhwz  = 0.5d0*hwz*tesla
    else
      if((hwt_count.ge.2).or.(hwp_count.ge.2)) exit
    end if
    dclimit_field_dependence(1) = hwa
    dclimit_field_dependence(2) = hwt
    dclimit_field_dependence(3) = hwp
!--------------- Allocating variables that depend on Npl ---------------
    neq = 4*Npl
!     lwa=neq*(3*neq+13)/2
    lwa=neq*(neq+1)/2
    allocate( sigmai2i(4,Npl),sigmaimunu2i(4,Npl,9,9),sigmaijmunu2i(4,Npl,Npl,9,9),eps1(Npl),sc_solu(neq), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1,sc_solu')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    allocate( diag(neq),qtf(neq),w(neq,4), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: diag,qtf,w')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    allocate( mx(Npl),my(Npl),mz(Npl),hdel(Npl),mp(Npl),hdelp(Npl),mm(Npl),hdelm(Npl),fvec(neq),jac(neq,neq),wa(lwa), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: mx,my,mz,hdel,mp,hdelp,mm,hdelm,fvec,jac,wa')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    allocate( mabs(Npl),mtheta(Npl),mphi(Npl),labs(Npl),ltheta(Npl),lphi(Npl),lpabs(Npl),lptheta(Npl),lpphi(Npl), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    allocate( mmlayer(Npl+2),layertype(Npl+2),U(Npl+2),mmlayermag(Npl+2),lambda(Npl+2),npart0(Npl+2), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: mmlayer,layertype,U,mmlayermag,lambda,npart0')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
!---------------------- Dimensions and identities ----------------------
    dimsigmaNpl = 4*Npl
    dim = dimsigmaNpl*9*9
    allocate( identt(dim,dim), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: identt')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    identt     = zero
    do i=1,dim
     identt(i,i) = zum
    end do
!----------- Creating bi-dimensional matrix of MPI processes  ----------
    if((itype.ge.5).and.(itype.le.8)) then ! Create matrix only when energy integration is involved
      call build_cartesian_grid()
    end if
    if(itype.eq.11) then ! Create matrix for dclimit - only one line of processes
      if(numprocs.le.pnt) then
        call build_cartesian_grid()
      else
        if(myrank.eq.0) write(*,"('[main] MPI not yet fully implemented for dc-limit calculation!')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if
!------------------------- Conversion arrays  --------------------------
    do nu=1,9 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4
      sigmaimunu2i(sigma,i,mu,nu)= (sigma-1)*Npl*9*9+(i-1)*9*9+(mu-1)*9+nu
      do j=1,Npl
        sigmaijmunu2i(sigma,i,j,mu,nu)= (sigma-1)*Npl*Npl*9*9+(i-1)*Npl*9*9+(j-1)*9*9+(mu-1)*9+nu
      end do
    end do ; end do ; end do ; end do
    do i=1,Npl ; do sigma=1,4
      sigmai2i(sigma,i) = (sigma-1)*Npl+i
    end do ; end do
!---------------------- Tight Binding parameters -----------------------
    call rs_hoppings()
!------------------------- Mounting U matrix ---------------------------
    allocate( Umatorb(dim,dim), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: Umatorb')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if

    Umatorb = zero
    if(Utype.ne.0) then
      do xi=5,9 ; do gamma=5,9 ; do nu=5,9 ; do mu=5,9 ; do i=1,Npl
        if((Utype.eq.1).and.(layertype(i+1).ne.2)) cycle
        if((mu.ne.nu).or.(gamma.ne.xi)) cycle
        Umatorb(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,i,gamma,xi)) = cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,i,gamma,xi)) = cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,i,gamma,xi)) = cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,i,gamma,xi)) = cmplx(U(i+1),0.d0)
      end do ; end do ; end do ; end do ; end do
      do xi=5,9 ; do gamma=5,9 ; do nu=5,9 ; do mu=5,9 ; do i=1,Npl
        if((Utype.eq.1).and.(layertype(i+1).ne.2)) cycle
        if((mu.ne.xi).or.(nu.ne.gamma)) cycle
        Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,i,gamma,xi)) = Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,i,gamma,xi))-cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,i,gamma,xi)) = Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,i,gamma,xi))-cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,i,gamma,xi)) = Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,i,gamma,xi))-cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,i,gamma,xi)) = Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,i,gamma,xi))-cmplx(U(i+1),0.d0)
      end do ; end do ; end do ; end do ; end do
    end if
!-----------------------------------------------------------------------
    select case (lattice)
    case("bcc110")
!------------------- Spin quantization direction -----------------------
      magnetization_axis_bcc110: select case (magaxis)
      case ("L")
    !   In-plane long axis: (x - short axis; y - out-of-plane; z - long axis)
        phi   = pi/4.d0 ! around z (counter-clockwise from the top X->Y)
        theta = pi/2.d0 ! around y' (counter-clockwise from the top Z->X')
      case ("P")
    !   Out-of-plane: (x - short axis; y - long axis; z - out-of-plane)
        phi   = -pi/4.d0 ! around z (counter-clockwise from the top X->Y)
        theta = pi/2.d0 ! around y' (counter-clockwise from the top Z->X')
      case ("S")
    !   In-plane short axis: (x - long axis; y - out-of-plane; z - short axis)
        phi   = pi/4.d0 ! around z (counter-clockwise from the top X->Y)
        theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
      case default
        if(myrank.eq.0) write(*,"('[main] Choose a correct magnetization axis!')")
        stop
    !     phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
    !     theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
      end select magnetization_axis_bcc110
!----------- Direction of applied in-plane electric field --------------
      direction_E_field_bcc110: select case (dirEfield)
      case ("L")
    !   In-plane long axis:
        dirEfieldvec = [hsq2,hsq2,0.d0]
      case ("S")
    !   In-plane short axis
        dirEfieldvec = [0.d0,0.d0,1.d0]
      case ("O")
    !   Other direction:
    !     dirEfieldvec = [1.d0/sq3,1.d0/sq3,1.d0/sq3] ! In the direction of the 1st n.n.
        dirEfieldvec = [1.d0/(sq2*sq3),1.d0/(sq2*sq3),-sq2/sq3] ! In-plane, Perpendicular to the 1st n.n.
      end select direction_E_field_bcc110
    case("fcc100")
!----------- Direction of applied in-plane electric field --------------
      read(unit=dirEfield,fmt=*) i
      direction_E_field_fcc100: select case (i)
      case (1:8)   !    In plane neighbors:
        dirEfieldvec = r0(i,:)
      case default !    Other direction:
        dirEfieldvec = [1.d0,0.d0,0.d0] ! In-plane
      end select direction_E_field_fcc100
!------------------- Spin quantization direction -----------------------
      read(unit=magaxis,fmt=*) j
      magnetization_axis_fcc100: select case (j)
      case (1:4)
    !   In-plane 1st n.n.:
        phi   = dble(2*j-1)*pi/4.d0 ! around z (counter-clockwise from the top X->Y)
        theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
      case (5:8)
    !   In-plane 2nd n.n.:
        phi   = dble(j-5)*pi/2.d0   ! around z (counter-clockwise from the top X->Y)
        theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
      case (9)
    !   Out-of-plane: (x in the direction of dirEfield)
        select case (i)             ! around z (counter-clockwise from the top X->Y)
        case (1:4)
          phi = dble(2*i-1)*pi/4.d0
        case (5:8)
          phi = dble(i-5)*pi/2.d0
        case default
          phi = 0.d0
        end select
        theta = 0.d0                ! around y' (counter-clockwise from the top Z->X')
      case default
        if(myrank.eq.0) write(*,"('[main] Choose a correct magnetization axis!')")
        stop
    !     phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
    !     theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
      end select magnetization_axis_fcc100
    case("fcc111")
!------------------- Spin quantization direction -----------------------
      read(unit=magaxis,fmt=*) j
      magnetization_axis_fcc111: select case (j)
!       case (1:4)
!     !   In-plane 1st n.n.:
!         phi   = dble(2*j-1)*pi/4.d0 ! around z (counter-clockwise from the top X->Y)
!         theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
!       case (5:8)
!     !   In-plane 2nd n.n.:
!         phi   = dble(j-5)*pi/2.d0   ! around z (counter-clockwise from the top X->Y)
!         theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
!       case (9)
!     !   Out-of-plane:
!         phi   = pi/4.d0             ! around z (counter-clockwise from the top X->Y)
!         theta = 0.d0                ! around y' (counter-clockwise from the top Z->X')
      case default
        if(myrank.eq.0) write(*,"('[main] Choose a correct magnetization axis!')")
        stop
    !     phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
    !     theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
      end select magnetization_axis_fcc111
!----------- Direction of applied in-plane electric field --------------
      read(unit=dirEfield,fmt=*) j
      direction_E_field_fcc111: select case (j)
!       case (1:8)   !    In plane neighbors:
!         dirEfieldvec = r0(j,:)
!       case default !    Other direction:
!         dirEfieldvec = [1.d0,0.d0,0.d0] ! In-plane
      end select direction_E_field_fcc111
    end select
!------- Writing parameters and data to be calculated on screen --------
    if(myrank.eq.0) call iowrite()
!---- L matrix in spin coordinates for given quantization direction ----
    call lp_matrix()
!------ Calculate L.S matrix for the given quantization direction ------
    call ls_matrix()
!------ Calculate L.B matrix for the given quantization direction ------
    call lb_matrix()
!------------------------ Calculate S.B matrix  ------------------------
    call sb_matrix()
!---------------------------- Debugging part ---------------------------
    if(ldebug) then
!       if(myrank.eq.0) then
      if(myrank.eq.0) write(*,"('Debugging...')")
      call debugging()
!       end if
      call MPI_Finalize(ierr)
      if (ierr.ne.0) then
        write(*,"('[main] ierr = ',i0,'. Something went wrong in the parallelization!')") ierr
      end if
      stop
    end if
!------------------- Only create files with headers --------------------
    if(lcreatefiles) then
      if(myrank.eq.0) then
        create_files: select case (itype)
        case (5)
          call openclose_chi_files(0)
          write(*,"('[main] Susceptibilities files created/overwritten!')")
        case (6)
          call openclose_chi_files(0)
          call openclose_disturbances_files(0)
          call openclose_beff_files(0)
          write(*,"('[main] Susceptibilities, disturbances and effective field files created/overwritten!')")
        case (7)
          call openclose_chi_files(0)
          call openclose_currents_files(0)
          write(*,"('[main] Susceptibilities and currents files created/overwritten!')")
        case (8)
          call openclose_chi_files(0)
          call openclose_disturbances_files(0)
          call openclose_currents_files(0)
          call openclose_beff_files(0)
          write(*,"('[main] Susceptibilities, disturbances, current and effective field files created/overwritten!')")
        case default
          write(*,"('[main] No files to create for selected option! (itype = ',i0,')')") itype
        end select create_files
      end if
      call MPI_Finalize(ierr)
      if (ierr.ne.0) then
        write(*,"('[main] ierr = ',i0,'. Something went wrong in the parallelization!')") ierr
      end if
      stop
    end if
!----------------- Check if files exist to add results -----------------
    if(laddresults) then
      check_files: select case (itype)
      case (5)
        call openclose_chi_files(1)
        call openclose_chi_files(2)
      case (6)
        call openclose_chi_files(1)
        call openclose_chi_files(2)
        call openclose_disturbances_files(1)
        call openclose_disturbances_files(2)
        call openclose_beff_files(1)
        call openclose_beff_files(2)
      case (7)
        call openclose_chi_files(1)
        call openclose_chi_files(2)
        call openclose_currents_files(1)
        call openclose_currents_files(2)
      case (8)
        if(.not.lhfresponses) then
          call openclose_chi_files(1)
          call openclose_chi_files(2)
        end if
        call openclose_disturbances_files(1)
        call openclose_disturbances_files(2)
        call openclose_currents_files(1)
        call openclose_currents_files(2)
        call openclose_beff_files(1)
        call openclose_beff_files(2)
        call openclose_torques_files(1)
        call openclose_torques_files(2)
      case (11)
        if(hwa_count*hwt_count*hwp_count.eq.1) then
          FIELD = .true.
          call openclose_dclimit_disturbances_files(1)
          call openclose_dclimit_disturbances_files(2)
          call openclose_dclimit_currents_files(1)
          call openclose_dclimit_currents_files(2)
          call openclose_dclimit_beff_files(1)
          call openclose_dclimit_beff_files(2)
          call openclose_dclimit_torques_files(1)
          call openclose_dclimit_torques_files(2)
        end if
      end select check_files
    end if
!-------------------------- Begin first test part ----------------------
    if((myrank.eq.0).and.(itype.eq.0)) then
      write(*,"('FIRST TEST PART')")

      allocate(ldosu(Npl,9),ldosd(Npl,9))
      allocate(trJij(nmaglayers,nmaglayers),Jij(nmaglayers,nmaglayers,3,3),Jijs(nmaglayers,nmaglayers,3,3),Jija(nmaglayers,nmaglayers,3,3))

      ! Parameters: center of band, magnetization, correlation functions
      eps1  = 0.d0
      mz   = 0.d0
      mp = zero
      ! Variables used in the hamiltonian
      do i=1,Npl
        hdel(i)   = 0.5d0*U(i+1)*mz(i)
        hdelp(i)  = 0.5d0*U(i+1)*mp(i)
      end do
      hdelp = zero
      hdelm = conjg(hdelp)

      emin = -2.d0  ! given in eV
      emax = 2.d0   ! given in eV
      npts = 400
      npt1 = npts+1
      test1_energy_loop2: do count=1,npt1
        e = emin + (count-1)*deltae
        e = (e/ry2ev)!+Ef ! Transform to Ry
        write(*,"(i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e

        ! Turning off SOC
!         lambda = 0.d0

        call ldos(e,ldosu,ldosd,Jij)

      end do test1_energy_loop2

      deallocate(ldosu,ldosd)
      deallocate(trJij,Jij,Jijs,Jija)

  !   Finalizing program
      if(myrank.eq.0) then
        call date_and_time(date, time, zone, values)
        write(*,"('[main] Finished on: ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") values(3),values(2),values(1),values(5),values(6),values(7)
        elapsed_time = MPI_Wtime() - start_program
        write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
      end if
      call MPI_Finalize(ierr)
      if (ierr.ne.0) then
        write(*,"('[main] ierr = ',i0,'. Something went wrong in the parallelization.')") ierr
      end if
      stop

    end if
!------------------------ Begin self-consistency -----------------------
  ! Reading previous shifts and mz from files
  ! (Try to read eps1 and mz if available - includes hdel, hdelp and hdelm calculations)
    iflag = 0
    call readwritesc(iflag,err)

    sc_file_status: select case (iflag)
    case(1)
      if(err.eq.0) then ! Same parameters
        if(skipsc) then ! Skip option ON
          if(myrank.eq.0) write(*,"('Existing results for the same parameters were read. Skipping self-consistency...')")
          selfcon = .false.
        else ! Skip option OFF
          if(myrank.eq.0) write(*,"('Existing results for the same parameters were read. Updating values...')")
          selfcon = .true.
        end if
      else ! Other parameters
        if(myrank.eq.0) write(*,"('Existing results for other parameters were read. Updating values...')")
        selfcon = .true.
      end if
      ! Putting read eps1 existing solutions into esp1_solu (first guess of the subroutine)
      if(selfcon) then
        sc_solu(1:Npl)         = eps1
        sc_solu(Npl+1:2*Npl)   = mx
        sc_solu(2*Npl+1:3*Npl) = my
        sc_solu(3*Npl+1:4*Npl) = mz
      end if
    case default !  If file doesn't exist
      if(myrank.eq.0) write(*,"('Self-consistency file does not exist.')")
      selfcon = .true.
      ! Parameters: center of band, magnetization, correlation functions
      eps1 = 0.d0
      mx = 0.d0
      my = 0.d0
      mz = 0.5d0
      do i=1,Npl
        if(layertype(i+1).eq.2) mz(i) = 2.d0
      end do
      mp = zero
      if(FIELD) then
        mx = mz*sin(hwt*pi)*cos(hwp*pi)
        my = mz*sin(hwt*pi)*sin(hwp*pi)
        mz = mz*cos(hwt*pi)
        mp = cmplx(mx,my,double)
      end if
      sc_solu(1:Npl)         = eps1
      sc_solu(Npl+1:2*Npl)   = mx
      sc_solu(2*Npl+1:3*Npl) = my
      sc_solu(3*Npl+1:4*Npl) = mz
      ! Variables used in the hamiltonian
      do i=1,Npl
        hdel(i)   = 0.5d0*U(i+1)*mz(i)
        hdelp(i)  = 0.5d0*U(i+1)*mp(i)
      end do
      hdelm = conjg(hdelp)
    end select sc_file_status

!     Rotate the magnetization to the direction of the field (useful for SOC=F)
!     do i=1,Npl
!       mabs(i) = sqrt(abs(mp(i))**2+(mz(i)**2))
!       mx(i)   = mabs(i)*sin(hwt*pi)*cos(hwp*pi)
!       my(i)   = mabs(i)*sin(hwt*pi)*sin(hwp*pi)
!       mz(i)   = mz(i)*cos(hwt*pi)
!       mp(i)   = cmplx(mx(i),my(i),double)
!     end do

!     ! Writing new eps1 and rotated mag to file (without self-consistency)
!     if(myrank.eq.0) then
!       iflag = 1
!       call readwritesc(iflag,err)
!     end if

!     call MPI_Finalize(ierr)
!     stop

    ! Self-consistency
    if(selfcon) then
!-------------------------------------------------------------------------------
      iter  = 1
      ifail = 0
      mpitag = (hwa_count-1)*hwp_npt1*hwt_npt1*Npl_f+(hwp_count-1)*hwt_npt1*Npl_f+(hwt_count-1)*Npl_f+Npl
      if(myrank.eq.0) write(*,"('Starting self-consistency:')")

      if(lslatec) then
        if(lnojac) then
          call dnsqe(selfconsistency,selfconsistencyjac,2,neq,sc_solu,fvec,tol,0,ifail,wa,lwa)
        else
          call dnsqe(selfconsistency,selfconsistencyjac,1,neq,sc_solu,fvec,tol,0,ifail,wa,lwa)
        end if
        ifail = ifail-1
      else
        if(lnojac) then
!           call c05nbf(selfconsistency,neq,sc_solu,fvec,tol,wa,lwa,ifail)
          maxfev = 200*(neq+1)
          ml = neq-1
          mr = neq-1
          epsfcn = 1.d-5
          mode = 1
          factor = 100.d0
          call c05ncf(selfconsistency,neq,sc_solu,fvec,tol,maxfev,ml,mr,epsfcn,diag,mode,factor,0,nfev,jac,4*Npl,wa,lwa,qtf,w,ifail)
        else
!           call c05pbf(selfconsistencyjacnag,neq,sc_solu,fvec,jac,neq,tol,wa,lwa,ifail)
          maxfev = 100*(neq+1)
          mode = 1
          diag = 1.d0
!           diag(Npl+1:4*Npl) = 100.d0
          factor = 100.d0
          call c05pcf(selfconsistencyjacnag,neq,sc_solu,fvec,jac,neq,tol,maxfev,diag,mode,factor,0,nfev,njev,wa,lwa,qtf,w,ifail)
        end if
      end if

      if(.not.lontheflysc) then
        ! Writing new eps1 and mz to file after self-consistency is done
        if(myrank.eq.0) then
          iflag = 1
          call readwritesc(iflag,err)
        end if
        call MPI_Bcast(scfile,len(scfile),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      end if
    end if

    if(lGSL) then
      allocate( lxm(Npl),lym(Npl),lzm(Npl),lxpm(Npl),lypm(Npl),lzpm(Npl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        if(myrank.eq.0) write(*,"('[main] Not enough memory for: lxm,lym,lzm,lxpm,lypm,lzpm')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      if(myrank.eq.0) write(*,"('[main] Calculating Orbital Angular Momentum ground state... ')")
      call L_gs()
    end if

    ! Calculating angles of GS magnetization and OAM
    do i = 1,Npl
      mabs(i)   = sqrt((mx(i)**2)+(my(i)**2)+(mz(i)**2))
      mtheta(i) = acos(mz(i)/mabs(i))
      mphi(i)   = atan2(my(i),mx(i))
      if(lGSL) then
        labs(i)   = sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
        ltheta(i) = acos(lzm(i)/sqrt(lxm(i)**2+lym(i)**2+lzm(i)**2))
        lphi(i)   = atan2(lym(i),lxm(i))
        lpabs(i)  = sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2))
        lptheta(i)= acos(lzpm(i)/sqrt(lxpm(i)**2+lypm(i)**2+lzpm(i)**2))
        lpphi(i)  = atan2(lypm(i),lxpm(i))
      end if
    end do
    if(sum(abs(mtheta)).gt.1.d-8) lrot = .true.

    if(myrank.eq.0) then
      write(*,"('|---------------------- Self-consistent ground state: ----------------------|')")
      write(*,"(11x,' *************** Center of d bands: ***************')")
      do i=1,Npl
        write(*,"(26x,'eps1(',i2.0,')=',f11.8)") i,eps1(i)
      end do
      write(*,"(11x,' *********** Magnetization components: **********')")
      do i=1,Npl
        write(*,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mx(i),i,my(i),i,mz(i)
        if(abs(mp(i)).ne.0) write(*,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") mtheta(i)/pi,mphi(i)/pi
      end do
      if(lGSL) then
        write(*,"(11x,' *** Orbital components in spin coordinates:  ***')")
        do i=1,Npl
          write(*,"(4x,'Lxp(',i2.0,')=',f11.8,4x,'Lyp(',i2.0,')=',f11.8,4x,'Lzp(',i2.0,')=',f11.8)") i,lxpm(i),i,lypm(i),i,lzpm(i)
          if(sqrt(lxpm(i)**2+lypm(i)**2).ne.0) write(*,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") lptheta(i)/pi,lpphi(i)/pi
        end do
        write(*,"(11x,' *** Orbital components in cubic coordinates: ***')")
        do i=1,Npl
          write(*,"(4x,'Lx (',i2.0,')=',f11.8,4x,'Ly (',i2.0,')=',f11.8,4x,'Lz (',i2.0,')=',f11.8)") i,lxm(i),i,lym(i),i,lzm(i)
        end do
        write(*,"(11x,' ******************** Total: ********************')")
        do i=1,Npl
          write(*,"(4x,'M (',i2.0,') =',f11.8,4x,'Lp (',i2.0,')=',f11.8,4x,'L (',i2.0,') =',f11.8)") i,mabs(i),i,lpabs(i),i,labs(i)
        end do
      else
        write(*,"(11x,' ******************** Total: ********************')")
        do i=1,Npl
          write(*,"(27x,'M (',i2.0,') =',f11.8)") i,mabs(i)
        end do
      end if
      write(*,"('|---------------------------------------------------------------------------|')")
      elapsed_time = MPI_Wtime() - start_program
      write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
    end if

!============================= MAIN PROGRAM ============================
    main_program: select case (itype)
    case (2)
!----------------------------- Begin test part -----------------------
!       if(myrank.eq.0) then
        if(myrank.eq.0) write(*,"('TESTING')")
        call eintidiag()

!         test2_energy_loop: do count=1,npt1
!           e = emin + (count-1)*deltae
!           write(*,"(i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e*ry2ev


!         end do test2_energy_loop
!       end if
!---------------------------- End test part ----------------------------


!-----------------------------------------------------------------------
    case (3)
      if(myrank.eq.0) then
        allocate(ldosu(Npl,9),ldosd(Npl,9))
        allocate(trJij(nmaglayers,nmaglayers),Jij(nmaglayers,nmaglayers,3,3),Jijs(nmaglayers,nmaglayers,3,3),Jija(nmaglayers,nmaglayers,3,3))

        write(*,"('CALCULATING LDOS AND EXCHANGE INTERACTIONS AS A FUNCTION OF ENERGY')")

        ! Opening files
        ! LDOS
        do i=1,Npl
          iw = 17+(i-1)*2
          write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/LDOS/ldosu_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,Npl,i,magaxis,socscale,ncp,eta,Utype,hwa,hwt,hwp
          open (unit=iw, file=varm,status='unknown')
          iw = iw+1
          write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/LDOS/ldosd_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,Npl,i,magaxis,socscale,ncp,eta,Utype,hwa,hwt,hwp
          open (unit=iw, file=varm,status='unknown')
        end do
        ! Exchange interactions
        do j=1,nmaglayers ; do i=1,nmaglayers
          iw = 99+(j-1)*nmaglayers*2+(i-1)*2
          if(i.eq.j) then
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/Jij/Jii_',i0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,Npl,mmlayermag(i)-1,magaxis,socscale,ncp,eta,Utype,hwa,hwt,hwp
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   energy      ,  Jii_xx           ,   Jii_yy  ')")
            iw = iw + 1
            ! TODO : Check how to write the anisotropy term here
          else
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/Jij/J_',i0,'_',i0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,Npl,mmlayermag(i)-1,mmlayermag(j)-1,magaxis,socscale,ncp,eta,Utype,hwa,hwt,hwp
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   energy      ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/Jij/Dz_',i0,'_',i0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,Npl,mmlayermag(i)-1,mmlayermag(j)-1,magaxis,socscale,ncp,eta,Utype,hwa,hwt,hwp
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   energy      , Dz = (Jxy - Jyx)/2       ')")
          end if
        end do ; end do

        ldos_energy_loop: do count=1,npt1
          e = emin + (count-1)*deltae
          write(*,"(i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e*ry2ev

          call ldos(e,ldosu,ldosd,Jij)

          do i=1,nmaglayers ; do j=1,nmaglayers
            trJij(i,j)    = 0.5d0*(Jij(i,j,1,1)+Jij(i,j,2,2))
            Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
            Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
            do mu=1,3
              Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
            end do
          end do ; end do

          ! Transform energy to eV if runoption is on
          e = e*ry2ev

          ! Writing into files
          ! LDOS
          ldos_writing_plane_loop: do i=1,Npl
              iw = 17+(i-1)*2
              write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
              iw = iw+1
              write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
          end do ldos_writing_plane_loop

          ! Exchange interactions
          jij_writing_loop: do j=1,nmaglayers ; do i=1,nmaglayers
            iw = 99+(j-1)*nmaglayers*2+(i-1)*2
            if(i.eq.j) then
              iw = iw + 1
              write(unit=iw,fmt="(3(es16.9,2x))") e,Jij(i,j,1,1),Jij(i,j,2,2)
              iw = iw + 1
            else
              iw = iw + 1
              write(unit=iw,fmt="(4(es16.9,2x))") e,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
              iw = iw + 1
              write(unit=iw,fmt="(2(es16.9,2x))") e,Jija(i,j,1,2)
            end if
          end do ; end do jij_writing_loop

        end do ldos_energy_loop

        ! Closing files
        do i=1,Npl
          iw = 17+(i-1)*2
          close (iw)
          iw = iw+1
          close (iw)
        end do
        do j=1,nmaglayers ; do i=1,nmaglayers
          iw = 99+(j-1)*nmaglayers*2+(i-1)*2
          if(i.eq.j) then
            iw = iw + 1
            close (iw)
            iw = iw + 1
          else
            iw = iw + 1
            close (iw)
            iw = iw + 1
            close (iw)
          end if
        end do ; end do

        deallocate(ldosu,ldosd)
        deallocate(trJij,Jij,Jijs,Jija)
      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (4)
      if(myrank.eq.0) then
        write(*,"('CALCULATING THE BAND STRUCTURE')")

        call bandstructure()

      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (5)
      call allocate_susceptibilities()
      if(myrank_row.eq.0) then
        allocate( temp(dim,dim), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          write(*,"('[main] Not enough memory for: temp')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        if(myrank_col.eq.0) then
          write(*,"('CALCULATING LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
          write(*,"('Qx = ',es10.3,', Qz = ',es10.3)") q(1),q(2)
          ! Creating files and writing headers
          if(.not.laddresults) then
            call openclose_chi_files(0)
          end if
        end if
      end if

      q     = [0.d0, 0.d0]
      if((myrank.eq.0).and.(skip_steps.gt.0)) write(*,"('[main] Skipping first ',i0,' step(s)...')") skip_steps

      chi_energy_loop: do count=1+skip_steps,MPIsteps
        mpitag = (hwa_count-1)*hwp_npt1*hwt_npt1*npt1+(hwp_count-1)*hwt_npt1*npt1+(hwt_count-1)*npt1+count
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',es10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e*ry2ev,myrank_col
        end if

        ! Start parallelized processes to calculate chiorb_hf and chiorbi0_hf for energy e
        call eintshechi(e)

        if(myrank_row.eq.0) then
          ! (1 + chi_hf*Umat)^-1
          temp = identt
          call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zum,temp,dim)
!           temp = identt + temp
          call invers(temp,dim)
          call zgemm('n','n',dim,dim,dim,zum,temp,dim,chiorb_hf,dim,zero,chiorb,dim)

          schi = zero
          schihf = zero
          ! Calculating RPA and HF susceptibilities
          calculate_susceptibility_chi: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
            schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
            schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
          end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_chi
          schi   = schi/ry2ev
          schihf = schihf/ry2ev
          ! Rotating susceptibilities to the magnetization direction
          if(lrot) then
            do i=1,Npl
              call build_rotation_matrices(mtheta(i),mphi(i),rottemp,1)
              rotmat_i(:,:,i) = rottemp
              call build_rotation_matrices(mtheta(i),mphi(i),rottemp,2)
              rotmat_j(:,:,i) = rottemp
            end do
            rotate_susceptibility_chi: do j=1,Npl ; do i=1,Npl
              rottemp  = rotmat_i(:,:,i)
              schitemp = schi(:,:,i,j)
              call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
              schi(:,:,i,j) = schitemp

              rottemp  = rotmat_i(:,:,i)
              schitemp = schihf(:,:,i,j)
              call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
              schihf(:,:,i,j) = schitemp
            end do ; end do rotate_susceptibility_chi
          end if

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            MPI_points_chi: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,1000,MPIComm_Col,stat,ierr)
                call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1100,MPIComm_Col,stat,ierr)
                call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1200,MPIComm_Col,stat,ierr)
              end if

              ! Transform energy to eV if runoption is on
              e = e*ry2ev

              ! DIAGONALIZING SUSCEPTIBILITY
              call diagonalize_susceptibilities()

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              ! Opening chi and diag files
              call openclose_chi_files(1)
              ! Writing susceptibilities
              call write_susceptibilities(e)
              ! Closing chi and diag files
              call openclose_chi_files(2)
            end do MPI_points_chi

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
              close(911)
              write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
              call system ('rm stop')
              write(*,"(a)") "[main] ('stop' file deleted!)"
              call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
            end if
          else
            call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,1000,MPIComm_Col,ierr)
            call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,1100,MPIComm_Col,ierr)
            call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,1200,MPIComm_Col,ierr)
          end if
        end if
      end do chi_energy_loop

      call deallocate_susceptibilities()
      if(myrank_row.eq.0) then
        deallocate(temp)
      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (6)
      q     = [0.d0, 0.d0]
      call allocate_susceptibilities()
      call allocate_disturbances()
      call allocate_beff()
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          write(*,"('CALCULATING DISTURBANCES AND LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
          write(*,"('Qx = ',es10.3,', Qz = ',es10.3)") q(1),q(2)
          ! Creating files and writing headers
          if(.not.laddresults) then
            call openclose_chi_files(0)
            call openclose_disturbances_files(0)
            call openclose_beff_files(0)
          end if
        end if
      end if

      if((myrank.eq.0).and.(skip_steps.gt.0)) write(*,"('[main] Skipping first ',i0,' step(s)...')") skip_steps

      disturbances_energy_loop: do count=1+skip_steps,MPIsteps
        mpitag = (hwa_count-1)*hwp_npt1*hwt_npt1*npt1+(hwp_count-1)*hwt_npt1*npt1+(hwt_count-1)*npt1+count
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',es10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e*ry2ev,myrank_col
        end if

        if(myrank.eq.0) write(*,"('[main] Calculating prefactor to use in disturbances calculation ')")
        call eintshechi(e)

        ! Broadcast chiorb_hf to all processors of the same row
        call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)

        ! prefactor = (1 + chi_hf*Umat)^-1
        prefactor = identt
        call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zum,prefactor,dim) !prefactor = 1+chi_hf*Umat
        call invers(prefactor,dim)
        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Calculated prefactor after: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        ! Start parallelized processes to calculate disturbances for energy e
        call eintshesd(e)

        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        if(myrank_row.eq.0) then
          ! Calculating the full matrix of RPA and HF susceptibilities for energy e
          call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hf,dim,zero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf

          chiinv = zero
          ! Calculating susceptibility to use on Beff calculation
          calculate_susceptibility_Beff_sd: do nu=1,9 ; do j=1,Npl ; do sigmap=1,4 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4
            chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))    ! +- , up- , down- , --
          end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_Beff_sd

          schi = zero
          schihf = zero
          ! Calculating RPA and HF susceptibilities
          calculate_susceptibility_sd: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
            schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
            schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
          end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_sd
          schi   = schi/ry2ev
          schihf = schihf/ry2ev
          ! Rotating susceptibilities to the magnetization direction
          if(lrot) then
            do i=1,Npl
              call build_rotation_matrices(mtheta(i),mphi(i),rottemp,1)
              rotmat_i(:,:,i) = rottemp
              call build_rotation_matrices(mtheta(i),mphi(i),rottemp,2)
              rotmat_j(:,:,i) = rottemp
            end do
            rotate_susceptibility_sd: do j=1,Npl ; do i=1,Npl
              rottemp  = rotmat_i(:,:,i)
              schitemp = schi(:,:,i,j)
              call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
              schi(:,:,i,j) = schitemp

              rottemp  = rotmat_i(:,:,i)
              schitemp = schihf(:,:,i,j)
              call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
              schihf(:,:,i,j) = schitemp
            end do ; end do rotate_susceptibility_sd
          end if

          disturbances = zero
          calculate_sd: do i=1,Npl
            ! Spin and charge disturbances
            do mu=1,9
              disturbances(1,i) = disturbances(1,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
              disturbances(2,i) = disturbances(2,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
              disturbances(3,i) = disturbances(3,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
              disturbances(4,i) = disturbances(4,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
            end do

            ! Spin disturbance matrix to calculate effective field
            sdmat(sigmai2i(1,i)) = disturbances(2,i) + zi*disturbances(3,i) ! +    = x + iy
            sdmat(sigmai2i(2,i)) = disturbances(1,i) + disturbances(4,i)    ! up   = 0 + z
            sdmat(sigmai2i(3,i)) = disturbances(1,i) - disturbances(4,i)    ! down = 0 - z
            sdmat(sigmai2i(4,i)) = disturbances(2,i) - zi*disturbances(3,i) ! -    = x - iy

            ! Orbital angular momentum disturbance
            do nu=1,9; do mu=1,9
              ldmat(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
              disturbances(5,i) = disturbances(5,i) + lxp(mu,nu)*ldmat(i,mu,nu)
              disturbances(6,i) = disturbances(6,i) + lyp(mu,nu)*ldmat(i,mu,nu)
              disturbances(7,i) = disturbances(7,i) + lzp(mu,nu)*ldmat(i,mu,nu)
            end do; end do
          end do calculate_sd
          disturbances(2:4,:) = 0.5d0*disturbances(2:4,:)
          disturbances(3,:)   = disturbances(3,:)/zi
          sdmat = 0.5d0*sdmat

          ! Effective field calculation
          call invers(chiinv,dimsigmaNpl) ! Inverse of the susceptibility chi^(-1)
          call zgemm('n','n',dimsigmaNpl,1,dimsigmaNpl,zum,chiinv,dimsigmaNpl,sdmat,dimsigmaNpl,zero,Beff,dimsigmaNpl) ! Beff = chi^(-1)*SD
          Beff = Beff*ry2ev

          plane_loop_effective_field_sd: do i=1,Npl
            Beff_cart(sigmai2i(1,i)) =          (Beff(sigmai2i(2,i)) + Beff(sigmai2i(3,i))) ! 0
            Beff_cart(sigmai2i(2,i)) = 0.5d0*   (Beff(sigmai2i(1,i)) + Beff(sigmai2i(4,i))) ! x
            Beff_cart(sigmai2i(3,i)) =-0.5d0*zi*(Beff(sigmai2i(1,i)) - Beff(sigmai2i(4,i))) ! y
            Beff_cart(sigmai2i(4,i)) = 0.5d0*   (Beff(sigmai2i(2,i)) - Beff(sigmai2i(3,i))) ! z
          end do plane_loop_effective_field_sd

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            MPI_points_sd: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,2000,MPIComm_Col,stat,ierr)
                call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2010,MPIComm_Col,stat,ierr)
                call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2020,MPIComm_Col,stat,ierr)
                call MPI_Recv(disturbances,7*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2030,MPIComm_Col,stat,ierr)
                call MPI_Recv(Beff_cart,dimsigmaNpl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2100,MPIComm_Col,stat,ierr)
              end if

              ! Transform energy to eV if runoption is on
              e = e*ry2ev

              ! DIAGONALIZING SUSCEPTIBILITY
              call diagonalize_susceptibilities()

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              ! Opening chi and diag files
              call openclose_chi_files(1)
              ! Writing susceptibilities
              call write_susceptibilities(e)
              ! Closing chi and diag files
              call openclose_chi_files(2)


              ! WRITING DISTURBANCES
              ! Opening disturbance files
              call openclose_disturbances_files(1)
              ! Writing disturbances
              call write_disturbances(e)
              ! Closing disturbance files
              call openclose_disturbances_files(2)

              ! Opening B effective files
              call openclose_beff_files(1)
              ! Writing effective fields
              call write_beff(e)
              ! Closing B effective files
              call openclose_beff_files(2)
            end do MPI_points_sd

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
              close(911)
              write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
              call system ('rm stop')
              write(*,"(a)") "[main] ('stop' file deleted!)"
              call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
            end if
          else
            call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,2000,MPIComm_Col,ierr)
            call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,2010,MPIComm_Col,ierr)
            call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,2020,MPIComm_Col,ierr)
            call MPI_Send(disturbances,7*Npl,MPI_DOUBLE_COMPLEX,0,2030,MPIComm_Col,ierr)
            call MPI_Send(Beff_cart,dimsigmaNpl,MPI_DOUBLE_COMPLEX,0,2100,MPIComm_Col,ierr)
          end if
        end if
      end do disturbances_energy_loop

      call deallocate_susceptibilities()
      call deallocate_disturbances()
      call deallocate_beff()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (7)
      call allocate_susceptibilities()
      call allocate_currents()
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          write(*,"('CALCULATING CURRENTS AND LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
          ! Creating files and writing headers
          if(.not.laddresults) then
            call openclose_chi_files(0)
            call openclose_currents_files(0)
          end if
        end if
      end if

      ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
      call allocate_prefactors()
      call OAM_curr_hopping_times_L()

      if((myrank.eq.0).and.(skip_steps.gt.0)) write(*,"('[main] Skipping first ',i0,' step(s)...')") skip_steps

      currents_energy_loop: do count=1+skip_steps,MPIsteps
        mpitag = (hwa_count-1)*hwp_npt1*hwt_npt1*npt1+(hwp_count-1)*hwt_npt1*npt1+(hwt_count-1)*npt1+count
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',es10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e*ry2ev,myrank_col
        end if

        if(myrank.eq.0) write(*,"('[main] Calculating prefactor to use in current calculation ')")
        call eintshechi(e)

        ! Broadcast chiorb_hf to all processors of the same row
        call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)

        ! prefactor = (1 + chi_hf*Umat)^-1
        prefactor = identt
        call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zum,prefactor,dim) !prefactor = 1+chi_hf*Umat
        call invers(prefactor,dim)
        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Calculated prefactor after: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        ! Start parallelized processes to calculate chiorb_hf and chiorbi0_hf for energy e
        call eintsheprllsc(e)

        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        if(myrank_row.eq.0) then
          ! Calculating the full matrix of RPA and HF susceptibilities for energy e
          call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hf,dim,zero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf

          schi = zero
          schihf = zero
          ! Calculating RPA and HF susceptibilities
          calculate_susceptibility_sc: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
            schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
            schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
          end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_sc
          schi   = schi/ry2ev
          schihf = schihf/ry2ev
          ! Rotating susceptibilities to the magnetization direction
          if(lrot) then
            do i=1,Npl
              call build_rotation_matrices(mtheta(i),mphi(i),rottemp,1)
              rotmat_i(:,:,i) = rottemp
              call build_rotation_matrices(mtheta(i),mphi(i),rottemp,2)
              rotmat_j(:,:,i) = rottemp
            end do
            rotate_susceptibility_sc: do j=1,Npl ; do i=1,Npl
              rottemp  = rotmat_i(:,:,i)
              schitemp = schi(:,:,i,j)
              call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
              schi(:,:,i,j) = schitemp

              rottemp  = rotmat_i(:,:,i)
              schitemp = schihf(:,:,i,j)
              call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
              schihf(:,:,i,j) = schitemp
            end do ; end do rotate_susceptibility_sc
          end if

          ! Calculating spin and charge current for each neighbor
          plane_loop_calculate_sc: do i=1,Npl
            neighbor_loop_calculate_sc: do neighbor=n0sc1,n0sc2
              ! Charge current
              currents(1,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)+ttchiorbiikl  (neighbor,sigmai2i(3,i),2)+ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
              ! Spin currents
              currents(2,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)+ttchiorbiikl  (neighbor,sigmai2i(4,i),2)+ttchiorbiikl  (neighbor,sigmai2i(4,i),3)
              currents(3,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)-ttchiorbiikl  (neighbor,sigmai2i(4,i),2)-ttchiorbiikl  (neighbor,sigmai2i(4,i),3)
              currents(4,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)-ttchiorbiikl  (neighbor,sigmai2i(3,i),2)-ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
              ! Orbital Angular Momentum currents
              currents(5,neighbor,i) = Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3)
              currents(6,neighbor,i) = Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3)
              currents(7,neighbor,i) = Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3)
            end do neighbor_loop_calculate_sc
          end do plane_loop_calculate_sc
          currents = currents*ry2ev
          currents(2:4,:,:) = -0.5d0*currents(2:4,:,:)
          currents(3,:,:) = currents(3,:,:)/zi
          ! Total currents for each neighbor direction (Sum of currents over all planes)
          total_currents = sum(currents,dim=3)

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            MPI_points_sc: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3000,MPIComm_Col,stat,ierr)
                call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3100,MPIComm_Col,stat,ierr)
                call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3200,MPIComm_Col,stat,ierr)
                call MPI_Recv(currents,7*n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3300,MPIComm_Col,stat,ierr)
                call MPI_Recv(total_currents,7*n0sc,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3400,MPIComm_Col,stat,ierr)
              end if

              ! Transform energy to eV if runoption is on
              e = e*ry2ev

              ! DIAGONALIZING SUSCEPTIBILITY
              call diagonalize_susceptibilities()

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              ! Opening chi and diag files
              call openclose_chi_files(1)
              ! Writing susceptibilities
              call write_susceptibilities(e)
              ! Closing chi and diag files
              call openclose_chi_files(2)

              ! Renormalizing disturbances and currents by the total charge current to neighbor renormnb
              if(renorm) then
                ! Obtaining current for renormalization
                Icabs  = abs(total_currents(1,renormnb))
                rcurrents = currents/Icabs
                rtotal_currents = total_currents/Icabs
              end if

              ! WRITING CURRENTS
              ! Opening current files
              call openclose_currents_files(1)
              ! Writing currents
              call write_currents(e)
              ! Closing current files
              call openclose_currents_files(2)
            end do MPI_points_sc

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
              close(911)
              write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
              call system ('rm stop')
              write(*,"(a)") "[main] ('stop' file deleted!)"
              call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
            end if
          else
            call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,3000,MPIComm_Col,ierr)
            call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,3100,MPIComm_Col,ierr)
            call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,3200,MPIComm_Col,ierr)
            call MPI_Send(currents,7*n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3300,MPIComm_Col,ierr)
            call MPI_Send(total_currents,7*n0sc,MPI_DOUBLE_COMPLEX,0,3400,MPIComm_Col,ierr)
          end if
        end if
      end do currents_energy_loop

      call deallocate_susceptibilities()
      call deallocate_currents()
      call deallocate_prefactors()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (8)
      call allocate_susceptibilities()
      call allocate_disturbances()
      call allocate_currents()
      call allocate_beff()
      call allocate_torques()
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          write(*,"('CALCULATING PARALLEL CURRENTS, DISTURBANCES, LOCAL SUSCEPTIBILITY, EFFECTIVE FIELDS AND SO-TORQUES AS A FUNCTION OF ENERGY')")
          ! Creating files and writing headers
          if(.not.laddresults) then
            if(.not.lhfresponses) call openclose_chi_files(0)
            call openclose_disturbances_files(0)
            call openclose_currents_files(0)
            call openclose_beff_files(0)
            call openclose_torques_files(0)
          end if
        end if
      end if

      ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
      call allocate_prefactors()
      call OAM_curr_hopping_times_L()

      if((myrank.eq.0).and.(skip_steps.gt.0)) write(*,"('[main] Skipping first ',i0,' step(s)...')") skip_steps

      all_energy_loop: do count=1+skip_steps,MPIsteps
        mpitag = (hwa_count-1)*hwp_npt1*hwt_npt1*npt1+(hwp_count-1)*hwt_npt1*npt1+(hwt_count-1)*npt1+count
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',es10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e*ry2ev,myrank_col
        end if

        if(lhfresponses) then
          if(myrank.eq.0) write(*,"('[main] No renormalization will be done. Setting prefactors to identity. ')")
          prefactor     = identt
          if(llinearsoc) prefactorlsoc = identt
        else
          if(myrank.eq.0) write(*,"('[main] Calculating prefactor to use in currents and disturbances calculation. ')")
          if(llinearsoc) then
            call eintshechilinearsoc(e) ! Note: chiorb_hflsoc = lambda*dchi_hf/dlambda(lambda=0)
            ! Broadcast chiorb_hflsoc to all processors of the same row
            call MPI_Bcast(chiorb_hflsoc,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)
          else
            call eintshechi(e)
          end if

          ! Broadcast chiorb_hf to all processors of the same row
          call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)

          ! prefactor = (1 + chi_hf*Umat)^-1
          prefactor     = identt
          call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zum,prefactor,dim) ! prefactor = 1+chi_hf*Umat
          call invers(prefactor,dim)
          if(llinearsoc) then
            prefactorlsoc = identt
            if (.not. allocated(chiorb)) allocate(chiorb(dim,dim))
            chiorb = chiorb_hf-chiorb_hflsoc ! the array chiorb will be used as a temporary array
            call zgemm('n','n',dim,dim,dim,zum,chiorb,dim,Umatorb,dim,zum,prefactorlsoc,dim) ! prefactorlsoc = 1+chiorb*Umat = 1+(chi_hf + dchi_hf/dlambda)*Umat
            call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,prefactorlsoc,dim,zero,chiorb,dim) ! chiorb = prefactor*prefactorlsoc
            call zgemm('n','n',dim,dim,dim,zum,chiorb,dim,prefactor,dim,zero,prefactorlsoc,dim) ! prefactorlsoc = chiorb*prefactor = prefactor*prefactorlsoc*prefactor
          end if
          if(myrank.eq.0) then
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Calculated prefactor after: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
          end if
        end if

        ! Start parallelized processes to calculate disturbances and currents for energy e
        call eintshe(e)

        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        if(myrank_row.eq.0) then

          if(.not.lhfresponses) then
            ! Calculating the full matrix of RPA and HF susceptibilities for energy e
            if(llinearsoc) then
              ! chiorb = prefactorlsoc*chi_hf + prefactor*chi_hflsoc
              call zgemm('n','n',dim,dim,dim,zum,prefactorlsoc,dim,chiorb_hf,dim,zero,chiorb,dim)
              call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hflsoc,dim,zum,chiorb,dim)
            else
              ! chiorb = prefactor*chi_hf
              call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hf,dim,zero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf
            end if

            chiinv = zero
            ! Calculating susceptibility to use on Beff calculation
            calculate_susceptibility_Beff_all: do nu=1,9 ; do j=1,Npl ; do sigmap=1,4 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4
              chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))    ! +- , up- , down- , --
            end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_Beff_all

            schi = zero
            schihf = zero
            ! Calculating RPA and HF susceptibilities
            calculate_susceptibility_all: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
              schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
              schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
            end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_all
            schi   = schi/ry2ev
            schihf = schihf/ry2ev
            ! Rotating susceptibilities to the magnetization direction
            if(lrot) then
              do i=1,Npl
                call build_rotation_matrices(mtheta(i),mphi(i),rottemp,1)
                rotmat_i(:,:,i) = rottemp
                call build_rotation_matrices(mtheta(i),mphi(i),rottemp,2)
                rotmat_j(:,:,i) = rottemp
              end do
              rotate_susceptibility_all: do j=1,Npl ; do i=1,Npl
                rottemp  = rotmat_i(:,:,i)
                schitemp = schi(:,:,i,j)
                call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
                rottemp  = rotmat_j(:,:,j)
                call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
                schi(:,:,i,j) = schitemp

                rottemp  = rotmat_i(:,:,i)
                schitemp = schihf(:,:,i,j)
                call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
                rottemp  = rotmat_j(:,:,j)
                call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
                schihf(:,:,i,j) = schitemp
              end do ; end do rotate_susceptibility_all
            end if
          end if ! lhfresponses

          disturbances = zero
          torques      = zero
          plane_loop_calculate_all: do i=1,Npl
            ! Spin and charge disturbances
            do mu=1,9
              disturbances(1,i) = disturbances(1,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
              disturbances(2,i) = disturbances(2,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
              disturbances(3,i) = disturbances(3,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))/zi
              disturbances(4,i) = disturbances(4,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
            end do

            ! Spin disturbance matrix to calculate effective field
            sdmat(sigmai2i(1,i)) = disturbances(2,i) + zi*disturbances(3,i) ! +    = x + iy
            sdmat(sigmai2i(2,i)) = disturbances(1,i) + disturbances(4,i)    ! up   = 0 + z
            sdmat(sigmai2i(3,i)) = disturbances(1,i) - disturbances(4,i)    ! down = 0 - z
            sdmat(sigmai2i(4,i)) = disturbances(2,i) - zi*disturbances(3,i) ! -    = x - iy

            ! Orbital angular momentum disturbance
            do nu=1,9; do mu=1,9
              ldmat(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
              disturbances(5,i) = disturbances(5,i) + lxp(mu,nu)*ldmat(i,mu,nu)
              disturbances(6,i) = disturbances(6,i) + lyp(mu,nu)*ldmat(i,mu,nu)
              disturbances(7,i) = disturbances(7,i) + lzp(mu,nu)*ldmat(i,mu,nu)
            end do; end do

            ! Spin-orbit torques
            do nu=1,9; do mu=1,9
              ! x component: Ly*Sz - Lz*Sy
              torques(1,1,i) = torques(1,1,i) + (   lyp(mu,nu)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3))) &
                                          + (zi*lzp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
              ! y component: Lz*Sx - Lx*Sz
              torques(1,2,i) = torques(1,2,i) + (   lzp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                          - (   lxp(mu,nu)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)))
              ! z component: Lx*Sy - Ly*Sx
              torques(1,3,i) = torques(1,3,i) - (zi*lxp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                          - (   lyp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
            end do; end do
            torques(1,:,i) = lambda(i+1)*torques(1,:,i)

            ! Exchange-correlation torques
            torques(2,1,i) = U(i+1)*(mz(i)*disturbances(3,i)-my(i)*disturbances(4,i))
            torques(2,2,i) = U(i+1)*(mx(i)*disturbances(4,i)-mz(i)*disturbances(2,i))
            torques(2,3,i) = U(i+1)*(my(i)*disturbances(2,i)-mx(i)*disturbances(3,i))

            ! Calculating spin and charge current for each neighbor
            neighbor_loop_calculate_all: do neighbor=n0sc1,n0sc2
              ! Charge current
              currents(1,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)+ttchiorbiikl  (neighbor,sigmai2i(3,i),2)+ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
              ! Spin currents
              currents(2,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)+ttchiorbiikl  (neighbor,sigmai2i(4,i),2)+ttchiorbiikl  (neighbor,sigmai2i(4,i),3)
              currents(3,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)-ttchiorbiikl  (neighbor,sigmai2i(4,i),2)-ttchiorbiikl  (neighbor,sigmai2i(4,i),3)
              currents(4,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)-ttchiorbiikl  (neighbor,sigmai2i(3,i),2)-ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
              ! Orbital Angular Momentum currents
              currents(5,neighbor,i) = Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3)
              currents(6,neighbor,i) = Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3)
              currents(7,neighbor,i) = Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3)
            end do neighbor_loop_calculate_all
          end do plane_loop_calculate_all
          disturbances(2:4,:) = 0.5d0*disturbances(2:4,:)
          torques  = 0.5d0*torques
          sdmat    = 0.5d0*sdmat
          currents = currents*ry2ev
          currents(2:4,:,:) = -0.5d0*currents(2:4,:,:)
          currents(3,:,:)   = currents(3,:,:)/zi
          ! Total currents for each neighbor direction (Sum of currents over all planes)
          total_currents = sum(currents,dim=3)

          if(.not.lhfresponses) then
            ! Effective field calculation
            call invers(chiinv,dimsigmaNpl) ! Inverse of the susceptibility chi^(-1)
            call zgemm('n','n',dimsigmaNpl,1,dimsigmaNpl,zum,chiinv,dimsigmaNpl,sdmat,dimsigmaNpl,zero,Beff,dimsigmaNpl) ! Beff = chi^(-1)*SD
            Beff = Beff*ry2ev

            plane_loop_effective_field_all: do i=1,Npl
              Beff_cart(sigmai2i(1,i)) =          (Beff(sigmai2i(2,i)) + Beff(sigmai2i(3,i))) ! 0
              Beff_cart(sigmai2i(2,i)) = 0.5d0*   (Beff(sigmai2i(1,i)) + Beff(sigmai2i(4,i))) ! x
              Beff_cart(sigmai2i(3,i)) =-0.5d0*zi*(Beff(sigmai2i(1,i)) - Beff(sigmai2i(4,i))) ! y
              Beff_cart(sigmai2i(4,i)) = 0.5d0*   (Beff(sigmai2i(2,i)) - Beff(sigmai2i(3,i))) ! z
            end do plane_loop_effective_field_all
          end if ! lhfresponses

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            MPI_points_all: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,4000,MPIComm_Col,stat,ierr)
                if(.not.lhfresponses) then
                  call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4010,MPIComm_Col,stat,ierr)
                  call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4020,MPIComm_Col,stat,ierr)
                  call MPI_Recv(Beff_cart,dimsigmaNpl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4170,MPIComm_Col,stat,ierr)
                end if
                call MPI_Recv(disturbances,7*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4030,MPIComm_Col,stat,ierr)
                call MPI_Recv(currents,7*n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4100,MPIComm_Col,stat,ierr)
                call MPI_Recv(total_currents,7*n0sc,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4110,MPIComm_Col,stat,ierr)
                call MPI_Recv(torques,6*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4120,MPIComm_Col,stat,ierr)
              end if

              ! Transform energy to eV if runoption is on
              e = e*ry2ev

              if(.not.lhfresponses) then
                ! DIAGONALIZING SUSCEPTIBILITY
                call diagonalize_susceptibilities()

                ! WRITING RPA AND HF SUSCEPTIBILITIES
                ! Opening chi and diag files
                call openclose_chi_files(1)
                ! Writing susceptibilities
                call write_susceptibilities(e)
                ! Closing chi and diag files
                call openclose_chi_files(2)
              end if

              ! Renormalizing disturbances and currents by the total charge current to neighbor renormnb
              if(renorm) then
                ! Obtaining current for renormalization
                Icabs  = abs(total_currents(1,renormnb))

                rdisturbances   = disturbances/Icabs
                rtorques        = torques/Icabs
                rcurrents       = currents/Icabs
                rtotal_currents = total_currents/Icabs
              end if

              ! WRITING DISTURBANCES
              ! Opening disturbance files
              call openclose_disturbances_files(1)
              ! Writing disturbances
              call write_disturbances(e)
              ! Closing disturbance files
              call openclose_disturbances_files(2)

              ! WRITING CURRENTS
              ! Opening current files
              call openclose_currents_files(1)
              ! Writing currents
              call write_currents(e)
              ! Closing current files
              call openclose_currents_files(2)

              if(.not.lhfresponses) then
                ! Opening B effective files
                call openclose_beff_files(1)
                ! Writing effective fields
                call write_beff(e)
                ! Closing B effective files
                call openclose_beff_files(2)
              end if

              ! WRITING TORQUES
              ! Opening torque files
              call openclose_torques_files(1)
              ! Writing torques
              call write_torques(e)
              ! Closing torque files
              call openclose_torques_files(2)

            end do MPI_points_all

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
              close(911)
              write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
              call system ('rm stop')
              write(*,"(a)") "[main] ('stop' file deleted!)"
              call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
            end if
          else
            call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,4000,MPIComm_Col,ierr)
            if(.not.lhfresponses) then
              call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,4010,MPIComm_Col,ierr)
              call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,4020,MPIComm_Col,ierr)
              call MPI_Send(Beff_cart,dimsigmaNpl,MPI_DOUBLE_COMPLEX,0,4170,MPIComm_Col,ierr)
            end if
            call MPI_Send(disturbances,7*Npl,MPI_DOUBLE_COMPLEX,0,4030,MPIComm_Col,ierr)
            call MPI_Send(currents,7*n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4100,MPIComm_Col,ierr)
            call MPI_Send(total_currents,7*n0sc,MPI_DOUBLE_COMPLEX,0,4110,MPIComm_Col,ierr)
            call MPI_Send(torques,6*Npl,MPI_DOUBLE_COMPLEX,0,4120,MPIComm_Col,ierr)
          end if
        end if
      end do all_energy_loop

      call deallocate_prefactors()
      call deallocate_susceptibilities()
      call deallocate_disturbances()
      call deallocate_currents()
      call deallocate_beff()
      call deallocate_torques()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (9)
      if(myrank.eq.0) then
        write(*,"('CALCULATING FERMI SURFACE')")

        call fermisurface(Ef)

      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (10)
      allocate(Jij(nmaglayers,nmaglayers,3,3),trJij(nmaglayers,nmaglayers),Jijs(nmaglayers,nmaglayers,3,3),Jija(nmaglayers,nmaglayers,3,3))

      if(myrank.eq.0) write(*,"('CALCULATING FULL TENSOR OF EXHANGE INTERACTIONS AND ANISOTROPIES AS A FUNCTION OF POSITION')")

      if(nmaglayers.eq.0) then
        if(myrank.eq.0) write(*,"(1x,'[main] No magnetic layers!')")
        call MPI_Finalize(ierr)
        stop
      end if

      ! Opening files for position dependence
      if((myrank.eq.0).and.(Npl.eq.Npl_i)) then
        ! Exchange interactions
        do j=1,nmaglayers ; do i=1,nmaglayers
          iw = 199+(j-1)*nmaglayers*2+(i-1)*2
          if(i.eq.j) then
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Jij/Jii_',i0,'_magaxis=',A,'_socscale=',f5.2,'_parts=',I0,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,i,magaxis,socscale,parts,ncp,eta,Utype,hwa,hwt,hwp
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#  Npl ,  Jii_xx           ,   Jii_yy  ')")
            iw = iw + 1
            ! TODO : Check how to write the anisotropy term here
          else
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Jij/J_',i0,'_',i0,'_magaxis=',A,'_socscale=',f5.2,'_parts=',I0,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,i,j,magaxis,socscale,parts,ncp,eta,Utype,hwa,hwt,hwp
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#  Npl ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Jij/Dz_',i0,'_',i0,'_magaxis=',A,'_socscale=',f5.2,'_parts=',I0,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,i,j,magaxis,socscale,parts,ncp,eta,Utype,hwa,hwt,hwp
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#  Npl , Dz = (Jxy - Jyx)/2       ')")
          end if
        end do ; end do
      end if

      call coupling(Jij)

      if(myrank.eq.0) then
        Jij = Jij*ry2ev

        do i=1,nmaglayers ; do j=1,nmaglayers
          trJij(i,j)    = 0.5d0*(Jij(i,j,1,1)+Jij(i,j,2,2))
          Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
          Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
          do mu=1,3
            Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
          end do
        end do ; end do

        ! Writing exchange couplings and anisotropies
        write(*,"('  ************************* Full tensor Jij:  *************************')")
        do i=1,nmaglayers ; do j=1,nmaglayers
        ! Writing on screen
        ! Writing original full tensor Jij
          if(i.eq.j) then
            write(*,"(3x,' ******** Magnetization components: (magaxis = ',a,') *******')") magaxis
            write(*,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mx(i),i,my(i),i,mz(i)
            write(*,"(' |--------------- i = ',i0,'   j = ',i0,': anisotropies ---------------|')") mmlayermag(i),mmlayermag(j)
          else
            write(*,"(' |----------- i = ',i0,'   j = ',i0,': exchange couplings -------------|')") mmlayermag(i),mmlayermag(j)
          end if
          write(*,"('             x                  y                  z')")
          write(*,"('  x  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,1,1),Jij(i,j,1,2),Jij(i,j,1,3)
          write(*,"('  y  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,2,1),Jij(i,j,2,2),Jij(i,j,2,3)
          write(*,"('  z  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,3,1),Jij(i,j,3,2),Jij(i,j,3,3)
        end do ; end do
        if(nmaglayers.gt.1) write(*,"('  *** Symmetric and antisymmetric exchange interactions:  ***')")
        do i=1,nmaglayers ; do j=1,nmaglayers
          if(i.eq.j) cycle
          write(*,"(' |--------------------- i = ',i0,'   j = ',i0,' -----------------------|')") mmlayermag(i),mmlayermag(j)
        ! Writing Heisenberg exchange interactions
          write(*,"('     Isotropic:     J     = ',es16.9)") trJij(i,j)
          write(*,"('   Anisotropic:     Js_xx = ',es16.9)") Jijs(i,j,1,1)
          write(*,"('                    Js_yy = ',es16.9)") Jijs(i,j,2,2)
          write(*,"('  DMI: Dz = (Jxy - Jyx)/2 = ',es16.9)") Jija(i,j,1,2)
          write(*,"(' --- z components of Jij (not physically correct) ---')")
          write(*,"('  Anisotropic:  Js_zz = ',es16.9)") Jijs(i,j,3,3)
          write(*,"('  DMI: Dy = (Jzx - Jxz)/2 = ',es16.9)") -Jija(i,j,1,3)
          write(*,"('  DMI: Dx = (Jyz - Jzy)/2 = ',es16.9)") Jija(i,j,2,3)
        end do ; end do

        ! Writing into files
        ! Exchange interactions
        exchange_writing_loop: do j=1,nmaglayers ; do i=1,nmaglayers
          iw = 199+(j-1)*nmaglayers*2+(i-1)*2
          if(i.eq.j) then
            iw = iw + 1
            write(unit=iw,fmt="(4x,i3,13x,2(es16.9,2x))") Npl,Jij(i,j,1,1),Jij(i,j,2,2)
            iw = iw + 1
          else
            iw = iw + 1
            write(unit=iw,fmt="(4x,i3,13x,3(es16.9,2x))") Npl,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
            iw = iw + 1
            write(unit=iw,fmt="(4x,i3,13x,es16.9,2x)") Npl,Jija(i,j,1,2)
          end if
        end do ; end do exchange_writing_loop

        ! Closing files
        if(Npl.eq.Npl_f) then
          ! Closing files
          do j=1,nmaglayers ; do i=1,nmaglayers
            iw = 199+(j-1)*nmaglayers*2+(i-1)*2
            if(i.eq.j) then
              iw = iw + 1
              close (iw)
              iw = iw + 1
            else
              iw = iw + 1
              close (iw)
              iw = iw + 1
              close (iw)
            end if
          end do ; end do
        end if
      end if

      deallocate(trJij,Jij,Jijs,Jija)
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (11)
      FIELD = .true. ! To use the correct file when hwa=0 (otherwise it does not write the "fieldpart" of the filename)
      call allocate_susceptibilities()
      call allocate_disturbances()
      call allocate_currents()
      call allocate_beff()
      call allocate_torques()
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          write(*,"('CALCULATING DC-LIMIT QUANTITIES AS A FUNCTION OF EXTERNAL FIELD (',a,')')") dcfield(dcfield_dependence)
          ! Creating files and writing headers
          if((.not.laddresults).and.(hwa_count*hwt_count*hwp_count.eq.1)) then
            call openclose_dclimit_disturbances_files(0)
            call openclose_dclimit_currents_files(0)
            call openclose_dclimit_beff_files(0)
            call openclose_dclimit_torques_files(0)
          end if
        end if
      end if

      ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
      call allocate_prefactors()
      call OAM_curr_hopping_times_L()

      e = emin

      if(lhfresponses) then
        if(myrank.eq.0) write(*,"('[main] No renormalization will be done. Setting prefactors to identity. ')")
        prefactor     = identt
        if(llinearsoc) prefactorlsoc = identt
      else
        if(myrank.eq.0) write(*,"('[main] Calculating prefactor to use in currents and disturbances calculation. ')")
        if(llinearsoc) then
          call eintshechilinearsoc(e) ! Note: chiorb_hflsoc = lambda*dchi_hf/dlambda(lambda=0)
          ! Broadcast chiorb_hflsoc to all processors of the same row
          call MPI_Bcast(chiorb_hflsoc,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)
        else
          call eintshechi(e)
        end if

        ! Broadcast chiorb_hf to all processors of the same row
        call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)

        ! prefactor = (1 + chi_hf*Umat)^-1
        prefactor     = identt
        call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zum,prefactor,dim) ! prefactor = 1+chi_hf*Umat
        call invers(prefactor,dim)
        if(llinearsoc) then
          prefactorlsoc = identt
          if (.not. allocated(chiorb)) allocate(chiorb(dim,dim))
          chiorb = chiorb_hf-chiorb_hflsoc ! the array chiorb will be used as a temporary array
          call zgemm('n','n',dim,dim,dim,zum,chiorb,dim,Umatorb,dim,zum,prefactorlsoc,dim) ! prefactorlsoc = 1+chiorb*Umat = 1+(chi_hf + dchi_hf/dlambda)*Umat
          call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,prefactorlsoc,dim,zero,chiorb,dim) ! chiorb = prefactor*prefactorlsoc
          call zgemm('n','n',dim,dim,dim,zum,chiorb,dim,prefactor,dim,zero,prefactorlsoc,dim) ! prefactorlsoc = chiorb*prefactor = prefactor*prefactorlsoc*prefactor
        end if
        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Calculated prefactor after: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if
      end if

      ! Start parallelized processes to calculate disturbances and currents for energy e
      call eintshe(e)

      if(myrank.eq.0) then
        elapsed_time = MPI_Wtime() - start_program
        write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
      end if

      if(myrank_row.eq.0) then

        if(.not.lhfresponses) then
          ! Calculating the full matrix of RPA and HF susceptibilities for energy e
          if(llinearsoc) then
            ! chiorb = prefactorlsoc*chi_hf + prefactor*chi_hflsoc
            call zgemm('n','n',dim,dim,dim,zum,prefactorlsoc,dim,chiorb_hf,dim,zero,chiorb,dim)
            call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hflsoc,dim,zum,chiorb,dim)
          else
            ! chiorb = prefactor*chi_hf
            call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hf,dim,zero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf
          end if

          chiinv = zero
          ! Calculating susceptibility to use on Beff calculation
          calculate_susceptibility_Beff_dclimit: do nu=1,9 ; do j=1,Npl ; do sigmap=1,4 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4
            chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))    ! +- , up- , down- , --
          end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_Beff_dclimit

!           schi = zero
!           schihf = zero
!           ! Calculating RPA and HF susceptibilities
!           calculate_susceptibility_dclimit: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
!             schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
!             schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
!           end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_dclimit
!           schi   = schi/ry2ev
!           schihf = schihf/ry2ev
!           ! Rotating susceptibilities to the magnetization direction
!           if(lrot) then
!             do i=1,Npl
!               call build_rotation_matrices(mtheta(i),mphi(i),rottemp,1)
!               rotmat_i(:,:,i) = rottemp
!               call build_rotation_matrices(mtheta(i),mphi(i),rottemp,2)
!               rotmat_j(:,:,i) = rottemp
!             end do
!             rotate_susceptibility_dclimit: do j=1,Npl ; do i=1,Npl
!               rottemp  = rotmat_i(:,:,i)
!               schitemp = schi(:,:,i,j)
!               call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
!               rottemp  = rotmat_j(:,:,j)
!               call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
!               schi(:,:,i,j) = schitemp

!               rottemp  = rotmat_i(:,:,i)
!               schitemp = schihf(:,:,i,j)
!               call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
!               rottemp  = rotmat_j(:,:,j)
!               call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
!               schihf(:,:,i,j) = schitemp
!             end do ; end do rotate_susceptibility_dclimit
!           end if
        end if ! lhfresponses

        disturbances = zero
        torques      = zero
        plane_loop_calculate_dclimit: do i=1,Npl
          ! Spin and charge disturbances
          do mu=1,9
            disturbances(1,i) = disturbances(1,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
            disturbances(2,i) = disturbances(2,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
            disturbances(3,i) = disturbances(3,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))/zi
            disturbances(4,i) = disturbances(4,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
          end do

          ! Spin disturbance matrix to calculate effective field
          sdmat(sigmai2i(1,i)) = disturbances(2,i) + zi*disturbances(3,i) ! +    = x + iy
          sdmat(sigmai2i(2,i)) = disturbances(1,i) + disturbances(4,i)    ! up   = 0 + z
          sdmat(sigmai2i(3,i)) = disturbances(1,i) - disturbances(4,i)    ! down = 0 - z
          sdmat(sigmai2i(4,i)) = disturbances(2,i) - zi*disturbances(3,i) ! -    = x - iy

          ! Orbital angular momentum disturbance
          do nu=1,9; do mu=1,9
            ldmat(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
            disturbances(5,i) = disturbances(5,i) + lxp(mu,nu)*ldmat(i,mu,nu)
            disturbances(6,i) = disturbances(6,i) + lyp(mu,nu)*ldmat(i,mu,nu)
            disturbances(7,i) = disturbances(7,i) + lzp(mu,nu)*ldmat(i,mu,nu)
          end do; end do

          ! Spin-orbit torques
          do nu=1,9; do mu=1,9
            ! x component: Ly*Sz - Lz*Sy
            torques(1,1,i) = torques(1,1,i) + (   lyp(mu,nu)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3))) &
                                        + (zi*lzp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
            ! y component: Lz*Sx - Lx*Sz
            torques(1,2,i) = torques(1,2,i) + (   lzp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                        - (   lxp(mu,nu)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)))
            ! z component: Lx*Sy - Ly*Sx
            torques(1,3,i) = torques(1,3,i) - (zi*lxp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                        - (   lyp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
          end do; end do
          torques(1,:,i) = lambda(i+1)*torques(1,:,i)

          ! Exchange-correlation torques
          torques(2,1,i) = U(i+1)*(mz(i)*disturbances(3,i)-my(i)*disturbances(4,i))
          torques(2,2,i) = U(i+1)*(mx(i)*disturbances(4,i)-mz(i)*disturbances(2,i))
          torques(2,3,i) = U(i+1)*(my(i)*disturbances(2,i)-mx(i)*disturbances(3,i))

          ! Calculating spin and charge current for each neighbor
          neighbor_loop_calculate_dclimit: do neighbor=n0sc1,n0sc2
            ! Charge current
            currents(1,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)+ttchiorbiikl  (neighbor,sigmai2i(3,i),2)+ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
            ! Spin currents
            currents(2,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)+ttchiorbiikl  (neighbor,sigmai2i(4,i),2)+ttchiorbiikl  (neighbor,sigmai2i(4,i),3)
            currents(3,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)-ttchiorbiikl  (neighbor,sigmai2i(4,i),2)-ttchiorbiikl  (neighbor,sigmai2i(4,i),3)/zi
            currents(4,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)-ttchiorbiikl  (neighbor,sigmai2i(3,i),2)-ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
            ! Orbital Angular Momentum currents
            currents(5,neighbor,i) = Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3)
            currents(6,neighbor,i) = Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3)
            currents(7,neighbor,i) = Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3)
          end do neighbor_loop_calculate_dclimit
        end do plane_loop_calculate_dclimit
        disturbances(2:4,:) = 0.5d0*disturbances(2:4,:)
        disturbances = disturbances/e
        torques  = 0.5d0*torques/e
        sdmat    = 0.5d0*sdmat/e
        currents = currents/e
        currents(2:4,:,:) = -0.5d0*currents(2:4,:,:)
        ! Total currents for each neighbor direction (Sum of currents over all planes)
        total_currents = sum(currents,dim=3)

        if(.not.lhfresponses) then
          ! Effective field calculation
          call invers(chiinv,dimsigmaNpl) ! Inverse of the susceptibility chi^(-1)
          call zgemm('n','n',dimsigmaNpl,1,dimsigmaNpl,zum,chiinv,dimsigmaNpl,sdmat,dimsigmaNpl,zero,Beff,dimsigmaNpl) ! Beff = chi^(-1)*SD
          Beff = Beff*ry2ev

          plane_loop_effective_field_dclimit: do i=1,Npl
            Beff_cart(sigmai2i(1,i)) =       (Beff(sigmai2i(2,i)) + Beff(sigmai2i(3,i)))    ! 0
            Beff_cart(sigmai2i(2,i)) = 0.5d0*(Beff(sigmai2i(1,i)) + Beff(sigmai2i(4,i)))    ! x
            Beff_cart(sigmai2i(3,i)) = 0.5d0*(Beff(sigmai2i(1,i)) - Beff(sigmai2i(4,i)))/zi ! y
            Beff_cart(sigmai2i(4,i)) = 0.5d0*(Beff(sigmai2i(2,i)) - Beff(sigmai2i(3,i)))    ! z
          end do plane_loop_effective_field_dclimit
        end if ! lhfresponses

        ! Sending results to myrank_row = myrank_col = 0 and writing on file
        if(myrank_col.eq.0) then
          MPI_points_dclimit: do mcount=1,MPIpts
            if (mcount.ne.1) then
              call MPI_Recv(dclimit_field_dependence,3,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,44000,MPIComm_Col,stat,ierr)
              if(.not.lhfresponses) then
!                 call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),44010,MPIComm_Col,stat,ierr)
!                 call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),44020,MPIComm_Col,stat,ierr)
                call MPI_Recv(Beff_cart,dimsigmaNpl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),44170,MPIComm_Col,stat,ierr)
              end if
              call MPI_Recv(disturbances,7*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),44030,MPIComm_Col,stat,ierr)
              call MPI_Recv(currents,7*n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),44100,MPIComm_Col,stat,ierr)
              call MPI_Recv(total_currents,7*n0sc,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),44110,MPIComm_Col,stat,ierr)
              call MPI_Recv(torques,6*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),44120,MPIComm_Col,stat,ierr)
            end if

!             if(.not.lhfresponses) then
!               ! DIAGONALIZING SUSCEPTIBILITY
!               call diagonalize_susceptibilities()

!               ! WRITING RPA AND HF SUSCEPTIBILITIES
!               ! Opening chi and diag files
!               call openclose_chi_files(1)
!               ! Writing susceptibilities
!               call write_susceptibilities(e)
!               ! Closing chi and diag files
!               call openclose_chi_files(2)
!             end if

            ! Renormalizing disturbances and currents by the total charge current to neighbor renormnb
            if(renorm) then
              ! Obtaining current for renormalization
              Icabs  = abs(total_currents(1,renormnb))

              rdisturbances   = disturbances/Icabs
              rtorques        = torques/Icabs
              rcurrents       = currents/Icabs
              rtotal_currents = total_currents/Icabs
            end if

            ! WRITING DISTURBANCES
            ! Opening disturbance files
            call openclose_dclimit_disturbances_files(1)
            ! Writing disturbances
            call write_dclimit_disturbances(dclimit_field_dependence(dcfield_dependence))
            ! Closing disturbance files
            call openclose_dclimit_disturbances_files(2)

            ! WRITING CURRENTS
            ! Opening current files
            call openclose_dclimit_currents_files(1)
            ! Writing currents
            call write_dclimit_currents(dclimit_field_dependence(dcfield_dependence))
            ! Closing current files
            call openclose_dclimit_currents_files(2)

            if(.not.lhfresponses) then
              ! Opening B effective files
              call openclose_dclimit_beff_files(1)
              ! Writing effective fields
              call write_dclimit_beff(dclimit_field_dependence(dcfield_dependence))
              ! Closing B effective files
              call openclose_dclimit_beff_files(2)
            end if

            ! WRITING TORQUES
            ! Opening torque files
            call openclose_dclimit_torques_files(1)
            ! Writing torques
            call write_dclimit_torques(dclimit_field_dependence(dcfield_dependence))
            ! Closing torque files
            call openclose_dclimit_torques_files(2)

          end do MPI_points_dclimit

          call date_and_time(date, time, zone, values)
          write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

          ! Emergency stop
          open(unit=911, file="stop", status='old', iostat=iw)
          if(iw.eq.0) then
            close(911)
            write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
            call system ('rm stop')
            write(*,"(a)") "[main] ('stop' file deleted!)"
            call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
          end if
        else
          call MPI_Send(dclimit_field_dependence,3,MPI_DOUBLE_PRECISION,0,44000,MPIComm_Col,ierr)
          if(.not.lhfresponses) then
!             call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,44010,MPIComm_Col,ierr)
!             call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,44020,MPIComm_Col,ierr)
            call MPI_Send(Beff_cart,dimsigmaNpl,MPI_DOUBLE_COMPLEX,0,44170,MPIComm_Col,ierr)
          end if
          call MPI_Send(disturbances,7*Npl,MPI_DOUBLE_COMPLEX,0,44030,MPIComm_Col,ierr)
          call MPI_Send(currents,7*n0sc*Npl,MPI_DOUBLE_COMPLEX,0,44100,MPIComm_Col,ierr)
          call MPI_Send(total_currents,7*n0sc,MPI_DOUBLE_COMPLEX,0,44110,MPIComm_Col,ierr)
          call MPI_Send(torques,6*Npl,MPI_DOUBLE_COMPLEX,0,44120,MPIComm_Col,ierr)
        end if
      end if

      call deallocate_prefactors()
      call deallocate_susceptibilities()
      call deallocate_disturbances()
      call deallocate_currents()
      call deallocate_beff()
      call deallocate_torques()
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!     case (12)
!       if(myrank.eq.0) then
!         allocate(sdl(Npl), STAT = AllocateStatus)
!         if (AllocateStatus.ne.0) then
!           call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
!         end if

!         write(*,"('CALCULATING PROBABILITY OF SPIN FLIP AS A FUNCTION OF POSITION')")
!         do i=1,Npl
!           write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/sdl_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'_magaxis=',A,'_socscale=',f5.2,'.dat')") SOC,Npl,i,parts,parts3,ncp,eta,Utype,hwa,hwt,hwp,magaxis,socscale
!           open (unit=555+i, file=varm,status='unknown')
!         end do

!         e = Ef
!         inn = 1
!         do pos=0,100
!           call spindifflength(e,pos,inn,sdl)
!           do i=1,Npl
!             sdl2 = abs(sdl(i))**2
!             ! Writing probability of spin-flip as a function os position
!             write(unit=555+i,fmt="(i0,2x,es16.9)") pos,sdl2
!           end do
!         end do
!         deallocate(sdl)
!       end if
!-----------------------------------------------------------------------


    end select main_program
!-------------- Deallocating variables that depend on Npl --------------
    deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1,sc_solu)
    deallocate(diag,qtf,w)
    deallocate(mx,my,mz,hdel,mp,hdelp,mm,hdelm,fvec,jac,wa)
    deallocate(mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi)
    if(lGSL) deallocate(lxm,lym,lzm,lxpm,lypm,lzpm)
    deallocate(mmlayer,layertype,U,mmlayermag,lambda,npart0)
    deallocate(identt,Umatorb)
    select case (plnn)
    case(1)
      deallocate(t00,t01)
    case(2)
      deallocate(t00,t01,t02)
    end select

    ! Emergency stop after the calculation for certain Npl is finished
    if(Npl_f.ne.Npl_i) then
      open(unit=911, file="stopNpl", status='old', iostat=iw)
      if(iw.eq.0) then
        close(911)
        write(*,"(a,i0,a)") "[main] Emergency 'stopNpl' file found! Stopping after Npl = ",Npl," ..."
        call system ('rm stopNpl')
        write(*,"(a)") "[main] ('stopNpl' file deleted!)"
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if
!---------------------------- Closing loops ----------------------------
  end do hw_phi ; end do hw_theta ; end do hw_intensity ; end do number_of_planes
!----------------------- Deallocating variables ------------------------
  deallocate(r0,c0,r1,c1,r2,c2)
  deallocate(kbz,wkbz,kbz2d)
!----------------------- Finalizing the program ------------------------
  if(myrank.eq.0) then
    call date_and_time(date, time, zone, values)
    write(*,"('[main] Finished on: ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") values(3),values(2),values(1),values(5),values(6),values(7)
    elapsed_time = MPI_Wtime() - start_program
    write(*,"('[main] Total elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
  end if
  call MPI_Finalize(ierr)
  if (ierr.ne.0) then
    write(*,"('[main] ierr = ',i0,'. Something went wrong in the parallelization!')") ierr
  end if
!=======================================================================
  stop
end program DHE