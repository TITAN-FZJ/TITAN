program SHE
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_magnet
  use mod_tight_binding
  use mod_prefactors
  use mod_generate_epoints
  use mod_generate_kpoints
  use mod_lattice
  use mod_progress
  use mod_mpi_pars
  use mod_dnsqe
  use MPI
  implicit none
  character(len=400)            :: varm
  character(len=8)              :: date
  character(len=10)             :: time
  character(len=5)              :: zone
  integer                       :: i,j,iw,pos,err,sigma,mu,nu,gamma,xi,inn,neighbor,iflag
  integer                       :: AllocateStatus,values(8),lwa,ifail=0
  real(double)                  :: e,sdl2,Icabs
  complex(double),allocatable   :: schi(:,:,:),schihf(:,:,:)
  complex(double),allocatable   :: sdx(:),sdy(:),sdz(:),chd(:),ldx(:),ldy(:),ldz(:)
  complex(double),allocatable   :: Ich(:,:),Isx(:,:),Isy(:,:),Isz(:,:),Ilx(:,:),Ily(:,:),Ilz(:,:)
  complex(double),allocatable   :: rsdx(:),rsdy(:),rsdz(:),rchd(:),rldx(:),rldy(:),rldz(:)
  complex(double),allocatable   :: rIch(:,:),rIsx(:,:),rIsy(:,:),rIsz(:,:),rIlx(:,:),rIly(:,:),rIlz(:,:)
  complex(double),allocatable   :: sdl(:)

  ! Exchange interaction
  real(double),dimension(:,:),allocatable     :: trJij
  real(double),dimension(:,:,:,:),allocatable :: Jij,Jijs,Jija

  ! LDOS
  real(double),dimension(:,:),allocatable :: ldosu,ldosd

  ! Self consistency variables
  logical                       :: selfcon
  real(double)                  :: mdif,sdif
  real(double),allocatable      :: npart(:),jac(:,:),wa(:),mag_0(:),eps1_solu(:)
  complex(double),allocatable   :: splus_0(:)

  ! Diagonalize susceptibilities
  integer                       :: lwork
  real(double),allocatable      :: rwork(:)
  complex(double),allocatable   :: eval(:),work(:)
  complex(double), dimension(:,:),allocatable       :: chimag,evecl,evecr
#ifdef _JUQUEEN
  integer :: ilo,ihi
  real(double) :: abnrm
  real(double),allocatable      :: dscale(:),rconde(:),rcondv(:)
#endif

  complex(double), dimension(:,:),   allocatable :: identt,temp,Umatorb
  complex(double), dimension(:,:),   allocatable :: chiorb_hf,chiorb,tchiorbiikl
  complex(double), dimension(:,:,:), allocatable :: templd
  complex(double), dimension(:,:,:), allocatable :: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl
!---------------------------- begin MPI vars ---------------------------
#ifndef _UFF
!$  integer :: provided
#endif
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^  end MPI vars ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  external nparticles,nparticlesjac,nparticlesjacnag
!------------------------ begin MPI initialization ---------------------
#ifdef _OPENMP
#ifdef _UFF
  call MPI_Init(ierr)
#else
  call MPI_Init_thread(MPI_THREAD_MULTIPLE,provided,ierr)
#endif
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
  if(myrank.eq.0) write(*,"('*** Running with openMP ***')")
#else
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
  if(myrank.eq.0) write(*,"('*** Running without openMP ***')")
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
!------------- Generating Cunningham's k points in 2D BZ ---------------
  select case (lattice)
  case("bcc110")
    call bcc110()
    call generate_kpoints_bcc110()
  case("fcc100")
    call fcc100()
    call generate_kpoints_fcc100()
  case default
    if(myrank.eq.0) write(*,"('[main] Lattice not defined: ',a,'!')") lattice
    stop
  end select
!---- Generating integration points of the complex energy integral -----
  call generate_imag_epoints()
!------------------------- Number of planes loop -----------------------
  number_of_planes: do Npl = Nplini,Nplfinal
!--------------- Allocating variables that depend on Npl ---------------
    lwa=Npl*(3*Npl+13)/2
    allocate( sigmai2i(4,Npl),sigmaimunu2i(4,Npl,9,9),sigmaijmunu2i(4,Npl,Npl,9,9),eps1(Npl),eps1_solu(Npl), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: sigmaimunu2i,sigmaijmunu2i,eps1,eps1_solu')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    allocate( mag(Npl),mag_0(Npl),hdel(Npl),splus(Npl),splus_0(Npl),usplus(Npl),sminus(Npl),usminus(Npl),npart(Npl),jac(Npl,Npl),wa(lwa), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: mag,mag_0,hdel,splus,splus_0,usplus,sminus,usminus,npart,jac,wa')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    allocate( mmlayer(Npl+2),layertype(Npl+2),U(Npl+2),mmlayermag(Npl+2),lambda(Npl+2),npart0(Npl+2), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[main] Not enough memory for: mmlayer,layertype,U,mmlayermag,lambda,npart0')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
!--------------------- Dimensions and identities- --------------------
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
!----------- Creating bi-dimensional matrix of MPI processes  -----------
    ! Bidimensional array should look like:
    !      0 1 2 3 ... pnt-1
    !    0
    !    1
    !    2
    !   ...
    ! MPIpts-1
    if((itype.ge.5).and.(itype.le.8)) then ! Create matrix only when energy integration is involved
      MPIpts = ceiling(dble(numprocs)/dble(pnt)) ! Number of rows to be used
      MPIdims = [MPIpts,pnt]
      if(numprocs.le.pnt) then  ! If number of processes is less than necessary for 1 energy integral
        MPIdims  = [MPIpts,numprocs]  ! Create only one array of processes, i.e., MPIpts = 1
        MPIsteps = npt1
      else
        if(mod(numprocs,pnt).ne.0) then
          if(myrank.eq.0) then
            write(*,"('[main] ERROR: number of processes not commensurable with total energy integral points!')")
            write(*,"('[main] Number of MPI processes: ',i0)") numprocs
            write(*,"('[main] Number of points required: ',i0)") npt1
            write(*,"('[main] Number of points in the energy integral: ',i0)") pnt
          end if
          stop
        end if
        MPIsteps = 1
        if(npt1*pnt.lt.numprocs) then ! If the number of processors is larger than complete calculation
          if(myrank.eq.0) then
            write(*,"('[main] ******************************** WARNING: ********************************')")
            write(*,"('[main]              Number of processes exceeds the total needed! ')")
            write(*,"('[main]     Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIpts-npt1
            write(*,"('[main] **************************************************************************')")
          end if
          npt1 = MPIpts
        else if(npt1*pnt.gt.numprocs) then ! If the number of processors is smaller than complete calculation, check commensurability
          MPIsteps = ceiling(dble(npt1)/dble(MPIpts))
          if(mod(npt1,MPIpts).ne.0) then
            if(myrank.eq.0) then
              write(*,"('[main] ******************************** WARNING: ********************************')")
              write(*,"('[main] Number of points to be calculated is not commensurable with processes used!')")
              write(*,"('[main]     Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIsteps*MPIpts-npt1
              write(*,"('[main] **************************************************************************')")
            end if
            npt1 = MPIsteps*MPIpts
          end if
        end if
      end if
      if(npt1.ne.1) then
        npts = npt1-1
      else
        npts = npt1
      end if
      ! Calculating variations of energy
      deltae = (emax - emin)/npts      ! variation of energy between each point
      MPIdelta = deltae*MPIpts         ! variation of energy between MPI steps

      ! Creating bidimensiontal Grid of tasks
      lperiodic = [.false.,.false.]
      lreorder  = .true.
      call MPI_Cart_create(MPI_COMM_WORLD,2,MPIdims,lperiodic,lreorder,MPIComm_Grid,ierr)

      ! Creating subarrays of rows and columns
      lrow = [.true.,.false.]
      call MPI_Cart_sub(MPIComm_Grid,lrow,MPIComm_Col,ierr) ! communicator inside a column (between rows)
      lcol = [.false.,.true.]
      call MPI_Cart_sub(MPIComm_Grid,lcol,MPIComm_Row,ierr) ! communicator inside a row (between columns)

      ! Obtaining process rank inside its row and column
      call MPI_Comm_rank(MPIComm_Row,myrank_row,ierr) ! Obtaining rank number inside its row
      call MPI_Comm_rank(MPIComm_Col,myrank_col,ierr) ! Obtaining rank number inside its column
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
        Umatorb(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,i,gamma,xi)) = cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,i,gamma,xi)) = cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,i,gamma,xi)) = cmplx(U(i+1),0.d0)
        Umatorb(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,i,gamma,xi)) = cmplx(U(i+1),0.d0)
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
    !   Out-of-plane:
        phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
        theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
      case default
        if(myrank.eq.0) write(*,"('[main] Choose a correct magnetization axis!')")
        stop
    !     phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
    !     theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
      end select magnetization_axis_fcc100
!----------- Direction of applied in-plane electric field --------------
      read(unit=dirEfield,fmt=*) j
      direction_E_field_fcc100: select case (j)
      case (1:8)   !    In plane neighbors:
        dirEfieldvec = r0(j,:)
      case default !    Other direction:
        dirEfieldvec = [1.d0,0.d0,0.d0] ! In-plane
      end select direction_E_field_fcc100
    end select
!-------------------- Writing parameters on screen ---------------------
    if(myrank.eq.0) then
      write(*,"('Running on ',i0,' MPI process(es)')") numprocs
      write(*,"('|------------- PARAMETERS: -------------|')")
      write(*,"(10x,'Npl = ',i0)") Npl
      write(*,"(1x,'DFT parameters: ',$)")
      dft_type: select case (dfttype)
      case ("T")
        write(*,"('Tight-binding basis')")
      case ("O")
        write(*,"('Orthogonal basis')")
      end select dft_type
      if(SOC) then
        write(*,"(1x,'Spin Orbit Coupling: ACTIVATED')")
      else
        write(*,"(1x,'Spin Orbit Coupling: DEACTIVATED')")
      end if
      write(*,"(12x,'Utype = ',i0)") Utype

      select case (lattice)
      case("bcc110")
        write(*,"(1x,'Spin quantization direction: ',$)")
        write_magnetization_axis_bcc110: select case (magaxis)
        case ("L")
          write(*,"('Long axis')")
        case ("S")
          write(*,"('Short axis')")
        case ("P")
          write(*,"('Perpendicular out-of-plane')")
        case default
          write(*,"(11x,'Other')")
        end select write_magnetization_axis_bcc110
        write(*,"(10x,'phi =',e9.2,'pi')") phi/pi
        write(*,"(8x,'theta =',e9.2,'pi')") theta/pi
        write(*,"(1x,'Electric field direction: ',$)")
        write_direction_E_field_bcc110: select case (dirEfield)
        case ("L")
          write(*,"('Long axis')")
        case ("S")
          write(*,"('Short axis')")
        case ("O")
          write(*,"('Other',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)
        end select write_direction_E_field_bcc110
      case("fcc100")
        write(*,"(1x,'Spin quantization direction: ',$)")
        read(unit=magaxis,fmt=*) j
        write_magnetization_axis_fcc100: select case (j)
        case (1:8)   !    In plane neighbors:
          write(*,"('Neighbor ',i0)") j
        case (9)
          write(*,"('Out-of-plane')")
        end select write_magnetization_axis_fcc100
        write(*,"(10x,'phi =',e9.2,'pi')") phi/pi
        write(*,"(8x,'theta =',e9.2,'pi')") theta/pi
        write(*,"(1x,'Electric field direction: ',$)")
        read(unit=dirEfield,fmt=*) j
        write_direction_E_field_fcc100: select case (j)
        case (1:8)   !    In plane neighbors:
          write(*,"('neighbor ',i0)") j
        case default
          write(*,"('Other',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)
        end select write_direction_E_field_fcc100
      end select
      if(renorm) then
        write(*,"(1x,'Current renormalization: ACTIVATED')")
        write(*,"(5x,'renormnb = ',i0)") renormnb
      else
        write(*,"(1x,'Current renormalization: DEACTIVATED')")
      end if
      write(*,"(10x,'ncp = ',i0)") ncp
      write(*,"(8x,'parts = ',i0,'x',i0)") parts,n1gl
      write(*,"(7x,'parts3 = ',i0,'x',i0)") parts3,n3gl
      write(*,"(10x,'eta =',e9.2)") eta
      write(*,"(10x,'hwx =',e9.2)") hwx
      write(*,"(10x,'hwy =',e9.2)") hwy
      write(*,"(10x,'hwz =',e9.2)") hwz
      if(runoptions.ne."") write(*,"(6x,'Activated options:',/,4x,a)") trim(runoptions)
    end if
!--------------- Writing data to be calculated on screen ---------------
    if(myrank.eq.0) then
      write(*,"('|------------ TO CALCULATE: ------------|')")
      write_itype: select case (itype)
      case (0)
        write(*,"(1x,'Test before SC')")
        write(*,"(8x,'n0sc1 = ',i0)") n0sc1
        write(*,"(8x,'n0sc2 = ',i0)") n0sc2
        write(*,"(9x,'emin =',e9.2)") emin
        write(*,"(9x,'emax =',e9.2)") emax
        write(*,"(1x,'Number of points to calculate: ',i0)") npt1
      case (1)
        write(*,"(1x,'Self-consistency only')")
      case (2)
        write(*,"(1x,'Test after SC')")
        write(*,"(8x,'n0sc1 = ',i0)") n0sc1
        write(*,"(8x,'n0sc2 = ',i0)") n0sc2
        write(*,"(9x,'emin =',e9.2)") emin
        write(*,"(9x,'emax =',e9.2)") emax
        write(*,"(1x,'Number of points to calculate: ',i0)") npt1
      case (3)
        write(*,"(1x,'LDOS and exchange interactions as a function of energy')")
        write(*,"(9x,'emin =',e9.2)") emin
        write(*,"(9x,'emax =',e9.2)") emax
        write(*,"(1x,'Number of points to calculate: ',i0)") npt1
      case (4)
        write(*,"(1x,'Band structure')")
        write(*,"(3x,'kdirection = ',a)") kdirection
        write(*,"(9x,'Number of points to calculate: ',i0)") npt1
      case (5)
        write(*,"(1x,'Local susceptibility as a function of energy')")
        write(*,"(9x,'emin =',e9.2)") emin
        write(*,"(9x,'emax =',e9.2)") emax
        write(*,"(5x,i0,' points')") npt1
        write(*,"(1x,'divided into ',i0,' steps')") MPIsteps
        write(*,"(1x,'of size ',e10.3)") MPIdelta
        write(*,"(1x,'each calculating ',i0,' points')") MPIpts
      case (6)
        write(*,"(1x,'Disturbances and local susceptibility as a function of energy')")
        write(*,"(9x,'emin =',e9.2)") emin
        write(*,"(9x,'emax =',e9.2)") emax
        write(*,"(5x,i0,' points')") npt1
        write(*,"(1x,'divided into ',i0,' steps')") MPIsteps
        write(*,"(1x,'of size ',e10.3)") MPIdelta
        write(*,"(1x,'each calculating ',i0,' points')") MPIpts
      case (7)
        write(*,"(1x,'Parallel currents and local susceptibility as a function of energy')")
        write(*,"(8x,'n0sc1 = ',i0)") n0sc1
        write(*,"(8x,'n0sc2 = ',i0)") n0sc2
        write(*,"(9x,'emin =',e9.2)") emin
        write(*,"(9x,'emax =',e9.2)") emax
        write(*,"(5x,i0,' points')") npt1
        write(*,"(1x,'divided into ',i0,' steps')") MPIsteps
        write(*,"(1x,'of size ',e10.3)") MPIdelta
        write(*,"(1x,'each calculating ',i0,' points')") MPIpts
      case (8)
        write(*,"(1x,'Parallel currents, disturbances and local susc. as a function of energy')")
        write(*,"(8x,'n0sc1 = ',i0)") n0sc1
        write(*,"(8x,'n0sc2 = ',i0)") n0sc2
        write(*,"(9x,'emin =',e9.2)") emin
        write(*,"(9x,'emax =',e9.2)") emax
        write(*,"(5x,i0,' points')") npt1
        write(*,"(1x,'divided into ',i0,' steps')") MPIsteps
        write(*,"(1x,'of size ',e10.3)") MPIdelta
        write(*,"(1x,'each calculating ',i0,' points')") MPIpts
      case (9)
        write(*,"(1x,'Fermi surface')")
      case (10)
        write(*,"(1x,'Exhange interactions and anisotropies (full tensor)')")
        if(nmaglayers.eq.0) then
          write(*,"(1x,'[main] No magnetic layers!')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        if(nmaglayers.eq.1) write(*,"(1x,'Only 1 magnetic layer: calculating only anisotropies')")
        write(*,"(8x,'from Npl = ',i0,' to ',i0)") Nplini,Nplfinal
      case (11)
        write(*,"(1x,'Probability of spin-flip as a function of position')")
      end select write_itype
      write(*,"('|---------------------------------------|')")
    end if
!---- L matrix in spin coordinates for given quantization direction ----
    call lp_matrix()
!------ Calculate L.S matrix for the given quantization direction ------
    call ls_matrix()
!------ Calculate L.B matrix for the given quantization direction ------
    call lb_matrix()
!------------------------ Calculate S.B matrix  ------------------------
    call sb_matrix()
!---------------------------- Debugging part ---------------------------
    if(index(runoptions,"idebug").gt.0) then
      if(myrank.eq.0) then
        write(*,"('Debugging...')")
        call debugging()
      end if
      call MPI_Finalize(ierr)
      if (ierr.ne.0) then
        write(*,"('[main] ierr = ',i0,'. Something went wrong in the parallelization!')") ierr
      end if
      stop
    end if
!------------------- Only create files with headers --------------------
    if(index(runoptions,"createfiles").ne.0) then
      if(myrank.eq.0) then
        create_files: select case (itype)
        case (5)
          call openclosechifiles(0)
          write(*,"('[main] Susceptibilities files created/overwritten!')")
        case (6)
          call openclosechifiles(0)
          call openclosesdfiles(0)
          write(*,"('[main] Susceptibilities and disturbances files created/overwritten!')")
        case (7)
          call openclosechifiles(0)
          call openclosescfiles(0)
          write(*,"('[main] Susceptibilities and currents files created/overwritten!')")
        case (8)
          call openclosechifiles(0)
          call openclosesdfiles(0)
          call openclosescfiles(0)
          write(*,"('[main] Susceptibilities, disturbances and current files created/overwritten!')")
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
!-------------------------- Begin first test part ----------------------
    if((myrank.eq.0).and.(itype.eq.0)) then
      write(*,"('FIRST TEST PART')")

      allocate(ldosu(Npl,9),ldosd(Npl,9))
      allocate(trJij(nmaglayers,nmaglayers),Jij(nmaglayers,nmaglayers,3,3),Jijs(nmaglayers,nmaglayers,3,3),Jija(nmaglayers,nmaglayers,3,3))

      ! Parameters: center of band, magnetization, correlation functions
      eps1  = 0.d0
      mag   = 0.d0
      splus = zero
      ! Variables used in the hamiltonian
      do i=1,Npl
        hdel(i)    = 0.5d0*U(i+1)*mag(i)
        usplus(i)  = U(i+1)*splus(i)
      end do
      usplus = zero
      usminus = conjg(usplus)

      emin = -2.d0  ! given in eV
      emax = 2.d0   ! given in eV
      npts = 400
      npt1 = npts+1
      test1_energy_loop2: do count=1,npt1
        e = emin + (count-1)*deltae
        e = (e/ry2ev)!+Ef ! Transform to Ry
        write(*,"(i0,' of ',i0,' points',', e = ',e10.3)") count,npt1,e

        ! Turning off SOC
        lambda = 0.d0

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
  ! Reading previous shifts and mag from files
  ! (Try to read eps1 and mag if available - includes hdel, usplus and usminus calculations)
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
      eps1_solu = eps1
    case default !  If file doesn't exist
      if(myrank.eq.0) write(*,"('Self-consistency file does not exist.')")
      selfcon = .true.
      ! Parameters: center of band, magnetization, correlation functions
      eps1_solu=0.d0
      mag = 0.5d0
      do i=1,Npl
        if(layertype(i+1).eq.2) mag(i) = 2.d0
      end do
      splus = zero
      ! Variables used in the hamiltonian
      do i=1,Npl
        hdel(i)    = 0.5d0*U(i+1)*mag(i)
        usplus(i)  = U(i+1)*splus(i)
      end do
      usminus = conjg(usplus)
    end select sc_file_status

    ! Self-consistency
    if(selfcon) then
      ifail = 0
      mdif  = 100.d0
      sdif  = 100.d0
      mag_0 = mag
      splus_0 = splus
      j     = 1
      if(myrank.eq.0) write(*,"('Starting self-consistency:')")
      self_consistency: do while ((mdif.gt.tol).or.(sdif.gt.tol))
        if(myrank.eq.0) write(*,"('|------------------- Iteration ',i0,' -------------------|')") j

        if(index(runoptions,"slatec").ne.0) then
          call dnsqe(nparticles,nparticlesjac,1,Npl,eps1_solu,npart,tol,0,ifail,wa,lwa)
          ifail = ifail-1
        else
          call c05pbf(nparticlesjacnag,Npl,eps1_solu,npart,jac,Npl,tol,wa,lwa,ifail)
  !         call c05nbf(nparticles,Npl,eps1_solu,npart,tol,wa,lwa,ifail)
        end if

        if (ifail.ne.0) then
          write(*,"('[main] Problem in self-consistency in rank ',i0,'(',a,')! ifail = ',i0)") myrank,trim(host),ifail
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        ! Using calculated values in the hamiltonian
  !       do i=1,Npl
  !         hdel(i)    = 0.25d0*U(i+1)*(mag(i)+mag_0(i))
  !         usplus(i)  = 0.5d0*U(i+1)*(splus(i)+splus_0(i))
  !       end do
        do i=1,Npl
          hdel(i)    = 0.5d0*U(i+1)*mag(i)
          usplus(i)  = U(i+1)*splus(i)
        end do
        usminus = conjg(usplus)
        ! Calculating differences between steps
        mdif  = sum(abs(mag-mag_0))
        sdif  = sum(abs(splus-splus_0))
        ! Redefining original correlations
        mag_0 = mag
        splus_0 = splus
        j   = j+1
        if(myrank.eq.0) then
          write(*,"('Differences to last step:')")
          write(*,"('mdif = ',e16.9)") mdif
          write(*,"('sdif = ',e16.9)") sdif
        end if

        ! Writing new eps1 and mag to file during self-consistency
        if(myrank.eq.0) then
          iflag = 1
          call readwritesc(iflag,err)
        end if

      end do self_consistency

      ! Writing new eps1 and mag to file after self-consistency
      if(myrank.eq.0) then
        iflag = 1
        call readwritesc(iflag,err)
      end if

    end if

    if(index(runoptions,"GSL").gt.0) then
      allocate( lxm(Npl),lym(Npl),lzm(Npl),lxpm(Npl),lypm(Npl),lzpm(Npl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        if(myrank.eq.0) write(*,"('[main] Not enough memory for: lxm,lym,lzm,lxpm,lypm,lzpm')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if

      if(myrank.eq.0) write(*,"('[main] Calculating Orbital Angular Momentum ground state... ')")
      call L_gs()
    end if

    if(myrank.eq.0) then
      write(*,"('|---------------------- Self-consistent ground state: ----------------------|')")
      write(*,"(11x,' *************** Center of d bands: ***************')")
      do i=1,Npl
        write(*,"(26x,'eps1(',I0,')=',e16.9)") i,eps1(i)
      end do
      write(*,"(11x,' *************** Spin components: ***************')")
      do i=1,Npl
        write(*,"('Sx (',I0,')=',e16.9,4x,'Sy (',I0,')=',e16.9,4x,'Sz (',I0,')=',e16.9)") i,real(splus(i)),i,aimag(splus(i)),i,mag(i)
      end do
      if(index(runoptions,"GSL").gt.0) then
        write(*,"(11x,' *** Orbital components in spin coordinates:  ***')")
        do i=1,Npl
          write(*,"('Lxp(',I0,')=',e16.9,4x,'Lyp(',I0,')=',e16.9,4x,'Lzp(',I0,')=',e16.9)") i,lxpm(i),i,lypm(i),i,lzpm(i)
        end do
        write(*,"(11x,' *** Orbital components in cubic coordinates: ***')")
        do i=1,Npl
          write(*,"('Lx (',I0,')=',e16.9,4x,'Ly (',I0,')=',e16.9,4x,'Lz (',I0,')=',e16.9)") i,lxm(i),i,lym(i),i,lzm(i)
        end do
        write(*,"(11x,' ******************** Total: ********************')")
        do i=1,Npl
          write(*,"('S (',I0,') =',e16.9,4x,'Lp (',I0,')=',e16.9,4x,'L (',I0,') =',e16.9)") i,sqrt((real(splus(i))**2)+(aimag(splus(i))**2)+(mag(i)**2)),i,sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2)),i,sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
        end do
      else
        write(*,"(11x,' ******************** Total: ********************')")
        do i=1,Npl
          write(*,"(27x,'S (',I0,') =',e16.9)") i,sqrt((real(splus(i))**2)+(aimag(splus(i))**2)+(mag(i)**2))
        end do
      end if
      write(*,"('|---------------------------------------------------------------------------|')")
      elapsed_time = MPI_Wtime() - start_program
      write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
    end if

    if(itype.eq.1) then
      ! Finalizing program
      if(myrank.eq.0) then
        call date_and_time(date, time, zone, values)
        write(*,"('[main] Finished on: ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") values(3),values(2),values(1),values(5),values(6),values(7)
        elapsed_time = MPI_Wtime() - start_program
        write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
      end if
      call MPI_Finalize(ierr)
      if (ierr.ne.0) then
        write(*,"('[main] ierr = ',i0,'. Something went wrong in the parallelization!')") ierr
      end if
      stop
    end if
!============================= MAIN PROGRAM ============================
    main_program: select case (itype)
    case (2)
!----------------------------- Begin test part -----------------------
      if(myrank.eq.0) then
        write(*,"('TESTING')")

        test2_energy_loop: do count=1,npt1
          e = emin + (count-1)*deltae
          write(*,"(i0,' of ',i0,' points',', e = ',e10.3)") count,npt1,e



        end do test2_energy_loop
      end if
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
          write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/LDOS/ldosu_layer',I0,'_magaxis=',A,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,i,magaxis,ncp,eta,Utype,hwx,hwy,hwz
          open (unit=iw, file=varm,status='unknown')
          iw = iw+1
          write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/LDOS/ldosd_layer',I0,'_magaxis=',A,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,i,magaxis,ncp,eta,Utype,hwx,hwy,hwz
          open (unit=iw, file=varm,status='unknown')
        end do
        ! Exchange interactions
        do j=1,nmaglayers ; do i=1,nmaglayers
          iw = 99+(j-1)*nmaglayers*2+(i-1)*2
          if(i.eq.j) then
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/Jij/Jii_',i0,'_magaxis=',A,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,mmlayermag(i)-1,magaxis,ncp,eta,Utype,hwx,hwy,hwz
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   energy      ,  Jii_xx           ,   Jii_yy  ')")
            iw = iw + 1
            ! TODO : Check how to write the anisotropy term here
          else
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/Jij/J_',i0,'_',i0,'_magaxis=',A,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,mmlayermag(i)-1,mmlayermag(j)-1,magaxis,ncp,eta,Utype,hwx,hwy,hwz
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   energy      ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/Jij/Dz_',i0,'_',i0,'_magaxis=',A,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,mmlayermag(i)-1,mmlayermag(j)-1,magaxis,ncp,eta,Utype,hwx,hwy,hwz
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   energy      , Dz = (Jxy - Jyx)/2       ')")
          end if
        end do ; end do

        ldos_energy_loop: do count=1,npt1
          e = emin + (count-1)*deltae
          write(*,"(i0,' of ',i0,' points',', e = ',e10.3)") count,npt1,e

          call ldos(e,ldosu,ldosd,Jij)

          do i=1,nmaglayers ; do j=1,nmaglayers
            trJij(i,j)    = 0.5d0*(Jij(i,j,1,1)+Jij(i,j,2,2))
            Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
            Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
            do mu=1,3
              Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
            end do
          end do ; end do

          ! Writing into files
          ! LDOS
          ldos_writing_plane_loop: do i=1,Npl
              iw = 17+(i-1)*2
              write(unit=iw,fmt="(5(e16.9,2x))") e*ry2ev,sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
              iw = iw+1
              write(unit=iw,fmt="(5(e16.9,2x))") e*ry2ev,sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
          end do ldos_writing_plane_loop

          ! Exchange interactions
          jij_writing_loop: do j=1,nmaglayers ; do i=1,nmaglayers
            iw = 99+(j-1)*nmaglayers*2+(i-1)*2
            if(i.eq.j) then
              iw = iw + 1
              write(unit=iw,fmt="(3(e16.9,2x))") e*ry2ev,Jij(i,j,1,1),Jij(i,j,2,2)
              iw = iw + 1
            else
              iw = iw + 1
              write(unit=iw,fmt="(4(e16.9,2x))") e*ry2ev,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
              iw = iw + 1
              write(unit=iw,fmt="(2(e16.9,2x))") e*ry2ev,Jija(i,j,1,2)
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
      if(myrank_row.eq.0) then
        allocate( schi(Npl,Npl,4),schihf(Npl,Npl,4), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          if(myrank.eq.0) write(*,"('[main] Not enough memory for: schi,schihf')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        allocate( temp(dim,dim),chiorb(dim,dim), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          write(*,"('[main] Not enough memory for: temp,chiorb')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if

        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            lwork = 33*nmaglayers
            allocate( chimag(nmaglayers,nmaglayers),rwork(2*nmaglayers),eval(nmaglayers),evecl(1,nmaglayers),evecr(nmaglayers,nmaglayers),work(lwork) )
#ifdef _JUQUEEN
            allocate( dscale(nmaglayers),rconde(nmaglayers),rcondv(nmaglayers) )
#endif
          end if
          write(*,"('CALCULATING LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
          write(*,"('Qx = ',e10.3,', Qz = ',e10.3)") q(1),q(2)
          ! Creating files and writing headers
          if(index(runoptions,"addresults").eq.0) then
            call openclosechifiles(0)
          end if
        end if
      end if
      allocate( chiorb_hf(dim,dim), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if

      q     = [0.d0, 0.d0]
      chi_energy_loop: do count=1,MPIsteps
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',e10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e,myrank_col
        end if

        ! Start parallelized processes to calculate chiorb_hf and chiorbi0_hf for energy e
        call eintshechi(e,chiorb_hf)

        if(myrank_row.eq.0) then
          ! (1 + chi_hf*Umat)^-1
          call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zero,temp,dim)
          temp = identt + temp
          call invers(temp,dim)
          call zgemm('n','n',dim,dim,dim,zum,temp,dim,chiorb_hf,dim,zero,chiorb,dim)

          schi = zero
          schihf = zero
          calculate_susceptibility: do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
            ! Calculating RPA and HF susceptibilities
            do nu=5,9 ; do mu=5,9
              schi(i,j,sigma) = schi(i,j,sigma) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
              schihf(i,j,sigma) = schihf(i,j,sigma) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
            end do ; end do
          end do ; end do ; end do calculate_susceptibility

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            ! Opening files for writing
            call openclosechifiles(1)
            MPI_points_chi: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,1000,MPIComm_Col,stat,ierr)
                call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1100,MPIComm_Col,stat,ierr)
                call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1200,MPIComm_Col,stat,ierr)
              end if

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              write(*,"(' #################  Susceptibilities:  #################')")
              write_susceptibility: do sigma=1,4 ; do j=1,Npl ;  do i=1,Npl
                iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i

                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schi(i,j,sigma)),aimag(schi(i,j,sigma)),abs(schi(i,j,sigma))
                if((sigma.eq.1).and.(i.eq.j)) write(*,"('e = ',e11.4,', plane: ',i0,' chi = (',e16.9,') + i(',e16.9,')')") e,i,real(schi(i,j,sigma)),aimag(schi(i,j,sigma))

                iw = iw+1000
                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schihf(i,j,sigma)),aimag(schihf(i,j,sigma)),abs(schihf(i,j,sigma))
              end do ; end do ; end do write_susceptibility

              ! DIAGONALIZING SUSCEPTIBILITY
              if(nmaglayers.gt.1) then
                do i=1,nmaglayers ; do j=1,nmaglayers
                  chimag(i,j) = schi(mmlayermag(i)-1,mmlayermag(j)-1,1)
                end do ; end do

#ifdef _JUQUEEN
                call zgeevx('N','N','V','N',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,ilo,ihi,dscale,abnrm,rconde,rcondv,work,lwork,rwork,ifail)
#else
                call zgeev('N','V',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,work,lwork,rwork,ifail)
#endif

                if(ifail.ne.0) then
                  write(*,*) '[main] Problem with diagonalization. ifail = ',ifail
                  stop
        !         else
        !           write(*,*) ' optimal lwork = ',work(1),' lwork = ',lwork
                end if
                write(varm,fmt="(a,i0,a)") '(',2*nmaglayers+1,'(e16.9,2x))'
                write(unit=1990,fmt=varm) e,(real(eval(i)),aimag(eval(i)),i=1,nmaglayers)
                do i=1,nmaglayers
                  write(unit=1990+i,fmt=varm) e,(real(evecr(j,i)),aimag(evecr(j,i)),j=1,nmaglayers)
                end do
              end if
            end do MPI_points_chi
            ! Closing files
            call openclosechifiles(2)

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
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

      deallocate(chiorb_hf)
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            deallocate(rwork,eval,evecl,evecr,work)
#ifdef _JUQUEEN
            deallocate(dscale,rconde,rcondv)
#endif
          end if
        end if
        deallocate(temp,chiorb)
        deallocate(schi,schihf)
      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (6)
      q     = [0.d0, 0.d0]
      if(myrank_row.eq.0) then
        allocate( schi(Npl,Npl,4),schihf(Npl,Npl,4),sdx(Npl),sdy(Npl),sdz(Npl),chd(Npl),ldx(Npl),ldy(Npl),ldz(Npl), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          if(myrank.eq.0) write(*,"('[main] Not enough memory for: schi,schihf,sdx,sdy,sdz,chd,ldx,ldy,ldz')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        allocate( templd(Npl,9,9),chiorb(dim,dim), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          write(*,"('[main] Not enough memory for: templd,chiorb')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            lwork = 33*nmaglayers
            allocate( chimag(nmaglayers,nmaglayers),rwork(2*nmaglayers),eval(nmaglayers),evecl(1,nmaglayers),evecr(nmaglayers,nmaglayers),work(lwork) )
#ifdef _JUQUEEN
            allocate( dscale(nmaglayers),rconde(nmaglayers),rcondv(nmaglayers) )
#endif
          end if
          write(*,"('CALCULATING DISTURBANCES AND LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
          write(*,"('Qx = ',e10.3,', Qz = ',e10.3)") q(1),q(2)
          ! Creating files and writing headers
          if(index(runoptions,"addresults").eq.0) then
            call openclosechifiles(0)
            call openclosesdfiles(0)
          end if
        end if
      end if
      allocate( chiorb_hf(dim,dim),prefactor(dim,dim),tchiorbiikl(dim,4), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(*,"('[main] Not enough memory for: chiorb_hf,prefactor,tchiorbiikl')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if

      sd_energy_loop: do count=1,MPIsteps
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',e10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e,myrank_col
        end if

        if(myrank.eq.0) write(*,"('[main] Calculating pre-factor to use in disturbances calculation ')")
        call eintshechi(e,chiorb_hf)

        ! Broadcast chiorb_hf to all processors of the same row
        call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)

        ! prefactor = (1 + chi_hf*Umat)^-1
        call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zero,prefactor,dim) !prefactor = chi_hf*Umat
        prefactor = identt + prefactor
        call invers(prefactor,dim)
        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Calculated prefactor after: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        ! Start parallelized processes to calculate disturbances for energy e
        call eintshesd(e,tchiorbiikl)

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
          calculate_susceptibility_sd: do nu=5,9 ; do mu=5,9 ; do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
            schi(i,j,sigma) = schi(i,j,sigma) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
            schihf(i,j,sigma) = schihf(i,j,sigma) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
          end do ; end do ; end do ; end do ; end do calculate_susceptibility_sd

          sdx = zero
          sdy = zero
          sdz = zero
          chd = zero
          ldx = zero
          ldy = zero
          ldz = zero
          calculate_sd: do i=1,Npl
            ! Spin and charge disturbances
            do mu=1,9
              sdx(i) = sdx(i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
              sdy(i) = sdy(i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
              sdz(i) = sdz(i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
              chd(i) = chd(i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
            end do

            ! Orbital angular momentum disturbance
            do nu=1,9; do mu=1,9
              templd(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
              ldx(i) = ldx(i) + lxp(mu,nu)*templd(i,mu,nu)
              ldy(i) = ldy(i) + lyp(mu,nu)*templd(i,mu,nu)
              ldz(i) = ldz(i) + lzp(mu,nu)*templd(i,mu,nu)
            end do; end do
          end do calculate_sd
          sdx = 0.5d0*sdx
          sdy = 0.5d0*sdy/zi
          sdz = 0.5d0*sdz

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            ! Opening files for writing
            call openclosechifiles(1)
            call openclosesdfiles(1)
            MPI_points_sd: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,2000,MPIComm_Col,stat,ierr)
                call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2100,MPIComm_Col,stat,ierr)
                call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2200,MPIComm_Col,stat,ierr)
                call MPI_Recv(sdx,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2300,MPIComm_Col,stat,ierr)
                call MPI_Recv(sdy,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2400,MPIComm_Col,stat,ierr)
                call MPI_Recv(sdz,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2500,MPIComm_Col,stat,ierr)
                call MPI_Recv(chd,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2600,MPIComm_Col,stat,ierr)
                call MPI_Recv(ldx,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2700,MPIComm_Col,stat,ierr)
                call MPI_Recv(ldy,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2800,MPIComm_Col,stat,ierr)
                call MPI_Recv(ldz,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),2900,MPIComm_Col,stat,ierr)
              end if

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              write_susceptibility_sd: do sigma=1,4 ; do j=1,Npl ;  do i=1,Npl
                iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i

                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schi(i,j,sigma)),aimag(schi(i,j,sigma)),abs(schi(i,j,sigma))

                iw = iw+1000
                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schihf(i,j,sigma)),aimag(schihf(i,j,sigma)),abs(schihf(i,j,sigma))
              end do ; end do ; end do write_susceptibility_sd

              ! DIAGONALIZING SUSCEPTIBILITY
              if(nmaglayers.gt.1) then
                do i=1,nmaglayers ; do j=1,nmaglayers
                  chimag(i,j) = schi(mmlayermag(i)-1,mmlayermag(j)-1,1)
                end do ; end do

#ifdef _JUQUEEN
                call zgeevx('N','N','V','N',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,ilo,ihi,dscale,abnrm,rconde,rcondv,work,lwork,rwork,ifail)
#else
                call zgeev('N','V',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,work,lwork,rwork,ifail)
#endif

                if(ifail.ne.0) then
                  write(*,*) '[main] Problem with diagonalization. ifail = ',ifail
                  stop
        !         else
        !           write(*,*) ' optimal lwork = ',work(1),' lwork = ',lwork
                end if
                write(varm,fmt="(a,i0,a)") '(',2*nmaglayers+1,'(e16.9,2x))'
                write(unit=1990,fmt=varm) e,(real(eval(i)),aimag(eval(i)),i=1,nmaglayers)
                do i=1,nmaglayers
                  write(unit=1990+i,fmt=varm) e,(real(evecr(j,i)),aimag(evecr(j,i)),j=1,nmaglayers)
                end do
              end if

              ! WRITING DISTURBANCES
              write_sd: do i=1,Npl
                write(*,"('|--------------- Plane: ',i0,' , Energy = ',e11.4,' ---------------|')") i,e
                write(*,"(' #################  Susceptibility:  #################')")
                write(*,"(' Chi+- = (',e16.9,') + i(',e16.9,')')") real(schi(i,i,1)),aimag(schi(i,i,1))

                ! Writing Spin, Charge and Orbital disturbances
                write(*,"(' ################# Charge disturbance: #################')")

                write(*,"('     Cd  = (',e16.9,') + i(',e16.9,')')") real(chd(i)),aimag(chd(i))
                write(*,"(' abs(Cd) = ',e16.9)") abs(chd(i))
                write(*,"('atan(Cd) = ',e16.9)") atan2(aimag(chd(i)),real(chd(i)))

                write(*,"(' ################# Spin disturbances:  #################')")

                write(*,"('     Sdx  = (',e16.9,') + i(',e16.9,')')") real(sdx(i)),aimag(sdx(i))
                write(*,"(' abs(Sdx) = ',e16.9)") abs(sdx(i))
                write(*,"('atan(Sdx) = ',e16.9)") atan2(aimag(sdx(i)),real(sdx(i)))

                write(*,"('     Sdy  = (',e16.9,') + i(',e16.9,')')") real(sdy(i)),aimag(sdy(i))
                write(*,"(' abs(Sdy) = ',e16.9)") abs(sdy(i))
                write(*,"('atan(Sdy) = ',e16.9)") atan2(aimag(sdy(i)),real(sdy(i)))

                write(*,"('     Sdz  = (',e16.9,') + i(',e16.9,')')") real(sdz(i)),aimag(sdz(i))
                write(*,"(' abs(Sdz) = ',e16.9)") abs(sdz(i))
                write(*,"('atan(Sdz) = ',e16.9)") atan2(aimag(sdz(i)),real(sdz(i)))

                ! Writing charge disturbance
                iw = 3000+(i-1)*7
                write(unit=iw+1,fmt="(7(e16.9,2x))") e,abs(chd(i)),real(chd(i)),aimag(chd(i)),atan2(aimag(chd(i)),real(chd(i))),real(chd(i))/abs(chd(i)),aimag(chd(i))/abs(chd(i))
                ! Writing x-component spin disturbance
                write(unit=iw+2,fmt="(7(e16.9,2x))") e,abs(sdx(i)),real(sdx(i)),aimag(sdx(i)),atan2(aimag(sdx(i)),real(sdx(i))),real(sdx(i))/abs(sdx(i)),aimag(sdx(i))/abs(sdx(i))
                ! Writing y-component spin disturbance
                write(unit=iw+3,fmt="(7(e16.9,2x))") e,abs(sdy(i)),real(sdy(i)),aimag(sdy(i)),atan2(aimag(sdy(i)),real(sdy(i))),real(sdy(i))/abs(sdy(i)),aimag(sdy(i))/abs(sdy(i))
                ! Writing z-component spin disturbance
                write(unit=iw+4,fmt="(7(e16.9,2x))") e,abs(sdz(i)),real(sdz(i)),aimag(sdz(i)),atan2(aimag(sdz(i)),real(sdz(i))),real(sdz(i))/abs(sdz(i)),aimag(sdz(i))/abs(sdz(i))

                write(*,"(' ################ Orbital disturbances: ################')")

                write(*,"('     Ldx  = (',e16.9,') + i(',e16.9,')')") real(ldx(i)),aimag(ldx(i))
                write(*,"(' abs(Ldx) = ',e16.9)") abs(ldx(i))
                write(*,"('atan(Ldx) = ',e16.9)") atan2(aimag(ldx(i)),real(ldx(i)))

                write(*,"('     Ldy  = (',e16.9,') + i(',e16.9,')')") real(ldy(i)),aimag(ldy(i))
                write(*,"(' abs(Ldy) = ',e16.9)") abs(ldy(i))
                write(*,"('atan(Ldy) = ',e16.9)") atan2(aimag(ldy(i)),real(ldy(i)))

                write(*,"('     Ldz  = (',e16.9,') + i(',e16.9,')')") real(ldz(i)),aimag(ldz(i))
                write(*,"(' abs(Ldz) = ',e16.9)") abs(ldz(i))
                write(*,"('atan(Ldz) = ',e16.9)") atan2(aimag(ldz(i)),real(ldz(i)))

                ! Writing x-component orbital disturbance
                write(unit=iw+5,fmt="(7(e16.9,2x))") e,abs(ldx(i)),real(ldx(i)),aimag(ldx(i)),atan2(aimag(ldx(i)),real(ldx(i))),real(ldx(i))/abs(ldx(i)),aimag(ldx(i))/abs(ldx(i))
                ! Writing y-component orbital disturbance
                write(unit=iw+6,fmt="(7(e16.9,2x))") e,abs(ldy(i)),real(ldy(i)),aimag(ldy(i)),atan2(aimag(ldy(i)),real(ldy(i))),real(ldy(i))/abs(ldy(i)),aimag(ldy(i))/abs(ldy(i))
                ! Writing z-component orbital disturbance
                write(unit=iw+7,fmt="(7(e16.9,2x))") e,abs(ldz(i)),real(ldz(i)),aimag(ldz(i)),atan2(aimag(ldz(i)),real(ldz(i))),real(ldz(i))/abs(ldz(i)),aimag(ldz(i))/abs(ldz(i))

              end do write_sd
            end do MPI_points_sd
            ! Closing files
            call openclosechifiles(2)
            call openclosesdfiles(2)

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
              write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
              call system ('rm stop')
              write(*,"(a)") "[main] ('stop' file deleted!)"
              call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
            end if
          else
            call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,2000,MPIComm_Col,ierr)
            call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,2100,MPIComm_Col,ierr)
            call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,2200,MPIComm_Col,ierr)
            call MPI_Send(sdx,Npl,MPI_DOUBLE_COMPLEX,0,2300,MPIComm_Col,ierr)
            call MPI_Send(sdy,Npl,MPI_DOUBLE_COMPLEX,0,2400,MPIComm_Col,ierr)
            call MPI_Send(sdz,Npl,MPI_DOUBLE_COMPLEX,0,2500,MPIComm_Col,ierr)
            call MPI_Send(chd,Npl,MPI_DOUBLE_COMPLEX,0,2600,MPIComm_Col,ierr)
            call MPI_Send(ldx,Npl,MPI_DOUBLE_COMPLEX,0,2700,MPIComm_Col,ierr)
            call MPI_Send(ldy,Npl,MPI_DOUBLE_COMPLEX,0,2800,MPIComm_Col,ierr)
            call MPI_Send(ldz,Npl,MPI_DOUBLE_COMPLEX,0,2900,MPIComm_Col,ierr)
          end if
        end if
      end do sd_energy_loop

      deallocate(chiorb_hf,prefactor,tchiorbiikl)
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            deallocate(rwork,eval,evecl,evecr,work)
#ifdef _JUQUEEN
            deallocate(dscale,rconde,rcondv)
#endif
          end if
        end if
        deallocate(templd,chiorb)
        deallocate(schi,schihf,sdx,sdy,sdz,chd,ldx,ldy,ldz)
      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (7)
      if(myrank_row.eq.0) then
        allocate( schi(Npl,Npl,4),schihf(Npl,Npl,4), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          if(myrank.eq.0) write(*,"('[main] Not enough memory for: schi,schihf')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        allocate( Ich(n0sc1:n0sc2,Npl),Isx(n0sc1:n0sc2,Npl),Isy(n0sc1:n0sc2,Npl),Isz(n0sc1:n0sc2,Npl),Ilx(n0sc1:n0sc2,Npl),Ily(n0sc1:n0sc2,Npl),Ilz(n0sc1:n0sc2,Npl), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          if(myrank.eq.0) write(*,"('[main] Not enough memory for: Ich,Isx,Isy,Isz,Ilx,Ily,Ilz')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        if(renorm) then
          allocate( rIch(n0sc1:n0sc2,Npl),rIsx(n0sc1:n0sc2,Npl),rIsy(n0sc1:n0sc2,Npl),rIsz(n0sc1:n0sc2,Npl),rIlx(n0sc1:n0sc2,Npl),rIly(n0sc1:n0sc2,Npl),rIlz(n0sc1:n0sc2,Npl), STAT = AllocateStatus )
          if (AllocateStatus.ne.0) then
            if(myrank.eq.0) write(*,"('[main] Not enough memory for: rIch,rIsx,rIsy,rIsz,rIlx,rIly,rIlz')")
            call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
          end if
        end if
        allocate( chiorb(dim,dim), STAT = AllocateStatus  )
        if (AllocateStatus.ne.0) then
          write(*,"('[main] Not enough memory for: chiorb')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            lwork = 33*nmaglayers
            allocate( chimag(nmaglayers,nmaglayers),rwork(2*nmaglayers),eval(nmaglayers),evecl(1,nmaglayers),evecr(nmaglayers,nmaglayers),work(lwork) )
#ifdef _JUQUEEN
            allocate( dscale(nmaglayers),rconde(nmaglayers),rcondv(nmaglayers) )
#endif
          end if
          write(*,"('CALCULATING CURRENTS AND LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
          ! Creating files and writing headers
          if(index(runoptions,"addresults").eq.0) then
            call openclosechifiles(0)
            call openclosescfiles(0)
          end if
        end if
      end if
      allocate( chiorb_hf(dim,dim),prefactor(dim,dim),ttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lxttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lyttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lzttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(*,"('[main] Not enough memory for: chiorb_hf,prefactor,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if

      ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
      call OAM_curr_hopping_times_L()

      sc_energy_loop: do count=1,MPIsteps
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',e10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e,myrank_col
        end if

        if(myrank.eq.0) write(*,"('[main] Calculating pre-factor to use in current calculation ')")
        call eintshechi(e,chiorb_hf)

        ! Broadcast chiorb_hf to all processors of the same row
        call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)

        ! prefactor = (1 + chi_hf*Umat)^-1
        call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zero,prefactor,dim) !prefactor = chi_hf*Umat
        prefactor = identt + prefactor
        call invers(prefactor,dim)
        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Calculated prefactor after: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        ! Start parallelized processes to calculate chiorb_hf and chiorbi0_hf for energy e
        call eintsheprllsc(e,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)

        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        if(myrank_row.eq.0) then
          ! Calculating the full matrix of RPA and HF susceptibilities for energy e
          call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hf,dim,zero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf

          schi = zero
          schihf = zero
          calculate_susceptibility_sc: do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
            ! Calculating RPA and HF susceptibilities
            do nu=5,9 ; do mu=5,9
              schi(i,j,sigma) = schi(i,j,sigma) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
              schihf(i,j,sigma) = schihf(i,j,sigma) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
            end do ; end do
          end do ; end do ; end do calculate_susceptibility_sc

          Ich = zero
          Isx = zero
          Isy = zero
          Isz = zero
          Ilx = zero
          Ily = zero
          Ilz = zero
          ! Calculating spin and charge current for each neighbor
          plane_loop_calculate_sc: do i=1,Npl
            neighbor_loop_calculate_sc: do neighbor=n0sc1,n0sc2
              ! Charge current
              Ich(neighbor,i) = Ich(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)+ttchiorbiikl  (neighbor,sigmai2i(3,i),2)+ttchiorbiikl  (neighbor,sigmai2i(3,i),3))
              ! Spin currents
              Isx(neighbor,i) = Isx(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)+ttchiorbiikl  (neighbor,sigmai2i(4,i),2)+ttchiorbiikl  (neighbor,sigmai2i(4,i),3))
              Isy(neighbor,i) = Isy(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)-ttchiorbiikl  (neighbor,sigmai2i(4,i),2)-ttchiorbiikl  (neighbor,sigmai2i(4,i),3))
              Isz(neighbor,i) = Isz(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)-ttchiorbiikl  (neighbor,sigmai2i(3,i),2)-ttchiorbiikl  (neighbor,sigmai2i(3,i),3))
              ! Orbital Angular Momentum currents
              Ilx(neighbor,i) = Ilx(neighbor,i) + (Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3))
              Ily(neighbor,i) = Ily(neighbor,i) + (Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3))
              Ilz(neighbor,i) = Ilz(neighbor,i) + (Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3))
            end do neighbor_loop_calculate_sc
          end do plane_loop_calculate_sc
  !         Ich = Ich
          Isx = -0.5d0*Isx
          Isy = -0.5d0*Isy/zi
          Isz = -0.5d0*Isz
  !         Ilx = Ilx
  !         Ily = Ily
  !         Ilz = Ilz

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            ! Opening files for writing
            call openclosechifiles(1)
            call openclosescfiles(1)
            MPI_points_sc: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3000,MPIComm_Col,stat,ierr)
                call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3100,MPIComm_Col,stat,ierr)
                call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3200,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ich,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3300,MPIComm_Col,stat,ierr)
                call MPI_Recv(Isx,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3400,MPIComm_Col,stat,ierr)
                call MPI_Recv(Isy,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3500,MPIComm_Col,stat,ierr)
                call MPI_Recv(Isz,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3600,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ilx,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3700,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ily,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3800,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ilz,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),3900,MPIComm_Col,stat,ierr)
              end if

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              write_susceptibility_sc: do sigma=1,4 ; do j=1,Npl ;  do i=1,Npl
                iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i

                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schi(i,j,sigma)),aimag(schi(i,j,sigma)),abs(schi(i,j,sigma))

                iw = iw+1000
                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schihf(i,j,sigma)),aimag(schihf(i,j,sigma)),abs(schihf(i,j,sigma))
              end do ; end do ; end do write_susceptibility_sc

              ! DIAGONALIZING SUSCEPTIBILITY
              if(nmaglayers.gt.1) then
                do i=1,nmaglayers ; do j=1,nmaglayers
                  chimag(i,j) = schi(mmlayermag(i)-1,mmlayermag(j)-1,1)
                end do ; end do

#ifdef _JUQUEEN
                call zgeevx('N','N','V','N',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,ilo,ihi,dscale,abnrm,rconde,rcondv,work,lwork,rwork,ifail)
#else
                call zgeev('N','V',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,work,lwork,rwork,ifail)
#endif

                if(ifail.ne.0) then
                  write(*,*) '[main] Problem with diagonalization. ifail = ',ifail
                  stop
        !         else
        !           write(*,*) ' optimal lwork = ',work(1),' lwork = ',lwork
                end if
                write(varm,fmt="(a,i0,a)") '(',2*nmaglayers+1,'(e16.9,2x))'
                write(unit=1990,fmt=varm) e,(real(eval(i)),aimag(eval(i)),i=1,nmaglayers)
                do i=1,nmaglayers
                  write(unit=1990+i,fmt=varm) e,(real(evecr(j,i)),aimag(evecr(j,i)),j=1,nmaglayers)
                end do
              end if

              ! Renormalizing currents by the charge current in plane 1
              if(renorm) then
                ! Obtaining current for renormalization
                Icabs  = abs(Ich(renormnb,1))

                rIch = Ich/Icabs
                rIsx = Isx/Icabs
                rIsy = Isy/Icabs
                rIsz = Isz/Icabs
                rIlx = Ilx/Icabs
                rIly = Ily/Icabs
                rIlz = Ilz/Icabs
              end if

              ! WRITING CURRENTS
              plane_loop_write_sc: do i=1,Npl
                write(*,"('|--------------- Plane: ',i0,' , Energy = ',e11.4,' ---------------|')") i,e
                write(*,"(' #################  Susceptibility:  #################')")
                write(*,"(' Chi+- = (',e16.9,') + i(',e16.9,')')") real(schi(i,i,1)),aimag(schi(i,i,1))

                neighbor_loop_write_sc: do neighbor=n0sc1,n0sc2
                  write(*,"('   ***************** Neighbor: ',i0,'  *****************')") neighbor

                  write(*,"('  Charge current:')")
                  write(*,"('     Ich = (',e16.9,') + i(',e16.9,')')") real(Ich(neighbor,i)),aimag(Ich(neighbor,i))
                  write(*,"(' abs(Ich) = ',e16.9)") abs(Ich(neighbor,i))
                  write(*,"('atan(Ich) = ',e16.9)") atan2(aimag(Ich(neighbor,i)),real(Ich(neighbor,i)))

                  write(*,"('  Spin currents:')")
                  write(*,"('     Isx  = (',e16.9,') + i(',e16.9,')')") real(Isx(neighbor,i)),aimag(Isx(neighbor,i))
                  write(*,"(' abs(Isx) = ',e16.9)") abs(Isx(neighbor,i))
                  write(*,"('atan(Isx) = ',e16.9)") atan2(aimag(Isx(neighbor,i)),real(Isx(neighbor,i)))

                  write(*,"('     Isy  = (',e16.9,') + i(',e16.9,')')") real(Isy(neighbor,i)),aimag(Isy(neighbor,i))
                  write(*,"(' abs(Isy) = ',e16.9)") abs(Isy(neighbor,i))
                  write(*,"('atan(Isy) = ',e16.9)") atan2(aimag(Isy(neighbor,i)),real(Isy(neighbor,i)))

                  write(*,"('     Isz  = (',e16.9,') + i(',e16.9,')')") real(Isz(neighbor,i)),aimag(Isz(neighbor,i))
                  write(*,"(' abs(Isz) = ',e16.9)") abs(Isz(neighbor,i))
                  write(*,"('atan(Isz) = ',e16.9)") atan2(aimag(Isz(neighbor,i)),real(Isz(neighbor,i)))

                  write(*,"('  Orbital Angular Momentum currents:')")
                  write(*,"('     Ilx  = (',e16.9,') + i(',e16.9,')')") real(Ilx(neighbor,i)),aimag(Ilx(neighbor,i))
                  write(*,"(' abs(Ilx) = ',e16.9)") abs(Ilx(neighbor,i))
                  write(*,"('atan(Ilx) = ',e16.9)") atan2(aimag(Ilx(neighbor,i)),real(Ilx(neighbor,i)))

                  write(*,"('     Ily  = (',e16.9,') + i(',e16.9,')')") real(Ily(neighbor,i)),aimag(Ily(neighbor,i))
                  write(*,"(' abs(Ily) = ',e16.9)") abs(Ily(neighbor,i))
                  write(*,"('atan(Ily) = ',e16.9)") atan2(aimag(Ily(neighbor,i)),real(Ily(neighbor,i)))

                  write(*,"('     Ilz  = (',e16.9,') + i(',e16.9,')')") real(Ilz(neighbor,i)),aimag(Ilz(neighbor,i))
                  write(*,"(' abs(Ilz) = ',e16.9)") abs(Ilz(neighbor,i))
                  write(*,"('atan(Ilz) = ',e16.9)") atan2(aimag(Ilz(neighbor,i)),real(Ilz(neighbor,i)))

                  ! CHARGE CURRENT
                  ! Writing charge current
                  iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7
                  write(unit=iw+1,fmt="(7(e16.9,2x))") e,abs(Ich(neighbor,i)),real(Ich(neighbor,i)),aimag(Ich(neighbor,i)),atan2(aimag(Ich(neighbor,i)),real(Ich(neighbor,i))),real(Ich(neighbor,i))/abs(Ich(neighbor,i)),aimag(Ich(neighbor,i))/abs(Ich(neighbor,i))

                  ! SPIN CURRENTS
                  ! Writing x-component spin current
                  write(unit=iw+2,fmt="(7(e16.9,2x))") e,abs(Isx(neighbor,i)),real(Isx(neighbor,i)),aimag(Isx(neighbor,i)),atan2(aimag(Isx(neighbor,i)),real(Isx(neighbor,i))),real(Isx(neighbor,i))/abs(Isx(neighbor,i)),aimag(Isx(neighbor,i))/abs(Isx(neighbor,i))
                  ! Writing y-component spin current
                  write(unit=iw+3,fmt="(7(e16.9,2x))") e,abs(Isy(neighbor,i)),real(Isy(neighbor,i)),aimag(Isy(neighbor,i)),atan2(aimag(Isy(neighbor,i)),real(Isy(neighbor,i))),real(Isy(neighbor,i))/abs(Isy(neighbor,i)),aimag(Isy(neighbor,i))/abs(Isy(neighbor,i))
                  ! Writing z-component spin current
                  write(unit=iw+4,fmt="(7(e16.9,2x))") e,abs(Isz(neighbor,i)),real(Isz(neighbor,i)),aimag(Isz(neighbor,i)),atan2(aimag(Isz(neighbor,i)),real(Isz(neighbor,i))),real(Isz(neighbor,i))/abs(Isz(neighbor,i)),aimag(Isz(neighbor,i))/abs(Isz(neighbor,i))

                  ! ORBITAL ANGULAR MOMENTUM CURRENTS
                  ! Writing x-component orbital angular momentum current
                  write(unit=iw+5,fmt="(7(e16.9,2x))") e,abs(Ilx(neighbor,i)),real(Ilx(neighbor,i)),aimag(Ilx(neighbor,i)),atan2(aimag(Ilx(neighbor,i)),real(Ilx(neighbor,i))),real(Ilx(neighbor,i))/abs(Ilx(neighbor,i)),aimag(Ilx(neighbor,i))/abs(Ilx(neighbor,i))
                  ! Writing y-component orbital angular momentum current
                  write(unit=iw+6,fmt="(7(e16.9,2x))") e,abs(Ily(neighbor,i)),real(Ily(neighbor,i)),aimag(Ily(neighbor,i)),atan2(aimag(Ily(neighbor,i)),real(Ily(neighbor,i))),real(Ily(neighbor,i))/abs(Ily(neighbor,i)),aimag(Ily(neighbor,i))/abs(Ily(neighbor,i))
                  ! Writing z-component orbital angular momentum current
                  write(unit=iw+7,fmt="(7(e16.9,2x))") e,abs(Ilz(neighbor,i)),real(Ilz(neighbor,i)),aimag(Ilz(neighbor,i)),atan2(aimag(Ilz(neighbor,i)),real(Ilz(neighbor,i))),real(Ilz(neighbor,i))/abs(Ilz(neighbor,i)),aimag(Ilz(neighbor,i))/abs(Ilz(neighbor,i))

                  ! Writing renormalized currents
                  if(renorm) then
                    ! CHARGE CURRENT
                    ! Writing renormalized charge current
                    write(unit=iw+1001,fmt="(7(e16.9,2x))") e,abs(rIch(neighbor,i)),real(rIch(neighbor,i)),aimag(rIch(neighbor,i)),atan2(aimag(rIch(neighbor,i)),real(rIch(neighbor,i))),real(rIch(neighbor,i))/abs(rIch(neighbor,i)),aimag(rIch(neighbor,i))/abs(rIch(neighbor,i))

                    ! SPIN CURRENTS
                    ! Writing renormalized x-component spin current
                    write(unit=iw+1002,fmt="(7(e16.9,2x))") e,abs(rIsx(neighbor,i)),real(rIsx(neighbor,i)),aimag(rIsx(neighbor,i)),atan2(aimag(rIsx(neighbor,i)),real(rIsx(neighbor,i))),real(rIsx(neighbor,i))/abs(rIsx(neighbor,i)),aimag(rIsx(neighbor,i))/abs(rIsx(neighbor,i))
                    ! Writing renormalized y-component spin current
                    write(unit=iw+1003,fmt="(7(e16.9,2x))") e,abs(rIsy(neighbor,i)),real(rIsy(neighbor,i)),aimag(rIsy(neighbor,i)),atan2(aimag(rIsy(neighbor,i)),real(rIsy(neighbor,i))),real(rIsy(neighbor,i))/abs(rIsy(neighbor,i)),aimag(rIsy(neighbor,i))/abs(rIsy(neighbor,i))
                    ! Writing renormalized z-component spin current
                    write(unit=iw+1004,fmt="(7(e16.9,2x))") e,abs(rIsz(neighbor,i)),real(rIsz(neighbor,i)),aimag(rIsz(neighbor,i)),atan2(aimag(rIsz(neighbor,i)),real(rIsz(neighbor,i))),real(rIsz(neighbor,i))/abs(rIsz(neighbor,i)),aimag(rIsz(neighbor,i))/abs(rIsz(neighbor,i))

                    ! ORBITAL ANGULAR MOMENTUM CURRENTS
                    ! Writing x-component orbital angular momentum current
                    write(unit=iw+1005,fmt="(7(e16.9,2x))") e,abs(rIlx(neighbor,i)),real(rIlx(neighbor,i)),aimag(rIlx(neighbor,i)),atan2(aimag(rIlx(neighbor,i)),real(rIlx(neighbor,i))),real(rIlx(neighbor,i))/abs(rIlx(neighbor,i)),aimag(rIlx(neighbor,i))/abs(rIlx(neighbor,i))
                    ! Writing y-component orbital angular momentum current
                    write(unit=iw+1006,fmt="(7(e16.9,2x))") e,abs(rIly(neighbor,i)),real(rIly(neighbor,i)),aimag(rIly(neighbor,i)),atan2(aimag(rIly(neighbor,i)),real(rIly(neighbor,i))),real(rIly(neighbor,i))/abs(rIly(neighbor,i)),aimag(rIly(neighbor,i))/abs(rIly(neighbor,i))
                    ! Writing z-component orbital angular momentum current
                    write(unit=iw+1007,fmt="(7(e16.9,2x))") e,abs(rIlz(neighbor,i)),real(rIlz(neighbor,i)),aimag(rIlz(neighbor,i)),atan2(aimag(rIlz(neighbor,i)),real(rIlz(neighbor,i))),real(rIlz(neighbor,i))/abs(rIlz(neighbor,i)),aimag(rIlz(neighbor,i))/abs(rIlz(neighbor,i))
                  end if
                end do neighbor_loop_write_sc
              end do plane_loop_write_sc
            end do MPI_points_sc
            ! Closing files
            call openclosechifiles(2)
            call openclosescfiles(2)

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
              write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
              call system ('rm stop')
              write(*,"(a)") "[main] ('stop' file deleted!)"
              call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
            end if
          else
            call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,3000,MPIComm_Col,ierr)
            call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,3100,MPIComm_Col,ierr)
            call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,3200,MPIComm_Col,ierr)
            call MPI_Send(Ich,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3300,MPIComm_Col,ierr)
            call MPI_Send(Isx,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3400,MPIComm_Col,ierr)
            call MPI_Send(Isy,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3500,MPIComm_Col,ierr)
            call MPI_Send(Isz,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3600,MPIComm_Col,ierr)
            call MPI_Send(Ilx,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3700,MPIComm_Col,ierr)
            call MPI_Send(Ily,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3800,MPIComm_Col,ierr)
            call MPI_Send(Ilz,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,3900,MPIComm_Col,ierr)
          end if
        end if
      end do sc_energy_loop

      deallocate(chiorb_hf,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            deallocate(rwork,eval,evecl,evecr,work)
#ifdef _JUQUEEN
            deallocate(dscale,rconde,rcondv)
#endif
          end if
        end if
        deallocate(chiorb)
        deallocate(schi,schihf)
        if(renorm) deallocate(rIch,rIsx,rIsy,rIsz,rIlx,rIly,rIlz)
        deallocate(Ich,Isx,Isy,Isz,Ilx,Ily,Ilz)
      end if
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
    case (8)
      if(myrank_row.eq.0) then
        allocate( schi(Npl,Npl,4),schihf(Npl,Npl,4),sdx(Npl),sdy(Npl),sdz(Npl),chd(Npl),ldx(Npl),ldy(Npl),ldz(Npl), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          if(myrank.eq.0) write(*,"('[main] Not enough memory for: schi,schihf,sdx,sdy,sdz,chd,ldx,ldy,ldz')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        allocate( Ich(n0sc1:n0sc2,Npl),Isx(n0sc1:n0sc2,Npl),Isy(n0sc1:n0sc2,Npl),Isz(n0sc1:n0sc2,Npl),Ilx(n0sc1:n0sc2,Npl),Ily(n0sc1:n0sc2,Npl),Ilz(n0sc1:n0sc2,Npl), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          if(myrank.eq.0) write(*,"('[main] Not enough memory for: Ich,Isx,Isy,Isz,Ilx,Ily,Ilz')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        if(renorm) then
          allocate( rsdx(Npl),rsdy(Npl),rsdz(Npl),rchd(Npl),rldx(Npl),rldy(Npl),rldz(Npl), STAT = AllocateStatus )
          if (AllocateStatus.ne.0) then
            if(myrank.eq.0) write(*,"('[main] Not enough memory for: rsdx,rsdy,rsdz,rchd,rldx,rldy,rldz')")
            call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
          end if
          allocate( rIch(n0sc1:n0sc2,Npl),rIsx(n0sc1:n0sc2,Npl),rIsy(n0sc1:n0sc2,Npl),rIsz(n0sc1:n0sc2,Npl),rIlx(n0sc1:n0sc2,Npl),rIly(n0sc1:n0sc2,Npl),rIlz(n0sc1:n0sc2,Npl), STAT = AllocateStatus )
          if (AllocateStatus.ne.0) then
            if(myrank.eq.0) write(*,"('[main] Not enough memory for: rIch,rIsx,rIsy,rIsz,rIlx,rIly,rIlz')")
            call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
          end if
        end if
        allocate( templd(Npl,9,9),chiorb(dim,dim), STAT = AllocateStatus  )
        if (AllocateStatus.ne.0) then
          write(*,"('[main] Not enough memory for: templd,chiorb')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            lwork = 33*nmaglayers
            allocate( chimag(nmaglayers,nmaglayers),rwork(2*nmaglayers),eval(nmaglayers),evecl(1,nmaglayers),evecr(nmaglayers,nmaglayers),work(lwork) )
#ifdef _JUQUEEN
            allocate( dscale(nmaglayers),rconde(nmaglayers),rcondv(nmaglayers) )
#endif
          end if
          write(*,"('CALCULATING PARALLEL CURRENTS, DISTURBANCES AND LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
          ! Creating files and writing headers
          if(index(runoptions,"addresults").eq.0) then
            call openclosechifiles(0)
            call openclosesdfiles(0)
            call openclosescfiles(0)
          end if
        end if
      end if
      allocate( chiorb_hf(dim,dim),prefactor(dim,dim),tchiorbiikl(dim,4),ttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lxttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lyttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lzttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4), STAT = AllocateStatus  )
      if (AllocateStatus.ne.0) then
        write(*,"('[main] Not enough memory for: chiorb_hf,prefactor,tchiorbiikl,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if

      ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
      call OAM_curr_hopping_times_L()

      all_energy_loop: do count=1,MPIsteps
        e = emin + deltae*myrank_col + MPIdelta*(count-1)
        if(myrank_row.eq.0) then
          write(*,"(i0,' of ',i0,' points',', e = ',e10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e,myrank_col
        end if

        if(myrank.eq.0) write(*,"('[main] Calculating pre-factor to use in currents and disturbances calculation ')")
        call eintshechi(e,chiorb_hf)

        ! Broadcast chiorb_hf to all processors of the same row
        call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPIComm_Row,ierr)

        ! prefactor = (1 + chi_hf*Umat)^-1
        call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zero,prefactor,dim) !prefactor = chi_hf*Umat
        prefactor = identt + prefactor
        call invers(prefactor,dim)
        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Calculated prefactor after: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        ! Start parallelized processes to calculate disturbances and currents for energy e
        call eintshe(e,tchiorbiikl,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)

        if(myrank.eq.0) then
          elapsed_time = MPI_Wtime() - start_program
          write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
        end if

        if(myrank_row.eq.0) then
          ! Calculating the full matrix of RPA and HF susceptibilities for energy e
          call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hf,dim,zero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf

          schi = zero
          schihf = zero
          calculate_susceptibility_all: do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
            ! Calculating RPA and HF susceptibilities
            do nu=5,9 ; do mu=5,9
              schi(i,j,sigma) = schi(i,j,sigma) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
              schihf(i,j,sigma) = schihf(i,j,sigma) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(4,j,nu,nu)) ! +- , up- , down- , --
            end do ; end do
          end do ; end do ; end do calculate_susceptibility_all

          sdx = zero
          sdy = zero
          sdz = zero
          chd = zero
          ldx = zero
          ldy = zero
          ldz = zero
          Ich = zero
          Isx = zero
          Isy = zero
          Isz = zero
          Ilx = zero
          Ily = zero
          Ilz = zero
          plane_loop_calculate_all: do i=1,Npl
            ! Spin and charge disturbances
            do mu=1,9
              sdx(i) = sdx(i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
              sdy(i) = sdy(i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
              sdz(i) = sdz(i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
              chd(i) = chd(i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
            end do

            ! Orbital angular momentum disturbance
            do nu=1,9; do mu=1,9
              templd(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
              ldx(i) = ldx(i) + lxp(mu,nu)*templd(i,mu,nu)
              ldy(i) = ldy(i) + lyp(mu,nu)*templd(i,mu,nu)
              ldz(i) = ldz(i) + lzp(mu,nu)*templd(i,mu,nu)
            end do; end do

            ! Calculating spin and charge current for each neighbor
            neighbor_loop_calculate_all: do neighbor=n0sc1,n0sc2
              ! Charge current
              Ich(neighbor,i) = Ich(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)+ttchiorbiikl  (neighbor,sigmai2i(3,i),2)+ttchiorbiikl  (neighbor,sigmai2i(3,i),3))
              ! Spin currents
              Isx(neighbor,i) = Isx(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)+ttchiorbiikl  (neighbor,sigmai2i(4,i),2)+ttchiorbiikl  (neighbor,sigmai2i(4,i),3))
              Isy(neighbor,i) = Isy(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)-ttchiorbiikl  (neighbor,sigmai2i(4,i),2)-ttchiorbiikl  (neighbor,sigmai2i(4,i),3))
              Isz(neighbor,i) = Isz(neighbor,i) + (ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)-ttchiorbiikl  (neighbor,sigmai2i(3,i),2)-ttchiorbiikl  (neighbor,sigmai2i(3,i),3))
              ! Orbital Angular Momentum currents
              Ilx(neighbor,i) = Ilx(neighbor,i) + (Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3))
              Ily(neighbor,i) = Ily(neighbor,i) + (Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3))
              Ilz(neighbor,i) = Ilz(neighbor,i) + (Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3))
            end do neighbor_loop_calculate_all
          end do plane_loop_calculate_all
          sdx = 0.5d0*sdx
          sdy = 0.5d0*sdy/zi
          sdz = 0.5d0*sdz
          Isx = -0.5d0*Isx
          Isy = -0.5d0*Isy/zi
          Isz = -0.5d0*Isz

          ! Sending results to myrank_row = myrank_col = 0 and writing on file
          if(myrank_col.eq.0) then
            ! Opening files for writing
            call openclosechifiles(1)
            call openclosesdfiles(1)
            call openclosescfiles(1)
            MPI_points_all: do mcount=1,MPIpts
              if (mcount.ne.1) then
                call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,4000,MPIComm_Col,stat,ierr)
                call MPI_Recv(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4010,MPIComm_Col,stat,ierr)
                call MPI_Recv(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4020,MPIComm_Col,stat,ierr)
                call MPI_Recv(sdx,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4030,MPIComm_Col,stat,ierr)
                call MPI_Recv(sdy,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4040,MPIComm_Col,stat,ierr)
                call MPI_Recv(sdz,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4050,MPIComm_Col,stat,ierr)
                call MPI_Recv(chd,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4060,MPIComm_Col,stat,ierr)
                call MPI_Recv(ldx,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4070,MPIComm_Col,stat,ierr)
                call MPI_Recv(ldy,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4080,MPIComm_Col,stat,ierr)
                call MPI_Recv(ldz,Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4090,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ich,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4100,MPIComm_Col,stat,ierr)
                call MPI_Recv(Isx,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4110,MPIComm_Col,stat,ierr)
                call MPI_Recv(Isy,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4120,MPIComm_Col,stat,ierr)
                call MPI_Recv(Isz,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4130,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ilx,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4140,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ily,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4150,MPIComm_Col,stat,ierr)
                call MPI_Recv(Ilz,n0sc*Npl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),4160,MPIComm_Col,stat,ierr)
              end if

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              write_susceptibility_all: do sigma=1,4 ; do j=1,Npl ;  do i=1,Npl
                iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i

                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schi(i,j,sigma)),aimag(schi(i,j,sigma)),abs(schi(i,j,sigma))

                iw = iw+1000
                write(unit=iw,fmt="(4(e16.9,2x))") e,real(schihf(i,j,sigma)),aimag(schihf(i,j,sigma)),abs(schihf(i,j,sigma))
              end do ; end do ; end do write_susceptibility_all

              ! DIAGONALIZING SUSCEPTIBILITY
              if(nmaglayers.gt.1) then
                do i=1,nmaglayers ; do j=1,nmaglayers
                  chimag(i,j) = schi(mmlayermag(i)-1,mmlayermag(j)-1,1)
                end do ; end do

#ifdef _JUQUEEN
                call zgeevx('N','N','V','N',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,ilo,ihi,dscale,abnrm,rconde,rcondv,work,lwork,rwork,ifail)
#else
                call zgeev('N','V',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,work,lwork,rwork,ifail)
#endif

                if(ifail.ne.0) then
                  write(*,*) '[main] Problem with diagonalization. ifail = ',ifail
                  stop
        !         else
        !           write(*,*) ' optimal lwork = ',work(1),' lwork = ',lwork
                end if
                write(varm,fmt="(a,i0,a)") '(',2*nmaglayers+1,'(e16.9,2x))'
                write(unit=1990,fmt=varm) e,(real(eval(i)),aimag(eval(i)),i=1,nmaglayers)
                do i=1,nmaglayers
                  write(unit=1990+i,fmt=varm) e,(real(evecr(j,i)),aimag(evecr(j,i)),j=1,nmaglayers)
                end do
              end if

              ! Renormalizing disturbances and currents by the charge current in plane 1
              if(renorm) then
                ! Obtaining current for renormalization
                Icabs  = abs(Ich(renormnb,1))

                rsdx = sdx/Icabs
                rsdy = sdy/Icabs
                rsdz = sdz/Icabs
                rchd = chd/Icabs
                rldx = ldx/Icabs
                rldy = ldy/Icabs
                rldz = ldz/Icabs

                rIch = Ich/Icabs
                rIsx = Isx/Icabs
                rIsy = Isy/Icabs
                rIsz = Isz/Icabs
                rIlx = Ilx/Icabs
                rIly = Ily/Icabs
                rIlz = Ilz/Icabs
              end if

              ! WRITING DISTURBANCES AND CURRENTS
              plane_loop_write_all: do i=1,Npl
                write(*,"('|--------------- Plane: ',i0,' , Energy = ',e11.4,' ---------------|')") i,e
                write(*,"(' #################  Susceptibility:  #################')")
                write(*,"(' Chi+- = (',e16.9,') + i(',e16.9,')')") real(schi(i,i,1)),aimag(schi(i,i,1))

                ! Writing Spin, Charge and Orbital disturbances
                write(*,"(' ################# Charge disturbance: #################')")

                write(*,"('     Cd  = (',e16.9,') + i(',e16.9,')')") real(chd(i)),aimag(chd(i))
                write(*,"(' abs(Cd) = ',e16.9)") abs(chd(i))
                write(*,"('atan(Cd) = ',e16.9)") atan2(aimag(chd(i)),real(chd(i)))

                write(*,"(' ################# Spin disturbances:  #################')")

                write(*,"('     Sdx  = (',e16.9,') + i(',e16.9,')')") real(sdx(i)),aimag(sdx(i))
                write(*,"(' abs(Sdx) = ',e16.9)") abs(sdx(i))
                write(*,"('atan(Sdx) = ',e16.9)") atan2(aimag(sdx(i)),real(sdx(i)))

                write(*,"('     Sdy  = (',e16.9,') + i(',e16.9,')')") real(sdy(i)),aimag(sdy(i))
                write(*,"(' abs(Sdy) = ',e16.9)") abs(sdy(i))
                write(*,"('atan(Sdy) = ',e16.9)") atan2(aimag(sdy(i)),real(sdy(i)))

                write(*,"('     Sdz  = (',e16.9,') + i(',e16.9,')')") real(sdz(i)),aimag(sdz(i))
                write(*,"(' abs(Sdz) = ',e16.9)") abs(sdz(i))
                write(*,"('atan(Sdz) = ',e16.9)") atan2(aimag(sdz(i)),real(sdz(i)))

                ! Writing charge disturbance
                iw = 3000+(i-1)*7
                write(unit=iw+1,fmt="(7(e16.9,2x))") e,abs(chd(i)),real(chd(i)),aimag(chd(i)),atan2(aimag(chd(i)),real(chd(i))),real(chd(i))/abs(chd(i)),aimag(chd(i))/abs(chd(i))
                ! Writing x-component spin disturbance
                write(unit=iw+2,fmt="(7(e16.9,2x))") e,abs(sdx(i)),real(sdx(i)),aimag(sdx(i)),atan2(aimag(sdx(i)),real(sdx(i))),real(sdx(i))/abs(sdx(i)),aimag(sdx(i))/abs(sdx(i))
                ! Writing y-component spin disturbance
                write(unit=iw+3,fmt="(7(e16.9,2x))") e,abs(sdy(i)),real(sdy(i)),aimag(sdy(i)),atan2(aimag(sdy(i)),real(sdy(i))),real(sdy(i))/abs(sdy(i)),aimag(sdy(i))/abs(sdy(i))
                ! Writing z-component spin disturbance
                write(unit=iw+4,fmt="(7(e16.9,2x))") e,abs(sdz(i)),real(sdz(i)),aimag(sdz(i)),atan2(aimag(sdz(i)),real(sdz(i))),real(sdz(i))/abs(sdz(i)),aimag(sdz(i))/abs(sdz(i))

                write(*,"(' ################ Orbital disturbances: ################')")

                write(*,"('     Ldx  = (',e16.9,') + i(',e16.9,')')") real(ldx(i)),aimag(ldx(i))
                write(*,"(' abs(Ldx) = ',e16.9)") abs(ldx(i))
                write(*,"('atan(Ldx) = ',e16.9)") atan2(aimag(ldx(i)),real(ldx(i)))

                write(*,"('     Ldy  = (',e16.9,') + i(',e16.9,')')") real(ldy(i)),aimag(ldy(i))
                write(*,"(' abs(Ldy) = ',e16.9)") abs(ldy(i))
                write(*,"('atan(Ldy) = ',e16.9)") atan2(aimag(ldy(i)),real(ldy(i)))

                write(*,"('     Ldz  = (',e16.9,') + i(',e16.9,')')") real(ldz(i)),aimag(ldz(i))
                write(*,"(' abs(Ldz) = ',e16.9)") abs(ldz(i))
                write(*,"('atan(Ldz) = ',e16.9)") atan2(aimag(ldz(i)),real(ldz(i)))

                ! Writing x-component orbital disturbance
                write(unit=iw+5,fmt="(7(e16.9,2x))") e,abs(ldx(i)),real(ldx(i)),aimag(ldx(i)),atan2(aimag(ldx(i)),real(ldx(i))),real(ldx(i))/abs(ldx(i)),aimag(ldx(i))/abs(ldx(i))
                ! Writing y-component orbital disturbance
                write(unit=iw+6,fmt="(7(e16.9,2x))") e,abs(ldy(i)),real(ldy(i)),aimag(ldy(i)),atan2(aimag(ldy(i)),real(ldy(i))),real(ldy(i))/abs(ldy(i)),aimag(ldy(i))/abs(ldy(i))
                ! Writing z-component orbital disturbance
                write(unit=iw+7,fmt="(7(e16.9,2x))") e,abs(ldz(i)),real(ldz(i)),aimag(ldz(i)),atan2(aimag(ldz(i)),real(ldz(i))),real(ldz(i))/abs(ldz(i)),aimag(ldz(i))/abs(ldz(i))

                ! Writing renormalized disturbances
                if(renorm) then
                  ! Writing renormalized charge disturbance
                  write(unit=iw+1001,fmt="(7(e16.9,2x))") e,abs(rchd(i)),real(rchd(i)),aimag(rchd(i)),atan2(aimag(rchd(i)),real(rchd(i))),real(rchd(i))/abs(rchd(i)),aimag(rchd(i))/abs(rchd(i))
                  ! Writing renormalized x-component spin disturbance
                  write(unit=iw+1002,fmt="(7(e16.9,2x))") e,abs(rsdx(i)),real(rsdx(i)),aimag(rsdx(i)),atan2(aimag(rsdx(i)),real(rsdx(i))),real(rsdx(i))/abs(rsdx(i)),aimag(rsdx(i))/abs(rsdx(i))
                  ! Writing renormalized y-component spin disturbance
                  write(unit=iw+1003,fmt="(7(e16.9,2x))") e,abs(rsdy(i)),real(rsdy(i)),aimag(rsdy(i)),atan2(aimag(rsdy(i)),real(rsdy(i))),real(rsdy(i))/abs(rsdy(i)),aimag(rsdy(i))/abs(rsdy(i))
                  ! Writing renormalized z-component spin disturbance
                  write(unit=iw+1004,fmt="(7(e16.9,2x))") e,abs(rsdz(i)),real(rsdz(i)),aimag(rsdz(i)),atan2(aimag(rsdz(i)),real(rsdz(i))),real(rsdz(i))/abs(rsdz(i)),aimag(rsdz(i))/abs(rsdz(i))

                  ! Writing renormalized x-component orbital disturbance
                  write(unit=iw+1005,fmt="(7(e16.9,2x))") e,abs(rldx(i)),real(rldx(i)),aimag(rldx(i)),atan2(aimag(rldx(i)),real(rldx(i))),real(rldx(i))/abs(rldx(i)),aimag(rldx(i))/abs(rldx(i))
                  ! Writing renormalized y-component orbital disturbance
                  write(unit=iw+1006,fmt="(7(e16.9,2x))") e,abs(rldy(i)),real(rldy(i)),aimag(rldy(i)),atan2(aimag(rldy(i)),real(rldy(i))),real(rldy(i))/abs(rldy(i)),aimag(rldy(i))/abs(rldy(i))
                  ! Writing renormalized z-component orbital disturbance
                  write(unit=iw+1007,fmt="(7(e16.9,2x))") e,abs(rldz(i)),real(rldz(i)),aimag(rldz(i)),atan2(aimag(rldz(i)),real(rldz(i))),real(rldz(i))/abs(rldz(i)),aimag(rldz(i))/abs(rldz(i))
                end if

                neighbor_loop_write_all: do neighbor=n0sc1,n0sc2
                  write(*,"('   ***************** Neighbor: ',i0,'  *****************')") neighbor

                  write(*,"('  Charge current:')")
                  write(*,"('     Ich = (',e16.9,') + i(',e16.9,')')") real(Ich(neighbor,i)),aimag(Ich(neighbor,i))
                  write(*,"(' abs(Ich) = ',e16.9)") abs(Ich(neighbor,i))
                  write(*,"('atan(Ich) = ',e16.9)") atan2(aimag(Ich(neighbor,i)),real(Ich(neighbor,i)))

                  write(*,"('  Spin currents:')")
                  write(*,"('     Isx  = (',e16.9,') + i(',e16.9,')')") real(Isx(neighbor,i)),aimag(Isx(neighbor,i))
                  write(*,"(' abs(Isx) = ',e16.9)") abs(Isx(neighbor,i))
                  write(*,"('atan(Isx) = ',e16.9)") atan2(aimag(Isx(neighbor,i)),real(Isx(neighbor,i)))

                  write(*,"('     Isy  = (',e16.9,') + i(',e16.9,')')") real(Isy(neighbor,i)),aimag(Isy(neighbor,i))
                  write(*,"(' abs(Isy) = ',e16.9)") abs(Isy(neighbor,i))
                  write(*,"('atan(Isy) = ',e16.9)") atan2(aimag(Isy(neighbor,i)),real(Isy(neighbor,i)))

                  write(*,"('     Isz  = (',e16.9,') + i(',e16.9,')')") real(Isz(neighbor,i)),aimag(Isz(neighbor,i))
                  write(*,"(' abs(Isz) = ',e16.9)") abs(Isz(neighbor,i))
                  write(*,"('atan(Isz) = ',e16.9)") atan2(aimag(Isz(neighbor,i)),real(Isz(neighbor,i)))

                  write(*,"('  Orbital Angular Momentum currents:')")
                  write(*,"('     Ilx  = (',e16.9,') + i(',e16.9,')')") real(Ilx(neighbor,i)),aimag(Ilx(neighbor,i))
                  write(*,"(' abs(Ilx) = ',e16.9)") abs(Ilx(neighbor,i))
                  write(*,"('atan(Ilx) = ',e16.9)") atan2(aimag(Ilx(neighbor,i)),real(Ilx(neighbor,i)))

                  write(*,"('     Ily  = (',e16.9,') + i(',e16.9,')')") real(Ily(neighbor,i)),aimag(Ily(neighbor,i))
                  write(*,"(' abs(Ily) = ',e16.9)") abs(Ily(neighbor,i))
                  write(*,"('atan(Ily) = ',e16.9)") atan2(aimag(Ily(neighbor,i)),real(Ily(neighbor,i)))

                  write(*,"('     Ilz  = (',e16.9,') + i(',e16.9,')')") real(Ilz(neighbor,i)),aimag(Ilz(neighbor,i))
                  write(*,"(' abs(Ilz) = ',e16.9)") abs(Ilz(neighbor,i))
                  write(*,"('atan(Ilz) = ',e16.9)") atan2(aimag(Ilz(neighbor,i)),real(Ilz(neighbor,i)))

                  ! CHARGE CURRENT
                  ! Writing charge current
                  iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7
                  write(unit=iw+1,fmt="(7(e16.9,2x))") e,abs(Ich(neighbor,i)),real(Ich(neighbor,i)),aimag(Ich(neighbor,i)),atan2(aimag(Ich(neighbor,i)),real(Ich(neighbor,i))),real(Ich(neighbor,i))/abs(Ich(neighbor,i)),aimag(Ich(neighbor,i))/abs(Ich(neighbor,i))

                  ! SPIN CURRENTS
                  ! Writing x-component spin current
                  write(unit=iw+2,fmt="(7(e16.9,2x))") e,abs(Isx(neighbor,i)),real(Isx(neighbor,i)),aimag(Isx(neighbor,i)),atan2(aimag(Isx(neighbor,i)),real(Isx(neighbor,i))),real(Isx(neighbor,i))/abs(Isx(neighbor,i)),aimag(Isx(neighbor,i))/abs(Isx(neighbor,i))
                  ! Writing y-component spin current
                  write(unit=iw+3,fmt="(7(e16.9,2x))") e,abs(Isy(neighbor,i)),real(Isy(neighbor,i)),aimag(Isy(neighbor,i)),atan2(aimag(Isy(neighbor,i)),real(Isy(neighbor,i))),real(Isy(neighbor,i))/abs(Isy(neighbor,i)),aimag(Isy(neighbor,i))/abs(Isy(neighbor,i))
                  ! Writing z-component spin current
                  write(unit=iw+4,fmt="(7(e16.9,2x))") e,abs(Isz(neighbor,i)),real(Isz(neighbor,i)),aimag(Isz(neighbor,i)),atan2(aimag(Isz(neighbor,i)),real(Isz(neighbor,i))),real(Isz(neighbor,i))/abs(Isz(neighbor,i)),aimag(Isz(neighbor,i))/abs(Isz(neighbor,i))

                  ! ORBITAL ANGULAR MOMENTUM CURRENTS
                  ! Writing x-component orbital angular momentum current
                  write(unit=iw+5,fmt="(7(e16.9,2x))") e,abs(Ilx(neighbor,i)),real(Ilx(neighbor,i)),aimag(Ilx(neighbor,i)),atan2(aimag(Ilx(neighbor,i)),real(Ilx(neighbor,i))),real(Ilx(neighbor,i))/abs(Ilx(neighbor,i)),aimag(Ilx(neighbor,i))/abs(Ilx(neighbor,i))
                  ! Writing y-component orbital angular momentum current
                  write(unit=iw+6,fmt="(7(e16.9,2x))") e,abs(Ily(neighbor,i)),real(Ily(neighbor,i)),aimag(Ily(neighbor,i)),atan2(aimag(Ily(neighbor,i)),real(Ily(neighbor,i))),real(Ily(neighbor,i))/abs(Ily(neighbor,i)),aimag(Ily(neighbor,i))/abs(Ily(neighbor,i))
                  ! Writing z-component orbital angular momentum current
                  write(unit=iw+7,fmt="(7(e16.9,2x))") e,abs(Ilz(neighbor,i)),real(Ilz(neighbor,i)),aimag(Ilz(neighbor,i)),atan2(aimag(Ilz(neighbor,i)),real(Ilz(neighbor,i))),real(Ilz(neighbor,i))/abs(Ilz(neighbor,i)),aimag(Ilz(neighbor,i))/abs(Ilz(neighbor,i))

                  ! Writing renormalized currents
                  if(renorm) then
                    ! CHARGE CURRENT
                    ! Writing renormalized charge current
                    write(unit=iw+1001,fmt="(7(e16.9,2x))") e,abs(rIch(neighbor,i)),real(rIch(neighbor,i)),aimag(rIch(neighbor,i)),atan2(aimag(rIch(neighbor,i)),real(rIch(neighbor,i))),real(rIch(neighbor,i))/abs(rIch(neighbor,i)),aimag(rIch(neighbor,i))/abs(rIch(neighbor,i))

                    ! SPIN CURRENTS
                    ! Writing renormalized x-component spin current
                    write(unit=iw+1002,fmt="(7(e16.9,2x))") e,abs(rIsx(neighbor,i)),real(rIsx(neighbor,i)),aimag(rIsx(neighbor,i)),atan2(aimag(rIsx(neighbor,i)),real(rIsx(neighbor,i))),real(rIsx(neighbor,i))/abs(rIsx(neighbor,i)),aimag(rIsx(neighbor,i))/abs(rIsx(neighbor,i))
                    ! Writing renormalized y-component spin current
                    write(unit=iw+1003,fmt="(7(e16.9,2x))") e,abs(rIsy(neighbor,i)),real(rIsy(neighbor,i)),aimag(rIsy(neighbor,i)),atan2(aimag(rIsy(neighbor,i)),real(rIsy(neighbor,i))),real(rIsy(neighbor,i))/abs(rIsy(neighbor,i)),aimag(rIsy(neighbor,i))/abs(rIsy(neighbor,i))
                    ! Writing renormalized z-component spin current
                    write(unit=iw+1004,fmt="(7(e16.9,2x))") e,abs(rIsz(neighbor,i)),real(rIsz(neighbor,i)),aimag(rIsz(neighbor,i)),atan2(aimag(rIsz(neighbor,i)),real(rIsz(neighbor,i))),real(rIsz(neighbor,i))/abs(rIsz(neighbor,i)),aimag(rIsz(neighbor,i))/abs(rIsz(neighbor,i))

                    ! ORBITAL ANGULAR MOMENTUM CURRENTS
                    ! Writing x-component orbital angular momentum current
                    write(unit=iw+1005,fmt="(7(e16.9,2x))") e,abs(rIlx(neighbor,i)),real(rIlx(neighbor,i)),aimag(rIlx(neighbor,i)),atan2(aimag(rIlx(neighbor,i)),real(rIlx(neighbor,i))),real(rIlx(neighbor,i))/abs(rIlx(neighbor,i)),aimag(rIlx(neighbor,i))/abs(rIlx(neighbor,i))
                    ! Writing y-component orbital angular momentum current
                    write(unit=iw+1006,fmt="(7(e16.9,2x))") e,abs(rIly(neighbor,i)),real(rIly(neighbor,i)),aimag(rIly(neighbor,i)),atan2(aimag(rIly(neighbor,i)),real(rIly(neighbor,i))),real(rIly(neighbor,i))/abs(rIly(neighbor,i)),aimag(rIly(neighbor,i))/abs(rIly(neighbor,i))
                    ! Writing z-component orbital angular momentum current
                    write(unit=iw+1007,fmt="(7(e16.9,2x))") e,abs(rIlz(neighbor,i)),real(rIlz(neighbor,i)),aimag(rIlz(neighbor,i)),atan2(aimag(rIlz(neighbor,i)),real(rIlz(neighbor,i))),real(rIlz(neighbor,i))/abs(rIlz(neighbor,i)),aimag(rIlz(neighbor,i))/abs(rIlz(neighbor,i))
                  end if
                end do neighbor_loop_write_all
              end do plane_loop_write_all
            end do MPI_points_all
            ! Closing files
            call openclosechifiles(2)
            call openclosesdfiles(2)
            call openclosescfiles(2)

            call date_and_time(date, time, zone, values)
            write(*,"('[main] Time after step ',i0,': ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") count,values(3),values(2),values(1),values(5),values(6),values(7)
            elapsed_time = MPI_Wtime() - start_program
            write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

            ! Emergency stop
            open(unit=911, file="stop", status='old', iostat=iw)
            if(iw.eq.0) then
              write(*,"(a,i0,a)") "[main] Emergency 'stop' file found! Stopping after step ",count," ..."
              call system ('rm stop')
              write(*,"(a)") "[main] ('stop' file deleted!)"
              call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
            end if
          else
            call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,4000,MPIComm_Col,ierr)
            call MPI_Send(schi,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,4010,MPIComm_Col,ierr)
            call MPI_Send(schihf,Npl*Npl*4,MPI_DOUBLE_COMPLEX,0,4020,MPIComm_Col,ierr)
            call MPI_Send(sdx,Npl,MPI_DOUBLE_COMPLEX,0,4030,MPIComm_Col,ierr)
            call MPI_Send(sdy,Npl,MPI_DOUBLE_COMPLEX,0,4040,MPIComm_Col,ierr)
            call MPI_Send(sdz,Npl,MPI_DOUBLE_COMPLEX,0,4050,MPIComm_Col,ierr)
            call MPI_Send(chd,Npl,MPI_DOUBLE_COMPLEX,0,4060,MPIComm_Col,ierr)
            call MPI_Send(ldx,Npl,MPI_DOUBLE_COMPLEX,0,4070,MPIComm_Col,ierr)
            call MPI_Send(ldy,Npl,MPI_DOUBLE_COMPLEX,0,4080,MPIComm_Col,ierr)
            call MPI_Send(ldz,Npl,MPI_DOUBLE_COMPLEX,0,4090,MPIComm_Col,ierr)
            call MPI_Send(Ich,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4100,MPIComm_Col,ierr)
            call MPI_Send(Isx,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4110,MPIComm_Col,ierr)
            call MPI_Send(Isy,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4120,MPIComm_Col,ierr)
            call MPI_Send(Isz,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4130,MPIComm_Col,ierr)
            call MPI_Send(Ilx,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4140,MPIComm_Col,ierr)
            call MPI_Send(Ily,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4150,MPIComm_Col,ierr)
            call MPI_Send(Ilz,n0sc*Npl,MPI_DOUBLE_COMPLEX,0,4160,MPIComm_Col,ierr)
          end if
        end if
      end do all_energy_loop

      deallocate(chiorb_hf,prefactor,tchiorbiikl,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)
      if(myrank_row.eq.0) then
        if(myrank_col.eq.0) then
          if(nmaglayers.gt.1) then
            deallocate(rwork,eval,evecl,evecr,work)
#ifdef _JUQUEEN
            deallocate(dscale,rconde,rcondv)
#endif
          end if
        end if
        deallocate(templd,chiorb)
        if(renorm) deallocate(rIch,rIsx,rIsy,rIsz,rIlx,rIly,rIlz,rsdx,rsdy,rsdz,rchd,rldx,rldy,rldz)
        deallocate(Ich,Isx,Isy,Isz,Ilx,Ily,Ilz)
        deallocate(schi,schihf,sdx,sdy,sdz,chd,ldx,ldy,ldz)
      end if
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

      ! Opening files for position dependence
      if((myrank.eq.0).and.(Npl.eq.Nplini)) then
        ! Exchange interactions
        do j=1,nmaglayers ; do i=1,nmaglayers
          iw = 199+(j-1)*nmaglayers*2+(i-1)*2
          if(i.eq.j) then
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Jij/Jii_',i0,'_magaxis=',A,'_parts=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,i,magaxis,parts,ncp,eta,Utype,hwx,hwy,hwz
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   Npl         ,  Jii_xx           ,   Jii_yy  ')")
            iw = iw + 1
            ! TODO : Check how to write the anisotropy term here
          else
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Jij/J_',i0,'_',i0,'_magaxis=',A,'_parts=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,i,j,magaxis,parts,ncp,eta,Utype,hwx,hwy,hwz
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   Npl         ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
            iw = iw + 1
            write(varm,"('./results/SOC=',L1,'/Jij/Dz_',i0,'_',i0,'_magaxis=',A,'_parts=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,i,j,magaxis,parts,ncp,eta,Utype,hwx,hwy,hwz
            open (unit=iw, file=varm,status='unknown')
            write(unit=iw, fmt="('#   Npl         , Dz = (Jxy - Jyx)/2       ')")
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
        write(*,"('  ******************** Full tensor Jij:  ********************')")
        do i=1,nmaglayers ; do j=1,nmaglayers
        ! Writing original full tensor Jij
          if(i.eq.j) then
            write(*,"(' |--------------- i = ',i0,'   j = ',i0,': anisotropies ---------------|')") mmlayermag(i),mmlayermag(j)
          else
            write(*,"(' |----------- i = ',i0,'   j = ',i0,': exchange couplings -------------|')") mmlayermag(i),mmlayermag(j)
          end if
          write(*,"('             x                  y                  z')")
          write(*,"('  x  (',e16.9,') (',e16.9,') (',e16.9,')')") Jij(i,j,1,1),Jij(i,j,1,2),Jij(i,j,1,3)
          write(*,"('  y  (',e16.9,') (',e16.9,') (',e16.9,')')") Jij(i,j,2,1),Jij(i,j,2,2),Jij(i,j,2,3)
          write(*,"('  z  (',e16.9,') (',e16.9,') (',e16.9,')')") Jij(i,j,3,1),Jij(i,j,3,2),Jij(i,j,3,3)
        end do ; end do
        write(*,"('  *** Symmetryc and antisymmetric exchange interactions:  ***')")
        do i=1,nmaglayers ; do j=1,nmaglayers
          if(i.eq.j) cycle
          write(*,"(' |--------------------- i = ',i0,'   j = ',i0,' -----------------------|')") mmlayermag(i),mmlayermag(j)
        ! Writing Heisenberg exchange interactions
          write(*,"('     Isotropic:     J     = ',e16.9)") trJij(i,j)
          write(*,"('   Anisotropic:     Js_xx = ',e16.9)") Jijs(i,j,1,1)
          write(*,"('                    Js_yy = ',e16.9)") Jijs(i,j,2,2)
          write(*,"('  DMI: Dz = (Jxy - Jyx)/2 = ',e16.9)") Jija(i,j,1,2)
          write(*,"(' --- z components of Jij (not physically correct) ---')")
          write(*,"('  Anisotropic:  Js_zz = ',e16.9)") Jijs(i,j,3,3)
          write(*,"('  DMI: Dy = (Jzx - Jxz)/2 = ',e16.9)") -Jija(i,j,1,3)
          write(*,"('  DMI: Dx = (Jyz - Jzy)/2 = ',e16.9)") Jija(i,j,2,3)
        end do ; end do

        ! Writing into files
        ! Exchange interactions
        exchange_writing_loop: do j=1,nmaglayers ; do i=1,nmaglayers
          iw = 199+(j-1)*nmaglayers*2+(i-1)*2
          if(i.eq.j) then
            iw = iw + 1
            write(unit=iw,fmt="(4x,i3,13x,2(e16.9,2x))") Npl,Jij(i,j,1,1),Jij(i,j,2,2)
            iw = iw + 1
          else
            iw = iw + 1
            write(unit=iw,fmt="(4x,i3,13x,3(e16.9,2x))") Npl,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
            iw = iw + 1
            write(unit=iw,fmt="(4x,i3,13x,e16.9,2x)") Npl,Jija(i,j,1,2)
          end if
        end do ; end do exchange_writing_loop

        ! Closing files
        if(Npl.eq.Nplfinal) then
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
      if(myrank.eq.0) then
        allocate(sdl(Npl), STAT = AllocateStatus)
        if (AllocateStatus.ne.0) then
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if

        write(*,"('CALCULATING PROBABILITY OF SPIN FLIP AS A FUNCTION OF POSITION')")
        do i=1,Npl
          write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/sdl_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'.dat')") SOC,Npl,i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis
          open (unit=555+i, file=varm,status='unknown')
        end do

        e = Ef
        inn = 1
        do pos=0,100
          call spindifflength(e,pos,inn,sdl)
          do i=1,Npl
            sdl2 = abs(sdl(i))**2
            ! Writing probability of spin-flip as a function os position
            write(unit=555+i,fmt="(i0,2x,e16.9)") pos,sdl2
          end do
        end do
        deallocate(sdl)
      end if
!-----------------------------------------------------------------------


    end select main_program
!-------------- Deallocating variables that depend on Npl --------------
    deallocate(sigmaimunu2i,sigmaijmunu2i,eps1,eps1_solu)
    deallocate(mag,mag_0,hdel,splus,splus_0,usplus,sminus,usminus,npart,jac,wa)
    if(index(runoptions,"GSL").gt.0) deallocate(lxm,lym,lzm,lxpm,lypm,lzpm)
    deallocate(mmlayer,layertype,U,mmlayermag,lambda,npart0)
    deallocate(identt,Umatorb)
    select case (plnn)
    case(1)
      deallocate(t00,t01)
    case(2)
      deallocate(t00,t01,t02)
    end select
  end do number_of_planes

  deallocate(r0,c0,r1,c1,r2,c2)
  deallocate(kbz,wkbz,kbz2d)
!----------------------- Finalizing the program ------------------------
  if(myrank.eq.0) then
    call date_and_time(date, time, zone, values)
    write(*,"('[main] Finished on: ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") values(3),values(2),values(1),values(5),values(6),values(7)
    elapsed_time = MPI_Wtime() - start_program
    write(*,"('[main] Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes / ',f7.4,' hours')") elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0
  end if
  call MPI_Finalize(ierr)
  if (ierr.ne.0) then
    write(*,"('[main] ierr = ',i0,'. Something went wrong in the parallelization!')") ierr
  end if
!=======================================================================
  stop
end program SHE