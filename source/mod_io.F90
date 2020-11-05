module mod_io
  ! Subroutines to read input file and variables containing these parameters
  implicit none
  logical :: log_unit = .false.
  character(len=12), parameter :: logfile = "parameter.in"

contains

  subroutine log_message(procedure, message)
    use mod_mpi_pars,   only: myrank
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
       if(log_unit) then
          write(output%unit, "('[',a,'] ',a,'')") procedure, trim(message)
       else
          write(*, "('[',a,'] ',a,'')") procedure, trim(message)
       end if
    end if
  end subroutine log_message

  subroutine log_error(procedure, message)
    use mod_mpi_pars, only: myrank,MPI_Abort,MPI_COMM_WORLD,errorcode,ierr
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
      if(log_unit) write(output%unit, "('[Error] [',a,'] ',a)") procedure, trim(message)
    end if
    write(*, "('[Error] [',a,'] ',a)") procedure, trim(message)
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    stop
  end subroutine log_error

  subroutine log_warning(procedure, message)
    use mod_mpi_pars,   only: myrank
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
       if(log_unit) then
          write(output%unit, "('[Warning] [',a,'] ',a)") procedure, trim(message)
       else
          write(*, "('[Warning] [',a,'] ',a)") procedure, trim(message)
       end if
    end if
  end subroutine log_warning

  subroutine get_parameters(filename,s)
    use mod_kind,              only: dp,int64
    use mod_input,             only: get_parameter,read_file,enable_input_logging,disable_input_logging
    use mod_parameters,        only: output,laddresults,lverbose,ldebug,lkpoints,&
                                     lpositions,lcreatefiles,lnolb,lhfresponses,&
                                     lnodiag,lsha,lcreatefolders,lwriteonscreen,runoptions,lsimplemix,&
                                     lcheckjac,llgtv,lsortfiles,leigenstates,lprintfieldonly,&
                                     itype,ry2ev,ltesla,eta,etap,dmax,emin,emax,&
                                     skip_steps,nEner,nEner1,nQvec,nQvec1,qbasis,renorm,renormnb,bands,band_cnt,&
                                     offset,dfttype,parField,parFreq,kptotal_in,kp_in,&
                                     nOrb,nOrb2,tbmode,fermi_layer
    use mod_superconductivity, only: lsuperCond,superCond
    use mod_self_consistency,  only: lontheflysc,lnojac,lforceoccup,lrotatemag,skipsc,scfile,magbasis,mag_tol
    use mod_system,            only: System_type,n0sc1,n0sc2
    use mod_SOC,               only: SOC,socscale,llinearsoc,llineargfsoc
    use mod_magnet,            only: lfield,tesla,hwa_i,hwa_f,hwa_npts,hwa_npt1,hwt_i,hwt_f,&
                                     hwt_npts,hwt_npt1,hwp_i,hwp_f,hwp_npts,hwp_npt1,hwx,hwy,&
                                     hwz,hwscale,hwtrotate,hwprotate,skip_steps_hw
    use ElectricField,         only: ElectricFieldMode,ElectricFieldVector,EFp,EFt
    use EnergyIntegration,     only: parts,parts3,pn1,pn2,pnt,n1gl,n3gl
    use mod_tools,             only: itos,rtos,vec_norm
    use adaptiveMesh,          only: minimumBZmesh
    use mod_fermi_surface,     only: lfs_loop,fs_energy_npts,fs_energy_npt1,fs_energy_i,fs_energy_f
    use mod_mpi_pars,          only: myrank,ierr
    use mod_imRK4_parameters,  only: integration_time,sc_tol,step,hE_0,hw1_m,hw_e,hw_m,tau_e,&
                                     polarization_e,polarization_m,polarization_vec_e,polarization_vec_m,&
                                     npulse_e,npulse_m,tau_m,delay_e,delay_m,lelectric,safe_factor,&
                                     lmagnetic,lpulse_e,lpulse_m,abs_tol,rel_tol
    use mod_expectation,       only: expectation_values,expectation_values_eigenstates,calc_GS_L_and_E,calc_GS_L_and_E_eigenstates
    use mod_greenfunction,     only: calc_green,greenlineargfsoc
    implicit none
    character(len=*),  intent(in)    :: filename
    type(System_type), intent(inout) :: s
    character(len=20), allocatable   :: s_vector(:)
    real(dp),          allocatable   :: vector(:)
    integer(int64),    allocatable   :: i_vector(:)
    integer :: i,cnt
    character(len=20) :: tmp_string

    external :: MPI_Finalize

    if(.not. read_file(filename)) &
      call log_error("get_parameters", "File " // trim(filename) // " not found!")

    if(myrank == 0) then
      if(.not. enable_input_logging(logfile)) &
        call log_warning("get_parameters", "couldn't enable logging.")
    end if
    if(.not. get_parameter("output", output%file)) &
      call log_error("get_parameters", "Output filename not given!")

    if(myrank==0) then
      !! Create folder for output file, if necessary
      call execute_command_line('mkdir -p $(dirname "'// trim(output%file) // '")')
      !! Open output file
      open(unit=output%unit, file=trim(output%file), status='replace')
    end if

    log_unit = .true.

! Print the Git version (VERSION is defined via CMake macro and defined with compiler flag -DVERSION='')
#if defined(VERSION)
    if(myrank==0) write(output%unit,"('Git commit: ',a)") VERSION
#else
    if(myrank==0) write(output%unit,"('Git commit: unknown')")
#endif

    if(myrank==0) &
      write(output%unit,"('[get_parameters] Reading parameters from ""',a,'"" file...')") trim(filename)
    !===============================================================================================
    !============= System configuration (Lattice + Reciprocal lattice) =============================
    !===============================================================================================
    if(.not. get_parameter("nn_stages", s%nStages,2)) &
      call log_warning("get_parameters","'nn_stages' missing. Using default value: 2")
    if(.not. get_parameter("relTol", s%relTol,0.05_dp)) &
      call log_warning("get_parameters","'relTol' missing. Using default value: 0.05")
    if(.not. get_parameter("sysdim", s%isysdim, 3)) &
      call log_warning("get_parameters", "'sysdim' missing. Using default value: 3")
    if(.not. get_parameter("nkpt", i_vector,cnt)) &
      call log_error("get_parameters","'nkpt' missing.")
    if(cnt == 1) then
      kptotal_in = int( i_vector(1), kind(kptotal_in) )
      select case(s%isysdim)
      case(3)
        kp_in(:)   = ceiling((dble(kptotal_in))**(0.333333333333333_dp), kind(kp_in(1)) )
        kptotal_in = int( kp_in(1) * kp_in(2) * kp_in(3), kind(kptotal_in) )
      case(2)
        kp_in(1:2) = ceiling((dble(kptotal_in))**(0.5_dp), kind(kp_in(1)) )
        kp_in(3)   = 1
        kptotal_in = int( kp_in(1) * kp_in(2), kind(kptotal_in) )
      case default
        kp_in(1)   = ceiling((dble(kptotal_in)), kind(kp_in(1)) )
        kp_in(2:3) = 1
        kptotal_in = int( kp_in(1), kind(kptotal_in) )
      end select

    else if(cnt == 3) then
      kp_in(1) = int( i_vector(1), kind(kp_in(1)) )
      kp_in(2) = int( i_vector(2), kind(kp_in(2)) )
      kp_in(3) = int( i_vector(3), kind(kp_in(3)) )
      kptotal_in = int( kp_in(1) * kp_in(2) * kp_in(3), kind(kptotal_in) )
    else
      call log_error("get_parameter", "'nkpt' has wrong size (expected 1 or 3).")
    end if
    if(.not. get_parameter("minimumBZmesh", minimumBZmesh, 1000)) &
      call log_warning("get_parameters", "'minimumBZmesh' missing. Using default value: 1000")
    if(.not. get_parameter("eta", eta)) &
      call log_error("get_parameters","'eta' missing.")
    if(.not. get_parameter("etap", etap, eta)) &
      call log_warning("get_parameters", "'etap' not found. Using default value: eta = " // trim(rtos(eta,"(es8.1)")) )
    !===============================================================================================
    !===============================================================================================
    !------------------------------------- Type of Calculation -------------------------------------
    if(.not. get_parameter("itype", itype)) &
      call log_error("get_parameters","'itype' missing.")
    if(.not. get_parameter("Options", s_vector, cnt)) &
      call log_warning("get_parameters","'Options' missing.")
    runoptions = ""
    do i = 1, cnt
      select case (s_vector(i))
      case ("ry2ev")
        ry2ev = 13.6_dp
      case ("tesla")
        tesla = 5.7883817555e-5_dp/13.6_dp ! Ry/T
        ltesla = .true.
      case ("verbose")
        lverbose = .true.
      case ("debug")
        ldebug = .true.
      case ("addresults")
        laddresults = .true.
      case ("createfiles")
        lcreatefiles = .true.
      case("createfolders")
        lcreatefolders = .true.
      case ("kpoints")
        lkpoints = .true.
      case ("positions")
        lpositions = .true.
      case ("lineargfsoc")
        llineargfsoc = .true.
        calc_green => greenlineargfsoc
      case ("linearsoc")
        llinearsoc = .true.
        calc_green => greenlineargfsoc
      case ("nojac")
        lnojac = .true.
      case ("hfresponses")
        lhfresponses  = .true.
      case ("ontheflysc")
        lontheflysc = .true.
      case ("rotatemag")
        lrotatemag = .true.
      case ("nolb")
        lnolb = .true.
      case ("nodiag")
        lnodiag = .true.
      case ("sha")
        lsha = .true.
      case ("writeonscreen")
        lwriteonscreen = .true.
      case ("sortfiles")
        lsortfiles = .true.
      case ("lgtv")
        llgtv = .true.
      case ("checkjac")
        lcheckjac = .true.
      case ("forceoccupation")
        lforceoccup = .true.
      case ("simplemix")
        lsimplemix = .true.
      case ("eigenstates")
        leigenstates = .true.
        expectation_values => expectation_values_eigenstates
        calc_GS_L_and_E => calc_GS_L_and_E_eigenstates
        lnojac = .true.
        call log_warning("get_parameters","eigenstates is used, jacobian deactivated (not implemented yet)")
      case ("printfieldonly")
        lprintfieldonly = .true.
      case("!")
        exit
      case default
        call log_warning("get_parameters","Runoption " // trim(s_vector(i)) // " not found!")
        cycle
      end select
      runoptions  = trim(runoptions) // " " // trim(s_vector(i))
    end do
    deallocate(s_vector)
    !-------------------------------------- In-Plane Currents --------------------------------------
    ! Delete?
    if(.not. get_parameter("n0sc1", n0sc1)) &
      call log_warning("get_parameters","'n0sc1' missing.")
    if(.not. get_parameter("n0sc2", n0sc2)) &
      call log_warning("get_parameters","'n0sc2' missing.")

    !------------------------------------- Fermi Surface -------------------------------------
    if(itype==5) then
      if(.not. get_parameter("fs_energy", vector, cnt)) then
        call log_warning("get_parameters","'fs_energy' missing. Using Fermi level only.")
      else
        if(cnt < 1) call log_error("get_parameters","'fs_energy' doesn't contain any parameter.")
        fs_energy_i = vector(1)
        if(cnt >= 2) fs_energy_f = vector(2)
        if(cnt >= 3) fs_energy_npts = int(vector(3))
        deallocate(vector)
        fs_energy_npt1 = fs_energy_npts + 1
        lfs_loop = .true.
      end if
    end if

    !------------------------------------- Spin-Orbit-Coupling -------------------------------------
    if(.not. get_parameter("SOC", SOC,.true.)) &
      call log_warning("get_parameters","'SOC' missing. Using default value: .true.")
    if(SOC) then
      if(.not. get_parameter("socscale", socscale, 1._dp)) &
        call log_warning("get_parameters","'socscale' missing. Using default value: 1._dp")
      if(llinearsoc) then
        output%SOCchar = "L"
      else
        output%SOCchar = "T"
      end if
      if(abs(socscale-1._dp)>1.e-6_dp) write(output%SOC,"('_socscale=',f5.2)") socscale
      if((llineargfsoc).or.(llinearsoc)) output%SOC = trim(output%SOC) // "_linearsoc"
    else
      output%SOCchar = "F"
    end if
    !---------------------------------------- Magnetization -----------------------------------------
    if(.not. get_parameter("magtol", mag_tol, 1.e-12_dp)) &
      call log_warning("get_parameters", "'magtol' not found. Using default value: 1.e-12_dp")
    if(.not. get_parameter("magbasis", magbasis)) &
      call log_warning("get_parameters","'magbasis' missing. Using default values for initial magnetization")
    !--------------------------------------- Electric Field ----------------------------------------
    if(.not. get_parameter("ebasis", tmp_string, "spherical")) &
      call log_warning("get_parameters","'ebasis' missing. Using default value: ""spherical""")
    select case (tmp_string)
    case("cartesian")
      if(.not. get_parameter("dirEfield", vector, cnt)) &
        call log_error("get_parameters","'dirEfield' missing.")
      if(cnt /= 3) call log_error("get_parameters","'dirEfield' has wrong size (size 3 required).")
      ElectricFieldMode = -1 ! TODO: Set it to a value if not determined otherwise?
      ElectricFieldVector(1:3) = vector(1:3)
      deallocate(vector)
    case("neighbor")
      if(.not. get_parameter("dirEfield", ElectricFieldMode)) &
        call log_error("get_parameters","'dirEfield' missing.")
    case("bravais")
      if(.not. get_parameter("dirEfield", i_vector, cnt)) &
        call log_error("get_parameters","'dirEfield' missing.")
      if(cnt /= 2) call log_error("get_parameters","'dirEfield' has wrong size (size 2 required).")
      ElectricFieldMode = -2 ! TODO: Add options to evaluate these values.
      ElectricFieldVector(1:2) = i_vector(1:2)
      deallocate(i_vector)
    case("spherical")
      if(.not. get_parameter("dirEfield", vector, cnt)) &
        call log_error("get_parameters", "'dirEfield' missing.")
      if(cnt /= 2) call log_error("get_parameters", "'dirEfield' has wrong size (size 2 required).")
      ElectricFieldMode = -3
      EFt = vector(1)
      EFp = vector(2)
      deallocate(vector)
    end select
    !------------------------------------- Static Magnetic Field -----------------------------------
    if(.not. get_parameter("FIELD", lfield, .false.)) &
      call log_warning("get_parameters","'lfield' missing. Using default value: .false.")
    if(lfield) then
      if(.not. get_parameter("hwa", vector, cnt)) &
        call log_error("get_parameters","'hwa' missing.")
      if(cnt < 1) call log_error("get_parameters","'hwa' doesn't contain any parameter.")
      hwa_i = vector(1)
      if(cnt >= 2) hwa_f = vector(2)
      if(cnt >= 3) hwa_npts = int(vector(3))
      deallocate(vector)
      hwa_npt1 = hwa_npts + 1

      if(.not. get_parameter("hwt", vector, cnt)) &
        call log_error("get_parameters","'hwt' missing.")
      if(cnt < 1) call log_error("get_parameters","'hwt' doesn't contain any parameter.")
      hwt_i = vector(1)
      if(cnt >= 2) hwt_f = vector(2)
      if(cnt >= 3) hwt_npts = int(vector(3))
      deallocate(vector)
      hwt_npt1 = hwt_npts + 1

      if(.not. get_parameter("hwp", vector, cnt)) &
        call log_error("get_parameters","'hwp' missing.")
      if(cnt < 1) call log_error("get_parameters","'hwp' doesn't contain any parameter.")
      hwp_i = vector(1)
      if(cnt >= 2) hwp_f = vector(2)
      if(cnt >= 3) hwp_npts = int(vector(3))
      deallocate(vector)
      hwp_npt1 = hwp_npts + 1

      if(abs(hwa_i) < 1.e-9_dp) then
        if(.not. get_parameter("hwx", hwx)) &
          call log_error("get_parameters","'hwx' missing.")
        if(.not. get_parameter("hwy", hwy)) &
          call log_error("get_parameters","'hwy' missing.")
        if(.not. get_parameter("hwz", hwz)) &
          call log_error("get_parameters","'hwz' missing.")
      end if
    end if
    if(.not. get_parameter("skip_steps_hw", skip_steps_hw, 0)) &
      call log_warning("get_parameters","'skip_steps_hw' missing. Using default value: 0")

    if(get_parameter("hwscale", vector, cnt)) then
       if(cnt < dmax)  hwscale(1:cnt)  = vector(1:cnt)
       if(cnt >= dmax) hwscale(1:dmax) = vector(1:dmax)
    end if
    if(allocated(vector)) deallocate(vector)

    if(get_parameter("hwtrotate", vector, cnt)) then
       if(cnt < dmax)  hwtrotate(1:cnt)  = vector(1:cnt)
       if(cnt >= dmax) hwtrotate(1:dmax) = vector(1:dmax)
    end if
    if(allocated(vector)) deallocate(vector)

    if(get_parameter("hwprotate", vector, cnt)) then
       if(cnt < dmax)  hwprotate(1:cnt)  = vector(1:cnt)
       if(cnt >= dmax) hwprotate(1:dmax) = vector(1:dmax)
    end if
    if(allocated(vector)) deallocate(vector)
    !------------------------------ Superconductivity Variables ------------------------------------
    if(.not. get_parameter("superCond", lsuperCond, .false.)) &
      call log_warning("get_parameters","'superCond' missing. Not using superconductivity.")
    superCond = merge(2,1,lsuperCond)
    !------------------------------------ Integration Variables ------------------------------------
    if(.not. get_parameter("parts", parts)) &
      call log_error("get_parameters","'parts' missing.")
    if(.not. get_parameter("parts3", parts3)) &
      call log_error("get_parameters","'parts3' missing.")
    write(output%Energy, "('_parts=',i0,'_parts3=',i0)") parts,parts3
    if(.not. get_parameter("n1gl", n1gl)) &
      call log_error("get_parameters","'n1gl' missing.")
    if(.not. get_parameter("n3gl", n3gl)) &
      call log_error("get_parameters","'n3gl' missing.")
    !------------------------------------ Loop Variables ------------------------------------
    ! Energy (frequency):
    if(.not. get_parameter("emin", emin)) &
      call log_error("get_parameters","'emin' missing.")
    if(.not. get_parameter("emax", emax)) &
      call log_error("get_parameters","'emax' missing.")
    if(.not. get_parameter("skip_steps", skip_steps, 0)) &
      call log_warning("get_parameters","'skip_steps' missing. Using default value: 0")
    if(.not. get_parameter("nEner", nEner)) &
      call log_error("get_parameters","'nEner' missing.")
    nEner1 = nEner + 1

    ! Wave vectors:
    if(.not. get_parameter("nQvec", nQvec, 0)) &
      call log_warning("get_parameters","'nQvec' not found. No wave vector loop will be done.")
    nQvec1 = nQvec + 1
    if(.not. get_parameter("renorm", renorm, .false.)) & ! Delete?
      call log_warning("get_parameters","'renorm' missing. Using default value .false.")
    if(renorm) then ! Delete?
      if(.not. get_parameter("renormnb", renormnb)) &
        call log_error("get_parameters","'renormnb' missing.")
    end if
    !----------------- Wave vector loop (band structure and susceptibility)  ----------------
    if((itype >= 7).and.(itype <= 9)) then
      ! Path or point to calculate the susceptibility
      if(.not. get_parameter("band", bands, band_cnt)) then
        call log_warning("get_parameters", "'band' missing. Using Gamma point only.")
        allocate(bands(1))
        bands = "G"
        band_cnt = 1
        nQvec1 = 1
      end if
      if((band_cnt == 1).and.(nQvec1 > 1)) then
        call log_warning("get_parameters", "Only one wave vector given. Using nQvec=1.")
        nQvec1 = 1
      end if
    endif
    if(itype == 4) then
      ! Path to calculate band structure (can't be single point in this case)
      if(.not. get_parameter("band", bands, band_cnt)) &
        call log_error("get_parameters", "'band' missing.")
      if((band_cnt < 2).or.(nQvec < 2)) call log_error("get_parameters", "Need at least two points for Band Structure")
    endif
    if(.not. get_parameter("qbasis", qbasis,"b")) &
      call log_warning("get_parameters","'qbasis' missing. Using reciprocal lattice.")
    !---------------------------------- Magnetic Self-consistency ----------------------------------
    if(.not. get_parameter("skipsc", skipsc, .false.)) &
      call log_warning("get_parameters","'skipsc' missing. Using default value: .false.")
    if(.not. get_parameter("scfile", scfile))&
      call log_warning("get_parameters","'scfile' missing. Using none.")
    !======================================== Tight-Binding ========================================

    if(.not. get_parameter("tbmode", tbmode, 1)) &
      call log_warning("get_parameters", "'tbmode' missing. Using default value 1 (SK parameters).")
    if(.not. get_parameter("nOrb",nOrb,9)) &
      call log_warning("get_parameters", "'nOrb' missing. Using default value: 9")
    nOrb2 = 2*nOrb
    !---------------------------------------- Slater-Koster ----------------------------------------
    if(tbmode == 1) then
      offset = 0
      dfttype = "S"

      ! if(.not. get_parameter("layers", layers, cnt)) call log_error("get_parameters", "'layers' missing.")
      ! if(cnt <= 0) call log_error("get_parameters", "'layers' No layers given.")
      ! Npl = cnt
      !
      ! if(get_parameter("Npl", i_vector, cnt)) then
      !   if(cnt < 1) call log_error("get_parameters","'Npl' doesn't contain any parameters.")
      !   Npl_i = i_vector(1)
      !   Npl_f = i_vector(1)
      !   if(cnt >= 2) Npl_f = i_vector(2)
      !   if(Npl_f < Npl_i) Npl_f = Npl_i
      !   deallocate(i_vector)
      !   if(Npl < Npl_f) call log_error("get_parameters", "'Npl' larger than amount of given layers")
      ! else
      !   call log_warning("get_parameters","'Npl' missing.")
      !   Npl_i = Npl
      !   Npl_f = Npl
      ! end if

      if(.not. get_parameter("fermi_layer", fermi_layer, 1)) &
        call log_warning("get_parameters", "'fermi_layer' not given. Using default value: fermi_layer = 1")

    !--------------------------------------------- DFT ---------------------------------------------
    else if(2 == tbmode) then
      stop "Not Implemented"
      !  offset = 1
      !  if(nstages /= 2) call log_error("get_parameters", "'tbmode' DFT Mode only supports nstages = 2")
       !
      !  if(.not. get_parameter("dfttype", dfttype)) call log_error("get_parameters","'dfttype' missing.")
       !
      !  if(.not. get_parameter("Npl", i_vector, cnt)) call log_error("get_parameters","'Npl' missing.")
      !  if(cnt < 1) call log_error("get_parameters","'Npl' doesn't contain any parameters.")
      !  Npl_i = i_vector(1)
      !  Npl_f = i_vector(1)
      !  if(cnt >= 2) Npl_f = i_vector(2)
      !  deallocate(i_vector)
       !
      !  ! Check number of planes
      !  if(Npl_f < Npl_i) Npl_f = Npl_i
       !
      !  if(.not. get_parameter("set1", set1)) call log_error("get_parameters","'set1' missing.")
       !
      !  if(.not. get_parameter("set2", set2)) call log_error("get_parameters","'set2' missing.")
       !
      !  if(get_parameter("addlayers", i_vector, cnt)) then
      !     if(cnt < 10) then
      !        addlayers(1:cnt) = i_vector(1:cnt)
      !        naddlayers = cnt
      !     else if(cnt >= 10) then
      !        addlayers(1:10) = i_vector(1:10)
      !        naddlayers = 10
      !     end if
      !  end if
      !  if(allocated(i_vector)) deallocate(i_vector)
       !
      !  ! Add 'naddlayers' to Npl
      !  if((naddlayers==1).and.(myrank==0)) write(outputunit,"('[get_parameters] WARNING: Added layers must include empty spheres! Only including one layer: naddlayers = ',i0)") naddlayers
      !  if((set1==9).or.(set2==9)) then
      !     naddlayers = 0
      !  end if
      !  if(naddlayers/=0) then
      !     Npl_i = Npl_i+naddlayers-1
      !     Npl_f = Npl_f+naddlayers-1
      !  end if
    else
       call log_error("get_parameters", "'tbmode' Unknown mode selected. (Choose either 1 or 2)")
    end if
    !==============================================================================================!
    ! REAL TIME PROPAGATION PARAMETERS
    if(itype==11) then
      if(.not. get_parameter("integration_time", integration_time)) &
        call log_error("get_parameters", "'integration_time' not found.")
      if(.not. get_parameter("step", step, integration_time/1.e4_dp )) &
        call log_warning("get_parameters", "'step' not found. Using default value: integration_time/1.d4")
      if(.not. get_parameter("sc_tol", sc_tol, 0.01_dp)) &
        call log_warning("get_parameters", "'sc_tol' not given. Using default value: sc_tol = 0.01_dp")
      if(.not. get_parameter("abs_tol", abs_tol, 1.e-3_dp)) &
        call log_warning("get_parameters", "'abs_tol' not given. Using default value: abs_tol = 0.001_dp")
      if(.not. get_parameter("rel_tol", rel_tol, 1.e-3_dp)) &
        call log_warning("get_parameters", "'rel_tol' not given. Using default value: rel_tol = 0.001_dp")
      if(.not. get_parameter("safe_factor", safe_factor, 0.9_dp)) &
        call log_warning("get_parameters", "'safe_factor' not given. Using default value: safe_factor = 0.9_dp")

      ! Reading ELECTRIC field variables
      if(.not. get_parameter("electric", lelectric,.false.)) &
        call log_warning("get_parameters", "'electric' not found. Electric field is not applied.")
      ! Electric field options:
      if(lelectric) then
        ! Pulse:
        if(.not. get_parameter("pulse_e", lpulse_e,.false.)) &
          call log_warning("get_parameters", "'pulse_e' not found. Oscillatory electric field is applied.")

        if(lpulse_e) then
          if(.not. get_parameter("npulse_e", npulse_e, 1)) &
            call log_warning("get_parameters","'npulse_e' missing. Using default npulse_e=1.")

          ! Allocating variables that depend on number of pulses
          allocate(polarization_e(npulse_e),polarization_vec_e(npulse_e,2,3),hE_0(npulse_e),hw_e(npulse_e),tau_e(npulse_e),delay_e(npulse_e))

          if(.not. get_parameter("polarization_e", s_vector, cnt)) then
            ! If 'polarization_e' is not given, try to find 'polarization_vec_ip_e' and/or 'polarization_vec_op_e'
            if(get_parameter("polarization_vec_ip_e", vector, cnt)) then
              if((cnt == 3).and.(npulse_e > 1)) then
                call log_warning("get_parameters","'polarization_vec_ip_e' has size 3. Using same in-phase polarization to all pulses.")
                forall(i=1:npulse_e) polarization_vec_e(i,1,:) = vector(1:3)
              else
                if(cnt < 3*npulse_e) call log_error("get_parameters","'polarization_vec_ip_e' has wrong size (size 3*npulse_e=" // trim(itos(3*npulse_e)) // " required).")
                forall(i=1:npulse_e) polarization_vec_e(i,1,:) = vector((i-1)*3+1:(i-1)*3+3)
              end if
              deallocate(vector)
            end if
            if(get_parameter("polarization_vec_op_e", vector, cnt)) then
              if((cnt == 3).and.(npulse_e > 1)) then
                call log_warning("get_parameters","'polarization_vec_op_e' has size 3. Using same out-of-phase polarization to all pulses.")
                forall(i=1:npulse_e) polarization_vec_e(i,2,:) = vector(1:3)
              else
                if(cnt < 3*npulse_e) call log_error("get_parameters","'polarization_vec_op_e' has wrong size (size 3*npulse_e=" // trim(itos(3*npulse_e)) // " required).")
                forall(i=1:npulse_e) polarization_vec_e(i,2,:) = vector((i-1)*3+1:(i-1)*3+3)
              end if
              deallocate(vector)
            end if
            do i=1,npulse_e
              if( sum(abs(polarization_vec_e(i,:,:))) < 1.e-15_dp ) &
                call log_error("get_parameters", "'polarization_vec_e' is zero for pulse " // trim(itos(i)) // ". Use polarization_vec_ip_e/polarization_vec_op_e or polarization_e to set the polarization.")
            end do
          else
            if(cnt < npulse_e) call log_error("get_parameters","'polarization_e' has wrong size (size npulse_e=" // trim(itos(npulse_e)) // " required).")
            forall(i=1:npulse_e) polarization_e(i) = trim(s_vector(i))
            deallocate(s_vector)
            do i=1,npulse_e
              select case(polarization_e(i))
              case("x")
                polarization_vec_e(i,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_e(i,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
              case("y")
                polarization_vec_e(i,1,:) = [ 0._dp, 1._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_e(i,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
              case("z")
                polarization_vec_e(i,1,:) = [ 0._dp, 0._dp, 1._dp] ! in-phase     (cos(wt))
                polarization_vec_e(i,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
              case("p")
                polarization_vec_e(i,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_e(i,2,:) = [ 0._dp, 1._dp, 0._dp] ! out-of-phase (sin(wt))
              case("m")
                polarization_vec_e(i,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_e(i,2,:) = [ 0._dp,-1._dp, 0._dp] ! out-of-phase (sin(wt))
              case default
                call log_error("get_parameters", "Electric polarization 'polarization_e = '"// trim(polarization_e(i)) //" not found.")
              end select
            end do
          end if

          if(.not. get_parameter("hE_0", vector, cnt)) &
            call log_error("get_parameters", "'hE_0' not found.")
          if(cnt < npulse_e) call log_error("get_parameters","'hE_0' has wrong size (size npulse_e=" // trim(itos(npulse_e)) // " required).")
          hE_0(:) = vector(1:npulse_e)
          deallocate(vector)
          if(.not. get_parameter("hw_e", vector, cnt)) &
            call log_error("get_parameters", "'hw_e' not found.")
          if(cnt < npulse_e) call log_error("get_parameters","'hw_e' has wrong size (size npulse_e=" // trim(itos(npulse_e)) // " required).")
          hw_e(:) = vector(1:npulse_e)
          deallocate(vector)
          if(.not. get_parameter("tau_e", vector, cnt)) &
            call log_error("get_parameters", "'tau_e' not found.")
          if(cnt < npulse_e) call log_error("get_parameters","'tau_e' has wrong size (size npulse_e=" // trim(itos(npulse_e)) // " required).")
          tau_e(:) = vector(1:npulse_e)
          deallocate(vector)

          if(.not. get_parameter("delay_e", vector, cnt)) then
            call log_warning("get_parameters", "'delay_e' not found. Center of the pulses is located at t=tau_e/2.")
            delay_e(:) = 0._dp ! No extra delay
          else
            if(cnt < npulse_e) call log_error("get_parameters","'delay_e' has wrong size (size npulse_e=" // trim(itos(npulse_e)) // " required).")
            delay_e(1) = vector(1)
            do i=2,npulse_e
              delay_e(i) = delay_e(i-1)+vector(i)
            end do
            deallocate(vector)
          end if
        else ! If not pulse:
          allocate(polarization_e(1),polarization_vec_e(1,2,3),hE_0(1),hw_e(1))
          if(.not. get_parameter("polarization_e", s_vector, cnt)) then
            ! If 'polarization_e' is not given, try to find 'polarization_vec_ip_e' and/or 'polarization_vec_op_e'
            if(get_parameter("polarization_vec_ip_e", vector, cnt)) then
              if(cnt < 3) call log_error("get_parameters","'polarization_vec_ip_e' has wrong size (size 3 required).")
              polarization_vec_e(1,1,:) = vector(1:3)
              deallocate(vector)
            end if
            if(get_parameter("polarization_vec_op_e", vector, cnt)) then
              if(cnt < 3) call log_error("get_parameters","'polarization_vec_op_e' has wrong size (size 3 required).")
              polarization_vec_e(1,2,:) = vector(1:3)
              deallocate(vector)
            end if
            if( sum(abs(polarization_vec_e(1,:,:))) < 1.e-15_dp ) &
              call log_error("get_parameters", "'polarization_vec_e' is zero for oscillatory field. Use polarization_vec_ip_e/polarization_vec_op_e or polarization_e to set the polarization.")
          else
            if(cnt > 1) call log_warning("get_parameters","Only first element of 'polarization_e' will be used.")
            polarization_e(1) = trim(s_vector(1))
            deallocate(s_vector)
            select case(polarization_e(1))
            case("x")
              polarization_vec_e(1,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_e(1,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
            case("y")
              polarization_vec_e(1,1,:) = [ 0._dp, 1._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_e(1,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
            case("z")
              polarization_vec_e(1,1,:) = [ 0._dp, 0._dp, 1._dp] ! in-phase     (cos(wt))
              polarization_vec_e(1,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
            case("p")
              polarization_vec_e(1,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_e(1,2,:) = [ 0._dp, 1._dp, 0._dp] ! out-of-phase (sin(wt))
            case("m")
              polarization_vec_e(1,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_e(1,2,:) = [ 0._dp,-1._dp, 0._dp] ! out-of-phase (sin(wt))
            case default
              call log_error("get_parameters", "Electric polarization 'polarization_e = '"// trim(polarization_e(1)) //" not found.")
            end select
          end if

          if(.not. get_parameter("hE_0", vector, cnt)) &
            call log_error("get_parameters", "'hE_0' not found.")
          if(cnt > 1) call log_warning("get_parameters","Only first element of 'hE_0' will be used.")
          hE_0(1) = vector(1)
          deallocate(vector)
          if(.not. get_parameter("hw_e", vector, cnt)) &
            call log_error("get_parameters", "'hw_e' not found.")
          if(cnt > 1) call log_warning("get_parameters","Only first element of 'hw_e' will be used.")
          hw_e(1) = vector(1)
          deallocate(vector)
        end if
      end if

      ! Reading MAGNETIC field variables
      if(.not. get_parameter("magnetic", lmagnetic,.false.)) &
        call log_warning("get_parameters", "'magnetic' not found. Magnetic field is not applied.")
      ! Electric field options:
      if(lmagnetic) then
        ! Pulse:
        if(.not. get_parameter("pulse_m", lpulse_m,.false.)) &
          call log_warning("get_parameters", "'pulse_m' not found. Oscillatory magnetic field is applied.")

        if(lpulse_m) then
          if(.not. get_parameter("npulse_m", npulse_m, 1)) &
            call log_warning("get_parameters","'npulse_m' missing. Using default npulse_m=1.")

          ! Allocating variables that depend on number of pulses
          allocate(polarization_m(npulse_m),polarization_vec_m(npulse_m,2,3),hw1_m(npulse_m),hw_m(npulse_m),tau_m(npulse_m),delay_m(npulse_m))

          if(.not. get_parameter("polarization_m", s_vector, cnt)) then
            ! If 'polarization_m' is not given, try to find 'polarization_vec_ip_m' and/or 'polarization_vec_op_m'
            if(get_parameter("polarization_vec_ip_m", vector, cnt)) then
              if((cnt == 3).and.(npulse_m > 1)) then
                call log_warning("get_parameters","'polarization_vec_ip_m' has size 3. Using same in-phase polarization to all pulses.")
                forall(i=1:npulse_m) polarization_vec_m(i,1,:) = vector(1:3)
              else
                if(cnt < 3*npulse_m) call log_error("get_parameters","'polarization_vec_ip_m' has wrong size (size 3*npulse_m=" // trim(itos(3*npulse_m)) // " required).")
                forall(i=1:npulse_m) polarization_vec_m(i,1,:) = vector((i-1)*3+1:(i-1)*3+3)
              end if
              deallocate(vector)
            end if
            if(get_parameter("polarization_vec_op_m", vector, cnt)) then
              if((cnt == 3).and.(npulse_m > 1)) then
                call log_warning("get_parameters","'polarization_vec_op_m' has size 3. Using same out-of-phase polarization to all pulses.")
                forall(i=1:npulse_m) polarization_vec_m(i,2,:) = vector(1:3)
              else
                if(cnt < 3*npulse_m) call log_error("get_parameters","'polarization_vec_op_m' has wrong size (size 3*npulse_m=" // trim(itos(3*npulse_m)) // " required).")
                forall(i=1:npulse_m) polarization_vec_m(i,2,:) = vector((i-1)*3+1:(i-1)*3+3)
              end if
              deallocate(vector)
            end if
            do i=1,npulse_m
              if( sum(abs(polarization_vec_m(i,:,:))) < 1.e-15_dp ) &
                call log_error("get_parameters", "'polarization_vec_m' is zero for pulse " // trim(itos(i)) // ". Use polarization_vec_m_ip/polarization_vec_m_op or polarization_m to set the polarization.")
            end do
          else
            if(cnt < npulse_m) call log_error("get_parameters","'polarization_m' has wrong size (size npulse_m=" // trim(itos(npulse_m)) // " required).")
            forall(i=1:npulse_m) polarization_m(i) = trim(s_vector(i))
            deallocate(s_vector)
            do i=1,npulse_m
              select case(polarization_m(i))
              case("x")
                polarization_vec_m(i,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_m(i,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
              case("y")
                polarization_vec_m(i,1,:) = [ 0._dp, 1._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_m(i,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
              case("z")
                polarization_vec_m(i,1,:) = [ 0._dp, 0._dp, 1._dp] ! in-phase     (cos(wt))
                polarization_vec_m(i,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
              case("p")
                polarization_vec_m(i,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_m(i,2,:) = [ 0._dp, 1._dp, 0._dp] ! out-of-phase (sin(wt))
              case("m")
                polarization_vec_m(i,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
                polarization_vec_m(i,2,:) = [ 0._dp,-1._dp, 0._dp] ! out-of-phase (sin(wt))
              case default
                call log_error("get_parameters", "Magnetic polarization 'polarization_m = '"// trim(polarization_m(i)) //" not found.")
              end select
            end do
          end if

          if(.not. get_parameter("hw1_m", vector, cnt)) &
            call log_error("get_parameters", "'hw1_m' not found.")
          if(cnt < npulse_m) call log_error("get_parameters","'hw1_m' has wrong size (size npulse_m=" // trim(itos(npulse_m)) // " required).")
          hE_0(:) = vector(1:npulse_m)
          deallocate(vector)
          if(.not. get_parameter("hw_m", vector, cnt)) &
            call log_error("get_parameters", "'hw_m' not found.")
          if(cnt < npulse_m) call log_error("get_parameters","'hw_m' has wrong size (size npulse_m=" // trim(itos(npulse_m)) // " required).")
          hw_m(:) = vector(1:npulse_m)
          deallocate(vector)
          if(.not. get_parameter("tau_m", vector, cnt)) &
            call log_error("get_parameters", "'tau_m' not found.")
          if(cnt < npulse_m) call log_error("get_parameters","'tau_m' has wrong size (size npulse_m=" // trim(itos(npulse_m)) // " required).")
          tau_m(:) = vector(1:npulse_m)
          deallocate(vector)

          if(.not. get_parameter("delay_m", vector, cnt)) then
            call log_warning("get_parameters", "'delay_m' not found. Center of the pulses is located at t=tau_m/2.")
            delay_m(:) = 0._dp ! No extra delay
          else
            if(cnt < npulse_m) call log_error("get_parameters","'delay_m' has wrong size (size npulse_m=" // trim(itos(npulse_m)) // " required).")
            delay_m(1) = vector(1)
            do i=2,npulse_m
              delay_m(i) = delay_m(i-1)+vector(i)
            end do
            deallocate(vector)
          end if
        else ! If not pulse:
          allocate(polarization_m(1),polarization_vec_m(1,2,3),hw1_m(1),hw_m(1))
          if(.not. get_parameter("polarization_m", s_vector, cnt)) then
            ! If 'polarization_m' is not given, try to find 'polarization_vec_ip_m' and/or 'polarization_vec_op_m'
            if(get_parameter("polarization_vec_ip_m", vector, cnt)) then
              if(cnt < 3) call log_error("get_parameters","'polarization_vec_ip_m' has wrong size (size 3 required).")
              polarization_vec_m(1,1,:) = vector(1:3)
              deallocate(vector)
            end if
            if(get_parameter("polarization_vec_op_m", vector, cnt)) then
              if(cnt < 3) call log_error("get_parameters","'polarization_vec_op_m' has wrong size (size 3 required).")
              polarization_vec_m(1,2,:) = vector(1:3)
              deallocate(vector)
            end if
            if( sum(abs(polarization_vec_m(1,:,:))) < 1.e-15_dp ) &
              call log_error("get_parameters", "'polarization_vec_m' is zero for oscillatory field. Use polarization_vec_m_ip/polarization_vec_m_op or polarization_m to set the polarization.")
          else
            if(cnt > 1) call log_warning("get_parameters","Only first element of 'polarization_m' will be used.")
            polarization_m(1) = trim(s_vector(1))
            deallocate(s_vector)
            select case(polarization_m(1))
            case("x")
              polarization_vec_m(1,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_m(1,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
            case("y")
              polarization_vec_m(1,1,:) = [ 0._dp, 1._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_m(1,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
            case("z")
              polarization_vec_m(1,1,:) = [ 0._dp, 0._dp, 1._dp] ! in-phase     (cos(wt))
              polarization_vec_m(1,2,:) = [ 0._dp, 0._dp, 0._dp] ! out-of-phase (sin(wt))
            case("p")
              polarization_vec_m(1,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_m(1,2,:) = [ 0._dp, 1._dp, 0._dp] ! out-of-phase (sin(wt))
            case("m")
              polarization_vec_m(1,1,:) = [ 1._dp, 0._dp, 0._dp] ! in-phase     (cos(wt))
              polarization_vec_m(1,2,:) = [ 0._dp,-1._dp, 0._dp] ! out-of-phase (sin(wt))
            case default
              call log_error("get_parameters", "Magnetic polarization 'polarization_m = '"// trim(polarization_m(1)) //" not found.")
            end select
          end if

          if(.not. get_parameter("hw1_m", vector, cnt)) &
            call log_error("get_parameters", "'hw1_m' not found.")
          if(cnt > 1) call log_warning("get_parameters","Only first element of 'hw1_m' will be used.")
          hw1_m(1) = vector(1)
          deallocate(vector)
          if(.not. get_parameter("hw_m", vector, cnt)) &
            call log_error("get_parameters", "'hw_m' not found.")
          if(cnt > 1) call log_warning("get_parameters","Only first element of 'hw_m' will be used.")
          hw_m(1) = vector(1)
          deallocate(vector)
        end if
      end if

      ! Setup output variable for filename
      if(lmagnetic) then
        output%time_field = "_magfield"
        if(lpulse_m) then
          output%time_field = trim(output%time_field) // "pulse"
        end if
      end if
      if(lelectric) then
        output%time_field = trim(output%time_field) // "_efield"
        if(lpulse_e) then
          output%time_field = trim(output%time_field) // "pulse"
        end if
      end if

    end if ! itype=11 - time_prop_input
    !==============================================================================================!
    if(.not. get_parameter("suffix", output%suffix)) &
      call log_warning("get_parameters","'suffix' missing. Using none.")
    if(.not. get_parameter("parField", parField, 1)) &
      call log_warning("get_parameters","'parField' missing. Using default value: 1")
    if(.not. get_parameter("parFreq", parFreq, 1)) &
      call log_warning("get_parameters","'parFreq' missing. Using default value: 1")
    if(myrank == 0) then
      if(.not. disable_input_logging()) &
          call log_warning("get_parameters", "Could not disable logging.")
    end if
    !==============================================================================================!
    if(myrank==0) &
      write(output%unit,"('[get_parameters] Finished reading from ""',a,'"" file')") trim(filename)


    !-------------------------------------------------------------------------------
    !*********** User manual additions / modifications in the input file **********!
    !     Npl_i  = 4
    !     Npl_f = 4
    !     nkpt = 6
    !     SOC = .true.
    !     magaxis = "5"
    !     runoptions = trim(runoptions)
    !     scfile = "results/selfconsistency/selfconsistency_Npl=4_dfttype=T_parts=2_U= 0.7E-01_hwa= 0.00E+00_hwt= 0.00E+00_hwp= 0.00E+00_nkpt=6_eta= 0.5E-03.dat"
    !-------------------------------------------------------------------------------
    ! Some consistency checks
    if((renorm).and.((renormnb<n0sc1).or.(renormnb>n0sc2))) then
        call log_error("get_parameters", "Invalid neighbor for renormalization: " // trim(itos(renormnb)) // ". Choose a value between " // trim(itos(n0sc1)) // " and " // trim(itos(n0sc2)) // ".")
    end if
    if(skip_steps<0) then
       if(myrank==0) write(output%unit,"('[get_parameters] Invalid number of energy steps to skip: ',i0)") skip_steps
       call MPI_Finalize(ierr)
       stop
    end if
    if(skip_steps_hw<0) then
       if(myrank==0) write(output%unit,"('[get_parameters] Invalid number of field steps to skip: ',i0)") skip_steps_hw
       call MPI_Finalize(ierr)
       stop
    end if
    if((lhfresponses).and.(itype==7).and.(myrank==0)) write(output%unit,"('[get_parameters] Susceptibility calculations already include HF responses. Ignoring ""hfresponses"" runoption')")
    ! Adjusting zeeman energy to Ry or eV
    tesla = tesla*ry2ev

    ! Turning off renormalization for non-current calculations
    if(itype/=8) renorm = .false.

    ! Set string for HF Responses
    if(lhfresponses) then
      output%hfr = "_HF"
    else
      output%hfr = ""
    end if

    ! Preparing dc-limit calculation
    ! if(itype==9) call prepare_dclimit() !TODO: Re-Include

    pn1=parts*n1gl
    pn2=parts3*n3gl
    pnt=pn1+pn2
  end subroutine get_parameters

  subroutine iowrite(s)
    use mod_mpi_pars
    use mod_parameters
    use mod_magnet
    use mod_System,            only: System_type, n0sc1, n0sc2
    use mod_BrillouinZone,     only: BZ => realBZ
    use mod_fermi_surface,     only: lfs_loop, fs_energy_i, fs_energy_f, fs_energy_npts
    use mod_SOC,               only: SOC, socscale
    use EnergyIntegration,     only: parts, parts3, n1gl, n3gl
    use ElectricField,         only: ElectricFieldMode, ElectricFieldVector, EFt, EFp
    use AdaptiveMesh,          only: minimumBZmesh
    use mod_superconductivity, only: lsupercond
    use mod_imRK4_parameters,  only: integration_time, sc_tol, hE_0, hw1_m, hw_e, hw_m, tau_e, &
                                    tau_m, delay_e, delay_m, lelectric, lmagnetic, lpulse_e, npulse_e, lpulse_m, npulse_m, &
                                    polarization_vec_e, polarization_vec_m, abs_tol, rel_tol, safe_factor
    !$ use omp_lib
    implicit none
    type(System_type), intent(in) :: s
    integer :: i,j
#ifdef _OPENMP
    write(output%unit_loop,"(10x,'Running on ',i0,' MPI process(es) WITH ',i0,' openMP')") numprocs, omp_get_max_threads()
#else
    write(output%unit_loop,"(10x,'Running on ',i0,' MPI process(es) WITHOUT openMP')") numprocs
#endif
    write(output%unit_loop,"('|------------------------------- PARAMETERS: -------------------------------|')")
    write(output%unit_loop,"(10x,'nAtoms = ',i0)") s%nAtoms
    ! write(output%unit_loop,"(1x,'DFT parameters: ')", advance='no')
    ! dft_type: select case (dfttype)
    ! case ("T")
    !    write(output%unit_loop,"('Tight-binding basis')")
    ! case ("O")
    !    write(output%unit_loop,"('Orthogonal basis')")
    ! end select dft_type
    if(SOC) then
      write(output%unit_loop,"(1x,'Spin Orbit Coupling: ACTIVATED')")
      write(output%unit_loop,"(5x,'socscale =',es9.2)") socscale
    else
      write(output%unit_loop,"(1x,'Spin Orbit Coupling: DEACTIVATED')")
    end if

    write(output%unit_loop,"(1x,'Electric field direction: ')", advance='no')
    select case(ElectricFieldMode)
    case(-3)
      write(output%unit_loop,"('Spherical theta=',f7.3,' phi=',f7.3)") EFt, EFp
    case(-2)
      write(output%unit_loop,"('Bravais ')")
    case(-1)
      write(output%unit_loop,"('Cartesian')")
    case(1:99)
      write(output%unit_loop, "('Neighbor ',i0)") ElectricFieldMode
    end select
    write(output%unit_loop,"(1x,'Direction: ')", advance='no')
    write(output%unit_loop,"('E = (',f6.3,',',f6.3,',',f6.3,')')") (ElectricFieldVector(i), i=1,3)
    if(renorm) then
      write(output%unit_loop,"(1x,'Current renormalization: ACTIVATED')")
      write(output%unit_loop,"(5x,'renormnb = ',i0)") renormnb
    else
      write(output%unit_loop,"(1x,'Current renormalization: DEACTIVATED')")
    end if
    write(output%unit_loop,"(9x,'nkpt = ',i0,' : ',i0,' x ',i0,' x ',i0)") BZ%nkpt, BZ%nkpt_x, BZ%nkpt_y, BZ%nkpt_z
    if(minimumBZmesh/=1000) write(output%unit_loop,"(9x,'minimumBZmesh = ',i0)") minimumBZmesh
    write(output%unit_loop,"(8x,'parts = ',i0,'x',i0)") parts,n1gl
    write(output%unit_loop,"(7x,'parts3 = ',i0,'x',i0)") parts3,n3gl
    write(output%unit_loop,"(10x,'eta =',es9.2)") eta
    if(lsuperCond) then
      write(output%unit_loop,"(1x,'Superconductivity: ACTIVATED')")
    else
      write(output%unit_loop,"(1x,'Superconductivity: DEACTIVATED')")
    end if

    if(lfield) then
      write(output%unit_loop,"(1x,'Static magnetic field: ACTIVATED')")
      write(output%unit_loop,"(10x,'hwx =',es9.2,5x,'|',5x,'hwa =',es9.2)") hwx,hw_list(hw_count,1)
      write(output%unit_loop,"(10x,'hwy =',es9.2,5x,'|',5x,'hwt =',f7.2)") hwy,hw_list(hw_count,2)
      write(output%unit_loop,"(10x,'hwz =',es9.2,5x,'|',5x,'hwp =',f7.2)") hwz,hw_list(hw_count,3)
    else
      write(output%unit_loop,"(1x,'Static magnetic field: DEACTIVATED')")
    end if
    if(runoptions/="") write(output%unit_loop,"(6x,'Activated options:',/,4x,a)") trim(runoptions)

    write(output%unit_loop,"('|------------------------------ TO CALCULATE: ------------------------------|')")
    write_itype: select case (itype)
    case (0)
      write(output%unit_loop,"(1x,'Test before SC')")
      write(output%unit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
      write(output%unit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
      write(output%unit_loop,"(9x,'emin =',es9.2)") emin
      write(output%unit_loop,"(9x,'emax =',es9.2)") emax
      write(output%unit_loop,"(1x,'Number of points to calculate: ',i0)") nEner1
    case (1)
      write(output%unit_loop,"(1x,'Self-consistency only')")
    case (2)
      write(output%unit_loop,"(1x,'Test after SC')")
      write(output%unit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
      write(output%unit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
      write(output%unit_loop,"(9x,'emin =',es9.2)") emin
      write(output%unit_loop,"(9x,'emax =',es9.2)") emax
      write(output%unit_loop,"(1x,'Number of points to calculate: ',i0)") nEner1
    case (3)
      write(output%unit_loop,"(1x,'LDOS and exchange interactions as a function of energy')")
      write(output%unit_loop,"(9x,'emin =',es9.2)") emin
      write(output%unit_loop,"(9x,'emax =',es9.2)") emax
      write(output%unit_loop,"(1x,'Number of points to calculate: ',i0)") nEner1
    case (4)
      write(output%unit_loop,"(1x,'Band structure')")
      write(output%unit_loop,"(2x,'Path along BZ: ',10(a,1x))") (trim(adjustl(bands(i))), i = 1,band_cnt)
      write(output%unit_loop,"(2x,'Number of wave vectors to calculate: ',i0)") nQvec1
    case (5)
      if(lfs_loop) then
        write(output%unit_loop,"(1x,'Charge, spin and orbital density at energy surface')")
        write(output%unit_loop,"(9x,'emin =',es9.2)") fs_energy_i
        write(output%unit_loop,"(9x,'emax =',es9.2)") fs_energy_f
        write(output%unit_loop,"(1x,'Number of points to calculate: ',i0)") fs_energy_npts
      else
        write(output%unit_loop,"(1x,'Charge, spin and orbital density at Fermi surface')")
      end if

    case (6)
      write(output%unit_loop,"(1x,'Exhange interactions and anisotropies (full tensor)')")
      if(s%nAtoms==1) write(output%unit_loop,"(1x,'Only 1 atom in the unit cell: calculating only anisotropies')")
      !write(outputunit_loop,"(8x,'from Npl = ',i0,' to ',i0)") Npl_i,Npl_f
    case (7)
      write(output%unit_loop,"(1x,'Local susceptibility as a function of energy')")
      write(output%unit_loop,"(2x,'Path along BZ: ',10(a,1x))") (trim(adjustl(bands(i))), i = 1,band_cnt)
      write(output%unit_loop,"(2x,'Number of wave vectors to calculate: ',i0)") nQvec1
      write(output%unit_loop,"(9x,'emin =',es9.2)") emin
      write(output%unit_loop,"(9x,'emax =',es9.2)") emax
      !write(output%unit_loop,"(1x,i0,' points divided into ',i0,' steps of size',es10.3,' each calculating ',i0,' points')") nEner1,MPIsteps,MPIdelta,MPIpts
    case (8)
      write(output%unit_loop,"(1x,'Parallel currents, disturbances and local susc. as a function of energy')")
      write(output%unit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
      write(output%unit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
      write(output%unit_loop,"(9x,'emin =',es9.2)") emin
      write(output%unit_loop,"(9x,'emax =',es9.2)") emax
      !write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps of energy size',es10.3,' each calculating ',i0,' points')") total_hw_npt1*nEner1,MPIsteps*MPIsteps_hw,MPIdelta,MPIpts_hw*MPIpts
    case (9)
      write(output%unit_loop,"(1x,'dc limit calculations as a function of ',a)") trim(dcfield(dcfield_dependence))
      write(output%unit_loop,"(1x,'e =',es9.2)") emin
      write(output%unit_loop,"(1x,'hwa_min =',es9.2)") hw_list(1,1)
      write(output%unit_loop,"(1x,'hwa_max =',es9.2)") hw_list(total_hw_npt1,1)
      write(output%unit_loop,"(1x,'hwt_min =',f7.2)") hw_list(1,2)
      write(output%unit_loop,"(1x,'hwt_max =',f7.2)") hw_list(total_hw_npt1,2)
      write(output%unit_loop,"(1x,'hwp_min =',f7.2)") hw_list(1,3)
      write(output%unit_loop,"(1x,'hwp_max =',f7.2)") hw_list(total_hw_npt1,3)
      !write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps, each calculating ',i0,' points')") total_hw_npt1*nEner1,MPIsteps*MPIsteps_hw,MPIpts_hw*MPIpts
    case (11)
      write(output%unit_loop, fmt="('Time propagation:')" )

      write(output%unit_loop,"(1x,'integration_time =',es9.2)") integration_time
      write(output%unit_loop,"(1x,'     sc_tol =',es9.2)") sc_tol
      write(output%unit_loop,"(1x,'    abs_tol =',es9.2)") abs_tol
      write(output%unit_loop,"(1x,'    rel_tol =',es9.2)") rel_tol
      write(output%unit_loop,"(1x,'safe_factor =',es9.2)") safe_factor
      if(lelectric) then
        write(output%unit_loop, fmt="('Electric field: ON')" )
        if(lpulse_e) then
          write(output%unit_loop, fmt="(i0, ' electric pulse(s):')" ) npulse_e
          do i = 1,npulse_e
            write(output%unit_loop,"(1x,'        hE_0 =',es9.2)") hE_0(i)
            write(output%unit_loop,"(1x,'        hw_e =',es9.2)") hw_e(i)
            write(output%unit_loop,"(1x,'Polarization = (',f6.3,',',f6.3,',',f6.3,') in-phase    ')") (polarization_vec_e(i,1,j), j=1,3)
            write(output%unit_loop,"(1x,'             = (',f6.3,',',f6.3,',',f6.3,') out-of-phase')") (polarization_vec_e(i,2,j), j=1,3)
            write(output%unit_loop,"(1x,'       tau_e =',es9.2)") tau_e(i)
            write(output%unit_loop,"(1x,'     delay_e =',es9.2)") delay_e(i)
          end do
        else
          write(output%unit_loop, fmt="(1x,'Oscillatory electric field:')" )
          write(output%unit_loop,"(1x,'        hE_0 =',es9.2)") hE_0(1)
          write(output%unit_loop,"(1x,'        hw_e =',es9.2)") hw_e(1)
          write(output%unit_loop,"(1x,'Polarization = (',f6.3,',',f6.3,',',f6.3,') in-phase    ')") (polarization_vec_e(1,1,j), j=1,3)
          write(output%unit_loop,"(1x,'             = (',f6.3,',',f6.3,',',f6.3,') out-of-phase')") (polarization_vec_e(1,2,j), j=1,3)
        end if
      end if
      if(lmagnetic) then
        write(output%unit_loop, fmt="('Magnetic field: ON')" )
        if(lpulse_m) then
          write(output%unit_loop, fmt="(i0, ' electric pulse(s):')" ) npulse_m
          do i = 1,npulse_m
            write(output%unit_loop,"(1x,'       hw1_m =',es9.2)") hw1_m(i)
            write(output%unit_loop,"(1x,'        hw_m =',es9.2)") hw_m(i)
            write(output%unit_loop,"(1x,'Polarization = (',f6.3,',',f6.3,',',f6.3,') in-phase    ')") (polarization_vec_m(i,1,j), j=1,3)
            write(output%unit_loop,"(1x,'             = (',f6.3,',',f6.3,',',f6.3,') out-of-phase')") (polarization_vec_m(i,2,j), j=1,3)
            write(output%unit_loop,"(1x,'       tau_m =',es9.2)") tau_m(i)
            write(output%unit_loop,"(1x,'     delay_m =',es9.2)") delay_m(i)
          end do
        else
          write(output%unit_loop, fmt="(1x,'Oscillatory electric field:')" )
          write(output%unit_loop,"(1x,'       hw1_m =',es9.2)") hw1_m(1)
          write(output%unit_loop,"(1x,'        hw_m =',es9.2)") hw_m(1)
          write(output%unit_loop,"(1x,'Polarization = (',f6.3,',',f6.3,',',f6.3,') in-phase    ')") (polarization_vec_m(1,1,j), j=1,3)
          write(output%unit_loop,"(1x,'             = (',f6.3,',',f6.3,',',f6.3,') out-of-phase')") (polarization_vec_m(1,2,j), j=1,3)
        end if
      end if
    end select write_itype
    write(output%unit_loop,"('|---------------------------------------------------------------------------|')")
  end subroutine iowrite


  ! Writing header for previously opened file of unit "unit"
  subroutine write_header(unit,title_line,Ef)
    use mod_kind, only: dp
    use mod_parameters, only: nQvec, nQvec1, bands, band_cnt, partial_length, itype
    implicit none
    integer,          intent(in)           :: unit
    character(len=*), intent(in)           :: title_line
    real(dp),     intent(in), optional :: Ef
    integer :: i

    ! LDOS header
    if(itype==2) then
      write(unit=unit, fmt="(a,2x,es16.9)") "# Ef ",Ef
    else ! Band structure and chi(q) header
      if(nQvec1/=1) then
        write(unit=unit, fmt="(a,2x,i0,2x,i0)") "# ", band_cnt, nQvec
        do i=1,band_cnt
          write(unit=unit, fmt="(a,2x,a,2x,es16.9)") "# ",trim(bands(i)), sum(partial_length(1:i))
        end do
        if(present(Ef)) write(unit=unit, fmt="(a,2x,es16.9)") "# Ef ",Ef
      end if
    end if

    write(unit=unit, fmt="(a)") trim(title_line)

  end subroutine write_header
end module mod_io
