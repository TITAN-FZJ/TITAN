module mod_io
  !! Subroutines to read input file and variables containing these parameters
  implicit none

contains



  function check_placeholders(initialFilename) result(modifiedFilename)
  !> This function will substitute pre-defined placeholders in the 'initialFilename'
  !> This is useful to output into different files when testing various values of parameters
    use mod_tools,      only: replaceStr
    use mod_parameters, only: kptotal_in,eta
    implicit none
    character(len=*)   :: initialFilename
    character(len=:), allocatable :: modifiedFilename
    character(len=100) :: temp

    modifiedFilename = initialFilename
    write (temp, "(i0)") kptotal_in
    modifiedFilename = replaceStr(modifiedFilename,"#nkpt",trim(temp))
    write (temp, "(es9.2)") eta
    modifiedFilename = replaceStr(modifiedFilename,"#eta",trim(adjustl(temp)))

  end function check_placeholders


  subroutine get_parameters(filename,s)
    use mod_kind,              only: dp,int64
    use mod_logging,           only: log_unit,logfile,log_store,log_message,log_warning,log_error
    use AtomTypes,             only: default_Orbs,default_nOrb
    use mod_input,             only: get_parameter,read_file,enable_input_logging,disable_input_logging
    use mod_parameters,        only: output,laddresults,lverbose,ldebug,lkpoints,&
                                     lpositions,lcreatefiles,lnolb,lhfresponses,&
                                     lnodiag,lsha,lcreatefolders,lwriteonscreen,runoptions,lsimplemix,&
                                     lcheckjac,llgtv,lsortfiles,leigenstates,lprintfieldonly,&
                                     itype,ry2ev,ltesla,eta,etap,addelectrons,dmax,emin,emax,&
                                     skip_steps,nEner,nEner1,nQvec,nQvec1,qbasis,renorm,renormnb,bands,band_cnt,&
                                     dfttype,parField,parFreq,kptotal_in,kp_in,&
                                     tbmode,fermi_layer,lfixEf,lEf_overwrite,Ef_overwrite,cluster_layers,qptotal_in,qp_in
    use mod_superconductivity, only: lsuperCond,superCond
    use mod_self_consistency,  only: lontheflysc,lnojac,lforceoccup,lrotatemag,skipsc,scfile,magbasis,mag_tol
    use mod_system,            only: System_type,n0sc1,n0sc2
    use mod_SOC,               only: SOC,socscale,llinearsoc,llineargfsoc
    use mod_magnet,            only: lfield,tesla,hwa_i,hwa_f,hwa_npts,hwa_npt1,hwt_i,hwt_f,&
                                     hwt_npts,hwt_npt1,hwp_i,hwp_f,hwp_npts,hwp_npt1,hwx,hwy,&
                                     hwz,hwscale,hwtrotate,hwprotate,skip_steps_hw,maxiter,cmix,lconstraining_field,constr_type
    use ElectricField,         only: ElectricFieldMode,ElectricFieldVector,EFp,EFt
    use EnergyIntegration,     only: parts,parts3,pn1,pn2,pnt,n1gl,n3gl
    use mod_tools,             only: itos,rtos,stoi,vec_norm,is_numeric
    use adaptiveMesh,          only: minimumBZmesh
    use mod_fermi_surface,     only: lfs_loop,fs_energy_npts,fs_energy_npt1,fs_energy_i,fs_energy_f
    use mod_mpi_pars,          only: myrank,ierr
!#ifndef _OLDMPI
!    use mod_mpi_pars,          only: MPI_Finalize
!#endif
    use mod_time_propagator,   only: integration_time,sc_tol,step,hE_0,hw1_m,hw_e,hw_m,tau_e,&
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
    character(len=20)  :: tmp_string
!#ifdef _OLDMPI
    external :: MPI_Finalize
!#endif

    intrinsic :: findloc

    if(.not. read_file(filename)) &
      call log_error("get_parameters", "File " // trim(filename) // " not found!")

    ! Logging of input variables in 'parameters.in' file
    if(myrank == 0) then
      if(.not. enable_input_logging(logfile)) &
        call log_warning("get_parameters", "couldn't enable logging.")
    end if
    !===============================================================================================
    !============= System configuration (Lattice + Reciprocal lattice) =============================
    !===============================================================================================
    ! Done before opening the output file to be able to use these variables in the output filename
    ! For this, the variable 'log_store' is used to temporarily store the output
    log_store = ""
    if(.not. get_parameter("nn_stages", s%nStages,2)) &
      call log_warning("get_parameters","'nn_stages' missing. Using default value: 2",log_store)
    if(.not. get_parameter("relTol", s%relTol,0.05_dp)) &
      call log_warning("get_parameters","'relTol' missing. Using default value: 0.05",log_store)
    if(.not. get_parameter("addelectrons", addelectrons, 0._dp)) &
      call log_warning("get_parameters","'addelectrons' missing. Using default value: 0.0",log_store)
    if(.not. get_parameter("sysdim", s%isysdim, 3)) &
      call log_warning("get_parameters", "'sysdim' missing. Using default value: 3",log_store)
    if(.not. get_parameter("nkpt", i_vector,cnt)) &
      call log_error("get_parameters","'nkpt' missing.")
    if(cnt == s%isysdim) then
      kp_in(1:s%isysdim) = int(i_vector(1:s%isysdim),kind(kp_in(1)))
      kptotal_in = product(kp_in(1:s%isysdim))
    else if(cnt == 1) then
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
    else
      call log_error("get_parameter", "'nkpt' has wrong size (expected 1 or isysdim).")
    end if

    if(.not. get_parameter("minimumBZmesh", i_vector,cnt)) then
      call log_warning("get_parameters", "'minimumBZmesh' missing. Using default value: 1000",log_store)
      minimumBZmesh = 1000
    else
      if(cnt == 1) then
        minimumBZmesh = int(i_vector(1),kind(minimumBZmesh))
      else
        call log_error("get_parameter", "'minimumBZmesh' has wrong size (expected 1).")
      end if
    end if
    if(minimumBZmesh>kptotal_in) then
      call log_warning("get_parameters", "'minimumBZmesh' is larger than the largest mesh. Using " // itos(kptotal_in))
      minimumBZmesh = kptotal_in
    end if
    if(.not. get_parameter("eta", eta)) &
      call log_error("get_parameters","'eta' missing.")
    if(.not. get_parameter("etap", etap, eta)) &
      call log_warning("get_parameters", "'etap' not found. Using default value: eta = " // trim(rtos(eta,"(es8.1)")),log_store )


    if(.not. get_parameter("output", output%file)) &
      call log_error("get_parameters", "Output filename not given!")

    ! Changing placeholders to variable values, if they exist
    output%file = trim(check_placeholders(output%file))

    if(myrank==0) then
      !! Create folder for output file, if necessary
      call execute_command_line('mkdir -p $(dirname "'// trim(output%file) // '")')
      !! Open output file
      open(unit=output%unit, file=trim(output%file), status='replace')
    end if

    log_unit = .true.

! Print the Git version (VERSION is defined via CMake macro and defined with compiler flag -DVERSION)
#if defined(VERSION)
    if(myrank==0) write(output%unit,"('Git commit: ',a)") VERSION
#else
    if(myrank==0) write(output%unit,"('Git commit: unknown')")
#endif

    if(myrank==0) &
      write(output%unit,"('[get_parameters] Reading parameters from ""',a,'"" file...')") trim(filename)

    ! Logging stored messages on 'log_store'
    call log_message("", log_store)

    !===============================================================================================
    !------------------------------------- Type of Calculation -------------------------------------
    if(.not. get_parameter("itype", itype)) &
      call log_error("get_parameters","'itype' missing.")
    if(.not. get_parameter("Options", s_vector, cnt)) then
      call log_warning("get_parameters","'Options' missing.")
    else
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
        case ("fixEf")
          lfixEf = .true.
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
    end if
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
      call log_warning("get_parameters","'FIELD' missing. Using default value: .false.")

    if(lfield) then
      if(.not. get_parameter("hwa", vector, cnt)) then
        call log_warning("get_parameters","Field is active but 'hwa' missing. Reading cartesian components...")
        if(.not. get_parameter("hwx", hwx)) &
          call log_error("get_parameters","'hwx' missing.")
        if(.not. get_parameter("hwy", hwy)) &
          call log_error("get_parameters","'hwy' missing.")
        if(.not. get_parameter("hwz", hwz)) &
          call log_error("get_parameters","'hwz' missing.")
      else
        select case (cnt)
        case(1)
          hwa_i = vector(1)
        case(3)
          hwa_f = vector(2)
          hwa_npts = int(vector(3))
        case default
          call log_error("get_parameters","'hwa' contain invalid number of parameters: " // itos(cnt) // ".")
        end select
        deallocate(vector)
        hwa_npt1 = hwa_npts + 1

        if(.not. get_parameter("hwt", vector, cnt)) &
          call log_error("get_parameters","'hwt' missing.")
        select case (cnt)
        case(1)
          hwt_i = vector(1)
        case(3)
          hwt_f = vector(2)
          hwt_npts = int(vector(3))
        case default
          call log_error("get_parameters","'hwt' contain invalid number of parameters: " // itos(cnt) // ".")
        end select
        deallocate(vector)
        hwt_npt1 = hwt_npts + 1

        if(.not. get_parameter("hwp", vector, cnt)) &
          call log_error("get_parameters","'hwp' missing.")
        select case (cnt)
        case(1)
          hwp_i = vector(1)
        case(3)
          hwp_f = vector(2)
          hwp_npts = int(vector(3))
        case default
          call log_error("get_parameters","'hwp' contain invalid number of parameters: " // itos(cnt) // ".")
        end select
        deallocate(vector)
        hwp_npt1 = hwp_npts + 1
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
    !--------------------------------------- Constraining Field ------------------------------------
    if(.not. get_parameter("constraining_field", lconstraining_field, .false.)) &
      call log_warning("get_parameters","'lconstraining_field' missing. Using default value: .false.")
    if(lconstraining_field) then
      if(magbasis=="") &
        call log_error("get_parameters","'magbasis' is required when constraining_field = T. ")
      if(.not. get_parameter("constr_type", constr_type)) &
        call log_error("get_parameters","'constr_type' is required when constraining_field = T.  Use one of the following: " // NEW_line('A') // &
                "1 - Transverse constraining field" // NEW_line('A') // &
                "2 - Full constranining field")
      select case(constr_type)
      case(1)
        call log_message("get_parameters","Transverse constraining field chosen (constr_type=1): a field is applied to keep the direction of m fixed (not its length).")
        if(.not. get_parameter("cmix", cmix, 1.e-2_dp)) &
          call log_warning("get_parameters","'cmix' missing. Using default value: 0.01")
        ! call read_initial_bconstr("initial_bconstr","c",s%nAtoms,initial_bconstr)
      case(2)
        call log_message("get_parameters","Full constraining field chosen (constr_type=2): the magnetic moment is completely fixed by the application of a constraining field.")
      case default
        call log_error("get_parameters","'constr_type' not recognized.  Use one of the following: " // NEW_line('A') // &
                "1 - Transverse constraining field (to fix direction of M)" // NEW_line('A') // &
                "2 - Full constraining field (to keep M fixed)")
      end select
    end if
    !------------------------------ Superconductivity Variables ------------------------------------
    if(.not. get_parameter("superCond", lsuperCond, .false.)) &
      call log_warning("get_parameters","'superCond' missing. Not using superconductivity.")
    superCond = merge(2,1,lsuperCond)

    if(itype==12) then
      if(.not. get_parameter("cluster_layers", cluster_layers)) then
        call log_warning("get_parameters","'cluster_layers' missing. Using default value 2.")
      end if
      if(.not. get_parameter("nqpt", i_vector,cnt)) then
        call log_warning("get_parameters","'nqpt' missing.  Using default value equal 'nkpt'.")
        qp_in(:) = kp_in(:)
        qptotal_in = kptotal_in
      else
        if(cnt == s%isysdim) then
          qp_in(1:s%isysdim) = int(i_vector(1:s%isysdim),kind(qp_in(1)))
          qptotal_in = product(qp_in(1:s%isysdim))
        else if(cnt == 1) then
          qptotal_in = int( i_vector(1), kind(qptotal_in) )
          select case(s%isysdim)
          case(3)
            qp_in(:)   = ceiling((dble(qptotal_in))**(0.333333333333333_dp), kind(qp_in(1)) )
            qptotal_in = int( qp_in(1) * qp_in(2) * qp_in(3), kind(qptotal_in) )
          case(2)
            qp_in(1:2) = ceiling((dble(qptotal_in))**(0.5_dp), kind(qp_in(1)) )
            qp_in(3)   = 1
            qptotal_in = int( qp_in(1) * qp_in(2), kind(qptotal_in) )
          case default
            qp_in(1)   = ceiling((dble(qptotal_in)), kind(qp_in(1)) )
            qp_in(2:3) = 1
            qptotal_in = int( qp_in(1), kind(qptotal_in) )
          end select
        else
          call log_error("get_parameter", "'nqpt' has wrong size (expected 1 or isysdim).")
        end if
      end if
    end if
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
    if((itype >= 6).and.(itype <= 9)) then
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
      ! Path to calculate band structure (cannot be single point in this case)
      if(.not. get_parameter("band", bands, band_cnt)) &
        call log_error("get_parameters", "'band' missing.")
      if((band_cnt < 2).or.(nQvec < 2)) call log_error("get_parameters", "Need at least two points for Band Structure")
    endif
    if(.not. get_parameter("qbasis", qbasis,"b")) &
      call log_warning("get_parameters","'qbasis' missing. Using reciprocal lattice.")
    !---------------------------------- Magnetic Self-consistency ----------------------------------
    if(.not. get_parameter("skipsc", skipsc, .false.)) &
      call log_warning("get_parameters","'skipsc' missing. Using default value: .false.")
    if(.not. get_parameter("maxiter", maxiter, 99999999)) &
      call log_warning("get_parameters","'maxiter' missing. Using default value: 99999999 (no limit)")
    if(.not. get_parameter("scfile", scfile))&
      call log_warning("get_parameters","'scfile' missing. Using none.")
    !======================================== Tight-Binding ========================================
    if(.not. get_parameter("tbmode", tbmode, 1)) &
      call log_warning("get_parameters", "'tbmode' missing. Using default value 1 (SK parameters).")

    if(.not.get_parameter("orbitals", s_vector, cnt)) then
      call log_warning("get_parameters","'orbitals' missing. Using the default 9 orbitals.")
      cnt = default_nOrb
      allocate(s_vector(cnt))
      s_vector(:) = default_Orbs(:)
    end if

    call get_orbitals("the whole system",s_vector,s%nOrb,s%nOrb2,s%nOrb2sc,s%nsOrb,s%npOrb,s%ndOrb,s%Orbs,s%sOrbs,s%pOrbs,s%dOrbs)
    deallocate(s_vector)

    if(.not. get_parameter("fermi_layer", fermi_layer, 1)) &
      call log_warning("get_parameters", "'fermi_layer' not given. Using default value: fermi_layer = 1")
    !---------------------------------------- Slater-Koster ----------------------------------------
    if(tbmode == 1) then
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

    !--------------------------------------------- DFT ---------------------------------------------
    else if(tbmode == 2) then
      dfttype = "D"
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
    if(get_parameter("Ef", Ef_overwrite)) then
      call log_warning("get_parameters", "Ef = " // trim(rtos(Ef_overwrite,"(es8.1)")) //" given. Overwriting (initial) Fermi energy. (use fixEf runoption to fix)")
      lEf_overwrite = .true.
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
            do i=1,npulse_e
              polarization_e(i) = trim(s_vector(i))
            end do
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
            do i=1,npulse_m
              polarization_m(i) = trim(s_vector(i))
            end do
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


  subroutine get_orbitals(string,input_orbitals,nOrb,nOrb2,nOrb2sc,nsOrb,npOrb,ndOrb,Orbs,sOrbs,pOrbs,dOrbs)
  !! Subroutine to read and parse the selected orbitals
    use mod_kind,    only: int32
    use mod_logging, only: log_message,log_error
    use mod_System,  only: System_type
    use mod_tools,   only: is_numeric,itos,stoi
    use AtomTypes,   only: default_Orbs
    use mod_superconductivity, only: superCond
    character(len=*)               , intent(in)    :: string
    character(len=20), dimension(:), intent(inout) :: input_orbitals
    integer(int32),                  intent(out)   :: nOrb,nOrb2,nOrb2sc,nsOrb,npOrb,ndOrb
    integer(int32),    dimension(:), allocatable, intent(out) :: Orbs,sOrbs,pOrbs,dOrbs

    integer :: i,iloc
    integer,  dimension(:), allocatable :: itmps_arr,itmpp_arr,itmpd_arr
    !! Variables to temporarily store the orbitals
    character(len=100) :: selected_orbitals,selected_sorbitals,selected_porbitals,selected_dorbitals
    !! Variables to print selected orbitals

    call log_message("get_orbitals","Selection of orbitals for " // trim(string) // ":")

    selected_orbitals = ""
    selected_sorbitals = ""
    selected_porbitals = ""
    selected_dorbitals = ""
    nOrb  = 0
    nsOrb = 0
    npOrb = 0
    ndOrb = 0

    do i = 1,size(input_orbitals)
      if(input_orbitals(i)(1:1) == "!") exit
      if(len_trim(input_orbitals(i)) == 0 .or. len_trim(input_orbitals(i)) == 4) cycle

      nOrb = nOrb + 1

      ! If the name of the orbital is given instead of a number, convert:
      if(.not.is_numeric( trim(input_orbitals(i)) )) then
        iloc = findloc( default_Orbs,input_orbitals(i)(1:3),dim=1 )
        if(iloc == 0) &
          call log_error("get_orbitals","Orbital not recognized: " // input_orbitals(i)(1:3) //". Use one of the following: " // NEW_line('A') // &
              "(1|s), (2|px), (3|py), (4|pz), (5|dxy), (6|dyz), (7|dzx), (8|dx2), (9|dz2)")
        input_orbitals(i) = itos( iloc )
      end if
    end do
    nOrb2 = nOrb*2
    nOrb2sc = superCond*nOrb2
    allocate(Orbs(nOrb),itmps_arr(nOrb),itmpp_arr(nOrb),itmpd_arr(nOrb))

    ! Looping over selected orbitals
    ! Transforms names to numbers
    ! and count s,p and d orbitals
    do i = 1,nOrb
      Orbs(i) = stoi( trim(input_orbitals(i)) )
      selected_orbitals = trim(selected_orbitals)  // " " // trim(default_Orbs(Orbs(i)))
      ! Checking orbital type, and storing information
      if(Orbs(i)==1) then
        nsOrb = nsOrb + 1
        itmps_arr(nsOrb) = i
        selected_sorbitals = trim(selected_sorbitals)  // " " // trim(default_Orbs(Orbs(i)))
      end if
      if((Orbs(i)>=2).and.(Orbs(i)<=4)) then
        npOrb = npOrb + 1
        itmpp_arr(npOrb) = i
        selected_porbitals = trim(selected_porbitals)  // " " // trim(default_Orbs(Orbs(i)))
      end if
      if((Orbs(i)>=5).and.(Orbs(i)<=9)) then
        ndOrb = ndOrb + 1
        itmpd_arr(ndOrb) = i
        selected_dorbitals = trim(selected_dorbitals)  // " " // trim(default_Orbs(Orbs(i)))
      end if
    end do

    call log_message("get_orbitals",trim(itos(nOrb)) // " orbitals selected:" // trim(selected_orbitals) // ", of which:")
    ! Storing position of each type of orbital
    ! s-orbitals
    if(nsOrb > 0) then
      allocate(sOrbs(nsOrb))
      sOrbs(1:nsOrb) = itmps_arr(1:nsOrb)
      call log_message("get_orbitals", trim(itos(nsOrb)) // " s orbitals:" // trim(selected_sorbitals) )
    end if
    ! p-orbitals
    if(npOrb > 0) then
      allocate(pOrbs(npOrb))
      pOrbs(1:npOrb) = itmpp_arr(1:npOrb)
      call log_message("get_orbitals", trim(itos(npOrb)) // " p orbitals:" // trim(selected_porbitals) )
    end if
    ! d-orbitals
    if(ndOrb > 0) then
      allocate(dOrbs(ndOrb))
      dOrbs(1:ndOrb) = itmpd_arr(1:ndOrb)
      call log_message("get_orbitals", trim(itos(ndOrb)) // " d orbitals:" // trim(selected_dorbitals) )
    end if

  end subroutine get_orbitals


  subroutine iowrite(s)
  !! Write read parameters on main output file
    use mod_mpi_pars,          only: numprocs
    use mod_parameters,        only: output,runoptions,bands,band_cnt,renorm,renormnb,eta,itype, &
                                     emin,emax,nener1,nqvec1
    use mod_magnet,            only: lfield,lconstraining_field,constr_type,hwx,hwy,hwz,hw_list,hw_count,dcfield,dcfield_dependence,total_hw_npt1
    use mod_System,            only: System_type,n0sc1,n0sc2
    use mod_BrillouinZone,     only: BZ => realBZ
    use mod_fermi_surface,     only: lfs_loop, fs_energy_i, fs_energy_f, fs_energy_npts
    use mod_SOC,               only: SOC, socscale
    use EnergyIntegration,     only: parts, parts3, n1gl, n3gl
    use ElectricField,         only: ElectricFieldMode, ElectricFieldVector, EFt, EFp
    use AdaptiveMesh,          only: minimumBZmesh
    use mod_superconductivity, only: lsupercond
    use mod_time_propagator,   only: integration_time, sc_tol, hE_0, hw1_m, hw_e, hw_m, tau_e, &
                                     tau_m, delay_e, delay_m, lelectric, lmagnetic, lpulse_e, npulse_e, lpulse_m, npulse_m, &
                                     polarization_vec_e, polarization_vec_m, abs_tol, rel_tol, safe_factor
#ifdef _GPU
    use mod_cuda,              only: num_gpus
#endif
    !$ use omp_lib
    implicit none
    type(System_type), intent(in) :: s
    integer :: i,j
#ifdef _OPENMP
    write(output%unit_loop,"(10x,'Running on ',i0,' MPI process(es) WITH ',i0,' openMP')") numprocs, omp_get_max_threads()
#else
    write(output%unit_loop,"(10x,'Running on ',i0,' MPI process(es) WITHOUT openMP')") numprocs
#endif
#ifdef _GPU
    write(output%unit_loop,"(10x,'      and ',i0,' GPUs per MPI process(es)')") num_gpus
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

    if(lconstraining_field) then
      write(output%unit_loop,"(1x,'Constranining field: ACTIVATED')")
      select case(constr_type)
      case(1)
        write(output%unit_loop,"(1x,'1 - Transverse constraining field (to fix direction of M)')")
      case(2)
        write(output%unit_loop,"(1x,'2 - Full constraining field (to keep M fixed)')")
      end select
    else
      write(output%unit_loop,"(1x,'Constranining field: DEACTIVATED')")
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
    if(trim(runoptions)/="") write(output%unit_loop,"(6x,'Activated options:',/,4x,a)") trim(runoptions)

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
