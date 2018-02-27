module mod_io
  ! Subroutines to read input file and variables containing these parameters
  implicit none
  logical :: log_unit = .false.
  character(len=12), parameter :: logfile = "parameter.in"

contains

  subroutine log_message(procedure, message)
    use mod_mpi_pars, only: myrank
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
    use mod_mpi_pars
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
       if(log_unit) then
          write(output%unit, "('[Error] [',a,'] ',a,'')") procedure, trim(message)
       else
          write(*, "('[Error] [',a,'] ',a,'')") procedure, trim(message)
       end if
    end if
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    stop
  end subroutine log_error

  subroutine log_warning(procedure, message)
    use mod_mpi_pars, only: myrank
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
       if(log_unit) then
          write(output%unit, "('[Warning] [',a,'] ',a,'')") procedure, trim(message)
       else
          write(*, "('[Warning] [',a,'] ',a,'')") procedure, trim(message)
       end if
    end if
  end subroutine log_warning

  subroutine get_parameters(filename, s)
    use mod_f90_kind, only: double
    use mod_mpi_pars
    use mod_input
    use mod_parameters, only: output, laddresults, lverbose, ldebug, lkpoints, &
                              lpositions, lcreatefiles, Utype, lnolb, lhfresponses, &
                              lnodiag, lsha, lcreatefolders, lwriteonscreen, runoptions, &
                              ltestcharge, llgtv, lsortfiles, magaxis, magaxisvec, &
                              itype, ry2ev, ltesla, eta, etap, dmax, emin, emax, deltae, &
                              skip_steps, npts, npt1, renorm, renormnb, bands, band_cnt, &
                              offset, dfttype, U, parField, parFreq, kptotal_in, kp_in
    use mod_self_consistency, only: lslatec, lontheflysc, lnojac, lGSL, lrotatemag, skipsc, scfile, mag_tol
    use mod_system, only: System, n0sc1, n0sc2
    use mod_SOC, only: SOC, socscale, llinearsoc, llineargfsoc
    use mod_magnet, only: lfield, tesla, hwa_i, hwa_f, hwa_npts, hwa_npt1, hwt_i, hwt_f, &
                          hwt_npts, hwt_npt1, hwp_i, hwp_f, hwp_npts, hwp_npt1, hwx, hwy, &
                          hwz, hwscale, hwtrotate, hwprotate, skip_steps_hw
    use TightBinding, only: tbmode, fermi_layer
    use ElectricField, only: ElectricFieldMode, ElectricFieldVector, EFp, EFt
    use EnergyIntegration, only: parts, parts3, pn1, pn2, pnt, n1gl, n3gl
    use mod_tools, only: itos
    use adaptiveMesh, only: minimumBZmesh
    implicit none
    character(len=*), intent(in) :: filename
    type(System), intent(inout) :: s
    character(len=20), allocatable :: s_vector(:)
    real(double), allocatable :: vector(:)
    integer, allocatable :: i_vector(:)
    integer :: i, cnt
    character(len=20) :: tmp_string
    if(.not. read_file(filename)) &
         call log_error("get_parameters", "File " // trim(filename) // " not found!")

    if(myrank == 0) then
       if(.not. enable_input_logging(logfile)) &
            call log_warning("get_parameters", "couldn't enable logging.")
    end if
    if(.not. get_parameter("output", output%file)) &
         call log_error("get_parameters", "Output filename not given!")

    if(myrank==0) open (unit=output%unit, file=trim(output%file), status='replace')
    log_unit = .true.

    if(myrank==0) &
         write(output%unit,"('[get_parameters] Reading parameters from ""',a,'"" file...')") trim(filename)

    !===============================================================================================
    !============= System configuration (Lattice + Reciprocal lattice) =============================
    !===============================================================================================
    if(.not. get_parameter("nn_stages", s%nStages,2)) call log_warning("get_parameters","'nn_stages' missing. Using default value of 2.")

    if(.not. get_parameter("bulk", s%lbulk, .true.)) call log_warning("get_parameters", "'bulk' missing. Using default value.")

    if(.not. get_parameter("nkpt", i_vector,cnt)) call log_error("get_parameters","'nkpt' missing.")
    if(cnt == 1) then
      kptotal_in = i_vector(1)
      if(s%lbulk) then
        kp_in(:) = ceiling((dble(kptotal_in))**(1.d0/3.d0))
        kptotal_in = kp_in(1) * kp_in(2) * kp_in(3)
      else
        kp_in(1:2) = ceiling((dble(kptotal_in))**(1.d0/2.d0))
        kp_in(3) = 0
        kptotal_in = kp_in(1) * kp_in(2)
      end if
    else if(cnt == 3) then
      kp_in(1) = i_vector(1)
      kp_in(2) = i_vector(2)
      kp_in(3) = i_vector(3)
      kptotal_in = kp_in(1) * kp_in(2) * kp_in(3)
    else
      call log_error("get_parameter", "'nkpt' has wrong size (expected 1 or 3)")
    end if

    if(.not. get_parameter("minimumBZmesh", minimumBZmesh, 1000)) call log_warning("get_parameters", "'minimumBZmesh' missing. Using default value.")
    !===============================================================================================
    !===============================================================================================

    !------------------------------------- Type of Calculation -------------------------------------
    if(.not. get_parameter("itype", itype)) call log_error("get_parameters","'itype' missing.")

    if(.not. get_parameter("Options", s_vector, cnt)) call log_error("get_parameters","'Options' missing.")
    runoptions = ""
    do i = 1, cnt
       select case (s_vector(i))
       case ("ry2ev")
          ry2ev = 13.6d0
       case ("tesla")
          tesla = 5.7883817555d-5/13.6d0 ! Ry/T
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
       case ("noUonall")
          if(Utype==1) then
             if(myrank==0) call log_warning("get_parameters","Runoption 'noUonNM' is already active.")
          else
             Utype = 0
          end if
       case ("noUonNM")
          if(Utype==0) then
             if(myrank==0) call log_warning("get_parameters","Runoption 'noUonall' is already active.")
          else
             Utype = 1
          end if
       case ("slatec")
          lslatec = .true.
       case ("GSL")
          lGSL = .true.
       case ("kpoints")
          lkpoints = .true.
       case ("positions")
          lpositions = .true.
       case ("lineargfsoc")
          llineargfsoc = .true.
       case ("linearsoc")
          llinearsoc = .true.
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
       case ("testcharge")
          ltestcharge = .true.
       case("!")
          exit
       case default
          call log_warning("get_parameters","Runoption " // trim(s_vector(i)) // " not found!")
          cycle
       end select
       runoptions  = trim(runoptions) // " " // trim(s_vector(i))
    end do


    !-------------------------------------- In-Plane Currents --------------------------------------
    if(.not. get_parameter("n0sc1", n0sc1)) call log_warning("get_parameters","'n0sc1' missing.")

    if(.not. get_parameter("n0sc2", n0sc2)) call log_warning("get_parameters","'n0sc2' missing.")


    !----------------------------------- Spin - Orbit - Coupling -----------------------------------
    if(.not. get_parameter("SOC", SOC)) &
        call log_warning("get_parameters","'SOC' missing. Using default value.")
    if(SOC) then
      if(.not. get_parameter("socscale", socscale, 1.d0)) &
          call log_warning("get_parameters","'socscale' missing. Using default value.")
      if(llinearsoc) then
        output%SOCchar = "L"
      else
        output%SOCchar = "T"
      end if
      if(abs(socscale-1.d0)>1.d-6) write(output%SOC,"('_socscale=',f5.2)") socscale
      if((llineargfsoc).or.(llinearsoc)) output%SOC = trim(output%SOC) // "_linearsoc"
    else
      output%SOCchar = "F"
    end if

    !---------------------------------------- Magnetization ----------------------------------------
    if(.not. get_parameter("magtol", mag_tol, 1.d-12)) call log_warning("get_parameters", "'magtol' not found. Using default value.")

    if(.not. get_parameter("magbasis", tmp_string)) then
        call log_warning("get_parameters","'magbasis' missing.")
      magaxis = 0
    else
      select case (tmp_string)
      case("cartesian")
        if(.not. get_parameter("magaxis", vector, cnt)) &
            call log_error("get_parameters","'magaxis' missing.")
        if(cnt /= 3) &
            call log_error("get_parameters","'magaxis' has wrong size (size 3 required).")
        magaxis = -1 ! TODO: Set it to a value if not determined otherwise?
        magaxisvec(1:3) = vector(1:3)
        deallocate(vector)
      case("neighbor")
        if(.not. get_parameter("magaxis", magaxis)) &
            call log_error("get_parameters","'magaxis' missing.")
      case("bravais")
        if(.not. get_parameter("magaxis", i_vector, cnt)) &
            call log_error("get_parameters","'magaxis' missing.")
        if(cnt /= 2) &
            call log_error("get_parameters","'magaxis' has wrong size (size 2 required).")
        magaxis = -2 ! TODO: Add options to evaluate these values.
        magaxisvec(1:2) = i_vector(1:2)
        deallocate(i_vector)
      case("spherical")
        if(.not. get_parameter("magaxis", vector, cnt)) &
            call log_error("get_parameters", "'magaxis' missing.")
        if(cnt /= 2) call log_error("get_parameters", "'magaxis' has wrong size (size 2 required).")
        magaxis = -3
        magaxisvec(1:2) = vector(1:2)
        deallocate(vector)
      end select
    end if

    !------------------------------------- Coulomb Interaction -------------------------------------
    if(.not. get_parameter("U", U, cnt)) then
      call log_warning("get_parameters","'U' not present. Using default value of 1eV for all sites.")
      allocate(U(1))
      U = 1.d0 / 13.6d0
    else if(cnt == 1) then
      call log_warning("get_parameters", "'U' has only single value. Using this on all sites.")
    else if(cnt < 0) then
      call log_error("get_parameters","'U' has size < 0.")
    end if

    !--------------------------------------- Electric Field ----------------------------------------
    if(.not. get_parameter("ebasis", tmp_string)) call log_error("get_parameters","'ebasis' missing.")
    select case (tmp_string)
    case("cartesian")
       if(.not. get_parameter("dirEfield", vector, cnt)) call log_error("get_parameters","'dirEfield' missing.")
       if(cnt /= 3) call log_error("get_parameters","'dirEfield' has wrong size (size 3 required).")
       ElectricFieldMode = -1 ! TODO: Set it to a value if not determined otherwise?
       ElectricFieldVector(1:3) = vector(1:3)
       deallocate(vector)
    case("neighbor")
       if(.not. get_parameter("dirEfield", ElectricFieldMode)) call log_error("get_parameters","'dirEfield' missing.")
    case("bravais")
       if(.not. get_parameter("dirEfield", i_vector, cnt)) call log_error("get_parameters","'dirEfield' missing.")
       if(cnt /= 2) call log_error("get_parameters","'dirEfield' has wrong size (size 2 required).")
       ElectricFieldMode = -2 ! TODO: Add options to evaluate these values.
       ElectricFieldVector(1:2) = i_vector(1:2)
       deallocate(i_vector)
    case("spherical")
       if(.not. get_parameter("dirEfield", vector, cnt)) call log_error("get_parameters", "'dirEfield' missing.")
       if(cnt /= 2) call log_error("get_parameters", "'dirEfield' has wrong size (size 2 required).")
       ElectricFieldMode = -3
       EFt = vector(1)
       EFp = vector(2)
       deallocate(vector)
    end select

    deallocate(s_vector)
    if(.not. get_parameter("eta", eta)) call log_error("get_parameters","'eta' missing.")
    if(.not. get_parameter("etap", etap, eta)) call log_warning("get_parameters", "'etap' not found. Using default value eta.")


    !------------------------------------- Static Magnetic Field -----------------------------------
    if(.not. get_parameter("FIELD", lfield, .false.)) call log_warning("get_parameters","'lfield' missing. Using default value.")
    if(lfield) then
       if(.not. get_parameter("hwa", vector, cnt)) call log_error("get_parameters","'hwa' missing.")
       if(cnt < 1) call log_error("get_parameters","'hwa' doesn't contain any parameters.")
       hwa_i = vector(1)
       if(cnt >= 2) hwa_f = vector(2)
       if(cnt >= 3) hwa_npts = vector(3)
       deallocate(vector)
       hwa_npt1 = hwa_npts + 1

       if(.not. get_parameter("hwt", vector, cnt)) call log_error("get_parameters","'hwt' missing.")
       if(cnt < 1) call log_error("get_parameters","'hwt' doesn't contain any parameters.")
       hwt_i = vector(1)
       if(cnt >= 2) hwt_f = vector(2)
       if(cnt >= 3) hwt_npts = vector(3)
       deallocate(vector)
       hwt_npt1 = hwt_npts + 1

       if(.not. get_parameter("hwp", vector, cnt)) call log_error("get_parameters","'hwp' missing.")
       if(cnt < 1) call log_error("get_parameters","'hwp' doesn't contain any parameters.")
       hwp_i = vector(1)
       if(cnt >= 2) hwp_f = vector(2)
       if(cnt >= 3) hwp_npts = vector(3)
       deallocate(vector)
       hwp_npt1 = hwp_npts + 1
       if(abs(hwa_i) < 1.d-9) then
         if(.not. get_parameter("hwx", hwx)) call log_error("get_parameters","'hwx' missing.")
         if(.not. get_parameter("hwy", hwy)) call log_error("get_parameters","'hwy' missing.")
         if(.not. get_parameter("hwz", hwz)) call log_error("get_parameters","'hwz' missing.")
       end if
    end if

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

    !------------------------------------ Integration Variables ------------------------------------


    if(.not. get_parameter("parts", parts)) call log_error("get_parameters","'parts' missing.")
    if(.not. get_parameter("parts3", parts3)) call log_error("get_parameters","'parts3' missing.")
    write(output%Energy, "('_parts=',i0,'_parts3=',i0)") parts,parts3


    if(.not. get_parameter("n1gl", n1gl)) call log_error("get_parameters","'n1gl' missing.")

    if(.not. get_parameter("n3gl", n3gl)) call log_error("get_parameters","'n3gl' missing.")

    if(.not. get_parameter("emin", emin)) call log_error("get_parameters","'emin' missing.")

    if(.not. get_parameter("emax", emax)) call log_error("get_parameters","'emax' missing.")

    if(.not. get_parameter("skip_steps", skip_steps, 0)) call log_warning("get_parameters","'skip_steps' missing.")

    if(.not. get_parameter("skip_steps_hw", skip_steps_hw, 0)) call log_warning("get_parameters","'skip_steps_hw' missing.")

    if(.not. get_parameter("npts", npts)) call log_error("get_parameters","'npts' missing.")
    npt1 = npts + 1

    if(.not. get_parameter("renorm", renorm)) call log_error("get_parameters","'renorm' missing.")
    if(renorm) then
       if(.not. get_parameter("renormnb", renormnb)) call log_error("get_parameters","'renormnb' missing.")
    end if



    !---------------------------------- Magnetic Self-consistency ----------------------------------
    if(.not. get_parameter("skipsc", skipsc, .false.)) call log_warning("get_parameters","'skipsc' missing. Using default value.")

    if(.not. get_parameter("scfile", scfile)) call log_warning("get_parameters","'scfile' missing.")



    !----------------------------------- Band Structure and LDOS -----------------------------------
    if(itype == 4) then
      if(.not. get_parameter("band", bands, band_cnt)) call log_error("get_parameters", "'band' missing.")
      if(band_cnt < 2) call log_error("get_parameters", "Need at least to Points for Band Structure")
    endif


    !======================================== Tight-Binding ========================================

    if(.not. get_parameter("tbmode", tbmode)) call log_error("get_parameters", "'tbmode' missing.")

    !---------------------------------------- Slater-Koster ----------------------------------------
    if(1 == tbmode) then
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

      if(.not. get_parameter("fermi_layer", fermi_layer, 1)) call log_warning("get_parameters", "'fermi_layer' not given. Using fermi_layer = 1")

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
    !==============================================================================================!

    if(.not. get_parameter("suffix", output%suffix)) call log_warning("get_parameters","'suffix' missing.")

    if(.not. get_parameter("parField", parField, 1)) call log_warning("get_parameters","'parField' missing. Using default value")
    if(.not. get_parameter("parFreq", parFreq, 1)) call log_warning("get_parameters","'parFreq' missing. Using default value")


    if(myrank == 0) then
      if(.not. disable_input_logging()) &
          call log_warning("get_parameters", "Could not disable logging.")
    end if


    if(myrank==0) &
        write(output%unit,"('[get_parameters] Finished reading from ""',a,'"" file')") trim(filename)


    !-------------------------------------------------------------------------------
    !*********** User manual additions / modifications in the input file **********!
    !     Npl_i  = 4
    !     Npl_f = 4
    !     nkpt = 6
    !     SOC = .true.
    !     magaxis = "5"
    !     runoptions = trim(runoptions) // " noUonNM"
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

    ! Energy loop step
    if(npts < 1) then
      deltae = 0
    else
      deltae = (emax - emin)/npts
    end if
    if(deltae<=1.d-14) npt1 = 1

    ! Preparing dc-limit calculation
    ! if(itype==9) call prepare_dclimit() !TODO: Re-Include


    pn1=parts*n1gl
    pn2=parts3*n3gl
    pnt=pn1+pn2
    return
  end subroutine get_parameters

  subroutine iowrite(s)
    use mod_mpi_pars
    use mod_parameters
    use mod_System, only: System, n0sc1, n0sc2
    use mod_BrillouinZone, only: BZ => realBZ
    use mod_SOC, only: SOC, socscale
    use mod_magnet
    use EnergyIntegration, only: parts, parts3, n1gl, n3gl
    use electricfield, only: ElectricFieldMode, ElectricFieldVector, EFt, EFp
    use adaptiveMesh, only: minimumBZmesh
    !$ use omp_lib
    implicit none
    type(System), intent(in) :: s
    integer :: i
#ifdef _OPENMP
    write(output%unit_loop,"(10x,'Running on ',i0,' MPI process(es) WITH ',i0,' openMP')") numprocs, omp_get_max_threads()
#else
    write(output%unit_loop,"(10x,'Running on ',i0,' MPI process(es) WITHOUT openMP')") numprocs
#endif
    write(output%unit_loop,"('|------------------------------- PARAMETERS: -------------------------------|')")
    write(output%unit_loop,"(10x,'nAtoms = ',i0)") s%nAtoms
    write(output%unit_loop,"(1x,'DFT parameters: ')", advance='no')
    dft_type: select case (dfttype)
    case ("T")
       write(output%unit_loop,"('Tight-binding basis')")
    case ("O")
       write(output%unit_loop,"('Orthogonal basis')")
    end select dft_type
    if(SOC) then
       write(output%unit_loop,"(1x,'Spin Orbit Coupling: ACTIVATED')")
       write(output%unit_loop,"(5x,'socscale =',es9.2)") socscale
    else
       write(output%unit_loop,"(1x,'Spin Orbit Coupling: DEACTIVATED')")
    end if
    write(output%unit_loop,"(8x,'Utype = ',i0)") Utype

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
       write(output%unit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (1)
       write(output%unit_loop,"(1x,'Self-consistency only')")
    case (2)
       write(output%unit_loop,"(1x,'Test after SC')")
       write(output%unit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
       write(output%unit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
       write(output%unit_loop,"(9x,'emin =',es9.2)") emin
       write(output%unit_loop,"(9x,'emax =',es9.2)") emax
       write(output%unit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (3)
       write(output%unit_loop,"(1x,'LDOS and exchange interactions as a function of energy')")
       write(output%unit_loop,"(9x,'emin =',es9.2)") emin
       write(output%unit_loop,"(9x,'emax =',es9.2)") emax
       write(output%unit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (4)
       write(output%unit_loop,"(1x,'Band structure')")
       write(output%unit_loop,"(9x,'Number of points to calculate: ',i0)") npt1
    case (5)
       write(output%unit_loop,"(1x,'Charge and spin density at Fermi surface')")
    case (6)
       write(output%unit_loop,"(1x,'Exhange interactions and anisotropies (full tensor)')")
       if(nmaglayers==1) write(output%unit_loop,"(1x,'Only 1 magnetic layer: calculating only anisotropies')")
       !write(outputunit_loop,"(8x,'from Npl = ',i0,' to ',i0)") Npl_i,Npl_f
    case (7)
       write(output%unit_loop,"(1x,'Local susceptibility as a function of energy')")
       write(output%unit_loop,"(9x,'emin =',es9.2)") emin
       write(output%unit_loop,"(9x,'emax =',es9.2)") emax
       !write(output%unit_loop,"(1x,i0,' points divided into ',i0,' steps of size',es10.3,' each calculating ',i0,' points')") npt1,MPIsteps,MPIdelta,MPIpts
    case (8)
       write(output%unit_loop,"(1x,'Parallel currents, disturbances and local susc. as a function of energy')")
       write(output%unit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
       write(output%unit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
       write(output%unit_loop,"(9x,'emin =',es9.2)") emin
       write(output%unit_loop,"(9x,'emax =',es9.2)") emax
       !write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps of energy size',es10.3,' each calculating ',i0,' points')") total_hw_npt1*npt1,MPIsteps*MPIsteps_hw,MPIdelta,MPIpts_hw*MPIpts
    case (9)
       write(output%unit_loop,"(1x,'dc limit calculations as a function of ',a)") trim(dcfield(dcfield_dependence))
       write(output%unit_loop,"(1x,'e =',es9.2)") emin
       write(output%unit_loop,"(1x,'hwa_min =',es9.2)") hw_list(1,1)
       write(output%unit_loop,"(1x,'hwa_max =',es9.2)") hw_list(total_hw_npt1,1)
       write(output%unit_loop,"(1x,'hwt_min =',f7.2)") hw_list(1,2)
       write(output%unit_loop,"(1x,'hwt_max =',f7.2)") hw_list(total_hw_npt1,2)
       write(output%unit_loop,"(1x,'hwp_min =',f7.2)") hw_list(1,3)
       write(output%unit_loop,"(1x,'hwp_max =',f7.2)") hw_list(total_hw_npt1,3)
       !write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps, each calculating ',i0,' points')") total_hw_npt1*npt1,MPIsteps*MPIsteps_hw,MPIpts_hw*MPIpts
    end select write_itype
    write(output%unit_loop,"('|---------------------------------------------------------------------------|')")
    return
  end subroutine iowrite

end module mod_io
