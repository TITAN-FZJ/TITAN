module mod_io
  use mod_f90_kind, only: double
  use mod_parameters
  implicit none
#ifdef _JUQUEEN
  character(len=200) :: filename="input_juqueen"
#else
  character(len=200) :: filename="input"
#endif
  logical :: log_unit = .false.

contains

  subroutine log_message(procedure, message)
    use mod_mpi_pars, only: myrank
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
       if(log_unit) then
          write(outputunit, "('[',a,'] ',a,'')") procedure, trim(message)
       else
          write(*, "('[',a,'] ',a,'')") procedure, trim(message)
       end if
    end if
  end subroutine log_message

  subroutine log_error(procedure, message)
    use mod_mpi_pars
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
       if(log_unit) then
          write(outputunit, "('[Error] [',a,'] ',a,'')") procedure, trim(message)
       else
          write(*, "('[Error] [',a,'] ',a,'')") procedure, trim(message)
       end if
    end if
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    stop
  end subroutine log_error

  subroutine log_warning(procedure, message)
    use mod_mpi_pars, only: myrank
    implicit none
    character(len=*), intent(in) :: procedure
    character(len=*), intent(in) :: message

    if(myrank == 0) then
       if(log_unit) then
          write(outputunit, "('[Warning] [',a,'] ',a,'')") procedure, trim(message)
       else
          write(*, "('[Warning] [',a,'] ',a,'')") procedure, trim(message)
       end if
    end if
  end subroutine log_warning

  subroutine get_parameters()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_input
    use mod_system, only: a0, a1, a2, a3, pln_normal, nkpt, nstages, n0sc1, n0sc2
    use mod_tight_binding, only: layers, tbmode, fermi_layer
    implicit none
    character(len=20), allocatable :: s_vector(:)
    real(double), allocatable :: vector(:)
    integer, allocatable :: i_vector(:)
    integer :: i, cnt
    character(len=20) :: tmp_string

    if(.not. read_file(filename)) &
         call log_error("get_parameters", "File " // trim(filename) // " not found!")

    if(myrank == 0) then
       if(.not. enable_input_logging("parameter.in")) &
            call log_warning("get_parameters", "couldn't enable logging.")
    end if
    if(.not. get_parameter("output", outputfile)) &
         call log_error("get_parameters", "Output filename not given!")

    if(myrank==0) open (unit=outputunit, file=trim(outputfile), status='replace')
    log_unit = .true.

    if(myrank==0) &
         write(outputunit,"('[get_parameters] Reading parameters from ""',a,'"" file...')") trim(filename)

    !===============================================================================================
    !============= System configuration (Lattice + Reciprocal lattice) =============================
    !===============================================================================================

    if(.not. get_parameter("a0", a0)) call log_error("get_parameters","'a0' missing.")

    if(.not. get_parameter("lattice", lattice)) call log_error("get_parameters","'lattice' missing.")

    select case(lattice)
    case("general")
       if(.not. get_parameter("a1", vector, cnt)) &
            call log_error("get_parameters"," 'a1' missing.")
       if(cnt /= 3) call log_error("get_parameters","'a1' has wrong size (size 3 required).")
       a1 = vector(1:3) * a0 / sqrt(dot_product(vector(1:3), vector(1:3)))
       deallocate(vector)

       if(.not. get_parameter("a2", vector, cnt)) call log_error("get_parameters","'a2' missing.")
       if(cnt /= 3) call log_error("get_parameters","'a3' has wrong size (size 3 required).")
       a2 = vector
       deallocate(vector)

       if(.not. get_parameter("a3", vector, cnt)) call log_error("get_parameters","'a3' missing.")
       if(cnt /= 3) call log_error("get_parameters","'a3' has wrong size (size 3 required).")
       a3 = vector
       deallocate(vector)
    case("bcc")
       a1 = [-0.5d0,  0.5d0,  0.5d0] * a0 !/ sqrt(0.75d0)
       a2 = [ 0.5d0, -0.5d0,  0.5d0] * a0 !/ sqrt(0.75d0)
       a3 = [ 0.5d0,  0.5d0, -0.5d0] * a0 !/ sqrt(0.75d0)
    case("fcc")
       a1 = [0.0d0, 0.5d0, 0.5d0] * a0 !/ sqrt(0.5d0)
       a2 = [0.5d0, 0.0d0, 0.5d0] * a0 !/ sqrt(0.5d0)
       a3 = [0.5d0, 0.5d0, 0.0d0] * a0 !/ sqrt(0.5d0)
    case("sc")
       a1 = [1.0d0, 0.0d0, 0.0d0] * a0
       a2 = [0.0d0, 1.0d0, 0.0d0] * a0
       a3 = [0.0d0, 0.0d0, 1.0d0] * a0
    case("hcp")
       a1 = [1.0d0, 0.0d0, 0.0d0] * a0
       a2 = [0.5d0, sqrt(3.0d0)*0.5d0, 0.0d0] * a0
       a3 = [0.0d0, 0.0d0, sqrt(8.0d0 / 3.0d0)] * a0 / sqrt( 8.0d0 / 3.0d0 )
    case default
       call log_error("get_parameters","Unknown 'lattice' option.")
    end select

    if(get_parameter("plane", vector, cnt)) then
      if(cnt /= 3) call log_error("get_parameters","'plane' has wrong size (size 3 required).")
      pln_normal = vector(1:3) / sqrt(dot_product(vector(1:3),vector(1:3)))
      deallocate(vector)
    else
      pln_normal = 0.d0
      call log_warning("get_parameters","'plane' missing.")
    end if

    if(.not. get_parameter("nkpt", nkpt)) call log_error("get_parameters","'nkpt' missing.")

    if(.not. get_parameter("nn_stages", nstages)) call log_warning("get_parameters","'nn_stages' missing.")

    if(.not. get_parameter("n0sc1", n0sc1)) call log_warning("get_parameters","'n0sc1' missing.")

    if(.not. get_parameter("n0sc2", n0sc2)) call log_warning("get_parameters","'n0sc2' missing.")

    !===============================================================================================
    !===============================================================================================

    if(.not. get_parameter("itype", itype)) call log_error("get_parameters","'itype' missing.")

    if(.not. get_parameter("SOC", SOC)) call log_error("get_parameters","'SOC' missing.")
    if(SOC) then
       if(.not. get_parameter("socscale", socscale)) call log_warning("get_parameters","'socscale' missing.")
    end if

    if(.not. get_parameter("magbasis", tmp_string)) then
       call log_warning("get_parameters","'magbasis' missing.")
       magaxis = 0
    else
       select case (tmp_string)
       case("cartesian")
          if(.not. get_parameter("magaxis", vector, cnt)) call log_error("get_parameters","'magaxis' missing.")
          if(cnt /= 3) call log_error("get_parameters","'magaxis' has wrong size (size 3 required).")
          magaxis = -1 ! TODO: Set it to a value if not determined otherwise?
          magaxisvec(1:3) = vector(1:3)
          deallocate(vector)
       case("neighbor")
          if(.not. get_parameter("magaxis", magaxis)) call log_error("get_parameters","'magaxis' missing.")
       case("bravais")
          if(.not. get_parameter("magaxis", i_vector, cnt)) call log_error("get_parameters","'magaxis' missing.")
          if(cnt /= 2) call log_error("get_parameters","'magaxis' has wrong size (size 2 required).")
          magaxis = -2 ! TODO: Add options to evaluate these values.
          magaxisvec(1:2) = i_vector(1:2)
          deallocate(i_vector)
        case("spherical")
           if(.not. get_parameter("magaxis", vector, cnt)) call log_error("get_parameters", "'magaxis' missing.")
           if(cnt /= 2) call log_error("get_parameters", "'magaxis' has wrong size (size 2 required).")
           magaxis = -3
           magaxisvec(1:2) = vector(1:2)
           deallocate(vector)
        end select
    end if

    if(.not. get_parameter("ebasis", tmp_string)) call log_error("get_parameters","'ebasis' missing.")
    select case (tmp_string)
    case("cartesian")
       if(.not. get_parameter("dirEfield", vector, cnt)) call log_error("get_parameters","'dirEfield' missing.")
       if(cnt /= 3) call log_error("get_parameters","'dirEfield' has wrong size (size 3 required).")
       dirEfield = -1 ! TODO: Set it to a value if not determined otherwise?
       dirEfieldvec(1:3) = vector(1:3)
       deallocate(vector)
    case("neighbor")
       if(.not. get_parameter("dirEfield", dirEfield)) call log_error("get_parameters","'dirEfield' missing.")
    case("bravais")
       if(.not. get_parameter("dirEfield", i_vector, cnt)) call log_error("get_parameters","'dirEfield' missing.")
       if(cnt /= 2) call log_error("get_parameters","'dirEfield' has wrong size (size 2 required).")
       dirEfield = -2 ! TODO: Add options to evaluate these values.
       dirEfieldvec(1:2) = i_vector(1:2)
       deallocate(i_vector)
    case("spherical")
       if(.not. get_parameter("dirEfield", vector, cnt)) call log_error("get_parameters", "'dirEfield' missing.")
       if(cnt /= 2) call log_error("get_parameters", "'dirEfield' has wrong size (size 2 required).")
       dirEfield = -3
       EFt = vector(1)
       EFp = vector(2)
       deallocate(vector)
    end select

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
    deallocate(s_vector)
    if(.not. get_parameter("eta", eta)) call log_error("get_parameters","'eta' missing.")

    if(.not. get_parameter("FIELD", lfield)) call log_error("get_parameters","'lfield' missing.")
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

       if(.not. get_parameter("hwx", hwx)) call log_error("get_parameters","'hwx' missing.")
       if(.not. get_parameter("hwy", hwy)) call log_error("get_parameters","'hwy' missing.")
       if(.not. get_parameter("hwz", hwz)) call log_error("get_parameters","'hwz' missing.")
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

    if(.not. get_parameter("parts", parts)) call log_error("get_parameters","'parts' missing.")

    if(.not. get_parameter("parts3", parts3)) call log_error("get_parameters","'parts3' missing.")

    if(.not. get_parameter("n1gl", n1gl)) call log_error("get_parameters","'n1gl' missing.")

    if(.not. get_parameter("n3gl", n3gl)) call log_error("get_parameters","'n3gl' missing.")

    if(.not. get_parameter("emin", emin)) call log_error("get_parameters","'emin' missing.")

    if(.not. get_parameter("emax", emax)) call log_error("get_parameters","'emax' missing.")

    if(.not. get_parameter("skip_steps", skip_steps)) call log_warning("get_parameters","'skip_steps' missing.")

    if(.not. get_parameter("skip_steps_hw", skip_steps_hw)) call log_warning("get_parameters","'skip_steps_hw' missing.")

    if(.not. get_parameter("npts", npts)) call log_error("get_parameters","'npts' missing.")
    npt1 = npts + 1

    if(.not. get_parameter("renorm", renorm)) call log_error("get_parameters","'renorm' missing.")
    if(renorm) then
       if(.not. get_parameter("renormnb", renormnb)) call log_error("get_parameters","'renormnb' missing.")
    end if

    if(.not. get_parameter("skipsc", skipsc)) call log_error("get_parameters","'skipsc' missing.")

    if(.not. get_parameter("scfile", scfile)) call log_warning("get_parameters","'scfile' missing.")

    if(itype == 4) then
      if(.not. get_parameter("band", bands, band_cnt)) call log_error("get_parameters", "'bands' missing.")
      if(band_cnt < 2) call log_error("get_parameters", "Need at least to Points for Band Structure")
    endif

    !==============================================================================================!
    !================================Tight-Binding Mode============================================!
    if(.not. get_parameter("tbmode", tbmode)) call log_error("get_parameters", "'tbmode' missing.")



    if(1 == tbmode) then
       !-----Slater-Koster--------------------------------------------------------------------------!
       dfttype = "S"

       if(.not. get_parameter("layers", layers, cnt)) call log_error("get_parameters", "'layers' missing.")
       if(cnt <= 0) call log_error("get_parameters", "'layers' No layers given.")
       Npl = cnt

       if(get_parameter("Npl", i_vector, cnt)) then
         if(cnt < 1) call log_error("get_parameters","'Npl' doesn't contain any parameters.")
         Npl_i = i_vector(1)
         Npl_f = i_vector(1)
         if(cnt >= 2) Npl_f = i_vector(2)
         if(Npl_f < Npl_i) Npl_f = Npl_i
         deallocate(i_vector)
         if(Npl < Npl_f) call log_error("get_parameters", "'Npl' larger than amount of given layers")
       else
         call log_warning("get_parameters","'Npl' missing.")
         Npl_i = Npl
         Npl_f = Npl
       end if
       if(.not. get_parameter("fermi_layer", fermi_layer)) call log_warning("get_parameters", "'fermi_layer' not given. Using fermi_layer = 1")
       if(fermi_layer <= 0 .or. fermi_layer > Npl) call log_error("get_parameters", "'fermi_layer' out of range.")

    else if(2 == tbmode) then
       !-----DFT------------------------------------------------------------------------------------!
       offset = 1
       if(nstages /= 2) call log_error("get_parameters", "'tbmode' DFT Mode only supports nstages = 2")

       if(.not. get_parameter("dfttype", dfttype)) call log_error("get_parameters","'dfttype' missing.")

       if(.not. get_parameter("Npl", i_vector, cnt)) call log_error("get_parameters","'Npl' missing.")
       if(cnt < 1) call log_error("get_parameters","'Npl' doesn't contain any parameters.")
       Npl_i = i_vector(1)
       Npl_f = i_vector(1)
       if(cnt >= 2) Npl_f = i_vector(2)
       deallocate(i_vector)

       ! Check number of planes
       if(Npl_f < Npl_i) Npl_f = Npl_i

       if(.not. get_parameter("set1", set1)) call log_error("get_parameters","'set1' missing.")

       if(.not. get_parameter("set2", set2)) call log_error("get_parameters","'set2' missing.")

       if(get_parameter("addlayers", i_vector, cnt)) then
          if(cnt < 10) then
             addlayers(1:cnt) = i_vector(1:cnt)
             naddlayers = cnt
          else if(cnt >= 10) then
             addlayers(1:10) = i_vector(1:10)
             naddlayers = 10
          end if
       end if
       if(allocated(i_vector)) deallocate(i_vector)

       ! Add 'naddlayers' to Npl
       if((naddlayers==1).and.(myrank==0)) write(outputunit,"('[get_parameters] WARNING: Added layers must include empty spheres! Only including one layer: naddlayers = ',i0)") naddlayers
       if((set1==9).or.(set2==9)) then
          naddlayers = 0
       end if
       if(naddlayers/=0) then
          Npl_i = Npl_i+naddlayers-1
          Npl_f = Npl_f+naddlayers-1
       end if


    else
       call log_error("get_parameters", "'tbmode' Unknown mode selected. (Choose either 1 or 2)0")
    end if
    !==============================================================================================!
    !==============================================================================================!

    if(.not. get_parameter("suffix", suffix)) call log_warning("get_parameters","'suffix' missing.")


    if(myrank == 0) then
       if(.not. disable_input_logging()) call log_warning("get_parameters", "couldn't disable logging.")
    end if


    if(myrank==0) write(outputunit,"('[get_parameters] Finished reading from ""',a,'"" file')") trim(filename)

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
       if(myrank==0) then
          write(outputunit,"('[get_parameters] Invalid neighbor for renormalization: ',i0,'!')") renormnb
          write(outputunit,"('[get_parameters] Choose a value between ',i0,' and ',i0,'.')") n0sc1,n0sc2
       end if
       call MPI_Finalize(ierr)
       stop
    end if
    if(skip_steps<0) then
       if(myrank==0) write(outputunit,"('[get_parameters] Invalid number of energy steps to skip: ',i0)") skip_steps
       call MPI_Finalize(ierr)
       stop
    end if
    if(skip_steps_hw<0) then
       if(myrank==0) write(outputunit,"('[get_parameters] Invalid number of field steps to skip: ',i0)") skip_steps_hw
       call MPI_Finalize(ierr)
       stop
    end if
    if((lhfresponses).and.(itype==7).and.(myrank==0)) write(outputunit,"('[get_parameters] Susceptibility calculations already include HF responses. Ignoring ""hfresponses"" runoption')")
    ! Adjusting zeeman energy to Ry or eV
    tesla = tesla*ry2ev

    ! Turning off renormalization for non-current calculations
    if(itype/=8) renorm = .false.

    ! Setting up external field variables and loops
    call prepare_field()

    ! Energy loop step
    deltae = (emax - emin)/npts
    if(deltae<=1.d-14) npt1 = 1

    ! Preparing dc-limit calculation
    if(itype==9) call prepare_dclimit()


    tol   = 1.d-10
    pn1=parts*n1gl
    pn2=parts3*n3gl
    pnt=pn1+pn2
    return
  end subroutine get_parameters

  subroutine iowrite()
    use mod_mpi_pars
    use mod_system, only: nkpt, n0sc1, n0sc2
    implicit none

#ifdef _OPENMP
    write(outputunit_loop,"(10x,'Running on ',i0,' MPI process(es) WITH openMP')") numprocs
#else
    write(outputunit_loop,"(10x,'Running on ',i0,' MPI process(es) WITHOUT openMP')") numprocs
#endif
    write(outputunit_loop,"('|------------------------------- PARAMETERS: -------------------------------|')")
    write(outputunit_loop,"(10x,'Npl = ',i0)") Npl
    write(outputunit_loop,"(1x,'DFT parameters: ')", advance='no')
    dft_type: select case (dfttype)
    case ("T")
       write(outputunit_loop,"('Tight-binding basis')")
    case ("O")
       write(outputunit_loop,"('Orthogonal basis')")
    end select dft_type
    if(SOC) then
       write(outputunit_loop,"(1x,'Spin Orbit Coupling: ACTIVATED')")
       write(outputunit_loop,"(5x,'socscale =',es9.2)") socscale
    else
       write(outputunit_loop,"(1x,'Spin Orbit Coupling: DEACTIVATED')")
    end if
    write(outputunit_loop,"(8x,'Utype = ',i0)") Utype

    write(outputunit_loop,"(1x,'Electric field direction: ')", advance='no')
    select case(dirEfield)
    case(-3)
       write(outputunit_loop,"('Spherical theta=',es8.1,' phi=',es8.1)") EFt, EFp
    case(-2)
       write(outputunit_loop,"('Bravais ')")
    case(-1)
       write(outputunit_loop,"('Cartesian')")
    case(1:99)
       write(outputunit_loop, "('Neighbor ',i0)") dirEfield
    end select
    write(outputunit_loop,"('Direction ',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)

    if(renorm) then
       write(outputunit_loop,"(1x,'Current renormalization: ACTIVATED')")
       write(outputunit_loop,"(5x,'renormnb = ',i0)") renormnb
    else
       write(outputunit_loop,"(1x,'Current renormalization: DEACTIVATED')")
    end if
    write(outputunit_loop,"(9x,'nkpt = ',i0)") nkpt
    write(outputunit_loop,"(8x,'parts = ',i0,'x',i0)") parts,n1gl
    write(outputunit_loop,"(7x,'parts3 = ',i0,'x',i0)") parts3,n3gl
    write(outputunit_loop,"(10x,'eta =',es9.2)") eta
    if(lfield) then
       write(outputunit_loop,"(1x,'Static magnetic field: ACTIVATED')")
       write(outputunit_loop,"(10x,'hwx =',es9.2,5x,'|',5x,'hwa =',es9.2)") hwx,hw_list(hw_count,1)
       write(outputunit_loop,"(10x,'hwy =',es9.2,5x,'|',5x,'hwt =',f6.3)") hwy,hw_list(hw_count,2)
       write(outputunit_loop,"(10x,'hwz =',es9.2,5x,'|',5x,'hwp =',f6.3)") hwz,hw_list(hw_count,3)
    else
       write(outputunit_loop,"(1x,'Static magnetic field: DEACTIVATED')")
    end if
    if(runoptions/="") write(outputunit_loop,"(6x,'Activated options:',/,4x,a)") trim(runoptions)

    write(outputunit_loop,"('|------------------------------ TO CALCULATE: ------------------------------|')")
    write_itype: select case (itype)
    case (0)
       write(outputunit_loop,"(1x,'Test before SC')")
       write(outputunit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
       write(outputunit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
       write(outputunit_loop,"(9x,'emin =',es9.2)") emin
       write(outputunit_loop,"(9x,'emax =',es9.2)") emax
       write(outputunit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (1)
       write(outputunit_loop,"(1x,'Self-consistency only')")
    case (2)
       write(outputunit_loop,"(1x,'Test after SC')")
       write(outputunit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
       write(outputunit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
       write(outputunit_loop,"(9x,'emin =',es9.2)") emin
       write(outputunit_loop,"(9x,'emax =',es9.2)") emax
       write(outputunit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (3)
       write(outputunit_loop,"(1x,'LDOS and exchange interactions as a function of energy')")
       write(outputunit_loop,"(9x,'emin =',es9.2)") emin
       write(outputunit_loop,"(9x,'emax =',es9.2)") emax
       write(outputunit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (4)
       write(outputunit_loop,"(1x,'Band structure')")
       write(outputunit_loop,"(9x,'Number of points to calculate: ',i0)") npt1
    case (5)
       write(outputunit_loop,"(1x,'Charge and spin density at Fermi surface')")
    case (6)
       write(outputunit_loop,"(1x,'Exhange interactions and anisotropies (full tensor)')")
       if(nmaglayers==1) write(outputunit_loop,"(1x,'Only 1 magnetic layer: calculating only anisotropies')")
       write(outputunit_loop,"(8x,'from Npl = ',i0,' to ',i0)") Npl_i,Npl_f
    case (7)
       write(outputunit_loop,"(1x,'Local susceptibility as a function of energy')")
       write(outputunit_loop,"(9x,'emin =',es9.2)") emin
       write(outputunit_loop,"(9x,'emax =',es9.2)") emax
       write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps of size',es10.3,' each calculating ',i0,' points')") npt1,MPIsteps,MPIdelta,MPIpts
    case (8)
       write(outputunit_loop,"(1x,'Parallel currents, disturbances and local susc. as a function of energy')")
       write(outputunit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
       write(outputunit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
       write(outputunit_loop,"(9x,'emin =',es9.2)") emin
       write(outputunit_loop,"(9x,'emax =',es9.2)") emax
       write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps of energy size',es10.3,' each calculating ',i0,' points')") total_hw_npt1*npt1,MPIsteps*MPIsteps_hw,MPIdelta,MPIpts_hw*MPIpts
    case (9)
       write(outputunit_loop,"(1x,'dc limit calculations as a function of ',a)") trim(dcfield(dcfield_dependence))
       write(outputunit_loop,"(1x,'e =',es9.2)") emin
       write(outputunit_loop,"(1x,'hwa_min =',es9.2)") hw_list(1,1)
       write(outputunit_loop,"(1x,'hwa_max =',es9.2)") hw_list(total_hw_npt1,1)
       write(outputunit_loop,"(1x,'hwt_min =',f6.3)") hw_list(1,2)
       write(outputunit_loop,"(1x,'hwt_max =',f6.3)") hw_list(total_hw_npt1,2)
       write(outputunit_loop,"(1x,'hwp_min =',f6.3)") hw_list(1,3)
       write(outputunit_loop,"(1x,'hwp_max =',f6.3)") hw_list(total_hw_npt1,3)
       write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps, each calculating ',i0,' points')") total_hw_npt1*npt1,MPIsteps*MPIsteps_hw,MPIpts_hw*MPIpts
    end select write_itype
    write(outputunit_loop,"('|---------------------------------------------------------------------------|')")
    return
  end subroutine iowrite

end module mod_io
