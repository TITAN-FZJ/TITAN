!> Main program:
program TITAN
  !! Initialize MPI, reads input, define lattice quantities then
  !! loops over number of planes and magnetic field to do the
  !! magnetic self-consistency and calculate either ground state
  !! quantities or response functions
  !use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_io
  use mod_system
  use mod_polyBasis, only: read_basis
  use Lattice
  use mod_BrillouinZone
  use TightBinding
  use mod_SOC
  use mod_magnet
  use ElectricField
  use adaptiveMesh
  use AtomTypes
  use mod_TCM
  use mod_self_consistency
  use EnergyIntegration
  use mod_progress
  use mod_mpi_pars
  use mod_Umatrix
  use mod_check_stop
  use mod_Atom_variables, only: allocate_Atom_variables, deallocate_Atom_variables
  use mod_tools, only: rtos
  use mod_initial_expectation, only: calc_initial_Uterms
  use mod_superconductivity
  use mod_time_propagator,   only: time_propagator
  !use mod_define_system TODO: Re-include
  !use mod_prefactors TODO: Re-include
  !use mod_lgtv_currents TODO: Re-include
  !use mod_sha TODO: Re-include
  !use mod_torques, only: ntypetorque
  implicit none

  character(len=500) :: arg = ""

  !------------------------ MPI initialization -------------------------
  call Initialize_MPI()

  !------------------------- Starting program --------------------------
  start_program = MPI_Wtime()

  !------------------------ Reading parameters -------------------------
  if(command_argument_count() >= 1) call get_command_argument(1, arg)
  if(len_trim(arg) > 0) then
    call get_parameters(arg, sys)
  else
    call get_parameters("input",sys)
  end if

  if(myrank == 0) call write_time(output%unit,'[main] Started on: ')

  !------------------- Useful constants and matrices -------------------
  call define_constants() ! TODO: Review

  !------------------- Define the lattice structure --------------------
  call read_basis("basis", sys)
  call initLattice(sys)
  write(output%Sites,fmt="(i0,'Sites')") sys%nAtoms
  ! Writing Positions into file
  if( lpositions .and. (myrank==0) ) call writeLattice(sys)

  !------------- Creating folders for current calculation ------------
  if(lcreatefolders) &
    call create_folder()

  !------------------ Set Loops and Integration Points -----------------
  call setMagneticLoopPoints()
  call setLoops(sys)

  !------------------ Creating grid of MPI processes  ------------------
  !call setup_MPI_grid(itype, pn1, nEner1, pnt,total_hw_npt1, nEner, deltae, emin, emax)
  call genMPIGrid(parField, total_hw_npt1, parFreq, nEner1 - skip_steps)

  !--- Generating integration points of the complex energy integral ----
  call allocate_energy_points()
  call generate_imag_epoints()

  !-------------------------- Debugging part -------------------------
  if(ldebug) then
    ! call initTightBinding(sys)
    ! call allocate_magnet_variables(sys%nAtoms, nOrb)
    ! call hamilt_sc(sys)
    !call hamilt_sc(s,rho,mp,mx,my,mz,sys,kp,hk)
    !if(myrank==0) then
    ! call debugging()
    !end if
    !MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    ! call deallocate_magnet_variables()
    ! stop
  end if

  !----- Calculating initial values in the hamiltonian with mag=0 ------
  if(.not.lsortfiles) call calc_initial_Uterms(sys)

  !----------- Generating k points for real axis integration -----------
  realBZ % nkpt_x = kp_in(1)
  realBZ % nkpt_y = kp_in(2)
  realBZ % nkpt_z = kp_in(3)
  call realBZ % count(sys)

  !---------------- Reading Tight Binding parameters -------------------
  call initTightBinding(sys)

  !---- Generating k meshes points for imaginary axis integration ------
  call generateAdaptiveMeshes(sys,pn1)

  !------------ Allocating variables that depend on nAtoms -------------
  call allocate_magnet_variables(sys%nAtoms, nOrb)
  call allocLS(sys%nAtoms,nOrb)
  call allocate_Atom_variables(sys%nAtoms,nOrb)

  !---------------------------- Dimensions -----------------------------
  dimH = sys%nAtoms*nOrb*2
  dimspinAtoms = 4 * sys%nAtoms
  dim = dimspinAtoms * nOrb * nOrb

  !---------- Conversion arrays for dynamical quantities ---------------
  call initConversionMatrices(sys%nAtoms,nOrb)

  !------- Initialize Stride Matrices for hamiltk and dtdksub ----------
  call initHamiltkStride(sys%nAtoms, nOrb)

  !--------------------- Lattice specific variables --------------------
  call initElectricField(sys%a1, sys%a2, sys%a3)

  !call setup_long_and_trans_current_neighbors(sys) !TODO: Not implemented TODO: Re-include

  !--------------------------- Build U array ---------------------------
  call build_U(sys)

  !-------------------------- Filename strings -------------------------
  write(output%info,"('_nkpt=',i0,'_eta=',a)") kptotal_in, trim(rtos(eta,"(es8.1)"))
  if(leigenstates) output%info = trim(output%info) // "_ev"

  !------------------------ MAGNETIC FIELD LOOP ------------------------
  if(myrank == 0 .and. skip_steps_hw > 0) &
  write(output%unit,"('[main] Skipping first ',i0,' field step(s)...')") skip_steps_hw

  do hw_count = startField, endField
    !---------- Opening files (general and one for each field) ---------
    if(rField == 0) then
      if(total_hw_npt1 == 1 .or. itype == 6 ) then
        output%file_loop = output%file
        output%unit_loop = output%unit
      else
        write(output%file_loop, fmt="(a,a,i0)") trim(output%file), '.', hw_count
        output%unit_loop = output%unit + hw_count
        open (unit=output%unit_loop, file=trim(output%file_loop), status='replace')
      end if
    end if

    !------------------- Initialize Magnetic Field ---------------------
    call initMagneticField(sys%nAtoms)
    if((llinearsoc) .or. (.not.SOC) .and. (.not.lfield) .and. (rField == 0)) &
      write(output%file_loop,"('[main] WARNING: No external magnetic field is applied and SOC is off/linear order: Goldstone mode is present!')")

    !----------------------- Update fieldpart --------------------------
    call set_fieldpart(hw_count)

    !----------------- Tests for coupling calculation ------------------
    if(itype == 6 .and. myrank == 0) then
      if(lfield) &
      call abortProgram("[main] Coupling calculation is valid for no external field!")
    end if

    !----- Writing parameters and data to be calculated on screen ------
    if(rField == 0) call iowrite(sys)

    !---- L matrix in global frame for given quantization direction ----
    call l_matrix()

    !---- Calculate L.S matrix for the given quantization direction ----
    if(SOC) call updateLS(sys,theta, phi)

    !---- Calculate L.B matrix for the given quantization direction ----
    call lb_matrix(sys%nAtoms, nOrb)

    !---------------------- Calculate S.B matrix  ----------------------
    call sb_matrix(sys%nAtoms, nOrb)

    !-------------------------- Debugging part -------------------------
    if(ldebug) then
      !if(myrank==0) then
      ! call debugging()
      !end if
    end if

    !----------------- Only create files with headers ------------------
    if(lcreatefiles) then
      if(rField == 0) call create_files()
      if((itype==7).or.(itype==8)) cycle
      call MPI_Finalize(ierr)
      call exit(0)
    end if
    !----------- Check if files exist to add results or sort -----------
    if( (lsortfiles .or. laddresults) .and. rField == 0 ) call check_files()

    !--------------------- Only sort existing files --------------------
    if(lsortfiles) then
      if(rField == 0) call sort_all_files()
      cycle
    end if

    !-- Only calculate long. and transv. currents from existing files --
    ! if(llgtv) then
    !   if(myrank==0) call read_calculate_lgtv_currents()
    !   if(itype==9) exit
    !   cycle
    ! end if

    !-------------- Only calculate SHA from existing files -------------
    ! if(lsha) then TODO: Re-include
    !   if(myrank==0) call read_currents_and_calculate_sha()
    !   if(itype==9) exit
    !   cycle
    ! end if

    !------------------------ Begin first test part --------------------
    if(rField == 0 .and. itype==0) then
      write(output%unit_loop,"('[main] FIRST TEST PART')")
      rho  = 0.d0
      mz  = 0.d0
      mp  = cZero
      ! Variables used in the hamiltonian

      call ldos()
      ! call debugging()
    end if

    !------------------------- Self-consistency ------------------------
    if(rField == 0) &
      call write_time(output%unit_loop,'[main] Time before self-consistency: ')

    call doSelfConsistency()

    !--------- Temporary execution to test this function ---------------

    ! call hamilt_sc(sys)

    ! if(rField == 0) &
    ! call band_structure(sys)

    !-------------------------------------------------------------------

    if(rField==0) &
      call write_time(output%unit_loop,'[main] Time after self-consistency: ')

    !-- Only calculate long. and transv. currents from existing files --
    ! if(llgtv) then TODO: Re-include
    !   if(myrank==0) call read_calculate_lgtv_currents()
    !   if(itype==9) exit
    !   cycle
    ! end if

    !========================== MAIN PROGRAM ===========================
    select case (itype)
    case(1)  ! Calculation of self-consistency only
      continue
    case (2) ! LDOS calculation
      call ldos()
    case (3) ! LDOS and coupling and a function of energy
      call ldos_and_coupling()
    case (4) ! Band structure (not parallelized)
      if(rField == 0) &
      call band_structure(sys)
    case (5) ! Iso-energy surface (for input = sys%Ef - Fermi surface)
      call fermi_surface(sys%Ef)
    case (6) ! Coupling tensor
      call coupling()
    case (7) ! Full magnetic susceptibility as a function of the frequency
      call calculate_chi()
    case (8) ! All responses to an electric field as a function of frequency
      call calculate_all()
    case (9) ! All responses to an electric field as a function of magnetic field
      call calculate_dc_limit()
    case(10) ! Calculation of Gilbert Damping by Kamberskys Torque Torque model
      call calculate_TCM()
    case(11) ! Propagation of ( H(k) + S.B(t) )
      call time_propagator(sys)
    end select
    !===================================================================

    ! Emergency stop after the calculation for a field is finished
    ! (not checked on last step)
    if(total_hw_npt1 /= 1 .and. hw_count < endField) &
      call check_stop("stopout",hw_count)
    !-------------------------- Closing files --------------------------
    if(total_hw_npt1 /= 1 .and.  itype /= 6 .and. rField == 0) &
      close(output%unit_loop)
  end do ! Ending of magnetic field loop

  !------------ Deallocating variables that depend on nAtoms------------
  call deallocate_Atom_variables()
  call deallocate_magnet_variables()

  !---------------------- Deallocating variables -----------------------
  !deallocate(r_nn, c_nn, l_nn)
  !deallocate(kbz,wkbz) !,kbz2d)

  !---------------------- Finalizing the program -----------------------
  if(myrank == 0) call write_time(output%unit,'[main] Finished on: ')
  if(myrank == 0) close(unit=output%unit)
  call MPI_Finalize(ierr)
  call exit(0)
end program TITAN
