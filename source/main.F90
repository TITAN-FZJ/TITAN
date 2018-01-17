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
  use mod_system         ! New
  use Lattice        ! New
  use mod_BrillouinZone  ! New
  use TightBinding   ! New
  use mod_SOC            ! New
  use mod_magnet
  use ElectricField
  use adaptiveMesh
  use AtomTypes
  use mod_gilbert_damping
  !use mod_define_system TODO: Re-include
  use mod_self_consistency
  !use mod_prefactors TODO: Re-include
  use EnergyIntegration
  use mod_progress
  use mod_mpi_pars
  use mod_Umatrix
  !use mod_lgtv_currents TODO: Re-include
  !use mod_sha TODO: Re-include
  !use mod_torques, only: ntypetorque
  implicit none

  integer           :: count_hw
  character(len=500) :: arg = ""

  !-------------------------- MPI initialization --------------------------
  call Initialize_MPI()

  !--------------------------- Starting program ---------------------------
  start_program = MPI_Wtime()

  !-------------------------- Reading parameters --------------------------
  if(command_argument_count() >= 1) call get_command_argument(1, arg)
  if(len_trim(arg) > 0) then
    call get_parameters(arg, sys)
  else
    call get_parameters("input",sys)
  end if

  if(myrank == 0) call write_time(outputunit,'[main] Started on: ')

  !-------------------- Useful constants and matrices --------------------
  call define_constants() ! TODO: Review

  !------------------- Set Loop and Integration Points -------------------
  call setMagneticLoopPoints()

  !----------- Creating bi-dimensional matrix of MPI processes  ----------
  !call setup_MPI_grid(itype, pn1, npt1, pnt,total_hw_npt1, npts, deltae, emin, emax)
  call genMPIGrid(parField,total_hw_npt1,parFreq,npt1 - skip_steps)
  !-------------------- Define the lattice structure ---------------------
  call initLattice(sys)
  ! Writing Positions into file
  if( lpositions .and. (myrank==0)) call writeLattice(sys)

  !-------------------- Generating k points in 2D BZ ---------------------
  call BZ % setup(sys%lbulk, sys%a1, sys%a2, sys%a3)
  ! Writing BZ points and weights into file
  if((lkpoints).and.(myrank==0)) call BZ % print()

  !---- Generating integration points of the complex energy integral -----
  call allocate_energy_points()
  call generate_imag_epoints()

  call generateAdaptiveMeshes()

  ! DFT tight binding mode has two additional vacuum layers
  !Npl_total = Npl
  !if(tbmode == 2) Npl_total = Npl + 2

  !--------------- Allocating variables that depend on Npl ---------------
  call allocate_magnet_variables(sys%nAtoms, nOrb)
  call allocLS(nOrb)
  call allocate_Npl_variables(sys%nAtoms) !TODO: Review

  !----------------------------- Dimensions ------------------------------
  dimsigmaNpl = 4 * sys%nAtoms
  dim = dimsigmaNpl * nOrb * nOrb

  !------------------------- Conversion arrays  --------------------------
  call initConversionMatrices(sys%nAtoms,nOrb)

  !------------------------- Defining the system -------------------------
  ! if(naddlayers/=0) then
  !   Npl_input = Npl-naddlayers+1 ! Npl is the total number of layers, including the added layers listed on inputcard
  ! else
  !   Npl_input = Npl
  ! end if
  write(strSites,fmt="(i0,'Sites')") sys%nAtoms
  ! if(tbmode == 2) call define_system()

  !---------------------- Tight Binding parameters -----------------------
  call initTightBinding(sys)

  !---- Initialize Stride Matrices for hamiltk and dtdksub --------------
  call initHamiltkStride(sys%nAtoms, nOrb)

  !---------------------- Lattice specific variables ---------------------
  call initElectricField(sys%a1, sys%a2, sys%a3)

  !call setup_long_and_trans_current_neighbors(sys) !TODO: Not implemented TODO: Re-include

  !------------------------- MAGNETIC FIELD LOOP -------------------------
  if(myrank == 0 .and. skip_steps_hw > 0) write(outputunit,"('[main] Skipping first ',i0,' field step(s)...')") skip_steps_hw

  do hw_count = startField, endField

    !------------ Opening files (general and one for each field) -----------
    if(rField == 0) then
      if(total_hw_npt1 == 1 .or. itype == 6 ) then
        outputfile_loop = outputfile
        outputunit_loop = outputunit
      else
        write(outputfile_loop, fmt="(a,a,i0)") trim(outputfile), '.', hw_count
        outputunit_loop = outputunit + hw_count
        open (unit=outputunit_loop, file=trim(outputfile_loop), status='replace')
      end if
    end if

    !------------------- Initialize Magnetic Field --------------------
    call initMagneticField(sys%nAtoms)

    !------------- Update fieldpart ----------------
    call set_fieldpart(hw_count)

    !------------------- Tests for coupling calculation --------------------
    if(itype == 6 .and. myrank == 0) then
      if(nmaglayers==0) call abortProgram("[main] No magnetic layers for coupling calculation!")
      if(lfield) call abortProgram("[main] Coupling calculation is valid for no external field!")
    end if

    !------- Writing parameters and data to be calculated on screen --------
    if(rField == 0) call iowrite(sys)

    !---- L matrix in global frame for given quantization direction ----
    call l_matrix()

    !------ Calculate L.S matrix for the given quantization direction ------
    call updateLS(theta, phi, nOrb)

    !------ Calculate L.B matrix for the given quantization direction ------
    call lb_matrix(sys%nAtoms, nOrb)

    !------------------------ Calculate S.B matrix  ------------------------
    call sb_matrix(sys%nAtoms, nOrb)

    !---------------------------- Debugging part ---------------------------
    ! if(ldebug) then TODO: Re-Include
    !   !if(myrank==0) then
    !   call debugging()
    !   !end if
    ! end if

    !------------------- Only create files with headers --------------------
    if(lcreatefiles) then
      if(rField == 0) call create_files()
      if((itype==7).or.(itype==8)) cycle
      call MPI_Finalize(ierr)
      call exit(0)
    end if
    !------------- Check if files exist to add results or sort -------------
    if( (lsortfiles .or. laddresults) .and. rField == 0 ) call check_files()

    !----------------------- Only sort existing files ----------------------
    if(lsortfiles) then
      if(rField == 0) call sort_all_files()
      cycle
    end if

    ! !---- Only calculate long. and transv. currents from existing files ----
    ! if(llgtv) then
    !   if(myrank==0) call read_calculate_lgtv_currents()
    !   if(itype==9) exit
    !   cycle
    ! end if

    !---------------- Only calculate SHA from existing files ---------------
    ! if(lsha) then TODO: Re-include
    !   if(myrank==0) call read_currents_and_calculate_sha()
    !   if(itype==9) exit
    !   cycle
    ! end if

    !-------------------------- Begin first test part ----------------------
    if(rField == 0 .and. itype==0) then
      write(outputunit_loop,"('[main] FIRST TEST PART')")
      n_t = 0.d0
      mz  = 0.d0
      mp  = cZero
      mm  = conjg(mp)
      ! Variables used in the hamiltonian

      call ldos()
      ! call debugging()
    end if

    if(lcreatefolders) call create_folder()


    ! Self-consistency

    ! Time now
    if(rField == 0)  call write_time(outputunit_loop,'[main] Time before self-consistency: ')

    call doSelfConsistency()

    ! Time now
    if(rField==0)  call write_time(outputunit_loop,'[main] Time after self-consistency: ')

    !---- Only calculate long. and transv. currents from existing files ----
    ! if(llgtv) then TODO: Re-include
    !   if(myrank==0) call read_calculate_lgtv_currents()
    !   if(itype==9) exit
    !   cycle
    ! end if

    !============================= MAIN PROGRAM =============================
    select case (itype)
    case(1) ! Calculation of self-consistency only
      continue
    case (2) !Debugging case
      !call debugging() TODO: Re-Include
      continue
    case (3) ! Ca
      call ldos_and_coupling()
    case (4) !
      if(rField == 0) call band_structure(sys)
    case (5) !
      call fermi_surface(sys%Ef)
    case (6) !
      call coupling()
    case (7) !
      call calculate_chi()
    case (8) !
      call calculate_all()
    case (9) !
      call calculate_dc_limit()
    case(10) ! Calculation of Gilbert Damping by Kamberskys Torque Torque model
      call calculate_gilbert_damping()
    end select
    !========================================================================

    ! Emergency stop after the calculation for a Npl or field is finished (don't check on last step)
    if(total_hw_npt1 /= 1 .and. count_hw < endField) call check_emergency_stop(sys%nAtoms, hw_list, hw_count)
    !---------------------------- Closing files ----------------------------
    if(total_hw_npt1 /= 1 .and.  itype /= 6 .and. rField == 0) close(outputunit_loop)
    !---------------------- Ending magnetic field loop ---------------------
  end do

  !-------------- Deallocating variables that depend on Npl --------------
  call deallocate_Npl_variables()

  !----------------------- Deallocating variables ------------------------
  !deallocate(r_nn, c_nn, l_nn)
  !deallocate(kbz,wkbz) !,kbz2d)

  !----------------------- Finalizing the program ------------------------
  if(myrank == 0) call write_time(outputunit,'[main] Finished on: ')
  if(myrank == 0) close(unit=outputunit)
  call MPI_Finalize(ierr)
  call exit(0)
end program TITAN
