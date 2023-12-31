#include "macros.fi"

!> Main program:
program TITAN
  !! Initialize MPI, reads input, define lattice quantities then
  !! loops over number of planes and magnetic field to do the
  !! magnetic self-consistency and calculate either ground state
  !! quantities or response functions
  use mod_kind,                only: dp
  use mod_constants,           only: cZero,define_constants
  use mod_parameters,          only: output,lpositions,lcreatefolders,parField,parFreq,nEner1,skip_steps,ldebug, &
                                     kp_in,kptotal_in,eta,leigenstates,itype,theta,phi, &
                                     laddresults,lsortfiles,lcreatefiles,arg,tbmode,lfixEf,addelectrons
  use mod_io,                  only: get_parameters,iowrite
  use mod_logging,             only: log_error,log_warning
  use Lattice,                 only: initLattice,writeLattice
  use mod_BrillouinZone,       only: realBZ,countBZ
  use mod_SOC,                 only: llinearsoc,SOC,updateLS
  use mod_magnet,              only: total_hw_npt1,skip_steps_hw,hw_count,lfield,rho,mz,mp,setMagneticLoopPoints,&
                                     allocate_magnet_variables,initMagneticField,set_fieldpart,l_matrix,lb_matrix,&
                                     sb_matrix,lconstraining_field
  use ElectricField,           only: initElectricField
  use adaptiveMesh,            only: generateAdaptiveMeshes,freeLocalEKMesh
  use mod_TCM,                 only: calculate_TCM
  use mod_ldos,                only: ldos,ldos_and_coupling
  use mod_band_structure,      only: band_structure
  use mod_self_consistency,    only: doSelfConsistency
  use EnergyIntegration,       only: pn1,allocate_energy_points,generate_imag_epoints
  use mod_progress,            only: start_program,write_time
  use mod_mpi_pars,            only: MPI_Wtime,MPI_COMM_WORLD,ierr,myrank,startField,endField,rField,Initialize_MPI,&
                                     genMPIGrid,abortProgram
  use mod_polyBasis,           only: read_basis
  use TightBinding,            only: initTightBinding
  use mod_fermi_surface,       only: fermi_surface
  use mod_check_stop,          only: check_stop
  use mod_Atom_variables,      only: allocate_Atom_variables
  use mod_tools,               only: rtos,itos
  use mod_init_expec,          only: calc_init_expec_SK,calc_init_expec_dft
  use mod_time_propagator,     only: time_propagator
  use mod_superconductivity,   only: lsuperCond,supercond,allocate_supercond_variables
  use mod_System,              only: s=>sys,init_Hamiltk_variables,initConversionMatrices,allocate_basis_variables
#ifdef _GPU
  use nvtx,                    only: nvtxStartRange,nvtxEndRange
  use mod_cuda,                only: num_gpus,result,create_handle,cudaGetDeviceCount,cudaSetDevice,cudaGetErrorString 
#endif
  !use mod_prefactors TODO: Re-include
  !use mod_lgtv_currents TODO: Re-include
  !use mod_sha TODO: Re-include
  !use mod_torques, only: ntypetorque
  implicit none

  external :: create_folder,create_files,check_files,sort_all_files
  external :: coupling,calculate_chi,calculate_all,calculate_dc_limit, real_coupling
  external :: setLoops,endTITAN,MPI_Barrier

  !------------------------ MPI initialization -------------------------
  call Initialize_MPI()

  !------------------------- Starting program --------------------------
  start_program = MPI_Wtime()

  !------------------------ Reading parameters -------------------------
  if(command_argument_count() >= 1) call get_command_argument(1, arg)
  if(len_trim(arg) > 0) then
    call get_parameters(arg, s)
  else
    call get_parameters("input",s)
  end if
  arg = ""

#ifdef _GPU
  !------------------------ GPU initialization -------------------------
  ! Getting number of GPUs seen by currrent MPI
  CUDA_CALL( cudaGetDeviceCount(num_gpus) )

  if(num_gpus == 0) &
    call log_error("main", "Something wrong in the setup: number of GPUS is zero!")
  ! Setting up the device for this process
  CUDA_CALL( cudaSetDevice(mod(myrank,num_gpus)) )
  ! Creating handle for cusolver
  call create_handle()
  !$acc init device_type(acc_device_nvidia)
  arg = trim(arg) // " GPU"
#endif

#ifdef DEBUG
  arg = trim(arg) // " DEBUG"
#endif

  !---------------- Writing time after initializations -----------------
  if(myrank == 0) call write_time('[main]' // trim(arg) // ' Started on: ',output%unit)

  !------------------- Define the lattice structure --------------------
  call read_basis("basis", s)
  call initLattice(s)
  write(output%Sites,fmt="(i0,'Sites')") s%nAtoms

  ! Writing Positions into file
  if( lpositions .and. (myrank==0) .and. (tbmode==1) ) call writeLattice(s)

  !------------- Creating folders for current calculation ------------
  if(lcreatefolders) then
    if (myrank == 0) call create_folder()
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
  end if

  !------------------ Set Loops and Integration Points -----------------
  call setMagneticLoopPoints()
  call setLoops(s)

  !------------------ Creating grid of MPI processes  ------------------
  !call setup_MPI_grid(itype, pn1, nEner1, pnt,total_hw_npt1, nEner, deltae, emin, emax)
  call genMPIGrid(parField, total_hw_npt1, parFreq, nEner1 - skip_steps)

  !--- Generating integration points of the complex energy integral ----
  call allocate_energy_points()
  call generate_imag_epoints()

  !-------------------------- Debugging part -------------------------
  if(ldebug) then
    ! call initTightBinding(s)
    ! call allocate_magnet_variables(s%nAtoms, s%nOrb)
    ! call hamilt_sc(s)
    !call hamilt_sc(s,rho,mp,mx,my,mz,s,kp,hk)
    !if(myrank==0) then
    ! call debugging()
    !end if
    !MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    ! call deallocate_magnet_variables()
    ! stop
  end if

  !---------------- Reading Tight Binding parameters -------------------
  call initTightBinding(s)

  !------------------- Useful constants and matrices -------------------
  call define_constants(s) ! TODO: Review

  !----------------- Tests for coupling calculation ------------------
  if((lconstraining_field).and.(sum(abs(s%Types(:)%Um))<1.e-8).and.(tbmode==1).and.(myrank == 0)) &
    call abortProgram("[main] Constraining fields need Um to induce a magnetic moment!")

  call flush(output%unit)

  ! Testing if add electrons is used with lfixEf or Un/=0
  if(lfixEf) then
    if(abs(addelectrons)>1.e-6) & 
      call log_warning("main", "addelectrons does not affect the calculation when lfixEf is used!" )
  else
    if((sum(abs(s%Types(:)%Un))>1.e-8).and.(abs(addelectrons)>1.e-8_dp)) &
      call log_error("main", "addelectrons is incompatible with Un! Either use Un=0 for all elements or addelectrons=0." )
  end if
  if(abs(addelectrons)>1.e-6_dp) then
    write(output%unit,"('[main] Electrons to add (or remove): ',es14.7)") addelectrons
    s%totalOccupation = s%totalOccupation + addelectrons
    output%suffix = "_add" // trim(rtos(addelectrons,"(es8.1)"))
  end if

  !-- Calculating initial values in the SK tight-binding hamiltonian ---
  if((tbmode==1).and.(.not.lsortfiles)) call calc_init_expec_SK(s)

  !----------- Generating k points for real axis integration -----------
  realBZ % nkpt_x = kp_in(1)
  realBZ % nkpt_y = kp_in(2)
  realBZ % nkpt_z = kp_in(3)
  call realBZ % countBZ(s)

  !---- Generating k meshes points for imaginary axis integration ------
  call generateAdaptiveMeshes(s,pn1)

  !------------ Allocating variables that depend on nAtoms -------------
  call allocate_basis_variables(s)
  call allocate_magnet_variables(s%nAtoms, s%nOrb)
  call allocate_supercond_variables(s%nAtoms, s%nOrb)
  call allocate_Atom_variables(s%nAtoms,s%nOrb)

  !------- Initialize Stride Matrices for hamiltk and dtdksub ----------
  call init_Hamiltk_variables(s,supercond)

  !---------- Conversion arrays for dynamical quantities ---------------
  call initConversionMatrices(s)

  !--------------------- Lattice specific variables --------------------
  call initElectricField(s%a1, s%a2, s%a3)

  !call setup_long_and_trans_current_neighbors(s) !TODO: Not implemented TODO: Re-include
  !-------------------------- Filename strings -------------------------
  write(output%info,"('_norb=',i0,'_nkpt=',i0,'_eta=',a)") s%nOrb,kptotal_in,trim(rtos(eta,"(es8.1)"))
  if(leigenstates) output%info = trim(output%info) // "_ev"
  if(lfixEf) output%info = trim(output%info) // "_fixEf"

  !-- Calculating initial values in the DFT tight-binding hamiltonian --
  if((tbmode==2).and.(.not.lsortfiles)) call calc_init_expec_dft(s)

  !--------------- Extra Token for Superconductivity files -------------
  if(lsuperCond) output%info = trim(output%info) // "_sc"

  !------------------------ MAGNETIC FIELD LOOP ------------------------
  if(myrank == 0 .and. skip_steps_hw > 0) &
    write(output%unit,"('[main] Skipping first ',i0,' field step(s)...')") skip_steps_hw

  call flush(output%unit)

  do hw_count = int(startField,4), int(endField,4)
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
    call initMagneticField(s)
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
    if(rField == 0) call iowrite(s)

    !---- L matrix in global frame for given quantization direction ----
    call l_matrix(s)

    !---- Calculate L.S matrix for the given quantization direction ----
    if(SOC) call updateLS(s,theta, phi)

    !---- Calculate L.B matrix for the given quantization direction ----
    call lb_matrix(s)

    !---------------------- Calculate S.B matrix  ----------------------
    call sb_matrix(s)

    !-------------------------- Debugging part -------------------------
    if(ldebug) then
      ! if(rField == 0) &
      !   call band_structure(s)
      ! call endTITAN()
      !if(myrank==0) then
      ! call debugging()
      !end if
    end if

    !----------------- Only create files with headers ------------------
    if(lcreatefiles) then
      if(rField == 0) call create_files()
      if((itype==7).or.(itype==8)) cycle
      call endTITAN()
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
    ! if(rField == 0 .and. itype==0) then
    !   write(output%unit_loop,"('[main] FIRST TEST PART')")
    !   rho  = 0._dp
    !   mz  = 0._dp
    !   mp  = cZero
    !   ! Variables used in the hamiltonian
    !
    !   call ldos()
    !   ! call debugging()
    ! end if

    !------------------------- Self-consistency ------------------------
    if(rField == 0) &
      call write_time('[main] Time before self-consistency: ',output%unit_loop)


#ifdef _GPU
    ! Starting marker of self-consistency for profiler
    call nvtxStartRange("Self-consistency",0)
#endif

    call doSelfConsistency()

#ifdef _GPU
    ! End of self-consistency marker
    call nvtxEndRange
#endif

    if(rField==0) &
      call write_time('[main] Time after self-consistency: ',output%unit_loop)

    call flush(output%unit)


    !--------- Temporary execution to test this function ---------------

    ! call hamilt_sc(s)

    ! if(rField == 0) &
    ! call band_structure(s)

    if(rField == 0 .and. itype==0) then
      write(output%unit_loop,"('[main] FIRST TEST PART')")
      rho  = 0._dp
      mz  = 0._dp
      mp  = cZero
      ! Variables used in the hamiltonian

      ! call expectation_values_eigenstates(s,rho,mp,mx,my,mz,delta_sc)

      ! call debugging()
    end if

    !-------------------------------------------------------------------

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
        call band_structure(s)
    case (5) ! Iso-energy surface (for input = s%Ef - Fermi surface)
      call fermi_surface()
    case (6) ! Coupling tensor
      call coupling()
    case (7) ! Full magnetic susceptibility as a function of the frequency
      call calculate_chi()
    case (8) ! All responses to an electric field as a function of frequency
      call calculate_all()
    case (9) ! All responses to an electric field as a function of magnetic field
      call calculate_dc_limit()
    case(10) ! Calculation of Gilbert Damping using Torque correlation models
      call calculate_TCM()
    case(11) ! Propagation of ( H(k) + S.B(t) )
      call time_propagator(s)
    case(12) ! Coupling tensor in real space
      call real_coupling()
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

  call endTITAN()
end program TITAN




subroutine endTITAN()
  use mod_parameters,          only: arg,output,leigenstates
  use mod_progress,            only: write_time
  use mod_BrillouinZone,       only: realBZ
  use mod_magnet,              only: deallocate_magnet_variables
  use adaptiveMesh,            only: deallocateAdaptiveMeshes
  use EnergyIntegration,       only: deallocate_energy_points
  use mod_Umatrix,             only: deallocate_Umatrix
  use mod_Atom_variables,      only: deallocate_Atom_variables
  use mod_hamiltonian,         only: deallocate_hamiltonian
  use mod_superconductivity,   only: deallocate_supercond_variables
  use mod_System,              only: deallocate_System_variables
  use mod_mpi_pars,            only: ierr,myrank
#ifdef _GPU
  use mod_cuda,                only: destroy_handle
#endif
  implicit none

  external :: MPI_Finalize,deallocateLoops

  !----------------------- Deallocating variables ----------------------
  call deallocate_Atom_variables()
  call deallocate_magnet_variables()
  call deallocate_System_variables()
  call deallocate_Umatrix()
  call deallocate_energy_points()
  call deallocateLoops()
  call deallocate_hamiltonian()
  if(leigenstates) call realBZ%free()
  call deallocateAdaptiveMeshes()
  call deallocate_supercond_variables()

  !---------------------- Deallocating variables -----------------------
  !deallocate(r_nn, c_nn, l_nn)
  !deallocate(kbz,wkbz) !,kbz2d)

  !---------------------- Finalizing the program -----------------------
#ifdef _GPU
  call destroy_handle()
#endif

  if(myrank == 0) call write_time('[endTITAN]' // trim(arg) // ' Finished on: ',output%unit)
  if(myrank == 0) close(unit=output%unit)
  call MPI_Finalize(ierr)
  call exit(0)
end subroutine endTITAN
