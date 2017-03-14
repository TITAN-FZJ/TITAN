program DHE
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_io
  use mod_magnet
  use mod_define_system
  use mod_tight_binding
  use mod_self_consistency
  use mod_prefactors
  use mod_generate_epoints
  use mod_generate_kpoints
  use mod_lattice
  use mod_progress
  use mod_mpi_pars
  use mod_lgtv_currents
  use mod_sha
  use mod_torques, only: ntypetorque
  implicit none
  integer           :: i,j,sigma,mu,nu,count_hw
  logical           :: lsuccess

#ifndef _UFF
!$  integer :: provided
#endif
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
!--------------------------- Starting program --------------------------
  start_program = MPI_Wtime()
  call read_input_file()
  call read_output_filename()
  if(myrank.eq.0) then
    open (unit=outputunit, file=trim(outputdhe), status='unknown')
    call write_time(outputunit,'[main] Started on: ')
  end if
!---------------------------- Getting hostname -------------------------
#ifdef _JUQUEEN
  call hostnm_(host)
#else
  call hostnm(host)
#endif
!-------------------- Useful constants and matrices --------------------
  call define_constants()
!------------------------- Reading parameters --------------------------
  call get_parameters()
!----------- Creating bi-dimensional matrix of MPI processes  ----------
  if((itype.eq.1).or.(itype.eq.6)) then ! Create column for field loop (no energy integration)
    call build_cartesian_grid_field(pn1)
  end if
  if((itype.ge.3).and.(itype.le.5)) then ! Create column for field loop (no energy integration)
    call build_cartesian_grid_field(1)
  end if
  if((itype.ge.7).and.(itype.le.8)) then ! Create matrix for energy dependence and integration
    call build_cartesian_grid()
    call build_cartesian_grid_field(npt1*pnt)
  end if
  if(itype.eq.9) then ! Create matrix for dclimit
    call build_cartesian_grid_field(pnt)
    call MPI_COMM_DUP(MPI_Comm_Col_hw,MPI_Comm_Col,ierr)
    call MPI_COMM_DUP(MPI_Comm_Row_hw,MPI_Comm_Row,ierr)
    myrank_col = myrank_col_hw
    myrank_row = myrank_row_hw
  end if
!-------------------- Define the lattice structure ---------------------
  call next_neighbour_init(a0*a1, a0*a2, a0*a3, pln_dir)
!-------------------- Generating k points in 2D BZ ---------------------
  call generate_kpoints(pln_dir)

  ! Writing BZ points and weights into files
  if((lkpoints).and.(myrank.eq.0)) then
    call write_kpoints_to_file()
  end if

  stop

!---- Generating integration points of the complex energy integral -----
  call generate_imag_epoints()
!------------------------ NUMBER OF PLANES LOOP ------------------------
  number_of_planes: do Npl = Npl_i,Npl_f
!--------------- Allocating variables that depend on Npl ---------------
    call allocate_Npl_variables()
!----------------------------- Dimensions ------------------------------
    dimsigmaNpl = 4*Npl
    dim = dimsigmaNpl*9*9
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
!------------------------- Defining the system -------------------------
    if(naddlayers.ne.0) then
      Npl_input = Npl-naddlayers+1 ! Npl is the total number of layers, including the added layers listed on inputcard
    else
      Npl_input = Npl
    end if
    call define_system()
!---------------------- Tight Binding parameters -----------------------
    call rs_hoppings()
!---------------------- Lattice specific variables ---------------------
    call lattice_definitions()
!------------------------- MAGNETIC FIELD LOOP -------------------------
    if((myrank.eq.0).and.(skip_steps_hw.gt.0)) write(outputunit,"('[main] Skipping first ',i0,' field step(s)...')") skip_steps_hw
    hw_loop: do count_hw=1+skip_steps_hw,MPIsteps_hw
      hw_count = 1 + myrank_col_hw + MPIpts_hw*(count_hw-1)
!------------ Opening files (general and one for each field) -----------
      if(myrank_row_hw.eq.0) then
        if((total_hw_npt1.eq.1).or.(itype.eq.6)) then
          outputdhe_loop  = outputdhe
          outputunit_loop = outputunit
        else
          write(outputdhe_loop,fmt="(a,a,i0)") trim(outputdhe),'.',hw_count
          outputunit_loop = outputunit+hw_count
          open (unit=outputunit_loop, file=trim(outputdhe_loop), status='unknown')
        end if
      end if
!---------------------- Turning off field for hwa=0 --------------------
      if(abs(hw_list(hw_count,1)).lt.1.d-8) then
        lfield = .false.
        ntypetorque = 2
        if((llinearsoc).or.(.not.SOC).and.(myrank_row_hw.eq.0)) write(outputdhe_loop,"('[main] WARNING: No external magnetic field is applied and SOC is off/linear order: Goldstone mode is present!')")
      else
        lfield = .true.
      end if
      sb   = zero
      lb   = zero
      hhwx = 0.d0
      hhwy = 0.d0
      hhwz = 0.d0
!--------------------- Defining the magnetic fields --------------------
      if(lfield) then
        ! Variables of the hamiltonian
        ! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
        ! to take into account the fact that we are considering negative
        ! external fields to get the peak at positive energies
        do i=1,Npl+2
          hhwx(i)  =-0.5d0*hwscale(i)*hw_list(hw_count,1)*sin((hw_list(hw_count,2)+hwtrotate(i))*pi)*cos((hw_list(hw_count,3)+hwprotate(i))*pi)*tesla
          hhwy(i)  =-0.5d0*hwscale(i)*hw_list(hw_count,1)*sin((hw_list(hw_count,2)+hwtrotate(i))*pi)*sin((hw_list(hw_count,3)+hwprotate(i))*pi)*tesla
          hhwz(i)  =-0.5d0*hwscale(i)*hw_list(hw_count,1)*cos((hw_list(hw_count,2)+hwtrotate(i))*pi)*tesla
          if(abs(hhwx(i)).lt.1.d-8) hhwx(i) = 0.d0
          if(abs(hhwy(i)).lt.1.d-8) hhwy(i) = 0.d0
          if(abs(hhwz(i)).lt.1.d-8) hhwz(i) = 0.d0
        end do
        ! Testing if hwscale is used
        lhwscale = any(abs(hwscale(1:Npl)-1.d0).gt.1.d-8)
        ! Testing if hwrotate is used
        lhwrotate = (any(abs(hwtrotate(1:Npl)).gt.1.d-8).or.any(abs(hwprotate(1:Npl)).gt.1.d-8))
      end if
!------------------- Tests for coupling calculation --------------------
      if(itype.eq.6) then
        if(nmaglayers.eq.0) then
          if(myrank.eq.0) write(outputunit,"('[main] No magnetic layers for coupling calculation!')")
          call MPI_Finalize(ierr)
          stop
        end if
        if(lfield) then
          if(myrank.eq.0) write(outputunit,"('[main] Coupling calculation is valid for no external field!')")
          call MPI_Finalize(ierr)
          stop
        end if
      end if
!------- Writing parameters and data to be calculated on screen --------
      if(myrank_row_hw.eq.0) call iowrite()
!---- L matrix in spin coordinates for given quantization direction ----
      call lp_matrix()
!------ Calculate L.S matrix for the given quantization direction ------
      call ls_matrix()
!------ Calculate L.B matrix for the given quantization direction ------
      if((lfield).and.(.not.lnolb)) call lb_matrix()
!------------------------ Calculate S.B matrix  ------------------------
      if(lfield) call sb_matrix()
!---------------------------- Debugging part ---------------------------
      if(ldebug) then
!       if(myrank.eq.0) then
        call debugging()
!       end if
      end if
!------------------- Only create files with headers --------------------
      if(lcreatefiles) then
        if(myrank.eq.0) then
          call create_files()
          if((itype.eq.7).or.(itype.eq.8)) cycle
        end if
        call MPI_Finalize(ierr)
        stop
      end if
!------------- Check if files exist to add results or sort -------------
      if(((lsortfiles).or.(laddresults)).and.(myrank.eq.0)) call check_files()
!----------------------- Only sort existing files ----------------------
      if(lsortfiles) then
        if(myrank.eq.0) call sort_all_files()
        cycle
      end if
! !---- Only calculate long. and transv. currents from existing files ----
!       if(llgtv) then
!         if(myrank.eq.0) call read_calculate_lgtv_currents()
!         if(itype.eq.9) exit
!         cycle
!       end if
!---------------- Only calculate SHA from existing files ---------------
      if(lsha) then
        if(myrank.eq.0) call read_currents_and_calculate_sha()
        if(itype.eq.9) exit
        cycle
      end if
!-------------------------- Begin first test part ----------------------
      if((myrank.eq.0).and.(itype.eq.0)) then
        write(outputunit,"('[main] FIRST TEST PART')")
        call debugging()
      end if
!--------------------------- Self-consistency --------------------------
      ! Trying to read previous shifts and m from files
      call read_previous_results(lsuccess)
      ! Rotate the magnetization to the direction of the field
      ! (useful for SOC=F)
      ! Only when some previous result was read from file (lsuccess=.true.)
      if((lrotatemag).and.(lsuccess)) then
        call rotate_magnetization_to_field()
        cycle
      end if
      ! Self-consistency
      if(lselfcon) call do_self_consistency()
      ! Writing new eps1 and mz to file after self-consistency is done
      if(.not.lontheflysc) call write_sc_results()
      ! Calculating ground state Orbital Angular Momentum
      if(lGSL) call L_gs()
      ! Writing self-consistency results on screen
      if(myrank_row_hw.eq.0)  call write_sc_results_on_screen()
      ! Time now
      if(myrank_row_hw.eq.0)  call write_time(outputunit_loop,'[main] Time after self-consistency: ')
!---- Only calculate long. and transv. currents from existing files ----
      if(llgtv) then
        if(myrank.eq.0) call read_calculate_lgtv_currents()
        if(itype.eq.9) exit
        cycle
      end if
!============================= MAIN PROGRAM ============================
      main_program: select case (itype)
      case (2)
        call debugging()
      case (3)
        if(myrank_row_hw.eq.0) call ldos_and_coupling()
      case (4)
        if(myrank_row_hw.eq.0) call band_structure()
      case (5)
        if(myrank_row_hw.eq.0) call fermi_surface(Ef)
      case (6)
        call coupling()
      case (7)
        call calculate_chi()
      case (8)
        call calculate_all()
      case (9)
        call calculate_dc_limit()
      end select main_program

      ! Emergency stop after the calculation for a Npl or field is finished (don't check on last step)
      if(((Npl_f.ne.Npl_i).or.(total_hw_npt1.ne.1)).and.((Npl.lt.Npl_f).or.(count_hw.lt.MPIsteps_hw))) call check_emergency_stop()
!---------------------------- Closing files ----------------------------
      if(((total_hw_npt1.ne.1).and.(itype.ne.6)).and.(myrank_row_hw.eq.0)) close(outputunit_loop)
!---------------------- Ending magnetic field loop ---------------------
    end do hw_loop
!-------------- Deallocating variables that depend on Npl --------------
    call deallocate_Npl_variables()
  end do number_of_planes
!----------------------- Deallocating variables ------------------------
  deallocate(r0,c0,r1,c1,r2,c2)
  deallocate(kbz,wkbz,kbz2d)
!----------------------- Finalizing the program ------------------------
  if(myrank.eq.0) call write_time(outputunit,'[main] Finished on: ')
  if(myrank.eq.0) close(unit=outputunit)
  call MPI_Finalize(ierr)
  if((ierr.ne.0).and.(myrank.eq.0)) write(outputunit,"('[main] Something went wrong in the parallelization! ierr = ',i0)") ierr
!=======================================================================
  stop
end program DHE
