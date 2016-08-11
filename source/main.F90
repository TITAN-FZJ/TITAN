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
  implicit none
  integer           :: i,j,iw,sigma,mu,nu
  integer           :: AllocateStatus
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
!     call date_and_time(date, time, zone, values)
!     write(outputunit,"('[main] Started on: ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s')") values(3),values(2),values(1),values(5),values(6),values(7)
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
    call build_cartesian_grid_field(numprocs)
  end if
  if(itype.eq.9) then ! Create matrix for dclimit
    call build_cartesian_grid_field(pnt)
    call MPI_COMM_DUP(MPIComm_Col_hw,MPIComm_Col,ierr)
    call MPI_COMM_DUP(MPIComm_Row_hw,MPIComm_Row,ierr)
    myrank_col = myrank_col_hw
    myrank_row = myrank_row_hw
  end if
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
    if(myrank.eq.0) write(outputunit,"('[main] Lattice not defined: ',a,'!')") lattice
    call MPI_Finalize(ierr)
    stop
  end select
  ! Writing BZ points and weights into files
  if((lkpoints).and.(myrank.eq.0)) then
    call write_kpoints_to_file()
  end if
!---- Generating integration points of the complex energy integral -----
  call generate_imag_epoints()
!--------------------------------- Loops -------------------------------
  number_of_planes: do Npl = Npl_i,Npl_f
    hw_loop: do hw_count=myrank_col_hw+1,total_hw_npt1,MPIpts_hw
!------------------------- Opening specific files ----------------------
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
        if((llinearsoc).or.(.not.SOC).and.(myrank_row_hw.eq.0)) write(outputdhe_loop,"('[main] WARNING: No external magnetic field is applied and SOC is off/linear order: Goldstone mode is present!')")
      else
        lfield = .true.
      end if
!--------------------- Defining the magnetic fields --------------------
      if(lfield) then
        hwx  = hw_list(hw_count,1)*sin(hw_list(hw_count,2)*pi)*cos(hw_list(hw_count,3)*pi)
        hwy  = hw_list(hw_count,1)*sin(hw_list(hw_count,2)*pi)*sin(hw_list(hw_count,3)*pi)
        hwz  = hw_list(hw_count,1)*cos(hw_list(hw_count,2)*pi)
        if(abs(hwx).lt.1.d-8) hwx = 0.d0
        if(abs(hwy).lt.1.d-8) hwy = 0.d0
        if(abs(hwz).lt.1.d-8) hwz = 0.d0
        ! Variables of the hamiltonian
        ! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
        ! to take into account the fact that we are considering negative
        ! external fields to get the peak in positive energies
        hhwx  =-0.5d0*hwx*tesla
        hhwy  =-0.5d0*hwy*tesla
        hhwz  =-0.5d0*hwz*tesla
!       else
!         if(hw_count.ge.2) exit hw_loop
      end if
!--------------- Allocating variables that depend on Npl ---------------
      allocate( sigmai2i(4,Npl),sigmaimunu2i(4,Npl,9,9),sigmaijmunu2i(4,Npl,Npl,9,9),eps1(Npl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(outputunit,"('[main] Not enough memory for: sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      allocate( mx(Npl),my(Npl),mz(Npl),mvec_cartesian(Npl,3),mvec_spherical(Npl,3),hdel(Npl),mp(Npl),hdelp(Npl),mm(Npl),hdelm(Npl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(outputunit,"('[main] Not enough memory for: mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      allocate( mabs(Npl),mtheta(Npl),mphi(Npl),labs(Npl),ltheta(Npl),lphi(Npl),lpabs(Npl),lptheta(Npl),lpphi(Npl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(outputunit,"('[main] Not enough memory for: mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      allocate( mmlayer(Npl+2),layertype(Npl+2),U(Npl+2),mmlayermag(Npl+2),lambda(Npl+2),npart0(Npl+2), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(outputunit,"('[main] Not enough memory for: mmlayer,layertype,U,mmlayermag,lambda,npart0')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
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
      call define_system()
!---------------------- Tight Binding parameters -----------------------
      call rs_hoppings()
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
!---------------------- Lattice specific variables ---------------------
      call lattice_definitions()
!------- Writing parameters and data to be calculated on screen --------
      if(myrank_row_hw.eq.0) call iowrite()
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
        call debugging()
!       end if
        call MPI_Finalize(ierr)
        if ((ierr.ne.0).and.(myrank.eq.0)) write(outputunit,"('[main] Something went wrong in the parallelization! ierr = ',i0)") ierr
        stop
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
!----------------- Check if files exist to add results -----------------
      if((laddresults).and.(myrank.eq.0)) call check_files()
!-------------------------- Begin first test part ----------------------
      if((myrank.eq.0).and.(itype.eq.0)) then
        write(outputunit,"('[main] FIRST TEST PART')")

!         ! Parameters: center of band, magnetization, correlation functions
!         eps1  = 0.d0
!         mz   = 0.d0
!         mp = zero
!         ! Variables used in the hamiltonian
!         do i=1,Npl
!           hdel(i)   = 0.5d0*U(i+1)*mz(i)
!           hdelp(i)  = 0.5d0*U(i+1)*mp(i)
!         end do
!         hdelp = zero
!         hdelm = conjg(hdelp)

!         emin = -2.d0  ! given in eV
!         emax = 2.d0   ! given in eV
!         npts = 400
!         npt1 = npts+1
!         test1_energy_loop2: do count=1,npt1
!           e = emin + (count-1)*deltae
!           write(outputunit,"(i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e

!           ! Turning off SOC
!   !         lambda = 0.d0

!           call ldos(e,ldosu,ldosd,Jij)

!         end do test1_energy_loop2


    !   Finalizing program
        if(myrank.eq.0) call write_time(outputunit,'[main] Finished on: ')
        call MPI_Finalize(ierr)
        stop

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

!============================= MAIN PROGRAM ============================
      main_program: select case (itype)
      case (2)
!----------------------------- Begin test part -----------------------
  !       if(myrank.eq.0) then
          if(myrank.eq.0) write(outputunit_loop,"('TESTING')")
          if(myrank.eq.0) write(outputunit_loop,"('Npl = ',i0)") Npl
  !         call diamagnetic_current()

  !         test2_energy_loop: do count=1,npt1
  !           e = emin + (count-1)*deltae
  !           write(outputunit,"(i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e


  !         end do test2_energy_loop
  !       end if
!---------------------------- End test part ----------------------------

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
!-------------- Deallocating variables that depend on Npl --------------
      deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1)
      deallocate(mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm)
      deallocate(mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi)
      if(lGSL) deallocate(lxm,lym,lzm,lxpm,lypm,lzpm)
      deallocate(mmlayer,layertype,U,mmlayermag,lambda,npart0)
      select case (plnn)
      case(1)
        deallocate(t00,t01)
      case(2)
        deallocate(t00,t01,t02)
      end select

      ! Emergency stop after the calculation for certain Npl is finished
      if((Npl_f.ne.Npl_i).and.(total_hw_npt1.ne.1)) then
        open(unit=911, file="stopout", status='old', iostat=iw)
        if(iw.eq.0) then
          close(911)
          write(outputunit,"('[main] Emergency ""stopout"" file found! Stopping after Npl = ',i0,', hwa = ',es9.2,' hwt=',f5.2,' hwp=',f5.2,'...')") Npl,hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
          call system ('rm stopout')
          write(outputunit,"('[main] (""stopout"" file deleted!)')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
!---------------------------- Closing loops ----------------------------
      if(((total_hw_npt1.ne.1).and.(itype.ne.6)).and.(myrank_row_hw.eq.0)) close(outputunit_loop)
    end do hw_loop
  end do number_of_planes
!----------------------- Deallocating variables ------------------------
  deallocate(r0,c0,r1,c1,r2,c2)
  deallocate(kbz,wkbz,kbz2d)
!----------------------- Finalizing the program ------------------------
  if(myrank.eq.0) call write_time(outputunit,'[main] Finished on: ')
  call MPI_Finalize(ierr)
  if((ierr.ne.0).and.(myrank.eq.0)) write(outputunit,"('[main] Something went wrong in the parallelization! ierr = ',i0)") ierr
  if(myrank.eq.0) close(unit=outputunit)
!=======================================================================
  stop
end program DHE