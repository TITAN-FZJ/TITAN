module mod_init_expec
!! Module with subroutines related to initial expectation values calculations

contains

  subroutine calc_init_expec_SK(s)
  !! This subroutine calculates initial expectation values with SK parameters
  !! where each element parameter comes from diferent calculations
    use mod_kind,              only: dp
    use mod_constants,         only: cZero
    use mod_System,            only: System_type,init_Hamiltk_variables,initConversionMatrices
    use TightBinding,          only: initTightBinding
    use mod_magnet,            only: l_matrix,lb,sb,allocate_magnet_variables,deallocate_magnet_variables,mzd0,mpd0,rho0,rhod0
    use mod_SOC,               only: ls,allocateLS
    use adaptiveMesh,          only: generateAdaptiveMeshes,genLocalEKMesh,freeLocalEKMesh
    use mod_parameters,        only: kp_in,kptotal_in,output,eta,leigenstates,lkpoints,dimH,dimHsc
    use mod_polyBasis,         only: read_basis
    use mod_mpi_pars,          only: myrank,FieldComm,rField,sField,rFreq,sFreq,FreqComm,abortProgram
    use Lattice,               only: initLattice
    use mod_progress,          only: write_time
    use mod_tools,             only: rtos
    use mod_Atom_variables,    only: allocate_Atom_variables,deallocate_Atom_variables
    use mod_BrillouinZone,     only: realBZ!,nkpt_x,nkpt_y,nkpt_z
    use EnergyIntegration,     only: pn1
    use mod_hamiltonian,       only: deallocate_hamiltonian
#ifdef _GPU
    use mod_superconductivity, only: lsupercond,supercond,allocate_supercond_variables,deallocate_supercond_variables,delta_sc,delta_sc_d
#else
    use mod_superconductivity, only: lsupercond,supercond,allocate_supercond_variables,deallocate_supercond_variables,delta_sc
#endif
    implicit none
    type(System_type), intent(inout) :: s

    type(System_type), allocatable   :: sys0(:)
    real(dp),    allocatable :: rho0_out(:,:),rhod0_out(:),mzd0_out(:)
    complex(dp), allocatable :: mpd0_out(:)
    logical :: lsupercond_temp,lfound
    integer :: i,j,mu,err,supercond_temp
    character(len=30)  :: formatvar

    if(myrank == 0) &
      call write_time('[calc_init_expec_SK] Obtaining initial densities: ',output%unit)

    allocate(sys0(s%nTypes))

    types_of_atoms: do i = 1, s%nTypes
      ! Checking if calculation is needed
      if(abs(s%Types(i)%Un)<1.e-8) then ! Since we only use paramagnetic SK parameters, only Un is checked here
        if(myrank == 0) &
          write(output%unit,"('[calc_init_expec_SK] Initial density for ""',a,'"" not needed (Un=0.0). Skipping...')") trim(s%Types(i)%Name)

        allocate( s%Types(i)%rho0(s%nOrb) )
        s%Types(i)%rho0(:) = 0._dp
        s%Types(i)%rhod0   = 0._dp
        s%Types(i)%mzd0    = 0._dp
        s%Types(i)%mpd0    = cZero
        cycle
      end if
      !------------------ Define the lattice structure -------------------
      call read_basis(trim(s%Types(i)%Name), sys0(i), .true.)

      !---------------------------- Dimensions -----------------------------
      dimH    = sys0(i)%nAtoms*sys0(i)%nOrb*2
      dimHsc  = dimH

      !------- Setting the subsystem properties from the general one -------
      sys0(i)%nStages = s%nStages
      sys0(i)%relTol  = s%relTol
      sys0(i)%nOrb    = s%nOrb
      sys0(i)%nOrb2   = s%nOrb2
      allocate(sys0(i)%Orbs(s%nOrb))
      sys0(i)%Orbs(:) = s%Orbs(:)
      sys0(i)%nsOrb   = s%nsOrb
      if(sys0(i)%nsOrb>0) then
        allocate(sys0(i)%sOrbs(sys0(i)%nsOrb))
        sys0(i)%sOrbs(:) = s%sOrbs(:)
      end if
      sys0(i)%npOrb   = s%npOrb
      if(sys0(i)%npOrb>0) then
        allocate(sys0(i)%pOrbs(sys0(i)%npOrb))
        sys0(i)%pOrbs(:) = s%pOrbs(:)
      end if
      sys0(i)%ndOrb   = s%ndOrb
      if(sys0(i)%ndOrb>0) then
        allocate(sys0(i)%dOrbs(sys0(i)%ndOrb))
        sys0(i)%dOrbs(:) = s%dOrbs(:)
      end if

      call initLattice(sys0(i))

      ! if(myrank==0) call writeLattice(sys0(i))
      ! stop
      !-------------------------- Filename strings -------------------------
      write(output%info,"('_norb=',i0,'_nkpt=',i0,'_eta=',a)") sys0(i)%nOrb,kptotal_in, trim(rtos(eta,"(es8.1)"))
      if(leigenstates) output%info = trim(output%info) // "_ev"

      !-------------------- Tight Binding parameters ---------------------
      call initTightBinding(sys0(i))

      !---------------- Reading from previous calculations -----------------
      !------- and calculating if file doesn't exist (or different U) ------
      if(.not.read_init_expecs(s%nOrb,s%Types(i),err)) then
        if(myrank == 0) write(output%unit,"('[calc_init_expec_SK] Calculating initial density for ""',a,'""...')") trim(s%Types(i)%Name)

        !----------- Allocating variables that depend on nAtoms ------------
        call allocate_magnet_variables(sys0(i)%nAtoms,sys0(i)%nOrb)
        call allocateLS(sys0(i)%nAtoms,sys0(i)%nOrb)
        call allocate_supercond_variables(sys0(i)%nAtoms,sys0(i)%nOrb)
        call allocate_Atom_variables(sys0(i)%nAtoms,sys0(i)%nOrb)
        allocate(rho0_out(sys0(i)%nOrb,sys0(i)%nAtoms),rhod0_out(sys0(i)%nAtoms),mpd0_out(sys0(i)%nAtoms),mzd0_out(sys0(i)%nAtoms))

        !---------- Generating k points for real axis integration ----------
        select case(sys0(i)%isysdim)
        case(3)
          realBZ % nkpt_x = ceiling((dble(kptotal_in))**(0.333333333333333_dp),kind(kp_in(1)) )
          realBZ % nkpt_y = realBZ % nkpt_x
          realBZ % nkpt_z = realBZ % nkpt_x
        case(2)
          realBZ % nkpt_x = ceiling((dble(kptotal_in))**(0.5_dp),kind(kp_in(1)) )
          realBZ % nkpt_y = realBZ % nkpt_x
          realBZ % nkpt_z = 1
        case default
          realBZ % nkpt_x = ceiling((dble(kptotal_in)), kind(kp_in(1)) )
          realBZ % nkpt_y = 1
          realBZ % nkpt_z = 1
        end select
        call realBZ % countBZ(sys0(i))

        !---- L matrix in global frame for given quantization direction ----
        call l_matrix(sys0(i)%nOrb,sys0(i)%Orbs)

        !----- Removing L.B, S.B and L.S matrices from the hamiltonian -----
        lb = cZero
        sb = cZero
        ls = cZero
        delta_sc = 0._dp
#ifdef _GPU
        delta_sc_d = delta_sc
#endif

        ! Turning off superconductivity
        lsuperCond_temp = lsuperCond
        supercond_temp  = supercond
        lsuperCond = .false.
        supercond = 1

        !------- Initialize Stride Matrices for hamiltk and dtdksub --------
        call init_Hamiltk_variables(sys0(i),supercond)

        !------------------------ Conversion arrays -------------------------
        call initConversionMatrices(sys0(i)%nAtoms,sys0(i)%nOrb)

        ! Distribute Energy Integration across all points available
        call realBZ % setup_fraction(sys0(i),rFreq(1), sFreq(1), FreqComm(1),lkpoints)
        if(.not.leigenstates) then
          !-- Generating k meshes points for imaginary axis integration ----
          call generateAdaptiveMeshes(sys0(i),pn1)
          !--- Distribute Energy Integration across all points available ---
          call genLocalEKMesh(sys0(i),rField,sField, FieldComm)
        end if

        !---------------- Calculating expectation values -------------------
        call calc_init_expec_values(sys0(i),rho0_out,rhod0_out,mpd0_out,mzd0_out)

        ! Storing initial densities
        lfound = .false.
        do j = 1,sys0(i)%nAtoms
          if(trim(sys0(i)%Types(s%Basis(j)%Material)%Name) == trim(s%Types(i)%Name)) then
            allocate( s%Types(i)%rho0(s%nOrb) )
            s%Types(i)%rho0(:) = rho0_out(:,j)
            s%Types(i)%rhod0   = rhod0_out(j)
            s%Types(i)%mzd0    = mzd0_out(j)
            s%Types(i)%mpd0    = mpd0_out(j)
            lfound = .true.
            exit
          end if
        end do
        if(.not.lfound) call abortProgram("[calc_init_expec_SK] Element " // trim(s%Types(i)%Name) // " not found in its elemental file!")

        !-------------- Writing initial densities to files -----------------
        if(myrank == 0) call write_init_expecs(s%nOrb,s%Types(i))

        !------------------------ Freeing memory ---------------------------
        if(.not.leigenstates) call freeLocalEKMesh()

        ! Recovering superconductivity
        lsuperCond = lsuperCond_temp
        supercond  = supercond_temp

        !-------------------- Deallocating variables -----------------------
        deallocate(rho0_out,rhod0_out,mpd0_out,mzd0_out)
        call deallocate_magnet_variables()
        call deallocate_Atom_variables()
        call deallocate_supercond_variables()
        call deallocate_hamiltonian()

      end if ! read_init_expecs

      !------------------------- Test printing ---------------------------
      if(myrank == 0) then
        write(output%unit,"('[calc_init_expec_SK] Initial densities from: ',a)") s%Types(i)%Name
        write(formatvar,fmt="(a,i0,a)") '(a,',s%nOrb,'(es14.7,2x))'
        write(unit=output%unit,fmt=formatvar) trim(s%Types(i)%Name), (s%Types(i)%rho0(mu),mu=1,s%nOrb)
        write(unit=output%unit,fmt="(a,4(2x,es14.7))") trim(s%Types(i)%Name), s%Types(i)%rhod0, s%Types(i)%mpd0%re, s%Types(i)%mpd0%im, s%Types(i)%mzd0
      end if
    end do types_of_atoms

    ! Transfering from occupations stored on Type to variables used in the hamiltonian
    if(.not.allocated(rho0 )) allocate(rho0(s%nOrb,s%nAtoms))
    if(.not.allocated(rhod0)) allocate(rhod0(s%nAtoms))
    if(.not.allocated(mzd0)) allocate(mzd0(s%nAtoms))
    if(.not.allocated(mpd0)) allocate(mpd0(s%nAtoms))
    do i=1,s%nAtoms
      rho0(:,i) = s%Types(s%Basis(i)%Material)%rho0(s%nOrb)
      rhod0(i)  = s%Types(s%Basis(i)%Material)%rhod0
      mzd0(i)   = s%Types(s%Basis(i)%Material)%mzd0
      mpd0(i)   = s%Types(s%Basis(i)%Material)%mpd0
    end do

    if(myrank == 0) &
      call write_time('[calc_init_expec_SK] Finished calculating initial density: ',output%unit)

  end subroutine calc_init_expec_SK


  subroutine calc_init_expec_dft(s)
  !! This subroutine calculates initial expectation values with parameters
  !! parameters from DFT calculations, where a single calculation is done
  !! for all elements (without any extra term added)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero
    use mod_System,            only: System_type
    use mod_magnet,            only: l_matrix,lb,sb,mzd0,mpd0,rho0,rhod0
    use mod_SOC,               only: ls
    use adaptiveMesh,          only: generateAdaptiveMeshes,genLocalEKMesh,freeLocalEKMesh
    use mod_parameters,        only: output,leigenstates,lkpoints
    ! use mod_parameters,        only: kp_in,kptotal_in,output,eta,leigenstates,lkpoints,dimH,dimHsc
    use mod_mpi_pars,          only: myrank,FieldComm,rField,sField,rFreq,sFreq,FreqComm,abortProgram
    use mod_progress,          only: write_time
    use mod_BrillouinZone,     only: realBZ
    use EnergyIntegration,     only: pn1
    use mod_hamiltonian,       only: deallocate_hamiltonian
#ifdef _GPU
    use mod_superconductivity, only: lsupercond,supercond,delta_sc,delta_sc_d
#else
    use mod_superconductivity, only: lsupercond,supercond,delta_sc
#endif
    implicit none
    type(System_type), intent(inout) :: s

    real(dp)    :: rho0_out(s%nOrb,s%nAtoms),rhod0_out(s%nAtoms),mzd0_out(s%nAtoms)
    complex(dp) :: mpd0_out(s%nAtoms)

    logical :: lsupercond_temp,lneed
    integer :: i,mu,err,supercond_temp
    character(len=30)  :: formatvar


    if(myrank == 0) &
      call write_time('[calc_init_expec_dft] Obtaining initial densities: ',output%unit)

    ! Checking if calculation is needed
    lneed = .false.
    types_of_atoms: do i = 1, s%nTypes
      if((abs(s%Types(i)%Un)<1.e-8).and.(abs(s%Types(i)%Um)<1.e-8)) then 
        ! Only not needed when both Un and Um are zero
        allocate( s%Types(i)%rho0(s%nOrb) )
        s%Types(i)%rho0(:) = 0._dp
        s%Types(i)%rhod0   = 0._dp
        s%Types(i)%mzd0    = 0._dp
        s%Types(i)%mpd0    = cZero
        cycle
      else 
        lneed = .true. ! If Un or Um is non zero for a single element, the calculation must be done
        exit           ! Since it's a single calculation, get all the other numbers too
      end if
    end do types_of_atoms

    ! If all values of Un and Um are zero, return
    if(.not.lneed) then
      if(myrank == 0) &
        write(output%unit,"('[calc_init_expec_dft] Initial densities for not needed (Un=Um=0.0). Skipping...')")
      ! Transfering from occupations stored on Type to variables used in the hamiltonian
      if(.not.allocated(rho0 )) allocate(rho0(s%nOrb,s%nAtoms))
      if(.not.allocated(rhod0)) allocate(rhod0(s%nAtoms))
      if(.not.allocated(mzd0)) allocate(mzd0(s%nAtoms))
      if(.not.allocated(mpd0)) allocate(mpd0(s%nAtoms))
      rho0  = 0._dp
      rhod0 = 0._dp
      mzd0  = 0._dp
      mpd0  = cZero
      return
    end if

    !---------------- Reading from previous calculations -----------------
    !------- and calculating if file doesn't exist (or different U) ------
    if(.not.read_init_expecs(s%nOrb,s%Types(i),err)) then
      if(myrank == 0) write(output%unit,"('[calc_init_expec_dft] Calculating initial densities for all elements...')") 

      !---- L matrix in global frame for given quantization direction ----
      call l_matrix(s%nOrb,s%Orbs)

      !----- Removing L.B, S.B and L.S matrices from the hamiltonian -----
      lb = cZero
      sb = cZero
      ls = cZero
      delta_sc = 0._dp
#ifdef _GPU
      delta_sc_d = delta_sc
#endif

      ! Turning off superconductivity
      lsuperCond_temp = lsuperCond
      supercond_temp  = supercond
      lsuperCond = .false.
      supercond = 1

      ! Distribute Energy Integration across all points available
      call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1),lkpoints)
      if(.not.leigenstates) then
        !-- Generating k meshes points for imaginary axis integration ----
        call generateAdaptiveMeshes(s,pn1)
        !--- Distribute Energy Integration across all points available ---
        call genLocalEKMesh(s,rField,sField, FieldComm)
      end if

      !---------------- Calculating expectation values -------------------
      call calc_init_expec_values(s,rho0_out,rhod0_out,mpd0_out,mzd0_out)

      ! Storing initial densities
      do i = 1,s%nAtoms
        if(.not.allocated(s%Types(i)%rho0)) allocate( s%Types(i)%rho0(s%nOrb) )
        s%Types(s%Basis(i)%Material)%rho0(:) = rho0_out(:,i)
        s%Types(s%Basis(i)%Material)%rhod0   = rhod0_out(i)
        s%Types(s%Basis(i)%Material)%mzd0    = mzd0_out(i)
        s%Types(s%Basis(i)%Material)%mpd0    = mpd0_out(i)
      end do

      !-------------- Writing initial densities to files -----------------
      if(myrank == 0) then
        do i = 1,s%nTypes
          call write_init_expecs(s%nOrb,s%Types(i))
        end do
      end if
      !------------------------ Freeing memory ---------------------------
      if(.not.leigenstates) call freeLocalEKMesh()

      ! Recovering superconductivity
      lsuperCond = lsuperCond_temp
      supercond  = supercond_temp
    end if ! read_init_expecs

    !------------------------- Test printing ---------------------------
    if(myrank == 0) then
      do i = 1,s%nTypes
        write(output%unit,"('[calc_init_expec_dft] Initial densities from: ',a)") s%Types(i)%Name
        write(formatvar,fmt="(a,i0,a)") '(a,',s%nOrb,'(es14.7,2x))'
        write(unit=output%unit,fmt=formatvar) trim(s%Types(i)%Name), (s%Types(i)%rho0(mu),mu=1,s%nOrb)
        write(unit=output%unit,fmt="(a,4(2x,es14.7))") trim(s%Types(i)%Name), s%Types(i)%rhod0, s%Types(i)%mpd0%re, s%Types(i)%mpd0%im, s%Types(i)%mzd0
      end do
    end if

    ! Transfering from occupations stored on Type to variables used in the hamiltonian
    if(.not.allocated(rho0 )) allocate(rho0(s%nOrb,s%nAtoms))
    if(.not.allocated(rhod0)) allocate(rhod0(s%nAtoms))
    if(.not.allocated(mzd0)) allocate(mzd0(s%nAtoms))
    if(.not.allocated(mpd0)) allocate(mpd0(s%nAtoms))
    do i=1,s%nAtoms
      rho0(:,i) = s%Types(s%Basis(i)%Material)%rho0(s%nOrb) ! Unit cell can have more than one atom (of one type)
      rhod0(i)  = s%Types(s%Basis(i)%Material)%rhod0
      mzd0(i)   = s%Types(s%Basis(i)%Material)%mzd0     ! but here only the occupation of first atom is used
      mpd0(i)   = s%Types(s%Basis(i)%Material)%mpd0     ! but here only the occupation of first atom is used
    end do

    call deallocate_hamiltonian()

    if(myrank == 0) &
      call write_time('[calc_init_expec_dft] Finished calculating initial density: ',output%unit)

  end subroutine calc_init_expec_dft


  subroutine calc_init_expec_values(s,rho0,rhod0,mpd0,mzd0)
    !! Calculates the expectation values of n_mu^s and n_i/2
    use mod_kind,        only: dp
    use mod_constants,   only: cZero
    use mod_System,      only: System_type
    use mod_expectation, only: expectation_values
    use mod_Umatrix,     only: allocate_Umatrix,hee
    implicit none
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: rho0
    real(dp),    dimension(s%nAtoms),        intent(out) :: rhod0,mzd0
    complex(dp), dimension(s%nAtoms),        intent(out) :: mpd0

    real(dp),    dimension(s%nOrb,s%nAtoms) :: mx,my,mz
    complex(dp), dimension(s%nOrb,s%nAtoms) :: mp
    real(dp),    dimension(s%nOrb,s%nAtoms) :: deltas
    integer      :: i,mu,mud

    call allocate_Umatrix(s%nAtoms,s%nOrb)
    hee = cZero

    call expectation_values(s,rho0,mp,mx,my,mz,deltas)

    rhod0 = 0._dp
    mzd0  = 0._dp
    mpd0  = cZero
    do i=1,s%nAtoms
      do mud=1,s%ndOrb
        mu = s%dOrbs(mud)
        rhod0(i) = rhod0(i) + rho0(mu,i)
        mzd0 (i) = mzd0 (i) + mz  (mu,i)
        mpd0 (i) = mpd0 (i) + cmplx(mx(mu,i),my(mu,i),dp)
      end do
    end do
  end subroutine calc_init_expec_values


  subroutine write_init_expecs(nOrb,material)
    !! Writes the initial orbital dependent densities (calculated with tight-binding hamiltonian only) into files
    use mod_parameters,    only: output,dfttype
    use EnergyIntegration, only: parts
    use AtomTypes,         only: AtomType
    implicit none
    character(len=500) :: filename
    character(len=30)  :: formatvar
    integer,        intent(in) :: nOrb
    type(AtomType), intent(in) :: material
    ! type(System_type), intent(in) :: sys0
    integer           :: mu,funit=98

    ! Defining and opening file:
    write(output%unit,"('[write_init_expecs] Writing initial densities of ""',a,'"" to file:')") trim(material%Name)
    write(filename,"('./results/FSOC/selfconsistency/initialrho_',a,'_dfttype=',a,'_parts=',i0,a,a,'.dat')") trim(material%Name),dfttype,parts,trim(output%info),trim(output%suffix)
    write(output%unit,"('[write_init_expecs] ',a)") trim(filename)
    open (unit=funit,status='replace',file=filename)

    ! Writing rho0(nOrb,nAtoms) and rhod0(nAtoms) to file
    write(formatvar,fmt="(a,i0,a)") '(',nOrb,'(es21.11,2x))'
    write(unit=funit,fmt=formatvar) (material%rho0(mu),mu=1,nOrb)
    write(unit=funit,fmt="(4(es21.11,2x))") material%rhod0,material%mpd0%re,material%mpd0%im,material%mzd0

    ! Writing U to file (to be compared in other calculations)
    write(unit=funit,fmt="(2(es21.11,2x))") material%Un, material%Um 

    close(unit=funit)
  end subroutine write_init_expecs


  function read_init_expecs(nOrb,material,err) result(success)
    !! Writes the initial orbital dependent densities (calculated with tight-binding hamiltonian only) into files
    use mod_kind,          only: dp
    use mod_parameters,    only: output,dfttype
    use EnergyIntegration, only: parts
    use AtomTypes,         only: AtomType
    use mod_mpi_pars,      only: rField,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,FieldComm,ierr
    implicit none
    integer,        intent(in)    :: nOrb
    type(AtomType), intent(inout) :: material
    integer,        intent(out)   :: err
    logical            :: success
    character(len=500) :: filename
    integer            :: mu,funit=97
    real(dp)           :: previous_rho0(nOrb),previous_rhod0,previous_mpd0(2),previous_mzd0,Un_tmp,Um_tmp

    external :: MPI_Bcast

    success = .false.

    write(filename,"('./results/FSOC/selfconsistency/initialrho_',a,'_dfttype=',a,'_parts=',i0,a,a,'.dat')") trim(material%Name),dfttype, parts,trim(output%info),trim(output%suffix)
    open(unit=funit,file=filename,status="old",iostat=err)
    if(err/=0) then
      if(rField==0) then
        write(output%unit,"('[read_init_expecs] Initial density file for ""',a,'"" does not exist:')") trim(material%Name)
        write(output%unit,"('[read_init_expecs] ',a)") trim(filename)
      end if
      return
    end if

    if(rField==0) then
      write(output%unit,"('[read_init_expecs] Initial density file for ""',a,'"" already exists. Reading it from file:')") trim(material%Name)
      write(output%unit,"(a)") trim(filename)
    end if

    read( unit=funit,fmt=*) (previous_rho0(mu),mu=1,nOrb)
    read( unit=funit,fmt=*) previous_rhod0,previous_mpd0(1),previous_mpd0(2),previous_mzd0
    read( unit=funit,fmt=*) Un_tmp, Um_tmp
    close(unit=funit)

    if((abs(Un_tmp - material%Un) > 1.e-15_dp).or.(abs(Um_tmp - material%Um) > 1.e-15_dp)) then
      if(rField==0) then
        write(output%unit,"('[read_init_expecs] Different value of Un, Um:')")
        write(output%unit,"('[read_init_expecs] Using for ',a,':', es16.9, es16.9,', Read from previous calculations: ', es16.9, es16.9)") trim(material%Name), material%Un, material%Um, Un_tmp, Um_tmp
        write(output%unit,"('[read_init_expecs] Recalculating expectation values...')")
      end if
      return
    end if

    call MPI_Bcast(previous_rho0 , nOrb  ,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
    call MPI_Bcast(previous_rhod0, 1     ,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
    call MPI_Bcast(previous_mpd0 , 1     ,MPI_DOUBLE_COMPLEX  ,0,FieldComm,ierr)
    call MPI_Bcast(previous_mzd0 , 1     ,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

    if(.not.allocated(material%rho0)) allocate(material%rho0(nOrb))
    material%rho0(:) = previous_rho0(:)
    material%rhod0   = previous_rhod0
    material%mzd0    = previous_mzd0
    material%mpd0    = cmplx(previous_mpd0(1),previous_mpd0(2),dp)

    success = .true.
  end function read_init_expecs

end module mod_init_expec
