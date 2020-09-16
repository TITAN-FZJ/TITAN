module mod_initial_expectation

contains

  subroutine calc_initial_Uterms(sys)
    use mod_kind,           only: dp
    use mod_constants,      only: cZero
    use mod_System,         only: System_type,initHamiltkStride
    use TightBinding,       only: initTightBinding
    use mod_magnet,         only: l_matrix,lb,sb,allocate_magnet_variables,deallocate_magnet_variables,rho0,rhod0
    use mod_SOC,            only: ls,allocateLS
    use adaptiveMesh,       only: generateAdaptiveMeshes,genLocalEKMesh,freeLocalEKMesh
    use mod_parameters,     only: nOrb,kp_in,kptotal_in,output,eta,leigenstates,lkpoints,dimH
    use mod_polyBasis,      only: read_basis
    use mod_mpi_pars,       only: myrank,FieldComm,rField,sField,rFreq,sFreq,FreqComm,abortProgram
    use Lattice,            only: initLattice
    use mod_progress,       only: write_time
    use mod_tools,          only: rtos, vec_norm
    use mod_Atom_variables, only: allocate_Atom_variables,deallocate_Atom_variables
    use mod_BrillouinZone,  only: realBZ!,nkpt_x,nkpt_y,nkpt_z
    use EnergyIntegration,  only: pn1
    use mod_hamiltonian,    only: deallocate_hamiltonian
    use mod_superconductivity, only: lsupercond,supercond,allocate_supercond_variables,deallocate_supercond_variables,singlet_coupling
    implicit none
    logical :: lsupercond_temp
    integer :: i,j,mu,err,supercond_temp
    type(System_type), intent(inout) :: sys
    type(System_type), allocatable   :: sys0(:)

    external :: build_U,initConversionMatrices
    if(myrank == 0) &
      call write_time(output%unit,'[calc_initial_Uterms] Obtaining initial densities: ')

    allocate(sys0(sys%nTypes))

    types_of_atoms: do i = 1, sys%nTypes
      !------------------ Define the lattice structure -------------------
      call read_basis(trim(sys%Types(i)%Name), sys0(i))
      if(sys0(i)%nTypes/=1) call abortProgram("[calc_initial_Uterms] Not implemented for parameter file with more than 1 type of atom!")

      !---------------------------- Dimensions -----------------------------
      dimH = sys0(i)%nAtoms*nOrb*2

      !------------- Setting the number of nearest neighbors -------------
      sys0(i)%nStages = sys%nStages
      sys0(i)%relTol  = sys%relTol

      call initLattice(sys0(i))

      ! if(myrank==0) call writeLattice(sys0(i))
      ! stop
      !-------------------------- Filename strings -------------------------
      write(output%info,"('_nkpt=',i0,'_eta=',a)") kptotal_in, trim(rtos(eta,"(es8.1)"))
      if(leigenstates) output%info = trim(output%info) // "_ev"

      !-------------------- Tight Binding parameters ---------------------
      call initTightBinding(sys0(i))

      !---------- Generating k points for real axis integration ----------
      select case(sys0(i)%Types(1)%isysdim)
      case(3)
        sys0(i)%isysdim = 3
        realBZ % nkpt_x = ceiling((dble(kptotal_in))**(0.333333333333333_dp),kind(kp_in(1)) )
        realBZ % nkpt_y = realBZ % nkpt_x
        realBZ % nkpt_z = realBZ % nkpt_x
      case(2)
        sys0(i)%isysdim = 2
        realBZ % nkpt_x = ceiling((dble(kptotal_in))**(0.5_dp),kind(kp_in(1)) )
        realBZ % nkpt_y = realBZ % nkpt_x
        realBZ % nkpt_z = 1
      case default
        sys0(i)%isysdim = 1
        realBZ % nkpt_x = ceiling((dble(kptotal_in)), kind(kp_in(1)) )
        realBZ % nkpt_y = 1
        realBZ % nkpt_z = 1
      end select
      call realBZ % countBZ(sys0(i))

      !---------------- Reading from previous calculations -----------------
      !------- and calculating if file doesn't exist (or different U) ------
      if(.not.read_initial_Uterms(sys0(i),err)) then
        if(myrank == 0) write(output%unit,"('[calc_initial_Uterms] Initial density file for ""',a,'"" does not exist. Calculating...')") trim(sys%Types(i)%Name)

        !----------- Allocating variables that depend on nAtoms ------------
        call allocate_magnet_variables(sys0(i)%nAtoms, nOrb)
        call allocateLS(sys0(i)%nAtoms,nOrb)
        call allocate_supercond_variables(sys0(i)%nAtoms, nOrb)
        call allocate_Atom_variables(sys0(i)%nAtoms,nOrb)

        !------- Initialize Stride Matrices for hamiltk and dtdksub --------
        call initHamiltkStride(sys0(i)%nAtoms, nOrb)

        !---- L matrix in global frame for given quantization direction ----
        call l_matrix()

        !-------------------------- Build U array --------------------------
        call build_U(sys0(i))

        !----- Removing L.B, S.B and L.S matrices from the hamiltonian -----
        lb = cZero
        sb = cZero
        ls = cZero
        singlet_coupling = 0._dp
        ! Turning off superconductivity
        lsuperCond_temp = lsuperCond
        supercond_temp  = supercond
        lsuperCond = .false.
        supercond = 1

        !------------------------ Conversion arrays -------------------------
        call initConversionMatrices(sys0(i)%nAtoms,nOrb)

        ! Distribute Energy Integration across all points available
        call realBZ % setup_fraction(sys0(i),rFreq(1), sFreq(1), FreqComm(1),lkpoints)
        if(.not.leigenstates) then
          !-- Generating k meshes points for imaginary axis integration ----
          call generateAdaptiveMeshes(sys0(i),pn1)
          !--- Distribute Energy Integration across all points available ---
          call genLocalEKMesh(sys0(i),rField,sField, FieldComm)
        end if

        !---------------- Calculating expectation values -------------------
        call calc_expectation_values(sys0(i))

        !-------------- Writing initial densities to files -----------------
        if(myrank == 0) call write_initial_Uterms(sys0(i))

        !------------------------ Freeing memory ---------------------------
        if(.not.leigenstates) call freeLocalEKMesh()

        ! Recovering superconductivity
        lsuperCond = lsuperCond_temp
        supercond  = supercond_temp

        !-------------------- Deallocating variables -----------------------
        call deallocate_magnet_variables()
        call deallocate_Atom_variables()
        call deallocate_supercond_variables()
        call deallocate_hamiltonian()

      end if ! read_initial_Uterms

      !------------------------- Test printing ---------------------------
      if(myrank == 0) then
        write(output%unit,"('[calc_initial_Uterms] Initial densities from: ',a)") trim(sys0(i)%Name)
        do j = 1,sys0(i)%nAtoms
          write(output%unit,"(a,i0,12(2x,es16.9))") trim(sys%Types(i)%Name), j, (sys0(i)%Types(1)%rho0(mu,j),mu=1,nOrb)
          write(output%unit,"(a,i0,12(2x,es16.9))") trim(sys%Types(i)%Name), j, sys0(i)%Types(1)%rhod0(j)
        end do
      end if
    end do types_of_atoms

    ! Transfering from occupations stored on Type to variables used in the hamiltonian
    if(.not.allocated(rho0 )) allocate(rho0(nOrb,sys%nAtoms))
    if(.not.allocated(rhod0)) allocate(rhod0(sys%nAtoms))
    do i=1,sys%nAtoms
      rho0(:,i) = sys0(sys%Basis(i)%Material)%Types(1)%rho0(nOrb,1) ! Unit cell can have more than one atom (of one type)
      rhod0(i)  = sys0(sys%Basis(i)%Material)%Types(1)%rhod0(1)     ! but here only the occupation of first atom is used
    end do

    if(myrank == 0) &
    call write_time(output%unit,'[calc_initial_Uterms] Finished calculating initial density: ')

  end subroutine calc_initial_Uterms


  subroutine calc_expectation_values(sys)
    !! Calculates the expectation values of n_mu^s and n_i/2
    use mod_kind, only: dp
    use mod_constants,     only: cZero
    use mod_System,        only: System_type
    use mod_parameters,    only: nOrb
    use mod_magnet,        only: rho,rhod,mp,mx,my,mz,mpd,mzd
    use mod_expectation,   only: expectation_values
    use mod_Umatrix
    use adaptiveMesh
    implicit none
    integer      :: i
    type(System_type),                   intent(inout) :: sys
    real(dp), dimension(:,:)    , allocatable :: rho0
    real(dp), dimension(:)      , allocatable :: rhod0
    real(dp), dimension(nOrb,sys%nAtoms)      :: deltas

    allocate( rho0(nOrb,sys%nAtoms),rhod0(sys%nAtoms) )

    mzd = 0._dp
    mpd = cZero
    do i = 1, sys%nAtoms
      rhod (i) = sys%Types(sys%Basis(i)%Material)%OccupationD
      rhod0(i) = sys%Types(sys%Basis(i)%Material)%OccupationD
    end do
    rho0 = 0._dp
    rho  = rho0
    call init_Umatrix(mzd,mpd,rhod,rhod0,rho,rho0,sys%nAtoms,nOrb)

    call expectation_values(sys,rho0,mp,mx,my,mz,deltas)

    do i = 1, sys%nAtoms
      rhod0(i)   = sum(rho0(5:9,i))
    end do

    ! The following lines only works for 1 type of atom in the unit cell
    allocate( sys%Types(1)%rho0(nOrb,sys%nAtoms),sys%Types(1)%rhod0(sys%nAtoms) )
    sys%Types(1)%rho0(:,:) = rho0(:,:)
    sys%Types(1)%rhod0(:)  = rhod0(:)

  end subroutine calc_expectation_values

  subroutine write_initial_Uterms(sys0)
    !! Writes the initial orbital dependent densities (calculated with tight-binding hamiltonian only) into files
    use mod_parameters,    only: nOrb, output, dfttype
    use EnergyIntegration, only: parts
    use mod_System,        only: System_type
    implicit none
    character(len=500) :: filename
    character(len=30) :: formatvar
    type(System_type), intent(in) :: sys0
    integer           :: j,mu

    ! Defining and opening file:
    write(output%unit,"('[write_initial_Uterms] Writing initial densities of ""',a,'"" to file:')") trim(sys0%Name)
    write(filename,"('./results/FSOC/selfconsistency/initialrho_',a,'_dfttype=',a,'_parts=',i0,a,a,'.dat')") trim(sys0%Types(1)%Name),dfttype,parts,trim(output%info),trim(output%suffix)
    write(output%unit,"('[write_initial_Uterms] ',a)") trim(filename)
    open (unit=98,status='replace',file=filename)

    ! Writing rho0(nOrb,nAtoms) and rhod0(nAtoms) to file
    write(formatvar,fmt="(a,i0,a)") '(',nOrb*sys0%nAtoms,'(es21.11,2x))'
    write(98,fmt=formatvar) ((sys0%Types(1)%rho0(mu,j),mu=1,nOrb),j=1,sys0%nAtoms)
    write(formatvar,fmt="(a,i0,a)") '(',sys0%nAtoms,'(es21.11,2x))'
    write(98,fmt=formatvar) (sys0%Types(1)%rhod0(j),j=1,sys0%nAtoms)

    ! Writing U to file (to be compared in other calculations)
    write(98,fmt="(es21.11)") sys0%Types(1)%Un, sys0%Types(1)%Um ! Only works for 1 atom in unit cell of elemental file

    close(98)
  end subroutine write_initial_Uterms


  function read_initial_Uterms(sys0,err) result(success)
    !! Writes the initial orbital dependent densities (calculated with tight-binding hamiltonian only) into files
    use mod_kind,          only: dp
    use mod_parameters,    only: nOrb, output, dfttype
    use EnergyIntegration, only: parts
    use mod_System,        only: System_type
    use mod_mpi_pars,      only: rField,MPI_DOUBLE_PRECISION,FieldComm,ierr
    implicit none
    type(System_type), intent(inout) :: sys0
    integer,      intent(out)   :: err
    logical            :: success
    character(len=500) :: filename
    integer            :: j,mu
    real(dp)           :: previous_results_rho0(nOrb,sys0%nAtoms),previous_results_rhod0(sys0%nAtoms),Un_tmp,Um_tmp

    external :: MPI_Bcast

    success = .false.

    write(filename,"('./results/FSOC/selfconsistency/initialrho_',a,'_dfttype=',a,'_parts=',i0,a,a,'.dat')") trim(sys0%Types(1)%Name),dfttype, parts,trim(output%info),trim(output%suffix)
    open(unit=97,file=filename,status="old",iostat=err)
    if(err/=0) then
      if(rField==0) then
        write(output%unit,"('[read_initial_Uterms] Initial density file for ""',a,'"" does not exist:')") trim(sys0%Types(1)%Name)
        write(output%unit,"('[read_initial_Uterms] ',a)") trim(filename)
      end if
      return
    end if

    if(rField==0) then
      write(output%unit,"('[read_initial_Uterms] Initial density file for ""',a,'"" already exists. Reading it from file:')") trim(sys0%Types(1)%Name)
      write(output%unit,"(a)") trim(filename)
    end if

    read(97,fmt=*) ((previous_results_rho0(mu,j),mu=1,nOrb),j=1,sys0%nAtoms)
    read(97,fmt=*) (previous_results_rhod0(j),j=1,sys0%nAtoms)
    read(97,fmt=*) Un_tmp, Um_tmp
    close(97)

    if((abs(Un_tmp - sys0%Types(1)%Un) > 1.e-15_dp).or.(abs(Um_tmp - sys0%Types(1)%Um) > 1.e-15_dp)) then ! Only works for 1 atom in unit cell of elemental file
      write(output%unit,"('[read_initial_Uterms] Different value of Un, Um:')")
      write(output%unit,"('[read_initial_Uterms] Using for ',a,':', es16.9, es16.9,', Read from previous calculations: ', es16.9, es16.9)") trim(sys0%Types(1)%Name), sys0%Types(1)%Un, sys0%Types(1)%Um, Un_tmp, Um_tmp
      write(output%unit,"('[read_initial_Uterms] Recalculating expectation values...')")
      return
    end if

    call MPI_Bcast(previous_results_rho0 ,nOrb*sys0%nAtoms,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
    call MPI_Bcast(previous_results_rhod0,sys0%nAtoms     ,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

    if(.not.allocated(sys0%Types(1)%rho0 )) allocate(sys0%Types(1)%rho0(nOrb,sys0%nAtoms))
    if(.not.allocated(sys0%Types(1)%rhod0)) allocate(sys0%Types(1)%rhod0(sys0%nAtoms))

    do j=1,sys0%nAtoms
      do mu=1,nOrb
        sys0%Types(1)%rho0(mu,j) = previous_results_rho0(mu,j)
      end do
      sys0%Types(1)%rhod0(j) = previous_results_rhod0(j)
    end do

    success = .true.
  end function read_initial_Uterms

end module mod_initial_expectation
