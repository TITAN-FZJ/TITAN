module TightBinding
  implicit none

contains

  subroutine initTightBinding(s)
    use mod_System,     only: System_type
    use mod_mpi_pars,   only: abortProgram
    use mod_parameters, only: tbmode,fermi_layer
    implicit none
    type(System_type), intent(inout) :: s
    integer :: i

    ! Getting information from all the elements of "basis" in the element files
    do i = 1, s%nTypes
      call readElementFile(s,i,tbmode)
    end do

    s%total_nOrb = s%Types(s%Basis(1)%Material)%nOrb
    do i = 2, s%nAtoms
      s%total_nOrb = s%total_nOrb + s%Types(s%Basis(i)%Material)%nOrb
    end do

    select case(tbmode)
    case(1)
      call get_SK_parameter(s,fermi_layer)
    case(2)
      call get_DFT_parameters(s,fermi_layer)
    case default
      call abortProgram("[initTightBinding] Only tbmode = 1 (SK parameters) or = 2 (parameters from DFT)")
    end select

  end subroutine initTightBinding


  subroutine get_DFT_parameters(s,fermi_layer)
  !! Gets hamiltonian from DFT input and parameters from elemental file
    use mod_system, only: System_type
    use mod_dft,    only: readHamiltonian
    implicit none
    type(System_type),   intent(inout) :: s
    integer,             intent(in)    :: fermi_layer
    integer                            :: i

    do i = 1, s%nAtoms
      ! Getting total number of electrons on the system
      s%totalOccupation = s%totalOccupation + s%Types(s%Basis(i)%Material)%Occupation
      ! Storing Hubbard e-e interaction per basis atom
      s%Basis(i)%Un = s%Types(s%Basis(i)%Material)%Un
      s%Basis(i)%Um = s%Types(s%Basis(i)%Material)%Um
    end do

    ! Read hamiltonian and dft parameters
    call readHamiltonian(s,s%Name)

    ! Setting Fermi energy
    s%Ef = s%Types(fermi_layer)%FermiLevel
  end subroutine get_DFT_parameters


  subroutine get_SK_parameter(s,fermi_layer)
  !! Gets parameters for SK two-center approximation and build hopping matrices
    use mod_kind,      only: dp
    use AtomTypes,     only: NeighborIndex
    use mod_system,    only: System_type
    use mod_constants, only: cZero
    implicit none
    type(System_type),   intent(inout) :: s
    integer,             intent(in)    :: fermi_layer
    integer                            :: i,j,k,l
    real(dp)                           :: scale_factor(2), mix(10,2), temp
    type(NeighborIndex), pointer       :: current
    real(dp), dimension(10), parameter :: expon = [1.0_dp,3.0_dp,3.0_dp,5.0_dp,5.0_dp,5.0_dp,2.0_dp,3.0_dp,4.0_dp,4.0_dp]

    nullify(current)

    ! Shifting all energies to a common Fermi level
    s%Ef = s%Types(fermi_layer)%FermiLevel
    do i = 1, s%nTypes
      do j = 1, s%Types(i)%nOrb
        s%Types(i)%onSite(j,j) = s%Types(i)%onSite(j,j) - s%Types(i)%FermiLevel + s%Ef
      end do
    end do

    ! Allocate & initialize Hopping variables
    do i = 1, s%nNeighbors
      allocate(s%Neighbors(i)%isHopping(s%nAtoms))
      s%Neighbors(i)%t0i = cZero
      s%Neighbors(i)%isHopping = .false.
    end do

    ! Setting up hopping matrices between all neighbors of the elements in basis file
    ! A mix (defaul: geometric; optional: simple average) between different elements is used
    !
    ! Scaling law by Andersen et al. O.K. Andersen, O. Jepsen, Physica 91B, 317 (1977); O.K. Andersen, W. Close. H. Nohl, Phys. Rev. B17, 1209 (1978)
    ! Distance dependence of tight binding matrix elements is given by V = C * d^(-[l+l'+1])
    ! e.g for ss hopping distance dependence is d^-1, for sp hopping d^-2
    ! do mu = 1, s%nOrb
    !   do nu = 1, s%nOrb
    !     s%Neighbors(current%index)%t0i(mu,nu,i) = s%Neighbors(current%index)%t0i(mu,nu,i) * scale_factor ** (mu + nu + 1)
    !   end do
    ! end do

    ! Looping over atoms in the unit cell
    do i = 1, s%nAtoms
      ! Getting total number of electrons on the system
      s%totalOccupation = s%totalOccupation + s%Types(s%Basis(i)%Material)%Occupation
      ! Storing Hubbard e-e interaction per basis atom
      s%Basis(i)%Un = s%Types(s%Basis(i)%Material)%Un
      s%Basis(i)%Um = s%Types(s%Basis(i)%Material)%Um
      ! Looping over neighbor stages and atoms in their unit cell and setting hopping parameters
      do j = 1, s%nAtoms
        do k = 1, s%nStages
          current => s%Basis(i)%NeighborList(k,j)%head
          do while(associated(current))
            ! Getting the smallest scale factor from the atoms in element file
            scale_factor = 0._dp
            do l=1,s%Types(s%Basis(i)%Material)%nAtoms
              if(.not.s%Types(s%Basis(i)%Material)%lelement(l)) cycle
              temp = s%Types(s%Basis(i)%Material)%stage(k,l) / s%Neighbors(current%index)%Distance(i)
              if(abs(temp-1._dp) < abs(scale_factor(1)-1._dp) ) scale_factor(1) = temp
            end do
            do l=1,s%Types(s%Basis(j)%Material)%nAtoms
              if(.not.s%Types(s%Basis(j)%Material)%lelement(l)) cycle
              temp = s%Types(s%Basis(j)%Material)%stage(k,l) / s%Neighbors(current%index)%Distance(i)
              if(abs(temp-1._dp) < abs(scale_factor(2)-1._dp) ) scale_factor(2) = temp
            end do

            do l = 1, 10
              mix(l,1) = s%Types(s%Basis(i)%Material)%Hopping(l,k) * scale_factor(1) ** expon(l)
              mix(l,2) = s%Types(s%Basis(j)%Material)%Hopping(l,k) * scale_factor(2) ** expon(l)
            end do
            if(.not.allocated(s%Neighbors(current%index)%t0i)) &
              allocate(s%Neighbors(current%index)%t0i(s%Types(s%Basis(j)%Material)%nOrb,s%Types(s%Basis(i)%Material)%nOrb,s%nAtoms))

            call set_hopping_matrix(s%Neighbors(current%index)%dirCos(:,i), &
                                    mix(:,1),mix(:,2),s%Types(s%Basis(i)%Material)%nOrb,s%Types(s%Basis(j)%Material)%nOrb, &
                                    s%Types(s%Basis(i)%Material)%Orbs,s%Types(s%Basis(j)%Material)%Orbs, &
                                    s%Neighbors(current%index)%t0i(:,:,i))
            s%Neighbors(current%index)%isHopping(i) = .true.
            current => current%next
          end do
        end do
      end do
    end do
  end subroutine get_SK_parameter

  subroutine set_hopping_matrix(dirCos,t1,t2,nOrb_i,nOrb_j,Orbs_i,Orbs_j,t0i)
    use mod_kind,       only: dp
    use mod_parameters, only: lsimplemix
    use AtomTypes,      only: default_orbitals
    implicit none
    real(dp),    dimension(3),             intent(in)    :: dirCos
    real(dp),    dimension(10),            intent(in)    :: t1, t2
    integer,                               intent(in)    :: nOrb_j,nOrb_i
    integer,     dimension(nOrb_j),        intent(in)    :: Orbs_j
    integer,     dimension(nOrb_i),        intent(in)    :: Orbs_i
    complex(dp), dimension(nOrb_j,nOrb_i), intent(inout) :: t0i
    integer                 :: i,j
    real(dp), dimension(10) :: mix
    real(dp), dimension(size(default_orbitals),size(default_orbitals)) :: bp

    do i = 1, 10
      if(lsimplemix) then
        mix(i) = 0.5_dp*(t1(i) + t2(i))
      else
        mix(i) = sign(sqrt(abs(t1(i)) * abs(t2(i))), t1(i) + t2(i))
      end if
    end do
    call intd(mix(1),mix(2),mix(3),mix(4),mix(5),mix(6),mix(7),mix(8),mix(9),mix(10),dirCos,bp)

    ! Selecting orbitals
    do j=1,nOrb_j
      do i=1,nOrb_i
        t0i(i,j) = cmplx(bp(Orbs_i(i),Orbs_j(j)),0.d0,dp)
      end do
    end do

  end subroutine set_hopping_matrix

  subroutine readElementFile(s,n,tbmode)
  !! Reading element file, including all the parameters
    use mod_kind,              only: dp
    use AtomTypes,             only: NeighborAtom,default_orbitals
    use mod_system,            only: System_type
    use mod_mpi_pars,          only: abortProgram
    use mod_tools,             only: ItoS,StoI,StoR,StoArray,next_line,vec_norm,vecDist,is_numeric
    use mod_io,                only: log_warning,log_error,log_message,get_orbitals
    use mod_input,             only: get_parameter
    use mod_superconductivity, only: lsuperCond
    use mod_SOC,               only: socscale
    implicit none
    type(System_type), intent(inout) :: s
    integer,           intent(in)    :: n
    integer,           intent(in)    :: tbmode
    integer :: nTypes
    integer :: nn_stages
    integer,  dimension(:), allocatable :: iMaterial
    integer,  dimension(:), allocatable :: type_count
    real(dp), dimension(3) :: dens
    integer            :: i,j,k,l,mu,nat,nCells,ios,line_count = 0
    integer, parameter :: word_length = 20, max_elements = 50
    character(len=word_length), dimension(max_elements) :: str_arr
    real(dp),                   dimension(max_elements) :: tmp_arr
    character(200) :: line
    character(50)  :: words(10)
    real(dp), dimension(3,3)  :: Bravais
    real(dp), dimension(size(default_orbitals)) :: on_site
    real(dp), dimension(:,:), allocatable :: position
    real(dp), dimension(:,:), allocatable :: localDistances
    type(NeighborAtom), dimension(:), allocatable :: list
    character(len=1)   :: coord_type
    character(len=100) :: Name
    integer :: f_unit = 995594, cnt = 0

    ! Opening file
    open(f_unit, file=trim(s%Types(n)%Name), status='old', iostat=ios)
    if(ios /= 0) call abortProgram("[readElementFile] Error occured when trying to read file " // trim(s%Types(n)%Name))
    line_count = 0

    ! Read parameter name
    read(f_unit, fmt='(A)', iostat=ios) Name
    Name = trim(adjustl(Name))

    if(tbmode==1) then
      ! Reading lattice parameter a0
      s%Types(n)%LatticeConstant = StoR(next_line("readElementFile",f_unit,"lattice parameter"))

      ! Bravais lattice
      do j = 1, 3
        Bravais(1:3,j) = StoR( next_line("readElementFile",f_unit,"Bravais lattice") ,3) 
      end do
      Bravais = Bravais * s%Types(n)%LatticeConstant
      s%Types(n)%a1 = Bravais(:,1)
      s%Types(n)%a2 = Bravais(:,2)
      s%Types(n)%a3 = Bravais(:,3)

      ! Read Different Elements in File
      str_arr(1:max_elements) = StoArray(next_line("readElementFile",f_unit,"different elements"),max_elements)
      nTypes = 0
      do i = 1, max_elements
        if(str_arr(i)(1:1) == "!") exit
        if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
        nTypes = nTypes + 1
      end do

      ! Reading name of elements
      allocate(s%Types(n)%Types(nTypes))
      j = 1
      do i = 1, max_elements
        if(str_arr(i)(1:1) == "!") exit
        if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
        s%Types(n)%Types(j) = str_arr(i)
        if(trim(s%Types(n)%Types(j)) == trim(s%Types(n)%Name)) s%Types(n)%nAtomEl = s%Types(n)%nAtomEl + 1

        j = j + 1
      end do

      select case(s%Types(n)%nAtomEl)
      case(0)
        call log_error("readElementFile", "Element " // trim(s%Types(n)%Name) // " not found in elemental file!")
      case(2:)
        call log_warning("readElementFile", "More than one instance of element " // trim(itos(i)) // " found in elemental file. Using the scale factor closest to 1.")
      end select

      ! Read amount of atoms per element
      allocate(type_count(nTypes))
      type_count(1:nTypes) = StoI(next_line("readElementFile",f_unit,"atoms per element"),nTypes)

      ! Count number of atoms
      s%Types(n)%nAtoms = sum(type_count(1:nTypes))
      if(s%Types(n)%nAtoms <= 0) call abortProgram("[readElementFile] No basis atoms given!")

      ! Read coordinate type
      str_arr(1:1) = StoArray(next_line("readElementFile",f_unit,"coordinate type"),1)
      coord_type = trim(str_arr(1))

      allocate(imaterial(s%Types(n)%nAtoms),s%Types(n)%lelement(s%Types(n)%nAtoms))
      s%Types(n)%lelement = .false.

      ! Read atom positions
      allocate(position(3,s%Types(n)%nAtoms))

      k = 0
      do i = 1, nTypes
        do j = 1, type_count(i)
          k = k + 1
          imaterial(k) = i
          if(trim(s%Types(n)%Types(i)) == trim(s%Types(n)%Name)) s%Types(n)%lelement(k) = .true.
          position(1:3,k) = StoR( next_line("readElementFile",f_unit,"basis atoms") ,3) 
          if(coord_type == 'C' .or. coord_type == 'c' .or. coord_type == 'K' .or. coord_type == 'k') then
            ! Position of atoms given in Cartesian coordinates
            position(:,k) = position(:,k) * s%Types(n)%LatticeConstant
          else
            ! Position of atoms given in Bravais (or Direct, Internal, Lattice) coordinates
            position(:,k) = position(1,k) * s%Types(n)%a1 + position(2,k) * s%Types(n)%a2 + position(3,k) * s%Types(n)%a3
          end if
        end do
      end do

      ! Read dimension of the system
      s%Types(n)%isysdim = StoI(next_line("readElementFile",f_unit,"system dimension"))

      do i = 1, s%Types(n)%isysdim
        if(vec_norm(Bravais(:,i),3) <= 1.e-9_dp) &
          call log_error("readElementFile", "Bravais vector a" // trim(itos(i)) // " not given.")
      end do
    end if

    ! Getting number of orbitals and Fermi energy
    line = next_line("readElementFile",f_unit,"orbitals or Fermi energy")

    ! Check if orbitals are given in elemental file
    orbital_selection: if(is_numeric( line(1:1) )) then
      ! If not, pass general selection of orbitals to each type
      s%Types(n)%nOrb  = s%nOrb
      allocate(s%Types(n)%Orbs(s%Types(n)%nOrb))
      s%Types(n)%Orbs(1:s%nOrb) = s%Orbs(1:s%nOrb)
      
      if(s%nsOrb>0) then
        s%Types(n)%nsOrb = s%nsOrb
        allocate(s%Types(n)%sOrbs(s%nsOrb))
        s%Types(n)%sOrbs(1:s%nsOrb) = s%sOrbs(1:s%nsOrb)
      end if

      if(s%npOrb>0) then
        s%Types(n)%npOrb = s%npOrb
        allocate(s%Types(n)%pOrbs(s%npOrb))
        s%Types(n)%pOrbs(1:s%npOrb) = s%pOrbs(1:s%npOrb)
      end if

      if(s%ndOrb>0) then
        s%Types(n)%ndOrb = s%ndOrb
        allocate(s%Types(n)%dOrbs(s%ndOrb))
        s%Types(n)%dOrbs(1:s%ndOrb) = s%dOrbs(1:s%ndOrb)
      end if

      ! call log_message("readElementFile",trim(itos(s%nOrb)) // " orbitals selected for " // trim(s%Types(n)%Name) //":" // trim(selected_orbitals))

      ! Read Fermi level
      s%Types(n)%FermiLevel = StoR(line)
    else orbital_selection
      ! Selecting orbitals for this type
      str_arr = ""
      str_arr(1:max_elements) = StoArray(line,max_elements)
      call get_orbitals(str_arr(2:),s%Types(n)%nOrb,s%Types(n)%nsOrb,s%Types(n)%npOrb,s%Types(n)%ndOrb,s%Types(n)%Orbs,s%Types(n)%sOrbs,s%Types(n)%pOrbs,s%Types(n)%dOrbs)
      ! selected_orbitals = ""
      ! selected_sorbitals = ""
      ! selected_porbitals = ""
      ! selected_dorbitals = ""
      ! s%nsOrb = 0
      ! s%npOrb = 0
      ! s%ndOrb = 0

      ! s%Types(n)%nOrb = 0
      ! str_arr = ""
      ! str_arr(1:max_elements) = StoArray(line,max_elements)
      ! do i = 2,max_elements
      !   if(str_arr(i)(1:1) == "!") exit
      !   if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle

      !   s%Types(n)%nOrb = s%Types(n)%nOrb + 1
      !   if(.not.is_numeric( trim(str_arr(i)) )) then
      !     iloc = findloc( default_orbitals,str_arr(i)(1:3),dim=1 )
      !     if(iloc == 0) &
      !       call log_error("readElementFile","Orbital not recognized: " // str_arr(i)(1:3) //". Use one of the following: " // NEW_line('A') // &
      !           "(1|s), (2|px), (3|py), (4|pz), (5|dxy), (6|dyz), (7|dzx), (8|dx2), (9|dz2)")
      !     str_arr(i) = itos( iloc )
      !   end if
      ! end do
      ! allocate(s%Types(n)%Orbs(s%Types(n)%nOrb),itmps_arr(s%Types(n)%nOrb),itmpp_arr(s%Types(n)%nOrb),itmpd_arr(s%Types(n)%nOrb))


      ! ! Looping over selected orbitals
      ! ! Transforms names to numbers
      ! ! and count d orbitals
      ! do i = 1,s%nOrb
      !   ! If the name of the orbital is given instead of a number, convert:
      !   if(.not.is_numeric( trim(s_vector(i)) )) then
      !     iloc = findloc( default_orbitals,s_vector(i)(1:3),dim=1 )
      !     if(iloc == 0) &
      !       call log_error("get_parameters","Orbital not recognized: " // s_vector(i)(1:3) //". Use one of the following: " // NEW_line('A') // &
      !           "(1|s), (2|px), (3|py), (4|pz), (5|dxy), (6|dyz), (7|dzx), (8|dx2), (9|dz2)")
      !     s_vector(i) = itos( iloc )
      !   end if
      !   s%Orbs(i) = stoi( trim(s_vector(i)) )
      !   selected_orbitals = trim(selected_orbitals)  // " " // trim(default_orbitals(s%Orbs(i)))
      !   ! Checking orbital type, and storing information
      !   if(s%Orbs(i)==1) then
      !     s%nsOrb = s%nsOrb + 1
      !     itmps_arr(s%nsOrb) = i
      !     selected_sorbitals = trim(selected_sorbitals)  // " " // trim(default_orbitals(s%Orbs(i)))
      !   end if
      !   if((s%Orbs(i)>=2).and.(s%Orbs(i)<=4)) then
      !     s%npOrb = s%npOrb + 1
      !     itmpp_arr(s%npOrb) = i
      !     selected_porbitals = trim(selected_porbitals)  // " " // trim(default_orbitals(s%Orbs(i)))
      !   end if
      !   if((s%Orbs(i)>=5).and.(s%Orbs(i)<=9)) then
      !     s%ndOrb = s%ndOrb + 1
      !     itmpd_arr(s%ndOrb) = i
      !     selected_dorbitals = trim(selected_dorbitals)  // " " // trim(default_orbitals(s%Orbs(i)))
      !   end if
      ! end do
      ! call log_message("get_parameters",trim(itos(s%nOrb)) // " orbitals selected:" // trim(selected_orbitals) // ", of which:")
      ! deallocate(s_vector)
      ! ! s-orbitals
      ! if(s%nsOrb > 0) then
      !   allocate(s%sOrbs(s%nsOrb))
      !   s%sOrbs(1:s%nsOrb) = itmps_arr(1:s%nsOrb)
      !   call log_message("get_parameters", trim(itos(s%nsOrb)) // " s orbitals:" // trim(selected_sorbitals) )
      ! end if
      ! ! p-orbitals
      ! if(s%npOrb > 0) then
      !   allocate(s%pOrbs(s%npOrb))
      !   s%pOrbs(1:s%npOrb) = itmpp_arr(1:s%npOrb)
      !   call log_message("get_parameters", trim(itos(s%npOrb)) // " p orbitals:" // trim(selected_porbitals) )
      ! end if
      ! ! d-orbitals
      ! if(s%ndOrb > 0) then
      !   allocate(s%dOrbs(s%ndOrb))
      !   s%dOrbs(1:s%ndOrb) = itmpd_arr(1:s%ndOrb)
      !   call log_message("get_parameters", trim(itos(s%ndOrb)) // " d orbitals:" // trim(selected_dorbitals) )
      ! end if













      ! s%Types(n)%nOrb = 0
      ! str_arr = ""
      ! str_arr(1:max_elements) = StoArray(line,max_elements)
      ! do i = 2,max_elements
      !   if(str_arr(i)(1:1) == "!") exit
      !   if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle

      !   s%Types(n)%nOrb = s%Types(n)%nOrb + 1
      !   if(.not.is_numeric( trim(str_arr(i)) )) then
      !     if(findloc( orbitals,trim(str_arr(i)),dim=1 ) == 0) &
      !       call log_error("readElementFile","Orbital not recognized: " // trim(str_arr(i)) //". Use one of the following: " // NEW_line('A') // &
      !            "(1 or s), (2|px), (3|py), (4|pz), (5|dxy), (6|dyz), (7|dzx), (8|dx2), (9|dz2)")
      !     str_arr(i) = itos( findloc(orbitals,trim(str_arr(i)),dim=1) )
      !   end if
      ! end do
      ! ! Selected orbitals
      ! allocate(s%Types(n)%Orbs(s%Types(n)%nOrb))
      ! selected_orbitals = ""
      ! selected_dorbitals = ""
      ! s%Types(n)%ndOrbs = 0
      ! do i=1,s%Types(n)%nOrb
      !   s%Types(n)%Orbs(i) = stoi( trim(str_arr(i+1)) )
      !   selected_orbitals = trim(selected_orbitals)  // " " // trim(orbitals(s%Types(n)%Orbs(i)))
      !   ! if it's d orbital:
      !   if((s%Types(n)%Orbs(i)>=5).and.(s%Types(n)%Orbs(i)<=9)) then
      !     s%Types(n)%ndOrbs = s%Types(n)%ndOrbs + 1
      !     itmp_arr(s%Types(n)%ndOrbs) = i
      !     selected_dorbitals = trim(selected_dorbitals)  // " " // trim(orbitals(s%Types(n)%Orbs(i)))
      !   end if
      ! end do
      ! call log_message("readElementFile",trim(itos(s%Types(n)%nOrb)) // " orbitals selected for " // trim(s%Types(n)%Name) // ":" // trim(selected_orbitals) )

      ! ! d-orbitals
      ! if(s%Types(n)%ndOrbs > 0) then
      !   allocate(s%Types(n)%dOrbs(s%Types(n)%ndOrbs))
      !   s%Types(n)%dOrbs(1:s%Types(n)%ndOrbs) = itmp_arr(1:s%Types(n)%ndOrbs)
      !   call log_message("readElementFile","of which " // trim(itos(s%Types(n)%nOrb)) // "d orbitals:" // trim(selected_dorbitals) )
      ! end if

      ! Read Fermi level
      s%Types(n)%FermiLevel = StoR(next_line("readElementFile",f_unit,"Fermi energy"))
    end if orbital_selection








    





    ! Read charge densitites for s p d
    ! line = next_line("readElementFile",f_unit)
    dens(1:3) = StoR(next_line("readElementFile",f_unit,"occupations"),3)
    s%Types(n)%OccupationS = dens(1)
    s%Types(n)%OccupationP = dens(2)
    s%Types(n)%OccupationD = dens(3)
    s%Types(n)%Occupation  = s%Types(n)%OccupationS+s%Types(n)%OccupationP+s%Types(n)%OccupationD

    ! Read Hubbard effective Coulomb interaction strength
    str_arr(1:2) = StoArray(next_line("readElementFile",f_unit,"Hubbard U"),2)
    cnt = 0
    do i = 1,2
      if(str_arr(i)(1:1) == "!") exit
      if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
      cnt = cnt + 1
      read(unit=str_arr(i), fmt=*,iostat=ios ) tmp_arr(cnt)
    end do    
    select case(cnt)
    case(1)
      s%Types(n)%Un = tmp_arr(1)
      s%Types(n)%Um = tmp_arr(1)
    case(2)
      s%Types(n)%Un = tmp_arr(1)
      s%Types(n)%Um = tmp_arr(2)
    case default
      call log_error("readElementFile","Something wrong in the definition of 'U'.")
    end select

    ! Read Spin-Orbit interaction strength for p and d
    tmp_arr(1:2) = StoR(next_line("readElementFile",f_unit,"SOI strength"),2)
    s%Types(n)%LambdaP = socscale*tmp_arr(1)
    s%Types(n)%LambdaD = socscale*tmp_arr(2)

    ! Read the superconducting parameter
    str_arr(1:max_elements) = StoArray(next_line("readElementFile",f_unit,"superconducting parameter"),max_elements)
    cnt = 0
    do i = 1, max_elements
      if(str_arr(i)(1:1) == "!") exit
      if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
      cnt = cnt + 1
      read(unit=str_arr(i), fmt=*,iostat=ios ) tmp_arr(cnt)
    end do

    allocate(s%Types(n)%lambda(s%Types(n)%nOrb))
    ! Sigle value for all orbitals
    if (cnt==1) then
      s%Types(n)%lambda(1:s%Types(n)%nOrb) = tmp_arr(1)
    ! One value per orbital
    else if (cnt==s%Types(n)%nOrb) then
      s%Types(n)%lambda(1:s%Types(n)%nOrb) = tmp_arr(1)
    ! One value per general orbital type (s,p,d)
    else if (cnt==3) then
      do mu=1,s%Types(n)%nOrb
        ! s-type orbital
        if(s%Types(n)%Orbs(mu)==1) then
          s%Types(n)%lambda(mu) = tmp_arr(1)
        end if
        ! p-type orbital
        if((s%Types(n)%Orbs(mu)>=2).and.(s%Types(n)%Orbs(mu)<=4)) then
          s%Types(n)%lambda(mu) = tmp_arr(2)
        end if
        ! d-type orbital
        if((s%Types(n)%Orbs(mu)>=5).and.(s%Types(n)%Orbs(mu)<=9)) then
          s%Types(n)%lambda(mu) = tmp_arr(3)
        end if
      end do
    ! One value per specific orbital type (s,px,py,pz,dxy,dyz,dzx,dx2,dz2)
    else if (cnt==9) then
      do mu=1,s%Types(n)%nOrb
        s%Types(n)%lambda(mu) = tmp_arr(s%Types(n)%Orbs(mu))
      end do
    else
      if(lsupercond) &
        call log_error("readElementFile","Something wrong in the definition of 'lambda'.")
    end if      

    ! Only continue for SK parameters
    if(tbmode/=1) then
      close(f_unit)
      return
    end if

    ! Read nearest neighbor stages
    nn_stages = StoI(next_line("readElementFile",f_unit,"nearest neighbor stages"))

    ! Read Hopping Parameter

    ! Read on-site terms
    ! s
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit=line, fmt=*, iostat=ios) (words(k), k=1,10)
    read(unit=words(3), fmt=*, iostat=ios) on_site(1)
    if(ios /= 0) &
      call log_error("readElementFile","Something wrong in the s tight-binding parameter of element " // trim(Name) //".")

    ! p
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit=line, fmt=*, iostat=ios) (words(k), k=1,10)
    if(words(4) == "2") then
      read(unit=words(3), fmt=*, iostat=ios) on_site(2)
      if(ios /= 0)  &
        call log_error("readElementFile","Something wrong in the p tight-binding parameters of element " // trim(Name) //".")
      on_site(3:4) = on_site(2)
    else
      do j=3,5
        read(unit=words(j), fmt=*, iostat=ios) on_site(j-1)
        if(ios /= 0)  &
          call log_error("readElementFile","Something wrong in the p tight-binding parameters of element " // trim(Name) //".")
      end do
    end if

    ! t2g
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit=line, fmt=*, iostat=ios) (words(k), k=1,10)
    if(words(4) == "3") then
      read(unit=words(3), fmt=*, iostat=ios) on_site(5)
      on_site(6:7) = on_site(5)
      if(ios /= 0)  &
        call log_error("readElementFile","Something wrong in the t2g tight-binding parameters of element " // trim(Name) //".")
    else
      do j=3,5
        read(unit=words(j), fmt=*, iostat=ios) on_site(j+2)
        if(ios /= 0)  &
          call log_error("readElementFile","Something wrong in the t2g tight-binding parameters of element " // trim(Name) //".")
      end do
    end if

    ! eg
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit=line, fmt=*, iostat=ios) (words(k), k=1,10)
    if(words(4) == "4") then
      read(unit=words(3), fmt=*, iostat=ios) on_site(8)
      on_site(9) = on_site(8)
      if(ios /= 0)  &
        call log_error("readElementFile","Something wrong in the eg tight-binding parameters of element " // trim(Name) //".")
    else
      do j=3,4
        read(unit=words(j), fmt=*, iostat=ios) on_site(j+5)
        if(ios /= 0)  &
          call log_error("readElementFile","Something wrong in the eg tight-binding parameters of element " // trim(Name) //".")
      end do
    end if

    ! Setting up on-site terms
    allocate(s%Types(n)%onSite(s%Types(n)%nOrb,s%Types(n)%nOrb))
    s%Types(n)%onSite = 0._dp
    do j=1,s%Types(n)%nOrb
      s%Types(n)%onSite(j,j) = on_site(s%Types(n)%Orbs(j))
    end do

    ! Reading two-center integrals
    allocate(s%Types(n)%Hopping(10,s%nStages))
    do j = 1, nn_stages
      do k = 1, 10
        read(f_unit, fmt='(A)', iostat = ios) line
        read(unit=line, fmt=*, iostat=ios) (words(l), l=1,10)
        if(j>s%nStages) exit
        read(unit=words(4), fmt=*, iostat=ios) s%Types(n)%Hopping(k,j)
        !s%Types(n)%Hopping(j,i) = s%Types(n)%Hopping(j,i) * (a0_corr ** expon(j)) ! Correction of hopping parameter by scaling law.
      end do
    end do

    close(f_unit)

    ! Number of unit cells to be generated along each dimensions
    ! For d dimensions it is (2*n+1)^d
    nCells = (2*s%nStages+1)**(s%Types(n)%isysdim)

    ! Allocate array for nn distances known to the system
    allocate(s%Types(n)%stage(s%nStages, s%Types(n)%nAtoms))
    ! Allocate array for all atoms in all unit cells
    allocate(list(nCells*s%Types(n)%nAtoms))

    allocate( localDistances(nCells * s%Types(n)%nAtoms, s%Types(n)%nAtoms))

    nat = 0

    localDistances = 1.e12_dp
    s%Types(n)%stage = 1.e12_dp
    do i = 1, nCells
      do j = 1, s%Types(n)%nAtoms
        nat = nat + 1

        ! Determine unit-cell indices
        select case(s%Types(n)%isysdim)
        case(3)
          list(nat)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(nat)%Cell(2) = mod( (i-1)/(2*s%nStages+1),(2*s%nStages+1) ) - s%nStages
          list(nat)%Cell(3) = mod( (i-1)/((2*s%nStages+1)*(2*s%nStages+1)),(2*s%nStages+1) ) - s%nStages
        case(2)
          list(nat)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(nat)%Cell(2) = mod( (i-1)/(2*s%nStages+1),(2*s%nStages+1) ) - s%nStages
          list(nat)%Cell(3) = 0
        case default
          list(nat)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(nat)%Cell(2) = 0
          list(nat)%Cell(3) = 0
        end select

        ! Cell position is R = i*a1 + j*a2 + k*a3
        list(nat)%CellVector = list(nat)%Cell(1) * s%Types(n)%a1 + list(nat)%Cell(2) * s%Types(n)%a2 + list(nat)%Cell(3) * s%Types(n)%a3
        ! Atom position is r = R + r_j
        list(nat)%Position = position(:,j) + list(nat)%CellVector

        ! Defining what kind of atom it is
        list(nat)%BasisIndex = j
        list(nat)%Material = imaterial(j)

        ! Allocate arrays for distances and directional cosines to all atoms in the central unit cell
        allocate(list(nat)%Distance(s%Types(n)%nAtoms))
        allocate(list(nat)%dirCos(3,s%Types(n)%nAtoms))

        ! Calculate distances and directional cosines
        do k = 1, s%Types(n)%nAtoms
          list(nat)%Distance(k) = vecDist(list(nat)%Position, position(:,k))
          list(nat)%dirCos(:,k) = 0._dp
          if(list(nat)%Distance(k) <= 1.e-9_dp) cycle
          list(nat)%dirCos(:,k) = (list(nat)%Position - position(:,k)) / list(nat)%Distance(k)

          ! Sort distances *new*
          localDistances(nat,k) = list(nat)%Distance(k)
          l = nat - 1
          do while(1 <= l)
            ! If distance of current atoms is larger than what is saved at position l, exit loop
            if(localDistances(l,k) - list(nat)%Distance(k) < 1.e-9_dp) exit
            localDistances(l+1,k) = localDistances(l,k)
            l = l - 1
          end do
          localDistances(l+1,k) = list(nat)%Distance(k)
        end do
      end do
    end do

    s%Types(n)%stage(1,:) = localDistances(1,:)
    do j = 1, s%Types(n)%nAtoms
      l = 1
      do i = 2, nat
        if(abs(localDistances(i,j) - s%Types(n)%stage(l,j)) < s%Types(n)%stage(1,j) * s%relTol) cycle
        l = l + 1
        if(l > s%nStages) exit
        s%Types(n)%stage(l,j) = localDistances(i,j)
      end do
    end do
    deallocate(localDistances)

  end subroutine readElementFile

  pure subroutine intd(sss,pps,ppp,dds,ddp,ddd,sps,sds,pds,pdp,w,b)
    use mod_kind, only: dp
    use mod_constants, only: sq3
    implicit none
    real(dp), intent(in)  :: sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd,w(3)
    real(dp), intent(out) :: b(9,9)
    real(dp)  :: x,y,z,xx,xy,yy,yz,zz,zx,xxyy,yyzz,zzxx
    real(dp)  :: aux,aux1,aux2,aux3,aux4,r3,f1,f2,f3,f4,f5,f8,g1,g2
    x=w(1)
    y=w(2)
    z=w(3)
    xx=x*x
    xy=x*y
    yy=y*y
    yz=y*z
    zz=z*z
    zx=z*x
    xxyy=xx*yy
    yyzz=yy*zz
    zzxx=zz*xx

    r3=sq3

    f1=xx+yy
    f2=xx-yy
    f3=zz-.5_dp*f1
    f8=3._dp*zz-1._dp
    f4=.5_dp*r3*f2*pds
    f5=.5_dp*f8*pds

    g1=1.5_dp*f2*dds
    g2=r3*f3*dds

    aux=pps-ppp
    aux1=r3*sds
    aux2=r3*xx*pds+(1._dp-2._dp*xx)*pdp
    aux3=(r3*yy*pds+(1._dp-2._dp*yy)*pdp)
    aux4=r3*zz*pds+(1._dp-2._dp*zz)*pdp

    b(1,1)=sss
    b(1,2)=x*sps
    b(1,3)=y*sps
    b(1,4)=z*sps
    b(1,5)=xy*aux1
    b(1,6)=yz*aux1
    b(1,7)=zx*aux1
    b(1,8)=.5_dp*f2*aux1
    b(1,9)=.5_dp*f8*sds

    b(2,1)=-b(1,2)
    b(2,2)=xx*pps+(1._dp-xx)*ppp
    b(2,3)=xy*aux
    b(2,4)=zx*aux
    b(2,5)=aux2*y
    b(2,6)=(r3*pds-2._dp*pdp)*xy*z
    b(2,7)=aux2*z
    b(2,8)=(f4+(1._dp-f2)*pdp)*x
    b(2,9)=(f5-r3*zz*pdp)*x

    b(3,1)=-b(1,3)
    b(3,2)= b(2,3)
    b(3,3)=yy*pps+(1._dp-yy)*ppp
    b(3,4)=yz*aux
    b(3,5)=aux3*x
    b(3,6)=aux3*z
    b(3,7)= b(2,6)
    b(3,8)=(f4-(1._dp+f2)*pdp)*y
    b(3,9)=(f5-r3*zz*pdp)*y

    b(4,1)=-b(1,4)
    b(4,2)= b(2,4)
    b(4,3)= b(3,4)
    b(4,4)=zz*pps+(1._dp-zz)*ppp
    b(4,5)= b(2,6)
    b(4,6)=aux4*y
    b(4,7)=aux4*x
    b(4,8)=(f4-f2*pdp)*z
    b(4,9)=(f5+r3*f1*pdp)*z

    b(5,1)= b(1,5)
    b(5,2)=-b(2,5)
    b(5,3)=-b(3,5)
    b(5,4)=-b(4,5)
    b(5,5)=3._dp*xxyy*dds+(f1-4._dp*xxyy)*ddp+(zz+xxyy)*ddd
    b(5,6)=(3._dp*yy*dds+(1._dp-4._dp*yy)*ddp+(yy-1._dp)*ddd)*zx
    b(5,7)=(3._dp*xx*dds+(1._dp-4._dp*xx)*ddp+(xx-1._dp)*ddd)*yz
    b(5,8)=(g1-2._dp*f2*ddp+.5_dp*f2*ddd)*xy
    b(5,9)=(g2-2._dp*r3*zz*ddp+.5_dp*r3*(1._dp+zz)*ddd)*xy

    b(6,1)= b(1,6)
    b(6,2)=-b(2,6)
    b(6,3)=-b(3,6)
    b(6,4)=-b(4,6)
    b(6,5)= b(5,6)
    b(6,6)=3._dp*yyzz*dds+(yy+zz-4._dp*yyzz)*ddp+(xx+yyzz)*ddd
    b(6,7)=(3._dp*zz*dds+(1._dp-4._dp*zz)*ddp+(zz-1._dp)*ddd)*xy
    b(6,8)=(g1-(1._dp+2._dp*f2)*ddp+(1._dp+.5_dp*f2)*ddd)*yz
    b(6,9)=(g2+r3*(f1-zz)*ddp-.5_dp*r3*f1*ddd)*yz

    b(7,1)= b(1,7)
    b(7,2)=-b(2,7)
    b(7,3)=-b(3,7)
    b(7,4)=-b(4,7)
    b(7,5)= b(5,7)
    b(7,6)= b(6,7)
    b(7,7)=3._dp*zzxx*dds+(zz+xx-4._dp*zzxx)*ddp+(yy+zzxx)*ddd
    b(7,8)=(g1+(1._dp-2._dp*f2)*ddp-(1._dp-.5_dp*f2)*ddd)*zx
    b(7,9)=(g2+r3*(f1-zz)*ddp-.5_dp*r3*f1*ddd)*zx

    b(8,1)= b(1,8)
    b(8,2)=-b(2,8)
    b(8,3)=-b(3,8)
    b(8,4)=-b(4,8)
    b(8,5)= b(5,8)
    b(8,6)= b(6,8)
    b(8,7)= b(7,8)
    b(8,8)=.75_dp*f2*f2*dds+(f1-f2*f2)*ddp+(zz+.25_dp*f2*f2)*ddd
    b(8,9)=.5_dp*f2*g2-r3*zz*f2*ddp+.25_dp*r3*(1._dp+zz)*f2*ddd

    b(9,1)= b(1,9)
    b(9,2)=-b(2,9)
    b(9,3)=-b(3,9)
    b(9,4)=-b(4,9)
    b(9,5)= b(5,9)
    b(9,6)= b(6,9)
    b(9,7)= b(7,9)
    b(9,8)= b(8,9)
    b(9,9)=f3*f3*dds+3._dp*zz*f1*ddp+.75_dp*f1*f1*ddd
  end subroutine intd

end module TightBinding
