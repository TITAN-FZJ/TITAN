!------------------------------------------------------------------------------------!
! TITAN - Time-dependent Transport and Angular momentum properties of Nanostructures !
!------------------------------------------------------------------------------------!
!
! MODULE: mod_Lattice
!
!> @author
!> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> Setup Real Space Lattice from given Lattice Vectors and Basis Atoms
!
! REVISION HISTORY:
! 05 July 2017 - Initial Version
!------------------------------------------------------------------------------------!

module Lattice
  use mod_kind,  only: dp
  use AtomTypes, only: NeighborAtom
  implicit none

  type(NeighborAtom), dimension(:), allocatable, private :: list
  !! Array containing all generated atoms

contains

  subroutine initLattice(s)
    use mod_system, only: System_type
    implicit none
    type(System_type), intent(inout) :: s
    integer :: isize

    call generateNeighbors(s, isize)
    call sortNeighbors(s, isize)

    deallocate(list)

  end subroutine initLattice

  subroutine generateNeighbors(s,isize)
  !! Generate list of all possible Neighbors to the atoms in unit cell 0
  !! (including positions of cells, position of atoms, distances to atoms at unit cell 0,
  !! atom types and directional cosines) and stores the distances of each 
  !! stage (within given tolerance)
    use mod_system, only: System_type
    use mod_tools,  only: vecDist
    implicit none
    type(System_type), intent(inout) :: s
    !! System for which lattice is generated
    integer,           intent(out)   :: isize
    !! Counter for number of atoms generated, return total number generated
    integer :: nCells
    !! Number of generated unit cells
    integer :: i,j,k,l
    !! Loop variables
    real(dp), allocatable, dimension(:,:) :: localDistances
    !! Array containing *all* distances

    ! Number of unit cells to be generated along each dimensions
    ! For d dimensions it is (2*n+1)^d
    nCells = (2*s%nStages+1)**(s%isysdim)

    ! Allocate array for nn distances known to the system
    allocate(s%Distances(s%nStages, s%nAtoms))
    ! Allocate array for all atoms in all unit cells
    allocate(list(nCells*s%nAtoms))

    allocate( localDistances(nCells * s%nAtoms, s%nAtoms))

    ! "isize" is added up until all atoms in all unit cells, nCells*s%nAtoms
    isize = 0
    ! Set known distances to value which is guaranteed to be bigger than any possible distance in the system
    localDistances = 1.e12_dp
    s%Distances = 1.e12_dp
    ! Looping over all the different possible cells
    do i = 1, nCells
      ! ...and atoms in the unit cell
      do j = 1, s%nAtoms
        ! "isize" is the current atom
        isize = isize + 1

        ! Determine the indices of the unit-cell containing atom "isize"
        select case(s%isysdim)
        case(3)
          list(isize)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(isize)%Cell(2) = mod( (i-1)/(2*s%nStages+1),(2*s%nStages+1) ) - s%nStages
          list(isize)%Cell(3) = mod( (i-1)/((2*s%nStages+1)*(2*s%nStages+1)),(2*s%nStages+1) ) - s%nStages
        case(2)
          list(isize)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(isize)%Cell(2) = mod( (i-1)/(2*s%nStages+1),(2*s%nStages+1) ) - s%nStages
          list(isize)%Cell(3) = 0
        case default
          list(isize)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(isize)%Cell(2) = 0
          list(isize)%Cell(3) = 0
        end select

        ! Cell position is R_i = i*a1 + j*a2 + k*a3
        list(isize)%CellVector = list(isize)%Cell(1) * s%a1 + list(isize)%Cell(2) * s%a2 + list(isize)%Cell(3) * s%a3
        ! Atom position is r = R_i + r_j
        list(isize)%Position = s%Basis(j)%Position + list(isize)%CellVector

        ! Defining what kind of atom it is (defined by which atom in the unit cell j)
        list(isize)%BasisIndex = j
        list(isize)%Material = s%Basis(j)%Material

        ! Allocate arrays for distances and directional cosines to all atoms in the central unit cell
        allocate(list(isize)%Distance(s%nAtoms))
        allocate(list(isize)%dirCos(3,s%nAtoms))

        ! Calculate distances to other atoms k in the unit cell 0, and respective directional cosines
        do k = 1, s%nAtoms
          list(isize)%Distance(k) = vecDist(list(isize)%Position, s%Basis(k)%Position)
          list(isize)%dirCos(:,k) = 0._dp
          ! If atom is the same, cycle
          if(list(isize)%Distance(k) <= 1.0e-9_dp) cycle
          ! Unit vector pointing from atom "k" to "isize"
          list(isize)%dirCos(:,k) = (list(isize)%Position - s%Basis(k)%Position) / list(isize)%Distance(k)

          ! Sort distances between all previous atoms "isize" and unit cell 0 atoms "k"
          ! The larger the first index, the larger the distance
          localDistances(isize,k) = list(isize)%Distance(k)
          l = isize - 1
          do while(l >= 1)
            ! If distance of current atoms is larger than what is saved 
            ! at a previous position (l), exit loop
            if( list(isize)%Distance(k) > localDistances(l,k)) exit
            ! Otherwise, move previous position l one index up (to open space for current one)
            localDistances(l+1,k) = localDistances(l,k)
            l = l - 1
          end do
          ! When the position is larger, fits the current one in place (l+1, due to l = l-1 at the end of the loop)
          localDistances(l+1,k) = list(isize)%Distance(k)
        end do
      end do
    end do

    ! Loop to store distances up to required stage (within the tolerance s%relTol)

    ! Getting the smallest distance for all atoms in unit cell 0 (first stage)
    s%Distances(1,:) = localDistances(1,:)
    ! Loop over all atoms in unit cell 0
    do j = 1, s%nAtoms
      l = 1
      ! Loop over all atoms (in all unit cells)
      do i = 2, isize
        ! If current distance is inside tolerance (relTol percent of first stage distance), cycle
        if(abs(localDistances(i,j) - s%Distances(l,j)) < s%Distances(1,j) * s%relTol) cycle

        ! Otherwise store the next stage, if desired
        l = l + 1
        if(l > s%nStages) exit
        s%Distances(l,j) = localDistances(i,j)
      end do
    end do
    deallocate(localDistances)
  end subroutine generateNeighbors

  subroutine sortNeighbors(s, isize)
  !! Sorts the neighbors that are inside the stages required
    use AtomTypes,  only: NeighborIndex
    use mod_system, only: System_type
    implicit none
    integer,           intent(in)    :: isize
    type(System_type), intent(inout) :: s
    integer :: nNeighbors

    integer :: i,j,k
    integer :: matchedNeighbor(s%nAtoms)
    logical :: found
    type(NeighborIndex), pointer :: local => null()

    ! Initialize Neighbor Lists for all Basis Atoms
    do i = 1, s%nAtoms
      allocate(s%Basis(i)%NeighborList(s%nStages, s%nAtoms))
      do j = 1, s%nAtoms
        do k = 1, s%nStages
          nullify(s%Basis(i)%NeighborList(k,j)%head)
        end do
      end do
    end do

    nNeighbors = 0
    ! Loop over all atoms of all cells
    do i = 1, isize
      matchedNeighbor = 0
      found = .false.
      ! Loop over all atoms of cell 0
      do j = 1, s%nAtoms
        ! Loop over all stages of its neighbors
        do k = 1, s%nStages
          ! If atom i is inside shells of atom j, found a neighbor
          if(abs(list(i)%Distance(j) - s%Distances(k,j)) < s%Distances(1,j) * s%relTol) then
            found = .true.
            matchedNeighbor(j) = k
            exit ! If neighbor is found, can leave the stages loop (as they are already inside)
          end if
        end do
      end do

      if(.not. found) cycle ! If current site "i" is not a neighbor of any atom "j" in the unit cell 0, go to next

      ! If found, site "i" is a neighbor of at least one site
      nNeighbors = nNeighbors + 1
      ! Move site i to ordered positions (increasing "nNeighbors")
      list(nNeighbors) = list(i)
      list(nNeighbors)%Distance = list(i)%Distance

      ! Ordering list sites in unit cell that are neighbors of site "i"
      do j = 1, s%nAtoms
        if(matchedNeighbor(j) == 0) cycle

        ! Storing the location of the current first element of the (basis list) 
        local => s%Basis(j)%NeighborList(matchedNeighbor(j),list(nNeighbors)%BasisIndex)%head

        ! Resetting the pointer of the (basis list) and pointing it to a newly created element
        nullify( s%Basis(j)%NeighborList(matchedNeighbor(j),list(nNeighbors)%BasisIndex)%head )
        allocate( s%Basis(j)%NeighborList(matchedNeighbor(j),list(nNeighbors)%BasisIndex)%head )
        ! At this point (basis list)%head points to a list with one element 

        ! Now make this new first element to be "nNeighbors" that is being counted outside, and the next element will be the previous head.
        s%Basis(j)%NeighborList(matchedNeighbor(j),list(nNeighbors)%BasisIndex)%head%index = nNeighbors
        s%Basis(j)%NeighborList(matchedNeighbor(j),list(nNeighbors)%BasisIndex)%head%next  => local

      end do
    end do

    s%nNeighbors = nNeighbors
    allocate(s%Neighbors(nNeighbors))
    s%Neighbors(1:nNeighbors) = list(1:nNeighbors)

  end subroutine sortNeighbors


  subroutine writeLattice(s)
  !! Writes out lattice information to file "Atoms, 
  !! including Bravais vectors, position of basis atoms 
  !! and their respective neighbors
    use AtomTypes,  only: NeighborIndex
    use mod_system, only: System_type
    implicit none
    type(System_type), intent(in) :: s
    integer :: i,j,k,ios
    integer :: out_unit = 99
    type(NeighborIndex), pointer :: current

    open(out_unit, file="Atoms", iostat = ios)
    if(ios /= 0) stop "[output_lattice] Something went wrong"

    write(out_unit, *) "# Lattice constant: "
    write(out_unit, *) s%a0
    write(out_unit, *) "# Lattice vectors: "
    write(out_unit, *) s%a1(1), s%a1(2), s%a1(3)
    write(out_unit, *) s%a2(1), s%a2(2), s%a2(3)
    write(out_unit, *) s%a3(1), s%a3(2), s%a3(3)
    write(out_unit, *) "# Basis Atoms: "
    do i = 1, s%nAtoms
      write(out_unit, *) trim(s%Types(s%Basis(i)%Material)%Name), s%Basis(i)%Position(1), s%Basis(i)%Position(2), s%Basis(i)%Position(3)
    end do

    write(out_unit, *) "========================================================================="
    write(out_unit, *) "#---------------------------- Neighbor Atoms ----------------------------"
    do i = 1, s%nAtoms
      write(out_unit, *) "# Atom: ", i
      do j = 1, s%nStages
        write(out_unit, *) "# Stage: ", j
        do k = 1, s%nAtoms
          current => s%Basis(i)%NeighborList(j,k)%head
          do while(associated(current))
            write(out_unit,*) trim(s%Types(s%Neighbors(current%index)%Material)%Name), &
                              s%Neighbors(current%index)%Position(1), s%Neighbors(current%index)%Position(2), s%Neighbors(current%index)%Position(3), &
                              s%Neighbors(current%index)%Distance(i)
            current => current%next
          end do
        end do
      end do
    end do
    close(out_unit)

  end subroutine writeLattice
end module Lattice
