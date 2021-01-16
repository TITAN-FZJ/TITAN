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
    integer :: size

    call generateNeighbors(s, size)
    call sortNeighbors(s, size)

    deallocate(list)

  end subroutine initLattice

  subroutine generateNeighbors(s,size)
    use mod_system, only: System_type
    use mod_tools,  only: vecDist
    implicit none
    type(System_type),                             intent(inout) :: s
    !! System for which lattice is generated
    integer,                                       intent(out)   :: size
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

    size = 0
    ! Set known distances to value which is guaranteed to be bigger than any possible distance in the system
    localDistances = 1.e12_dp
    s%Distances = 1.e12_dp
    do i = 1, nCells
      do j = 1, s%nAtoms
        size = size + 1

        ! Determine unit-cell indices
        select case(s%isysdim)
        case(3)
          list(size)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(size)%Cell(2) = mod( (i-1)/(2*s%nStages+1),(2*s%nStages+1) ) - s%nStages
          list(size)%Cell(3) = mod( (i-1)/((2*s%nStages+1)*(2*s%nStages+1)),(2*s%nStages+1) ) - s%nStages
        case(2)
          list(size)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(size)%Cell(2) = mod( (i-1)/(2*s%nStages+1),(2*s%nStages+1) ) - s%nStages
          list(size)%Cell(3) = 0
        case default
          list(size)%Cell(1) = mod( (i-1),(2*s%nStages+1) ) - s%nStages
          list(size)%Cell(2) = 0
          list(size)%Cell(3) = 0
        end select

        ! Cell position is R = i*a1 + j*a2 + k*a3
        list(size)%CellVector = list(size)%Cell(1) * s%a1 + list(size)%Cell(2) * s%a2 + list(size)%Cell(3) * s%a3
        ! Atom position is r = R + r_j
        list(size)%Position = s%Basis(j)%Position + list(size)%CellVector

        ! Defining what kind of atom it is
        list(size)%BasisIndex = j
        list(size)%Material = s%Basis(j)%Material

        ! Allocate arrays for distances and directional cosines to all atoms in the central unit cell
        allocate(list(size)%Distance(s%nAtoms))
        allocate(list(size)%dirCos(3,s%nAtoms))

        ! Calculate distances and directional cosines
        do k = 1, s%nAtoms
          list(size)%Distance(k) = vecDist(list(size)%Position, s%Basis(k)%Position)
          list(size)%dirCos(:,k) = 0._dp
          if(list(size)%Distance(k) <= 1.0e-9_dp) cycle
          list(size)%dirCos(:,k) = (list(size)%Position - s%Basis(k)%Position) / list(size)%Distance(k)

          ! Sort distances *new*
          localDistances(size,k) = list(size)%Distance(k)
          l = size - 1
          do while(1 <= l)
            ! If distance of current atoms is larger than what is saved at position l, exit loop
            if(localDistances(l,k) - list(size)%Distance(k) < 1.e-9_dp) exit
            localDistances(l+1,k) = localDistances(l,k)
            l = l - 1
          end do
          localDistances(l+1,k) = list(size)%Distance(k)

          ! ! Sort distances
          ! l = s%nStages
          ! do while(1 <= l)
          !   ! If distance of current atoms is larger than what is saved for the current nn stage, exit loop
          !   if(s%Distances(l,k) - list(size)%Distance(k) < 1.e-9_dp) exit
          !
          !   ! When a larger stage exists, move current stage l one up and set new value
          !   if(l < s%nStages) s%Distances(l+1,k) = s%Distances(l,k)
          !   l = l - 1
          ! end do
          ! if(l == 0) then
          !   s%Distances(l+1,k) = list(size)%Distance(k)
          ! elseif(abs(s%Distances(l,k)-list(size)%Distance(k)) >= 1.e-9_dp .and. l < s%nStages) then
          !   s%Distances(l+1,k) = list(size)%Distance(k)
          ! end if
        end do
      end do
    end do

    s%Distances(1,:) = localDistances(1,:)
    do j = 1, s%nAtoms
      l = 1
      do i = 2, size
        if(abs(localDistances(i,j) - s%Distances(l,j)) < s%Distances(1,j) * s%relTol) cycle
        l = l + 1
        if(l > s%nStages) exit
        s%Distances(l,j) = localDistances(i,j)
      end do
    end do
    deallocate(localDistances)
  end subroutine generateNeighbors

  subroutine sortNeighbors(s, size)
    use AtomTypes,  only: NeighborIndex
    use mod_system, only: System_type
    implicit none
    integer,      intent(in)    :: size
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

    do i = 1, size
      matchedNeighbor = 0
      found = .false.
      do j = 1, s%nAtoms
        do k = 1, s%nStages
          if(abs(list(i)%Distance(j) - s%Distances(k,j)) < s%Distances(1,j) * s%relTol) then
            found = .true.
            matchedNeighbor(j) = k
            exit
          end if
        end do
      end do

      if(.not. found) cycle

      nNeighbors = nNeighbors + 1
      list(nNeighbors) = list(i)
      list(nNeighbors)%Distance = list(i)%Distance

      do j = 1, s%nAtoms
        if(matchedNeighbor(j) == 0) cycle

        ! Storing the location of the current first element of the (basis list) 
        local => s%Basis(j)%NeighborList(matchedNeighbor(j),list(nNeighbors)%BasisIndex)%head

        ! Creates resetting the pointer of the (basis list) and pointing it to a newly created element
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
