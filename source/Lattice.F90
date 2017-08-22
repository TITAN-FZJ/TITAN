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
  use mod_f90_kind, only: double
  implicit none

  integer :: LatticeMode = 2

  contains

    subroutine initLattice(s)
      use mod_monoBasis, only: monoBasis => create_basis
      use mod_polyBasis, only: polyBasis => read_basis
      use mod_system, only: System
      use AtomTypes, only: NeighborAtom
      implicit none
      type(System), intent(inout) :: s
      type(NeighborAtom), dimension(:), allocatable :: list
      integer :: size

      
      if(LatticeMode == 1) then
        call monoBasis()
      elseif(LatticeMode == 2) then
        call polyBasis("basis", s)
      end if

      call generateNeighbors(s, list, size)
      call sortNeighbors(s, list,size)

      return
    end subroutine initLattice

    subroutine generateNeighbors(s,list,size)
      use AtomTypes, only: NeighborAtom
      use mod_system, only: System
      implicit none
      type(System), intent(inout) :: s
      type(NeighborAtom), intent(out), dimension(:), allocatable :: list
      integer,intent(out) :: size
      integer :: nCells
      integer :: i,j,k,l

      nCells = (2*s%nStages+1)*(2*s%nStages+1)
      if(s%lbulk) nCells = nCells * (2*s%nStages+1)
      allocate(s%Distances(s%nStages, s%nAtoms))
      allocate(list(nCells*s%nAtoms))

      size = 0
      s%Distances = 1.d+12
      do i = 1, nCells
        do j = 1, s%nAtoms
          size = size + 1

          list(size)%Cell = [mod((i-1)/(2*s%nStages+1),(2*s%nStages+1))-s%nStages, mod((i-1),(2*s%nStages+1))-s%nStages, 0]
          if(s%lbulk) list(size)%Cell(3) = mod((i-1)/((2*s%nStages+1)*(2*s%nStages+1)), (2*s%nStages+1)) - s%nStages

          list(size)%CellVector = list(size)%Cell(1) * s%a1 + list(size)%Cell(2) * s%a2 + list(size)%Cell(3) * s%a3
          list(size)%Position = s%Basis(j)%Position + list(size)%CellVector

          list(size)%BasisIndex = j
          list(size)%Material = s%Basis(j)%Material

          allocate(list(size)%Distance(s%nAtoms))
          allocate(list(size)%dirCos(3,s%nAtoms))
          do k = 1, s%nAtoms
            list(size)%Distance(k) = sqrt(dot_product(list(size)%Position-s%Basis(k)%Position, list(size)%Position-s%Basis(k)%Position))
            list(size)%dirCos(:,k) = 0.d0
            if(list(size)%Distance(k) <= 1.0d-9) cycle
            list(size)%dirCos(:,k) = (list(size)%Position - s%Basis(k)%Position) / sqrt(dot_product(list(size)%Position - s%Basis(k)%Position,list(size)%Position - s%Basis(k)%Position))
            l = s%nStages
            do while(1 <= l)
              if(s%Distances(l,k) - list(size)%Distance(k) < 1.d-9) exit
              if(l < s%nStages) s%Distances(l+1,k) = s%Distances(l,k)
              l = l - 1
            end do
            if(l == 0) then
              s%Distances(l+1,k) = list(size)%Distance(k)
            elseif(abs(s%Distances(l,k)-list(size)%Distance(k)) >= 1.d-9 .and. l < s%nStages) then
              s%Distances(l+1,k) = list(size)%Distance(k)
            end if
          end do
        end do
      end do
    end subroutine generateNeighbors

    subroutine sortNeighbors(s,list, size)
      use AtomTypes, only: NeighborAtom, add_elem
      use mod_system, only: System
      implicit none
      type(System), intent(inout) :: s
      integer, intent(in) :: size
      integer :: nNeighbors
      type(NeighborAtom), dimension(size), intent(inout) :: list

      integer :: i,j,k
      integer :: matchedNeighbor(s%nAtoms)
      logical :: found


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
            if(abs(list(i)%Distance(j) - s%Distances(k,j)) < 1.d-9) then
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
          call add_elem(s%Basis(j)%NeighborList(matchedNeighbor(j), list(nNeighbors)%BasisIndex), nNeighbors)
        end do
      end do

      s%nNeighbors = nNeighbors
      allocate(s%Neighbors(nNeighbors))
      s%Neighbors(1:nNeighbors) = list(1:nNeighbors)

      return
    end subroutine sortNeighbors


    subroutine writeLattice(s)
      use AtomTypes, only: NeighborIndex
      use mod_system, only: System
      implicit none
      type(System), intent(in) :: s
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
      return

    end subroutine writeLattice
end module Lattice
