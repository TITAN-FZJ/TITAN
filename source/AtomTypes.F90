module AtomTypes
use mod_f90_kind, only: double

type NeighborIndex
  integer :: index
  type(NeighborIndex), pointer :: next => null()
end type NeighborIndex

type NeighborHead
  type(NeighborIndex), pointer :: head => null()
end type

type :: Atom
  real(double), dimension(3) :: Position
  integer :: Material
end type

type, extends(Atom) :: BasisAtom
  type(NeighborHead), dimension(:,:), allocatable :: NeighborList
end type BasisAtom

type, extends(Atom) :: NeighborAtom
  integer :: BasisIndex
  integer, dimension(3) :: Cell
  real(double), dimension(3) :: CellVector
  real(double), dimension(:), allocatable :: Distance
  real(double), dimension(:,:), allocatable :: dirCos
  real(double), dimension(:,:,:), allocatable :: t0i ! nOrb, nOrb, nAtoms
  logical, dimension(:), allocatable :: isHopping
end type NeighborAtom

type :: AtomKind
  character(len=50) :: Name
  real(double), dimension(:,:), allocatable :: onSite ! on site matrix; size (nOrb, nOrb)
  real(double), dimension(:,:), allocatable :: Hopping ! hopping matrix; size (10, nStages)
  real(double) :: Lambda ! SOC Coupling constant
  real(double) :: FermiLevel
  real(double) :: Occupation
  real(double) :: LatticeConstant
end type AtomKind


contains
  subroutine add_elem(list, index)
    implicit none
    type(NeighborHead), intent(inout) :: list
    integer :: index
    type(NeighborIndex), pointer :: local => null()

    local => list%head
    nullify(list%head)
    allocate(list%head)

    list%head%index = index
    list%head%next  => local

    return
  end subroutine add_elem

end module AtomTypes
