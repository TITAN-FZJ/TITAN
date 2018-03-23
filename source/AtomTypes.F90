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
  !! Distance to basis atoms
  real(double), dimension(:,:), allocatable :: dirCos
  !! Directional cosine to basis atoms
  real(double), dimension(:,:,:), allocatable :: t0i
  !! Hopping to basis atoms; size (nOrb, nOrb, nAtoms)
  logical, dimension(:), allocatable :: isHopping
  !! Flags saying that the atom is in neighbor range of basis atoms; size (nAtoms)
end type NeighborAtom

type :: AtomType
  !! Store information for each element type
  character(len=50) :: Name
  !! Element name
  real(double), dimension(:,:), allocatable :: onSite
  !! on site matrix; size (nOrb, nOrb)
  real(double), dimension(:,:), allocatable :: Hopping
  !! hopping matrix; size (10, nStages)
  real(double) :: LambdaP=0.d0, LambdaD=0.d0
  !! SOC Coupling constants
  real(double) :: FermiLevel
  !! Fermi level
  real(double) :: Occupation, OccupationS, OccupationP, OccupationD
  !! Total, s, p and d occupations
  real(double) :: LatticeConstant
  !! Lattice constant
  real(double), dimension(3) :: a1,a2,a3
  !! Lattice vectors
  real(double), dimension(:), allocatable :: Stage
  !! Neighbor distances
  real(double) :: U
  !! Coulomb strength
  real(double), dimension(:,:), allocatable :: rho0
  !! Initial occupation per orbital obtained from hopping parameters only
  real(double), dimension(:), allocatable   :: rhod0
  !! Initial occupation of d orbitals obtained from hopping parameters only
end type AtomType


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

      end subroutine add_elem

end module AtomTypes
