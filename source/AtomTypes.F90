module AtomTypes
  use mod_kind, only: dp

  type NeighborIndex
    integer :: index
    type(NeighborIndex), pointer :: next => null()
  end type NeighborIndex

  type NeighborHead
    type(NeighborIndex), pointer :: head => null()
  end type

  type :: Atom
    real(dp), dimension(3) :: Position
    integer :: Material
  end type

  type, extends(Atom) :: BasisAtom
    type(NeighborHead), dimension(:,:), allocatable :: NeighborList
  end type BasisAtom

  type, extends(Atom) :: NeighborAtom
    integer :: BasisIndex
    integer, dimension(3) :: Cell
    real(dp), dimension(3) :: CellVector
    real(dp), dimension(:), allocatable :: Distance
    !! Distance to basis atoms
    real(dp), dimension(:,:), allocatable :: dirCos
    !! Directional cosine to basis atoms
    real(dp), dimension(:,:,:), allocatable :: t0i
    !! Hopping to basis atoms; size (nOrb, nOrb, nAtoms)
    logical, dimension(:), allocatable :: isHopping
    !! Flags saying that the atom is in neighbor range of basis atoms; size (nAtoms)
  end type NeighborAtom

  type :: AtomType
    !! Store information for each element type
    character(len=50) :: Name
    !! Element name
    integer :: nAtoms
    !! Number of atoms in elemental file
    integer :: nAtomEl=0
    !! Number of atoms of a given element in elemental file
    character(len=50), dimension(:), allocatable :: Types
    !! Element Types inside elemental file
    logical, dimension(:), allocatable :: lelement
    !! Logical variable storing where are the elements of a given element in elemental file
    real(dp), dimension(:,:), allocatable :: onSite
    !! on site matrix; size (nOrb, nOrb)
    real(dp), dimension(:,:), allocatable :: Hopping
    !! hopping matrix; size (10, nStages)
    real(dp) :: LambdaP=0._dp, LambdaD=0._dp
    !! SOC Coupling constants
    real(dp) :: FermiLevel
    !! Fermi level
    real(dp) :: Occupation, OccupationS, OccupationP, OccupationD
    !! Total, s, p and d occupations
    real(dp) :: LatticeConstant
    !! Lattice constant
    integer      :: isysdim
    !! Dimension of the system (1D, 2D, 3D)
    real(dp), dimension(3) :: a1,a2,a3
    !! Lattice vectors
    real(dp), dimension(:,:), allocatable :: Stage
    !! Neighbor distances
    real(dp) :: Un, Um
    !! Effective Coulomb interaction strength - Hubbard U (charge and magnetic)
    real(dp), dimension(9) :: lambda
    !! Superconducting coupling strength
    real(dp), dimension(:,:), allocatable :: rho0
    !! Initial occupation per orbital obtained from hopping parameters only
    real(dp), dimension(:), allocatable   :: rhod0
    !! Initial occupation of d orbitals obtained from hopping parameters only
  end type AtomType

end module AtomTypes
