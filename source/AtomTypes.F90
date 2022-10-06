module AtomTypes
  use mod_kind, only: dp,int32
  implicit none
  
  integer(int32) :: default_nOrb=9
  !! Default number of orbitals
  character(len=3),dimension(9) :: default_Orbs = ["s  ","px ","py ","pz ","dxy","dyz","dzx","dx2","dz2"]
  !! Orbitals that are implemented and recognized on TITAN
  integer(int32) :: default_nsOrb=1, default_sOrbs(1) = [1]
  !! Default s Orbitals 
  integer(int32) :: default_npOrb=3, default_pOrbs(3) = [2,3,4]
  !! Default p Orbitals 
  integer(int32) :: default_ndOrb=5,default_dOrbs(5) = [5,6,7,8,9]
  !! Default d Orbitals 

  type NeighborIndex
    integer :: index
    type(NeighborIndex), pointer :: next => null()
  end type NeighborIndex

  type NeighborHead
    type(NeighborIndex), pointer :: head => null()
  end type

  type :: Atom
    integer  :: Material
    !! Material index of neighboring atom
    real(dp), dimension(3) :: Position
    !! Position of the Atom
  end type

  type, extends(Atom) :: BasisAtom
    type(NeighborHead), dimension(:,:), allocatable :: NeighborList
    !! List of Neighboring atoms
    real(dp) :: Un=0._dp, Um=0._dp
    !! Effective Coulomb interaction strength - Hubbard U (charge and magnetic)
    complex(dp), dimension(:,:), allocatable :: lb, sb
#ifdef _GPU
    complex(dp), dimension(:,:), allocatable, device :: lb_d, sb_d
#endif
    !! Zeeman OAM and spin interaction matrices (on the CPU and on the GPU)
    complex(dp), dimension(:,:), allocatable :: ls
#ifdef _GPU
    complex(dp), dimension(:,:), allocatable, device :: ls_d
#endif
    !! Spin-orbit coupling interaction matrix (on the CPU and on the GPU)
    complex(dp), dimension(:,:,:), allocatable :: lpvec
    !! Angular momentum vector matrices in local frame
    complex(dp), dimension(:,:), allocatable :: onSite
    !! on site matrix for each site in the unit cell; size (nOrb, nOrb)
  end type BasisAtom

  type, extends(Atom) :: NeighborAtom
    integer :: BasisIndex
    !! Index to store atom of the Basis
    integer, dimension(3) :: Cell
    !! Index of the cell (3D)
    real(dp), dimension(3) :: CellVector
    !! Vector position of the cell
    real(dp), dimension(:), allocatable :: Distance
    !! Distance to basis atoms
    real(dp), dimension(:,:), allocatable :: dirCos
    !! Directional cosine to basis atoms
    complex(dp), dimension(:,:,:), allocatable :: t0i
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
    complex(dp), dimension(:,:), allocatable :: onSite
    !! on site matrix for each Atom type; size (nOrb, nOrb)
    real(dp), dimension(:,:), allocatable :: Hopping
    !! hopping matrix; size (10, nStages)
    real(dp) :: LambdaP=0._dp, LambdaD=0._dp
    !! SOC Coupling constants
    real(dp) :: FermiLevel
    !! Fermi level
    
    integer(int32), dimension(:), allocatable :: Orbs
    !! Types of selected orbitals 
    integer(int32), dimension(:), allocatable :: sOrbs,pOrbs,dOrbs
    !! Indices of selected s,p,d orbitals
    integer(int32) :: nOrb,nsOrb,npOrb,ndOrb,nOrb2,nOrb2sc
    !! Number of orbitals, number of s,p,d orbitals
    !! 2*(number of orbitals) (for spin) and supercond*2*nOrb
    real(dp) :: Occupation, OccupationS, OccupationP, OccupationD
    !! Total, s, p and d occupations
    real(dp) :: LatticeConstant
    !! Lattice constant
    integer  :: isysdim
    !! Dimension of the system (1D, 2D, 3D)
    real(dp), dimension(3) :: a1,a2,a3
    !! Lattice vectors
    real(dp), dimension(:,:), allocatable :: Stage
    !! Neighbor distances
    real(dp) :: Un=0._dp, Um=0._dp
    !! Effective Coulomb interaction strength - Hubbard U (charge and magnetic)
    real(dp), dimension(:), allocatable :: lambda
    !! Superconducting coupling strength
    real(dp), dimension(:), allocatable :: rho0
    !! Initial occupation per orbital obtained from hopping parameters only
    real(dp)    :: rhod0,mzd0
    !! Initial occupation and z-component of magnetization of d orbitals obtained from hopping parameters only
    complex(dp) :: mpd0
    !! Initial +-component of magnetization of d orbitals obtained from hopping parameters only
    complex(dp), dimension(:,:,:), allocatable :: lvec
    !! Angular momentum vector matrices in global frame

    complex(dp), dimension(:,:), allocatable :: ident_norb
    !! Identity in orbital space
    complex(dp), dimension(:,:), allocatable :: ident_norb2
    !! Identity in spin and orbital space
    complex(dp), dimension(:,:), allocatable :: ident_dorb
    !! Identity in spin and orbital space (only non-zero on d-orbitals)
    complex(dp), dimension(:,:,:), allocatable :: pauli_orb
    !! Pauli matrices in spin and orbital space (x,y,z,+,-)
    complex(dp), dimension(:,:,:), allocatable :: pauli_dorb
    !! Pauli matrices in spin and orbital space (x,y,z,+,-) (only non-zero on d-orbitals)
  end type AtomType

end module AtomTypes
