!-------------------------------------------------------------------------------
! TITAN - Time-dependent description of Itinerant electrons: Transport and Angular momentum properties of Nanostructures
!-------------------------------------------------------------------------------
!
! MODULE: mod_system
!
!> @author
!> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> Describes the system in real and reciprocal space.
!> Contains function to generate all neighbors and the Brillouin Zone for a layered system or bulk.
!
! REVISION HISTORY:
! 28 April 2017 - Initial Version
!-------------------------------------------------------------------------------
module mod_system
  use mod_kind, only: dp
  use AtomTypes, only: BasisAtom, NeighborAtom, AtomType
  implicit none

  type :: System
    !! Variables for the system to be calculated
    character(len=200) :: Name = ""
    !! Name of the system
    real(dp) :: a0
    !! Lattice parameter
    real(dp), dimension(3) :: a1, a2, a3
    !! Lattice vectors
    real(dp), dimension(3) :: b1, b2, b3
    !! Reciprocal vectors
    real(dp) :: vol
    !! Volume of Brillouin Zone
    integer :: isysdim = 3
    !! Dimension of the system: 3D (default), 2D or 1D

    integer :: nAtoms = 0
    !! Number of atoms in the system
    type(BasisAtom), dimension(:), allocatable :: Basis
    integer :: nNeighbors
    real(dp)  :: Ef
    !! Fermi energy
    real(dp)  :: totalOccupation = 0
    !! Total occupation of the system
    type(NeighborAtom), dimension(:), allocatable :: Neighbors
    integer :: nStages = 0
    !! Number of nearest neighbors
    real(dp) :: relTol
    !! Tolerance for shell radius
    real(dp), dimension(:,:), allocatable :: Distances
    !! List of all distances in nnstages range; size (nStages, nAtoms)
    integer :: nTypes = 0
    !! Number of different atom types
    type(AtomType), dimension(:), allocatable :: Types
    !! List of types
  end type System

  type(System) :: sys

  integer :: n0sc1 !< first neighbor to calculate the in-plane spin and charge current
  integer :: n0sc2 !< last neighbor to calculate the in-plane spin and charge current
  integer :: n0sc  !< Number of neighbors to calculate currents
  real(dp), dimension(3) :: pln_normal

  integer, dimension(:,:), allocatable :: ia
  integer, dimension(:,:), allocatable :: ia_sc

contains

  subroutine initHamiltkStride(nAtoms, nOrb)
    implicit none
    integer, intent(in) :: nAtoms, nOrb
    integer :: i, offsetParameter

    offsetParameter = nAtoms*nOrb*2

    if(allocated(ia)) deallocate(ia)
    if(allocated(ia_sc)) deallocate(ia_sc)
    allocate(ia(4,nAtoms))
    allocate(ia_sc(4,nAtoms))
    do i = 1, nAtoms
      ia(1,i) = (i-1) * 2 * nOrb + 1
      ia(2,i) = ia(1,i) + nOrb - 1
      ia(3,i) = ia(1,i) + nOrb
      ia(4,i) = ia(3,i) + nOrb - 1
      !superconductivity block has doubled dimensions in each spin
      ia_sc(1,i) = (i-1) * 2 * nOrb + 1
      ia_sc(2,i) = ia_sc(1,i) + nOrb*2 - 1
      ia_sc(3,i) = ia_sc(1,i) + offsetParameter
      ia_sc(4,i) = ia_sc(3,i) + nOrb*2 - 1
    end do
  end subroutine initHamiltkStride

  subroutine deallocate_System_variables()
  !! This subroutine deallocates the system variables
    implicit none

    deallocate( sys%Basis )
    deallocate( sys%Neighbors )
    deallocate( sys%Distances )
    deallocate( sys%Types )
    deallocate( ia )
    deallocate( ia_sc )

  end subroutine deallocate_System_variables

end module mod_system
