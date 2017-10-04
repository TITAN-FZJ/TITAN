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
  use mod_f90_kind, only: double
  use AtomTypes, only: BasisAtom, NeighborAtom, AtomType
  implicit none

  type :: System
    character(len=200) :: Name = ""
    real(double), dimension(3) :: a1, a2, a3
    real(double) :: a0
    logical :: lbulk = .false.

    integer :: nAtoms = 0
    type(BasisAtom), dimension(:), allocatable :: Basis
    integer :: nNeighbors
    type(NeighborAtom), dimension(:), allocatable :: Neighbors
    integer :: nStages = 0
    real(double), dimension(:,:), allocatable :: Distances
    !! List of all distances in nnstages range; size (nStages, nAtoms)
    integer :: nTypes = 0
    !! Number of different atom types
    type(AtomType), dimension(:), allocatable :: Types
  end type System

  type(System) :: sys

  integer :: n0sc1 !< first neighbor to calculate the in-plane spin and charge current
  integer :: n0sc2 !< last neighbor to calculate the in-plane spin and charge current
  integer :: n0sc  !< Number of neighbors to calculate currents
  real(double), dimension(3) :: pln_normal

  integer, dimension(:,:), allocatable :: ia



 contains

   subroutine initHamiltkStride(nAtoms, nOrb)
     implicit none
     integer, intent(in) :: nAtoms, nOrb
     integer :: i
     allocate(ia(4,nAtoms))
     do i = 1, nAtoms
        ia(1,i) = (i-1) * 2 * nOrb + 1
        ia(2,i) = ia(1,i) + nOrb - 1
        ia(3,i) = ia(1,i) + nOrb
        ia(4,i) = ia(3,i) + nOrb - 1
     end do
     return
   end subroutine initHamiltkStride
end module mod_system