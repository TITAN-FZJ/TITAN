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

  type :: System_type
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
  end type System_type

  type(System_type) :: sys

  integer :: n0sc1 !< first neighbor to calculate the in-plane spin and charge current
  integer :: n0sc2 !< last neighbor to calculate the in-plane spin and charge current
  integer :: n0sc  !< Number of neighbors to calculate currents
  real(dp), dimension(3) :: pln_normal

  integer, dimension(:,:), allocatable :: ia
  integer, dimension(:,:), allocatable :: ia_sc

contains

  subroutine initHamiltkStride(s, superCond)
    use AtomTypes,      only: lorbital_selection,orbitals
    use mod_parameters, only: dimH,dimspinAtoms,dimens,dimHsc
    implicit none
    type(System_type), intent(in) :: s
    integer,           intent(in) :: superCond
    integer :: i,nOrb,offsetParameter

    if(.not.lorbital_selection) then
      nOrb = size(orbitals)
      dimH = s%nAtoms*nOrb*2
      dimHsc  = dimH*superCond
      dimspinAtoms = 4 * s%nAtoms
      dimens = 4 * s%nAtoms * nOrb * nOrb

      offsetParameter = s%nAtoms*nOrb*2

      if(allocated(ia)) deallocate(ia)
      if(allocated(ia_sc)) deallocate(ia_sc)
      allocate(ia(4,s%nAtoms))
      allocate(ia_sc(4,s%nAtoms))
      do i = 1, s%nAtoms
        ia(1,i) = (i-1) * 2 * nOrb + 1   ! Begin up
        ia(2,i) = ia(1,i) + nOrb - 1     ! End up
        ia(3,i) = ia(2,i) + 1            ! Begin down
        ia(4,i) = ia(3,i) + nOrb - 1     ! End down
        ! Superconductivity block has doubled dimensions in each spin
        ia_sc(1,i) = (i-1) * 2 * nOrb + 1           ! Begin first block (electrons) 1 to 2*nOrb
        ia_sc(2,i) = ia_sc(1,i) + nOrb*2 - 1        ! End first block (electrons)
        ia_sc(3,i) = ia_sc(1,i) + dimH              ! Begin second block (holes) 1 to 2*nOrb + dimH
        ia_sc(4,i) = ia_sc(3,i) + nOrb*2 - 1        ! End second block (holes)
      end do      
    end if
  end subroutine initHamiltkStride


  subroutine initConversionMatrices(nAtoms, nOrbs)
  !! This subroutine mounts the conversion matrices from 4 to 2 ranks
    use mod_parameters, only: sigmai2i,sigmaimunu2i,sigmaijmunu2i,isigmamu2n,n2isigmamu
    implicit none
    integer, intent(in) :: nAtoms, nOrbs
    integer :: nu, mu, i, sigma, j, kount

    !------------------------- Conversion arrays  --------------------------
    do nu = 1, nOrbs
      do mu = 1, nOrbs
        do i = 1, nAtoms
          do sigma = 1, 4
            sigmaimunu2i(sigma,i,mu,nu) = (sigma-1)*nAtoms*nOrbs*nOrbs + (i-1)*nOrbs*nOrbs + (mu-1)*nOrbs + nu
            do j = 1, nAtoms
              sigmaijmunu2i(sigma,i,j,mu,nu) = (sigma-1)*nAtoms*nAtoms*nOrbs*nOrbs + (i-1)*nAtoms*nOrbs*nOrbs + (j-1)*nOrbs*nOrbs + (mu-1)*nOrbs + nu
            end do
          end do
        end do
      end do
    end do

    do i = 1, nAtoms
      do sigma = 1, 4
        sigmai2i(sigma,i) = (sigma-1)*nAtoms + i
      end do
    end do

    ! Conversion array from local atomic orbitals to (i,sigma,mu) to eigenvectors n
    do i = 1, nAtoms
      do sigma = 1, 2
        do mu = 1, nOrbs
          kount = (i-1)*2*nOrbs + (sigma-1)*nOrbs + mu
          isigmamu2n(i,sigma,mu) = kount
          n2isigmamu(kount,1) = i
          n2isigmamu(kount,2) = sigma
          n2isigmamu(kount,3) = mu
        end do
      end do
    end do

  end subroutine initConversionMatrices


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
