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
  use mod_kind,  only: dp,int32
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
    integer(int32) :: isysdim = 3
    !! Dimension of the system: 3D (default), 2D or 1D

    integer(int32) :: nAtoms = 0
    !! Number of atoms in the system
    integer :: ndAtoms = 0
    !! Number of atoms with d-orbitals in the system
    integer(int32), dimension(:), allocatable :: Orbs
    !! Types of selected orbitals
    integer(int32), dimension(:), allocatable :: sOrbs,pOrbs,dOrbs
    !! Indices of selected s,p,d orbitals
    integer(int32) :: nOrb,nOrb2,nOrb2sc
    !! Number of orbitals per atom (when all types use the same, otherwise, max number) 
    !! 2*(number of orbitals) (for spin) and supercond*2*nOrb
    integer(int32) :: total_nOrb
    !! Total number of orbitals
    integer(int32) :: nsOrb,npOrb,ndOrb
    !! Number of s,p,d orbitals
    type(BasisAtom), dimension(:), allocatable :: Basis
    !! Information of the basis
    integer(int32) :: nNeighbors
    !! Number of neighbors to be considered
    real(dp)  :: Ef
    !! Fermi energy
    real(dp)  :: totalOccupation = 0._dp
    !! Total occupation of the system
    type(NeighborAtom), dimension(:), allocatable :: Neighbors
    !! Information of the neighbors
    integer(int32) :: nStages = 0
    !! Number of nearest neighbors
    real(dp) :: relTol
    !! Tolerance for shell radius
    real(dp), dimension(:,:), allocatable :: Distances
    !! List of all distances in nnstages range; size (nStages, nAtoms)
    integer(int32) :: nTypes = 0
    !! Number of different atom types
    type(AtomType), dimension(:), allocatable :: Types
    !! List of types
  end type System_type

  type(System_type) :: sys
  !! System variables

  integer(int32) :: n0sc1 !< first neighbor to calculate the in-plane spin and charge current
  integer(int32) :: n0sc2 !< last neighbor to calculate the in-plane spin and charge current
  integer(int32) :: n0sc  !< Number of neighbors to calculate currents
  real(dp), dimension(3) :: pln_normal

  integer(int32), dimension(:,:), allocatable :: ia
#ifdef _GPU
  integer(int32), dimension(:,:), allocatable, device :: ia_d
  integer, allocatable, device :: nOrb_d(:)
#endif
  integer(int32), dimension(:,:), allocatable :: ia_sc

contains

  subroutine init_Hamiltk_variables(s,superCond)
  !! This subroutine builds hamiltonian variables (dimensions and conversions)
    use mod_parameters, only: dimH,dimspinAtoms,dimens,dimHsc
    implicit none
    type(System_type), intent(in) :: s
    integer(int32),    intent(in) :: superCond
    integer(int32) :: i,offsetParameter

    if(allocated(ia)) deallocate(ia)
    if(allocated(ia_sc)) deallocate(ia_sc)
    allocate(ia(4,s%nAtoms))
    allocate(ia_sc(4,s%nAtoms))

    dimH = s%Types(s%Basis(1)%Material)%nOrb
    dimens = s%Types(s%Basis(1)%Material)%nOrb * s%Types(s%Basis(1)%Material)%nOrb
  
    ia(1,1) = 1
    ia(2,1) = ia(1,1) + s%Types(s%Basis(1)%Material)%nOrb - 1
    ia(3,1) = ia(2,1) + 1
    ia(4,1) = ia(3,1) + s%Types(s%Basis(1)%Material)%nOrb - 1
    do i = 2, s%nAtoms
      dimH = dimH + s%Types(s%Basis(i)%Material)%nOrb
      dimens = dimens + s%Types(s%Basis(i)%Material)%nOrb * s%Types(s%Basis(i)%Material)%nOrb

      ia(1,i) = ia(4,i-1) + 1                                    ! Begin up 
      ia(2,i) = ia(1,i) + s%Types(s%Basis(i)%Material)%nOrb - 1  ! End up
      ia(3,i) = ia(2,i) + 1                                      ! Begin down 
      ia(4,i) = ia(3,i) + s%Types(s%Basis(i)%Material)%nOrb - 1  ! End down
    end do
    offsetParameter = dimH
    dimH = dimH*2
    dimens = dimens*4
    dimHsc  = dimH*superCond
    dimspinAtoms = 4 * s%nAtoms

    ! Superconductivity strides (need dimH)
    ia_sc(1,1) = 1
    ia_sc(2,1) = ia_sc(1,1) + 2*s%Types(s%Basis(1)%Material)%nOrb - 1
    ia_sc(3,1) = ia_sc(1,1) + dimH
    ia_sc(4,1) = ia_sc(3,1) + 2*s%Types(s%Basis(1)%Material)%nOrb - 1
    do i = 2, s%nAtoms
      ! Superconductivity block has doubled dimensions in each spin
      ia_sc(1,i) = ia_sc(2,i-1) + 1                                        ! Begin first block (electrons) 1 to 2*nOrb
      ia_sc(2,i) = ia_sc(1,i) + 2*s%Types(s%Basis(i)%Material)%nOrb - 1    ! End first block (electrons)
      ia_sc(3,i) = ia_sc(1,i) + dimH                                       ! Begin second block (holes) 1 to 2*nOrb + dimH
      ia_sc(4,i) = ia_sc(3,i) + 2*s%Types(s%Basis(i)%Material)%nOrb - 1    ! End second block (holes)
    end do

#ifdef _GPU
    if(allocated(ia_d)) deallocate(ia_d)
    allocate(ia_d(4,s%nAtoms))
    ia_d = ia
#endif

  end subroutine init_Hamiltk_variables


  subroutine initConversionMatrices(s)
  !> This subroutine mounts the conversion matrices from 4 to 2 ranks
#ifdef _GPU
    use mod_parameters, only: sigmai2i,sigmaimunu2i,sigmaijmunu2i,isigmamu2n,isigmamu2n_d,n2isigmamu
#else
    use mod_parameters, only: sigmai2i,sigmaimunu2i,sigmaijmunu2i,isigmamu2n,n2isigmamu
#endif
    implicit none
    type(System_type), intent(in) :: s
    integer :: nu, mu, i, sigma, j, kount

    do i = 1, s%nAtoms
      do sigma = 1, 4
        sigmai2i(sigma,i) = (sigma-1)*s%nAtoms + i
      end do
    end do

    !------------------------- Conversion arrays  --------------------------
    kount = 0
    do i = 1, s%nAtoms
      do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
        do nu = 1, s%Types(s%Basis(i)%Material)%nOrb
          do sigma = 1, 4
            ! kount = (sigma-1)*nAtoms*nOrb*nOrb + (i-1)*nOrb*nOrb + (mu-1)*nOrb + nu
            kount = kount + 1
            sigmaimunu2i(sigma,i,mu,nu) = kount
          end do
        end do
      end do
    end do

    kount = 0
    do i = 1, s%nAtoms
      do mu = 1, s%Types(s%Basis(i)%Material)%nOrb 
        do j = 1, s%nAtoms
          do nu = 1, s%Types(s%Basis(j)%Material)%nOrb
            do sigma = 1, 4
              ! kount = (sigma-1)*nAtoms*nAtoms*nOrb*nOrb + (i-1)*nAtoms*nOrb*nOrb + (j-1)*nOrb*nOrb + (mu-1)*nOrb + nu
              kount = kount + 1
              sigmaijmunu2i(sigma,i,j,mu,nu) = kount
            end do
          end do
        end do
      end do
    end do

    ! Conversion array from local atomic orbitals to (i,sigma,mu) to eigenvectors n
    kount = 0
    do i = 1, s%nAtoms
      do sigma = 1, 2
        do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
          ! kount = (i-1)*2*nOrb + (sigma-1)*nOrb + mu
          kount = kount + 1
          isigmamu2n(i,sigma,mu) = kount
          n2isigmamu(kount,1) = i
          n2isigmamu(kount,2) = sigma
          n2isigmamu(kount,3) = mu
        end do
      end do
    end do

#ifdef _GPU
    isigmamu2n_d = isigmamu2n
#endif

  end subroutine initConversionMatrices

  subroutine allocate_basis_variables(s)
  !! Allocates variables inside Basis derived datatype
    implicit none
    type(System_type), intent(inout) :: s
    integer :: i

    do i = 1,s%nAtoms

      if(.not.allocated(s%Basis(i)%sb)) allocate( s%Basis(i)%sb(s%Types(s%Basis(i)%Material)%nOrb2,s%Types(s%Basis(i)%Material)%nOrb2))
      s%Basis(i)%sb = cmplx(0._dp,0._dp,dp)
      if(.not.allocated(s%Basis(i)%lb)) allocate( s%Basis(i)%lb(s%Types(s%Basis(i)%Material)%nOrb2,s%Types(s%Basis(i)%Material)%nOrb2))
      s%Basis(i)%lb = cmplx(0._dp,0._dp,dp)
      if(.not.allocated(s%Basis(i)%ls)) allocate( s%Basis(i)%ls(s%Types(s%Basis(i)%Material)%nOrb2,s%Types(s%Basis(i)%Material)%nOrb2))
      s%Basis(i)%ls = cmplx(0._dp,0._dp,dp)
#ifdef _GPU
      if(.not.allocated(s%Basis(i)%sb_d)) allocate( s%Basis(i)%sb_d(s%Types(s%Basis(i)%Material)%nOrb2,s%Types(s%Basis(i)%Material)%nOrb2))
      if(.not.allocated(s%Basis(i)%lb_d)) allocate( s%Basis(i)%lb_d(s%Types(s%Basis(i)%Material)%nOrb2,s%Types(s%Basis(i)%Material)%nOrb2))
      if(.not.allocated(s%Basis(i)%ls_d)) allocate( s%Basis(i)%ls_d(s%Types(s%Basis(i)%Material)%nOrb2,s%Types(s%Basis(i)%Material)%nOrb2))
#endif
      if(.not.allocated(s%Basis(i)%lpvec)) allocate( s%Basis(i)%lpvec(s%Types(s%Basis(i)%Material)%nOrb,s%Types(s%Basis(i)%Material)%nOrb,3) )

    end do

  end subroutine allocate_basis_variables


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
