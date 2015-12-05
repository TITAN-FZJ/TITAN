module mod_parameters
  use mod_f90_kind
  implicit none
  integer       :: ncp
  real(double)  :: eta,hwx,hwy,hwz
  real(double)  :: Ef,q(2)
  integer       :: dimsigmaNpl,dim
!========================================================================================!
! Effective intra-site electron electron interaction
  real(double),allocatable  :: U(:)
  integer       :: Utype
!========================================================================================!
! Lattice and surface direction
  character(len=6) :: lattice
  integer          :: Npl,Nplini,Nplfinal
! Lattice parameter (define the units of distance in the program)
  real(double)     :: a0
!========================================================================================!
! Equilibrium magnetization and direction of in-plane applied electric field:
  character(len=1)   :: magaxis,dirEfield
  real(double)       :: theta,phi       ! Euler Angles for the magnetization direction
  real(double)       :: dirEfieldvec(3) ! Direction vector of the electric field
!========================================================================================!
! type of calculation - defined in input file 'inputdhe'
  integer :: itype
!========================================================================================!
! Turn on/off SOC
  logical :: SOC
! Rescale of SOC parameter
  real(double) :: socscale
!========================================================================================!
! Number of parts to divide energy integral I1+I2 and I3
  integer      :: pn1,pn2,pnt
  integer      :: parts,parts3
  integer      :: n1gl,n3gl
  real(double) :: tol
!========================================================================================!
! Band structure
  character(len=2)  :: kdirection
!========================================================================================!
! Number of points and interval of energy/wave vector/position calculations
  integer      :: npts,npt1
  real(double) :: emin,emax,deltae
  real(double) :: qxmin,qxmax,qzmin,qzmax
!========================================================================================!
! Conversion arrays
  integer,allocatable :: sigmaimunu2i(:,:,:,:),sigmaijmunu2i(:,:,:,:,:),sigmai2i(:,:)
!========================================================================================!
! Run options are stored in this string
  character(len=100)          :: runoptions
  real(double)                :: ry2ev              ! Optional conversion of ry to eV
!========================================================================================!
! Activate debug options
  logical :: verbose
  logical :: idebug
!========================================================================================!
! n0sc1 - first neighbor to calculate the in-plane spin and charge current
! n0sc2 - last neighbor to calculate the in-plane spin and charge current
! n0sc  - Number of neighbors to calculate currents
  integer :: n0sc1,n0sc2,n0sc
!========================================================================================!
! Current renormalization
  logical :: renorm
  integer :: renormnb
!========================================================================================!
! Skip self-consistency
  logical :: skipsc
! Give a file to start self-consistency
  character(len=200) :: scfile
  integer            :: scfileini
!========================================================================================!
! Choose between tight-binding (T) or orthogonal (O) DFT parameters
  character(len=1)  :: dfttype
!========================================================================================!
! Layer conversion
  integer, dimension(:), allocatable         :: mmlayer
! Number and list of magnetic layers
  integer :: nmaglayers
  integer, dimension(:), allocatable       :: mmlayermag
! Layer type: 1 - Empty Sphere; 2 - Magnetic ; 3 - Bulk ; 0 - Other
  integer, dimension(:), allocatable       :: layertype
!========================================================================================!
! Hostname of rank
  character(len=50) :: host
!========================================================================================!
! Set of tight-binding parameters to be used
! in the first half (set1) and second half (set2) of the slab
! NOTE: the Fermi energy is obtained from the first half.
  integer :: set1,set2
!========================================================================================!
end module mod_parameters