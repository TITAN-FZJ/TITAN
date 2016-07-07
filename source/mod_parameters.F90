module mod_parameters
  use mod_f90_kind
  implicit none
  integer       :: ncp
  real(double)  :: eta
  real(double)  :: Ef,q(3)
  integer       :: dimsigmaNpl,dim
!========================================================================================!
! Effective intra-site electron electron interaction
  real(double),allocatable  :: U(:)
  integer       :: Utype = 2
! Use HF susceptibilities to calculate currents, disturbances and accumulations (don't renormalize)
  logical       :: lhfresponses = .false.
!========================================================================================!
! Lattice and surface direction
  character(len=6) :: lattice
  integer          :: Npl,Npl_i,Npl_f
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
! Linear SOC
  logical :: llineargfsoc = .false.,llinearsoc = .false.
! Rescale of SOC parameter
  real(double) :: socscale
!========================================================================================!
! Turn on/off static magnetic field, option to give in magnetic field in tesla
  logical :: lfield
! Values of magnetic field in cartesian or spherical coordinates
  character(len=9),dimension(7)  :: dcfield = ["hwa      ","hwt      ","hwp      ","hwahwt   ","hwahwp   ","hwthwp   ","hwahwthwp"]
  character(len=60)              :: dc_header,dcprefix
  character(len=60),allocatable  :: dc_fields(:)
  real(double)     ,allocatable  :: hw_list(:,:)
  integer       :: dcfield_dependence=0
  real(double)  :: hwx=0.d0,hwy=0.d0,hwz=0.d0,tesla=1.d0
  real(double)  :: hwa,hwa_i=0.d0,hwa_f=0.d0,hwa_s
  integer       :: hwa_npts,hwa_npt1=1,hwa_count
  real(double)  :: hwt,hwt_i=0.d0,hwt_f=0.d0,hwt_s
  integer       :: hwt_npts,hwt_npt1=1,hwt_count
  real(double)  :: hwp,hwp_i=0.d0,hwp_f=0.d0,hwp_s
  integer       :: hwp_npts,hwp_npt1=1,hwp_count
  integer       :: hw_count,total_hw_npt1
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
! Number of steps to skip from the beginning
! (useful to get the same points after a calculation has stopped)
  integer :: skip_steps = 0
!========================================================================================!
! Conversion arrays
  integer,allocatable :: sigmaimunu2i(:,:,:,:),sigmaijmunu2i(:,:,:,:,:),sigmai2i(:,:)
!========================================================================================!
! Run options are stored in this string
  character(len=200)          :: runoptions
  real(double)                :: ry2ev=1.d0              ! Optional conversion of ry to eV
!========================================================================================!
! Logical variables for runoptions
  logical :: lkpoints     = .false.
  logical :: ltesla       = .false.
  logical :: lcreatefiles = .false.
  logical :: laddresults  = .false.
  logical :: lGSL         = .false.
  logical :: lslatec      = .false.
  logical :: lnojac       = .false.
  logical :: lontheflysc  = .false.
  logical :: lrotatemag   = .false.
!========================================================================================!
! Activate debug options
  logical :: lverbose = .false.
  logical :: ldebug = .false.
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
  logical :: skipsc,lselfcon
! Give a file to start self-consistency
  character(len=200) :: scfile
!========================================================================================!
! Output file
  integer            :: outputunit=123456789,outputunit_loop
  character(len=200) :: outputdhe,outputdhe_loop
!========================================================================================!
! Choose between tight-binding (T) or orthogonal (O) DFT parameters
  character(len=1)  :: dfttype
!========================================================================================!
! Layer conversion
  integer, dimension(:), allocatable       :: mmlayer
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