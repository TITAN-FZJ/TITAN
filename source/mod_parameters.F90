module mod_parameters
  use mod_f90_kind, only: double
  implicit none
  !========================================================================================!
  integer,parameter :: dmax=20
  !! Maximum dimension for Npl-dependent quantities that must be read from input (obsolete?)
  real(double)  :: eta
  !! Small imaginary part included in the energy of GF or susceptibility z = w + i.eta
  real(double)  :: etap
  !! Small imaginary part included in the energy of GF  z = E + i.eta'
  real(double)  :: q(3)
  !! q-vector for the dependence of response functions (not used yet)
  ! Dimension variables:
  integer       :: dimspinAtoms
  !! Dimension: 4 spins x number of atoms in the unit cell
  integer       :: dim
  !! Dimension: 4 spins x number of atoms in the unit cell x number of orbitals^2

  integer, dimension(3) :: kp_in
  !! Number of k-points in each direction
  integer :: kptotal_in
  !! Total number of k-points

  !========================================================================================!
  real(double), allocatable  :: U(:)
  !! Effective intra-site electron electron interaction
  logical       :: lhfresponses = .false.
  !! Use HF susceptibilities to calculate currents, disturbances and accumulations (don't renormalize)
  !========================================================================================!
  ! Lattice and surface direction
  character(len=6) :: latticeName
  !! Lattice description; general or bcc, fcc, hcp, cubic
  !integer          :: Npl,Npl_i,Npl_f,Npl_input
  !integer          :: Npl_total
  !! Obsolete?
  logical          :: bulk = .false.
  !! Flag turning on/off bulk calculations, default: .false., not used yet
  !========================================================================================!
  character(len=50)  :: magbasis = ""
  !! Basis to give initial magnetization in 'initialmag' file
  real(double),allocatable       :: initialmag(:,:)
  !! Initial guess for magnetization
  real(double)       :: theta=0.d0,phi=0.d0
  !! Euler Angles for the magnetization frame of reference
  !========================================================================================!
  integer :: itype
  !! type of calculation - defined in input file 'input'

  !========================================================================================!
  ! Number of parts to divide energy integral I1+I2 and I3
  !========================================================================================!
  character(len=5), dimension(:), allocatable :: bands
  !! Band structure
  integer :: band_cnt
  !========================================================================================!
  ! Number of points and interval of energy/wave vector/position calculations
  integer      :: npts,npt1,count
  !! Number of energy points (+1) and counter
  real(double) :: emin,emax,deltae
  !! Minimum energy, maximum energy, and step size
  real(double) :: qxmin,qxmax,qzmin,qzmax
  !! Minimum and maximum q-vector
  integer :: skip_steps = 0
  !! Number of steps to skip from the beginning
  !! (useful to get the same points after a calculation has stopped)
  !========================================================================================!
  integer,allocatable :: sigmaimunu2i(:,:,:,:),sigmaijmunu2i(:,:,:,:,:),sigmai2i(:,:)
  !! Conversion arrays
  !========================================================================================!
  character(len=200)          :: runoptions
  !! Run options are stored in this string
  real(double)                :: ry2ev=1.d0
  !! Optional conversion of ry to eV
  !========================================================================================!
  ! Logical variables for runoptions
  logical :: lkpoints       = .false.
  logical :: lpositions     = .false.
  logical :: ltesla         = .false.
  logical :: lcreatefiles   = .false.
  logical :: lcreatefolders = .false.
  logical :: laddresults    = .false.
  logical :: lnolb          = .false.
  logical :: lnodiag        = .false.
  logical :: lwriteonscreen = .false.
  logical :: lsortfiles     = .false.
  logical :: lsha           = .false.
  logical :: llgtv          = .false.
  logical :: lcheckjac      = .false.
  !========================================================================================!
  ! Activate debug options
  logical :: lverbose = .false.
  logical :: ldebug   = .false.
  !========================================================================================!
  ! Longitudinal and transverse, and Spin Hall Angle calculation
  integer, dimension(:),  allocatable :: sha_longitudinal,sha_transverse
  !! In-plane longitudinal and transverse neighbors
  real(double), dimension(:),  allocatable :: long_cos(:),transv_cos(:)
  !! In-plane longitudinal and transverse cosines
  integer :: longitudinal_neighbors
  !! Number of longitudinal neighbors
  integer :: transverse_neighbors
  !! Number of transverse neighbors

  !! String for HF Responses
  !========================================================================================!
  ! n0sc1 - first neighbor to calculate the in-plane spin and charge current
  ! n0sc2 - last neighbor to calculate the in-plane spin and charge current
  ! n0sc  - Number of neighbors to calculate currents
  ! integer :: n0sc1,n0sc2,n0sc
  !========================================================================================!
  ! Current renormalization (obsolete)
  logical :: renorm
  integer :: renormnb
  !========================================================================================!
  character(len=500)  :: missing_files=""
  !! Variable to store missing filenames
  !========================================================================================!
  ! Suffix to use on filenames (to avoid overwriting while comparing different results)
  !========================================================================================!


  !========================================================================================!
  character(len=1)  :: dfttype
  !! Choose between tight-binding (T) or orthogonal (O) DFT parameters (obsolete?)
  !========================================================================================!
  ! Set of tight-binding parameters to be used
  ! in the first half (set1) and second half (set2) of the slab
  ! NOTE: the Fermi energy is obtained from the first half.
  integer :: set1,set2
  ! Layers to add after set2 (maximum of 10) - Must include empty spheres on the list
  integer :: addlayers(10),naddlayers=0
  !========================================================================================!
  integer :: offset = 0

  integer :: parField = 1
  integer :: parFreq = 1

  type :: Filename
     integer :: unit=123456789
     character(len=200) :: file

     integer :: unit_loop
     character(len=200) :: file_loop

     character(len=50) :: BField   = "" ! hwa, hwt, hwp, ltesla, lnolb, lhwscale, lhwrotate,
     character(len=50) :: dcBField = "" ! dc special case of fieldpart

     character(len=50) :: EField = ""   ! EFt, EFp

     character(len=50) :: suffix = ""
     character(len=50) :: Sites = ""

     character(len=50) :: SOC = ""
     character(len=1)  :: SOCchar = ""

     character(len=50) :: Energy = "" ! parts, parts3

     character(len=3)  :: hfr = ""

     character(len=50) :: info = ""
  end type Filename
  type(Filename) :: output
end module mod_parameters
