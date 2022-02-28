module mod_parameters
  use mod_kind, only: dp, int32, int64
  implicit none
  !========================================================================================!
  integer(int32),parameter :: dmax=20
  !! Maximum dimension for quantities that must be read from input
  real(dp)  :: eta
  !! Small imaginary part included in the energy of GF or susceptibility z = w + i.eta
  real(dp)  :: etap
  !! Small imaginary part included in the energy of GF  z = E + i.eta'
  real(dp) :: addelectrons
  !! Electrons to add (or remove, when negative) for the whole system (default = 0 set in mod_io)
  !! Only for the case of Un==0 (otherwise, Un requires site- and orbital-dependent quantity)
  real(dp)  :: q(3)
  !! q-vector for the dependence of response functions
  ! Dimension variables:
  integer(int32) :: dimH
  !! Dimension of the Hamiltonian: 2 spins x number of atoms in the unit cell
  integer(int32) :: dimHsc
  !! Dimension of the Hamiltonian taking into account superconductivity (may be dimH or 2*dimH)
  integer(int32) :: dimspinAtoms
  !! Dimension: 4 spins x number of atoms in the unit cell
  integer(int32) :: dimens
  !! Dimension: 4 spins x number of atoms in the unit cell x number of orbitals^2

  integer(int32), dimension(3) :: kp_in, qp_in 
  !! Number of k-points in each direction
  integer(int64) :: kptotal_in, qptotal_in
  !! Total number of k-points

  integer(int32) :: tbmode
  !! Tight-binding mode: (1) Slater-Koster, (2) DFT
  integer(int32) :: fermi_layer
  !! Which site will be used for fermi level (Maybe remove it and read from outside?)
  real(dp) :: Ef_overwrite
  logical  :: lEf_overwrite = .false.
  !! Overwrite Fermi energy variables (that is set immediately before self-consistency)
  character(len=20), dimension(:), allocatable :: layers
  !! Number of layers (Obsolete?)
  !========================================================================================!
  integer(int32) :: cluster_layers = 2
  !! Number of cells around the origin for the cluster generation of the real space Jij
  !========================================================================================!
  logical        :: lhfresponses = .false.
  !! Use HF susceptibilities to calculate currents, disturbances and accumulations (don't renormalize)
  !========================================================================================!
  real(dp)       :: theta=0._dp,phi=0._dp
  !! Euler Angles for the magnetization frame of reference
  !========================================================================================!
  integer(int32) :: itype
  !! type of calculation - defined in input file 'input'
  !========================================================================================!
  character(len=10), dimension(:), allocatable :: bands
  !! Band structure points
  integer(int32)    :: band_cnt
  !! Number of points along the loop path
  character(len=20) :: kdirection
  !! Path of symmetry points followed in the Brillouin Zone
  real(dp), allocatable :: partial_length(:)
  !! Length of each segment on the Brillouin zone
  !========================================================================================!
  integer(int32)    :: nEner,nEner1
  ! Number of points of energy (frequency) loops
  integer(int64)    :: kount
  ! Counter for frequency loops
  real(dp) :: emin,emax,deltae
  !! Minimum and maximum energy (frequency), and step size
  integer(int32)    :: skip_steps = 0
  !! Number of steps to skip from the beginning
  !! (useful to get the same points after a calculation has stopped)
  integer(int32)    :: nQvec,nQvec1
  !! Number of points of wave vector loops, and step size (band structure and susceptibility)
  real(dp) :: deltak
  !! Step size of wave vector loops (band structure and susceptibility)
  real(dp), dimension(:,:), allocatable :: band_points
  !! Band points used in wave vector loop
  real(dp), allocatable :: kpoints(:,:)
  !! Kpoints in the wave vector loop
  character(len=50)  :: qbasis = ""
  !! Basis to use on kpoints given in kbands file. Default: reciprocal lattice vectors
  character(len=400) :: bsfile, wsfile
  !! Filename for band structure and weights calculations
  !========================================================================================!
  integer(int32), allocatable :: sigmaimunu2i(:,:,:,:),sigmaijmunu2i(:,:,:,:,:),sigmai2i(:,:),isigmamu2n(:,:,:), n2isigmamu(:,:)
  !! Conversion arrays
#ifdef _GPU
  integer(int32), allocatable, device :: isigmamu2n_d(:,:,:)
  !! Conversion arrays on GPUs
#endif
  !========================================================================================!
  character(len=200)          :: runoptions
  !! Run options are stored in this string
  real(dp)                :: ry2ev=1._dp
  !! Optional conversion of ry to eV
  !========================================================================================!
  ! Logical variables for runoptions
  logical :: lkpoints        = .false.
  logical :: lpositions      = .false.
  logical :: ltesla          = .false.
  logical :: lcreatefiles    = .false.
  logical :: lcreatefolders  = .false.
  logical :: laddresults     = .false.
  logical :: lnolb           = .false.
  logical :: lnodiag         = .false.
  logical :: lwriteonscreen  = .false.
  logical :: lsortfiles      = .false.
  logical :: lsha            = .false.
  logical :: llgtv           = .false.
  logical :: lcheckjac       = .false.
  logical :: lsimplemix      = .false.
  logical :: leigenstates    = .false.
  logical :: lprintfieldonly = .false.
  logical :: lfixEf          = .false.
  !========================================================================================!
  ! Activate debug options
  logical :: lverbose = .false.
  logical :: ldebug   = .false.
  !========================================================================================!
  ! Longitudinal and transverse, and Spin Hall Angle calculation
  integer(int32), dimension(:),  allocatable :: sha_longitudinal,sha_transverse
  !! In-plane longitudinal and transverse neighbors
  real(dp), dimension(:),  allocatable :: long_cos(:),transv_cos(:)
  !! In-plane longitudinal and transverse cosines
  integer(int32) :: longitudinal_neighbors
  !! Number of longitudinal neighbors
  integer(int32) :: transverse_neighbors
  !! Number of transverse neighbors

  !! String for HF Responses
  !========================================================================================!
  ! n0sc1 - first neighbor to calculate the in-plane spin and charge current
  ! n0sc2 - last neighbor to calculate the in-plane spin and charge current
  ! n0sc  - Number of neighbors to calculate currents
  ! integer :: n0sc1,n0sc2,n0sc
  !========================================================================================!
  ! Current renormalization (obsolete)
  logical   :: renorm
  integer(int32) :: renormnb
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
  integer(int32) :: set1,set2
  ! Layers to add after set2 (maximum of 10) - Must include empty spheres on the list
  integer(int32) :: addlayers(10),naddlayers=0
  !========================================================================================!

  integer(int32) :: parField = 1
  integer(int32) :: parFreq = 1

  type :: Filename
     integer :: unit=123456789
     character(len=200) :: file

     integer :: unit_loop
     character(len=200) :: file_loop

     character(len=50) :: BField   = "" ! hwa, hwt, hwp, ltesla, lnolb, lhwscale, lhwrotate,
     character(len=50) :: dcBField = "" ! dc special case of fieldpart

     character(len=50) :: EField   = ""   ! EFt, EFp

     character(len=50) :: suffix = ""
     character(len=50) :: Sites  = ""

     character(len=50) :: SOC = ""
     character(len=1)  :: SOCchar = ""

     character(len=50) :: Energy = "" ! parts, parts3

     character(len=3)  :: hfr = ""

     character(len=50) :: info = ""
     ! Strings for time-dependent calculations
     character(len=100) :: time_field = ""
     character(len=50),dimension(:), allocatable :: observable
  end type Filename
  type(Filename) :: output
  !! Strings to output filenames

  character(len=500) :: arg = ""
  !! String to store type of calculation (DEBUG, GPU)

end module mod_parameters
