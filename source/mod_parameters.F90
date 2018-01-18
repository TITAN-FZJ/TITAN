module mod_parameters
  use mod_f90_kind, only: double
  implicit none
  !========================================================================================!
  integer,parameter :: dmax=20 !< Maximum dimension for Npl-dependent quantities that must be read from input
!  integer       :: nkpt         !< Number of k-point generation vectors (S. Cunningham, Phys. Rev. B 10, 4988 (1974))
  !XXX: integer       :: nkpt        !< Number of generated k-points
  real(double)  :: eta         !< Small imaginary part included in the energy z = E + i.eta
  integer :: total_nkpt
  real(double)  :: q(3)        !< q-vector for the dependence of response functions (not used yet)
  ! Dimensions
  integer       :: dimsigmaNpl !< Dimensions (Maybe enhance description?)
  integer       :: dim         !< Dimensions (Maybe enhance description?)
  !========================================================================================!
  real(double), allocatable  :: U(:)       !< Effective intra-site electron electron interaction
  integer       :: Utype = 2              !< Description missing
  logical       :: lhfresponses = .false. !< Use HF susceptibilities to calculate currents, disturbances and accumulations (don't renormalize)
  !========================================================================================!
  ! Lattice and surface direction
  character(len=6) :: latticeName                   !< Lattice description; general or bcc, fcc, hcp, cubic
  !integer          :: Npl,Npl_i,Npl_f,Npl_input !< Description missing.
  !integer          :: Npl_total
  logical          :: bulk = .false.            !< Flag turning on/off bulk calculations, default: .false., not used yet
  !========================================================================================!
  integer            :: magaxis                   !< Initial guess for magnetization
  real(double)       :: magaxisvec(3)
  real(double)       :: theta=0.d0,phi=0.d0       !< Euler Angles for the magnetization frame of reference
  !========================================================================================!
  integer :: itype  !< type of calculation - defined in input file 'input'

  !========================================================================================!
  ! Number of parts to divide energy integral I1+I2 and I3
  !========================================================================================!
  ! Band structure
  character(len=5), dimension(:), allocatable :: bands
  integer :: band_cnt
  !========================================================================================!
  ! Number of points and interval of energy/wave vector/position calculations
  integer      :: npts,npt1,count
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
  logical :: ltestcharge    = .false.
  !========================================================================================!
  ! Activate debug options
  logical :: lverbose = .false.
  logical :: ldebug   = .false.
  !========================================================================================!
  ! Longitudinal and transverse, and Spin Hall Angle calculation
  integer, dimension(:),  allocatable :: sha_longitudinal,sha_transverse ! In-plane longitudinal and transverse neighbors
  real(double), dimension(:),  allocatable :: long_cos(:),transv_cos(:)  ! In-plane longitudinal and transverse cosines
  integer :: longitudinal_neighbors
  integer :: transverse_neighbors

  character(len=3)    :: hfr = ""
  !! String for HF Responses
  !========================================================================================!
  ! n0sc1 - first neighbor to calculate the in-plane spin and charge current
  ! n0sc2 - last neighbor to calculate the in-plane spin and charge current
  ! n0sc  - Number of neighbors to calculate currents
  ! integer :: n0sc1,n0sc2,n0sc
  !========================================================================================!
  ! Current renormalization
  logical :: renorm
  integer :: renormnb
  !========================================================================================!
  ! Variable to store missing filenames
  character(len=500)  :: missing_files=""
  !========================================================================================!
  ! Suffix to use on filenames (to avoid overwriting while comparing different results)
  character(len=50)  :: suffix=""
  !========================================================================================!
  ! Output file
  integer            :: outputunit=123456789,outputunit_loop
  character(len=200) :: outputfile, outputfile_loop

  character(len=50) :: fieldpart   = "" ! hwa, hwt, hwp, ltesla, lnolb, lhwscale, lhwrotate,
  character(len=50) :: dcfieldpart = "" ! dc special case of fieldpart

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
  ! Set of tight-binding parameters to be used
  ! in the first half (set1) and second half (set2) of the slab
  ! NOTE: the Fermi energy is obtained from the first half.
  integer :: set1,set2
  ! Layers to add after set2 (maximum of 10) - Must include empty spheres on the list
  integer :: addlayers(10),naddlayers=0
  character(len=50) :: strSites
  !========================================================================================!
  integer :: offset = 0

  integer :: parField = 1
  integer :: parFreq = 1

end module mod_parameters
