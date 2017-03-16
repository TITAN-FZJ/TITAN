module mod_parameters
  use mod_f90_kind, only: double
  implicit none
  !========================================================================================!
  integer,parameter :: dmax=20 !< Maximum dimension for Npl-dependent quantities that must be read from input
  integer       :: ncp         !< Number of k-point generation vectors (S. Cunningham, Phys. Rev. B 10, 4988 (1974))
  integer       :: nkpt        !< Approximate number of k-points
  real(double)  :: eta         !< Small imaginary part included in the energy z = E + i.eta
  real(double)  :: Ef          !< Fermi energy (read from mod_tight_binding)
  real(double)  :: q(3)        !< q-vector for the dependence of response functions (not used yet)
  ! Dimensions
  integer       :: dimsigmaNpl !< Dimensions (Maybe enhance description?)
  integer       :: dim         !< Dimensions (Maybe enhance description?)
  !========================================================================================!
  real(double),allocatable  :: U(:)       !< Effective intra-site electron electron interaction
  integer       :: Utype = 2              !< Description missing
  logical       :: lhfresponses = .false. !< Use HF susceptibilities to calculate currents, disturbances and accumulations (don't renormalize)
  !========================================================================================!
  ! Lattice and surface direction
  character(len=6) :: lattice                   !< Lattice description; general(not implemented) or bcc, fcc, hcp
  integer          :: Npl,Npl_i,Npl_f,Npl_input !< Description missing.
  real(double)     :: a1(3), a2(3), a3(3)       !< Lattice unit vectors
  real(double)     :: pln_dir(3)                !< Plane vector describing the layer in the unit cell described by a1,a2,a3
  real(double)     :: a0                        !< Lattice parameter (define the units of distance in the program)
  !========================================================================================!
  character(len=1)   :: magaxis         !< Equilibrium magnetization
  character(len=1)   :: dirEfield       !< Direction of in-plane applied electric field
  real(double)       :: theta,phi       !< Euler Angles for the magnetization direction
  real(double)       :: dirEfieldvec(3) !< Direction vector of the electric field
  !========================================================================================!
  integer :: itype  !< type of calculation - defined in input file 'inputdhe'
  !========================================================================================!
  logical :: SOC                                          !< Turn on/off SOC
  logical :: llineargfsoc = .false.,llinearsoc = .false.  !< Linear SOC
  real(double) :: socscale                                !< Rescale of SOC parameter
  !========================================================================================!
  logical :: lfield !< Turn on/off static magnetic field, option to give in magnetic field in tesla
  ! Values of magnetic field in cartesian or spherical coordinates
  character(len=9),dimension(7)  :: dcfield = ["hwa      ","hwt      ","hwp      ","hwahwt   ","hwahwp   ","hwthwp   ","hwahwthwp"]
  character(len=60)              :: dc_header
  character(len=60),allocatable  :: dc_fields(:),dcprefix(:)
  real(double)     ,allocatable  :: hw_list(:,:)
  integer       :: dcfield_dependence=0,dc_count=0
  real(double)  :: hwx=0.d0,hwy=0.d0,hwz=0.d0,tesla=1.d0
  real(double)  :: hwa,hwa_i=0.d0,hwa_f=0.d0,hwa_s
  integer       :: hwa_npts=0,hwa_npt1=1,hwa_count
  real(double)  :: hwt,hwt_i=0.d0,hwt_f=0.d0,hwt_s
  integer       :: hwt_npts=0,hwt_npt1=1,hwt_count
  real(double)  :: hwp,hwp_i=0.d0,hwp_f=0.d0,hwp_s
  integer       :: hwp_npts=0,hwp_npt1=1,hwp_count
  integer       :: hw_count,total_hw_npt1
  ! Layer-resolved scale of magnetic field (including empty spheres)
  logical       :: lhwscale        = .false.
  real(double)  :: hwscale(dmax)   = 1.d0
  logical       :: lhwrotate       = .false.
  real(double)  :: hwtrotate(dmax) = 0.d0, hwprotate(dmax) = 0.d0
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
  integer      :: npts,npt1,count
  real(double) :: emin,emax,deltae
  real(double) :: qxmin,qxmax,qzmin,qzmax
  ! Number of steps to skip from the beginning
  ! (useful to get the same points after a calculation has stopped)
  integer :: skip_steps = 0, skip_steps_hw = 0
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
  logical :: ltesla         = .false.
  logical :: lcreatefiles   = .false.
  logical :: laddresults    = .false.
  logical :: lGSL           = .false.
  logical :: lslatec        = .false.
  logical :: lnojac         = .false.
  logical :: lontheflysc    = .false.
  logical :: lrotatemag     = .false.
  logical :: lnolb          = .false.
  logical :: lnodiag        = .false.
  logical :: lwriteonscreen = .false.
  logical :: lsortfiles     = .false.
  logical :: lsha           = .false.
  logical :: llgtv          = .false.
  !========================================================================================!
  ! Activate debug options
  logical :: lverbose = .false.
  logical :: ldebug   = .false.
  !========================================================================================!
  ! Longitudinal and transverse, and Spin Hall Angle calculation
  integer, dimension(:),  allocatable :: sha_longitudinal,sha_transverse ! In-plane longitudinal and transverse neighbors
  real(double), dimension(:),  allocatable :: long_cos(:),transv_cos(:)  ! In-plane longitudinal and transverse cosines
  !========================================================================================!
  ! DC voltage calculations
  logical :: lvdc = .false.
  integer :: vdcneighbor(2)=0 ! Longitudinal (1) and transverse (2) neighbor to calculate Vdc
  integer :: mvdcvector(3)    ! Mapping of the magnetization direction for the selected applied field/magnetization frame of reference
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
  ! Suffix to use on filenames (to avoid overwriting while comparing different results)
  character(len=50)  :: suffix=""
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
  ! Layers to add after set2 (maximum of 10) - Must include empty spheres on the list
  integer :: addlayers(10),naddlayers=0
  character(len=50) :: Npl_folder
  !========================================================================================!
end module mod_parameters
