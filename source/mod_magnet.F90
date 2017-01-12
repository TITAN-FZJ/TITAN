module mod_magnet
  use mod_f90_kind
  implicit none
  integer                     :: iter                                   ! self-consistency iteration
  real(double),allocatable    :: mx(:),my(:),mz(:),mvec_cartesian(:,:),hdel(:)    ! Magnetization and exchange split delta/2
  real(double),allocatable    :: lxm(:),lym(:),lzm(:)                   ! Orbital angular momentum in lattice frame of reference
  real(double),allocatable    :: lxpm(:),lypm(:),lzpm(:)                ! Orbital angular momentum in spin frame of reference
  complex(double),allocatable :: mp(:),hdelp(:)
  complex(double),allocatable :: mm(:),hdelm(:)
  real(double),allocatable    :: mabs(:),mtheta(:),mphi(:),mvec_spherical(:,:)
  real(double),allocatable    :: labs(:),ltheta(:),lphi(:)
  real(double),allocatable    :: lpabs(:),lptheta(:),lpphi(:)
  real(double),allocatable    :: eps1(:)                                ! Center of the bands for each l - eps(Npl)
  real(double),allocatable    :: hhwx(:),hhwy(:),hhwz(:)                ! Static magnetic fields in each direction
  complex(double),allocatable :: lb(:,:,:),sb(:,:,:)                    ! Zeeman matrices
  complex(double), dimension(9,9)   :: lxp,lyp,lzp                      ! Angular momentum matrices in spin coordinate system
end module mod_magnet