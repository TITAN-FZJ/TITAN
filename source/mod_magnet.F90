module mod_magnet
	use mod_f90_kind
	implicit none
	integer											:: iter 							! self-consistency iteration
	real(double),allocatable 		:: mz(:),hdel(:)			! Magnetization and exchange split delta/2
	real(double),allocatable 		:: lxm(:),lym(:),lzm(:),lxpm(:),lypm(:),lzpm(:)
	complex(double),allocatable :: mp(:),hdelp(:)
	complex(double),allocatable :: mm(:),hdelm(:)
	real(double),allocatable    :: mabs(:),mtheta(:),mphi(:)
	real(double),allocatable    :: labs(:),ltheta(:),lphi(:)
	real(double),allocatable    :: lpabs(:),lptheta(:),lpphi(:)
	real(double),allocatable 		:: eps1(:) 				  	! Center of the bands for each l - eps(Npl)
	real(double)								:: hhwx,hhwy,hhwz 		! Magnetic fields in each direction
	complex(double), dimension(18,18) :: lb,sb 				! Zeeman matrices
	complex(double), dimension(9,9) 	:: lxp,lyp,lzp 	! Angular momentum matrices in spin coordinate system
end module mod_magnet