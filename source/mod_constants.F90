module mod_constants
	use mod_f90_kind
	implicit none
	real(double) :: pi, tpi, sq2, hsq2, sq3
	complex(double), parameter :: zero=(0.d0,0.d0), zum=(1.d0,0.d0), zi=(0.d0,1.d0)
	complex(double)	:: identorb18(18,18),identorb9(9,9)
end module mod_constants