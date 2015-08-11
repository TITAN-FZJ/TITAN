module mod_f90_kind
	implicit none

	intrinsic kind,selected_int_kind,selected_real_kind	! we use these,
	private kind,selected_int_kind,selected_real_kind	! but do not force
														! them on the user.
!
!	Indicator that the KIND= is not available for this compiler/host
	integer, parameter :: not_available = -1
!
!	Real and Complex numbers
!		Single precision
	integer, parameter :: single  = kind(0.0)
!		Double precision
	integer, parameter :: double  = kind(0.0d0)
!		Quadruple precision
	integer, parameter :: quad    = selected_real_kind(p=30)
!
!	Integers numbers
!		Single byte integer
	integer, parameter :: int8    = selected_int_kind(2)
!		Two byte integer
	integer, parameter :: int16   = selected_int_kind(4)
!		Four byte integer
	integer, parameter :: int32   = selected_int_kind(9)
!		Eight byte integer
	integer, parameter :: int64   = selected_int_kind(18)
!
!	Logical values
!		Single byte logical
	integer, parameter :: byte    = 1
!		Two byte logical
	integer, parameter :: twobyte = 2
!		Four byte logical
	integer, parameter :: word    = kind(.TRUE.)
!
!		Character itype
!		Normal single byte character (ASCII sequence)
	integer, parameter :: ascii   = kind('x')
end module mod_f90_kind