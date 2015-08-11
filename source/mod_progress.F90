module mod_progress
	use mod_f90_kind
	implicit none
	character(len=130) :: progbar
	real(double)			 :: start_time,elapsed_time
	integer            :: prog
end module mod_progress
