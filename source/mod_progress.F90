module mod_progress
	use mod_f90_kind
	implicit none
	character(len=130) :: progbar
	character(len=1)	 :: spiner(4) = ['|','/','-','\']
	real(double)			 :: start_program,start_time,elapsed_time
	integer            :: prog
end module mod_progress
