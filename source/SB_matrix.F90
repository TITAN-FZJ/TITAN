subroutine sb_matrix()
	use mod_f90_kind
	use mod_constants
	use mod_magnet
	implicit none
	integer	:: mu,nu

	sb = zero
	do mu=1,9
		nu=mu+9
		sb(mu,mu) = -hhwz
		sb(nu,nu) =  hhwz
		sb(mu,nu) = -(hhwx-zi*hhwy)
		sb(nu,mu) = -(hhwx+zi*hhwy)
	end do

	return
end subroutine sb_matrix