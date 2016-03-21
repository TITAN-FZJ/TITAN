! Spin Zeeman hamiltonian
subroutine sb_matrix()
	use mod_f90_kind
	use mod_constants
	use mod_magnet
	implicit none
	integer	:: mu,nu

!	The minus sign takes into account the fact that we are considering negative
!	external fields to get the peak in positive energies
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