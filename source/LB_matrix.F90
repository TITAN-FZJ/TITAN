! Orbital Zeeman hamiltonian
subroutine lb_matrix()
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_magnet
	implicit none

! 	complex(double), dimension(18,18),intent(out) :: lb
	complex(double), dimension(9,9) :: lbsigma

	lb = zero

	if(lnolb) return

!	There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
! to take into account the fact that we are considering negative
!	external fields to get the peak in positive energies
	lbsigma = 0.5d0*(lxp*hhwx + lyp*hhwy + lzp*hhwz)

	lb( 1: 9, 1: 9) = lbsigma(:,:)
	lb(10:18,10:18) = lbsigma(:,:)

	return
end subroutine lb_matrix
