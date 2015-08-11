! subroutine lb_matrix()
! 	use mod_f90_kind
! 	use mod_constants
! 	use mod_parameters
! 	use mod_magnet
! 	implicit none

! ! 	real(double),intent(in)	:: hhwx,hhwy,hhwz
! ! 	real(double),intent(in)	:: theta,phi
! 	real(double)			:: hhwxp,hhwyp,hhwzp
! ! 	complex(double), dimension(18,18),intent(out) :: lb
! 	complex(double), dimension(9,9) :: Lx,Ly,Lz,lbsigma

! 	call L_matrix(Lx,Ly,Lz)

! 	hhwxp = (hhwx*cos(theta)*cos(phi)) - (hhwy*sin(phi)) + (hhwz*sin(theta)*cos(phi))
! 	hhwyp = (hhwx*cos(theta)*sin(phi)) + (hhwy*cos(phi)) + (hhwz*sin(theta)*sin(phi))
! 	hhwzp = -(hhwx*sin(theta)) + (hhwz*cos(theta))

! !	The minus sign takes into account the fact that we are considering negative
! !	external fields to get the peak in positive energies
! 	lbsigma = - 0.5d0*(Lx*hhwxp + Ly*hhwyp + Lz*hhwzp)

! 	lb = zero
! 	lb( 1: 9, 1: 9) = lbsigma(:,:)
! 	lb(10:18,10:18) = lbsigma(:,:)

! 	return
! end subroutine lb_matrix

subroutine lb_matrix()
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_magnet
	implicit none

! 	complex(double), dimension(18,18),intent(out) :: lb
	complex(double), dimension(9,9) :: lbsigma

! 	call lp_matrix()

!	The minus sign takes into account the fact that we are considering negative
!	external fields to get the peak in positive energies
	lbsigma = - 0.5d0*(lxp*hhwx + lyp*hhwy + lzp*hhwz)

	lb = zero
	lb( 1: 9, 1: 9) = lbsigma(:,:)
	lb(10:18,10:18) = lbsigma(:,:)

	return
end subroutine lb_matrix