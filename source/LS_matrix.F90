subroutine ls_matrix()
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_tight_binding, only: ls
	implicit none
!	complex(double), dimension(18,18), intent(out) :: ls
!	complex(double), dimension(5,5) :: l_x,l_y,l_z

	! the spin-orbit matrix
	ls = zero

	! diagonal in spin
	! p-block
! 	ls( 2, 3) = -0.5d0*zi*cos(theta)
! 	ls( 2, 4) =  0.5d0*zi*sin(theta)*sin(phi)
! 	ls( 3, 2) =  conjg(ls(2,3))
! 	ls( 3, 4) = -0.5d0*zi*sin(theta)*cos(phi)
! 	ls( 4, 2) =  conjg(ls(2,4))
! 	ls( 4, 3) =  conjg(ls(3,4))

! 	ls(11:13,11:13) = -ls(2:4,2:4)

	! d-block
	ls( 5, 6) =  0.5d0*zi*sin(theta)*sin(phi)
	ls( 5, 7) = -0.5d0*zi*sin(theta)*cos(phi)
	ls( 5, 8) =  zi*cos(theta)
	ls( 6, 5) =  conjg(ls(5,6))
	ls( 6, 7) =  0.5d0*zi*cos(theta)
	ls( 6, 8) = -0.5d0*zi*sin(theta)*cos(phi)
	ls( 6, 9) = -0.5d0*zi*sq3*sin(theta)*cos(phi)
	ls( 7, 5) =  conjg(ls(5,7))
	ls( 7, 6) =  conjg(ls(6,7))
	ls( 7, 8) = -0.5d0*zi*sin(theta)*sin(phi)
	ls( 7, 9) =  0.5d0*zi*sq3*sin(theta)*sin(phi)
	ls( 8, 5) =  conjg(ls(5,8))
	ls( 8, 6) =  conjg(ls(6,8))
	ls( 8, 7) =  conjg(ls(7,8))
	ls( 9, 6) =  conjg(ls(6,9))
	ls( 9, 7) =  conjg(ls(7,9))

	ls(14:18,14:18) = -ls(5:9,5:9)

	! off-diagonal in spin
	! p-block
! 	ls( 2,12) =  0.5d0*zi*sin(theta)
! 	ls( 2,13) =  0.5d0*(cos(phi)+zi*cos(theta)*sin(phi))
! 	ls( 3,11) = -0.5d0*zi*sin(theta)
! 	ls( 3,13) =  0.5d0*(sin(phi)-zi*cos(theta)*cos(phi))
! 	ls( 4,11) = -0.5d0*(cos(phi)+zi*cos(theta)*sin(phi))
! 	ls( 4,12) = -0.5d0*(sin(phi)-zi*cos(theta)*cos(phi))

! 	ls(11:13,2:4) = transpose(conjg(ls(2:4,11:13)))

	! d-block
	ls( 5,15) =  0.5d0*(cos(phi)+zi*cos(theta)*sin(phi))
	ls( 5,16) =  0.5d0*(sin(phi)-zi*cos(theta)*cos(phi))
	ls( 5,17) = -zi*sin(theta)
	ls( 6,14) = -0.5d0*(cos(phi)+zi*cos(theta)*sin(phi))
	ls( 6,16) = -0.5d0*zi*sin(theta)
	ls( 6,17) =  0.5d0*(sin(phi)-zi*cos(theta)*cos(phi))
	ls( 6,18) =  0.5d0*sq3*(sin(phi)-zi*cos(theta)*cos(phi))
	ls( 7,14) =  -0.5d0*(sin(phi)-zi*cos(theta)*cos(phi))
	ls( 7,15) =  0.5d0*zi*sin(theta)
	ls( 7,17) = -0.5d0*(cos(phi)+zi*cos(theta)*sin(phi))
	ls( 7,18) =  0.5d0*sq3*(cos(phi)+zi*cos(theta)*sin(phi))
	ls( 8,14) =  zi*sin(theta)
	ls( 8,15) = -0.5d0*(sin(phi)-zi*cos(theta)*cos(phi))
	ls( 8,16) =  0.5d0*(cos(phi)+zi*cos(theta)*sin(phi))
	ls( 9,15) = -0.5d0*sq3*(sin(phi)-zi*cos(theta)*cos(phi))
	ls( 9,16) = -0.5d0*sq3*(cos(phi)+zi*cos(theta)*sin(phi))

	ls(14:18,5:9) = transpose(conjg(ls(5:9,14:18)))

! 	l_x = ls(14:18, 5: 9)+ls( 5: 9,14:18)
! 	l_y = -zi*(ls(14:18, 5: 9)-ls( 5: 9,14:18))
! 	l_z = 2.d0*ls( 5: 9, 5: 9)

	return
end subroutine ls_matrix