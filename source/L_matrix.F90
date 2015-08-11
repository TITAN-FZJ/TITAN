! This subroutine calculate the orbital angular momentum matrix in the cubic system of coordinates
subroutine l_matrix(Lx,Ly,Lz)
	use mod_f90_kind
	use mod_constants
	implicit none
	complex(double), dimension(9,9), intent(out) :: Lx,Ly,Lz
	complex(double), dimension(9,9) :: Lp,Lm

	Lz = zero

	Lz(2,3) = -zi
	Lz(3,2) = zi

	Lz(5,8) = 2.d0*zi
	Lz(8,5) = -2.d0*zi
	Lz(6,7) = zi
	Lz(7,6) = -zi

	Lp = zero
	Lm = zero

	Lp(2,4) = -zum
	Lp(3,4) = -zi
	Lp(4,2) = zum
	Lp(4,3) = zi

	Lp(5,6) = -zum
	Lp(5,7) = -zi

	Lp(6,5) = zum
	Lp(6,8) = -zi
	Lp(6,9) = -sq3*zi

	Lp(7,5) = zi
	Lp(7,8) = zum
	Lp(7,9) = -sq3

	Lp(8,6) = zi
	Lp(8,7) = -zum

	Lp(9,6) = sq3*zi
	Lp(9,7) = sq3

	Lm = transpose(conjg(Lp))

	Lx = 0.5d0*(Lp+Lm)
	Ly = -0.5d0*zi*(Lp-Lm)

	return
end subroutine l_matrix

! This subroutine calculate the orbital angular momentum matrix in the spin system of coordinates
subroutine lp_matrix()
	use mod_f90_kind
	use mod_parameters, only: theta, phi
	use mod_constants
	use mod_magnet
	implicit none
	complex(double), dimension(9,9) :: Lp,Lm,Lx,Ly,Lz

	Lz = zero

	Lz(2,3) = -zi
	Lz(3,2) = zi

	Lz(5,8) = 2.d0*zi
	Lz(8,5) = -2.d0*zi
	Lz(6,7) = zi
	Lz(7,6) = -zi

	Lp = zero

	Lp(2,4) = -zum
	Lp(3,4) = -zi
	Lp(4,2) = zum
	Lp(4,3) = zi

	Lp(5,6) = -zum
	Lp(5,7) = -zi

	Lp(6,5) = zum
	Lp(6,8) = -zi
	Lp(6,9) = -sq3*zi

	Lp(7,5) = zi
	Lp(7,8) = zum
	Lp(7,9) = -sq3

	Lp(8,6) = zi
	Lp(8,7) = -zum

	Lp(9,6) = sq3*zi
	Lp(9,7) = sq3

	Lm = transpose(conjg(Lp))

	Lx = 0.5d0*(Lp+Lm)
	Ly = -0.5d0*zi*(Lp-Lm)

	lxp = (Lx*cos(theta)*cos(phi))+(Ly*cos(theta)*sin(phi))-(Lz*sin(theta))
	lyp =-(Lx*sin(phi))+(Ly*cos(phi))
	lzp = (Lx*sin(theta)*cos(phi))+(Ly*sin(theta)*sin(phi))+(Lz*cos(theta))

	return
end subroutine lp_matrix