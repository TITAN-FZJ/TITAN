module mod_generate_epoints
	use mod_f90_kind
	use mod_parameters, only: parts, parts3, eta, Ef, n1gl, n3gl, pn1, pn2, pnt
	implicit none
!========================================================================================!
!	Energy integration variables
!	Number of parts to divide energy integral I1+I2 and I3
	integer :: nepoints
	real(double),allocatable      :: x2(:),p2(:)
	real(double),allocatable			:: y(:),wght(:)

contains
	!-- Generating integration points of the complex energy integral --
	subroutine generate_imag_epoints()
		implicit none
		integer												:: i,j,k,AllocateStatus
		real(double)									:: e1,e2,et1,xx
		real(double),allocatable			:: x1(:),p1(:)

		allocate( x1(pn1),p1(pn1),y(pn1),wght(pn1), STAT = AllocateStatus )
		do k=parts,1,-1
			if(parts.ne.1)then
				if(k.eq.1) then
					e1 = 10.d0**(-2)
					e2 = 1.d0
				else if(k.ne.parts) then
					e1 = 10.d0**(-k-1)
					e2 = 10.d0**(-k)
				else
					e1 = 0.d0
					e2 = 10.d0**(-k)
				end if
			else
				e1 = 0.d0
				e2 = 1.d0
			end if
			i = (k-1)*n1gl+1
			j = k*n1gl
			call gauleg(e1,e2,x1(i:j),p1(i:j),n1gl)
		end do
		et1 = 1.d0+eta
		do k=1,pn1
			y(k)	= (x1(k)+eta)/(1.d0-x1(k))
			xx		= (1.d0-x1(k))*(1.d0-x1(k))
			wght(k)	= et1*p1(k)/xx
		end do
		deallocate(x1,p1)
		return
	end subroutine generate_imag_epoints

	! Generating energy points in the real axis for third integration
	subroutine generate_real_epoints(e)
		implicit none
		integer												:: i,j,k,AllocateStatus
		real(double)									:: e1,e2
		real(double),intent(in)				:: e

		allocate( x2(pn2),p2(pn2), STAT = AllocateStatus )

		if (abs(e).ge.1.d-10) then
			nepoints = pnt
			do k=1,parts3
				e1 = Ef-((parts3-(k-1))*e/parts3)
				e2 = Ef-((parts3-k)*e/parts3)

				i = (k-1)*n3gl+1
				j = k*n3gl

				call gauleg(e1,e2,x2(i:j),p2(i:j),n3gl)
			end do
		else
			nepoints = pn1
		end if
		return
	end subroutine generate_real_epoints

	!	This subroutine returns the N points in the abcissae x and its weight w using Gauss-Legendre
	subroutine gauleg(x1,x2,x,w,n)
	  implicit none
	  integer, 			intent(in) :: n
	  real(double), intent(in) :: x1, x2
	  real(double), dimension(n), intent(out) :: x, w
	  integer :: i, j, m
	  real(double) :: p1, p2, p3, pp, xl, xm, z, z1
	  real(double), parameter :: eps=3.d-14

	  m = (n+1)/2
	  xm = 0.5d0*(x2+x1)
	  xl = 0.5d0*(x2-x1)
	  do i=1,m
	    z = cos(3.141592654d0*(i-0.25d0)/(n+0.5d0))
	    z1 = 0.0
	    do while(abs(z-z1).gt. eps)
	      p1 = 1.0d0
	      p2 = 0.0d0
	      do j=1,n
	        p3 = p2
	        p2 = p1
	        p1 = ((2.0d0*j-1.0d0)*z*p2-(j-1.0d0)*p3)/j
	      end do
	      pp = n*(z*p1-p2)/(z*z-1.0d0)
	      z1 = z
	      z = z1 - p1/pp
	    end do
	    x(i) = xm - xl*z
	    x(n+1-i) = xm + xl*z
	    w(i) = (2.0d0*xl)/((1.0d0-z*z)*pp*pp)
	    w(n+1-i) = w(i)
	  end do
	  return
	end subroutine gauleg

end module mod_generate_epoints