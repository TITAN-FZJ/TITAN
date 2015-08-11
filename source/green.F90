! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine green(er,ei,kp,gf)
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_magnet
	implicit none
	integer			:: i,j,i0,i1,j0,j1
	real(double), intent(in)	:: er,ei,kp(3)
	complex(double)	:: ec
	complex(double),dimension((Npl+2)*18,(Npl+2)*18)	:: gslab,identes,hk
	complex(double),dimension(Npl,Npl,18,18)	:: gf

	ec		= cmplx(er,ei,double)
	identes     = zero
	do i=1,(Npl+2)*18
	 identes(i,i) = zum
	end do

	gslab	= zero

	call hamiltk(kp,hk)

	gslab = ec*identes - hk

	call invers(gslab,(Npl+2)*18)

	! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
	do i=1,Npl
		do j=1,Npl
			i0 = i*18+1
			i1 = i0+17
			j0 = j*18+1
			j1 = j0+17
			gf(i,j,:,:) = gslab(i0:i1,j0:j1)
		end do
	end do

	return
end subroutine green

! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine green_es(er,ei,kp,gf)
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_magnet
	implicit none
	integer			:: i,j,i0,i1,j0,j1
	real(double), intent(in)	:: er,ei,kp(3)
	complex(double)	:: ec
	complex(double),dimension((Npl+2)*18,(Npl+2)*18)	:: gslab,identes,hk
	complex(double),dimension(Npl+2,Npl+2,18,18)	:: gf

	ec		= cmplx(er,ei,double)
	identes     = zero
	do i=1,(Npl+2)*18
	 identes(i,i) = zum
	end do

	gslab	= zero

	call hamiltk(kp,hk)

	gslab = ec*identes - hk

	call invers(gslab,(Npl+2)*18)

	! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
	do i=1,Npl+2
		do j=1,Npl+2
			i0 = (i-1)*18+1
			i1 = i0+17
			j0 = (j-1)*18+1
			j1 = j0+17
			gf(i,j,:,:) = gslab(i0:i1,j0:j1)
		end do
	end do

	return
end subroutine green_es
