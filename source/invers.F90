!	This subroutine calculates the inverse of NNxNN matrix 'matriz'

subroutine invers(matriz,nn)
	use mod_f90_kind
	implicit none
	integer :: nn,info
	integer :: lwork
	integer, dimension(nn) :: ipiv
	complex(double), dimension(nn,nn) :: matriz
	complex(double), dimension(nn*4) :: work

	lwork = 4*nn
	info = 0
	call zgetrf(nn,nn,matriz,nn,ipiv,info)
	if (info/=0) then
		write(*,*)'ifail=',info
	end if
	call zgetri(nn,matriz,nn,ipiv,work,lwork,info)

	if (info/=0) then
		write(*,*)'ifail=',info
		stop '[invers] Singular matrix!'
	end if

	return
end subroutine invers