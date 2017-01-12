! This subroutine calculates the inverse of NNxNN matrix 'matriz'
subroutine invers(matriz,nn)
  use mod_f90_kind
  use mod_mpi_pars
  use mod_parameters, only: outputunit
  implicit none
  integer :: nn,info
  integer :: lwork
  integer, dimension(nn) :: ipiv
  complex(double), dimension(nn,nn) :: matriz
  complex(double), dimension(nn*4) :: work

  lwork = 4*nn
  info = 0
  call zgetrf(nn,nn,matriz,nn,ipiv,info)
  if (info.gt.0) then
    write(outputunit,"('[invers] Singular matrix! info = ',i0,' on rank ',i0)") info,myrank
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  if (info.lt.0) then
    write(outputunit,"('[invers] Illegal value of argument ',i0,' on rank ',i0)") -info,myrank
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  call zgetri(nn,matriz,nn,ipiv,work,lwork,info)

  if (info.ne.0) then
    write(outputunit,"('[invers] Singular matrix! info = ',i0,' on rank ',i0)") info,myrank
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

  return
end subroutine invers