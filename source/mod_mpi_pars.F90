module mod_mpi_pars
	use mod_f90_kind
	use MPI
	implicit none
	integer :: myrank,numprocs,errorcode,ierr,count,mcount
	integer :: myrank_row,myrank_col
	integer :: MPIdims(2),MPIsteps,MPIpts
	real(double) :: MPIdelta
	integer :: MPIComm_Grid,MPIComm_Row,MPIComm_Col
	logical :: lreorder,lperiodic(2),lrow(2),lcol(2)
	integer, dimension(MPI_STATUS_SIZE) :: stat
end module mod_mpi_pars