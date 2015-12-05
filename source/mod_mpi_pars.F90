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

contains

  ! Bidimensional array should look like:
  !      0 1 2 3 ... pnt-1
  !    0
  !    1
  !    2
  !   ...
  ! MPIpts-1
	subroutine build_cartesian_grid()
		use mod_parameters, only: pnt, npt1, npts, deltae, emin, emax
		implicit none

    MPIpts = ceiling(dble(numprocs)/dble(pnt)) ! Number of rows to be used
    MPIdims = [MPIpts,pnt]
    if(numprocs.le.pnt) then  ! If number of processes is less than necessary for 1 energy integral
      MPIdims  = [MPIpts,numprocs]  ! Create only one array of processes, i.e., MPIpts = 1
      MPIsteps = npt1
    else
      if(mod(numprocs,pnt).ne.0) then
        if(myrank.eq.0) then
          write(*,"('[main] ERROR: number of processes not commensurable with total energy integral points!')")
          write(*,"('[main] Number of MPI processes: ',i0)") numprocs
          write(*,"('[main] Number of points required: ',i0)") npt1
          write(*,"('[main] Number of points in the energy integral: ',i0)") pnt
        end if
        stop
      end if
      MPIsteps = 1
      if(npt1*pnt.lt.numprocs) then ! If the number of processors is larger than complete calculation
        if(myrank.eq.0) then
          write(*,"('[main] ******************************** WARNING: ********************************')")
          write(*,"('[main]              Number of processes exceeds the total needed! ')")
          write(*,"('[main]     Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIpts-npt1
          write(*,"('[main] **************************************************************************')")
        end if
        npt1 = MPIpts
      else if(npt1*pnt.gt.numprocs) then ! If the number of processors is smaller than complete calculation, check commensurability
        MPIsteps = ceiling(dble(npt1)/dble(MPIpts))
        if(mod(npt1,MPIpts).ne.0) then
          if(myrank.eq.0) then
            write(*,"('[main] ******************************** WARNING: ********************************')")
            write(*,"('[main] Number of points to be calculated is not commensurable with processes used!')")
            write(*,"('[main]     Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIsteps*MPIpts-npt1
            write(*,"('[main] **************************************************************************')")
          end if
          npt1 = MPIsteps*MPIpts
        end if
      end if
    end if
    if(npt1.ne.1) then
      npts = npt1-1
    else
      npts = npt1
    end if
    ! Calculating variations of energy
    deltae = (emax - emin)/npts      ! variation of energy between each point
    MPIdelta = deltae*MPIpts         ! variation of energy between MPI steps

    ! Creating bidimensiontal Grid of tasks
    lperiodic = [.false.,.false.]
    lreorder  = .true.
    call MPI_Cart_create(MPI_COMM_WORLD,2,MPIdims,lperiodic,lreorder,MPIComm_Grid,ierr)

    ! Creating subarrays of rows and columns
    lrow = [.true.,.false.]
    call MPI_Cart_sub(MPIComm_Grid,lrow,MPIComm_Col,ierr) ! communicator inside a column (between rows)
    lcol = [.false.,.true.]
    call MPI_Cart_sub(MPIComm_Grid,lcol,MPIComm_Row,ierr) ! communicator inside a row (between columns)

    ! Obtaining process rank inside its row and column
    call MPI_Comm_rank(MPIComm_Row,myrank_row,ierr) ! Obtaining rank number inside its row
    call MPI_Comm_rank(MPIComm_Col,myrank_col,ierr) ! Obtaining rank number inside its column
	end subroutine build_cartesian_grid
end module mod_mpi_pars