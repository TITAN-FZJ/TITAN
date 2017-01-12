module mod_mpi_pars
  use mod_f90_kind
  use MPI
  implicit none
  integer :: myrank,numprocs,errorcode,ierr,mcount,mpitag
  integer, dimension(MPI_STATUS_SIZE) :: stat
  ! Grid for energy loop calculations
  integer :: myrank_row,myrank_col
  integer :: MPIsteps,MPIpts
  integer :: MPIComm_Grid,MPI_Comm_Row,MPI_Comm_Col,MPI_Comm_Energy
  real(double) :: MPIdelta
  ! Grid for field loop calculation
  integer :: myrank_row_hw,myrank_col_hw
  integer :: MPIsteps_hw,MPIpts_hw=1
  integer :: MPIComm_Grid_hw,MPI_Comm_Row_hw,MPI_Comm_Col_hw


contains

  ! Bidimensional array should look like:
  !    myrank_col: / myrank_row: 0 1 2 3 ... pnt-1 (number inside a row)
  !        0
  !        1
  !        2
  !       ...
  !     MPIpts-1
  !(number inside a col)
  subroutine build_cartesian_grid()
    use mod_parameters, only: total_hw_npt1, Npl_i, Npl_f, pnt, npt1, npts, deltae, emin, emax, outputunit
    implicit none
    logical :: lreorder,lperiodic(2),lrow(2),lcol(2)
    integer :: MPIdims(2)

    MPIpts  = ceiling(dble(numprocs)/dble(pnt)) ! Number of rows to be used
    MPIdims = [MPIpts,pnt]
    if(numprocs.le.pnt) then  ! If number of processes is less than necessary for 1 energy integral
      MPIdims  = [MPIpts,numprocs]  ! Create only one array of processes, i.e., MPIpts = 1
      MPIsteps = npt1
    else
      if(mod(numprocs,pnt).ne.0) then
        if(myrank.eq.0) then
          write(outputunit,"('[build_cartesian_grid] ************************************** ERROR: **************************************')")
          write(outputunit,"('[build_cartesian_grid]      Number of processes not commensurable with total energy integral points!')")
          write(outputunit,"('[build_cartesian_grid]    Number of MPI processes: ',i0)") numprocs
          write(outputunit,"('[build_cartesian_grid]    Number of points required: ',i0)") npt1
          write(outputunit,"('[build_cartesian_grid]    Number of points in the energy integral: ',i0)") pnt
          write(outputunit,"('[build_cartesian_grid] ************************************************************************************')")
        end if
        call MPI_Finalize(ierr)
        stop
      end if
      MPIsteps = 1
      if(total_hw_npt1.eq.1) then ! If there's no loop on field, complete number of points
        if(npt1*pnt.lt.numprocs) then ! If the number of processors is larger (but commensurable) than complete calculation
          if(myrank.eq.0) then
            write(outputunit,"('[build_cartesian_grid] ************************************* WARNING: *************************************')")
            write(outputunit,"('[build_cartesian_grid]              Number of processes exceeds the total needed! ')")
            write(outputunit,"('[build_cartesian_grid]     Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIpts-npt1
            write(outputunit,"('[build_cartesian_grid] ************************************************************************************')")
          end if
          npt1 = MPIpts
        else if(npt1*pnt.gt.numprocs) then ! If the number of processors is smaller than complete calculation, check commensurability
          MPIsteps = ceiling(dble(npt1)/dble(MPIpts))
          if(mod(npt1,MPIpts).ne.0) then
            if(myrank.eq.0) then
              write(outputunit,"('[build_cartesian_grid] ************************************* WARNING: *************************************')")
              write(outputunit,"('[build_cartesian_grid]    Number of points to be calculated is not commensurable with processes used!')")
              write(outputunit,"('[build_cartesian_grid]   Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIsteps*MPIpts-npt1
              write(outputunit,"('[build_cartesian_grid] ************************************************************************************')")
            end if
            npt1 = MPIsteps*MPIpts
          end if
        end if
      else
        MPIpts  = npt1
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
    call MPI_Cart_sub(MPIComm_Grid,lrow,MPI_Comm_Col,ierr) ! communicator inside a column (between rows)
    lcol = [.false.,.true.]
    call MPI_Cart_sub(MPIComm_Grid,lcol,MPI_Comm_Row,ierr) ! communicator inside a row (between columns)

    ! Obtaining process rank inside its row and column
    call MPI_Comm_rank(MPI_Comm_Row,myrank_row,ierr) ! Obtaining rank number inside its row
    call MPI_Comm_rank(MPI_Comm_Col,myrank_col,ierr) ! Obtaining rank number inside its column

    if(myrank.eq.0) write(outputunit,"('[build_cartesian_grid] Created grid with ',i0,' rows (myrank_col) x ',i0,' columns (myrank_row)')") MPIdims(1),MPIdims(2)

    return
  end subroutine build_cartesian_grid


  subroutine build_cartesian_grid_field(elements_in_a_row)
    use mod_parameters, only: total_hw_npt1,outputunit,npt1
    implicit none
    logical :: lreorder,lperiodic(2),lrow(2),lcol(2)
    integer :: MPIdims(2)
    integer, intent(in) :: elements_in_a_row

    MPIpts_hw  = ceiling(dble(numprocs)/dble(elements_in_a_row)) ! Number of rows to be used
    MPIdims    = [MPIpts_hw,elements_in_a_row]
    if(numprocs.le.elements_in_a_row) then  ! If number of processes is less than necessary for 1 energy integral
      MPIdims     = [MPIpts_hw,numprocs]  ! Create only one array of processes, i.e., MPIpts_hw = 1
      MPIsteps_hw = total_hw_npt1
    else
      if(mod(numprocs,elements_in_a_row).ne.0) then
        if(myrank.eq.0) then
          write(outputunit,"('[build_cartesian_grid_field] ************************************** ERROR: **************************************')")
          write(outputunit,"('[build_cartesian_grid_field]      Number of processes not commensurable with total energy integral points!')")
          write(outputunit,"('[build_cartesian_grid_field]    Number of MPI processes: ',i0)") numprocs
          write(outputunit,"('[build_cartesian_grid_field]    Number of points required: ',i0)") total_hw_npt1
          write(outputunit,"('[build_cartesian_grid_field]    Number of elements in a row: ',i0)") elements_in_a_row
          write(outputunit,"('[build_cartesian_grid_field] ************************************************************************************')")
        end if
        call MPI_Finalize(ierr)
        stop
      end if
      MPIsteps_hw = 1
      if(npt1*total_hw_npt1*elements_in_a_row.lt.numprocs) then ! If the number of processors is larger (but commensurable) than complete calculation
        if(myrank.eq.0) then
          write(outputunit,"('[build_cartesian_grid_field] ************************************** ERROR: **************************************')")
          write(outputunit,"('[build_cartesian_grid_field]              Number of processes exceeds the total needed! ')")
          write(outputunit,"('[build_cartesian_grid_field]     Complete the ',i0,' points with more ',i0,' points not to waste computing time. ')") total_hw_npt1,MPIpts_hw-total_hw_npt1
          write(outputunit,"('[build_cartesian_grid_field] ************************************************************************************')")
        end if
        call MPI_Finalize(ierr)
        stop
      else if(total_hw_npt1*elements_in_a_row.gt.numprocs) then ! If the number of processors is smaller than complete calculation, check commensurability
        MPIsteps_hw = ceiling(dble(total_hw_npt1)/dble(MPIpts_hw))
        if(mod(npt1*total_hw_npt1,MPIpts_hw).ne.0) then
          if(myrank.eq.0) then
            write(outputunit,"('[build_cartesian_grid_field] ************************************** ERROR: **************************************')")
            write(outputunit,"('[build_cartesian_grid_field]    Number of points to be calculated is not commensurable with processes used!')")
            write(outputunit,"('[build_cartesian_grid_field]    Complete ',i0,' points with more ',i0,' points not to waste computing time. ')") total_hw_npt1,MPIsteps_hw*MPIpts_hw-total_hw_npt1
            write(outputunit,"('[build_cartesian_grid_field] ************************************************************************************')")
          end if
          call MPI_Finalize(ierr)
          stop
        end if
      end if
    end if

    ! Creating bidimensiontal Grid of tasks
    lperiodic = [.false.,.false.]
    lreorder  = .true.
    call MPI_Cart_create(MPI_COMM_WORLD,2,MPIdims,lperiodic,lreorder,MPIComm_Grid_hw,ierr)

    ! Creating subarrays of rows and columns
    lrow = [.true.,.false.]
    call MPI_Cart_sub(MPIComm_Grid_hw,lrow,MPI_Comm_Col_hw,ierr) ! communicator inside a column (between rows)
    lcol = [.false.,.true.]
    call MPI_Cart_sub(MPIComm_Grid_hw,lcol,MPI_Comm_Row_hw,ierr) ! communicator inside a row (between columns)

    ! Obtaining process rank inside its row and column
    call MPI_Comm_rank(MPI_Comm_Row_hw,myrank_row_hw,ierr) ! Obtaining rank number inside its row
    call MPI_Comm_rank(MPI_Comm_Col_hw,myrank_col_hw,ierr) ! Obtaining rank number inside its column

    if(myrank.eq.0) write(outputunit,"('[build_cartesian_grid_field] Created grid with ',i0,' rows (myrank_col_hw) x ',i0,' columns (myrank_row_hw)')") MPIdims(1),MPIdims(2)

    return
  end subroutine build_cartesian_grid_field
end module mod_mpi_pars