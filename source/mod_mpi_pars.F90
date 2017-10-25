module mod_mpi_pars
  use mod_f90_kind, only: double
  use MPI
  implicit none
  integer :: ierr, errorcode
  integer :: myrank,numprocs !,mcount,mpitag
  ! integer :: numprocs_row, numprocs_col
  ! integer :: numprocs_row_hw, numprocs_col_hw
  integer, dimension(MPI_STATUS_SIZE) :: stat
  ! Grid for energy loop calculations
  ! integer :: myrank_row,myrank_col
  ! integer :: MPIsteps,MPIpts
  ! integer :: MPIComm_Grid,MPI_Comm_Row,MPI_Comm_Col,MPI_Comm_Energy
  ! real(double) :: MPIdelta
  ! Grid for field loop calculation
  ! integer :: myrank_row_hw,myrank_col_hw
  ! integer :: MPIsteps_hw,MPIpts_hw=1
  ! integer :: MPIComm_Grid_hw,MPI_Comm_Row_hw,MPI_Comm_Col_hw


  ! New Stuff
  integer :: FieldComm
  !! Magnetic Field Communicator
  integer :: sField
  !! Size of Magnetic Field Communicator
  integer :: rField
  !! Rank in Magnetic Field Communicator

  integer, dimension(2) :: FreqComm
  !! Frequency Communicators, (1) Row, (2) Col (for sending to 1 Process)
  integer, dimension(2) :: sFreq
  !! Size of Frequency Communicators
  integer, dimension(2) :: rFreq
  !! Rank in Frequency Communicators

  integer :: startField
  !! First Point of Magnetic Field workload
  integer :: endField
  !! Last Point of Magnetic Field workload
  integer :: startFreq
  !! First Point of Frequency workload
  integer :: endFreq
  !! Last Point of Frequency workload

contains

  subroutine Initialize_MPI()
    implicit none

 !#ifndef _UFF
 !#endif
    integer :: provided
    !integer :: ierr

#ifdef _OPENMP
#ifdef _UFF
      call MPI_Init(ierr)
#else
      call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
#endif
      call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
#else
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
#endif
    return
  end subroutine

  ! Wrapper function for MPI Abort
  subroutine abortProgram(str)
    use mod_parameters, only: outputunit
    implicit none
    character(len=*), intent(in) :: str
    !! Abort Message

    write(outputunit,"(a)") str
    close(outputunit)
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)

    return
  end subroutine abortProgram

  !
  ! subroutine setup_MPI_grid(itype, pn1, npt1, pnt, total_hw_npt1, npts, deltae, emin, emax)
  !   !! Chooses how to setup the MPI Grids depending on the calculation to be performed and the amount of available nodes
  !   use mod_f90_kind, only: double
  !   implicit none
  !   integer, intent(in) :: itype
  !   !! Type of calculation
  !   integer, intent(in) :: pn1
  !   !! Energy integration points along the imaginary axis (eta -> infinity)
  !   integer, intent(in) :: pnt
  !   !! Total number of energy integration points along imaginary and real axis (eta -> infinity + Ef-e -> Ef)
  !   integer, intent(in) :: total_hw_npt1
  !   !! Number of magnetic field points
  !   real(double), intent(in) :: emin, emax
  !   !! Range of frequency points
  !   real(double), intent(inout) :: deltae
  !   !! Delta inbetween frequency points
  !   integer, intent(inout) :: npts
  !   !! Number of frequency points
  !   integer, intent(inout) :: npt1
  !   !! Number of frequency points + 1
  !   integer :: ierr
  !
  !   if((itype==1).or.(itype==6) .or. (itype==10)) then ! Create column for field loop (only energy integration along imaginary axis)
  !     call build_cartesian_grid_field(pn1, total_hw_npt1, npt1)
  !   end if
  !   if( ((itype>=3).and.(itype<=5)) .or. itype==0) then ! Create column for field loop (no energy integration)
  !     call build_cartesian_grid_field(1, total_hw_npt1, npt1)
  !   end if
  !   if((itype>=7).and.(itype<=8)) then ! Create matrix for energy dependence and integration
  !     call build_cartesian_grid(total_hw_npt1, pnt, npt1, npts, deltae, emin, emax)
  !     call build_cartesian_grid_field(npt1*pnt, total_hw_npt1, npt1)
  !   end if
  !   if(itype==9) then ! Create matrix for dclimit
  !     call build_cartesian_grid_field(pnt, total_hw_npt1, npt1)
  !     call MPI_COMM_DUP(MPI_Comm_Col_hw,MPI_Comm_Col,ierr)
  !     call MPI_COMM_DUP(MPI_Comm_Row_hw,MPI_Comm_Row,ierr)
  !     myrank_col = myrank_col_hw
  !     myrank_row = myrank_row_hw
  !   end if
  !
  ! end subroutine setup_MPI_grid
  !
  ! ! Bidimensional array should look like:
  ! !    myrank_col: / myrank_row: 0 1 2 3 ... pnt-1 (number inside a row)
  ! !        0
  ! !        1
  ! !        2
  ! !       ...
  ! !     MPIpts-1
  ! !(number inside a col)
  ! subroutine build_cartesian_grid(total_hw_npt1, pnt, npt1, npts, deltae, emin, emax)
  !   use mod_f90_kind, only: double
  !   use mod_parameters, only: outputunit
  !   implicit none
  !   integer, intent(in) :: total_hw_npt1
  !   !! Number of magnetic field points
  !   integer, intent(in) :: pnt
  !   !! Total number of points for Energy Integration (i.e. real + imaginary axis)
  !   real(double), intent(in) :: emin, emax
  !   !! Frequency boundaries
  !   integer, intent(inout) :: npts
  !   !! Total number of frequency points
  !   integer, intent(inout) :: npt1
  !   !! Total number of frequency points + 1
  !   real(double), intent(inout) :: deltae
  !   !! Delta between two frequency points
  !
  !   logical :: lreorder,lperiodic(2),lrow(2),lcol(2)
  !   integer :: MPIdims(2)
  !
  !
  !   MPIpts  = ceiling(dble(numprocs)/dble(pnt)) ! Number of rows to be used
  !   MPIdims = [MPIpts,pnt]
  !
  !   if(numprocs <= pnt) then  ! If number of processes is less than necessary for 1 energy integral
  !     MPIdims  = [MPIpts,numprocs]  ! Create only one array of processes, i.e., MPIpts = 1 (one row with numprocs columns)
  !     MPIsteps = npt1 ! perform npt1 iterations
  !   else ! If number of processes exceeds the number necessary for a single energy integral
  !     if(mod(numprocs,pnt)/=0) then ! If number of processes is not divisible by the number of points per integral -> abort
  !       if(myrank==0) then
  !         write(outputunit,"('[build_cartesian_grid] ************************************** ERROR: **************************************')")
  !         write(outputunit,"('[build_cartesian_grid]      Number of processes not commensurable with total energy integral points!')")
  !         write(outputunit,"('[build_cartesian_grid]    Number of MPI processes: ',i0)") numprocs
  !         write(outputunit,"('[build_cartesian_grid]    Number of points required: ',i0)") npt1
  !         write(outputunit,"('[build_cartesian_grid]    Number of points in the energy integral: ',i0)") pnt
  !         write(outputunit,"('[build_cartesian_grid] ************************************************************************************')")
  !       end if
  !       call MPI_Finalize(ierr)
  !       stop
  !     end if
  !
  !     MPIsteps = 1
  !     if(total_hw_npt1==1) then ! If there's no loop on field, complete number of points
  !       if(numprocs > npt1*pnt) then ! If the number of processors is larger (but commensurable) than complete calculation (add points to calculation)
  !         if(myrank==0) then
  !           write(outputunit,"('[build_cartesian_grid] ************************************* WARNING: *************************************')")
  !           write(outputunit,"('[build_cartesian_grid]              Number of processes exceeds the total needed! ')")
  !           write(outputunit,"('[build_cartesian_grid]     Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIpts-npt1
  !           write(outputunit,"('[build_cartesian_grid] ************************************************************************************')")
  !         end if
  !         npt1 = MPIpts
  !       else if(numprocs < npt1*pnt) then ! If the number of processors is smaller than complete calculation, check commensurability
  !         MPIsteps = ceiling(dble(npt1)/dble(MPIpts))
  !         if(mod(npt1,MPIpts)/=0) then ! If number of frequency points is not divisible by number of rows add points to calculation
  !           if(myrank==0) then
  !             write(outputunit,"('[build_cartesian_grid] ************************************* WARNING: *************************************')")
  !             write(outputunit,"('[build_cartesian_grid]    Number of points to be calculated is not commensurable with processes used!')")
  !             write(outputunit,"('[build_cartesian_grid]   Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIsteps*MPIpts-npt1
  !             write(outputunit,"('[build_cartesian_grid] ************************************************************************************')")
  !           end if
  !           npt1 = MPIsteps*MPIpts
  !         end if
  !       end if
  !     else ! If there is a loop on the field,
  !       if(numprocs > npt1*pnt) then ! If number of processes exceeds number of points necessary for energy integration and frequencies only use first npt1*pnt nodes
  !         MPIpts = npt1
  !       else
  !         MPIsteps = ceiling(dble(npt1)/dble(MPIpts))
  !         if(mod(npt1,MPIpts)/=0) then ! If number of frequencies is not devisible by number of mpi points, add points to calculation
  !           if(myrank==0) then
  !             write(outputunit,"('[build_cartesian_grid] ************************************* WARNING: *************************************')")
  !             write(outputunit,"('[build_cartesian_grid]    Number of points to be calculated is not commensurable with processes used!')")
  !             write(outputunit,"('[build_cartesian_grid]   Completing ',i0,' points with more ',i0,' points not to waste computing time. ')") npt1,MPIsteps*MPIpts-npt1
  !             write(outputunit,"('[build_cartesian_grid] ************************************************************************************')")
  !           end if
  !           npt1 = MPIsteps*MPIpts
  !         end if
  !       end if
  !     end if
  !   end if
  !   if(npt1/=1) then
  !     npts = npt1-1
  !   else
  !     npts = npt1
  !   end if
  !   ! Calculating variations of energy
  !   deltae = (emax - emin)/npts      ! variation of energy between each point
  !   MPIdelta = deltae*MPIpts         ! variation of energy between MPI steps
  !
  !   ! Creating bidimensiontal Grid of tasks
  !   lperiodic = [.false.,.false.]
  !   lreorder  = .true.
  !   call MPI_Cart_create(MPI_COMM_WORLD,2,MPIdims,lperiodic,lreorder,MPIComm_Grid,ierr)
  !
  !   ! Creating subarrays of rows and columns
  !   lrow = [.true.,.false.]
  !   call MPI_Cart_sub(MPIComm_Grid,lrow,MPI_Comm_Col,ierr) ! communicator inside a column (between rows)
  !   lcol = [.false.,.true.]
  !   call MPI_Cart_sub(MPIComm_Grid,lcol,MPI_Comm_Row,ierr) ! communicator inside a row (between columns)
  !
  !   ! Obtaining process rank inside its row and column
  !   call MPI_Comm_rank(MPI_Comm_Row, myrank_row, ierr) ! Obtaining rank number inside its row
  !   call MPI_Comm_size(MPI_Comm_Row, numprocs_row, ierr) ! Obtain size of row
  !   call MPI_Comm_rank(MPI_Comm_Col, myrank_col, ierr) ! Obtaining rank number inside its column
  !   call MPI_Comm_size(MPI_Comm_Col, numprocs_col, ierr) ! Obtain size of row
  !
  !   if(myrank==0) write(outputunit,"('[build_cartesian_grid] Created grid with ',i0,' rows (myrank_col) x ',i0,' columns (myrank_row)')") MPIdims(1),MPIdims(2)
  !
  !   return
  ! end subroutine build_cartesian_grid
  !
  !
  ! subroutine build_cartesian_grid_field(elements_in_a_row, total_hw_npt1, npt1)
  ! !! Create communicator for magnetic field loop
  !   use mod_parameters, only: outputunit
  !   implicit none
  !   integer, intent(in) :: elements_in_a_row
  !   !! Number of elements that are to be in a row
  !   integer, intent(in) :: total_hw_npt1
  !   !! Number of magnetic field points (number of rows)
  !   integer, intent(in) :: npt1
  !   !! Number of frequency points + 1
  !   logical :: lreorder,lperiodic(2),lrow(2),lcol(2)
  !   integer :: MPIdims(2)
  !
  !   MPIpts_hw  = ceiling(dble(numprocs)/dble(elements_in_a_row)) ! Number of rows to be used
  !   MPIdims    = [MPIpts_hw,elements_in_a_row]
  !   if(numprocs<=elements_in_a_row) then  ! If number of processes is less than necessary for 1 energy integral
  !     MPIdims     = [MPIpts_hw,numprocs]  ! Create only one array of processes, i.e., MPIpts_hw = 1
  !     MPIsteps_hw = total_hw_npt1
  !   else
  !     if(mod(numprocs,elements_in_a_row)/=0) then
  !       if(myrank==0) then
  !         write(outputunit,"('[build_cartesian_grid_field] ************************************** ERROR: **************************************')")
  !         write(outputunit,"('[build_cartesian_grid_field]      Number of processes not commensurable with total energy integral points!')")
  !         write(outputunit,"('[build_cartesian_grid_field]    Number of MPI processes: ',i0)") numprocs
  !         write(outputunit,"('[build_cartesian_grid_field]    Number of points required: ',i0)") total_hw_npt1
  !         write(outputunit,"('[build_cartesian_grid_field]    Number of elements in a row: ',i0)") elements_in_a_row
  !         write(outputunit,"('[build_cartesian_grid_field] ************************************************************************************')")
  !       end if
  !       call MPI_Finalize(ierr)
  !       stop
  !     end if
  !     MPIsteps_hw = 1
  !     if(npt1*total_hw_npt1*elements_in_a_row<numprocs) then ! If the number of processors is larger (but commensurable) than complete calculation
  !       if(myrank==0) then
  !         write(outputunit,"('[build_cartesian_grid_field] ************************************** ERROR: **************************************')")
  !         write(outputunit,"('[build_cartesian_grid_field]              Number of processes exceeds the total needed! ')")
  !         write(outputunit,"('[build_cartesian_grid_field]     Complete the ',i0,' points with more ',i0,' points not to waste computing time. ')") total_hw_npt1,MPIpts_hw-total_hw_npt1
  !         write(outputunit,"('[build_cartesian_grid_field] ************************************************************************************')")
  !       end if
  !       call MPI_Finalize(ierr)
  !       stop
  !     else if(total_hw_npt1*elements_in_a_row>numprocs) then ! If the number of processors is smaller than complete calculation, check commensurability
  !       MPIsteps_hw = ceiling(dble(total_hw_npt1)/dble(MPIpts_hw))
  !       if(mod(npt1*total_hw_npt1,MPIpts_hw)/=0) then
  !         if(myrank==0) then
  !           write(outputunit,"('[build_cartesian_grid_field] ************************************** ERROR: **************************************')")
  !           write(outputunit,"('[build_cartesian_grid_field]    Number of points to be calculated is not commensurable with processes used!')")
  !           write(outputunit,"('[build_cartesian_grid_field]    Complete ',i0,' points with more ',i0,' points not to waste computing time. ')") total_hw_npt1,MPIsteps_hw*MPIpts_hw-total_hw_npt1
  !           write(outputunit,"('[build_cartesian_grid_field] ************************************************************************************')")
  !         end if
  !         call MPI_Finalize(ierr)
  !         stop
  !       end if
  !     end if
  !   end if
  !
  !   ! Creating bidimensiontal Grid of tasks
  !   lperiodic = [.false.,.false.]
  !   lreorder  = .true.
  !
  !   call MPI_Cart_create(MPI_COMM_WORLD,2,MPIdims,lperiodic,lreorder,MPIComm_Grid_hw,ierr)
  !
  !   ! Creating subarrays of rows and columns
  !   lrow = [.true.,.false.]
  !   call MPI_Cart_sub(MPIComm_Grid_hw,lrow,MPI_Comm_Col_hw,ierr) ! communicator inside a column (between rows)
  !   lcol = [.false.,.true.]
  !   call MPI_Cart_sub(MPIComm_Grid_hw,lcol,MPI_Comm_Row_hw,ierr) ! communicator inside a row (between columns)
  !
  !   ! Obtaining process rank inside its row and column
  !   call MPI_Comm_rank(MPI_Comm_Row_hw,myrank_row_hw,ierr) ! Obtaining rank number inside its row
  !   call MPI_Comm_size(MPI_Comm_Row_hw, numprocs_row_hw,ierr)      ! Obtain size of row
  !   call MPI_Comm_rank(MPI_Comm_Col_hw,myrank_col_hw,ierr) ! Obtaining rank number inside its column
  !   call MPI_Comm_size(MPI_Comm_Col_hw, numprocs_col_hw,ierr)      ! Obtain size of row
  !
  !   if(myrank==0) write(outputunit,"('[build_cartesian_grid_field] Created grid with ',i0,' rows (myrank_col_hw) x ',i0,' columns (myrank_row_hw)')") MPIdims(1),MPIdims(2)
  !
  !   return
  ! end subroutine build_cartesian_grid_field

  integer function genFieldComm(nFields, FieldID, comm)
     implicit none
     integer, intent(in) :: nFields
     !! Number of field points to be calculated in parallel
     integer, intent(out) :: FieldID
     !! Communicator ID
     integer, intent(in) :: comm
     !! MPI Communicator to split up
     integer :: rank, procs, ierr
     integer :: group
     call MPI_Comm_rank(comm, rank, ierr)
     call MPI_Comm_size(comm, procs, ierr)
     if(mod(procs,nFields) /= 0) call abortProgram("[genFieldComm] Number of MPI processes have to be divisible by nFields")
     group = procs / nFields
     FieldID = rank / group
     call MPI_Comm_split(comm, FieldID, rank, genFieldComm, ierr)
     return
  end function genFieldComm

  function genFreqComm(nFreq, FreqID, comm)
     implicit none
     integer, intent(in) :: nFreq
     !! Number of frequency points to be calculated in parallel
     integer, intent(out) :: FreqID
     !! Communicator ID
     integer, intent(in) :: comm
     !! MPI Communicator to split up
     integer, dimension(2) :: genFreqComm

     integer :: rank, rank_row, procs, ierr
     integer :: group

     call MPI_Comm_rank(comm, rank, ierr)
     call MPI_Comm_size(comm, procs,ierr)
     if(mod(procs,nFreq) /= 0) call abortProgram("[genFreqComm] Number of MPI processes have to be divisble by nFreq")
     group = procs / nFreq
     FreqID = rank / group
     call MPI_Comm_split(comm, FreqID, rank, genFreqComm(1), ierr)

     call MPI_Comm_rank(genFreqComm(1), rank_row, ierr)
     call MPI_Comm_split(comm, rank_row, rank, genFreqComm(2), ierr)

     return
  end function genFreqComm

  subroutine genMPIGrid(nFields, nFieldPoints, nFreq, nFreqPoints)
     implicit none
     integer, intent(in) :: nFields
     integer, intent(in) :: nFreq
     integer, intent(in) :: nFieldPoints
     integer, intent(in) :: nFreqPoints
     integer :: FreqID, FieldID

     integer :: i, ierr

     FieldComm = genFieldComm(nFields, FieldID, MPI_COMM_WORLD)
     FreqComm = genFreqComm(nFreq, FreqID, FieldComm)
     call MPI_Comm_rank(FieldComm, rField, ierr)
     call MPI_Comm_size(FieldComm, sField, ierr)

     do i = 1,2
        call MPI_Comm_rank(FreqComm(i), rFreq(i), ierr)
        call MPI_Comm_size(FreqComm(i), sFreq(i), ierr)
     end do

     ! Calculate Field workload for each Field Communicator
     call calcWorkload(nFieldPoints,nFields,FieldID,startField,endField)
     ! Calculate Field workload for each Frequency Communicator
     call calcWorkload(nFreqPoints,nFreq,FreqID,startFreq,endFreq)

     return
  end subroutine genMPIGrid

  subroutine calcWorkload(points, procs, rank, firstPoint, lastPoint)
     implicit none
     integer, intent(in) :: points
     !! Number of points to split up
     integer, intent(in) :: procs
     !! Number of processes
     integer, intent(in) :: rank
     !! MPI Rank
     integer, intent(out) :: firstPoint
     !! First point of workload
     integer, intent(out) :: lastPoint
     !! Last point of workload
     integer :: work, remainder

     remainder = mod(points, procs)
     if(rank < remainder) then
        work = ceiling(dble(points) / dble(procs))
        firstPoint = rank * work + 1
        lastPoint = (rank + 1) * work
     else
        work = floor(dble(points) / dble(procs))
        firstPoint = rank * work + 1 + remainder
        lastPoint = (rank + 1) * work + remainder
     end if
     return
  end subroutine calcWorkload

end module mod_mpi_pars
