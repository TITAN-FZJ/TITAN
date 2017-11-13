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

  integer*8 :: startField
  !! First Point of Magnetic Field workload
  integer*8 :: endField
  !! Last Point of Magnetic Field workload
  integer*8 :: startFreq
  !! First Point of Frequency workload
  integer*8 :: endFreq
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
     call calcWorkload(int(nFieldPoints,8),nFields,FieldID,startField,endField)
     ! Calculate Field workload for each Frequency Communicator
     call calcWorkload(int(nFreqPoints,8),nFreq,FreqID,startFreq,endFreq)
     return
  end subroutine genMPIGrid

  subroutine calcWorkload(points, procs, rank, firstPoint, lastPoint)
     implicit none
     integer*8, intent(in) :: points
     !! Number of points to split up
     integer, intent(in) :: procs
     !! Number of processes
     integer, intent(in) :: rank
     !! MPI Rank
     integer*8, intent(out) :: firstPoint
     !! First point of workload
     integer*8, intent(out) :: lastPoint
     !! Last point of workload
     integer*8 :: work, remainder

     remainder = mod(points, int(procs,8))
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
