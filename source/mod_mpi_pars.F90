module mod_mpi_pars
  use mod_kind, only: dp, int32, int64
  use MPI
  ! use MPI_f08
  implicit none
  integer(int32) :: ierr, errorcode
  integer(int32) :: myrank,numprocs
  integer(int32), dimension(MPI_STATUS_SIZE) :: stat
  ! MPI_f08:
  ! TYPE(MPI_Status) :: stat

  integer(int32) :: FieldComm
  ! MPI_f08:
  ! type(MPI_Comm) :: FieldComm
  !! Magnetic Field Communicator
  integer(int32) :: sField
  !! Size of Magnetic Field Communicator
  integer(int32) :: rField
  !! Rank in Magnetic Field Communicator

  integer(int32), dimension(2) :: FreqComm
  ! MPI_f08:
  ! type(MPI_Comm), dimension(2) :: FreqComm
  !! Frequency Communicators, (1) Row, (2) Col (for sending to 1 Process)
  integer(int32), dimension(2) :: sFreq
  !! Size of Frequency Communicators
  integer(int32), dimension(2) :: rFreq
  !! Rank in Frequency Communicators

  integer(int64) :: startField
  !! First Point of Magnetic Field workload
  integer(int64) :: endField
  !! Last Point of Magnetic Field workload
  integer(int64) :: startFreq
  !! First Point of Frequency workload
  integer(int64) :: endFreq
  !! Last Point of Frequency workload

contains

  subroutine Initialize_MPI()
    implicit none

    integer(int32) :: provided

#ifdef _OPENMP
      call MPI_Init_thread(MPI_THREAD_FUNNELED,provided,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
#else
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,numprocs,ierr)
#endif
  end subroutine

  subroutine abortProgram(str)
    implicit none
    character(len=*), intent(in) :: str
    character(len=15) :: filename
    integer           :: unit=123456787
    !! Abort Message

    write(filename, "('error.',i0)") myrank
    open (unit=unit, file=trim(filename), status='replace', form='formatted')
    write(unit=unit,fmt="(a)") str
    close(unit=unit)
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)

  end subroutine abortProgram

  function genFieldComm(nFields, FieldID, comm)
    implicit none
    integer(int32), intent(in)      :: nFields
    !! Number of field points to be calculated in parallel
    integer(int32), intent(out)     :: FieldID
    !! Communicator ID
    integer(int32), intent(in) :: comm
    ! MPI_f08:
    ! type(MPI_Comm), intent(in) :: comm
    !! MPI Communicator to split up

    integer(int32) :: genFieldComm
    integer(int32) :: rank, procs, ierror
    integer(int32) :: group

    ! MPI_f08:
    ! type(MPI_Comm) :: genFieldComm
    ! integer(int32)      :: rank, procs, ierror
    ! integer(int32)      :: group

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs, ierror)
    if(mod(procs,nFields) /= 0) call abortProgram("[genFieldComm] Number of MPI processes have to be divisible by nFields")
    group = procs / nFields
    FieldID = rank / group
    call MPI_Comm_split(comm, FieldID, rank, genFieldComm, ierror)
  end function genFieldComm

  function genFreqComm(nFreq, FreqID, comm)
    implicit none
    integer(int32),        intent(in) :: nFreq
    !! Number of frequency points to be calculated in parallel
    integer(int32),       intent(out) :: FreqID
    !! Communicator ID
    integer(int32), intent(in) :: comm
    ! MPI_f08:
    ! type(MPI_Comm),   intent(in) :: comm
    
    !! MPI Communicator to split up
    integer(int32), dimension(2) :: genFreqComm
    ! MPI_f08:
    ! type(MPI_Comm), dimension(2) :: genFreqComm

    integer(int32) :: rank, rank_row, procs, ierror
    integer(int32) :: group

    call MPI_Comm_rank(comm, rank, ierror)
    call MPI_Comm_size(comm, procs,ierror)
    if(mod(procs,nFreq) /= 0) call abortProgram("[genFreqComm] Number of MPI processes have to be divisble by nFreq")
    group = procs / nFreq
    FreqID = rank / group
    call MPI_Comm_split(comm, FreqID, rank, genFreqComm(1), ierror)

    call MPI_Comm_rank(genFreqComm(1), rank_row, ierror)
    call MPI_Comm_split(comm, rank_row, rank, genFreqComm(2), ierror)
  end function genFreqComm

  subroutine genMPIGrid(nFields, nFieldPoints, nFreq, nFreqPoints)
    implicit none
    integer(int32), intent(in) :: nFields
    integer(int32), intent(in) :: nFreq
    integer(int32), intent(in) :: nFieldPoints
    integer(int32), intent(in) :: nFreqPoints
    integer(int32) :: FreqID, FieldID

    integer(int32) :: i, ierror

    FieldComm = genFieldComm(int(nFields,4), FieldID, MPI_COMM_WORLD)
    FreqComm = genFreqComm(int(nFreq,4), FreqID, FieldComm)
    call MPI_Comm_rank(FieldComm, rField, ierror)
    call MPI_Comm_size(FieldComm, sField, ierror)

    do i = 1,2
      call MPI_Comm_rank(FreqComm(i), rFreq(i), ierror)
      call MPI_Comm_size(FreqComm(i), sFreq(i), ierror)
    end do

    ! Calculate Field workload for each Field Communicator
    call calcWorkload(int(nFieldPoints,8),int(nFields,4),FieldID,startField,endField)
    ! Calculate Field workload for each Frequency Communicator
    call calcWorkload(int(nFreqPoints,8),int(nFreq,4),FreqID,startFreq,endFreq)
  end subroutine genMPIGrid

  subroutine calcWorkload(points, procs, rank, firstPoint, lastPoint)
    implicit none
    integer(int64), intent(in) :: points
    !! Number of points to split up
    integer(int32), intent(in) :: procs
    !! Number of processes
    integer(int32), intent(in) :: rank
    !! MPI Rank
    integer(int64), intent(out) :: firstPoint
    !! First point of workload
    integer(int64), intent(out) :: lastPoint
    !! Last point of workload
    integer(int64) :: work, remainder

    remainder = mod( points, int(procs,kind(points)) )
    if( int( rank,kind(remainder) ) < remainder ) then
      work = ceiling( dble(points)/dble(procs) , kind(work) )
      firstPoint = int(rank,kind(firstPoint)) * work + 1
      lastPoint = (int(rank,kind(lastPoint)) + 1) * work
    else
      work = floor( dble(points)/dble(procs) , kind(work) )
      firstPoint = int(rank,kind(firstPoint)) * work + 1 + remainder
      lastPoint = (int(rank,kind(lastPoint)) + 1) * work + remainder
    end if
  end subroutine calcWorkload

end module mod_mpi_pars
