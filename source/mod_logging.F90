module mod_logging
  !! Subroutines to for logging variables and procedures
  implicit none
  logical :: log_unit = .false.
  character(len=12), parameter :: logfile = "parameter.in"
  character(len=:), allocatable :: log_store

contains

  subroutine log_message(proced, message, var)
    !! This subroutine logs a 'message' containing the procedure where it was called (passed via 'proced')
    !! if the string variable 'var' is present, the message is logged in it
    use mod_mpi_pars,   only: myrank
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: proced
    character(len=*), intent(in) :: message
    character(len=:), allocatable, intent(inout), optional :: var


    ! To store initial messages before output file is open
    if (present(var)) then
      var = trim(var) // "[Warning] [" // proced // "] "// trim(message) // NEW_LINE('a')
      return
    end if
    
    if(myrank == 0) then
      if(log_unit) then
        ! To write the stored messages without the square brackeds and skipping a line
        if (proced == "") then
          write(output%unit, "(a)", advance="no") trim(message)
        else
          write(output%unit, "('[',a,'] ',a)") proced, trim(message)
        end if
      else
        write(*, "('[',a,'] ',a)") proced, trim(message)
      end if
    end if
  end subroutine log_message


  subroutine log_warning(proced, message, var)
    !! This subroutine logs a warning 'message' containing the procedure where it was called (passed via 'proced')
    !! if the string variable 'var' is present, the message is logged in it
    use mod_mpi_pars,   only: myrank
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: proced
    character(len=*), intent(in) :: message
    character(len=:), allocatable, intent(inout), optional :: var

    ! To store initial messages before output file is open
    if (present(var)) then
      var = trim(var) // "[Warning] [" // proced // "] "// trim(message) // NEW_LINE('a')
      return
    end if

    if(myrank == 0) then
      if(log_unit) then
        write(output%unit, "('[Warning] [',a,'] ',a)") proced, trim(message)
      else
        write(*, "('[Warning] [',a,'] ',a)") proced, trim(message)
      end if
    end if
  end subroutine log_warning


  subroutine log_error(proced, message)
    !! This subroutine logs an error 'message' containing the procedure where it was called (passed via 'proced')
    !! and then aborts the run
    use mod_mpi_pars,   only: myrank,MPI_Abort,MPI_COMM_WORLD,errorcode,ierr
    use mod_parameters, only: output
    implicit none
    character(len=*), intent(in) :: proced
    character(len=*), intent(in) :: message

    if(myrank == 0) then
      if(log_unit) write(output%unit, "('[Error] [',a,'] ',a)") proced, trim(message)
    end if
    write(*, "('[Error] [',a,'] ',a)") proced, trim(message)
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    stop
  end subroutine log_error

end module mod_logging