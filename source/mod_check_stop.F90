module mod_check_stop
  implicit none
  character(len=30) :: stopfilename=""    !! Filename for the new stop file

contains

! This subroutine check for stop file 'filename'
! and create a unique subfile when more than one iteration is required.
! This file counts down the remaining number of iterations to be made. 
! Since the unique file is read at every iteration,
! it can be editted while the program is running.
! This is useful to run different jobs with the same loop
! (in conjunction with variable skip_steps_hw)
  subroutine check_stop(filename,hw_count,e)
    use mod_kind,       only: dp
    use mod_parameters, only: output
    use mod_magnet,     only: hw_list
    use mod_mpi_pars,   only: abortProgram,myrank,MPI_CHARACTER,MPI_COMM_WORLD,ierr
    implicit none
    integer,          intent(in) :: hw_count     !! Counter for magnetic field
    real(dp),         intent(in), optional :: e  !! Current frequency 
    character(len=*), intent(in) :: filename     !! Initial filename to check
    character(len=8)  :: cdate
    character(len=10) :: ctime
    character(len=5)  :: zone
    integer           :: values(8)
    integer           :: ios,istop
  
    external :: MPI_Bcast

    ! Check if stop file exists (if hasn't been read before)
    if(trim(stopfilename)=="") then
      open(unit=911, file=trim(filename), status='old', iostat=ios)
      if(ios/=0) return ! If it doesn't exist, return

      ! If exists, generate a unique filename at myrank=0 and broadcast to all
      ! (to be able to run more than one version of the program at the same time)
      call date_and_time(cdate, ctime, zone, values)
      if(myrank==0) write(stopfilename,fmt="(a,'_',i0,i0,i0,i2.2,i2.2,i2.2)") trim(stopfilename),values(3),values(2),values(1),values(5),values(6),values(7)
      call MPI_Bcast(stopfilename,30,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    else
      open(unit=911, file=trim(stopfilename), status='old', iostat=ios)
    end if

    ! Read number of iterations to be made before stopping
    read(unit=911,fmt=*,iostat=ios) istop 
    close(911)
    if((ios/=0).or.(istop<=1)) then
      ! If no number is given or if it reaches 0, the program stops
      if(myrank==0) then
        if(present(e)) then
          write(output%unit,"('[check_stop] Emergency ""',a,'"" file found! Stopping after e = ',es9.2,' was calculated...')") trim(filename),e
        else
          write(output%unit,"('[check_stop] Emergency ""',a,'"" file found! Stopping after hwa = ',es9.2,' hwt=',f7.2,' hwp=',f7.2,' was calculated...')") trim(filename),hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
        end if
        ! open(unit=911, file=trim(stopfilename), status='replace')
        ! write(911,"(i0)") istop-1 ! Removing one iteration of the new file
        ! close(911)
        ! call system ('rm stopout')
        ! write(outputunit,"('[check_stop] (""stopout"" file deleted!)')")
        call abortProgram("[check_stop] REMEMBER TO DELETE" // trim(filename) // "FILES!")
      end if

    else 
      ! If a number > 1 is given, counts down the number of iterations:
      if(myrank==0) then
        write(output%unit,"('[check_stop] Emergency ""',a,'"" file found! ',i0,' more iterations before stopping...')") trim(filename),istop-1
        open(unit=911, file=trim(stopfilename), status='replace')
        write(911,"(i0)") istop-1 ! Removing one iteration of the file
        close(911)
      end if
    end if
  end subroutine check_stop

end module mod_check_stop