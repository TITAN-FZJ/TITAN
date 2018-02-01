! This subroutine check for stop file
! and create unique subfiles when more than one iteration is required
! This is useful to run different jobs with the same loop
! (in conjunction with variable skip_steps_hw)
subroutine check_emergency_stop(nAtoms, hw_list, hw_count)
  use mod_parameters, only: output
  use mod_mpi_pars
  implicit none
  integer, intent(in) :: nAtoms, hw_count
  real(double), dimension(hw_count,3), intent(in) :: hw_list
  character(len=8)  :: date
  character(len=10) :: time
  character(len=5)  :: zone
  integer           :: values(8)
  character(len=30) :: stopfilename="stopout"
  integer           :: ios,istop

  open(unit=911, file=trim(stopfilename), status='old', iostat=ios)
  if(ios==0) then
    read(unit=911,fmt=*,iostat=ios) istop ! Checking number of iterations to be made before stopping
    close(911)
    if((istop<=1).or.(ios/=0)) then
      if(myrank==0) then
        write(output%unit,"('[main] Emergency ""stopout"" file found! Stopping after Npl = ',i0,', hwa = ',es9.2,' hwt=',f7.2,' hwp=',f7.2,'...')") nAtoms,hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
        open(unit=911, file=trim(stopfilename), status='replace')
        write(911,"(i0)") istop-1 ! Removing one iteration of the file
        close(911)
        ! call system ('rm stopout')
        ! write(outputunit,"('[main] (""stopout"" file deleted!)')")
        write(*,"('[main] REMEMBER TO DELETE FILE ""',a,'""!')") trim(stopfilename)
      end if
      call MPI_Finalize(ierr)
      stop
    else
      ! Generating a unique filename at myrank=0 (to be able to run more than one version of the program at the same time)
      call date_and_time(date, time, zone, values)
      if(trim(stopfilename)=="stopout") then
        if(myrank==0) write(stopfilename,fmt="(a,'_',i0,i0,i0,i2.2,i2.2,i2.2)") trim(stopfilename),values(3),values(2),values(1),values(5),values(6),values(7)
        call MPI_Bcast(stopfilename,30,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      end if

      if(myrank==0) then
        write(output%unit,"('[main] Emergency ""stopout"" file found! ',i0,' more iterations before stopping...')") istop-1
        open(unit=911, file=trim(stopfilename), status='replace')
        write(911,"(i0)") istop-1 ! Removing one iteration of the file
        close(911)
      end if
    end if
  end if
end subroutine check_emergency_stop
