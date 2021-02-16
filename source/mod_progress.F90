module mod_progress
  use mod_kind, only: dp
  implicit none
  real(dp) :: start_program,start_time,elapsed_time

  interface progress_bar
     module procedure progress_bar_int, &
                      progress_bar_double_int
  end interface progress_bar

contains

  subroutine progress_bar_int(message,current_point,total_points,unit)
    use mod_mpi_pars, only: MPI_Wtime
    implicit none
    integer         , intent(in) :: current_point,total_points
    character(len=*), intent(in) :: message
    integer         , intent(in), optional :: unit
    character(len=2)   :: spiner(4) = ['| ','/ ','- ','\ ']
    integer            :: porcent

    porcent = floor(current_point*100._dp/total_points)
    elapsed_time = MPI_Wtime() - start_program
    !write(unit,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1,$)") spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600._dp),int(mod(elapsed_time,3600._dp)/60._dp),int(mod(mod(elapsed_time,3600._dp),60._dp)),char(13)
    if(present(unit)) then
      write(unit,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1)", advance='no') spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600._dp),int(mod(elapsed_time,3600._dp)/60._dp),int(mod(mod(elapsed_time,3600._dp),60._dp)),char(13)
    else
      write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1)", advance='no') spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600._dp),int(mod(elapsed_time,3600._dp)/60._dp),int(mod(mod(elapsed_time,3600._dp),60._dp)),char(13)
    end if

  end subroutine progress_bar_int


  subroutine progress_bar_double_int(message,current_point,total_points,unit)
    use mod_kind,     only: int64
    use mod_mpi_pars, only: MPI_Wtime
    implicit none
    integer(int64)  , intent(in) :: current_point,total_points
    character(len=*), intent(in) :: message
    integer         , intent(in), optional :: unit
    character(len=2)   :: spiner(4) = ['| ','/ ','- ','\ ']
    integer            :: porcent

    porcent = floor(current_point*100._dp/total_points)
    elapsed_time = MPI_Wtime() - start_program
    !write(unit,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1,$)") spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600._dp),int(mod(elapsed_time,3600._dp)/60._dp),int(mod(mod(elapsed_time,3600._dp),60._dp)),char(13)
    if(present(unit)) then
      write(unit,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1)", advance='no') spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600._dp),int(mod(elapsed_time,3600._dp)/60._dp),int(mod(mod(elapsed_time,3600._dp),60._dp)),char(13)
    else
      write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1)", advance='no') spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600._dp),int(mod(elapsed_time,3600._dp)/60._dp),int(mod(mod(elapsed_time,3600._dp),60._dp)),char(13)
    end if

  end subroutine progress_bar_double_int


  subroutine write_time(message,unit)
    use mod_mpi_pars, only: MPI_Wtime
    implicit none
    character(len=8)  :: cdate
    character(len=10) :: ctime
    character(len=5)  :: zone
    integer           :: values(8)
    character(len=*), intent(in) :: message
    integer         , intent(in), optional :: unit

    call date_and_time(cdate,ctime,zone,values)
    elapsed_time = MPI_Wtime() - start_program
    if(present(unit)) then
      write(unit,"(a,' ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s (Elapsed time: ',f8.1,'s / ',f7.2,'m / ',f6.3,'h)')") trim(message),values(3),values(2),values(1),values(5),values(6),values(7),elapsed_time,elapsed_time/60._dp,elapsed_time/3600._dp
    else
      write(*,"(a,' ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s (Elapsed time: ',f8.1,'s / ',f7.2,'m / ',f6.3,'h)')") trim(message),values(3),values(2),values(1),values(5),values(6),values(7),elapsed_time,elapsed_time/60._dp,elapsed_time/3600._dp
    end if

  end subroutine write_time

end module mod_progress
