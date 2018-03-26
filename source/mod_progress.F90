module mod_progress
  use mod_f90_kind
  implicit none
  real(double)             :: start_program,start_time,elapsed_time
contains

  subroutine progress_bar(unit,message,current_point,total_points)
    use mod_mpi_pars
    implicit none
    integer         ,intent(in) :: unit,current_point,total_points
    character(len=*),intent(in) :: message
    character(len=2)   :: spiner(4) = ['| ','/ ','- ','\ ']
    integer            :: porcent

    porcent = floor(current_point*100.d0/total_points)
    elapsed_time = MPI_Wtime() - start_program
    !write(unit,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1,$)") spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),char(13)
    write(unit,"(a1,2x,i3,'% (',i0,'/',i0,') of ',a,' done. Total time: ',i2,'h:',i2,'m:',i2,'s',a1)", advance='no') spiner(mod(current_point,4)+1),porcent,current_point,total_points,trim(message),int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),char(13)

  end subroutine progress_bar

  subroutine write_time(unit,message)
    use mod_mpi_pars
    implicit none
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    integer           :: values(8)
    character(len=*),intent(in) :: message
    integer         ,intent(in) :: unit

    call date_and_time(date, time, zone, values)
    elapsed_time = MPI_Wtime() - start_program
    write(unit,"(a,' ',i0,'/',i0,'/',i0,' at ',i2.2,'h',i2.2,'m',i2.2,'s (Elapsed time: ',f8.1,'s / ',f7.2,'m / ',f6.3,'h)')") trim(message),values(3),values(2),values(1),values(5),values(6),values(7),elapsed_time,elapsed_time/60.d0,elapsed_time/3600.d0

  end subroutine write_time

end module mod_progress
