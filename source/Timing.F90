module Timing
use mod_f90_kind, only: double
enum, bind(C)
enumerator :: INIT_TIME = 1
enumerator :: SC_TIME = 2
enumerator :: MAG_TIME = 3
enumerator :: JAC_TIME = 4
enumerator :: CHI_TIME = 5
enumerator :: TOTAL_TIME = 6
end enum

real(double), dimension(6), PRIVATE :: t0 = 0.d0, time = 0.d0

contains

  subroutine start_time(index)
    implicit none
    integer, intent(in) :: index

    call cpu_time(t0(index))
    return
  end subroutine start_time

  subroutine end_time(index)
    use mod_f90_kind, only: double
    implicit none
    integer, intent(in) :: index
    real(double) :: t

    call cpu_time(t)
    time(index) = t - t0(index)
    return
  end subroutine end_time

  subroutine print_time()
    implicit none

    print *, "-----------------------------Timings------------------------"
    print *, "Time for initialization: ", time(INIT_TIME)
    print *, "Time for Self-Consistency: ", time(SC_TIME)
    print *, "Time for Magnetization: ", time(MAG_TIME)
    print *, "Time for Jacobian: ", time(JAC_TIME)
    print *, "Time for Susceptibility: ", time(CHI_TIME)
    print *, "Total time: ", time(TOTAL_TIME)

    return
  end subroutine print_time
end module Timing
