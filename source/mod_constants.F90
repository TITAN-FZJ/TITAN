module mod_constants
  use mod_f90_kind
  implicit none
  real(double) :: pi, tpi, sq2, hsq2, sq3
  complex(double), parameter :: zero=(0.d0,0.d0), zum=(1.d0,0.d0), zi=(0.d0,1.d0)
  complex(double) :: identorb18(18,18),identorb9(9,9)

contains
  subroutine define_constants()

    integer :: i

    pi  = 4.d0*atan(1.d0)
    tpi = 2.d0*pi
    sq2 = sqrt(2.d0)
    hsq2 = 0.5d0*sq2
    sq3 = sqrt(3.d0)
    identorb9  = zero
    identorb18 = zero
    do i=1,9
      identorb9(i,i) = zum
    end do
    do i=1,18
      identorb18(i,i) = zum
    end do

    return
  end subroutine define_constants
end module mod_constants