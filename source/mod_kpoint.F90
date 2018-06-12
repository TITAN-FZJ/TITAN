module mod_kpoint
  use mod_f90_kinds, only: double
  implicit none
  type :: kpoint_set
    integer :: size
    real(double), dimension(:), allocatable :: w
    real(double), dimension(:,3), allocatable :: p
  end type
  
end module