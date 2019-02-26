!! Fermi-Dirac distribution
function fd_dist(ef,beta,ek)
  use mod_f90_kind
  implicit none
  real(double), intent(in)  :: ef
!! Fermi energy
  real(double), intent(in)  :: beta
!! Inverse temperature
  real(double), intent(in)  :: ek
!!  Energy
  real(double)              :: fd_dist
!! Fermi-Dirac distribution
! ----------------------------------------------------------------------
  real(double) :: x

  x = 0.5d0*beta*(ef-ek)
  if (abs(x) < 10.d0) then
    fd_dist = 0.5d0*(1.d0 + tanh(x))
  else if (x > 0.d0) then
    fd_dist = 1.d0
  else
    fd_dist = 0.d0
  end if
end function fd_dist