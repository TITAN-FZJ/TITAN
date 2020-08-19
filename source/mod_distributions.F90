!! This module contains the Fermi-Dirac distribution
module mod_distributions

contains

  !! Fermi-Dirac distribution
  pure elemental function fd_dist(ef,beta,ek)
    !$OMP DECLARE SIMD
    use mod_kind, only: dp
    implicit none
    real(dp), intent(in)  :: ef
  !! Fermi energy
    real(dp), intent(in)  :: beta
  !! Inverse temperature
    real(dp), intent(in)  :: ek
  !!  Energy
    real(dp)              :: fd_dist
  !! Fermi-Dirac distribution
  ! ----------------------------------------------------------------------
    real(dp) :: x

    x = 0.5_dp*beta*(ef-ek)
    if (abs(x) < 10._dp) then
      fd_dist = 0.5_dp*(1._dp + tanh(x))
    else if (x > 0._dp) then
      fd_dist = 1._dp
    else
      fd_dist = 0._dp
    end if
  end function fd_dist

end module mod_distributions