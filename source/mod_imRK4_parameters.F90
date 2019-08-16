!> module for system parameters.
module mod_imRK4_parameters
  use mod_f90_kind, only: double
  implicit none
  real(double) :: hE_0, hw1_m
  !! Intensity of electric and magnetic fields
  real(double) :: hw_e, hw_m
  !! Frequency (hwt) of electric and magnetic fields
  real(double) :: tau_e, tau_m
  !! Pulse length of electric and magnetic fields
  real(double) :: delay_e, delay_m
  !! Time delay in electric and magnetic fields pulses
  logical      :: lelectric, lmagnetic, lpulse_e, lpulse_m
  !! Logical variables for choosing between (oscilatory or laser fields)
  real(double) :: omega
  !! Mix of frequencies
  real(double) :: integration_time
  !! Real integration time 
  real(double) :: step
  !! Step size
  real(double) :: sc_tol
  !! Time propagation self consistency tolerence 
  real(double) :: abs_tol, rel_tol, Delta
  !! Step size control error(ERR) tolerence
  integer      :: time
  !! Integer integration time
  integer      :: dimH2
  !! Dimension: 2*dimension of the Hamiltonian (dimH)
  real(double) :: ERR
  !! Error for the calculation of the step size in time propagation
  real(double) :: field_direction_m(3), field_direction_e(3)
  !! Direction of the magnetic and electric field pulses
  real(double) :: time_conv = 6.582d-7
  !! Conversion of time units to picosecond
end module mod_imRK4_parameters
