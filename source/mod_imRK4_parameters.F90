!> module for system parameters.
module mod_imRK4_parameters
  use mod_f90_kind, only: double
  implicit none
  logical      :: lelectric, lmagnetic, lpulse_e, lpulse_m
  !! Logical variables for choosing which field is applied
  logical      :: npulse_e, npulse_m
  !! Number of pulses (electric and magnetic)
  real(double) :: hE_0, hw1_m
  !! Intensity of electric and magnetic fields
  real(double) :: hw_e, hw_m
  !! Frequency (hwt) of electric and magnetic fields
  real(double) :: tau_e, tau_m
  !! Pulse length of electric and magnetic fields
  real(double) :: delay_e, delay_m
  !! Time delay in electric and magnetic fields pulses
  real(double) :: integration_time
  !! Real integration time 
  real(double) :: step
  !! Step size
  real(double) :: sc_tol
  !! Time propagation self consistency tolerence 
  real(double) :: abs_tol, rel_tol, safe_factor
  !! Step size control error(ERR) tolerance
  integer      :: time
  !! Integer integration time
  integer      :: dimH2
  !! Dimension: 2*dimension of the Hamiltonian (dimH)
  real(double) :: ERR
  !! Error for the calculation of the step size in time propagation
  real(double) :: field_direction_m(3)
  !! Direction of the magnetic field
  real(double) :: time_conv = 6.582d-7
  !! Conversion of time units to picosecond
  real(double)         :: polarization_vector_e(3)
  !! Polarization vector of electric field 
  character(len = 20)  :: shape_e
  !! Shape of the electric pulse
  character(len = 20)  ::  polarization_type_e 
  !! Electric field polarization type 
  real(double)         :: helicity
  !! Helicity of the circular electric field

end module mod_imRK4_parameters
