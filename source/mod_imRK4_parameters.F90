!> module for system parameters.
module mod_imRK4_parameters
  use mod_f90_kind, only: double
  implicit none
  logical      :: lelectric, lmagnetic, lpulse_e, lpulse_m
  !! Logical variables for choosing which field is applied

  integer      :: npulse_e
  !! Number of electric pulses
  real(double),     dimension(:) , allocatable :: hE_0
  !! Intensity of electric fields
  real(double),     dimension(:) , allocatable :: hw_e
  !! Frequency (hwt) of electric fields
  real(double),     dimension(:) , allocatable :: tau_e
  !! Pulse length of electric fields
  real(double),     dimension(:) , allocatable :: delay_e
  !! Time delay in electric fields pulses
  character(len=1), dimension(:) , allocatable ::  polarization_e 
  !! Polarization of electric field
  real(double), dimension(:,:,:) , allocatable ::  polarization_vec_e
  !! Polarization vector, inphase (cos) and out-of-phase (sin) of electric field

  integer      :: npulse_m
  !! Number of magnetic pulses
  real(double),     dimension(:) , allocatable :: hw1_m
  !! Intensity of magnetic fields
  real(double),     dimension(:) , allocatable :: hw_m
  !! Frequency (w.t) of magnetic fields
  real(double),     dimension(:) , allocatable :: tau_m
  !! Pulse length of magnetic fields
  real(double),     dimension(:) , allocatable :: delay_m
  !! Time delay in magnetic fields pulses
  character(len=1), dimension(:) , allocatable ::  polarization_m
  !! Polarization of magnetic field
  real(double), dimension(:,:,:) , allocatable ::  polarization_vec_m
  !! Polarization vector, inphase (cos) and out-of-phase (sin) of magnetic field

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
  real(double) :: time_conv = 4.84d-5
  !! Conversion of time units to picosecond

end module mod_imRK4_parameters
