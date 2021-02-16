!> module for system parameters.
module mod_imRK4_parameters
  use mod_kind, only: dp
  implicit none
  logical      :: lelectric, lmagnetic, lpulse_e, lpulse_m
  !! Logical variables for choosing which field is applied

  integer      :: npulse_e
  !! Number of electric pulses
  real(dp),     dimension(:) , allocatable :: hE_0
  !! Intensity of electric fields
  real(dp),     dimension(:) , allocatable :: hw_e
  !! Frequency (hwt) of electric fields
  real(dp),     dimension(:) , allocatable :: tau_e
  !! Pulse length of electric fields
  real(dp),     dimension(:) , allocatable :: delay_e
  !! Time delay in electric fields pulses
  character(len=1), dimension(:) , allocatable ::  polarization_e 
  !! Polarization of electric field
  real(dp), dimension(:,:,:) , allocatable ::  polarization_vec_e
  !! Polarization vector, inphase (cos) and out-of-phase (sin) of electric field

  integer      :: npulse_m
  !! Number of magnetic pulses
  real(dp),     dimension(:) , allocatable :: hw1_m
  !! Intensity of magnetic fields
  real(dp),     dimension(:) , allocatable :: hw_m
  !! Frequency (w.t) of magnetic fields
  real(dp),     dimension(:) , allocatable :: tau_m
  !! Pulse length of magnetic fields
  real(dp),     dimension(:) , allocatable :: delay_m
  !! Time delay in magnetic fields pulses
  character(len=1), dimension(:) , allocatable ::  polarization_m
  !! Polarization of magnetic field
  real(dp), dimension(:,:,:) , allocatable ::  polarization_vec_m
  !! Polarization vector, inphase (cos) and out-of-phase (sin) of magnetic field

  real(dp) :: integration_time
  !! Real integration time 
  real(dp) :: step
  !! Step size
  real(dp) :: sc_tol
  !! Time propagation self consistency tolerence 
  real(dp) :: abs_tol, rel_tol, safe_factor
  !! Step size control error(ERR) tolerance
  integer  :: dimH2
  !! Dimension: 2*dimension of the Hamiltonian (2*dimHsc)
  real(dp) :: ERR
  !! Error for the calculation of the step size in time propagation
  real(dp) :: time_conv = 4.84e-5_dp
  !! Conversion of time units to picosecond

end module mod_imRK4_parameters
