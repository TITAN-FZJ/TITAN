!> module for system parameters.
module mod_imRK4_parameters
  use mod_f90_kind, only: double
  use mod_constants
  implicit none
  real(double) :: hw1, hw, omega, integration_time, step, sc_tol
  integer      :: N, time
  integer      :: dimH2
  ! Dimension: 2*dimension of the Hamiltonian (dimH)

  ! namelist /tp_parameters/ hw1, hw, N, integration_time, sc_tol
  
  contains
  
  subroutine initialize()
    use mod_RK_matrices, only: A, build_identity
    use mod_parameters,  only: dimH
    use mod_magnet,      only: hwa
    implicit none
    integer              :: dim_I
    integer              :: ios

    omega= sqrt(hw1**2+ (hwa-0.5*hw)**2) 
    step= (2*pi)/(N*max(hw, omega))
    time= int(integration_time/step)


    ! Building identity
    dim_I = size(A,1)*dimH
    call build_identity(dim_I) 
  end subroutine initialize

end module mod_imRK4_parameters
