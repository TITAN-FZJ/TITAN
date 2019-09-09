!> module for Runge-Kutta method parameters and Butcher tableu, Pauli matricies,
!> A_inverse matrix, M1 matrix, identity matricies, coefficients b's, d's, c's and c_avg.
module mod_RK_matrices 
  use mod_f90_kind, only: double
  implicit none
  complex(double), dimension(:,:), allocatable :: id(:,:),id2(:,:)
  ! real(double),    dimension(4,4), parameter   :: M1= reshape([ 0.25d0,0.d0,-0.03867513d0, 0.d0,0.d0,0.25d0,0.d0,-0.03867513d0,0.53867513d0,0.d0,0.25d0,0.d0,0.d0,0.53867513d0,0.d0,0.25d0 ],[ 4,4 ],Order=[ 2,1 ])
  complex(double), dimension(:,:), allocatable :: M1
  real(double),    dimension(2,2), parameter   :: A_inverse= reshape([ 3.d0, 0.46410162d0, -6.46410162d0, 3.d0 ],[ 2,2 ],Order=[ 2,1 ])
  complex(double), dimension(2,2), parameter   :: A= reshape([ 0.25d0, -0.03867513d0, 0.53867513d0, 0.25d0 ],[ 2,2 ],Order=[ 2,1 ])
  real(double),    dimension(2)  , parameter   :: b= [ 0.5d0, 0.5d0 ]
  real(double)                   , parameter   :: d1= -1.73205081d0
  real(double)                   , parameter   :: d2=  1.73205081d0
  real(double)                   , parameter   :: d1_hat= 6.4641018d0
  real(double)                   , parameter   :: d2_hat= -0.46410171d0
  real(double)                   , parameter   :: c1= 0.21132486540518713d0
  real(double)                   , parameter   :: c2= 0.5288675134594812d0
  real(double)                   , parameter   :: c_avg= (c1+c2)/2.d0
contains

  ! subroutine for identity matricies of dimension dim_I. 
  subroutine build_identity(dim_I,ident)
    use mod_f90_kind, only: double
    implicit none
    integer,                                 intent(in)  :: dim_I
    complex(double), dimension(dim_I,dim_I), intent(out) :: ident
    integer :: n

    ident = 0.d0
    do n= 1, dim_I
      ident(n,n)= 1.d0
    end do    

  end subroutine build_identity

end module mod_RK_matrices
