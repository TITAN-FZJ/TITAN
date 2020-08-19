!> module for Runge-Kutta method parameters and Butcher tableu, Pauli matricies,
!> A_inverse matrix, M1 matrix, identity matricies, coefficients b's, d's, c's and c_avg.
module mod_RK_matrices 
  use mod_kind, only: dp
  implicit none
  complex(dp), dimension(:,:), allocatable :: id(:,:),id2(:,:)
  ! real(dp),    dimension(4,4), parameter   :: M1= reshape([ 0.25_dp,0._dp,-0.03867513_dp, 0._dp,0._dp,0.25_dp,0._dp,-0.03867513_dp,0.53867513_dp,0._dp,0.25_dp,0._dp,0._dp,0.53867513_dp,0._dp,0.25_dp ],[ 4,4 ],Order=[ 2,1 ])
  complex(dp), dimension(:,:), allocatable :: M1
  real(dp),    dimension(2,2), parameter   :: A_inverse= reshape([ 3._dp, 0.46410162_dp, -6.46410162_dp, 3._dp ],[ 2,2 ],Order=[ 2,1 ])
  complex(dp), dimension(2,2), parameter   :: A= reshape([ 0.25_dp, -0.03867513_dp, 0.53867513_dp, 0.25_dp ],[ 2,2 ],Order=[ 2,1 ])
  real(dp),    dimension(2)  , parameter   :: b= [ 0.5_dp, 0.5_dp ]
  real(dp)                   , parameter   :: d1= -1.73205081_dp
  real(dp)                   , parameter   :: d2=  1.73205081_dp
  real(dp)                   , parameter   :: d1_hat= 6.4641018_dp
  real(dp)                   , parameter   :: d2_hat= -0.46410171_dp
  real(dp)                   , parameter   :: c1= 0.21132486540518713_dp
  real(dp)                   , parameter   :: c2= 0.5288675134594812_dp
  real(dp)                   , parameter   :: c_avg= (c1+c2)/2._dp
contains

  ! subroutine for identity matricies of dimension dim_I. 
  subroutine build_identity(dim_I,ident)
    use mod_kind, only: dp
    implicit none
    integer,                                 intent(in)  :: dim_I
    complex(dp), dimension(dim_I,dim_I), intent(out) :: ident
    integer :: n

    ident = 0._dp
    do n= 1, dim_I
      ident(n,n)= 1._dp
    end do    

  end subroutine build_identity

end module mod_RK_matrices
