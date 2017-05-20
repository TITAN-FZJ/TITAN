! Orbital Zeeman hamiltonian
subroutine lb_matrix()
  use mod_parameters, only: Npl_total
  use mod_constants,  only: zero
  use mod_f90_kind,   only: double
  use mod_magnet,     only: lxp, lyp, lzp, hhwx, hhwy, hhwz, lb
  implicit none
  integer :: i
  complex(double), dimension(9,9,Npl_total) :: lbsigma

! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
! to take into account the fact that we are considering negative
! external fields to get the peak at positive energies
  lb = zero
  do i=1,Npl_total
    lbsigma( 1:9,   1:9,  i) = 0.5d0*(lxp*hhwx(i) + lyp*hhwy(i) + lzp*hhwz(i))
    lb(      1:9,   1:9,  i) = lbsigma(:,:,i)
    lb(     10:18, 10:18, i) = lbsigma(:,:,i)
  end do

  return
end subroutine lb_matrix
