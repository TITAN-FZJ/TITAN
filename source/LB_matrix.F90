! Orbital Zeeman hamiltonian
subroutine lb_matrix()
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_magnet
  implicit none
  integer :: i
  complex(double), dimension(Npl+2,9,9) :: lbsigma

! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
! to take into account the fact that we are considering negative
! external fields to get the peak at positive energies
  lb = zero
  do i=1,Npl+2
    lbsigma(i,:,:) = 0.5d0*(lxp*hhwx(i) + lyp*hhwy(i) + lzp*hhwz(i))
    lb(i, 1: 9, 1: 9) = lbsigma(i,:,:)
    lb(i,10:18,10:18) = lbsigma(i,:,:)
  end do

  return
end subroutine lb_matrix
