! Spin Zeeman hamiltonian
subroutine sb_matrix()
  use mod_f90_kind
  use mod_constants
  use mod_parameters, only: Npl_total
  use mod_magnet
  implicit none
  integer :: i,mu,nu

! There is an extra  minus sign in the definition of hhwx,hhwy,hhwz
! to take into account the fact that we are considering negative
! external fields to get the peak at positive energies
  sb = zero
  do i=1,Npl_total ; do mu=1,9
    nu=mu+9
    sb(i,mu,mu) = hhwz(i)
    sb(i,nu,nu) =-hhwz(i)
    sb(i,mu,nu) = hhwx(i)-zi*hhwy(i)
    sb(i,nu,mu) = hhwx(i)+zi*hhwy(i)
  end do ; end do

  return
end subroutine sb_matrix
