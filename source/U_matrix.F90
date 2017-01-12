! Mounts the effective electron-electron interaction Hamiltonian
subroutine U_matrix(hee)
  use mod_f90_kind
  use mod_constants
  use mod_magnet
  use mod_parameters
  implicit none
  integer :: i,mu,nu
  complex(double), dimension(Npl+2,18,18), intent(out) :: hee

  ! the effective electronic interaction matrix
  hee = zero

  ! Diagonal terms (in orbital)
  do i=1,Npl
    do mu=5,9
      nu=mu+9
      hee(i+1,mu,mu) = eps1(i)-hdel(i)
      hee(i+1,nu,nu) = eps1(i)+hdel(i)
      hee(i+1,mu,nu) = -hdelm(i)
      hee(i+1,nu,mu) = -hdelp(i)
    end do
  end do

  return
end subroutine U_matrix