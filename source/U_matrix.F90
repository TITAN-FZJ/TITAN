! Mounts the effective electron-electron interaction Hamiltonian
subroutine U_matrix(hee)
  use mod_f90_kind,      only: double
  use mod_constants,     only: zero
  use mod_magnet,        only: eps1, hdel, hdelm, hdelp
  use mod_parameters,    only: Npl, Npl_total, offset
  implicit none
  integer :: i,mu,nu
  complex(double), dimension(Npl_total,18,18), intent(out) :: hee

  ! the effective electronic interaction matrix
  hee = zero

  ! Diagonal terms (in orbital)
    do i=1,Npl
      do mu=5,9
        nu=mu+9
        hee(i+offset,mu,mu) = eps1(i)-hdel(i)
        hee(i+offset,nu,nu) = eps1(i)+hdel(i)
        hee(i+offset,mu,nu) = -hdelm(i)
        hee(i+offset,nu,mu) = -hdelp(i)
      end do
    end do
  return
end subroutine U_matrix
