! Mounts the effective electron-electron interaction Hamiltonian
subroutine U_matrix(hee)
  use mod_f90_kind,      only: double
  use mod_constants,     only: zero
  use mod_magnet,        only: eps1, hdel, hdelm, hdelp
  use mod_parameters,    only: Npl, Npl_total
  use mod_tight_binding, only: tbmode
  implicit none
  integer :: i,mu,nu
  complex(double), dimension(Npl_total,18,18), intent(out) :: hee

  ! the effective electronic interaction matrix
  hee = zero

  ! Diagonal terms (in orbital)
  if(tbmode == 2) then
    do i=1,Npl
      do mu=5,9
        nu=mu+9
        hee(i+1,mu,mu) = eps1(i)-hdel(i)
        hee(i+1,nu,nu) = eps1(i)+hdel(i)
        hee(i+1,mu,nu) = -hdelm(i)
        hee(i+1,nu,mu) = -hdelp(i)
      end do
    end do
  else
    do i=1,Npl
      do mu=5,9
        nu=mu+9
        hee(i,mu,mu) = eps1(i)-hdel(i)
        hee(i,nu,nu) = eps1(i)+hdel(i)
        hee(i,mu,nu) = -hdelm(i)
        hee(i,nu,mu) = -hdelp(i)
      end do
    end do
  end if
  return
end subroutine U_matrix
