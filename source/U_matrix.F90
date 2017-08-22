! Mounts the effective electron-electron interaction Hamiltonian
subroutine U_matrix(hee, nAtoms, nOrb)
  use mod_f90_kind,      only: double
  use mod_constants,     only: zero
  use mod_magnet,        only: eps1, hdel, hdelm, hdelp
  use mod_parameters, only: offset
  implicit none
  integer, intent(in) :: nAtoms, nOrb
  integer :: i,mu,nu
  complex(double), dimension(2*nOrb,2*nOrb,nAtoms), intent(out) :: hee

  ! the effective electronic interaction matrix
  hee = zero
  if(2*nOrb /= 18) stop "U_matrix only implemented for nOrb = 9"
  ! Diagonal terms (in orbital)
    do i=1,nAtoms
      do mu=5, 9
        nu=mu+9
        hee(mu,mu,i+offset) = eps1(i)-hdel(i)
        hee(nu,nu,i+offset) = eps1(i)+hdel(i)
        hee(mu,nu,i+offset) = -hdelm(i)
        hee(nu,mu,i+offset) = -hdelp(i)
      end do
    end do
  return
end subroutine U_matrix
