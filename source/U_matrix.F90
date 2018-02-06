module mod_Umatrix
  ! Mounts the effective electron-electron interaction Hamiltonian
  use mod_f90_kind, only: double
  implicit none

  complex(double), dimension(:,:,:), allocatable :: hee

contains

  subroutine allocate_Umatrix(nAtoms, nOrb)
    use mod_mpi_pars, only: AbortProgram
    implicit none

    integer, intent(in) :: nAtoms, nOrb
    integer :: AllocateStatus
    if(allocated(hee)) deallocate(hee)
    allocate(hee(2*nOrb,2*nOrb,nAtoms), stat=AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocate_Umatrix] Failed to allocate 'hee'.")

    return
  end subroutine allocate_Umatrix

  subroutine update_Umatrix(mz, mp, rho, Un0, Ucc0, nAtoms, nOrb)
    use mod_f90_kind, only: double
    use mod_constants, only: cZero
    use mod_parameters, only: offset, U
    use mod_parameters, only: U
    implicit none

    integer, intent(in) :: nAtoms, nOrb
    real(double), dimension(nAtoms), intent(in) :: mz
    complex(double), dimension(nAtoms), intent(in) :: mp
    real(double), dimension(nAtoms), intent(in) :: rho
    complex(double), dimension(nOrb,nOrb,nAtoms) :: Ucc0
    real(double), dimension(nAtoms)              :: Un0

    integer :: i, mu, nu

    hee = cZero
    do i=1,nAtoms
      ! Off-diagonal terms (in orbital)
      hee(     5:  nOrb,     5:  nOrb,i+offset) = Ucc0(5:nOrb,5:nOrb,i)
      hee(nOrb+5:2*nOrb,nOrb+5:2*nOrb,i+offset) = Ucc0(5:nOrb,5:nOrb,i)
      ! Diagonal terms (in orbital)
      do mu=5, nOrb
        nu=mu+nOrb
        hee(mu,mu,i+offset) = hee(mu,mu,i+offset)+0.5d0*U(i+offset) * (-mz(i) - rho(i) ) + Un0(i)
        hee(nu,nu,i+offset) = hee(nu,nu,i+offset)+0.5d0*U(i+offset) * ( mz(i) - rho(i) ) + Un0(i)
        hee(mu,nu,i+offset) = hee(mu,nu,i+offset)-0.5d0*U(i+offset)* conjg(mp(i))
        hee(nu,mu,i+offset) = hee(nu,mu,i+offset)-0.5d0*U(i+offset)* mp(i)
      end do
    end do
    return
  end subroutine update_Umatrix

  subroutine init_Umatrix(mz, mp, rho, Un0, Ucc0, nAtoms, nOrb)
    use mod_f90_kind, only: double
    use mod_mpi_pars, only: abortProgram
    implicit none

    integer, intent(in) :: nAtoms, nOrb
    real(double), dimension(nAtoms), intent(in) :: mz, rho
    complex(double), dimension(nAtoms), intent(in) :: mp
    complex(double), dimension(nOrb,nOrb,nAtoms) :: Ucc0
    real(double), dimension(nAtoms)              :: Un0

    if(nOrb /= 9) call abortProgram("[init_Umatrix] Umatrix only implemented for nOrb = 9")

    call allocate_Umatrix(nAtoms,nOrb)
    call update_Umatrix(mz,mp,rho,Un0,Ucc0,nAtoms,nOrb)

    return
  end subroutine

end module mod_Umatrix
