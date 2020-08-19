module mod_Umatrix
  ! Mounts the effective electron-electron interaction Hamiltonian
  use mod_kind, only: dp
  implicit none

  complex(dp), dimension(:,:,:), allocatable :: hee

contains

  subroutine allocate_Umatrix(nAtoms,nOrb)
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms, nOrb
    integer             :: AllocateStatus

    if(allocated(hee)) deallocate(hee)
    allocate(hee(2*nOrb,2*nOrb,nAtoms), stat=AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocate_Umatrix] Failed to allocate 'hee'.")

  end subroutine allocate_Umatrix

  subroutine deallocate_Umatrix()
    implicit none

    if(allocated(hee)) deallocate(hee)

  end subroutine deallocate_Umatrix

  subroutine update_Umatrix(mz,mp,rhod,rhod0,rho,rho0,nAtoms,nOrb)
    use mod_kind, only: dp
    use mod_constants,  only: cZero
    use mod_parameters, only: Un, Um
    implicit none
    integer,                              intent(in) :: nAtoms, nOrb
    real(dp),    dimension(nAtoms),   intent(in) :: mz,rhod,rhod0
    complex(dp), dimension(nAtoms),   intent(in) :: mp
    real(dp), dimension(nOrb,nAtoms), intent(in) :: rho,rho0
    integer :: i, mu, nu

    hee = cZero
    do i=1,nAtoms
      do mu=5, nOrb
        nu=mu+nOrb
        hee(mu,mu,i) = Un(i)*(rho(mu,i) - rho0(mu,i)) + 0.5_dp*Un(i)*( - rhod(i) + rhod0(i)) &
                       + 0.5_dp*Um(i)*(-mz(i))
        hee(nu,nu,i) = Un(i)*(rho(mu,i) - rho0(mu,i)) + 0.5_dp*Un(i)*( - rhod(i) + rhod0(i)) &
                       + 0.5_dp*Um(i)*( mz(i))
        hee(mu,nu,i) = - 0.5_dp*Um(i)*conjg(mp(i))
        hee(nu,mu,i) = - 0.5_dp*Um(i)*mp(i)
      end do
    end do
  end subroutine update_Umatrix

  subroutine init_Umatrix(mz,mp,rhod,rhod0,rho,rho0,nAtoms,nOrb)
    use mod_kind, only: dp
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer,                              intent(in) :: nAtoms,nOrb
    real(dp),    dimension(nAtoms),   intent(in) :: mz,rhod,rhod0
    complex(dp), dimension(nAtoms),   intent(in) :: mp
    real(dp), dimension(nOrb,nAtoms), intent(in) :: rho,rho0

    if(nOrb /= 9) call abortProgram("[init_Umatrix] Umatrix only implemented for nOrb = 9")

    call allocate_Umatrix(nAtoms,nOrb)
    call update_Umatrix(mz,mp,rhod,rhod0,rho,rho0,nAtoms,nOrb)

  end subroutine init_Umatrix
end module mod_Umatrix
