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

    allocate(hee(2*nOrb,2*nOrb,nAtoms), stat=AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocate_Umatrix] Failed to allocate 'hee'.")

    return
  end subroutine allocate_Umatrix

  subroutine update_Umatrix(eps1, hdel, hdelm, hdelp, ndel, nAtoms, nOrb)
    use mod_f90_kind, only: double
    use mod_constants, only: cZero
    use mod_parameters, only: offset
    use mod_System, only: s => sys
    use mod_parameters, only: U
    implicit none

    integer, intent(in) :: nAtoms, nOrb
    real(double), dimension(nAtoms), intent(in) :: eps1, hdel
    complex(double), dimension(nAtoms), intent(in) :: hdelm, hdelp
    real(double), dimension(nAtoms), intent(in) :: ndel

    integer :: i, mu, nu

    hee = cZero
    ! Diagonal terms (in orbital)
    do i=1,nAtoms
      do mu=5, nOrb
        nu=mu+nOrb
        hee(mu,mu,i+offset) = eps1(i)-hdel(i) + ndel(i) - 0.5d0*U(i)*s%Types(s%Basis(i)%Material)%Occupation
        hee(nu,nu,i+offset) = eps1(i)+hdel(i) + ndel(i) - 0.5d0*U(i)*s%Types(s%Basis(i)%Material)%Occupation
        hee(mu,nu,i+offset) = -hdelm(i)
        hee(nu,mu,i+offset) = -hdelp(i)
      end do
    end do
    return
  end subroutine update_Umatrix

  subroutine init_Umatrix(eps1, hdel, hdelm, hdelp, ndel, nAtoms, nOrb)
    use mod_f90_kind, only: double
    use mod_mpi_pars, only: abortProgram
    implicit none

    integer, intent(in) :: nAtoms, nOrb
    real(double), dimension(nAtoms), intent(in) :: eps1, hdel
    complex(double), dimension(nAtoms), intent(in) :: hdelm, hdelp
    real(double), dimension(nAtoms), intent(in) :: ndel

    if(nOrb /= 9) call abortProgram("[init_Umatrix] Umatrix only implemented for nOrb = 9")

    call allocate_Umatrix(nAtoms,nOrb)
    call update_Umatrix(eps1,hdel,hdelm,hdelp,ndel,nAtoms,nOrb)

    return
  end subroutine

end module mod_Umatrix
