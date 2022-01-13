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

  subroutine update_Umatrix(mzd,mzd0,mpd,mpd0,rhod,rhod0,rho,rho0,s)
  !! This subroutine receives the magnetic moment and occupation per orbital
  !! and per site and updates the Hubbard term of the Hamiltonian
    use mod_kind,       only: dp
    use mod_constants,  only: cZero
    use mod_System,     only: System_type
    implicit none
    type(System_type),                       intent(in) :: s
    real(dp),    dimension(s%nAtoms),        intent(in) :: mzd,mzd0,rhod,rhod0
    complex(dp), dimension(s%nAtoms),        intent(in) :: mpd,mpd0
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(in) :: rho,rho0
    integer :: i,mu,nud,mud

    hee = cZero
    do i=1,s%nAtoms
      do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
        mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
        nud=mud+s%Types(s%Basis(i)%Material)%nOrb
        hee(mud,mud,i) = s%Basis(i)%Un*(rho(mud,i) - rho0(mud,i)) + 0.5_dp*s%Basis(i)%Un*( - rhod(i) + rhod0(i)) &
                        +0.5_dp*s%Basis(i)%Um*(-mzd(i)+mzd0(i))
        hee(nud,nud,i) = s%Basis(i)%Un*(rho(mud,i) - rho0(mud,i)) + 0.5_dp*s%Basis(i)%Un*( - rhod(i) + rhod0(i)) &
                        +0.5_dp*s%Basis(i)%Um*( mzd(i)-mzd0(i))
        hee(mud,nud,i) =-0.5_dp*s%Basis(i)%Um*conjg(mpd(i)-mpd0(i))
        hee(nud,mud,i) =-0.5_dp*s%Basis(i)%Um*     (mpd(i)-mpd0(i))
      end do
    end do
  end subroutine update_Umatrix

  subroutine init_Umatrix(mzd,mzd0,mpd,mpd0,rhod,rhod0,rho,rho0,s)
  !! This subroutine receives the magnetic moment and occupation per orbital
  !! and per site and initializes the Hubbard term of the Hamiltonian
    use mod_kind,   only: dp
    use mod_System, only: System_type
    implicit none
    type(System_type),                       intent(in) :: s
    real(dp),    dimension(s%nAtoms),        intent(in) :: mzd,mzd0,rhod,rhod0
    complex(dp), dimension(s%nAtoms),        intent(in) :: mpd,mpd0
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(in) :: rho,rho0

    call allocate_Umatrix(s%nAtoms,s%nOrb)
    call update_Umatrix(mzd,mzd0,mpd,mpd0,rhod,rhod0,rho,rho0,s)

  end subroutine init_Umatrix
end module mod_Umatrix
