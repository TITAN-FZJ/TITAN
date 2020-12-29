! This module implements BCS like superconductivity using the Bogoliuvob-de Gennes
! method (check Bogoliubov-de Gennes Method and Its Applications - Jian-Xin Zhu
! Springer Verlag, it is not implemented exactly like there, but the theory holds)
!
! The superconducting gap is given by \Delta = \lambda*delta_sc
! lambda is supposed to hold the constant value that multiplies the
! expected value of the the cooper channels (delta_sc)
! delta_sc holds the expected value of c\up c\down for each orbital
! delta_sc is expected to be real within BCS
! check Uriel's thesis
!
! lambda is not supposed to change during the execution
! singlet coupling is expected to change at every step until it converges,
! therefore it is defined outside the subroutines so it can be called and
! modified from other modules

module mod_superconductivity
  use mod_kind, only: dp
  implicit none
  logical :: lsuperCond = .false.
  integer :: superCond
  real(dp), dimension(:,:), allocatable  :: delta_sc
#ifdef _GPU
  real(dp), dimension(:,:), allocatable,device  :: delta_sc_d
#endif
  integer :: flag = 0

contains

  subroutine allocate_supercond_variables(nAtoms,nOrbs)
    use mod_mpi_pars,  only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms
    integer, intent(in) :: nOrbs

    allocate( delta_sc(nOrbs,nAtoms))
    delta_sc = 0._dp
#ifdef _GPU
    allocate( delta_sc_d(nOrbs,nAtoms))
#endif

  end subroutine allocate_supercond_variables

  subroutine deallocate_supercond_variables()
    implicit none

    if(allocated(delta_sc)) deallocate(delta_sc)
#ifdef _GPU
    if(allocated(delta_sc_d)) deallocate(delta_sc_d)
#endif
  end subroutine deallocate_supercond_variables


  ! TODO: Either remove this function or write to a file
  subroutine print_hamilt(sys,hk)
    use mod_kind,       only: dp
    use mod_parameters, only: nOrb2
    use mod_system,     only: System_type, initHamiltkStride
    use mod_parameters, only: dimH
    implicit none

    type(System_type),                              intent(in) :: sys
    complex(dp), dimension(2*dimH,2*dimH), intent(in) :: hk
    integer :: i,j

    do j = 1, sys%nAtoms*nOrb2*2
      do i = 1, sys%nAtoms*nOrb2*2
        write(*,*) real(hk(i,j)), imag(hk(i,j))
      end do
    end do

  end subroutine print_hamilt


  ! TODO: Either remove this function or write to a file
  subroutine print_hamilt_entry(hk,i,j)
    use mod_kind,       only: dp
    use mod_system,     only: System_type, initHamiltkStride
    use mod_parameters, only: dimH
    implicit none
    complex(dp), dimension(2*dimH,2*dimH), intent(in) :: hk
    integer,                               intent(in):: i,j

    write(*,*) real(hk(i,j)), imag(hk(i,j))

  end subroutine print_hamilt_entry

  subroutine update_delta_sc(sys,couplings)
    use mod_kind,       only: dp
    use mod_system,     only: System_type, initHamiltkStride
    use mod_parameters, only: nOrb
    implicit none
    type(System_type),        intent(in)  :: sys
    real(dp), dimension(nOrb,sys%nAtoms)  :: couplings

    delta_sc = couplings
#ifdef _GPU
    delta_sc_d = delta_sc
#endif

  end subroutine update_delta_sc

  subroutine bcs_pairing(sys,delta,hk_sc)
    use mod_kind,       only: dp
    use mod_system,     only: System_type
    use mod_parameters, only: nOrb, isigmamu2n, dimH
    implicit none

    type(System_type),                       intent(in)    :: sys
    real(dp),    dimension(nOrb,sys%nAtoms), intent(in)    :: delta
    complex(dp), dimension(2*dimH,2*dimH),   intent(inout) :: hk_sc
    integer :: i, mu

    ! Populate the entries for the singlet pairing of the s-orbitals
    ! Assuming that the order to populate the hamiltonian hk was
    ! First all up spins, then all down, from s to d
    ! This is, s^ px^ py^ ... d^ s* px* ... d*, where ^ (*) means spin up (down)

    ! The row and column of the s^ electron is 1
    ! The row and column of the s* electron is sys%nAtoms*nOrb2/2
    ! The row and column of the h^ hole is sys%nAtoms*nOrb2 + 1
    ! The row and column of the h* hole is sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2

    ! s^ and h* couple with -delta_s
    ! s* and h^ couple with delta_s
    ! h^ and s* couple with delta_s*
    ! h* and s^ couple with -delta_s*

    do i = 1,sys%nAtoms
      do mu = 1,nOrb
        hk_sc(isigmamu2n(i,1,mu)     ,isigmamu2n(i,2,mu)+dimH) = - cmplx(delta(mu,i),0._dp,dp)
        hk_sc(isigmamu2n(i,2,mu)     ,isigmamu2n(i,1,mu)+dimH) =   cmplx(delta(mu,i),0._dp,dp)
        hk_sc(isigmamu2n(i,2,mu)+dimH,isigmamu2n(i,1,mu)     ) = - cmplx(delta(mu,i),0._dp,dp)
        hk_sc(isigmamu2n(i,1,mu)+dimH,isigmamu2n(i,2,mu)     ) =   cmplx(delta(mu,i),0._dp,dp)
      end do
    end do


    ! do i = 1,sys%nAtoms
    !   do mu = 1,nOrb
    !     hk_sc(isigmamu2n(i,1,mu)     ,isigmamu2n(i,2,mu)+dimH) = - delta(mu,i)
    !     hk_sc(isigmamu2n(i,2,mu)     ,isigmamu2n(i,1,mu)+dimH) =   delta(mu,i)
    !     hk_sc(isigmamu2n(i,2,mu)+dimH,isigmamu2n(i,1,mu)     ) = - conjg(delta(mu,i))
    !     hk_sc(isigmamu2n(i,1,mu)+dimH,isigmamu2n(i,2,mu)     ) =   conjg(delta(mu,i))
    !     ! write(*,*) i, " ", mu, " ", delta(mu,i)
    !   end do
    ! end do

  end subroutine bcs_pairing


end module mod_superconductivity