! This module implements BCS like superconductivity using the Bogoliuvob-de Gennes
! method (check Bogoliubov-de Gennes Method and Its Applications - Jian-Xin Zhu
! Springer Verlag, it is not implemented exactly like there, but the theory holds)
!
! The superconducting gap is given by \Delta = \lambda*singlet_coupling
! lambda is supposed to hold the constant value that multiplies the
! expected value of the the cooper channels (singlet_coupling)
! singlet_coupling holds the expected value of c\up c\down for each orbital
! singlet_coupling is expected to be real within BCS
! check Uriel's thesis
!
! lambda is not supposed to change during the execution
! singlet coupling is expected to change at every step until it converges,
! therefore it is defined outside the subroutines so it can be called and
! modified from other modules

module mod_superconductivity
  use mod_f90_kind,       only: double
  implicit none
  logical :: lsuperCond = .false.
  integer :: superCond
  real(double), dimension(:,:), allocatable  :: singlet_coupling
  integer :: flag = 0

contains

  subroutine allocate_supercond_variables(nAtoms,nOrbs)
    use mod_mpi_pars,  only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms
    integer, intent(in) :: nOrbs
    integer :: AllocateStatus

    allocate( singlet_coupling(nOrbs,nAtoms), stat = AllocateStatus)
    singlet_coupling = 0.d0
    if(AllocateStatus /= 0) call abortProgram("[allocate_supercond_variables] Not enough memory for: singlet_coupling")

  end subroutine allocate_supercond_variables

  subroutine deallocate_supercond_variables()
    implicit none

    if(allocated(singlet_coupling)) deallocate(singlet_coupling)
  end subroutine deallocate_supercond_variables

  subroutine hamiltk_sc(sys,kp,hk_sc)
    use mod_f90_kind,       only: double
    use mod_parameters,     only: dimH
    use mod_system,         only: System, initHamiltkStride
    use mod_constants,      only: cZero,cOne
    implicit none
    type(System), intent(in)  :: sys
    real(double), intent(in)  :: kp(3)
    complex(double),dimension(2*dimH,2*dimH), intent(inout) :: hk_sc
    complex(double),dimension(  dimH,  dimH) :: hk, dummy, hk2
    integer      :: i, j
    real(double) :: mkp(3)

    ! The form of the superconducting hamiltonian depends on a series of decisions,
    ! such as how to choose the basis after the Bogoliuvob-de Gennes transformation
    ! and how do we define this transformation. In our particular case, we choose to
    ! have a basis (u_up, u_down, v_up, v_down), and the BdG transformation we perform
    ! is c_{i\sigma} = sum_n u_{in\sigma}\gamma_n + v_{in\sigma}\gamma^{n*}\gamma_n^\dagger
    ! The final form of the hamiltonian is something like
    ! | H(k) - E_f        Delta     |
    ! |   Delta^*   -(H(-k) - E_f)* |
    ! roughly. Look at this paper 10.1103/RevModPhys.87.1037 , and to Uriel's thesis to
    ! get a better idea of how to construct this operator

    !!This part is neccessary since hamiltk returns a hamiltonian, hermitian
    !!up to 1e-16 precision. This interfered with some expectation value calculations
    !!causing an artificial anisotropy
    !!With this I'm just making sure that diagonal blocks of the superconducting
    !!Hamiltonian are hermitian, and therefore the total hamiltonian is hermitian.

    ! call hamiltk(sys,kp,hk)

    call hamiltk(sys,kp,dummy)

    do concurrent (i = 1:dimH , j = 1:dimH)
      hk(i,j) = (dummy(i,j) + conjg(dummy(j,i)))/2.d0
    end do

    dummy = cZero

    mkp = -kp
    call hamiltk(sys,mkp,dummy)

    do concurrent (i = 1:dimH , j = 1:dimH)
      hk2(i,j) = (dummy(i,j) + conjg(dummy(j,i)))/2.d0
    end do

    ! Populate the diagonal blocks of the hamiltonian. i.e. electron-electron
    ! and hole-hole interactions
    hk_sc = cZero
    hk_sc(     1:  dimH,     1:  dimH) = hk
    hk_sc(dimH+1:2*dimH,dimH+1:2*dimH) = -conjg(hk2)
    ! The diagonal terms involve also the Fermi Energy/chemical potential, as we can see below
    ! Check any superconductivity reference for this detail
    do concurrent (i = 1:dimH)
      hk_sc(     i,     i) = hk_sc(     i,     i) - sys%Ef*cOne
      hk_sc(dimH+i,dimH+i) = hk_sc(dimH+i,dimH+i) + sys%Ef*cOne
    end do

    ! Once this is done the remaining part is to populate the non-diagonal blocks
    ! of the hamiltonian. There are several ways to do it.

    call bcs_pairing(sys, singlet_coupling, hk_sc)

  end subroutine hamiltk_sc

  ! TODO: Either remove this function or write to a file
  subroutine print_hamilt(sys,hk)
    use mod_f90_kind,       only: double
    use mod_parameters,     only: nOrb2
    use mod_system,         only: System, initHamiltkStride
    use mod_parameters,     only: dimH
    implicit none

    type(System),                              intent(in) :: sys
    complex(double), dimension(2*dimH,2*dimH), intent(in) :: hk
    integer :: i,j

    do j = 1, sys%nAtoms*nOrb2*2
      do i = 1, sys%nAtoms*nOrb2*2
        write(*,*) real(hk(i,j)), imag(hk(i,j))
      end do
    end do

  end subroutine print_hamilt

  ! TODO: Either remove this function or write to a file
  subroutine print_hamilt_entry(hk,i,j)
    use mod_f90_kind,       only: double
    use mod_system,         only: System, initHamiltkStride
    use mod_parameters,     only: dimH
    implicit none
    complex(double), dimension(2*dimH,2*dimH), intent(in) :: hk
    integer, intent(in):: i,j

    write(*,*) real(hk(i,j)), imag(hk(i,j))

  end subroutine print_hamilt_entry

  subroutine update_singlet_couplings(sys,couplings)
    use mod_f90_kind,       only: double
    use mod_system,         only: System, initHamiltkStride
    use mod_parameters,     only: nOrb
    implicit none

    type(System),   intent(in)  :: sys

    real(double), dimension(nOrb,sys%nAtoms)  :: couplings
    integer :: i,mu

    do concurrent (i = 1:sys%nAtoms, mu = 1:nOrb)
      singlet_coupling(mu,i) = couplings(mu,i)
    end do

  end subroutine update_singlet_couplings

  subroutine bcs_pairing(sys,delta, hk_sc)
    use mod_f90_kind,       only: double
    use mod_system,         only: System, initHamiltkStride
    use mod_parameters,     only: nOrb, isigmamu2n, dimH
    implicit none

    type(System),                                 intent(in) :: sys
    real(double),    dimension(nOrb,sys%nAtoms),  intent(in) :: delta
    complex(double), dimension(2*dimH,2*dimH), intent(inout) :: hk_sc
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

    do concurrent (i = 1:sys%nAtoms, mu = 1:nOrb)
      hk_sc(isigmamu2n(i,1,mu)     ,isigmamu2n(i,2,mu)+dimH) = - cmplx(delta(mu,i),0.d0,double)
      hk_sc(isigmamu2n(i,2,mu)     ,isigmamu2n(i,1,mu)+dimH) =   cmplx(delta(mu,i),0.d0,double)
      hk_sc(isigmamu2n(i,2,mu)+dimH,isigmamu2n(i,1,mu)     ) = - cmplx(delta(mu,i),0.d0,double)
      hk_sc(isigmamu2n(i,1,mu)+dimH,isigmamu2n(i,2,mu)     ) =   cmplx(delta(mu,i),0.d0,double)
      ! write(*,*) i, " ", mu, " ", delta(mu,i)
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

  subroutine green_sc(er,ei,sys,kp,gf)
    use mod_f90_kind,   only: double
    use mod_constants,  only: cZero,cOne
    use mod_System,     only: ia_sc, System
    use mod_parameters, only: nOrb2, dimH
    implicit none
    integer     :: i,j,d

    real(double), intent(in) :: er,ei,kp(3)
    type(System), intent(in) :: sys
    complex(double),dimension(nOrb2*superCond, nOrb2*superCond, sys%nAtoms, sys%nAtoms), intent(out)  :: gf

    complex(double) :: ec
    complex(double),dimension(dimH*superCond, dimH*superCond) :: gslab,hk

    ! Dimension of the matrices (they are square)
    d = dimH * superCond

    ec    = cmplx(er,ei,double)

    gslab = cZero
    do concurrent (i = 1:d)
      gslab(i,i) = ec
    end do

    call hamiltk_sc(sys,kp,hk)

    !zaxpy performs gslab = gslab + (-cOne)*hk
    !d*d si the dimension of the array
    !both 1's are for the strides to use on each array
    !for example we can add each second or third entry of the array
    !on this case we just use every component
    call zaxpy(d*d,-cOne,hk,1,gslab,1)

    call invers(gslab, d)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Block to print
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if(flag == 1) then
    !     do i =1,dimH*superCond
    !         do j = 1, dimH*superCond
    !             write(*,*) real(gslab(i,j)), aimag(gslab(i,j))
    !         end do
    !     end do
    !     stop
    ! end if
    !
    ! flag = flag + 1

    ! Put the slab Green's function [A(nAtoms*36,nAtoms*36)] in the A(i,j,mu,nu) form
    do concurrent (j = 1:sys%nAtoms, i = 1:sys%nAtoms)
      gf(      1:  nOrb2,      1:  nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(1,j):ia_sc(2,j))
      gf(      1:  nOrb2,nOrb2+1:2*nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(3,j):ia_sc(4,j))
      gf(nOrb2+1:2*nOrb2,      1:  nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(1,j):ia_sc(2,j))
      gf(nOrb2+1:2*nOrb2,nOrb2+1:2*nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(3,j):ia_sc(4,j))
    end do

  end subroutine green_sc

end module mod_superconductivity
