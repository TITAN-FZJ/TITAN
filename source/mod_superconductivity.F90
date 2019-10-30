! This module implements BCS like superconductivity using the Bogoliuvob-de Gennes
! method (check Bogoliubov-de Gennes Method and Its Applications - Jian-Xin Zhu
! Springer Verlag, it is not implemented exactly like there, but the theory holds)
!
! The superconducting gap is given by \Delta = \lambda*singlet_coupling
! lambda is supposed to hold the constant value that multiplies the
! expected value of the the cooper channels (singlet_coupling)
! singlet_coupling holds the expected value of c\up c\down for each orbital
! check Uriel's thesis
!
! lambda is not supposed to change during the execution
! singlet coupling is expected to change at every step until it converges,
! therefore it is defined outside the subroutines so it can be called and
! modified from other modules
!
! IMPORTANT: For the moment some things are hard-coded, for example in the arrays
! lambda, singlet_coupling, it is assumed that only one atom is present, a therefore
! just nine orbitals in total are needed. Sometimes this is made explicit, for
! example when calculating the expected values
module mod_superconductivity
  use mod_f90_kind,       only: double
   implicit none
  logical :: lsuperCond = .false.
  integer :: superCond
  complex(double), dimension(1:9)  :: lambda, singlet_coupling
contains

  subroutine hamiltk_sc(sys,kp,hk_sc)
    use mod_f90_kind,       only: double
    use mod_parameters,     only: output, kpoints, nOrb2
    use mod_system,         only: System, initHamiltkStride
    use mod_constants,      only: cZero,cOne
    use mod_parameters,     only: offset
    implicit none
    type(System),   intent(in)  :: sys
    real(double), intent(in)  :: kp(3)
    complex(double),dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2) :: hk
    complex(double),dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2), intent(out) :: hk_sc
    integer :: j, i

    ! Initialize hk, this will be used to calculate the superconducting
    ! counterpart, namely hk_sc
    hk = cZero

    ! The form of the superconducting hamiltonian depends on a series of decisions,
    ! such as how to choose the basis after the Bogoliuvob-de Gennes transformation
    ! and how do we define this transformation. In our particular case, we choose to
    ! have a basis (u_up, u_down, v_up, v_down), and the BdG transformation we perform
    ! is c_{i\sigma} = sum_n u_{in\sigma}\gamma_n + v_{in\sigma}\gamma^{n*}\gamma_n^\dagger
    ! The final form of the hamiltonian is something like
    ! | H - E_f   Delta     |
    ! | Delta^* -(H - E_f)* |
    ! roughly. Look at this paper 10.1103/RevModPhys.87.1037 , and to Uriel's thesis to
    ! get a better idea of how to construct this operator

    call hamiltk(sys,kp,hk)

    ! Populate the diagonal blocks of the hamiltonian. i.e. electron-electron
    ! and hole-hole interactions
    hk_sc = cZero
    hk_sc(1:sys%nAtoms*nOrb2,1:sys%nAtoms*nOrb2) = hk
    hk_sc(sys%nAtoms*nOrb2+1:sys%nAtoms*nOrb2*2,sys%nAtoms*nOrb2+1:sys%nAtoms*nOrb2*2) = -conjg(hk)
    ! The diagonal terms involve also the Fermi Energy/chemical potential, as we can see below
    ! Check any superconductivity reference for this detail
    do i = 1, sys%nAtoms*nOrb2
        hk_sc(i,i) = hk_sc(i,i) - sys%Ef*cOne
        hk_sc(sys%nAtoms*nOrb2+i,sys%nAtoms*nOrb2+i) = hk_sc(sys%nAtoms*nOrb2+i,sys%nAtoms*nOrb2+i) + sys%Ef*cOne
    end do

    ! Once this is done the remaining part is to populate the non-diagonal blocks
    ! of the hamiltonian. There are several ways to do it.

    ! Later we can add conditional clauses that call or not this functions
    call bcs_s_pairing(sys, singlet_coupling(1),hk_sc)

    call bcs_p_pairing(sys, 0, singlet_coupling(2), hk_sc)
    call bcs_p_pairing(sys, 1, singlet_coupling(3), hk_sc)
    call bcs_p_pairing(sys, 2, singlet_coupling(4), hk_sc)

    call bcs_d_pairing(sys, 0, singlet_coupling(5), hk_sc)
    call bcs_d_pairing(sys, 1, singlet_coupling(6), hk_sc)
    call bcs_d_pairing(sys, 2, singlet_coupling(7), hk_sc)
    call bcs_d_pairing(sys, 3, singlet_coupling(8), hk_sc)
    call bcs_d_pairing(sys, 4, singlet_coupling(9), hk_sc)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Block used to print the hamiltonian
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !Prints to check the shape of the matrix
    ! do i = 1, sys%nAtoms*nOrb2*2
    !   do j = 1, sys%nAtoms*nOrb2*2
    !     write(*,*) real(hk_sc(i,j)), imag(hk_sc(i,j))
    !   end do
    ! end do

  end subroutine hamiltk_sc

  subroutine update_singlet_couplings(couplings)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints, nOrb2
      use mod_system,         only: System, initHamiltkStride
      use mod_constants,      only: cZero,cOne
      use mod_parameters,     only: offset
      implicit none

      complex(double), dimension(1:9)  :: couplings
      integer :: i

      do i = 1,9
          singlet_coupling(i) = lambda(i)*cOne*couplings(i)
      end do

      ! write(*,*) "Couplings: ", singlet_coupling

  end subroutine update_singlet_couplings

  subroutine bcs_s_pairing(sys,delta_s, hk_sc)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints, nOrb2
      use mod_system,         only: System, initHamiltkStride
      use mod_constants,      only: cZero,cOne
      use mod_parameters,     only: offset
      implicit none

      type(System),                              intent(in)  :: sys
      complex(double), intent(in) :: delta_s
      complex(double), dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2), intent(inout) :: hk_sc

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

      hk_sc(1,sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2 + 1) = -delta_s
      hk_sc(sys%nAtoms*nOrb2/2 + 1, sys%nAtoms*nOrb2 + 1 ) = delta_s
      hk_sc(sys%nAtoms*nOrb2 + 1, sys%nAtoms*nOrb2/2 + 1) = conjg(delta_s)
      hk_sc(sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2 + 1,1) = -conjg(delta_s)

  end subroutine bcs_s_pairing

  subroutine bcs_p_pairing(sys, label, delta_p, hk_sc)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints, nOrb2
      use mod_system,         only: System, initHamiltkStride
      use mod_constants,  only: cZero,cOne
      use mod_parameters, only: offset
      implicit none

      ! TODO put a restriction so label can only be in {0,1,2}

      type(System),                              intent(in)  :: sys
      ! "label" is a parameter to pick the specific p orbital (one of the 3) to
      ! couple, it can be 0, 1, or 2
      integer :: label, elecOrbs, indexJump
      complex(double), intent(in) :: delta_p
      complex(double), dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2), intent(inout) :: hk_sc

      elecOrbs = sys%nAtoms*nOrb2 ! Number of up and down electron orbitals

      ! Same idea to populate as in the s-pairing section above.

      indexJump = 1 + label ! local variable, to find the p orbitals respect to s

      hk_sc(1 + indexJump,elecOrbs + elecOrbs/2 + indexJump + 1 ) = -delta_p
      hk_sc(elecOrbs/2 + indexJump + 1, elecOrbs + 1 + indexJump) = delta_p
      hk_sc(elecOrbs + 1 + indexJump, elecOrbs/2 + indexJump + 1) = conjg(delta_p)
      hk_sc(elecOrbs + elecOrbs/2 + indexJump + 1 ,1 + indexJump) = -conjg(delta_p)

  end subroutine bcs_p_pairing

  subroutine bcs_d_pairing(sys, label, delta_d, hk_sc)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints, nOrb2
      use mod_system,         only: System, initHamiltkStride
      use mod_constants,  only: cZero,cOne
      use mod_parameters, only: offset
      implicit none

      ! TODO put a restriction so label can only be in {0,1,2,3,4}

      type(System),                              intent(in)  :: sys
      ! "label" is a parameter to pick the specific p orbital to couple, it can
      ! be in {0,1,2,3,4}
      integer :: label, elecOrbs, indexJump
      complex(double), intent(in) :: delta_d
      complex(double), dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2), intent(inout) :: hk_sc

      elecOrbs = sys%nAtoms*nOrb2 ! Number of up and down electron orbitals

      ! Same idea to populate as in the s-pairing section above.

      indexJump = 4 + label ! local variable, to find the p orbitals respect to s

      hk_sc(1 + indexJump,elecOrbs + elecOrbs/2 + indexJump+ 1) = -delta_d
      hk_sc(elecOrbs/2 + indexJump + 1, elecOrbs + 1 + indexJump) = delta_d
      hk_sc(elecOrbs + 1 + indexJump, elecOrbs/2 + indexJump + 1) = conjg(delta_d)
      hk_sc(elecOrbs + elecOrbs/2 + indexJump+ 1,1 + indexJump) = -conjg(delta_d)

  end subroutine bcs_d_pairing

end module mod_superconductivity
