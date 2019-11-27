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
! lambda, singlet_coupling, it is assumed that only one atom is present, and therefore
! just nine orbitals in total are needed. Sometimes this is made explicit, for
! example when calculating the expected values
module mod_superconductivity
  use mod_f90_kind,       only: double
   implicit none
  logical :: lsuperCond = .false.
  integer :: superCond
  complex(double), dimension(:,:), allocatable  :: singlet_coupling
  integer :: flag = 0

contains

  subroutine allocate_super_variables(nAtoms,nOrbs)
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms
    integer, intent(in) :: nOrbs
    integer :: AllocateStatus

    allocate( singlet_coupling(nOrbs,nAtoms), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocate_super_variables] Not enough memory for: singlet_coupling")

end subroutine allocate_super_variables

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

    call bcs_pairing(sys, singlet_coupling,hk_sc)

    ! call bcs_p_pairing(sys, 0, singlet_coupling(2), hk_sc)
    ! call bcs_p_pairing(sys, 1, singlet_coupling(3), hk_sc)
    ! call bcs_p_pairing(sys, 2, singlet_coupling(4), hk_sc)
    !
    ! call bcs_d_pairing(sys, 0, singlet_coupling(5), hk_sc)
    ! call bcs_d_pairing(sys, 1, singlet_coupling(6), hk_sc)
    ! call bcs_d_pairing(sys, 2, singlet_coupling(7), hk_sc)
    ! call bcs_d_pairing(sys, 3, singlet_coupling(8), hk_sc)
    ! call bcs_d_pairing(sys, 4, singlet_coupling(9), hk_sc)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Block used to print the hamiltonian
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Prints to check the shape of the matrix
    ! do i = 1, sys%nAtoms*nOrb2*2
    !   do j = 1, sys%nAtoms*nOrb2*2
    !     write(*,*) real(hk_sc(i,j)), imag(hk_sc(i,j))
    !   end do
    ! end do

    stop

  end subroutine hamiltk_sc

  subroutine update_singlet_couplings(sys,couplings)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints, nOrb2
      use mod_system,         only: System, initHamiltkStride
      use mod_constants,      only: cZero,cOne
      use mod_parameters,     only: offset, nOrb
      implicit none

      type(System),   intent(in)  :: sys

      complex(double), dimension(nOrb,sys%nAtoms)  :: couplings
      integer :: i,mu

      ! sys%Types(i)%lambda(1:9)
      ! do i=sys%nAtoms

      do i = 1,sys%nAtoms
          do mu = 1,nOrb
              singlet_coupling(mu,i) = sys%Types(i)%lambda(mu)*cOne*couplings(mu,i)
          end do
      end do

      ! write(*,*) "Couplings: ", singlet_coupling

  end subroutine update_singlet_couplings

  subroutine bcs_pairing(sys,delta, hk_sc)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints, nOrb2
      use mod_system,         only: System, initHamiltkStride
      use mod_constants,      only: cZero,cOne
      use mod_parameters,     only: offset, nOrb, isigmamu2n
      implicit none

      type(System),                                intent(in) :: sys
      complex(double), dimension(nOrb,sys%nAtoms), intent(in) :: delta
      complex(double), dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2), intent(inout) :: hk_sc
      integer :: dummy, i, mu

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

      dummy = nOrb2*sys%nAtoms

      do i = 1,sys%nAtoms
          do mu = 1,nOrb
              hk_sc(isigmamu2n(i,1,mu),isigmamu2n(i,2,mu)+dummy) = - delta(mu,i)
              hk_sc(isigmamu2n(i,2,mu),isigmamu2n(i,1,mu)+dummy) = delta(mu,i)
              hk_sc(isigmamu2n(i,2,mu)+dummy,isigmamu2n(i,1,mu)) = - conjg(delta(mu,i))
              hk_sc(isigmamu2n(i,1,mu)+dummy,isigmamu2n(i,2,mu)) = conjg(delta(mu,i))
              write(*,*) i, " ", mu, " ", delta(mu,i)
          end do
      end do



  end subroutine bcs_pairing

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

  ! subroutine green_sc(er,ei,sys,kp,gf)
  subroutine green_sc(sys,kp)
    use mod_f90_kind,   only: double
    use mod_constants,  only: cZero,cOne
    use mod_System,     only: ia, System
    use mod_parameters, only: nOrb2, offset
    implicit none
    integer     :: i,j,d
    real(double) :: er,ei
    real(double), intent(in) :: kp(3)
    type(System), intent(in) :: sys
    complex(double) :: ec
    complex(double),dimension(sys%nAtoms*nOrb2*2.0, sys%nAtoms*nOrb2*2.0) :: gslab,hk
    complex(double),dimension(nOrb2, nOrb2, sys%nAtoms, sys%nAtoms)  :: gf

    ! real(double), intent(in) :: er,ei,kp(3)
    ! type(System), intent(in) :: sys
    ! complex(double) :: ec
    ! complex(double),dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2) :: gslab,hk
    ! complex(double),dimension(nOrb2, nOrb2, sys%nAtoms, sys%nAtoms), intent(out)  :: gf

    ! Dimension of the matrices (they are square)
    d = sys%nAtoms * nOrb2 * 2.0

    if (flag < 1) then
        write(*,*) "d = ", d
        write(*,*) "sys%nAtoms = ", sys%nAtoms
        write(*,*) "nOrb2 = ", nOrb2
    end if

    ec    = cmplx(er,ei,double)

    if (flag < 1) then
        write(*,*) "er = ", er
        write(*,*) "ei = ", ei
    end if

    gslab = cZero
    do i = 1, d
      gslab(i,i) = ec
    end do

    call hamiltk_sc(sys,kp,hk)

    call zaxpy(d*d,-cOne,hk,1,gslab,1)
    ! !gslab(i,j) = gslab(i,j) - hk(i,j)
    call invers(gslab, d)

    if (flag < 1) then
        !Prints to check the shape of the matrix
        do i = 1, sys%nAtoms*nOrb2*2
          do j = 1, sys%nAtoms*nOrb2*2
            write(*,*) real(gslab(i,j)), imag(gslab(i,j))
          end do
        end do
    end if
    !
    ! ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
    ! !dir$ ivdep:loop
    ! do j = 1, sys%nAtoms
    !   !dir$ ivdep:loop
    !   do i = 1, sys%nAtoms
    !     gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    !   end do
    ! end do

    flag = flag + 1

end subroutine green_sc

end module mod_superconductivity
