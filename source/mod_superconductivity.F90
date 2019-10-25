module mod_superconductivity
  use mod_f90_kind,       only: double
   implicit none
  logical :: lsuperCond = .false.
  integer :: superCond
  complex(double), dimension(1:9)  :: lambda, singlet_coupling
contains

  subroutine hamiltk_sc(sys,kp,hk_sc)
    use mod_f90_kind,       only: double
    use mod_parameters,     only: output, kpoints
    use mod_system,         only: System, initHamiltkStride
    use TightBinding,       only: nOrb2, initTightBinding
    use mod_constants,      only: cZero,cOne
    use mod_parameters,     only: offset
    implicit none
    type(System),   intent(in)  :: sys
    real(double), intent(in)  :: kp(3)
    complex(double),dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2) :: hk
    complex(double),dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2), intent(out) :: hk_sc
    integer :: j, i
    ! Coupling constant fot the s-orbital pairing
    complex(double) :: delta_s, delta_p, delta_d

    ! singlet_coupling = singlet_coupling*lambda

    ! write(*,*) singlet_coupling(:)

    ! delta_s = lambda(1)*cOne
    ! delta_p = lambda(2)*cOne !+ conjg(lambda(2)*cOne)
    ! delta_d = lambda(5)*cOne

    ! singlet_coupling(:) = lambda(:)*cOne

    ! Initialize hk, this will be used to calculate the superconducting
    ! counterpart, namely hk_sc
    hk = cZero

    call hamiltk(sys,kp,hk)

    ! Populate the diagonal blocks of the hamiltonian. i.e. electron-electron
    ! and hole-hole interactions
    hk_sc = cZero
    hk_sc(1:sys%nAtoms*nOrb2,1:sys%nAtoms*nOrb2) = hk
    hk_sc(sys%nAtoms*nOrb2+1:sys%nAtoms*nOrb2*2,sys%nAtoms*nOrb2+1:sys%nAtoms*nOrb2*2) = -conjg(hk)

    do i = 1, sys%nAtoms*nOrb2
        hk_sc(i,i) = hk_sc(i,i) - sys%Ef*cOne
        hk_sc(sys%nAtoms*nOrb2+i,sys%nAtoms*nOrb2+i) = hk_sc(sys%nAtoms*nOrb2+i,sys%nAtoms*nOrb2+i) + sys%Ef*cOne
    end do

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

    ! ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
    ! !call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)

  end subroutine hamiltk_sc

  subroutine update_singlet_couplings(couplings)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints
      use mod_system,         only: System, initHamiltkStride
      use TightBinding,       only: nOrb2, initTightBinding
      use mod_constants,      only: cZero,cOne
      use mod_parameters,     only: offset
      implicit none

      complex(double), dimension(1:9)  :: couplings
      integer :: i

      do i = 1,9
          singlet_coupling(i) = lambda(i)*cOne*couplings(i)
      end do

  end subroutine update_singlet_couplings

  subroutine bcs_s_pairing(sys,delta_s, hk_sc)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints
      use mod_system,         only: System, initHamiltkStride
      use TightBinding,       only: nOrb2, initTightBinding
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

      ! write(*,*) sys%nAtoms*nOrb2 + 1, sys%nAtoms*nOrb2/2 +1
      ! write(*,*) 1,sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2 + 1
      ! hk_sc(1,sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2) = -delta_s
      ! hk_sc(sys%nAtoms*nOrb2/2, sys%nAtoms*nOrb2 + 1 ) = delta_s
      ! hk_sc(sys%nAtoms*nOrb2 + 1, sys%nAtoms*nOrb2/2) = conjg(delta_s)
      ! hk_sc(sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2,1) = -conjg(delta_s)

      hk_sc(1,sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2 + 1) = -delta_s
      hk_sc(sys%nAtoms*nOrb2/2 + 1, sys%nAtoms*nOrb2 + 1 ) = delta_s
      hk_sc(sys%nAtoms*nOrb2 + 1, sys%nAtoms*nOrb2/2 + 1) = conjg(delta_s)
      hk_sc(sys%nAtoms*nOrb2 + sys%nAtoms*nOrb2/2 + 1,1) = -conjg(delta_s)

  end subroutine bcs_s_pairing

  subroutine bcs_p_pairing(sys, label, delta_p, hk_sc)
      use mod_f90_kind,       only: double
      use mod_parameters,     only: output, kpoints
      use mod_system,         only: System, initHamiltkStride
      use TightBinding,       only: nOrb2, initTightBinding
      use mod_constants,  only: cZero,cOne
      use mod_parameters, only: offset
      implicit none

      ! TODO put a restriction so label can only be in {0,1,2}

      type(System),                              intent(in)  :: sys
      ! "label" is a parameter to pick the specific p orbital to couple, it can
      ! be 0, 1, or 2
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
      use mod_parameters,     only: output, kpoints
      use mod_system,         only: System, initHamiltkStride
      use TightBinding,       only: nOrb2, initTightBinding
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
  !   use mod_f90_kind,   only: double
  !   use mod_constants,  only: cZero,cOne
  !   use mod_System,     only: ia, System
  !   use mod_parameters, only: offset
  !   use TightBinding,   only: nOrb2
  !   implicit none
  !   integer     :: i,j,d
  !   real(double), intent(in) :: er,ei,kp(3)
  !   type(System), intent(in) :: sys
  !   complex(double) :: ec
  !   complex(double),dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2) :: gslab,hk
  !   complex(double),dimension(nOrb2, nOrb2, sys%nAtoms, sys%nAtoms), intent(out)  :: gf
  !   !
  !   ! d = sys%nAtoms * nOrb2 *2
  !   !
  !   ! ec    = cmplx(er,ei,double)
  !   !
  !   ! gslab = cZero
  !   ! do i = 1, d
  !   !   gslab(i,i) = ec
  !   ! end do
  !   !
  !   ! call hamiltk_sc(sys,kp,hk)
  !   !
  !   ! call zaxpy(d*d,-cOne,hk,1,gslab,1)
  !   ! !gslab(i,j) = gslab(i,j) - hk(i,j)
  !   ! call invers(gslab, d)
  !   !
  !   ! ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
  !   ! !dir$ ivdep:loop
  !   ! do j = 1, sys%nAtoms
  !   !   !dir$ ivdep:loop
  !   !   do i = 1, sys%nAtoms
  !   !     gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
  !   !   end do
  !   ! end do
  !
  ! end subroutine green_sc

end module mod_superconductivity
