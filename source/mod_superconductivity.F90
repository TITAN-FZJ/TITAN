module mod_superconductivity
  implicit none
  logical :: lsuperCond = .false.
  integer :: superCond
contains

  subroutine hamiltk_sc(sys,kp,hk_sc)
  ! subroutine hamilt_sc(sys,rho,mp,mx,my,mz,sys,kp,hk)
    use mod_f90_kind,       only: double
    use mod_BrillouinZone,  only: realBZ
    use mod_parameters,     only: output, kpoints
    use mod_system,         only: System, initHamiltkStride, ia
    use TightBinding,       only: nOrb,nOrb2, initTightBinding
    use ElectricField,      only: EshiftBZ,ElectricFieldVector
    use mod_constants,  only: cZero,cOne
    use mod_parameters, only: offset
    implicit none
    type(System),                              intent(in)  :: sys
  !   real(double),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
  !   complex(double), dimension(nOrb,s%nAtoms), intent(out) :: mp
    real(double), intent(in)  :: kp(3)
    complex(double),dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2) :: hk
    complex(double),dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2), intent(out) :: hk_sc
    integer :: j, i

    hk = cZero

    !kp = EshiftBZ*ElectricFieldVector

    ! write(*,*) kp
    ! kp = realBZ%kp(1:3,1)

    ! write(*,*) kp
    ! allocate( hk(dimH*2,dimH*2))
    !
    ! !$omp parallel default(none) &
    ! !$omp& firstprivate(lwork,dimH) &
    ! !$omp& private(iz,kp,weight,hk,eval,work,rwork,info,expec_0, expec_p, expec_z) &
    ! !$omp& shared(s,output,nOrb2,realBZ,rho,mp,mz,EshiftBZ,ElectricFieldVector)
    !
    ! rho = 0.d0
    ! mp  = 0.d0
    ! mz  = 0.d0
    !
    ! call initHamiltkStride(sys%nAtoms, nOrb)
    call hamiltk(sys,kp,hk)

    hk_sc = cZero
    hk_sc(1:sys%nAtoms*nOrb2,1:sys%nAtoms*nOrb2) = hk
    hk_sc(sys%nAtoms*nOrb2+1:sys%nAtoms*nOrb2*2,sys%nAtoms*nOrb2+1:sys%nAtoms*nOrb2*2) = -conjg(hk)

    ! do i = 1, sys%nAtoms*nOrb2*2
    !   do j = 1, sys%nAtoms*nOrb2*2
    !     write(*,*) real(hk_sc(i,j)), imag(hk_sc(i,j))
    !   end do
    ! end do

    ! ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
    ! !call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)
    !
    ! deallocate(rwork,eval,work)

    ! write(*,*) "You are an amazing man Uriel"
    !
    ! write(*,*) "Also very handsome I may say"

  end subroutine hamiltk_sc

  subroutine green_sc(er,ei,sys,kp,gf)
    use mod_f90_kind,   only: double
    use mod_constants,  only: cZero,cOne
    use mod_System,     only: ia, System
    use mod_parameters, only: offset
    use TightBinding,   only: nOrb2
    implicit none
    integer     :: i,j,d
    real(double), intent(in) :: er,ei,kp(3)
    type(System), intent(in) :: sys
    complex(double) :: ec
    complex(double),dimension(sys%nAtoms*nOrb2*2, sys%nAtoms*nOrb2*2) :: gslab,hk
    complex(double),dimension(nOrb2, nOrb2, sys%nAtoms, sys%nAtoms), intent(out)  :: gf
    !
    ! d = sys%nAtoms * nOrb2 *2
    !
    ! ec    = cmplx(er,ei,double)
    !
    ! gslab = cZero
    ! do i = 1, d
    !   gslab(i,i) = ec
    ! end do
    !
    ! call hamiltk_sc(sys,kp,hk)
    !
    ! call zaxpy(d*d,-cOne,hk,1,gslab,1)
    ! !gslab(i,j) = gslab(i,j) - hk(i,j)
    ! call invers(gslab, d)
    !
    ! ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
    ! !dir$ ivdep:loop
    ! do j = 1, sys%nAtoms
    !   !dir$ ivdep:loop
    !   do i = 1, sys%nAtoms
    !     gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    !   end do
    ! end do

  end subroutine green_sc

end module mod_superconductivity
