module mod_superconductivity
  implicit none

contains

  subroutine hamilt_sc(s)
  ! subroutine hamilt_sc(s,rho,mp,mx,my,mz,sys,kp,hk)
    use mod_f90_kind,       only: double
    use mod_BrillouinZone,  only: realBZ
  !   use mod_parameters,     only: output
    use mod_system,         only: System
    use TightBinding,       only: nOrb,nOrb2
  !   use ElectricField,      only: EshiftBZ,ElectricFieldVector
    implicit none
    type(System),                              intent(in)  :: s
  !   real(double),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
  !   complex(double), dimension(nOrb,s%nAtoms), intent(out) :: mp
    ! real(double), intent(in) :: kp(3)
    real(double) :: kp(3)

  !   type(System), intent(in) :: sys
    complex(double),dimension(s%nAtoms*nOrb2, s%nAtoms*nOrb2) :: hk
    !
    ! integer                      :: iz, info !, mu,i
    ! integer                      :: lwork,dimH
    ! real(double)                 :: weight, kp(3)
    ! real(double),    dimension(nOrb,s%nAtoms)   :: expec_0, expec_z
    ! complex(double), dimension(nOrb,s%nAtoms)   :: expec_p
    ! real(double),    dimension(:),  allocatable :: rwork(:), eval(:)
    ! complex(double),                allocatable :: work(:), hk(:,:), dummyhk(:,:)
    !
    ! dimH  = (s%nAtoms)*nOrb2
    ! lwork = 21*dimH
    !
    ! allocate( hk(dimH*2,dimH*2),dummyhk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork) )
    ! integer                      :: lwork,dimH

    write(*,*) kp(1:3)

    kp(1) = 3.0
    kp(2) = 3.0
    kp(3) = 3.0

    write(*,*) kp(1:3)
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
    ! call hamiltk(s,kp,hk)
    !
    ! hk(1:dimH,1:dimH) = dummyhk
    ! hk(dimH+1:dimH*2,dimH+1:dimH*2) = -dummyhk
    !
    ! ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
    ! !call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)
    !
    ! deallocate(hk,rwork,eval,work)
    write(*,*) "You are an amazing man Uriel"

    write(*,*) "Also very handsome I may say"

  end subroutine hamilt_sc

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
    ! d = sys%nAtoms * nOrb2
    !
    ! ec    = cmplx(er,ei,double)
    !
    ! gslab = cZero
    ! do i = 1, d
    !   gslab(i,i) = ec
    ! end do
    !
    ! call hamiltk(sys,kp,hk)
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
