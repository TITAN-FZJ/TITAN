!   Calculates spin-resolved LDOS and energy-dependence of exchange interactions
subroutine ldos_energy(e,ldosu,ldosd)
  use mod_f90_kind, only: double
  use mod_constants, only: pi
  use mod_parameters, only: lverbose, eta, outputunit_loop
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use mod_progress
  !use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer           :: nthreads,mythread
  integer             :: i,mu,nu,iz
  real(double)        :: kp(3)
  real(double),intent(in) :: e
  real(double),intent(out) :: ldosu(s%nAtoms,nOrb),ldosd(s%nAtoms,nOrb)
  complex(double),dimension(s%nAtoms, s%nAtoms, 2*nOrb, 2*nOrb) :: gf
  complex(double),dimension(s%nAtoms, nOrb) :: gfdiagu,gfdiagd

  ldosu = 0.d0
  ldosd = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,mu,nu,gfdiagu,gfdiagd) &
!$omp& shared(lverbose,s,e,eta,ldosu,ldosd,nthreads,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if(mythread.eq.0) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[ldos_jij_energy] Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:ldosu,ldosd)
  do iz=1,s%nkpt
!$  if((mythread.eq.0)) then
      if(lverbose) call progress_bar(outputunit_loop,"kpoints",iz,s%nkpt)
!$  end if
    kp = s%kbz(:,iz)

    ! Green function on energy E + ieta, and wave vector kp
    call green(e,eta,kp,gf)

    ! Density of states
    do mu=1,nOrb; do i=1,s%nAtoms
      nu=mu+nOrb
      gfdiagu(i,mu) = - aimag(gf(i,i,mu,mu))*s%wkbz(iz)
      gfdiagd(i,mu) = - aimag(gf(i,i,nu,nu))*s%wkbz(iz)
    end do ; end do

    ldosu = ldosu + gfdiagu
    ldosd = ldosd + gfdiagd

  end do
!$omp end do
!$omp end parallel

  ldosu  = ldosu/pi
  ldosd  = ldosd/pi

  return
end subroutine ldos_energy
