!   Calculates spin-resolved LDOS and energy-dependence of exchange interactions
subroutine ldos_energy(e,ldosu,ldosd)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_system, only: kbz, wkbz, nkpt
  use mod_progress
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer           :: nthreads,mythread
  integer             :: i,mu,nu,iz
  real(double)        :: kp(3)
  real(double),intent(in)     :: e
  real(double),intent(out)    :: ldosu(Npl,9),ldosd(Npl,9)
  complex(double),dimension(Npl,Npl,18,18)    :: gf
  complex(double),dimension(Npl,9)            :: gfdiagu,gfdiagd

  ldosu = 0.d0
  ldosd = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,mu,nu,gfdiagu,gfdiagd) &
!$omp& shared(lverbose,kbz,nkpt,wkbz,e,eta,Npl,ldosu,ldosd,nthreads,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if(mythread.eq.0) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[ldos_jij_energy] Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:ldosu,ldosd)
  kpoints: do iz=1,nkpt
!$  if((mythread.eq.0)) then
      if(lverbose) call progress_bar(outputunit_loop,"kpoints",iz,nkpt)
!$  end if
    kp = kbz(:,iz)

    ! Green function on energy E + ieta, and wave vector kp
    call green(e,eta,kp,gf)

    ! Density of states
    do mu=1,9; do i=1,Npl
      nu=mu+9
      gfdiagu(i,mu) = - aimag(gf(i,i,mu,mu))*wkbz(iz)
      gfdiagd(i,mu) = - aimag(gf(i,i,nu,nu))*wkbz(iz)
    end do ; end do

    ldosu = ldosu + gfdiagu
    ldosd = ldosd + gfdiagd

  end do kpoints
!$omp end do
!$omp end parallel

  ldosu  = ldosu/pi
  ldosd  = ldosd/pi

  return
end subroutine ldos_energy
