!   Calculates spin-resolved LDOS and energy-dependence of exchange interactions
subroutine ldos_energy(e,ldosu,ldosd)
  use mod_f90_kind,      only: double
  use mod_constants,     only: pi
  use mod_parameters,    only: eta
  use mod_system,        only: s => sys
  use mod_BrillouinZone, only: realBZ
  use TightBinding,      only: nOrb,nOrb2
  use mod_mpi_pars
  implicit none
  real(double), intent(in) :: e
  real(double), dimension(s%nAtoms, nOrb), intent(out) :: ldosu, ldosd
  complex(double), dimension(nOrb2, nOrb2, s%nAtoms, s%nAtoms) :: gf
  complex(double), dimension(s%nAtoms, nOrb) :: gfdiagu,gfdiagd
  real(double), dimension(3) :: kp
  real(double) :: weight
  integer :: i,mu,nu
  integer*8 :: iz

  ldosu = 0.d0
  ldosd = 0.d0

!$omp parallel default(none) &
!$omp& private(iz,kp,weight,gf,i,mu,nu,gfdiagu,gfdiagd) &
!$omp& shared(s,realBZ,e,eta,ldosu,ldosd)

!$omp do reduction(+:ldosu,ldosd)
  do iz = 1,realBZ%workload
    kp = realBZ%kp(1:3,iz)
    weight = realBZ%w(iz)
    ! Green function on energy E + ieta, and wave vector kp
    call green(e,eta,s,kp,gf)

    ! Density of states
    do mu=1,nOrb
      do i=1,s%nAtoms
         nu = mu + nOrb
         gfdiagu(i,mu) = - aimag(gf(mu,mu,i,i)) * weight
         gfdiagd(i,mu) = - aimag(gf(nu,nu,i,i)) * weight
      end do
   end do

    ldosu = ldosu + gfdiagu
    ldosd = ldosd + gfdiagd

  end do
!$omp end do
!$omp end parallel

  ldosu  = ldosu/pi
  ldosd  = ldosd/pi

    if(rFreq(1) == 0) then
       call MPI_Reduce(MPI_IN_PLACE, ldosu , s%nAtoms*nOrb, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
       call MPI_Reduce(MPI_IN_PLACE, ldosd , s%nAtoms*nOrb, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    else
       call MPI_Reduce(ldosu , ldosu , s%nAtoms*nOrb, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
       call MPI_Reduce(ldosd , ldosd , s%nAtoms*nOrb, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    end if

end subroutine ldos_energy
