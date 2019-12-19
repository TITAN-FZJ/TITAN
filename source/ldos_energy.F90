!   Calculates spin-resolved LDOS and energy-dependence of exchange interactions
subroutine ldos_energy(e,ldosu,ldosd)
  use mod_f90_kind,      only: double
  use mod_constants,     only: pi
  use mod_parameters,    only: nOrb, nOrb2, eta
  use mod_system,        only: s => sys
  use mod_BrillouinZone, only: realBZ
  use mod_mpi_pars
  use mod_superconductivity, only: lsupercond, green_sc, superCond
  implicit none
  real(double), intent(in) :: e
  real(double), dimension(s%nAtoms, nOrb*superCond), intent(out) :: ldosu, ldosd
  complex(double), dimension(nOrb2*superCond, nOrb2*superCond, s%nAtoms, s%nAtoms) :: gf
  complex(double), dimension(s%nAtoms, nOrb*superCond) :: gfdiagu,gfdiagd
  real(double), dimension(3) :: kp
  real(double) :: weight
  integer :: i,mu,nu
  integer*8 :: iz

  ldosu = 0.d0
  ldosd = 0.d0

!$omp parallel default(none) &
!$omp& private(iz,kp,weight,gf,i,mu,nu,gfdiagu,gfdiagd) &
!$omp& shared(s,realBZ,e,nOrb,eta,ldosu,ldosd,lsupercond,nOrb2)

!$omp do reduction(+:ldosu,ldosd)
  do iz = 1,realBZ%workload
    kp = realBZ%kp(1:3,iz)
    weight = realBZ%w(iz)
    ! Green function on energy E + ieta, and wave vector kp
    if(lsupercond == .true.) then
        call green_sc(e,eta,s,kp,gf)
    else
        call green(e,eta,s,kp,gf)
    end if

    ! Density of states
    do mu=1,nOrb
      do i=1,s%nAtoms
         nu = mu + nOrb
         gfdiagu(i,mu) = - aimag(gf(mu,mu,i,i)) * weight
         gfdiagd(i,mu) = - aimag(gf(nu,nu,i,i)) * weight
         if(lsupercond) then
             gfdiagu(i,mu+nOrb) = - aimag(gf(mu+nOrb2,mu+nOrb2,i,i)) * weight
             gfdiagd(i,mu+nOrb) = - aimag(gf(nu+nOrb2,nu+nOrb2,i,i)) * weight
         end if
      end do
   end do

    ldosu = ldosu + gfdiagu
    ldosd = ldosd + gfdiagd

  end do
!$omp end do
!$omp end parallel

  ldosu  = ldosu/pi
  ldosd  = ldosd/pi

  call MPI_Allreduce(MPI_IN_PLACE, ldosu , s%nAtoms*nOrb*superCond, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1), ierr)
  call MPI_Allreduce(MPI_IN_PLACE, ldosd , s%nAtoms*nOrb*superCond, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1), ierr)


end subroutine ldos_energy
