!   Calculates spin-resolved LDOS and energy-dependence of exchange interactions
subroutine ldos_jij_energy(e,ldosu,ldosd,Jijint)
  use mod_f90_kind, only: double
  use mod_constants, only: pi, cZero, cOne, pauli_dorb
  use mod_parameters
  use mod_progress
  use mod_system, only: s => sys
  use TightBinding, only: nOrb,nOrb2
  use mod_magnet, only: hdel,mz
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer           :: nthreads,mythread
  integer             :: i,j,mu,nu,iz,alpha
  real(double)        :: kp(3),Jijkan(nmaglayers,3,3),Jijk(nmaglayers,nmaglayers,3,3)
  real(double),intent(in)     :: e
  real(double),intent(out)    :: ldosu(s%nAtoms, nOrb),ldosd(s%nAtoms, nOrb)
  real(double),intent(out)    :: Jijint(nmaglayers,nmaglayers,3,3)
  complex(double)     :: paulimatan(3,3,nOrb2, nOrb2)
  complex(double),dimension(nOrb2, nOrb2, s%nAtoms, s%nAtoms) :: gf
  complex(double),dimension(s%nAtoms, nOrb) :: gfdiagu,gfdiagd
  complex(double),dimension(nOrb2, nOrb2) :: gij,gji,temp1,temp2,paulia,paulib
  real(double),dimension(:,:),allocatable    :: ldosu_loc,ldosd_loc
  real(double),dimension(:,:,:,:),allocatable    :: Jijint_loc

! (x,y,z)-tensor formed by Pauli matrices to calculate anisotropy term (when i=j)
  paulimatan = cZero
  paulimatan(1,1,:,:) = -pauli_dorb(3,:,:)
  paulimatan(2,2,:,:) = -pauli_dorb(3,:,:)
  paulimatan(1,3,:,:) = -pauli_dorb(1,:,:)
  paulimatan(3,1,:,:) = -pauli_dorb(1,:,:)
  paulimatan(2,3,:,:) = -pauli_dorb(2,:,:)
  paulimatan(3,2,:,:) = -pauli_dorb(2,:,:)

  ldosu = 0.d0
  ldosd = 0.d0
  Jijint = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,gij,gji,paulia,paulib,i,j,mu,nu,alpha,gfdiagu,gfdiagd,Jijk,Jijkan,temp1,temp2,ldosu_loc,ldosd_loc,Jijint_loc) &
!$omp& shared(lverbose,s,e,eta,hdel,mz,nmaglayers,mmlayermag,pauli_dorb,paulimatan,ldosu,ldosd,Jijint,nthreads,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if(mythread==0) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[ldos_jij_energy] Number of threads: ',i0)") nthreads
!$  end if
allocate(ldosu_loc(s%nAtoms, nOrb), ldosd_loc(s%nAtoms, nOrb), Jijint_loc(nmaglayers,nmaglayers,3,3))
ldosu_loc = 0.d0
ldosd_loc = 0.d0
Jijint_loc = 0.d0
!reduction(+:ldosu,ldosd,Jijint)
!$omp do
do iz=1,s%nkpt
!$  if((mythread==0)) then
      if(lverbose) call progress_bar(outputunit_loop,"kpoints",iz,s%nkpt)
!$  end if
    kp = s%kbz(:,iz)

    ! Green function on energy E + ieta, and wave vector kp
    call green(e,eta,kp,gf)

    ! Exchange interaction tensor
    Jijk   = 0.d0
    Jijkan = 0.d0
    do nu = 1,3 ; do mu = 1,3 ; do j = 1,nmaglayers ; do i = 1,nmaglayers
      paulia = pauli_dorb(mu,:,:)
      gij = gf(:,:,mmlayermag(i)-1,mmlayermag(j)-1)
      paulib = pauli_dorb(nu,:,:)
      gji = gf(:,:,mmlayermag(j)-1,mmlayermag(i)-1)
      call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,paulia,nOrb2,gij,   nOrb2,cZero,temp1,nOrb2)
      call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp1, nOrb2,paulib,nOrb2,cZero,temp2,nOrb2)
      call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp2, nOrb2,gji,   nOrb2,cZero,temp1,nOrb2)
      do alpha = 1, nOrb2
        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
      end do

      Jijk(i,j,mu,nu) = -hdel(mmlayermag(i)-1)*hdel(mmlayermag(j)-1)*Jijk(i,j,mu,nu)*s%wkbz(iz)/(mz(mmlayermag(i)-1)*mz(mmlayermag(j)-1))

      ! Anisotropy (on-site) term
      if(i==j) then
        gij = gf(:,:,mmlayermag(i)-1,mmlayermag(i)-1)
        paulia = paulimatan(mu,nu,:,:)
        call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,gij,nOrb2,paulia,nOrb2,cZero,temp1,nOrb2)

        do alpha = 1,nOrb2
          Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
        end do

        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + hdel(mmlayermag(i)-1)*Jijkan(i,mu,nu)*s%wkbz(iz)/(mz(mmlayermag(i)-1)**2)
      end if
    end do ; end do ; end do ; end do

    ! Density of states
    do mu=1,nOrb; do i=1,s%nAtoms
      nu=mu+nOrb
      gfdiagu(i,mu) = - aimag(gf(mu,mu,i,i))*s%wkbz(iz)
      gfdiagd(i,mu) = - aimag(gf(nu,nu,i,i))*s%wkbz(iz)
    end do ; end do

    ldosu_loc = ldosu_loc + gfdiagu
    ldosd_loc = ldosd_loc + gfdiagd
    Jijint_loc = Jijint_loc + Jijk

  end do
!$omp end do
!$omp critical
ldosu = ldosu + ldosu_loc
ldosd = ldosd + ldosd_loc
Jijint = Jijint + Jijint_loc
!$omp end critical
!$omp end parallel

  ldosu  = ldosu/pi
  ldosd  = ldosd/pi
  Jijint = Jijint/pi

  return
end subroutine ldos_jij_energy
