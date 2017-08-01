! Integration of Green functions over k values to calculate the number of particles
subroutine sumk_jij(er,ei,Jijint)
  use mod_f90_kind, only: double
  use mod_constants, only: pi, zum, zero, pauli_dorb
  use mod_parameters, only: mmlayermag, U, lverbose, q, mmlayermag, outputunit, nmaglayers
  use mod_magnet, only: mx,my,mz,mabs
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use mod_mpi_pars
  use mod_progress
!$  use omp_lib
  implicit none
!$  integer     :: nthreads,mythread
  integer       :: iz,i,j,mu,nu,alpha
  real(double),intent(in)   :: er,ei
  real(double),intent(out)  :: Jijint(nmaglayers,nmaglayers,3,3)
  real(double)  :: kp(3),kminusq(3),evec(3,nmaglayers),Jijk(nmaglayers,nmaglayers,3,3),Jijkan(nmaglayers,3,3)
  complex(double) :: dbxcdm(nmaglayers,3,2*nOrb,2*nOrb),d2bxcdm2(nmaglayers,3,3,2*nOrb,2*nOrb),paulievec(nmaglayers,2*nOrb,2*nOrb)
  complex(double),dimension(2*nOrb,2*nOrb)             :: gij,gji,temp1,temp2,paulia,paulib
  complex(double),dimension(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb)     :: gf,gfq

  do iz=1,nmaglayers
    ! Unit vector along the direction of the magnetization of each magnetic plane
    evec(:,iz) = [ mx(mmlayermag(iz)-1), my(mmlayermag(iz)-1), mz(mmlayermag(iz)-1) ]/mabs(mmlayermag(iz)-1)
    ! Inner product of pauli matrix in spin and orbital space and unit vector evec
    paulievec(iz,:,:) = pauli_dorb(1,:,:)*evec(1,iz)+pauli_dorb(2,:,:)*evec(2,iz)+pauli_dorb(3,:,:)*evec(3,iz)
  end do

! Derivative of Bxc*sigma*evec w.r.t. m_i (Bxc = -U.m/2)
  do i=1,3 ; do iz=1,nmaglayers
    dbxcdm(iz,i,:,:) = -0.5d0*U(mmlayermag(iz))*(pauli_dorb(i,:,:)-(paulievec(iz,:,:))*evec(i,iz))
  end do ; end do

! Second derivative of Bxc w.r.t. m_i (Bxc = -U.m/2)
  do j=1,3 ; do i = 1,3 ; do iz=1,nmaglayers
    d2bxcdm2(iz,i,j,:,:) = evec(i,iz)*pauli_dorb(j,:,:) + pauli_dorb(i,:,:)*evec(j,iz) - 3*paulievec(iz,:,:)*evec(i,iz)*evec(j,iz)
    if(i==j) d2bxcdm2(iz,i,j,:,:) = d2bxcdm2(iz,i,j,:,:) + paulievec(iz,:,:)
    d2bxcdm2(iz,i,j,:,:) = 0.5d0*U(mmlayermag(iz))*d2bxcdm2(iz,i,j,:,:)/(mabs(mmlayermag(iz)-1))
  end do ; end do ; end do

  Jijint = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,gfq,gij,gji,paulia,paulib,i,j,mu,nu,alpha,kminusq,Jijk,Jijkan,temp1,temp2) &
!$omp& shared(lverbose,s,q,er,ei,Jijint,dbxcdm,d2bxcdm2,mz,nmaglayers,mmlayermag,myrank,nthreads,outputunit)
!$  mythread = omp_get_thread_num()
!$  if((mythread==0).and.(myrank==0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit,"('[jij_energy] Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:Jijint)
  kpoints: do iz=1,s%nkpt
!$  if((mythread==0)) then
      if((myrank==0).and.(lverbose)) call progress_bar(outputunit,"kpoints",iz,s%nkpt)
!$   end if

    kp = s%kbz(:,iz)

    ! Green function on energy Ef + iy, and wave vector kp
    call green(er,ei,kp,gf)

    ! Green function on energy Ef + iy, and wave vector kp-q
    kminusq = kp-q
    call green(er,ei,kminusq,gfq)

    Jijk   = 0.d0
    Jijkan = 0.d0
    do nu = 1,3 ; do mu = 1,3 ; do j = 1,nmaglayers ; do i = 1,nmaglayers
      paulia = dbxcdm(i,mu,:,:)
      gij = gf(mmlayermag(i)-1,mmlayermag(j)-1,:,:)
      paulib = dbxcdm(j,nu,:,:)
      gji = gfq(mmlayermag(j)-1,mmlayermag(i)-1,:,:)
      call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,zum,paulia,2*nOrb,gij,   2*nOrb,zero,temp1,2*nOrb)
      call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,zum,temp1, 2*nOrb,paulib,2*nOrb,zero,temp2,2*nOrb)
      call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,zum,temp2, 2*nOrb,gji,   2*nOrb,zero,temp1,2*nOrb)
      ! Trace over orbitals and spins
      do alpha = 1,2*nOrb
        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
      end do

      ! Anisotropy (on-site) term
      if(i==j) then
        gij = gf(mmlayermag(i)-1,mmlayermag(i)-1,:,:)
        paulia = d2bxcdm2(i,mu,nu,:,:)
        call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,zum,gij,2*nOrb,paulia,2*nOrb,zero,temp1,2*nOrb)
        ! Trace over orbitals and spins
        do alpha = 1,2*nOrb
          Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
        end do

        Jijk(i,i,mu,nu) = Jijk(i,i,mu,nu) + Jijkan(i,mu,nu)
      end if
    end do ; end do ; end do ; end do

    Jijk = Jijk*s%wkbz(iz)

    Jijint = Jijint + Jijk

  end do kpoints
!$omp end do
!$omp end parallel

  Jijint = -Jijint/pi

!   write(outputunit,"('  ******************** *******  ********************')")
!   do i=1,nmaglayers ; do j=1,nmaglayers
!   ! Writing original full tensor Jij
!     if(i==j) then
!       write(outputunit,"(' |--------------- i = ',i0,'   j = ',i0,': anisotropies ---------------|')") i,j
!     else
!       write(outputunit,"(' |----------- i = ',i0,'   j = ',i0,': exchange couplings -------------|')") i,j
!     end if
!     write(outputunit,"('             x                  y                  z')")
!     write(outputunit,"('  x  (',es16.9,') (',es16.9,') (',es16.9,')')") Jijint(i,j,1,1),Jijint(i,j,1,2),Jijint(i,j,1,3)
!     write(outputunit,"('  y  (',es16.9,') (',es16.9,') (',es16.9,')')") Jijint(i,j,2,1),Jijint(i,j,2,2),Jijint(i,j,2,3)
!     write(outputunit,"('  z  (',es16.9,') (',es16.9,') (',es16.9,')')") Jijint(i,j,3,1),Jijint(i,j,3,2),Jijint(i,j,3,3)
!   end do ; end do

  return
end subroutine sumk_jij
