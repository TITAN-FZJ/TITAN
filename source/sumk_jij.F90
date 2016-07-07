! Integration of Green functions over k values to calculate the number of particles
subroutine sumk_jij(er,ei,Jijint)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_magnet, only: mx,my,mz,mabs
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer         :: nthreads,mythread
  integer       :: iz,i,j,mu,nu,alpha
  real(double),intent(in)   :: er,ei
  real(double),intent(out)  :: Jijint(nmaglayers,nmaglayers,3,3)
  real(double)  :: kp(3),evec(nmaglayers,3),Jijk(nmaglayers,nmaglayers,3,3),Jijkan(nmaglayers,3,3)
  complex(double) :: pauli(3,18,18),dbxcdm(nmaglayers,3,18,18),d2bxcdm2(nmaglayers,3,3,18,18),paulievec(nmaglayers,18,18)
  complex(double),dimension(18,18)             :: gij,gji,temp1,temp2,paulia,paulib
  complex(double),dimension(Npl,Npl,18,18)     :: gf,gfq

#ifndef _JUQUEEN
  open(outputunit,carriagecontrol ='fortran')
#endif

! Pauli matrices in spin and orbital space
  pauli = 0.d0
  do mu = 5,9
    nu = mu+9
    ! pauli matrix x
    pauli(1,mu,nu) = zum
    pauli(1,nu,mu) = zum
    ! pauli matrix y
    pauli(2,mu,nu) = -zi
    pauli(2,nu,mu) = zi
    ! pauli matrix z
    pauli(3,mu,mu) = zum
    pauli(3,nu,nu) = -zum
  end do

  do iz=1,nmaglayers
    ! Unit vector along the direction of the magnetization of each magnetic plane
    evec(iz,:) = [ mx(mmlayermag(iz)-1), my(mmlayermag(iz)-1), mz(mmlayermag(iz)-1) ]/mabs(mmlayermag(iz)-1)
    ! Inner product of pauli matrix in spin and orbital space and unit vector evec
    paulievec(iz,:,:) = pauli(1,:,:)*evec(iz,1)+pauli(2,:,:)*evec(iz,2)+pauli(3,:,:)*evec(iz,3)
  end do

! Derivative of Bxc*sigma*evec w.r.t. m_i (Bxc = -U.m/2)
  do i=1,3 ; do iz=1,nmaglayers
    dbxcdm(iz,i,:,:) = -0.5d0*U(mmlayermag(iz))*(pauli(i,:,:)-(paulievec(iz,:,:))*evec(iz,i))
  end do ; end do

! Second derivative of Bxc w.r.t. m_i (Bxc = -U.m/2)
  do j=1,3 ; do i = 1,3 ; do iz=1,nmaglayers
    d2bxcdm2(iz,i,j,:,:) = evec(iz,i)*pauli(j,:,:) + pauli(i,:,:)*evec(iz,j) - 3*paulievec(iz,:,:)*evec(iz,i)*evec(iz,j)
    if(i.eq.j) d2bxcdm2(iz,i,j,:,:) = d2bxcdm2(iz,i,j,:,:) + paulievec(iz,:,:)
    d2bxcdm2(iz,i,j,:,:) = 0.5d0*U(mmlayermag(iz))*d2bxcdm2(iz,i,j,:,:)/(mabs(mmlayermag(iz)-1))
  end do ; end do ; end do

  Jijint = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,gfq,gij,gji,paulia,paulib,i,j,mu,nu,alpha,Jijk,Jijkan,temp1,temp2) &
!$omp& shared(prog,spiner,elapsed_time,start_program,progbar,kbz,q,wkbz,nkpoints,er,ei,Jijint,dbxcdm,d2bxcdm2,pauli,Npl,mz,nmaglayers,mmlayermag,myrank,nthreads,pi,outputunit)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit,"('[jij_energy] Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:Jijint)
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if(myrank.eq.0) then
        ! Progress bar
        prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
        write(outputunit,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
#else
        elapsed_time = MPI_Wtime() - start_program
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
        write(outputunit,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
      end if
!$   end if

    kp = kbz(iz,:)

    ! Green function on energy Ef + iy, and wave vector kp
    call green(er,ei,kp,gf)

    ! Green function on energy Ef + iy, and wave vector kp-q
    call green(er,ei,kp-q,gfq)

    Jijk   = 0.d0
    Jijkan = 0.d0
    do nu = 1,3 ; do mu = 1,3 ; do j = 1,nmaglayers ; do i = 1,nmaglayers
      paulia = dbxcdm(i,mu,:,:)
      gij = gf(mmlayermag(i)-1,mmlayermag(j)-1,:,:)
      paulib = dbxcdm(j,nu,:,:)
      gji = gfq(mmlayermag(j)-1,mmlayermag(i)-1,:,:)
      call zgemm('n','n',18,18,18,zum,paulia,18,gij,18,zero,temp1,18)
      call zgemm('n','n',18,18,18,zum,temp1,18,paulib,18,zero,temp2,18)
      call zgemm('n','n',18,18,18,zum,temp2,18,gji,18,zero,temp1,18)
      ! Trace over orbitals and spins
      do alpha = 1,18
        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
      end do

      ! Anisotropy (on-site) term
      if(i.eq.j) then
        gij = gf(mmlayermag(i)-1,mmlayermag(i)-1,:,:)
        paulia = d2bxcdm2(i,mu,nu,:,:)
        call zgemm('n','n',18,18,18,zum,gij,18,paulia,18,zero,temp1,18)
        ! Trace over orbitals and spins
        do alpha = 1,18
          Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
        end do

        Jijk(i,i,mu,nu) = Jijk(i,i,mu,nu) + Jijkan(i,mu,nu)
      end if
    end do ; end do ; end do ; end do

    Jijk = Jijk*wkbz(iz)

    Jijint = Jijint + Jijk

  end do kpoints
!$omp end do
!$omp end parallel

  Jijint = -Jijint/pi

!   write(outputunit,"('  ******************** *******  ********************')")
!   do i=1,nmaglayers ; do j=1,nmaglayers
!   ! Writing original full tensor Jij
!     if(i.eq.j) then
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