! Integration of Green functions over k values to calculate the number of particles
subroutine sumk_jij(er,ei,Jijint)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_magnet, only: hdel,mag
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer         :: nthreads,mythread
  integer       :: iz,i,j,mu,nu,alpha
  real(double),intent(in)   :: er,ei
  real(double),intent(out)  :: Jijint(nmaglayers,nmaglayers,3,3)
  real(double)  :: kp(3),Jijk(nmaglayers,nmaglayers,3,3),Jijkan(nmaglayers,3,3)
  complex(double) :: pauli(3,18,18),paulimatan(3,3,18,18)
  complex(double),dimension(18,18)             :: gij,gji,temp1,temp2,paulia,paulib
  complex(double),dimension(Npl,Npl,18,18)     :: gf

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
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
! (x,y,z)-tensor formed by Pauli matrices to calculate anisotropy term (when i=j)
  paulimatan = zero
  paulimatan(1,1,:,:) = -pauli(3,:,:)
  paulimatan(2,2,:,:) = -pauli(3,:,:)
  paulimatan(1,3,:,:) = -pauli(1,:,:)
  paulimatan(3,1,:,:) = -pauli(1,:,:)
  paulimatan(2,3,:,:) = -pauli(2,:,:)
  paulimatan(3,2,:,:) = -pauli(2,:,:)

  Jijint = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,gij,gji,paulia,paulib,i,j,mu,nu,alpha,Jijk,Jijkan,temp1,temp2) &
!$omp& shared(prog,spiner,elapsed_time,start_time,progbar,kbz,wkbz,nkpoints,er,ei,Jijint,pauli,paulimatan,Npl,hdel,mag,nmaglayers,mmlayermag,myrank,nthreads)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if(myrank.eq.0) then
        ! Progress bar
        prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
#else
        elapsed_time = MPI_Wtime() - start_time
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
        write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
      end if
!$   end if

    kp = kbz(iz,:)

    ! Green function on energy Ef + iy, and wave vector kp
    call green(er,ei,kp,gf)

    Jijk   = 0.d0
    Jijkan = 0.d0
    do nu = 1,3 ; do mu = 1,3 ; do j = 1,nmaglayers ; do i = 1,nmaglayers
      paulia = pauli(mu,:,:)
      gij = gf(mmlayermag(i)-1,mmlayermag(j)-1,:,:)
      paulib = pauli(nu,:,:)
      gji = gf(mmlayermag(j)-1,mmlayermag(i)-1,:,:)
      call zgemm('n','n',18,18,18,zum,paulia,18,gij,18,zero,temp1,18)
      call zgemm('n','n',18,18,18,zum,temp1,18,paulib,18,zero,temp2,18)
      call zgemm('n','n',18,18,18,zum,temp2,18,gji,18,zero,temp1,18)
      do alpha = 1,18
        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
      end do

      Jijk(i,j,mu,nu) = -hdel(mmlayermag(i)-1)*hdel(mmlayermag(j)-1)*Jijk(i,j,mu,nu)*wkbz(iz)/(mag(mmlayermag(i)-1)*mag(mmlayermag(j)-1))

      ! Anisotropy (on-site) term
      if(i.eq.j) then
        gij = gf(mmlayermag(i)-1,mmlayermag(i)-1,:,:)
        paulia = paulimatan(mu,nu,:,:)
        call zgemm('n','n',18,18,18,zum,gij,18,paulia,18,zero,temp1,18)

        do alpha = 1,18
          Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
        end do

        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + hdel(mmlayermag(i)-1)*Jijkan(i,mu,nu)*wkbz(iz)/(mag(mmlayermag(i)-1)**2)
      end if
    end do ; end do ; end do ; end do

    !$omp critical
    Jijint = Jijint + Jijk
    !$omp end critical
  end do kpoints
!$omp end do
!$omp end parallel

!   write(*,"('  ******************** *******  ********************')")
!   do i=1,nmaglayers ; do j=1,nmaglayers
!   ! Writing original full tensor Jij
!     if(i.eq.j) then
!       write(*,"(' |--------------- i = ',i0,'   j = ',i0,': anisotropies ---------------|')") i,j
!     else
!       write(*,"(' |----------- i = ',i0,'   j = ',i0,': exchange couplings -------------|')") i,j
!     end if
!     write(*,"('             x                  y                  z')")
!     write(*,"('  x  (',e16.9,') (',e16.9,') (',e16.9,')')") Jijint(i,j,1,1),Jijint(i,j,1,2),Jijint(i,j,1,3)
!     write(*,"('  y  (',e16.9,') (',e16.9,') (',e16.9,')')") Jijint(i,j,2,1),Jijint(i,j,2,2),Jijint(i,j,2,3)
!     write(*,"('  z  (',e16.9,') (',e16.9,') (',e16.9,')')") Jijint(i,j,3,1),Jijint(i,j,3,2),Jijint(i,j,3,3)
!   end do ; end do

  return
end subroutine sumk_jij