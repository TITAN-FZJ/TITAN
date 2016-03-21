!   Calculates spin-resolved LDOS and energy-dependence of exchange interactions
subroutine ldos(e,ldosu,ldosd,Jijint)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_magnet, only: hdel,mz
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer       :: nthreads,mythread
  integer             :: i,j,mu,nu,iz,alpha
  real(double)        :: kp(3),Jijkan(nmaglayers,3,3),Jijk(nmaglayers,nmaglayers,3,3)
  real(double),intent(in)     :: e
  real(double),intent(out)    :: ldosu(Npl,9),ldosd(Npl,9)
  real(double),intent(out)    :: Jijint(nmaglayers,nmaglayers,3,3)
  complex(double)     :: pauli(3,18,18),paulimatan(3,3,18,18)
  complex(double),dimension(Npl,Npl,18,18)    :: gf
  complex(double),dimension(Npl,9)            :: gfdiagu,gfdiagd
  complex(double),dimension(18,18)            :: gij,gji,temp1,temp2,paulia,paulib

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

  ldosu = 0.d0
  ldosd = 0.d0
  Jijint = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,gij,gji,paulia,paulib,i,j,mu,nu,alpha,gfdiagu,gfdiagd,Jijk,Jijkan,temp1,temp2) &
!$omp& shared(prog,spiner,elapsed_time,start_time,progbar,kbz,nkpoints,wkbz,e,eta,Npl,hdel,mz,nmaglayers,mmlayermag,pauli,paulimatan,ldosu,ldosd,Jijint,myrank,nthreads)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:ldosu,ldosd,Jijint)
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      ! Progress bar
      prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
      write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum ',a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,char(13)
#else
      elapsed_time = MPI_Wtime() - start_time
      write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
!$   end if
    kp = kbz(iz,:)

    ! Green function on energy E + ieta, and wave vector kp
    call green(e,eta,kp,gf)

    ! Exchange interaction tensor
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

      Jijk(i,j,mu,nu) = -hdel(mmlayermag(i)-1)*hdel(mmlayermag(j)-1)*Jijk(i,j,mu,nu)*wkbz(iz)/(mz(mmlayermag(i)-1)*mz(mmlayermag(j)-1))

      ! Anisotropy (on-site) term
      if(i.eq.j) then
        gij = gf(mmlayermag(i)-1,mmlayermag(i)-1,:,:)
        paulia = paulimatan(mu,nu,:,:)
        call zgemm('n','n',18,18,18,zum,gij,18,paulia,18,zero,temp1,18)

        do alpha = 1,18
          Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
        end do

        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + hdel(mmlayermag(i)-1)*Jijkan(i,mu,nu)*wkbz(iz)/(mz(mmlayermag(i)-1)**2)
      end if
    end do ; end do ; end do ; end do

    ! Density of states
    do mu=1,9; do i=1,Npl
      nu=mu+9
      gfdiagu(i,mu) = - aimag(gf(i,i,mu,mu))*wkbz(iz)
      gfdiagd(i,mu) = - aimag(gf(i,i,nu,nu))*wkbz(iz)
    end do ; end do

    ldosu = ldosu + gfdiagu
    ldosd = ldosd + gfdiagd
    Jijint = Jijint + Jijk

  end do kpoints
!$omp end do
!$omp end parallel

  ldosu  = ldosu/(pi*ry2ev)
  ldosd  = ldosd/(pi*ry2ev)
  Jijint = Jijint*ry2ev/pi

  return
end subroutine ldos

!   Calculates spin-resolved LDOS (including empty spheres)
subroutine ldos_es(e)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer       :: nthreads,mythread
  character(len=400) :: varm
  integer         :: i,j,mu,nu,iz
  real(double)    :: ldosu(Npl+2,9),ldosd(Npl+2,9)
  real(double)    :: e,kp(3)
  complex(double),dimension(Npl,9)                :: gfdiagu,gfdiagd
  complex(double),dimension(Npl+2,Npl+2,18,18)    :: gf

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

  ldosu = 0.d0
  ldosd = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,j,mu,nu,nthreads) &
!$omp& shared(prog,spiner,elapsed_time,start_time,progbar,kbz,nkpoints,wkbz,e,eta,Npl,gfdiagu,gfdiagd,ldosu,ldosd,myrank)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:ldosu,ldosd)
  kpoints: do iz=1,nkpoints
    ! Progress bar
!$  if((mythread.eq.0)) then
      prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
      write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum ',a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,char(13)
#else
      elapsed_time = MPI_Wtime() - start_time
      write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
!$   end if

    kp = kbz(iz,:)

    ! Green function on energy E + ieta, and wave vector kp
    call green_es(e,eta,kp,gf)

    do mu=1,9; do i=1,Npl+2
      nu=mu+9
      gfdiagu(i,mu) = - aimag(gf(i,i,mu,mu))*wkbz(iz)
      gfdiagd(i,mu) = - aimag(gf(i,i,nu,nu))*wkbz(iz)
    end do ; end do

    ldosu = ldosu + gfdiagu
    ldosd = ldosd + gfdiagd

  end do kpoints
!$omp end do
!$omp end parallel

  ldosu  = ldosu/(pi*ry2ev)
  ldosd  = ldosd/(pi*ry2ev)

  ! Transform energy to eV if runoption is on
  e = e*ry2ev

  !Writing on files
  do i=1,Npl+2
    write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/LDOS/ldosu_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,i,magaxis,socscale,ncp,eta,Utype,hwx,hwy,hwz
    open (unit=117+i, file=varm,status='unknown')
    write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/LDOS/ldosd_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,i,magaxis,socscale,ncp,eta,Utype,hwx,hwy,hwz
    open (unit=517+i, file=varm,status='unknown')
  end do
  ldos_writing_plane_loop: do i=1,Npl+2
      write(unit=117+i,fmt="(5(e16.9,2x))") e,sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
      write(unit=517+i,fmt="(5(e16.9,2x))") e,sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
  end do ldos_writing_plane_loop

  do i=1,Npl+2
    close (117+i)
    close (517+i)
  end do

  return
end subroutine ldos_es