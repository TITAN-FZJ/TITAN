!   Calculates iso-energy surface (e=Ef for Fermi surface)
subroutine fermisurface(e)
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_generate_kpoints
	use mod_progress
	use mod_mpi_pars
!$ 	use omp_lib
	implicit none
!$	integer					 :: nthreads,mythread
  character(len=400) :: varm
	integer         	 :: i,j,mu,nu,iz
	real(double)    	 :: fsu_layer(Npl,nkpoints),fsd_layer(Npl,nkpoints),fsu_total(nkpoints),fsd_total(nkpoints)
  real(double),intent(in)    :: e
	real(double)    	 :: kp(3)
	complex(double),dimension(Npl,Npl,18,18)    :: gf

#ifndef _JUQUEEN
      open(6,carriagecontrol ='fortran')
#endif

  fsu_layer = 0.d0
  fsd_layer = 0.d0
  fsu_total = 0.d0
  fsd_total = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,j,mu,nu) &
!$omp& shared(prog,spiner,elapsed_time,start_time,progbar,e,kbz,kbz2d,nkpoints,Ef,eta,Npl,nthreads,pi,fsu_layer,fsd_layer,fsu_total,fsd_total)
!$  mythread = omp_get_thread_num()
!$  if(mythread.eq.0) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do
  fermi_surface_kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
    ! Progress bar
    prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
    write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of kpoints done ',a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,char(13)
#else
    elapsed_time = MPI_Wtime() - start_time
    write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
    write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
!$   end if

    kp = kbz(iz,:)

    ! Green function on energy Ef + ieta, and wave vector kp
    call green(e,eta,kp,gf)

    do mu=1,9; do i=1,Npl
      nu=mu+9
      fsu_layer(i,iz) = fsu_layer(i,iz) - aimag(gf(i,i,mu,mu))/pi
      fsd_layer(i,iz) = fsd_layer(i,iz) - aimag(gf(i,i,nu,nu))/pi
    end do ; end do

    do i=1,Npl
      fsu_total(iz) = fsu_total(iz) + fsu_layer(i,iz)
      fsd_total(iz) = fsd_total(iz) + fsd_layer(i,iz)
    end do

  end do fermi_surface_kpoints
!$omp end do
!$omp end parallel

  ! Writing on files
  do i=1,Npl
    write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/FS/fsu_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,i,magaxis,socscale,ncp,eta,Utype,hwx,hwy,hwz
    open (unit=17+i, file=varm,status='unknown')
    write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/FS/fsd_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,i,magaxis,socscale,ncp,eta,Utype,hwx,hwy,hwz
    open (unit=57+i, file=varm,status='unknown')
  end do
  write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/FS/fsu_total_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,magaxis,socscale,ncp,eta,Utype,hwx,hwy,hwz
  open (unit=97, file=varm,status='unknown')
  write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/FS/fsd_total_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'.dat')") SOC,Npl,magaxis,socscale,ncp,eta,Utype,hwx,hwy,hwz
  open (unit=98, file=varm,status='unknown')

  writing_fermi_surface: do iz=1,nkpoints
    write_plane_loop_fs: do i=1,Npl
      write(unit=17+i,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsu_layer(i,iz)/ry2ev
      write(unit=57+i,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsd_layer(i,iz)/ry2ev
    end do write_plane_loop_fs

    write(unit=97,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsu_total(iz)/ry2ev
    write(unit=98,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsd_total(iz)/ry2ev
  end do writing_fermi_surface

  do i=1,Npl
    close (17+i)
    close (57+i)
  end do
  close(97)
  close(98)

	return
end subroutine fermisurface