!   Calculates iso-energy surface (e=Ef for Fermi surface)
subroutine fermi_surface(e)
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
  character(len=50)  :: fieldpart,socpart
  character(len=1)   :: SOCc
	integer         	 :: i,mu,nu,iz
	real(double)    	 :: fsu_layer(Npl,nkpoints),fsd_layer(Npl,nkpoints),fsu_total(nkpoints),fsd_total(nkpoints)
  real(double),intent(in)    :: e
	real(double)    	 :: kp(3)
	complex(double),dimension(Npl,Npl,18,18)    :: gf

  write(outputunit_loop,"('CALCULATING FERMI SURFACE')")

! #ifndef _JUQUEEN
!   open(outputunit_loop,carriagecontrol ='fortran')
! #endif

  fsu_layer = 0.d0
  fsd_layer = 0.d0
  fsu_total = 0.d0
  fsd_total = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,mu,nu) &
!$omp& shared(prog,spiner,elapsed_time,start_program,progbar,e,kbz,kbz2d,nkpoints,Ef,eta,Npl,nthreads,pi,fsu_layer,fsd_layer,fsu_total,fsd_total,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if(mythread.eq.0) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[fermi_surface] Number of threads: ',i0)") nthreads
!$  end if

!$omp do
  fermi_surface_kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
    ! Progress bar
    prog = floor(iz*100.d0/nkpoints)
! #ifdef _JUQUEEN
    write(outputunit_loop,"(a1,2x,i3,'% (',i0,'/',i0,') of kpoints done ',a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,char(13)
! #else
!     elapsed_time = MPI_Wtime() - start_program
!     write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
!     write(outputunit_loop,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
! #endif
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

  fieldpart = ""
  socpart   = ""
  if(SOC) then
    if(llinearsoc) then
      SOCc = "L"
    else
      SOCc = "T"
    end if
    write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
  else
    SOCc = "F"
  end if
  if(lfield) then
    write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
    if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
  end if

  ! Writing on files
  do i=1,Npl
    write(varm,"('./results/',a1,'SOC/',i0,'Npl/FS/fsu_layer',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,i,ncp,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=17+(mpitag-1)*Npl+i, file=varm,status='unknown')
    write(varm,"('./results/',a1,'SOC/',i0,'Npl/FS/fsd_layer',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,i,ncp,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=57+(mpitag-1)*Npl+i, file=varm,status='unknown')
  end do
  write(varm,"('./results/',a1,'SOC/',i0,'Npl/FS/fsu_total_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,ncp,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=97+mpitag, file=varm,status='unknown')
  write(varm,"('./results/',a1,'SOC/',i0,'Npl/FS/fsd_total_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,ncp,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=98+mpitag, file=varm,status='unknown')

  writing_fermi_surface: do iz=1,nkpoints
    write_plane_loop_fs: do i=1,Npl
      write(unit=17+(mpitag-1)*Npl+i,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsu_layer(i,iz)
      write(unit=57+(mpitag-1)*Npl+i,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsd_layer(i,iz)
    end do write_plane_loop_fs

    write(unit=97+mpitag,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsu_total(iz)
    write(unit=98+mpitag,fmt="(3(f16.11,2x))") kbz2d(iz,1),kbz2d(iz,2),fsd_total(iz)
  end do writing_fermi_surface

  do i=1,Npl
    close (17+(mpitag-1)*Npl+i)
    close (57+(mpitag-1)*Npl+i)
  end do
  close(97+mpitag)
  close(98+mpitag)

	return
end subroutine fermi_surface