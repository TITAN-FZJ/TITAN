!   Calculates magnetic LDOS
subroutine ldos(e,ldosu,ldosd)
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_generate_kpoints
	use mod_progress
	use mod_mpi_pars
!$ 	use omp_lib
	implicit none
!$	integer				:: nthreads,mythread
	integer         :: i,j,mu,nu,iz
	real(double)    :: ldosu(Npl,9),ldosd(Npl,9)
	real(double)    :: e,kp(3)
	complex(double),dimension(Npl,Npl,18,18)    :: gf

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

	ldosu = 0.d0
	ldosd = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,j,mu,nu) &
!$omp& shared(prog,elapsed_time,start_time,progbar,kbz,nkpoints,wkbz,e,eta,Npl,ldosu,ldosd,myrank,nthreads)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:ldosu,ldosd)
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if(myrank.eq.0) then
        ! Progress bar
        prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
        progress_bar: select case (mod(iz,4))
        case(0)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '|',prog,char(13)
        case(1)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '/',prog,char(13)
        case(2)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '-',prog,char(13)
        case(3)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '\',prog,char(13)
        end select progress_bar
#else
        elapsed_time = MPI_Wtime() - start_time
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
        write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
      end if
!$   end if
		kp = kbz(iz,:)

    ! Green function on energy E + ieta, and wave vector kp
		call green(e,eta,kp,gf)
		do mu=1,9; do i=1,Npl
			nu=mu+9
			ldosu(i,mu) = ldosu(i,mu) - aimag(gf(i,i,mu,mu))*wkbz(iz)
			ldosd(i,mu) = ldosd(i,mu) - aimag(gf(i,i,nu,nu))*wkbz(iz)
		end do ; end do
	end do kpoints
!$omp end do
!$omp end parallel

	ldosu  = ldosu/pi
	ldosd  = ldosd/pi

	return
end subroutine ldos

!   Calculates magnetic LDOS
subroutine ldos_es(e,ldosu,ldosd)
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_generate_kpoints
	use mod_progress
	use mod_mpi_pars
!$ 	use omp_lib
	implicit none
!$	integer				:: nthreads,mythread
	integer         :: i,j,mu,nu,iz
	real(double)    :: ldosu(Npl+2,9),ldosd(Npl+2,9)
	real(double)    :: e,kp(3)
	complex(double),dimension(Npl+2,Npl+2,18,18)    :: gf

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

	ldosu = 0.d0
	ldosd = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,j,mu,nu,nthreads) &
!$omp& shared(prog,elapsed_time,start_time,progbar,kbz,nkpoints,wkbz,e,eta,Npl,ldosu,ldosd,myrank)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:ldosu,ldosd)
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if(myrank.eq.0) then
        ! Progress bar
        prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
        progress_bar: select case (mod(iz,4))
        case(0)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '|',prog,char(13)
        case(1)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '/',prog,char(13)
        case(2)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '-',prog,char(13)
        case(3)
          write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '\',prog,char(13)
        end select progress_bar
#else
        elapsed_time = MPI_Wtime() - start_time
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
        write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
      end if
!$   end if
		kp = kbz(iz,:)

    ! Green function on energy E + ieta, and wave vector kp
		call green_es(e,eta,kp,gf)
		do mu=1,9; do i=1,Npl+2
			nu=mu+9
			ldosu(i,mu) = ldosu(i,mu) - aimag(gf(i,i,mu,mu))*wkbz(iz)
			ldosd(i,mu) = ldosd(i,mu) - aimag(gf(i,i,nu,nu))*wkbz(iz)
		end do ; end do
	end do kpoints
!$omp end do
!$omp end parallel

	ldosu  = ldosu/pi
	ldosd  = ldosd/pi

	return
end subroutine ldos_es