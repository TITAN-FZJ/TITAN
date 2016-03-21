! Integration of Green functions over k values to calculate the number of particles
subroutine sumk_npart(er,ei,gdiaguur,gdiagddr,gdiagud,gdiagdu)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_mpi_pars
  use MPI
!$  use omp_lib
  implicit none
!$  integer         :: nthreads,mythread
  integer       :: iz,i,j,mu,mup
  real(double)  :: kp(3)
  real(double),intent(in)  :: er,ei
  real(double),dimension(Npl,9),intent(out)    :: gdiaguur,gdiagddr
  complex(double),dimension(Npl,9),intent(out) :: gdiagud,gdiagdu
  complex(double),dimension(Npl,Npl,18,18)     :: gf

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

  gdiaguur= 0.d0
  gdiagddr= 0.d0
  gdiagud = zero
  gdiagdu = zero

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,j,mu,mup) &
!$omp& shared(llineargfsoc,llinearsoc,prog,spiner,elapsed_time,start_program,progbar,lverbose,kbz,wkbz,nkpoints,er,ei,gdiaguur,gdiagddr,gdiagud,gdiagdu,Npl,myrank,nthreads)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if((myrank.eq.0).and.(lverbose)) then
        ! Progress bar
        prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of nparticles k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
#else
        elapsed_time = MPI_Wtime() - start_program
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
        write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
      end if
!$   end if

    kp = kbz(iz,:)

    ! Green function on energy Ef + iy, and wave vector kp
    if((llineargfsoc).or.(llinearsoc)) then
      call greenlineargfsoc(er,ei,kp,gf)
    else
      call green(er,ei,kp,gf)
    end if
    !$omp critical
    planes: do i=1,Npl
      orbital_index: do mu=1,9
        mup = mu+9
        gdiaguur(i,mu) = gdiaguur(i,mu) + real(gf(i,i,mu,mu)*wkbz(iz))
        gdiagddr(i,mu) = gdiagddr(i,mu) + real(gf(i,i,mup,mup)*wkbz(iz))
        gdiagud(i,mu) = gdiagud(i,mu) + (gf(i,i,mu,mup)*wkbz(iz))
        gdiagdu(i,mu) = gdiagdu(i,mu) + (gf(i,i,mup,mu)*wkbz(iz))
      end do orbital_index
    end do planes
    !$omp end critical
  end do kpoints
!$omp end do
!$omp end parallel
  return
end subroutine sumk_npart

! Integration of Green functions over k values to calculate the number of particles
subroutine sumk_npartjac(er,ei,ggr)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_mpi_pars
  use MPI
!$  use omp_lib
  implicit none
!$  integer         :: nthreads,mythread
  integer     :: iz,i,j,mu,nu
  real(double)  :: kp(3)
  real(double),intent(in)  :: er,ei
  real(double),dimension(Npl,Npl),intent(out)  :: ggr
  complex(double),dimension(Npl,Npl,18,18)     :: gf,gvg

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

  ggr    = 0.d0
!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,gvg,i,j,mu,nu) &
!$omp& shared(llineargfsoc,llinearsoc,prog,spiner,elapsed_time,start_program,progbar,lverbose,kbz,wkbz,nkpoints,er,ei,ggr,Npl,myrank,nthreads)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if((myrank.eq.0).and.(lverbose)) then
        ! Progress bar
        prog = floor(iz*100.d0/nkpoints)
#ifdef _JUQUEEN
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of jacobian k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
#else
        elapsed_time = MPI_Wtime() - start_program
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(iz+1)*20/nkpoints, "a,' ',i0,'%')"
        write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(iz+1)*20/nkpoints),100*(iz+1)/nkpoints
#endif
      end if
!$   end if

    kp = kbz(iz,:)

    ! Green function on energy Ef + iy, and wave vector kp
    if((llineargfsoc).or.(llinearsoc)) then
      call greenlinearsoc(er,ei,kp,gf,gvg)
      gf = gf + gvg
    else
      call green(er,ei,kp,gf)
    end if
    !$omp critical
    planes_i: do i=1,Npl
      planes_j: do j=1,Npl
        first_orbital_index: do mu=1,18
          second_orbital_index: do nu=1,18
            if ((mod(nu-1,9)+1).lt.5) cycle
            ggr(i,j) = ggr(i,j) + real(gf(i,j,mu,nu)*gf(j,i,nu,mu)*wkbz(iz))
            if((llineargfsoc).or.(llinearsoc)) ggr(i,j) = ggr(i,j) - real(gvg(i,j,mu,nu)*gvg(j,i,nu,mu)*wkbz(iz))
          end do second_orbital_index
        end do first_orbital_index
      end do planes_j
    end do planes_i
    !$omp end critical
  end do kpoints
!$omp end do
!$omp end parallel

  return
end subroutine sumk_npartjac