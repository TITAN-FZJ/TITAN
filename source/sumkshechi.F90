! ----------- Sum over wave vectors to calculate spin disturbance -----------
subroutine sumkshechi(e,ep,Fint,iflag)
	use mod_f90_kind
	use mod_parameters
	use mod_constants
	use mod_generate_kpoints
	use mod_progress
 	use mod_mpi_pars
 	use MPI
!$ 	use omp_lib
	implicit none
!$	integer					:: nthreads,mythread
	integer 				:: AllocateStatus
	integer         :: i,j,l,mu,nu,gamma,xi,iz
	integer, intent(in) :: iflag
	real(double),intent(in)			:: e,ep
	real(double)								:: kp(3)
	complex(double),dimension(:,:,:,:),allocatable 		:: gf
	complex(double),dimension(:,:,:,:,:),allocatable  :: gfmat,gfuu,gfud,gfdu,gfdd
	complex(double),dimension(dim,dim),intent(out)    :: Fint
	complex(double),dimension(:,:),allocatable  			:: df1

	Fint      = zero

!$omp parallel default(none) &
!$omp& private(errorcode,ierr,mythread,AllocateStatus,iz,kp,gf,gfmat,l,gfuu,gfud,gfdu,gfdd,i,j,mu,nu,gamma,xi,df1) &
!$omp& shared(prog,spiner,runoptions,kbz,wkbz,e,ep,iflag,Fint,nkpoints,Ef,eta,nthreads,myrank,Npl,dim,sigmaimunu2i,sigmaijmunu2i)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if
	allocate(df1(dim,dim),gf(Npl,Npl,18,18),gfuu(2,Npl,Npl,9,9),gfud(2,Npl,Npl,9,9),gfdu(2,Npl,Npl,9,9),gfdd(2,Npl,Npl,9,9),gfmat(2,Npl,Npl,18,18), STAT = AllocateStatus  )
  if (AllocateStatus.ne.0) then
    write(*,"('Not enough memory for: df1,gf,gfuu,gfud,gfdu,gfdd,gfmat')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do schedule(static)
	do iz=1,nkpoints
!    Progress bar
!$  if((mythread.eq.0)) then
      if((myrank.eq.0).and.(index(runoptions,"verbose").gt.0)) then
        prog = floor(iz*100.d0/nkpoints)
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
      end if
!$   end if

		kp = kbz(iz,:)

		if(iflag.eq.0)then
			df1 = zero
			gfmat = zero

      ! Green function at (k+q,E_F+E+iy)
      call green(Ef+e,ep,kp,gf)
      gfmat(1,:,:,:,:) = gf

      ! Green function at (k,E_F+iy)
      call green(Ef,ep,kp,gf)
      gfmat(2,:,:,:,:) = gf

      do l=1,2
        gfuu(l,:,:,:,:) = gfmat(l,:,:, 1: 9, 1: 9)
        gfud(l,:,:,:,:) = gfmat(l,:,:, 1: 9,10:18)
        gfdu(l,:,:,:,:) = gfmat(l,:,:,10:18, 1: 9)
        gfdd(l,:,:,:,:) = gfmat(l,:,:,10:18,10:18)
      end do

      do xi=1,9 ; do gamma=1,9 ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu)))*wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu)))*wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu)))*wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu)))*wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do

			!$omp critical
#ifdef _JUQUEEN
      call zgeadd(Fint,dim,'N',df1,dim,'N',Fint,dim,dim,dim)
#else
      call ZAXPY(dim*dim,zum,df1,1,Fint,1)              !       Fint      = Fint + df1
#endif
			!$omp end critical
		else
			df1 = zero
			gfmat = zero

			! Green function at (k+q,E_F+E+iy)
			call  green(ep+e,eta,kp,gf)
			gfmat(1,:,:,:,:) = gf

			! Green function at (k,E_F+iy)
			call green(ep,eta,kp,gf)
			gfmat(2,:,:,:,:) = gf

			do l=1,2
				gfuu(l,:,:,:,:) = gfmat(l,:,:, 1: 9, 1: 9)
				gfud(l,:,:,:,:) = gfmat(l,:,:, 1: 9,10:18)
				gfdu(l,:,:,:,:) = gfmat(l,:,:,10:18, 1: 9)
				gfdd(l,:,:,:,:) = gfmat(l,:,:,10:18,10:18)
			end do

      do xi=1,9 ; do gamma=1,9 ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))*wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))*wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))*wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))*wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do

			!$omp critical
#ifdef _JUQUEEN
      call zgeadd(Fint,dim,'N',df1,dim,'N',Fint,dim,dim,dim)
#else
      call ZAXPY(dim*dim,zum,df1,1,Fint,1)              !       Fint      = Fint + df1
#endif
			!$omp end critical
		end if
	end do
!$omp end do
	deallocate(df1)
	deallocate(gf,gfuu,gfud,gfdu,gfdd,gfmat)
!$omp end parallel

	Fint     = Fint/tpi

	return
end subroutine sumkshechi