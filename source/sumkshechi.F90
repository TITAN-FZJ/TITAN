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
	integer         :: i,j,mu,nu,gamma,xi,iz
	integer, intent(in) :: iflag
	real(double),intent(in)			:: e,ep
	real(double)								:: kp(3)
	complex(double),dimension(:,:,:,:),allocatable 		:: gf
	complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd
	complex(double),dimension(dim,dim),intent(out)    :: Fint
	complex(double),dimension(:,:),allocatable  			:: df1

	Fint      = zero

!$omp parallel default(none) &
!$omp& private(errorcode,ierr,mythread,AllocateStatus,iz,kp,gf,gfuu,gfud,gfdu,gfdd,i,j,mu,nu,gamma,xi,df1) &
!$omp& shared(llineargfsoc,prog,spiner,lverbose,kbz,wkbz,e,ep,iflag,Fint,nkpoints,Ef,eta,nthreads,myrank,Npl,dim,sigmaimunu2i,sigmaijmunu2i)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if
	allocate(df1(dim,dim),gf(Npl,Npl,18,18),gfuu(Npl,Npl,9,9,2),gfud(Npl,Npl,9,9,2),gfdu(Npl,Npl,9,9,2),gfdd(Npl,Npl,9,9,2), STAT = AllocateStatus  )
  if (AllocateStatus.ne.0) then
    write(*,"('[sumkshechi] Not enough memory for: df1,gf,gfuu,gfud,gfdu,gfdd')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do schedule(auto)
	kpoints: do iz=1,nkpoints
    ! Progress bar
!$  if((mythread.eq.0)) then
      if((myrank.eq.0).and.(lverbose)) then
        prog = floor(iz*100.d0/nkpoints)
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
      end if
!$   end if

		kp = kbz(iz,:)

		if(iflag.eq.0)then
      ! Green function at (k+q,E_F+E+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef+e,ep,kp,gf)
      else
        call green(Ef+e,ep,kp,gf)
      end if
      gfuu(:,:,:,:,1) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,1) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,1) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,1) = gf(:,:,10:18,10:18)

      ! Green function at (k,E_F+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef,ep,kp,gf)
      else
        call green(Ef,ep,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,2) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,2) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,2) = gf(:,:,10:18,10:18)

!dir$ simd
      do xi=1,9 ; do gamma=1,9 ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do

		else
			! Green function at (k+q,E_F+E+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(ep+e,eta,kp,gf)
      else
        call green(ep+e,eta,kp,gf)
      end if
      gfuu(:,:,:,:,1) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,1) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,1) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,1) = gf(:,:,10:18,10:18)

			! Green function at (k,E_F+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(ep,eta,kp,gf)
      else
        call green(ep,eta,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,2) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,2) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,2) = gf(:,:,10:18,10:18)

!dir$ simd
      do xi=1,9 ; do gamma=1,9 ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do
		end if

    !$omp critical
#ifdef _JUQUEEN
    call zgeadd(Fint,dim,'N',df1,dim,'N',Fint,dim,dim,dim)
#else
    call ZAXPY(dim*dim,zum,df1,1,Fint,1)              !       Fint      = Fint + df1
!     Fint      = Fint + df1
#endif
    !$omp end critical
	end do kpoints
!$omp end do
	deallocate(df1)
	deallocate(gf,gfuu,gfud,gfdu,gfdd)
!$omp end parallel

	Fint     = Fint/tpi

	return
end subroutine sumkshechi


! ----------- Sum over wave vectors to calculate spin susceptibility -----------
! -------------- to be used in the calculation of linear SOC chi ---------------
subroutine sumkshechilinearsoc(e,ep,Fint,Fintlsoc,iflag)
  use mod_f90_kind
  use mod_parameters
  use mod_constants
  use mod_generate_kpoints
  use mod_progress
  use mod_mpi_pars
  use MPI
!$  use omp_lib
  implicit none
!$  integer         :: nthreads,mythread
  integer         :: AllocateStatus
  integer         :: i,j,mu,nu,gamma,xi,iz
  integer, intent(in) :: iflag
  real(double),intent(in)     :: e,ep
  real(double)                :: kp(3)
  complex(double),dimension(:,:,:,:),allocatable    :: gf,gvg
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd,gvguu,gvgud,gvgdu,gvgdd
  complex(double),dimension(dim,dim),intent(out)    :: Fint,Fintlsoc
  complex(double),dimension(:,:),allocatable        :: df1,df1lsoc

  Fint      = zero
  Fintlsoc  = zero

!$omp parallel default(none) &
!$omp& private(errorcode,ierr,mythread,AllocateStatus,iz,kp,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,i,j,mu,nu,gamma,xi,df1,df1lsoc) &
!$omp& shared(llineargfsoc,prog,spiner,lverbose,kbz,wkbz,e,ep,iflag,Fint,Fintlsoc,nkpoints,Ef,eta,nthreads,myrank,Npl,dim,sigmaimunu2i,sigmaijmunu2i)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if
  allocate(df1(dim,dim),gf(Npl,Npl,18,18),gfuu(Npl,Npl,9,9,2),gfud(Npl,Npl,9,9,2),gfdu(Npl,Npl,9,9,2),gfdd(Npl,Npl,9,9,2), STAT = AllocateStatus  )
  if (AllocateStatus.ne.0) then
    write(*,"('[sumkshechilinearsoc] Not enough memory for: df1,gf,gfuu,gfud,gfdu,gfdd')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  allocate( df1lsoc(dim,dim),gvg(Npl,Npl,18,18),gvguu(Npl,Npl,9,9,2),gvgud(Npl,Npl,9,9,2),gvgdu(Npl,Npl,9,9,2),gvgdd(Npl,Npl,9,9,2), STAT = AllocateStatus  )
  if (AllocateStatus.ne.0) then
    write(*,"('[sumkshechilinearsoc] Not enough memory for: df1lsoc,gvg,gvguu,gvgud,gvgdu,gvgdd')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do schedule(auto)
  kpoints: do iz=1,nkpoints
    ! Progress bar
!$  if((mythread.eq.0)) then
      if((myrank.eq.0).and.(lverbose)) then
        prog = floor(iz*100.d0/nkpoints)
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
      end if
!$   end if

    kp = kbz(iz,:)

    if(iflag.eq.0)then
      ! Green function at (k+q,E_F+E+iy)
      call greenlinearsoc(Ef+e,ep,kp,gf,gvg)
      gfuu(:,:,:,:,1) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,1) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,1) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,1) = gf(:,:,10:18,10:18)
      gvguu(:,:,:,:,1) = gvg(:,:, 1: 9, 1: 9)
      gvgud(:,:,:,:,1) = gvg(:,:, 1: 9,10:18)
      gvgdu(:,:,:,:,1) = gvg(:,:,10:18, 1: 9)
      gvgdd(:,:,:,:,1) = gvg(:,:,10:18,10:18)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(Ef,ep,kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,2) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,2) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,2) = gf(:,:,10:18,10:18)
      gvguu(:,:,:,:,2) = gvg(:,:, 1: 9, 1: 9)
      gvgud(:,:,:,:,2) = gvg(:,:, 1: 9,10:18)
      gvgdu(:,:,:,:,2) = gvg(:,:,10:18, 1: 9)
      gvgdd(:,:,:,:,2) = gvg(:,:,10:18,10:18)

!dir$ simd
      do xi=1,9 ; do gamma=1,9 ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*wkbz(iz)


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*wkbz(iz)

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*wkbz(iz)

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*wkbz(iz)

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do

    else
      ! Green function at (k+q,E_F+E+iy)
      call greenlinearsoc(ep+e,eta,kp,gf,gvg)
      gfuu(:,:,:,:,1) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,1) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,1) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,1) = gf(:,:,10:18,10:18)
      gvguu(:,:,:,:,1) = gvg(:,:, 1: 9, 1: 9)
      gvgud(:,:,:,:,1) = gvg(:,:, 1: 9,10:18)
      gvgdu(:,:,:,:,1) = gvg(:,:,10:18, 1: 9)
      gvgdd(:,:,:,:,1) = gvg(:,:,10:18,10:18)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(ep,eta,kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(:,:, 1: 9, 1: 9)
      gfud(:,:,:,:,2) = gf(:,:, 1: 9,10:18)
      gfdu(:,:,:,:,2) = gf(:,:,10:18, 1: 9)
      gfdd(:,:,:,:,2) = gf(:,:,10:18,10:18)
      gvguu(:,:,:,:,2) = gvg(:,:, 1: 9, 1: 9)
      gvgud(:,:,:,:,2) = gvg(:,:, 1: 9,10:18)
      gvgdu(:,:,:,:,2) = gvg(:,:,10:18, 1: 9)
      gvgdd(:,:,:,:,2) = gvg(:,:,10:18,10:18)

!dir$ simd
      do xi=1,9 ; do gamma=1,9 ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*wkbz(iz)


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*wkbz(iz)

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*wkbz(iz)

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*wkbz(iz)

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do
    end if

    !$omp critical
#ifdef _JUQUEEN
    call zgeadd(Fint,dim,'N',df1,dim,'N',Fint,dim,dim,dim)
    call zgeadd(Fintlsoc,dim,'N',df1lsoc,dim,'N',Fintlsoc,dim,dim,dim)
#else
    call ZAXPY(dim*dim,zum,df1,1,Fint,1)              !       Fint      = Fint + df1
    call ZAXPY(dim*dim,zum,df1lsoc,1,Fintlsoc,1)      !       Fintlsoc  = Fintlsoc + df1lsoc
!     Fint      = Fint + df1
#endif
    !$omp end critical
  end do kpoints
!$omp end do
  deallocate(df1)
  deallocate(gvg,gvguu,gvgud,gvgdu,gvgdd)
  deallocate(gf,gfuu,gfud,gfdu,gfdd)
!$omp end parallel

  Fint     = Fint/tpi
  Fintlsoc = Fintlsoc/tpi

  return
end subroutine sumkshechilinearsoc