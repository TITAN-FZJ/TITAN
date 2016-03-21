! -------- sum over wave vectors to calculate parallel spin current --------
subroutine sumkprllsc(e,ep,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,iflag)
  use mod_f90_kind
  use mod_parameters
  use mod_constants
  use mod_lattice
  use mod_generate_kpoints
  use mod_tight_binding, only: t00
  use mod_prefactors
  use mod_progress
  use mod_mpi_pars
  use MPI
!$  use omp_lib
  implicit none
!$  integer       :: nthreads,mythread
  integer         :: AllocateStatus
  integer         :: i,j,l,mu,nu,gamma,xi,sigma,sigmap
  integer         :: iz,neighbor
  integer,      intent(in) :: iflag
  real(double), intent(in) :: e,ep
  real(double)             :: kp(3)
  complex(double)          :: expikr(n0sc1:n0sc2),wkbzc
  complex(double),dimension(:,:,:,:),allocatable    :: gf,dtdk
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfmat,gfuu,gfud,gfdu,gfdd
  complex(double),dimension(:,:),allocatable        :: df1iikl,pfdf1iikl
  complex(double),dimension(n0sc1:n0sc2,Npl,9,9)    :: prett,preLxtt,preLytt,preLztt
  complex(double),dimension(n0sc1:n0sc2,dimsigmaNpl,4),  intent(out) :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl

  ttFintiikl   = zero
  LxttFintiikl = zero
  LyttFintiikl = zero
  LzttFintiikl = zero

!$omp parallel default(none) &
!$omp& private(errorcode,ierr,mythread,AllocateStatus,iz,wkbzc,kp,df1iikl,pfdf1iikl,prett,preLxtt,preLytt,preLztt,dtdk,gf,gfmat,expikr,gfuu,gfud,gfdu,gfdd,sigma,sigmap,i,j,l,mu,nu,gamma,xi,neighbor) &
!$omp& shared(ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,prefactor,prog,spiner,lverbose,myrank,kbz,wkbz,iflag,e,ep,nkpoints,r0,Ef,eta,nthreads,sigmai2i,sigmaimunu2i,sigmaijmunu2i,dim,Npl,n0sc1,n0sc2,plnn,t00,lxpt,lypt,lzpt,tlxp,tlyp,tlzp)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if
  allocate( df1iikl(dim,4),pfdf1iikl(dim,4),gf(Npl,Npl,18,18),dtdk(Npl,Npl,9,9),gfuu(2,Npl,Npl,9,9),gfud(2,Npl,Npl,9,9),gfdu(2,Npl,Npl,9,9),gfdd(2,Npl,Npl,9,9),gfmat(2,Npl,Npl,18,18), STAT = AllocateStatus  )
  if (AllocateStatus.ne.0) then
    write(*,"('[sumkprllsc] Not enough memory for: df1iikl,,pfdf1iikl,gf,dtdk,gfuu,gfud,gfdu,gfdd,gfmat')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do schedule(static)
  do iz=1,nkpoints
!    Progress bar
!$  if((mythread.eq.0)) then
      if((myrank.eq.0).and.(lverbose)) then
        prog = floor(iz*100.d0/nkpoints)
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of k-sum on rank ',i0,a1,$)") spiner(mod(iz,4)+1),prog,iz,nkpoints,myrank,char(13)
      end if
!$   end if

    kp = kbz(iz,:)

    df1iikl = zero
    gfmat   = zero

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(kp,dtdk)

    do neighbor=n0sc1,n0sc2
      expikr(neighbor) = exp(zi*(kp(1)*r0(neighbor,1)+kp(2)*r0(neighbor,2)+kp(3)*r0(neighbor,3)))
    end do

    ! Calculating the prefactor (L).t.exp - t.(L).exp
    do nu=1,9 ; do mu=1,9 ; do neighbor=n0sc1,n0sc2 ; do i=1,Npl
      prett  (neighbor,i,mu,nu) = t00(i+1,neighbor,mu,nu)*expikr(neighbor)-(t00(i+1,neighbor,nu,mu)*conjg(expikr(neighbor)))
      preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    end do ; end do ; end do ; end do

    if(iflag.eq.0)then
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

      do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        if(abs(j-l).gt.plnn) cycle
        df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1) + (gfdd(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2) + (gfdu(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3) + (gfdd(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4) + (gfdu(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1) + (gfud(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2) + (gfuu(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3) + (gfud(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4) + (gfuu(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1) + (gfdd(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2) + (gfdu(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3) + (gfdd(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4) + (gfdu(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1) + (gfud(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2) + (gfuu(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3) + (gfud(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4) + (gfuu(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)
      end do ; end do ; end do ; end do ; end do ; end do ; end do

    else

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

      do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        if(abs(j-l).gt.plnn) cycle
        df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1)-zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2)-zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3)-zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4)-zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1)-zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2)-zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3)-zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4)-zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1)-zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2)-zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3)-zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4)-zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1)-zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2)-zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3)-zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4)-zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)
      end do ; end do ; end do ; end do ; end do ; end do ; end do

    end if

    ! Multiplying the prefactor by the susceptibility and the k-point weight
    ! pfdf1iikl = prefactor*df1iikl*wkbz(iz)
    wkbzc = cmplx(wkbz(iz),0.d0)
    call zgemm('n','n',dim,4,dim,wkbzc,prefactor,dim,df1iikl,dim,zero,pfdf1iikl,dim)

    !$omp critical
    do sigmap=1,4 ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4 ; do neighbor=n0sc1,n0sc2
      ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) = ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) + (prett  (neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
      ! Orbital angular momentum currents
      LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLxtt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
      LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLytt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
      LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLztt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    end do ; end do ; end do ; end do ; end do ; end do
    !$omp end critical

  end do ! iz
!$omp end do
  deallocate(df1iikl,pfdf1iikl)
  deallocate(gf,dtdk,gfuu,gfud,gfdu,gfdd,gfmat)
!$omp end parallel

  ttFintiikl   = ttFintiikl/tpi
  LxttFintiikl = LxttFintiikl/tpi
  LyttFintiikl = LyttFintiikl/tpi
  LzttFintiikl = LzttFintiikl/tpi
  write(*,*) myrank,sum(abs(ttFintiikl))
  write(*,*) myrank,sum(abs(LxttFintiikl))
  write(*,*) myrank,sum(abs(LyttFintiikl))
  write(*,*) myrank,sum(abs(LzttFintiikl))

  return
end subroutine sumkprllsc