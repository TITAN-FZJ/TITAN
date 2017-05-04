! -------- sum over wave vectors to calculate parallel spin current --------
subroutine sumk(e,ep,tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,iflag)
  use mod_f90_kind
  use mod_parameters
  use mod_constants
  use mod_system, only: l_nn, r_nn, npln, nkpt, kbz, wkbz
  use mod_tight_binding, only: t00
  use mod_prefactors
  use mod_progress
  use mod_mpi_pars
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
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd
  complex(double),dimension(:,:),allocatable        :: df1iikl,pfdf1iikl
  complex(double),dimension(n0sc1:n0sc2,Npl,9,9)    :: prett,preLxtt,preLytt,preLztt
  complex(double),dimension(n0sc1:n0sc2,dimsigmaNpl,4), intent(out) :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl
  complex(double),dimension(dim,4),                     intent(out) :: tFintiikl

  tFintiikl    = zero
  ttFintiikl   = zero
  LxttFintiikl = zero
  LyttFintiikl = zero
  LzttFintiikl = zero

!$omp parallel default(none) &
!$omp& private(errorcode,ierr,mythread,AllocateStatus,iz,wkbzc,kp,df1iikl,pfdf1iikl,prett,preLxtt,preLytt,preLztt,dtdk,gf,expikr,gfuu,gfud,gfdu,gfdd,sigma,sigmap,i,j,l,mu,nu,gamma,xi,neighbor) &
!$omp& shared(llineargfsoc,tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,prefactor,lverbose,myrank_row_hw,kbz,wkbz,iflag,e,ep,nkpoints,r0,Ef,eta,nthreads,sigmai2i,sigmaimunu2i,sigmaijmunu2i,dim,Npl,n0sc1,n0sc2,npln,t00,lxpt,lypt,lzpt,tlxp,tlyp,tlzp,outputunit,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if((mythread==0).and.(myrank_row_hw==0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[sumk] Number of threads: ',i0)") nthreads
!$  end if
  allocate( df1iikl(dim,4),pfdf1iikl(dim,4),gf(Npl,Npl,18,18),dtdk(Npl,Npl,9,9),gfuu(Npl,Npl,9,9,2),gfud(Npl,Npl,9,9,2),gfdu(Npl,Npl,9,9,2),gfdd(Npl,Npl,9,9,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) then
    write(outputunit,"('[sumk] Not enough memory for: df1iikl,pfdf1iikl,gf,dtdk,gfuu,gfud,gfdu,gfdd')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do schedule(auto)
  do iz=1,nkpt
!    Progress bar
!$  if((mythread==0)) then
      if((myrank_row_hw==0).and.(lverbose)) call progress_bar(outputunit_loop,"kpoints",iz,nkpt)
!$   end if

    kp = kbz(:,iz)

    df1iikl = zero

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(kp,dtdk)

    do neighbor=n0sc1,n0sc2
      expikr(neighbor) = exp(zi * dot_product(kp, r_nn(:,neighbor)))
    end do

    ! Calculating the prefactor (L).t.exp - t.(L).exp
    do nu=1,9 ; do mu=1,9 ; do neighbor=n0sc1,n0sc2 ; do i=1,Npl
      prett  (neighbor,i,mu,nu) = t00(i+1,neighbor,mu,nu)*expikr(neighbor)-(t00(i+1,neighbor,nu,mu)*conjg(expikr(neighbor)))
      preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    end do ; end do ; end do ; end do

    if(iflag==0)then
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

      do nu=1,9 ; do mu=1,9 ; do i=1,Npl; do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl
        if(abs(j-l) > npln) cycle
        df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1) + (gfdd(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2) + (gfdu(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3) + (gfdd(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4) + (gfdu(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1) + (gfud(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2) + (gfuu(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3) + (gfud(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4) + (gfuu(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1) + (gfdd(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2) + (gfdu(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3) + (gfdd(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4) + (gfdu(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1) + (gfud(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2) + (gfuu(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3) + (gfud(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4) + (gfuu(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
      end do ; end do ; end do ; end do ; end do ; end do ; end do

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

      do nu=1,9 ; do gamma=1,9 ; do mu=1,9 ; do i=1,Npl ; do xi=1,9 ; do l=1,Npl ; do j=1,Npl
        if(abs(j-l) > npln) cycle
        df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
      end do ; end do ; end do ; end do ; end do ; end do ; end do

    end if

    ! Multiplying the prefactor by the susceptibility and the k-point weight
    ! pfdf1iikl = prefactor*df1iikl
    wkbzc = cmplx(wkbz(iz),0.d0)
    call zgemm('n','n',dim,4,dim,wkbzc,prefactor,dim,df1iikl,dim,zero,pfdf1iikl,dim)

    !$omp critical

    ! Integrating matrices

    ! tFintiikl = tFintiikl + pfdf1iikl
#ifdef _JUQUEEN
    call zgeadd(tFintiikl,dim,'N',pfdf1iikl,dim,'N',tFintiikl,dim,dim,4)
#else
    call ZAXPY(dim*4,zum,pfdf1iikl,1,tFintiikl,1)
!       tFintiikl = tFintiikl + pfdf1iikl
#endif
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
  deallocate(gf,dtdk,gfuu,gfud,gfdu,gfdd)
!$omp end parallel

  tFintiikl    = tFintiikl/tpi
  ttFintiikl   = ttFintiikl/tpi
  LxttFintiikl = LxttFintiikl/tpi
  LyttFintiikl = LyttFintiikl/tpi
  LzttFintiikl = LzttFintiikl/tpi

  return
end subroutine sumk


! -------- sum over wave vectors to calculate parallel spin current --------
subroutine sumklinearsoc(e,ep,tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,iflag)
  use mod_f90_kind
  use mod_parameters
  use mod_constants
  use mod_system, only: l_nn, r_nn, npln, nkpt, kbz, wkbz
  use mod_tight_binding, only: t00
  use mod_prefactors
  use mod_progress
  use mod_mpi_pars, only: myrank_row_hw,errorcode,ierr
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
  complex(double),dimension(:,:,:,:),allocatable    :: gf,dtdk,gvg
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd,gvguu,gvgud,gvgdu,gvgdd
  complex(double),dimension(:,:),allocatable        :: df1iikl,pfdf1iikl,df1lsoc
  complex(double),dimension(n0sc1:n0sc2,Npl,9,9)    :: prett,preLxtt,preLytt,preLztt
  complex(double),dimension(n0sc1:n0sc2,dimsigmaNpl,4), intent(out) :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl
  complex(double),dimension(dim,4),                     intent(out) :: tFintiikl

  tFintiikl    = zero
  ttFintiikl   = zero
  LxttFintiikl = zero
  LyttFintiikl = zero
  LzttFintiikl = zero

!$omp parallel default(none) &
!$omp& private(errorcode,ierr,mythread,AllocateStatus,iz,wkbzc,kp,df1iikl,pfdf1iikl,df1lsoc,prett,preLxtt,preLytt,preLztt,dtdk,expikr,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,sigma,sigmap,i,j,l,mu,nu,gamma,xi,neighbor) &
!$omp& shared(tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,prefactor,prefactorlsoc,lverbose,myrank_row_hw,kbz,wkbz,iflag,e,ep,nkpoints,r0,Ef,eta,nthreads,sigmai2i,sigmaimunu2i,sigmaijmunu2i,dim,Npl,n0sc1,n0sc2,npln,t00,lxpt,lypt,lzpt,tlxp,tlyp,tlzp,outputunit,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if((mythread==0).and.(myrank_row_hw==0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[sumklinearsoc] Number of threads: ',i0)") nthreads
!$  end if
  allocate( df1iikl(dim,4),pfdf1iikl(dim,4),df1lsoc(dim,4),gf(Npl,Npl,18,18),dtdk(Npl,Npl,9,9),gfuu(Npl,Npl,9,9,2),gfud(Npl,Npl,9,9,2),gfdu(Npl,Npl,9,9,2),gfdd(Npl,Npl,9,9,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) then
    write(outputunit,"('[sumklinearsoc] Not enough memory for: df1iikl,pfdf1iikl,df1lsoc,gf,dtdk,gfuu,gfud,gfdu,gfdd')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  allocate( gvg(Npl,Npl,18,18),gvguu(Npl,Npl,9,9,2),gvgud(Npl,Npl,9,9,2),gvgdu(Npl,Npl,9,9,2),gvgdd(Npl,Npl,9,9,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) then
    write(outputunit,"('[sumklinearsoc] Not enough memory for: gvg,gvguu,gvgud,gvgdu,gvgdd')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do schedule(auto)
  do iz=1,nkpt
!    Progress bar
!$  if((mythread==0)) then
      if((myrank_row_hw==0).and.(lverbose)) call progress_bar(outputunit_loop,"kpoints",iz,nkpt)
!$   end if

    kp = kbz(:,iz)

    df1iikl = zero
    df1lsoc = zero

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(kp,dtdk)

    do neighbor=n0sc1,n0sc2
      expikr(neighbor) = exp(zi*dot_product(kp, r_nn(:,neighbor)))
    end do

    ! Calculating the prefactor (L).t.exp - t.(L).exp
    do nu=1,9 ; do mu=1,9 ; do neighbor=n0sc1,n0sc2 ; do i=1,Npl
      prett  (neighbor,i,mu,nu) = t00(i+1,neighbor,mu,nu)*expikr(neighbor)-(t00(i+1,neighbor,nu,mu)*conjg(expikr(neighbor)))
      preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    end do ; end do ; end do ; end do

    if(iflag==0)then
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

      do nu=1,9 ; do mu=1,9 ; do i=1,Npl; do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl
        if(abs(j-l) > npln) cycle
        df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1) + (gfdd(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2) + (gfdu(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3) + (gfdd(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4) + (gfdu(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1) + (gfud(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2) + (gfuu(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3) + (gfud(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4) + (gfuu(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1) + (gfdd(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2) + (gfdu(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3) + (gfdd(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4) + (gfdu(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1) + (gfud(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2) + (gfuu(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3) + (gfud(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4) + (gfuu(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)


        df1lsoc(sigmaimunu2i(1,i,mu,nu),1) = df1lsoc(sigmaimunu2i(1,i,mu,nu),1) + (gvgdd(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gvguu(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvguu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),2) = df1lsoc(sigmaimunu2i(1,i,mu,nu),2) + (gvgdu(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gvguu(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvguu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),3) = df1lsoc(sigmaimunu2i(1,i,mu,nu),3) + (gvgdd(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gvgud(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),4) = df1lsoc(sigmaimunu2i(1,i,mu,nu),4) + (gvgdu(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gvgud(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1lsoc(sigmaimunu2i(2,i,mu,nu),1) = df1lsoc(sigmaimunu2i(2,i,mu,nu),1) + (gvgud(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gvguu(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvguu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),2) = df1lsoc(sigmaimunu2i(2,i,mu,nu),2) + (gvguu(i,j,nu,gamma,1)*gfuu(l,i,xi,mu,2) + conjg(gvguu(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvguu(l,i,xi,mu,2) + conjg(gfuu(i,l,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),3) = df1lsoc(sigmaimunu2i(2,i,mu,nu),3) + (gvgud(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gvgud(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),4) = df1lsoc(sigmaimunu2i(2,i,mu,nu),4) + (gvguu(i,j,nu,gamma,1)*gfdu(l,i,xi,mu,2) + conjg(gvgud(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdu(l,i,xi,mu,2) + conjg(gfud(i,l,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1lsoc(sigmaimunu2i(3,i,mu,nu),1) = df1lsoc(sigmaimunu2i(3,i,mu,nu),1) + (gvgdd(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gvgdu(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),2) = df1lsoc(sigmaimunu2i(3,i,mu,nu),2) + (gvgdu(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gvgdu(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),3) = df1lsoc(sigmaimunu2i(3,i,mu,nu),3) + (gvgdd(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gvgdd(i,l,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),4) = df1lsoc(sigmaimunu2i(3,i,mu,nu),4) + (gvgdu(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gvgdd(i,l,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

        df1lsoc(sigmaimunu2i(4,i,mu,nu),1) = df1lsoc(sigmaimunu2i(4,i,mu,nu),1) + (gvgud(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gvgdu(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),2) = df1lsoc(sigmaimunu2i(4,i,mu,nu),2) + (gvguu(i,j,nu,gamma,1)*gfud(l,i,xi,mu,2) + conjg(gvgdu(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgud(l,i,xi,mu,2) + conjg(gfdu(i,l,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),3) = df1lsoc(sigmaimunu2i(4,i,mu,nu),3) + (gvgud(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gvgdd(i,l,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),4) = df1lsoc(sigmaimunu2i(4,i,mu,nu),4) + (gvguu(i,j,nu,gamma,1)*gfdd(l,i,xi,mu,2) + conjg(gvgdd(i,l,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdd(l,i,xi,mu,2) + conjg(gfdd(i,l,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*dtdk(j,l,gamma,xi)

      end do ; end do ; end do ; end do ; end do ; end do ; end do

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

      do nu=1,9 ; do gamma=1,9 ; do mu=1,9 ; do i=1,Npl ; do xi=1,9 ; do l=1,Npl ; do j=1,Npl
        if(abs(j-l) > npln) cycle
        df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3)-zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4)-zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)

        df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3)-zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)
        df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4)-zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))*dtdk(j,l,gamma,xi)


        df1lsoc(sigmaimunu2i(1,i,mu,nu),1) = df1lsoc(sigmaimunu2i(1,i,mu,nu),1)-zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvguu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),2) = df1lsoc(sigmaimunu2i(1,i,mu,nu),2)-zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvguu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),3) = df1lsoc(sigmaimunu2i(1,i,mu,nu),3)-zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgud(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),4) = df1lsoc(sigmaimunu2i(1,i,mu,nu),4)-zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgud(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)

        df1lsoc(sigmaimunu2i(2,i,mu,nu),1) = df1lsoc(sigmaimunu2i(2,i,mu,nu),1)-zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvguu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),2) = df1lsoc(sigmaimunu2i(2,i,mu,nu),2)-zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfuu(i,l,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvguu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),3) = df1lsoc(sigmaimunu2i(2,i,mu,nu),3)-zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgud(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),4) = df1lsoc(sigmaimunu2i(2,i,mu,nu),4)-zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfud(i,l,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgud(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)

        df1lsoc(sigmaimunu2i(3,i,mu,nu),1) = df1lsoc(sigmaimunu2i(3,i,mu,nu),1)-zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),2) = df1lsoc(sigmaimunu2i(3,i,mu,nu),2)-zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),3) = df1lsoc(sigmaimunu2i(3,i,mu,nu),3)-zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdd(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),4) = df1lsoc(sigmaimunu2i(3,i,mu,nu),4)-zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdd(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)

        df1lsoc(sigmaimunu2i(4,i,mu,nu),1) = df1lsoc(sigmaimunu2i(4,i,mu,nu),1)-zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),2) = df1lsoc(sigmaimunu2i(4,i,mu,nu),2)-zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdu(i,l,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdu(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),3) = df1lsoc(sigmaimunu2i(4,i,mu,nu),3)-zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdd(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),4) = df1lsoc(sigmaimunu2i(4,i,mu,nu),4)-zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdd(i,l,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdd(i,l,mu,xi,2)))*dtdk(j,l,gamma,xi)
      end do ; end do ; end do ; end do ; end do ; end do ; end do

    end if

    ! Multiplying the prefactor by the susceptibility and the k-point weight
    ! pfdf1iikl = prefactorlsoc*df1iikl
    wkbzc = cmplx(wkbz(iz),0.d0)
    call zgemm('n','n',dim,4,dim,wkbzc,prefactorlsoc,dim,df1iikl,dim,zero,pfdf1iikl,dim)
    ! pfdf1lsoc = prefactorlsoc*df1iikl + prefactor*df1lsoc
    call zgemm('n','n',dim,4,dim,wkbzc,prefactor,dim,df1lsoc,dim,zum,pfdf1iikl,dim)

    !$omp critical

    ! Integrating matrices

    ! tFintiikl = tFintiikl + pfdf1iikl
#ifdef _JUQUEEN
    call zgeadd(tFintiikl,dim,'N',pfdf1iikl,dim,'N',tFintiikl,dim,dim,4)
#else
    call ZAXPY(dim*4,zum,pfdf1iikl,1,tFintiikl,1)
!       tFintiikl = tFintiikl + pfdf1iikl
#endif
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
  deallocate(df1iikl,pfdf1iikl,df1lsoc)
  deallocate(gvg,gvguu,gvgud,gvgdu,gvgdd)
  deallocate(gf,dtdk,gfuu,gfud,gfdu,gfdd)
!$omp end parallel

  tFintiikl    = tFintiikl/tpi
  ttFintiikl   = ttFintiikl/tpi
  LxttFintiikl = LxttFintiikl/tpi
  LyttFintiikl = LyttFintiikl/tpi
  LzttFintiikl = LzttFintiikl/tpi

  return
end subroutine sumklinearsoc
