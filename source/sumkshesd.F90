! ----------- Sum over wave vectors to calculate spin disturbance -----------
subroutine sumkshesd(e,ep,Fint,tFintiikl,iflag)
  use mod_f90_kind
  use mod_parameters
  use mod_constants
  use mod_lattice, only: plnn
  use mod_generate_kpoints
  use mod_progress
  use mod_mpi_pars
  use MPI
!$  use omp_lib
  implicit none
!$  integer         :: nthreads,mythread
  integer         :: AllocateStatus
  integer         :: i,j,l,mu,nu,gamma,xi,sigma,sigmap,iz
  integer, intent(in) :: iflag
  real(double),intent(in)     :: e,ep
  real(double)                :: kp(3)
  complex(double),dimension(:,:,:,:),allocatable    :: gf,dtdk
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfmat,gfuu,gfud,gfdu,gfdd
  complex(double),dimension(dim,dim),intent(out)    :: Fint
  complex(double),dimension(dim,dimNpl),intent(out) :: tFintiikl
  complex(double),dimension(:,:),allocatable        :: df1
  complex(double),dimension(:,:),allocatable        :: df1iikl

  Fint      = zero
  tFintiikl = zero

!$omp parallel default(none) &
!$omp& private(errorcode,ierr,mythread,AllocateStatus,iz,kp,df1iikl,dtdk,gf,gfmat,l,gfuu,gfud,gfdu,gfdd,sigma,sigmap,i,j,mu,nu,gamma,xi,df1) &
!$omp& shared(prog,runoptions,kbz,wkbz,e,ep,iflag,Fint,tFintiikl,nkpoints,Ef,eta,nthreads,myrank,Npl,plnn,dim,dimNpl,sigmaimunu2i,sigmaijmunu2i)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if
  allocate(df1(dim,dim),df1iikl(dim,dimNpl),gf(Npl,Npl,18,18),dtdk(Npl,Npl,9,9),gfuu(2,Npl,Npl,9,9),gfud(2,Npl,Npl,9,9),gfdu(2,Npl,Npl,9,9),gfdd(2,Npl,Npl,9,9),gfmat(2,Npl,Npl,18,18), STAT = AllocateStatus  )
  if (AllocateStatus.ne.0) then
    write(*,"('[sumkshesd] Not enough memory for: df1,df1iikl,gf,dtdk,gfuu,gfud,gfdu,gfdd,gfmat')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do schedule(static)
  do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if((myrank.eq.0).and.(index(runoptions,"verbose").gt.0)) then
        prog = floor(iz*100.d0/nkpoints)
        progress_bar: select case (mod(iz,4))
        case(0)
          write(*,"(a1,2x,i3,'% of k-sum on rank 0',a1,$)") '|',prog,char(13)
        case(1)
          write(*,"(a1,2x,i3,'% of k-sum on rank 0',a1,$)") '/',prog,char(13)
        case(2)
          write(*,"(a1,2x,i3,'% of k-sum on rank 0',a1,$)") '-',prog,char(13)
        case(3)
          write(*,"(a1,2x,i3,'% of k-sum on rank 0',a1,$)") '\',prog,char(13)
        end select progress_bar
      end if
!$  end if

    kp = kbz(iz,:)

    if(iflag.eq.0)then
      df1 = zero
      df1iikl = zero
      gfmat = zero

      ! Calculating derivative of in-plane and n.n. inter-plane hoppings
      call dtdksub(kp,dtdk)

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
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = gfdu(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = gfdu(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = gfdd(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = gfdd(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = gfuu(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = gfuu(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = gfud(1,i,j,nu,gamma)*gfdu(2,j,i,xi,mu) + conjg(gfud(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = gfud(1,i,j,nu,gamma)*gfuu(2,j,i,xi,mu) + conjg(gfuu(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = gfdu(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = gfdu(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfud(1,j,i,gamma,nu))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = gfdd(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = gfdd(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfdd(1,j,i,gamma,nu))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = gfuu(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = gfuu(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfuu(1,j,i,gamma,nu))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = gfud(1,i,j,nu,gamma)*gfdd(2,j,i,xi,mu) + conjg(gfdd(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = gfud(1,i,j,nu,gamma)*gfud(2,j,i,xi,mu) + conjg(gfdu(2,i,j,mu,xi)*gfdu(1,j,i,gamma,nu))
      end do ; end do ; end do ; end do ; end do ; end do

      do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        if(abs(j-l).gt.plnn) cycle
        ! Diagonal values in j,l are the same calculated for chi_ij
        if(l.eq.j) then
          do sigmap=1,4 ; do sigma=1,4
            df1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,l,gamma,xi)) = df1(sigmaimunu2i(sigma,i,mu,nu),sigmaimunu2i(sigmap,j,gamma,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          end do ; end do
        else
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)

          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfdu(2,l,i,xi,mu) + conjg(gfud(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfuu(2,l,i,xi,mu) + conjg(gfuu(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)

          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = (gfdu(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfud(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = (gfdd(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfdd(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)

          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = (gfuu(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfuu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfdd(2,l,i,xi,mu) + conjg(gfdd(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = (gfud(1,i,j,nu,gamma)*gfud(2,l,i,xi,mu) + conjg(gfdu(2,i,l,mu,xi)*gfdu(1,j,i,gamma,nu)))*dtdk(j,l,gamma,xi)*wkbz(iz)
        end if
      end do ; end do ; end do ; end do ; end do ; end do ; end do

      !$omp critical
      Fint      = Fint + (df1*wkbz(iz))
      tFintiikl = tFintiikl + df1iikl
!       do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl ; do sigmap=1,4 ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4
!         if(abs(j-l).gt.plnn) cycle
!         tFintiikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,l,gamma,xi)) = tFintiikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,l,gamma,xi)) + df1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,l,gamma,xi))
!       end do ; end do ; end do ; end do ; end do ; end do ; end do ; end do ; end do
      !$omp end critical
    else
      df1 = zero
      df1iikl = zero
      gfmat = zero

      ! Calculating derivative of in-plane and n.n. inter-plane hoppings
      call dtdksub(kp,dtdk)

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
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfud(2,i,j,mu,xi))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,j,mu,xi))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,j,mu,xi))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,j,mu,xi))
      end do ; end do ; end do ; end do ; end do ; end do

      do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl
        if(abs(j-l).gt.plnn) cycle
        ! Diagonal values in j,l are the same calculated for chi_ij
        if(l.eq.j) then
          do sigmap=1,4 ; do sigma=1,4
            df1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,j,gamma,xi)) = -zi*df1(sigmaimunu2i(sigma,i,mu,nu),sigmaimunu2i(sigmap,j,gamma,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          end do ; end do
        else
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(1,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)

          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfud(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(2,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfuu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)

          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = -zi*(gfdu(1,i,j,nu,gamma)-conjg(gfud(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(3,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = -zi*(gfdd(1,i,j,nu,gamma)-conjg(gfdd(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)

          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(1,j,l,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(2,j,l,gamma,xi)) = -zi*(gfuu(1,i,j,nu,gamma)-conjg(gfuu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(3,j,l,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdd(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
          df1iikl(sigmaimunu2i(4,i,mu,nu),sigmaijmunu2i(4,j,l,gamma,xi)) = -zi*(gfud(1,i,j,nu,gamma)-conjg(gfdu(1,j,i,gamma,nu)))*conjg(gfdu(2,i,l,mu,xi))*dtdk(j,l,gamma,xi)*wkbz(iz)
        end if
      end do ; end do ; end do ; end do ; end do ; end do ; end do

      !$omp critical
      Fint     = Fint - (zi*df1*wkbz(iz))
      tFintiikl = tFintiikl + df1iikl
!       do xi=1,9 ; do gamma=1,9 ; do l=1,Npl ; do j=1,Npl ; do sigmap=1,4 ; do nu=1,9 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4
!         if(abs(j-l).gt.plnn) cycle
!         tFintiikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,l,gamma,xi)) = tFintiikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,l,gamma,xi)) + df1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmaijmunu2i(sigmap,j,l,gamma,xi))
!       end do ; end do ; end do ; end do ; end do ; end do ; end do ; end do ; end do
      !$omp end critical
    end if
  end do
!$omp end do
  deallocate(df1,df1iikl)
  deallocate(gf,dtdk,gfuu,gfud,gfdu,gfdd,gfmat)
!$omp end parallel

  Fint     = Fint/tpi
  tFintiikl = tFintiikl/tpi

  return
end subroutine sumkshesd