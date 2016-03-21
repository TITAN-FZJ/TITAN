! Integration of Green functions over k vectors to calculate
! the jacobian of self-consistency system
subroutine sumk_selfconjac(er,ei,ggr)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_mpi_pars
  use MPI
!$  use omp_lib
  implicit none
!$  integer                :: nthreads,mythread
  integer                  :: AllocateStatus
  integer                  :: iz,i,j,i0,j0,mu,nu,sigma,sigmap
  real(double)             :: kp(3)
  real(double),intent(in)  :: er,ei
  real(double),dimension(4*Npl,4*Npl),intent(out) :: ggr
  complex(double)                                 :: mhalfU(4,Npl),wkbzc
  complex(double),dimension(4,18,18)              :: pauli,paulid,temp1,temp2
  complex(double),dimension(18,18)                :: gij,gji,temp,paulitemp
  complex(double),dimension(:,:,:,:,:,:),allocatable    :: gdHdxg,gvgdHdxgvg
  complex(double),dimension(Npl,Npl,18,18)        :: gf,gvg

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

! Pauli matrices in spin and orbital space
  pauli  = zero
  paulid = zero
  do mu = 1,9
    nu = mu+9
    ! identity
    pauli(1,mu,mu) = zum
    pauli(1,nu,nu) = zum
    if (mu.lt.5) cycle     ! Pauli matrices for d orbitals only
    ! paulid matrix x
    pauli(2,mu,nu) = zum
    pauli(2,nu,mu) = zum
    ! paulid matrix y
    pauli(3,mu,nu) = -zi
    pauli(3,nu,mu) = zi
    ! paulid matrix z
    pauli(4,mu,mu) = zum
    pauli(4,nu,nu) = -zum

    ! identity
    paulid(1,mu,mu) = zum
    paulid(1,nu,nu) = zum
  end do
  paulid(2:4,:,:) = pauli(2:4,:,:)

! Prefactor -U/2 in dH/dm and 1 in dH/deps1
  do j=1,Npl
    mhalfU(1,j) = zum
    mhalfU(2:4,j) = -0.5d0*U(j+1)
  end do

  ggr    = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,AllocateStatus,errorcode,ierr,iz,kp,wkbzc,gf,gvg,temp,temp1,temp2,gij,gji,paulitemp,gdHdxg,gvgdHdxgvg,i,j,i0,j0,mu,nu,sigma) &
!$omp& shared(llineargfsoc,llinearsoc,prog,spiner,elapsed_time,start_program,progbar,lverbose,kbz,wkbz,nkpoints,er,ei,ggr,mhalfU,pauli,paulid,Npl,myrank,nthreads)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if
  allocate( gdHdxg(4,4,Npl,Npl,18,18),gvgdHdxgvg(4,4,Npl,Npl,18,18) , STAT = AllocateStatus  )
  if (AllocateStatus.ne.0) then
    write(*,"('[sumk_selfconjac] Not enough memory for: gdHdxg,gvgdHdxgvg')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!$omp do reduction(+:ggr)
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
    wkbzc = cmplx(wkbz(iz),0.d0)

    ! Green function on energy Ef + iy, and wave vector kp
    if((llineargfsoc).or.(llinearsoc)) then
      call greenlinearsoc(er,ei,kp,gf,gvg)
      gf = gf + gvg
    else
      call green(er,ei,kp,gf)
    end if

    do j=1,Npl ; do i=1,Npl
      gij = gf(i,j,:,:)
      gji = gf(j,i,:,:)

      do sigma = 1,4
        ! temp1 =  pauli*g_ij
        paulitemp = pauli(sigma,:,:)
        call zgemm('n','n',18,18,18,zum,paulitemp,18,gij,18,zero,temp,18)
        temp1(sigma,:,:) = temp
      end do

      do sigmap = 1,4
        ! temp2 = (-U/2) * sigma* g_ji
        paulitemp = paulid(sigmap,:,:)
        call zgemm('n','n',18,18,18,mhalfU(sigmap,j),paulitemp,18,gji,18,zero,temp,18)
        temp2(sigmap,:,:) = temp
      end do

      do sigmap = 1,4 ; do sigma = 1,4
        ! gdHdxg = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji
        gij = temp1(sigma,:,:)
        gji = temp2(sigmap,:,:)
        call zgemm('n','n',18,18,18,wkbzc,gij,18,gji,18,zero,temp,18)
        gdHdxg(sigma,sigmap,i,j,:,:) = temp
      end do ; end do


      if((llineargfsoc).or.(llinearsoc)) then ! non-linear term
        gij = gvg(i,j,:,:)
        gji = gvg(j,i,:,:)

        do sigma = 1,4
          ! temp1 = wkbz* pauli*gvg_ij
          paulitemp = pauli(sigma,:,:)
          call zgemm('n','n',18,18,18,zum,paulitemp,18,gij,18,zero,temp,18)
          temp1(sigma,:,:) = temp
        end do

        do sigmap = 1,4
          ! temp2 = (-U/2) * sigma* gvg_ji
          paulitemp = paulid(sigmap,:,:)
          call zgemm('n','n',18,18,18,mhalfU(sigmap,j),paulitemp,18,gji,18,zero,temp,18)
          temp2(sigmap,:,:) = temp
        end do

        do sigmap = 1,4 ; do sigma = 1,4
          ! gdHdxg = temp1*temp2 = wkbz* pauli*gvg_ij*(-U/2)*sigma* gvg_ji
          gij = temp1(sigma,:,:)
          gji = temp2(sigmap,:,:)
          call zgemm('n','n',18,18,18,wkbzc,gij,18,gji,18,zero,temp,18)
          gvgdHdxgvg(sigma,sigmap,i,j,:,:) = temp
        end do ; end do
      end if
    end do ; end do

    ! removing non-linear SOC term
    if((llineargfsoc).or.(llinearsoc)) gdHdxg = gdHdxg - gvgdHdxgvg

!    !$omp critical
    do mu=1,18 ; do j=1,Npl ; do i=1,Npl ; do sigmap=1,4 ; do sigma=1,4
      i0 = (sigma-1)*Npl + i
      j0 = (sigmap-1)*Npl + j
      ! Trace over orbitals and spins of the real part
      ggr(i0,j0) = ggr(i0,j0) + real(gdHdxg(sigma,sigmap,i,j,mu,mu))
    end do ; end do ; end do ; end do ; end do
!    !$omp end critical
  end do kpoints
!$omp end do
  deallocate(gdHdxg,gvgdHdxgvg)
!$omp end parallel

  return
end subroutine sumk_selfconjac