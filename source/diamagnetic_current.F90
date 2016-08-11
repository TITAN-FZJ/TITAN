! ---------- Diamagnetic current: Energy integration ---------
subroutine diamagnetic_current()
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_epoints
  use mod_progress
  use mod_mpi_pars
  use MPI
  implicit none
  integer           :: i,j
  real(double),dimension(n0sc1:n0sc2,Npl) :: Idia,Idia_total
!--------------------- begin MPI vars --------------------
  integer :: ix,itask
  integer :: ncount
  ncount=n0sc*Npl
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ix = myrank+1
  itask = numprocs ! Number of tasks done initially

  ! Calculating the number of particles for each spin and orbital using a complex integral
  if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
    call sumk_idia(Ef,y(ix),Idia)
    Idia_total = wght(ix)*Idia

    if(lverbose) write(outputunit,"('[diamagnetic_current] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
    do i=2,pn1
      if(lverbose) call progress_bar(outputunit,"energy points",i,pn1)

      call MPI_Recv(Idia,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,1604,MPI_COMM_WORLD,stat,ierr)

      Idia_total = Idia_total + Idia

      ! If the number of processors is less than the total number of points, sends
      ! the rest of the points to the ones that finish first
      if (itask.lt.pn1) then
        itask = itask + 1
        call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)
      else
        call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
      end if
    end do
  else
    ! Other processors calculate each point of the integral and waits for new points
    do
      if(ix.gt.pn1) exit

      ! First and second integrations (in the complex plane)
      call sumk_idia(Ef,y(ix),Idia)
      Idia = wght(ix)*Idia

!       if(lverbose) write(outputunit,"('[diamagnetic_current] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
      ! Sending results to process 0
      call MPI_Send(Idia,ncount,MPI_DOUBLE_PRECISION,0,1604,MPI_COMM_WORLD,ierr)
      ! Receiving new point or signal to exit
      call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
      if(ix.eq.0) exit
    end do
  end if

  if(myrank.eq.0) then
    do j=n0sc1,n0sc2
      write(outputunit,*) (Idia_total(j,i),i=1,Npl)
    end do
  end if
  return
end subroutine diamagnetic_current

!   !   Calculates k-sum for diamagnetic current
!   subroutine sumk_diamag(e,ep,Idia)
!     use mod_f90_kind
!     use mod_constants
!     use mod_parameters
!     use mod_lattice
!     use mod_generate_kpoints
!     use mod_mpi_pars
! !$  use omp_lib
!     implicit none
! !$  integer       :: nthreads,mythread
!     integer         :: neighbor,mu,nu,iz,i,mup,nup
!     real(double)    :: kp(3)
!     real(double),intent(in)    :: e,ep
!     real(double),dimension(n0sc1:n0sc2,Npl),intent(out) :: Idia
!     complex(double) :: expikr(n0sc1:n0sc2)
!     complex(double),dimension(Npl,Npl,18,18)            :: gf

!     Idia = 0.d0

! !$omp parallel default(none) &
! !$omp& private(mythread,iz,kp,i,mu,nu,mup,nup,gf,expikr) &
! !$omp& shared(kbz,nkpoints,wkbz,Npl,myrank,nthreads,llineargfsoc,n0sc1,n0sc2,Ef,e,ep,r0,pji,Idia)
! !$  mythread = omp_get_thread_num()
! !$  if((mythread.eq.0).and.(myrank.eq.0)) then
! !$    nthreads = omp_get_num_threads()
! !$    write(outputunit,"('[sumk_diamag] Number of threads: ',i0)") nthreads
! !$  end if

! !$omp do reduction(+:Idia)
!     kpoints: do iz=1,nkpoints
!       kp = kbz(iz,:)

!       do neighbor=n0sc1,n0sc2
!         expikr(neighbor) = exp(zi*(kp(1)*r0(neighbor,1)+kp(2)*r0(neighbor,2)+kp(3)*r0(neighbor,3)))
!       end do

!       ! Green function at (k+q,E_F+E+iy)
!       if(llineargfsoc) then
!         call greenlineargfsoc(e,ep,kp,gf)
!       else
!         call green(e,ep,kp,gf)
!       end if

!       do neighbor=n0sc1,n0sc2 ; do i=1,Npl ; do mu=1,9; do nu=1,9
!         mup=mu+9
!         nup=nu+9
!         Idia(neighbor,i) = Idia(neighbor,i) + real((gf(i,i,mu,nu)+gf(i,i,mup,nup))*(pji(neighbor,i,i,nu,mu)*expikr(neighbor)-pji(neighbor,i,i,mu,nu)*conjg(expikr(neighbor))))
!       end do ; end do ; end do ; end do

!     end do kpoints
! !$omp end do
! !$omp end parallel

!     return
!   end subroutine sumk_diamag

!   Calculates the momentum matrix in real space
subroutine sumk_idia(e,ep,Idia)
  use mod_f90_kind
  use mod_parameters, only: Npl,n0sc1,n0sc2,llineargfsoc
  use mod_lattice
  use mod_generate_kpoints
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer       :: nthreads,mythread
  integer         :: neighbor,iz,i,mu,nu,mup,nup
  real(double)    :: kp(3)
  complex(double) :: expikr(n0sc1:n0sc2)
  real(double),intent(in)    :: e,ep
  real(double),dimension(n0sc1:n0sc2,Npl),intent(out) :: Idia
  complex(double),dimension(Npl,Npl,18,18)            :: gf
  complex(double),dimension(Npl,Npl,9,9)              :: dtdk
  complex(double),dimension(n0sc1:n0sc2,Npl,9,9)      :: pji
  complex(double),dimension(n0sc1:n0sc2,Npl,18,18)    :: gij,gji

  pji = zero
  gji = zero
  gij = zero

!$omp parallel default(none) &
!$omp& private(mythread,neighbor,iz,kp,i,dtdk,expikr,gf) &
!$omp& shared(kbz,nkpoints,wkbz,myrank,nthreads,n0sc1,n0sc2,llineargfsoc,Npl,r0,e,ep,pji,gij,gji,outputunit)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit,"('[sumk_idia] Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:pji,gij,gji)
  kpoints: do iz=1,nkpoints
    kp = kbz(iz,:)

    do neighbor=n0sc1,n0sc2
      expikr(neighbor) = exp(-zi*(kp(1)*r0(neighbor,1)+kp(2)*r0(neighbor,2)+kp(3)*r0(neighbor,3)))
    end do

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(kp,dtdk)

    ! Green function at (k+q,E_F+E+iy)
    if(llineargfsoc) then
      call greenlineargfsoc(e,ep,kp,gf)
    else
      call green(e,ep,kp,gf)
    end if

    do neighbor=n0sc1,n0sc2 ; do i=1,Npl
      pji(neighbor,i,:,:) = pji(neighbor,i,:,:) + expikr(neighbor)*dtdk(i,i,:,:)*wkbz(iz)
      gij(neighbor,i,:,:) = gij(neighbor,i,:,:) + conjg(expikr(neighbor))*gf(i,i,:,:)*wkbz(iz)
      gji(neighbor,i,:,:) = gji(neighbor,i,:,:) + expikr(neighbor)*gf(i,i,:,:)*wkbz(iz)
    end do ; end do

!       !$omp critical
!       do neighbor=n0sc1,n0sc2 ; do i=1,Npl ; do mu=1,9 ; do nu=1,9
!         pji(neighbor,i,mu,nu) = pji(neighbor,i,mu,nu) + expikr(neighbor)*dtdk(i,i,mu,nu)*wkbz(iz)
!       end do ; end do ; end do ; end do
!       do neighbor=n0sc1,n0sc2 ; do i=1,Npl ; do mu=1,18 ; do nu=1,18
!         gij(neighbor,i,mu,nu) = gij(neighbor,i,mu,nu) + conjg(expikr(neighbor))*gf(i,i,mu,nu)*wkbz(iz)
!         gji(neighbor,i,mu,nu) = gji(neighbor,i,mu,nu) + expikr(neighbor)*gf(i,i,mu,nu)*wkbz(iz)
!       end do ; end do ; end do ; end do
!       !$omp end critical

  end do kpoints
!$omp end do
!$omp end parallel

!     if(myrank.eq.0) then
!       write(outputunit,*) "pij"
!       do neighbor=n0sc1,n0sc2
!         write(outputunit,*) (sum(pji(neighbor,iz,:,:)),iz=1,Npl)
!       end do

!       write(outputunit,*) "gij"
!       do neighbor=n0sc1,n0sc2
!         write(outputunit,*) (sum(gij(neighbor,iz,:,:)),iz=1,Npl)
!       end do

!       write(outputunit,*) "gji"
!       do neighbor=n0sc1,n0sc2
!         write(outputunit,*) (sum(gij(neighbor,iz,:,:)),iz=1,Npl)
!       end do
!     end if
!     call MPI_Finalize(ierr)
!     stop

  Idia = 0.d0
  do neighbor=n0sc1,n0sc2 ; do i=1,Npl ; do mu=1,9 ; do nu=1,9
    mup = mu+9
    nup = nu+9
    Idia(neighbor,i) = Idia(neighbor,i) + zi*pji(neighbor,i,mu,nu)*aimag(gij(neighbor,i,nu,mu)+gji(neighbor,i,mu,nu)+gij(neighbor,i,nup,mup)+gji(neighbor,i,mup,nup))
  end do ; end do ; end do ; end do

  return
end subroutine sumk_idia