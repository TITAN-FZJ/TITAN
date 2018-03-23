module mod_diamagnetic_current
  use mod_f90_kind, only: double
  implicit none

  ! Diamagnetic current for each neighbor and each plane
  real(double),dimension(:,:), allocatable :: Idia_total
contains
  ! This subroutine allocates variables related to the diamagnetic current
  subroutine allocate_idia()
    use mod_mpi_pars,   only: myrank
    use mod_parameters, only: Npl
    use mod_system,     only: n0sc1, n0sc2
    implicit none

    if(myrank==0) allocate( Idia_total(n0sc1:n0sc2,Npl) )

      end subroutine allocate_idia

  ! This subroutine deallocates variables related to the diamagnetic current
  subroutine deallocate_idia()
    use mod_mpi_pars, only: myrank
    implicit none
    if(myrank==0) deallocate( Idia_total )
      end subroutine deallocate_idia

  ! ---------- Diamagnetic current: Energy integration ---------
  subroutine calculate_idia()
    use mod_f90_kind,         only: double
    use MPI,                  only: MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SOURCE, MPI_ANY_TAG, MPI_ANY_SOURCE, MPI_COMM_WORLD
    use mod_mpi_pars,         only: myrank, numprocs, stat, ierr
    use mod_generate_epoints, only: y, wght
    use mod_progress,         only: progress_bar
    use mod_system,           only: n0sc1, n0sc2, n0sc, s => sys
    use mod_parameters

    implicit none

    integer           :: i,j
    real(double),dimension(n0sc1:n0sc2,Npl) :: Idia
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount
    ncount=n0sc*Npl
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    if(myrank==0) write(outputunit,"('[diamagnetic_current] Calculating diamagnetic current... ')", advance='no')
    ix = myrank+1
    itask = numprocs ! Number of tasks done initially

    ! Calculating the number of particles for each spin and orbital using a complex integral
    if (myrank==0) then ! Process 0 receives all results and send new tasks if necessary
      call sumk_idia(s%Ef,y(ix)+eta,Idia)
      Idia_total = wght(ix)*Idia

      if(lverbose) write(outputunit,"(' Finished point ',i0,' in rank ',i0)") ix,myrank
      do i=2,pn1
        if(lverbose) call progress_bar(outputunit,"energy points",i,pn1)

        call MPI_Recv(Idia,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,1604,MPI_COMM_WORLD,stat,ierr)

        Idia_total = Idia_total + Idia

        ! If the number of processors is less than the total number of points, sends
        ! the rest of the points to the ones that finish first
        if (itask<pn1) then
          itask = itask + 1
          call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)
        else
          call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
        end if
      end do
    else
      ! Other processors calculate each point of the integral and waits for new points
      do
        if(ix>pn1) exit

        ! First and second integrations (in the complex plane)
        call sumk_idia(s%Ef,y(ix),Idia)
        Idia = wght(ix)*Idia

  !       if(lverbose) write(outputunit,"('[diamagnetic_current] Finished point ',i0,' in rank ',i0)") ix,myrank
        ! Sending results to process 0
        call MPI_Send(Idia,ncount,MPI_DOUBLE_PRECISION,0,1604,MPI_COMM_WORLD,ierr)
        ! Receiving new point or signal to exit
        call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
        if(ix==0) exit
      end do
    end if

    if(myrank==0) then
      write(outputunit,"('Finished!')")
      do i=1,Npl
        write(outputunit,fmt="(10(es11.4,2x))") (Idia_total(j,i),j=n0sc1,n0sc2)
      end do
    end if
      end subroutine calculate_idia

  !   Calculates the momentum matrix in real space
  subroutine sumk_idia(e,ep,Idia)
    use mod_f90_kind
    use mod_constants
    use TightBinding
    use mod_parameters, only: Npl,llineargfsoc,outputunit
    use mod_system,     only: n0sc1,n0sc2,r_nn,nkpt,kbz,wkbz, s=>sys
    !use mod_generate_kpoints
    use mod_mpi_pars
!$  use omp_lib
    implicit none
!$  integer       :: nthreads,mythread
    integer         :: neighbor,iz,i,mu,nu,mup,nup
    real(double)    :: kp(3)
    complex(double) :: expikr(n0sc1:n0sc2)
    real(double),intent(in)    :: e,ep
    real(double),dimension(n0sc1:n0sc2,Npl),intent(out) :: Idia
    complex(double),dimension(Npl,Npl,nOrb2,nOrb2)            :: gf
    complex(double),dimension(nOrb,nOrb,Npl,Npl)              :: dtdk
    complex(double),dimension(n0sc1:n0sc2,Npl,nOrb,nOrb)      :: pij
    complex(double),dimension(n0sc1:n0sc2,Npl,nOrb2,nOrb2)    :: gij,gji

    pij = cZero
    gji = cZero
    gij = cZero

!$omp parallel default(none) &
!$omp& private(mythread,neighbor,iz,kp,i,dtdk,expikr,gf) &
!$omp& shared(kbz,nkpt,wkbz,myrank,nthreads,n0sc1,n0sc2,llineargfsoc,Npl,r_nn,e,ep,pij,gij,gji, outputunit)
!$  mythread = omp_get_thread_num()
!$  if((mythread==0).and.(myrank==0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit,"('[sumk_idia] Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:pij,gij,gji)
    kpoints: do iz=1,nkpt
      kp = kbz(:,iz)

      do neighbor=n0sc1,n0sc2
        expikr(neighbor) = exp(-cI*dot_product(kp, r_nn(:,neighbor)))
      end do

      ! Calculating derivative of in-plane and n.n. inter-plane hoppings
      call dtdksub(kp,dtdk)

      ! Green function at (k+q,E_F+E+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(e,ep,s,kp,gf)
      else
        call green(e,ep,s,kp,gf)
      end if

      do neighbor=n0sc1,n0sc2
        do i=1,Npl
          pij(neighbor,i,:,:) = pij(neighbor,i,:,:) + expikr(neighbor)*dtdk(:,:,i,i)*wkbz(iz)
          gij(neighbor,i,:,:) = gij(neighbor,i,:,:) + expikr(neighbor)*gf(:,:,i,i)*wkbz(iz)
          gji(neighbor,i,:,:) = gji(neighbor,i,:,:) + conjg(expikr(neighbor))*gf(:,:,i,i)*wkbz(iz)
        end do
      end do

  !       !$omp critical
  !       do neighbor=n0sc1,n0sc2 ; do i=1,Npl ; do mu=1,9 ; do nu=1,9
  !         pij(neighbor,i,mu,nu) = pij(neighbor,i,mu,nu) + expikr(neighbor)*dtdk(mu,nu,i,i)*wkbz(iz)
  !       end do ; end do ; end do ; end do
  !       do neighbor=n0sc1,n0sc2 ; do i=1,Npl ; do mu=1,18 ; do nu=1,18
  !         gij(neighbor,i,mu,nu) = gij(neighbor,i,mu,nu) + conjg(expikr(neighbor))*gf(mu,nu,i,i)*wkbz(iz)
  !         gji(neighbor,i,mu,nu) = gji(neighbor,i,mu,nu) + expikr(neighbor)*gf(mu,nu,i,i)*wkbz(iz)
  !       end do ; end do ; end do ; end do
  !       !$omp end critical

    end do kpoints
!$omp end do
!$omp end parallel

      ! if(myrank==0) then
      !   write(outputunit,*) "pij"
      !   do neighbor=n0sc1,n0sc2
      !     write(outputunit,*) (sum(pij(neighbor,iz,:,:)),iz=1,Npl)
      !   end do

      !   write(outputunit,*) "gij"
      !   do neighbor=n0sc1,n0sc2
      !     write(outputunit,*) (sum(gij(neighbor,iz,:,:)),iz=1,Npl)
      !   end do

      !   write(outputunit,*) "gji"
      !   do neighbor=n0sc1,n0sc2
      !     write(outputunit,*) (sum(gij(neighbor,iz,:,:)),iz=1,Npl)
      !   end do
      ! end if
      ! call MPI_Finalize(ierr)
      ! stop

    Idia = 0.d0
    do neighbor=n0sc1,n0sc2
      do i=1,Npl
        do mu=1,9
          do nu=1,9
            mup = mu+9
            nup = nu+9
            ! Idia(neighbor,i) = Idia(neighbor,i) + cI*pij(neighbor,i,mu,nu)*aimag(gij(neighbor,i,nu,mu)+gji(neighbor,i,mu,nu)+gij(neighbor,i,nup,mup)+gji(neighbor,i,mup,nup))
            Idia(neighbor,i) = Idia(neighbor,i) + real(pij(neighbor,i,mu,nu)*( (gji(neighbor,i,nu,mu)+gij(neighbor,i,mu,nu)) + (gji(neighbor,i,nup,mup)+gij(neighbor,i,mup,nup)) ))/pi
          end do
        end do
      end do
    end do

      end subroutine sumk_idia
end module mod_diamagnetic_current
