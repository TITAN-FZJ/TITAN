! Calculates the orbital angular momentum ground state
subroutine L_gs(s)
  use mod_f90_kind, only: double
  use mod_constants, only: zero,pi
  use mod_System, only: System
  use TightBinding, only: nOrb
  use mod_parameters, only: outputunit_loop, Ef, lverbose, host
  use mod_magnet
  use EnergyIntegration, only: y, wght, pn1
  use mod_mpi_pars
  implicit none
  type(System), intent(inout) :: s
  integer     :: i,mu,nu
  complex(double), dimension(nOrb, nOrb)   :: lx,ly,lz
  complex(double),dimension(s%nAtoms, nOrb, nOrb)  :: gupgd,gupgdint
!--------------------- begin MPI vars --------------------
  integer :: ix,itask
  integer :: ncount
  ncount=s%nAtoms*nOrb*nOrb
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  allocate( lxm(s%nAtoms),lym(s%nAtoms),lzm(s%nAtoms),lxpm(s%nAtoms),lypm(s%nAtoms),lzpm(s%nAtoms) )
  if(myrank_row_hw==0) write(outputunit_loop,"('[L_gs] Calculating Orbital Angular Momentum ground state... ')")


  ix = myrank_row_hw+1
  itask = numprocs ! Number of tasks done initially

!   Calculating the number of particles for each spin and orbital using a complex integral

  gupgdint  = zero
  do while(ix <= pn1)
    call sumk_L_gs(Ef,y(ix),gupgd, s%kbz, s%wkbz, s%nkpt, s%nAtoms)
    gupgdint = gupgdint + gupgd*wght(ix)
    ix = ix + numprocs_row
  end do
  call MPI_Allreduce(MPI_IN_PLACE, gupgdint, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)

  !
  ! if(myrank_row_hw/=0) then
  !   return
  ! end if

  call l_matrix(lx,ly,lz)

  lxpm = 0.d0
  lypm = 0.d0
  lzpm = 0.d0
  lxm = 0.d0
  lym = 0.d0
  lzm = 0.d0
  do nu=5,9 ; do mu=5,9 ; do i=1,s%nAtoms
    lxpm(i) = lxpm(i) + real(lxp(mu,nu)*gupgdint(i,nu,mu))
    lypm(i) = lypm(i) + real(lyp(mu,nu)*gupgdint(i,nu,mu))
    lzpm(i) = lzpm(i) + real(lzp(mu,nu)*gupgdint(i,nu,mu))
    lxm(i)  = lxm(i)  + real(lx (mu,nu)*gupgdint(i,nu,mu))
    lym(i)  = lym(i)  + real(ly (mu,nu)*gupgdint(i,nu,mu))
    lzm(i)  = lzm(i)  + real(lz (mu,nu)*gupgdint(i,nu,mu))
  end do ; end do ; end do

  ! Calculating angles of GS OAM (in units of pi)
  do i = 1,s%nAtoms
    labs(i)   = sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
    ltheta(i) = acos(lzm(i)/sqrt(lxm(i)**2+lym(i)**2+lzm(i)**2))/pi
    lphi(i)   = atan2(lym(i),lxm(i))/pi
    lpabs(i)  = sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2))
    lptheta(i)= acos(lzpm(i)/sqrt(lxpm(i)**2+lypm(i)**2+lzpm(i)**2))/pi
    lpphi(i)  = atan2(lypm(i),lxpm(i))/pi
  end do

  return
end subroutine L_gs



subroutine sumk_L_gs(e,ep,gupgd, kbz, wkbz, nkpt, nAtoms)
  use mod_f90_kind
  use mod_constants, only: pi, zero
!$  use mod_parameters, only: outputunit_loop
  use TightBinding, only: nOrb
!$  use mod_mpi_pars, only: myrank_row_hw
!$  use omp_lib
  implicit none
  integer, intent(in) :: nkpt, nAtoms
  real(double), dimension(3,nkpt), intent(in) :: kbz
  real(double), dimension(nkpt), intent(in) :: wkbz
!$  integer       :: nthreads,mythread
  integer       :: i,mu,nu,mup,nup,iz
  real(double)  :: kp(3)
  real(double),intent(in) :: e,ep
  complex(double), dimension(nAtoms, nOrb, nOrb), intent(out) :: gupgd
  complex(double), dimension(nAtoms,nAtoms,2*nOrb, 2*nOrb) :: gf

  gupgd   = zero

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,mu,nu,mup,nup) &
!$omp& shared(kbz,nkpt,wkbz,e,ep,nAtoms,gupgd,myrank_row_hw,nthreads,outputunit_loop)
!$  mythread = omp_get_thread_num()
!!$  if((mythread==0).and.(myrank_row_hw==0)) then
!!$    nthreads = omp_get_num_threads()
!!$    write(outputunit_loop,"('[L_gs] Number of threads: ',i0)") nthreads
!!$  end if

!$omp do reduction(+:gupgd)
  kpoints: do iz=1,nkpt
    kp = kbz(:,iz)

    !Green function on energy Ef + iy, and wave vector kp
    call green(e,ep,kp,gf)

    do i=1,nAtoms
      do mu=1,nOrb
        mup = mu+nOrb
        do nu=1,nOrb
          nup = nu+nOrb
          gupgd(i,mu,nu) = gupgd(i,mu,nu) + ((gf(i,i,mu,nu)+gf(i,i,mup,nup))*wkbz(iz))
        end do
      end do
    end do
  end do kpoints
!$omp end do
!$omp end parallel

  gupgd = gupgd/pi

  return
end subroutine sumk_L_gs
