! ---------- Spin disturbance: Energy integration ---------
subroutine eintshechi(e, count)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, zum, zi, tpi
  use mod_parameters, only: eta, ef, dim, sigmaijmunu2i, sigmaimunu2i
  use EnergyIntegration, only: generate_real_epoints, y, wght, x2, p2, nepoints, pn1
  use mod_susceptibilities, only: chiorb_hf
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use mod_SOC, only: llineargfsoc
  use mod_mpi_pars
  !$  use omp_lib
  implicit none
  real(double), intent(in)    :: e
  integer, intent(in) :: count

  integer :: AllocateStatus
  complex(double), dimension(:,:),allocatable :: Fint



  integer         :: i,j,mu,nu,gamma,xi,iz
  real(double)                :: kp(3)
  complex(double),dimension(:,:,:,:),allocatable    :: gf
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd
  complex(double),dimension(:,:),allocatable        :: df1



!--------------------- begin MPI vars --------------------
  integer :: ix,ix2
  integer :: start, start1, start2, end, end1, end2, work, remainder
  integer :: ncount
  ncount=dim*dim
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  ! Calculate workload for each MPI process
  remainder = mod(nepoints,numprocs_row)
  if(myrank_row < remainder) then
    work = ceiling(dble(nepoints) / dble(numprocs_row))
    start = myrank_row*work + 1
    end = (myrank_row+1) * work
  else
    work = floor(dble(nepoints) / dble(numprocs_row))
    start = myrank_row*work + 1 + remainder
    end = (myrank_row+1) * work + remainder
  end if

  start1 = 1
  end1 = pn1
  start2 = pn1 + 1
  end2 = nepoints
  if(start <= pn1) then
    start1 = start
  else
    start1 = pn1 + 1
    start2 = start
  end if

  if(end <= pn1) then
    end1 = end
    end2 = 0
  else
    end2 = end
  end if
  start2 = start2 - pn1
  end2 = end2 - pn1

  ! Starting to calculate energy integral
  chiorb_hf = zero

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,iz,ix,ix2,i,j,mu,nu,gamma,xi,kp,gf,gfuu,gfud,gfdu,gfdd,df1,Fint) &
  !$omp& shared(llineargfsoc,start1,start2,end1,end2,s,e,y,x2,wght,p2,Ef,eta,dim,sigmaimunu2i,sigmaijmunu2i,chiorb_hf)
  allocate(df1(dim,dim), Fint(dim,dim), &
           gf  (s%nAtoms,s%nAtoms,2*nOrb,2*nOrb), &
           gfuu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
           gfud(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
           gfdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
           gfdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[eintshechi] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")
  Fint = zero

  !$omp do schedule(static) collapse(2)
  do ix = start1, end1
    do iz = 1, s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E_F+E+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef+e,y(ix),kp,gf)
      else
        call green(Ef+e,y(ix),kp,gf)
      end if
      gfuu(:,:,:,:,1) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,1) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E_F+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef,y(ix),kp,gf)
      else
        call green(Ef,y(ix),kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,2) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))
      end do ; end do ; end do ; end do ; end do ; end do

      Fint = Fint + df1*wght(ix)*s%wkbz(iz)
    end do !end nkpt loop
  end do !end ix <= pn1 loop
  !$omp end do nowait

  !$omp do schedule(static) collapse(2)
  do ix2 = start2, end2 ! Third integration (on the real axis)
    do iz = 1, s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E'+E+i.eta)
      if(llineargfsoc) then
        call greenlineargfsoc(x2(ix2)+e,eta,kp,gf)
      else
        call green(x2(ix2)+e,eta,kp,gf)
      end if
      gfuu(:,:,:,:,1) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,1) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E'+i.eta)
      if(llineargfsoc) then
        call greenlineargfsoc(x2(ix2),eta,kp,gf)
      else
        call green(x2(ix2),eta,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,2) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))
      end do ; end do ; end do ; end do ; end do ; end do

      ! Locally add up df1
      Fint = Fint + df1*p2(ix2)*s%wkbz(iz)
    end do !end nkpt loop
  end do !end pn1+1, nepoints loop
  !$omp end do nowait

  Fint = Fint / tpi
  !$omp critical
    chiorb_hf = chiorb_hf + Fint
  !$omp end critical

  deallocate(df1, Fint)
  deallocate(gf,gfuu,gfud,gfdu,gfdd)
  !$omp end parallel
  if(myrank_row == 0) then
    call MPI_Reduce(MPI_IN_PLACE, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
  else
    call MPI_Reduce(chiorb_hf, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
  end if
  return
end subroutine eintshechi

! -------------------- Spin disturbance: Energy integration --------------------
! -------------- to be used in the calculation of linear SOC chi ---------------
subroutine eintshechilinearsoc(e, count)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, zum, zi, tpi
  use mod_parameters, only: eta, ef, dim, sigmaijmunu2i, sigmaimunu2i
  use EnergyIntegration, only: generate_real_epoints,y, wght, x2, p2, nepoints, pn1
  use mod_susceptibilities, only: chiorb_hf,chiorb_hflsoc
  use mod_mpi_pars
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use mod_SOC, only: llineargfsoc
  !$  use omp_lib

  implicit none
  integer, intent(in) :: count
  real(double), intent(in) :: e

  integer :: AllocateStatus
  integer :: ix,ix2, iz
  integer :: i,j,mu,nu,gamma,xi
  real(double) :: kp(3)
  complex(double), dimension(:,:,:,:), allocatable :: gf,gvg
  complex(double), dimension(:,:,:,:,:), allocatable :: gfuu,gfud,gfdu,gfdd
  complex(double), dimension(:,:,:,:,:), allocatable :: gvguu,gvgud,gvgdu,gvgdd
  complex(double), dimension(:,:), allocatable :: Fint,Fintlsoc
  complex(double), dimension(:,:), allocatable :: df1,df1lsoc

  !--------------------- begin MPI vars --------------------
  integer :: start, start1, start2, end, end1, end2, work, remainder
  integer :: ncount
  ncount=dim*dim
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  ! Calculate workload for each MPI process
  remainder = mod(nepoints,numprocs_row)
  if(myrank_row < remainder) then
    work = ceiling(dble(nepoints) / dble(numprocs_row))
    start = myrank_row*work + 1
    end = (myrank_row+1) * work
  else
    work = floor(dble(nepoints) / dble(numprocs_row))
    start = myrank_row*work + 1 + remainder
    end = (myrank_row+1) * work + remainder
  end if

  start1 = 1
  end1 = pn1
  start2 = pn1 + 1
  end2 = nepoints
  if(start <= pn1) then
    start1 = start
  else
    start1 = pn1 + 1
    start2 = start
  end if

  if(end <= pn1) then
    end1 = end
    end2 = 0
  else
    end2 = end
  end if
  start2 = start2 - pn1
  end2 = end2 - pn1

  !------------------------------------------------------

  ! Starting to calculate energy integral
  chiorb_hf     = zero
  chiorb_hflsoc = zero

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,ix2,iz,i,j,mu,nu,gamma,xi,kp,Fint,Fintlsoc,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,df1,df1lsoc) &
  !$omp& shared(llineargfsoc,start1,end1,start2,end2,s,e,y,wght,x2,p2,Ef,eta,dim,sigmaimunu2i,sigmaijmunu2i,chiorb_hf,chiorb_hflsoc)
  allocate(df1(dim,dim), Fint(dim,dim), &
           gf(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb), &
           gfuu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
           gfud(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
           gfdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
           gfdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[eintshechilinearsoc] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")

  allocate( df1lsoc(dim,dim), Fintlsoc(dim,dim), &
            gvg(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb), &
            gvguu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
            gvgud(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
            gvgdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
            gvgdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[eintshechilinearsoc] Not enough memory for: df1lsoc,Fintsloc,gvg,gvguu,gvgud,gvgdu,gvgdd")

  Fint      = zero
  Fintlsoc  = zero

  ! Starting to calculate energy integral
  !$omp do schedule(static) collapse(2)
  do ix = start1, end1
    do iz=1,s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E_F+E+iy)
      call greenlinearsoc(Ef+e,y(ix),kp,gf,gvg)
      gfuu(:,:,:,:,1) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,1) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,1) = gvg(:,:, 1: nOrb, 1: nOrb)
      gvgud(:,:,:,:,1) = gvg(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gvgdd(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(Ef,y(ix),kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(:,:,     1:  nOrb,     1:  nOrb)
      gfud(:,:,:,:,2) = gf(:,:,     1:  nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,     1:  nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,2) = gvg(:,:,     1:  nOrb,     1:  nOrb)
      gvgud(:,:,:,:,2) = gvg(:,:,     1:  nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb,     1:  nOrb)
      gvgdd(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))
      end do ; end do ; end do ; end do ; end do ; end do

      Fint = Fint + df1 * s%wkbz(iz) * wght(ix)
      Fintlsoc = Fintlsoc + df1lsoc * s%wkbz(iz) * wght(ix)

    end do
  end do
  !$omp end do nowait

  !$omp do schedule(static) collapse(2)
  do ix2 = start2, end2
    do iz=1,s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E_F+E+iy)
      call greenlinearsoc(x2(ix2)+e,eta,kp,gf,gvg)
      gfuu(:,:,:,:,1) = gf(:,:,     1:  nOrb,     1:  nOrb)
      gfud(:,:,:,:,1) = gf(:,:,     1:  nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,     1:  nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,1) = gvg(:,:,     1:  nOrb,     1:  nOrb)
      gvgud(:,:,:,:,1) = gvg(:,:,     1:  nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb,     1:  nOrb)
      gvgdd(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(x2(ix2),eta,kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(:,:,     1:  nOrb,     1: nOrb)
      gfud(:,:,:,:,2) = gf(:,:,     1:  nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,     1: nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,2) = gvg(:,:,     1:  nOrb,     1:  nOrb)
      gvgud(:,:,:,:,2) = gvg(:,:,     1:  nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb,     1:  nOrb)
      gvgdd(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))
      end do ; end do ; end do ; end do ; end do ; end do

      Fint = Fint + df1 * s%wkbz(iz) * p2(ix2)
      Fintlsoc = Fintlsoc + df1lsoc * s%wkbz(iz) * p2(ix2)
    end do
  end do
  !$omp end do nowait

  Fint = Fint / tpi
  Fintlsoc = Fintlsoc / tpi

  !$omp critical
    chiorb_hf = chiorb_hf + Fint
    chiorb_hflsoc = chiorb_hflsoc + Fintlsoc
  !$omp end critical

  deallocate(df1, df1lsoc)
  deallocate(Fint, Fintlsoc)
  deallocate(gvg,gvguu,gvgud,gvgdu,gvgdd)
  deallocate(gf,gfuu,gfud,gfdu,gfdd)
  !$omp end parallel
  if(myrank_row == 0) then
    call MPI_Reduce(MPI_IN_PLACE, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
    call MPI_Reduce(MPI_IN_PLACE, chiorb_hflsoc, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
  else
    call MPI_Reduce(chiorb_hf, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
    call MPI_Reduce(chiorb_hflsoc, chiorb_hflsoc, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
  end if

  return
end subroutine eintshechilinearsoc
