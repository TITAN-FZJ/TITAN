! ---------- Spin disturbance: Energy integration ---------
subroutine eintshechi(e)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, zum, zi, tpi
  use mod_parameters, only: eta, ef, dim, sigmaijmunu2i, sigmaimunu2i
  use EnergyIntegration, only: generate_real_epoints, y, wght, x2, p2, nepoints, pn1
  use mod_susceptibilities, only: chiorb_hf
  use mod_system, only: s => sys
  use TightBinding, only: nOrb,nOrb2
  use mod_SOC, only: llineargfsoc
  use mod_mpi_pars
  !$  use omp_lib
  implicit none
  real(double), intent(in)    :: e

  integer :: AllocateStatus
  complex(double), dimension(:,:),allocatable :: Fint



  integer         :: i,j,mu,nu,gamma,xi,iz
  real(double)                :: kp(3)
  complex(double),dimension(:,:,:,:),allocatable    :: gf
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd
  complex(double),dimension(:,:),allocatable        :: df1
  real(double) :: weight


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
  !$omp& private(AllocateStatus,iz,ix,ix2,i,j,mu,nu,gamma,xi,kp,weight,gf,gfuu,gfud,gfdu,gfdd,df1,Fint) &
  !$omp& shared(llineargfsoc,start1,start2,end1,end2,s,e,y,x2,wght,p2,Ef,eta,dim,sigmaimunu2i,sigmaijmunu2i,chiorb_hf)
  allocate(df1(dim,dim), Fint(dim,dim), &
           gf  (nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
           gfuu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfud(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdd(nOrb,nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[eintshechi] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")
  Fint = zero

  !$omp do schedule(static) collapse(2)
  do ix = start1, end1
    do iz = 1, s%nkpt
      kp = s%kbz(:,iz)
      weight = wght(ix)*s%wkbz(iz)
      ! Green function at (k+q,E_F+E+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef+e,y(ix),kp,gf)
      else
        call green(Ef+e,y(ix),kp,gf)
      end if
      gfuu(:,:,:,:,1) = gf(     1:  nOrb,     1:  nOrb, :,:)
      gfud(:,:,:,:,1) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gfdu(:,:,:,:,1) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gfdd(:,:,:,:,1) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

      ! Green function at (k,E_F+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef,y(ix),kp,gf)
      else
        call green(Ef,y(ix),kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(     1:  nOrb,     1:  nOrb, :,:)
      gfud(:,:,:,:,2) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gfdu(:,:,:,:,2) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gfdd(:,:,:,:,2) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

      !dir$ simd
      do xi=1,nOrb
        do nu=1,nOrb
          do gamma=1,nOrb
            do mu=1,nOrb
              do j=1,s%nAtoms
                do i=1,s%nAtoms
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))

                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))

                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))

                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
                end do
              end do
            end do
          end do
        end do
      end do

      Fint = Fint + df1*weight
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
      gfuu(:,:,:,:,1) = gf(     1:  nOrb,     1:  nOrb, :,:)
      gfud(:,:,:,:,1) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gfdu(:,:,:,:,1) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gfdd(:,:,:,:,1) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

      ! Green function at (k,E'+i.eta)
      if(llineargfsoc) then
        call greenlineargfsoc(x2(ix2),eta,kp,gf)
      else
        call green(x2(ix2),eta,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(     1:  nOrb,     1:  nOrb, :,:)
      gfud(:,:,:,:,2) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gfdu(:,:,:,:,2) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gfdd(:,:,:,:,2) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

      !dir$ simd
      do xi=1,nOrb
        do nu=1,nOrb
          do gamma=1,nOrb
            do mu=1,nOrb
              do j=1,s%nAtoms
                do i=1,s%nAtoms
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))

                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
                end do
              end do
            end do
          end do
        end do
      end do

      ! Locally add up df1
      Fint = Fint + df1*(p2(ix2)*s%wkbz(iz))
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
subroutine eintshechilinearsoc(e)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, zum, zi, tpi
  use mod_parameters, only: eta, ef, dim, sigmaijmunu2i, sigmaimunu2i
  use EnergyIntegration, only: generate_real_epoints,y, wght, x2, p2, nepoints, pn1
  use mod_susceptibilities, only: chiorb_hf,chiorb_hflsoc
  use mod_mpi_pars
  use mod_system, only: s => sys
  use TightBinding, only: nOrb,nOrb2
  use mod_SOC, only: llineargfsoc
  !$  use omp_lib

  implicit none
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
           gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
           gfuu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfud(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdd(nOrb,nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[eintshechilinearsoc] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")

  allocate( df1lsoc(dim,dim), Fintlsoc(dim,dim), &
            gvg(nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
            gvguu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gvgud(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gvgdu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gvgdd(nOrb,nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
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
      gfuu(:,:,:,:,1) = gf(     1:  nOrb,      1: nOrb, :,:)
      gfud(:,:,:,:,1) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gfdu(:,:,:,:,1) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gfdd(:,:,:,:,1) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)
      gvguu(:,:,:,:,1) = gvg(     1:  nOrb,     1:  nOrb, :,:)
      gvgud(:,:,:,:,1) = gvg(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gvgdu(:,:,:,:,1) = gvg(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gvgdd(:,:,:,:,1) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(Ef,y(ix),kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(:,:,     1:  nOrb,     1:  nOrb)
      gfud(:,:,:,:,2) = gf(:,:,     1:  nOrb,nOrb+1:nOrb2)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:nOrb2,     1:  nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:nOrb2,nOrb+1:nOrb2)
      gvguu(:,:,:,:,2) = gvg(:,:,     1:  nOrb,     1:  nOrb)
      gvgud(:,:,:,:,2) = gvg(:,:,     1:  nOrb,nOrb+1:nOrb2)
      gvgdu(:,:,:,:,2) = gvg(:,:,nOrb+1:nOrb2,     1:  nOrb)
      gvgdd(:,:,:,:,2) = gvg(:,:,nOrb+1:nOrb2,nOrb+1:nOrb2)

      !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvgdd(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvgud(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvgdd(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvgud(gamma,nu,j,i,1)))

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvgdu(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvguu(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvgdu(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvguu(gamma,nu,j,i,1)))

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvgdd(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvgud(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvgdd(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvgud(gamma,nu,j,i,1)))

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvgdu(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvguu(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvgdu(gamma,nu,j,i,1)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvguu(gamma,nu,j,i,1)))
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
      gfuu(:,:,:,:,1) = gf(     1:  nOrb,     1:  nOrb, :,:)
      gfud(:,:,:,:,1) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gfdu(:,:,:,:,1) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gfdd(:,:,:,:,1) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)
      gvguu(:,:,:,:,1) = gvg(     1:  nOrb,     1:  nOrb, :,:)
      gvgud(:,:,:,:,1) = gvg(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gvgdu(:,:,:,:,1) = gvg(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gvgdd(:,:,:,:,1) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(x2(ix2),eta,kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(     1:  nOrb,     1:  nOrb, :,:)
      gfud(:,:,:,:,2) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gfdu(:,:,:,:,2) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gfdd(:,:,:,:,2) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)
      gvguu(:,:,:,:,2) = gvg(     1:  nOrb,     1:  nOrb, :,:)
      gvgud(:,:,:,:,2) = gvg(     1:  nOrb,nOrb+1:nOrb2, :,:)
      gvgdu(:,:,:,:,2) = gvg(nOrb+1:nOrb2,     1:  nOrb, :,:)
      gvgdd(:,:,:,:,2) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

      !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
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
