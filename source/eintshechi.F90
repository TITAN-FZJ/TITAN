! ---------- Spin disturbance: Energy integration ---------
subroutine eintshechi(q,e)
  use mod_f90_kind,         only: double
  use mod_constants,        only: cZero, cI, tpi
  use mod_parameters,       only: nOrb, nOrb2, eta, etap, dim, sigmaimunu2i
  use EnergyIntegration,    only: generate_real_epoints, y, wght, x2, p2, pn2
  use mod_susceptibilities, only: chiorb_hf
  use mod_system,           only: s => sys
  use mod_BrillouinZone,    only: realBZ
  use mod_SOC,              only: llineargfsoc
  use mod_mpi_pars
  use adaptiveMesh
  implicit none
  real(double), intent(in)    :: e, q(3)

  integer*4 :: AllocateStatus
  complex(double), dimension(:,:),allocatable :: Fint

  integer*4       :: i,j,mu,nu,gamma,xi
  real(double)    :: kp(3),kpq(3)
  complex(double),dimension(:,:,:,:),allocatable    :: gf
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd
  real(double) :: weight, ep
  integer*4, dimension(4) :: index1, index2

!--------------------- begin MPI vars --------------------
  integer*8 :: ix
  integer*8 :: ix2, nep,nkp
  integer*8 :: real_points
  integer*4 :: ncount
  ncount=dim*dim
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  real_points = 0
  if(abs(e) >= 1.d-15) real_points = int(pn2*realBZ%workload,8)

  ! Starting to calculate energy integral

  allocate(Fint(dim,dim), STAT = AllocateStatus  )
  if (AllocateStatus/=0) &
  call abortProgram("[eintshechi] Not enough memory for: Fint")
  Fint = cZero
  chiorb_hf = cZero

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,ix2,i,j,mu,nu,gamma,xi,nep,nkp,ep,kp,kpq,weight,gf,gfuu,gfud,gfdu,gfdd,index1, index2) &
  !$omp& shared(llineargfsoc,bzs,s,nOrb,nOrb2,realBZ,local_points,q,e,y,wght,x2,p2,pn2,real_points,E_k_imag_mesh,eta,etap,dim,sigmaimunu2i,Fint,chiorb_hf)
  allocate(gf  (nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
           gfuu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfud(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdd(nOrb,nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) &
  call abortProgram("[eintshechi] Not enough memory for: gf,gfuu,gfud,gfdu,gfdd")

  !$omp do schedule(static) reduction(+:chiorb_hf)
  do ix = 1, local_points
      ep = y( E_k_imag_mesh(1,ix) )
      kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
      weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix))
      ! Green function at (k+q,E_F+E+iy)
      kpq = kp+q
      if(llineargfsoc) then
        call greenlineargfsoc(s%Ef+e,ep+eta,s,kpq,gf)
      else
        call green(s%Ef+e,ep+eta,s,kpq,gf)
      end if
      gfuu(:,:,:,:,1) = gf(     1:nOrb ,     1:nOrb ,:,:)
      gfud(:,:,:,:,1) = gf(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu(:,:,:,:,1) = gf(nOrb+1:nOrb2,     1:nOrb ,:,:)
      gfdd(:,:,:,:,1) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

      ! Green function at (k,E_F+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(s%Ef,ep+etap,s,kp,gf)
      else
        call green(s%Ef,ep+etap,s,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(     1:nOrb ,     1:nOrb ,:,:)
      gfud(:,:,:,:,2) = gf(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu(:,:,:,:,2) = gf(nOrb+1:nOrb2,     1:nOrb ,:,:)
      gfdd(:,:,:,:,2) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

      !dir$ ivdep:loop
      do j=1,s%nAtoms
        !dir$ ivdep:loop
        do i=1,s%nAtoms
          !dir$ ivdep:loop
          do gamma=1,nOrb
            !dir$ ivdep:loop
            do mu=1,nOrb
              !dir$ ivdep:loop
              do xi=1,nOrb
                !dir$ ivdep:loop
                do nu=1,nOrb
                  index1(1:4) = sigmaimunu2i(1:4,i,mu,nu)
                  index2(1:4) = sigmaimunu2i(1:4,j,gamma,xi)
                  chiorb_hf(index1(1),index2(1)) = chiorb_hf(index1(1),index2(1)) + weight*(gfdd(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  chiorb_hf(index1(1),index2(2)) = chiorb_hf(index1(1),index2(2)) + weight*(gfdu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))
                  chiorb_hf(index1(1),index2(3)) = chiorb_hf(index1(1),index2(3)) + weight*(gfdd(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  chiorb_hf(index1(1),index2(4)) = chiorb_hf(index1(1),index2(4)) + weight*(gfdu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))

                  chiorb_hf(index1(2),index2(1)) = chiorb_hf(index1(2),index2(1)) + weight*(gfud(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  chiorb_hf(index1(2),index2(2)) = chiorb_hf(index1(2),index2(2)) + weight*(gfuu(nu,gamma,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
                  chiorb_hf(index1(2),index2(3)) = chiorb_hf(index1(2),index2(3)) + weight*(gfud(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  chiorb_hf(index1(2),index2(4)) = chiorb_hf(index1(2),index2(4)) + weight*(gfuu(nu,gamma,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))

                  chiorb_hf(index1(3),index2(1)) = chiorb_hf(index1(3),index2(1)) + weight*(gfdd(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  chiorb_hf(index1(3),index2(2)) = chiorb_hf(index1(3),index2(2)) + weight*(gfdu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))
                  chiorb_hf(index1(3),index2(3)) = chiorb_hf(index1(3),index2(3)) + weight*(gfdd(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdd(gamma,nu,j,i,1)))
                  chiorb_hf(index1(3),index2(4)) = chiorb_hf(index1(3),index2(4)) + weight*(gfdu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfud(gamma,nu,j,i,1)))

                  chiorb_hf(index1(4),index2(1)) = chiorb_hf(index1(4),index2(1)) + weight*(gfud(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  chiorb_hf(index1(4),index2(2)) = chiorb_hf(index1(4),index2(2)) + weight*(gfuu(nu,gamma,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
                  chiorb_hf(index1(4),index2(3)) = chiorb_hf(index1(4),index2(3)) + weight*(gfud(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdu(gamma,nu,j,i,1)))
                  chiorb_hf(index1(4),index2(4)) = chiorb_hf(index1(4),index2(4)) + weight*(gfuu(nu,gamma,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfuu(gamma,nu,j,i,1)))
                end do
              end do
            end do
          end do
        end do
      end do
      !call zaxpy(dim*dim,weight,df1,1,Fint,1)
      !Fint = Fint + df1*weight
  end do !end ix <= pn1 loop
  !$omp end do nowait

  !$omp do schedule(static) reduction(+:Fint)
  do ix2 = 1, real_points ! Third integration (on the real axis)
      nep = (ix2-1) / int(realBZ % workload,8) + 1
      nkp = mod(ix2-1, int(realBZ % workload,8)) + 1
      ep = x2(nep)
      kp = realBZ % kp(:,nkp)
      weight = p2(nep) * realBZ % w(nkp)
      ! Green function at (k+q,E'+E+i.eta)
      kpq = kp+q
      if(llineargfsoc) then
        call greenlineargfsoc(ep+e,eta,s,kpq,gf)
      else
        call green(ep+e,eta,s,kpq,gf)
      end if
      gfuu(:,:,:,:,1) = gf(     1:nOrb ,     1:nOrb ,:,:)
      gfud(:,:,:,:,1) = gf(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu(:,:,:,:,1) = gf(nOrb+1:nOrb2,     1:nOrb ,:,:)
      gfdd(:,:,:,:,1) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

      ! Green function at (k,E'+i.eta)
      if(llineargfsoc) then
        call greenlineargfsoc(ep,etap,s,kp,gf)
      else
        call green(ep,etap,s,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(     1:nOrb ,     1:nOrb ,:,:)
      gfud(:,:,:,:,2) = gf(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu(:,:,:,:,2) = gf(nOrb+1:nOrb2,     1:nOrb ,:,:)
      gfdd(:,:,:,:,2) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

      !dir$ ivdep:loop
      do j=1,s%nAtoms
        !dir$ ivdep:loop
        do i=1,s%nAtoms
          !dir$ ivdep:loop
          do gamma=1,nOrb
            !dir$ ivdep:loop
            do mu=1,nOrb
              !dir$ ivdep:loop
              do xi=1,nOrb
                !dir$ ivdep:loop
                do nu=1,nOrb
                  index1(1:4) = sigmaimunu2i(1:4,i,mu,nu)
                  index2(1:4) = sigmaimunu2i(1:4,j,gamma,xi)
                  Fint(index1(1),index2(1)) = Fint(index1(1),index2(1)) - weight * cI*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  Fint(index1(1),index2(2)) = Fint(index1(1),index2(2)) - weight * cI*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  Fint(index1(1),index2(3)) = Fint(index1(1),index2(3)) - weight * cI*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
                  Fint(index1(1),index2(4)) = Fint(index1(1),index2(4)) - weight * cI*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

                  Fint(index1(2),index2(1)) = Fint(index1(2),index2(1)) - weight * cI*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  Fint(index1(2),index2(2)) = Fint(index1(2),index2(2)) - weight * cI*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
                  Fint(index1(2),index2(3)) = Fint(index1(2),index2(3)) - weight * cI*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
                  Fint(index1(2),index2(4)) = Fint(index1(2),index2(4)) - weight * cI*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

                  Fint(index1(3),index2(1)) = Fint(index1(3),index2(1)) - weight * cI*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  Fint(index1(3),index2(2)) = Fint(index1(3),index2(2)) - weight * cI*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  Fint(index1(3),index2(3)) = Fint(index1(3),index2(3)) - weight * cI*(gfdd(nu,gamma,i,j,1) - conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
                  Fint(index1(3),index2(4)) = Fint(index1(3),index2(4)) - weight * cI*(gfdu(nu,gamma,i,j,1) - conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))

                  Fint(index1(4),index2(1)) = Fint(index1(4),index2(1)) - weight * cI*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  Fint(index1(4),index2(2)) = Fint(index1(4),index2(2)) - weight * cI*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
                  Fint(index1(4),index2(3)) = Fint(index1(4),index2(3)) - weight * cI*(gfud(nu,gamma,i,j,1) - conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
                  Fint(index1(4),index2(4)) = Fint(index1(4),index2(4)) - weight * cI*(gfuu(nu,gamma,i,j,1) - conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
                end do
              end do
            end do
          end do
        end do
      end do

      ! Locally add up df1
      ! call zaxpy(dim*dim,(p2(ix2)*s%wkbz(iz)),df1,1,Fint,1)
      !Fint = Fint + df1*(p2(ix2)*s%wkbz(iz))
  end do !end pn1+1, nepoints loop
  !$omp end do nowait

  deallocate(gf,gfuu,gfud,gfdu,gfdd)
  !$omp end parallel

  chiorb_hf = (chiorb_hf + Fint) / tpi

  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1), ierr)

  deallocate(Fint)
end subroutine eintshechi

! -------------------- Spin disturbance: Energy integration --------------------
! -------------- to be used in the calculation of linear SOC chi ---------------
subroutine eintshechilinearsoc(q,e)
  use mod_f90_kind,         only: double
  use mod_constants,        only: cZero, cI, tpi
  use mod_parameters,       only: nOrb, nOrb2, eta, etap, dim, sigmaimunu2i
  use EnergyIntegration,    only: generate_real_epoints,y, wght, x2, p2, pn2
  use mod_susceptibilities, only: chiorb_hf,chiorb_hflsoc
  use mod_system,           only: s => sys
  use mod_BrillouinZone,    only: realBZ
  use mod_mpi_pars
  use adaptiveMesh
  implicit none
  real(double), intent(in) :: e,q(3)

  integer*4 :: AllocateStatus
  integer*4 :: i,j,mu,nu,gamma,xi
  integer*8 :: ix,ix2, nep, nkp
  real(double) :: kp(3), kpq(3), ep
  real(double) :: weight
  complex(double), dimension(:,:,:,:), allocatable :: gf,gvg
  complex(double), dimension(:,:,:,:,:), allocatable :: gfuu,gfud,gfdu,gfdd
  complex(double), dimension(:,:,:,:,:), allocatable :: gvguu,gvgud,gvgdu,gvgdd
  complex(double), dimension(:,:), allocatable :: Fint,Fintlsoc
  complex(double), dimension(:,:), allocatable :: df1,df1lsoc

  !--------------------- begin MPI vars --------------------
  integer*4 :: ncount
  integer*8 :: real_points
  ncount=dim*dim
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  real_points = 0
  if(abs(e) >= 1.d-15) real_points = int(pn2*realBZ%workload,8)

  !------------------------------------------------------

  ! Starting to calculate energy integral
  chiorb_hf     = cZero
  chiorb_hflsoc = cZero

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,ix2,i,j,mu,nu,nep,nkp,gamma,xi,kp,kpq,ep,weight,Fint,Fintlsoc,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,df1,df1lsoc) &
  !$omp& shared(local_points,s,nOrb,nOrb2,bzs,realBZ,real_points,E_k_imag_mesh,q,e,y,wght,x2,p2,eta,etap,dim,sigmaimunu2i,chiorb_hf,chiorb_hflsoc)
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

  Fint      = cZero
  Fintlsoc  = cZero

  ! Starting to calculate energy integral
  !$omp do schedule(static)
  do ix = 1, local_points
      ep = y( E_k_imag_mesh(1,ix) )
      kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
      weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix))
      ! Green function at (k+q,E_F+E+iy)
      kpq = kp+q
      call greenlinearsoc(s%Ef+e,ep+eta,s,kpq,gf,gvg)
      gfuu (:,:,:,:,1) = gf (     1:nOrb ,     1:nOrb ,:,:)
      gfud (:,:,:,:,1) = gf (     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu (:,:,:,:,1) = gf (nOrb+1:nOrb2,     1:nOrb ,:,:)
      gfdd (:,:,:,:,1) = gf (nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)
      gvguu(:,:,:,:,1) = gvg(     1:nOrb ,     1:nOrb ,:,:)
      gvgud(:,:,:,:,1) = gvg(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gvgdu(:,:,:,:,1) = gvg(nOrb+1:nOrb2,     1:nOrb ,:,:)
      gvgdd(:,:,:,:,1) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(s%Ef,ep+etap,s,kp,gf,gvg)
      gfuu (:,:,:,:,2) = gf (     1:nOrb ,     1: nOrb,:,:)
      gfud (:,:,:,:,2) = gf (     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu (:,:,:,:,2) = gf (nOrb+1:nOrb2,     1: nOrb,:,:)
      gfdd (:,:,:,:,2) = gf (nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)
      gvguu(:,:,:,:,2) = gvg(     1:nOrb ,     1: nOrb,:,:)
      gvgud(:,:,:,:,2) = gvg(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gvgdu(:,:,:,:,2) = gvg(nOrb+1:nOrb2,     1: nOrb,:,:)
      gvgdd(:,:,:,:,2) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

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

      Fint = Fint + df1 * weight
      Fintlsoc = Fintlsoc + df1lsoc * weight

  end do
  !$omp end do nowait

  !$omp do schedule(static)
  do ix2 = 1, realBZ%workload ! Third integration (on the real axis)
      nep = (ix2-1) / int(realBZ % workload,8) + 1
      nkp = mod(ix2-1, int(realBZ % workload,8))+1
      ep = x2(nep)
      kp = realBZ % kp(:,nkp)
      weight = p2(nep) * realBZ % w(nkp)

      ! Green function at (k+q,E_F+E+iy)
      kpq = kp+q
      call greenlinearsoc(ep+e,eta,s,kpq,gf,gvg)
      gfuu (:,:,:,:,1) = gf (     1:nOrb ,     1:nOrb ,:,:)
      gfud (:,:,:,:,1) = gf (     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu (:,:,:,:,1) = gf (nOrb+1:nOrb2,     1:nOrb ,:,:)
      gfdd (:,:,:,:,1) = gf (nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)
      gvguu(:,:,:,:,1) = gvg(     1:nOrb ,     1:nOrb ,:,:)
      gvgud(:,:,:,:,1) = gvg(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gvgdu(:,:,:,:,1) = gvg(nOrb+1:nOrb2,     1:nOrb ,:,:)
      gvgdd(:,:,:,:,1) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(ep,etap,s,kp,gf,gvg)
      gfuu (:,:,:,:,2) = gf (     1:nOrb ,     1: nOrb,:,:)
      gfud (:,:,:,:,2) = gf (     1:nOrb ,nOrb+1:nOrb2,:,:)
      gfdu (:,:,:,:,2) = gf (nOrb+1:nOrb2,     1: nOrb,:,:)
      gfdd (:,:,:,:,2) = gf (nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)
      gvguu(:,:,:,:,2) = gvg(     1:nOrb ,     1: nOrb,:,:)
      gvgud(:,:,:,:,2) = gvg(     1:nOrb ,nOrb+1:nOrb2,:,:)
      gvgdu(:,:,:,:,2) = gvg(nOrb+1:nOrb2,     1: nOrb,:,:)
      gvgdd(:,:,:,:,2) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2,:,:)

      !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
      end do ; end do ; end do ; end do ; end do ; end do

      Fint = Fint + df1 * weight
      Fintlsoc = Fintlsoc + df1lsoc * weight
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

  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hf    , ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1), ierr)
  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hflsoc, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1), ierr)

end subroutine eintshechilinearsoc
