! ---------- Spin disturbance: Energy integration ---------
subroutine eintshechi(q,e)
  use mod_kind,             only: dp,int32,int64
  use mod_constants,        only: cZero,cI,tpi
  use mod_parameters,       only: eta,etap,dimens,sigmaimunu2i
  use EnergyIntegration,    only: generate_real_epoints,y,wght,x2,p2,pn2
  use mod_susceptibilities, only: chiorb_hf
  use mod_system,           only: s => sys
  use mod_BrillouinZone,    only: realBZ
  use mod_SOC,              only: llineargfsoc
  use mod_hamiltonian,      only: hamilt_local
  use mod_greenfunction,    only: calc_green
  use adaptiveMesh,         only: bzs,E_k_imag_mesh,local_points
  use mod_mpi_pars,         only: MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr,abortProgram
  implicit none
  real(dp), intent(in)    :: e, q(3)

  integer(int32) :: AllocateStatus
  complex(dp), dimension(:,:),allocatable :: Fint

  integer(int32), dimension(4) :: index1, index2
  integer(int32) :: i,j,mu,nu,gamma,xi
  real(dp)       :: kp(3),kpq(3)
  real(dp)       :: weight, ep
  integer(int64) :: ix
  integer(int64) :: ix2, nep,nkp
  integer(int64) :: real_points
  integer(int32) :: ncount
  complex(dp), dimension(:,:,:,:),   allocatable  :: gf
  complex(dp), dimension(:,:,:,:,:), allocatable  :: gfuu,gfud,gfdu,gfdd

  external :: MPI_Allreduce

  ncount=dimens*dimens

  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  real_points = 0
  if(abs(e) >= 1.e-15_dp) real_points = int(pn2*realBZ%workload,8)

  ! Build local hamiltonian
  if(.not.llineargfsoc) call hamilt_local(s)

  ! Starting to calculate energy integral
  allocate(Fint(dimens,dimens), STAT = AllocateStatus  )
  if (AllocateStatus/=0) &
    call abortProgram("[eintshechi] Not enough memory for: Fint")
  Fint = cZero
  chiorb_hf = cZero

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,ix2,i,j,mu,nu,gamma,xi,nep,nkp,ep,kp,kpq,weight,gf,gfuu,gfud,gfdu,gfdd,index1, index2) &
  !$omp& shared(bzs,s,calc_green,realBZ,local_points,q,e,y,wght,x2,p2,pn2,real_points,E_k_imag_mesh,eta,etap,dimens,sigmaimunu2i,Fint,chiorb_hf)
  allocate(gf  (s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms), &
           gfuu(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
           gfud(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
           gfdu(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
           gfdd(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) &
    call abortProgram("[eintshechi] Not enough memory for: gf,gfuu,gfud,gfdu,gfdd")

  !$omp do schedule(static) reduction(+:chiorb_hf)
  do ix = 1, local_points
    ep = y( E_k_imag_mesh(1,ix) )
    kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
    weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix))
    ! Green function at (k+q,E_F+E+iy)
    kpq = kp+q
    call calc_green(s%Ef+e,ep+eta,s,kpq,gf)
    gfuu(:,:,:,:,1) = gf(       1:s%nOrb ,       1:s%nOrb ,:,:)
    gfud(:,:,:,:,1) = gf(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu(:,:,:,:,1) = gf(s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gfdd(:,:,:,:,1) = gf(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    ! Green function at (k,E_F+iy)
    call calc_green(s%Ef,ep+etap,s,kp,gf)
    gfuu(:,:,:,:,2) = gf(       1:s%nOrb ,       1:s%nOrb ,:,:)
    gfud(:,:,:,:,2) = gf(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu(:,:,:,:,2) = gf(s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gfdd(:,:,:,:,2) = gf(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    do j=1,s%nAtoms; do i=1,s%nAtoms; do gamma=1,s%nOrb; do mu=1,s%nOrb; do xi=1,s%nOrb; do nu=1,s%nOrb
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
    end do; end do; end do; end do; end do; end do
    !call zaxpy(dimens*dimens,weight,df1,1,Fint,1)
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
    call calc_green(ep+e,eta,s,kpq,gf)
    gfuu(:,:,:,:,1) = gf(       1:s%nOrb ,       1:s%nOrb ,:,:)
    gfud(:,:,:,:,1) = gf(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu(:,:,:,:,1) = gf(s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gfdd(:,:,:,:,1) = gf(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    ! Green function at (k,E'+i.eta)
    call calc_green(ep,etap,s,kp,gf)
    gfuu(:,:,:,:,2) = gf(       1:s%nOrb ,       1:s%nOrb ,:,:)
    gfud(:,:,:,:,2) = gf(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu(:,:,:,:,2) = gf(s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gfdd(:,:,:,:,2) = gf(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    do j=1,s%nAtoms; do i=1,s%nAtoms; do gamma=1,s%nOrb; do mu=1,s%nOrb; do xi=1,s%nOrb; do nu=1,s%nOrb
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
    end do; end do; end do; end do; end do; end do

    ! Locally add up df1
    ! call zaxpy(dimens*dimens,(p2(ix2)*s%wkbz(iz)),df1,1,Fint,1)
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
  use mod_kind,             only: dp,int32,int64
  use mod_constants,        only: cZero, cI, tpi
  use mod_parameters,       only: eta,etap,dimens,sigmaimunu2i
  use EnergyIntegration,    only: generate_real_epoints,y,wght,x2,p2,pn2
  use mod_susceptibilities, only: chiorb_hf,chiorb_hflsoc
  use mod_system,           only: s => sys
  use mod_BrillouinZone,    only: realBZ
  use mod_greenfunction,    only: greenlinearsoc
  use adaptiveMesh,         only: bzs,E_k_imag_mesh,local_points
  use mod_mpi_pars,         only: MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr,abortProgram
  implicit none
  real(dp), intent(in) :: e,q(3)

  integer(int32) :: AllocateStatus
  integer(int32) :: i,j,mu,nu,gamma,xi
  integer(int64) :: ix,ix2, nep, nkp
  real(dp) :: kp(3), kpq(3), ep
  real(dp) :: weight
  complex(dp), dimension(:,:,:,:), allocatable :: gf,gvg
  complex(dp), dimension(:,:,:,:,:), allocatable :: gfuu,gfud,gfdu,gfdd
  complex(dp), dimension(:,:,:,:,:), allocatable :: gvguu,gvgud,gvgdu,gvgdd
  complex(dp), dimension(:,:), allocatable :: Fint,Fintlsoc
  complex(dp), dimension(:,:), allocatable :: df1,df1lsoc
  integer(int32) :: ncount
  integer(int64) :: real_points

  external :: MPI_Allreduce

  ncount=dimens*dimens


  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  real_points = 0
  if(abs(e) >= 1.e-15_dp) real_points = int(pn2*realBZ%workload,8)

  !------------------------------------------------------

  ! Starting to calculate energy integral
  chiorb_hf     = cZero
  chiorb_hflsoc = cZero

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,ix2,i,j,mu,nu,nep,nkp,gamma,xi,kp,kpq,ep,weight,Fint,Fintlsoc,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,df1,df1lsoc) &
  !$omp& shared(local_points,s,bzs,realBZ,real_points,E_k_imag_mesh,q,e,y,wght,x2,p2,eta,etap,dimens,sigmaimunu2i,chiorb_hf,chiorb_hflsoc)
  allocate(df1(dimens,dimens), Fint(dimens,dimens), &
           gf(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms), &
           gfuu(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
           gfud(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
           gfdu(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
           gfdd(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[eintshechilinearsoc] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")

  allocate( df1lsoc(dimens,dimens), Fintlsoc(dimens,dimens), &
            gvg(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms), &
            gvguu(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
            gvgud(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
            gvgdu(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), &
            gvgdd(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
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
    gfuu (:,:,:,:,1) = gf (       1:s%nOrb ,       1:s%nOrb ,:,:)
    gfud (:,:,:,:,1) = gf (       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu (:,:,:,:,1) = gf (s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gfdd (:,:,:,:,1) = gf (s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)
    gvguu(:,:,:,:,1) = gvg(       1:s%nOrb ,       1:s%nOrb ,:,:)
    gvgud(:,:,:,:,1) = gvg(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gvgdu(:,:,:,:,1) = gvg(s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gvgdd(:,:,:,:,1) = gvg(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    ! Green function at (k,E_F+iy)
    call greenlinearsoc(s%Ef,ep+etap,s,kp,gf,gvg)
    gfuu (:,:,:,:,2) = gf (       1:s%nOrb ,       1: s%nOrb,:,:)
    gfud (:,:,:,:,2) = gf (       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu (:,:,:,:,2) = gf (s%nOrb+1:s%nOrb2,       1: s%nOrb,:,:)
    gfdd (:,:,:,:,2) = gf (s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)
    gvguu(:,:,:,:,2) = gvg(       1:s%nOrb ,       1: s%nOrb,:,:)
    gvgud(:,:,:,:,2) = gvg(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gvgdu(:,:,:,:,2) = gvg(s%nOrb+1:s%nOrb2,       1: s%nOrb,:,:)
    gvgdd(:,:,:,:,2) = gvg(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    do j=1,s%nAtoms; do i=1,s%nAtoms; do gamma=1,s%nOrb; do mu=1,s%nOrb; do xi=1,s%nOrb; do nu=1,s%nOrb
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
    end do; end do; end do; end do; end do; end do

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
    gfuu (:,:,:,:,1) = gf (       1:s%nOrb ,       1:s%nOrb ,:,:)
    gfud (:,:,:,:,1) = gf (       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu (:,:,:,:,1) = gf (s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gfdd (:,:,:,:,1) = gf (s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)
    gvguu(:,:,:,:,1) = gvg(       1:s%nOrb ,       1:s%nOrb ,:,:)
    gvgud(:,:,:,:,1) = gvg(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gvgdu(:,:,:,:,1) = gvg(s%nOrb+1:s%nOrb2,       1:s%nOrb ,:,:)
    gvgdd(:,:,:,:,1) = gvg(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    ! Green function at (k,E_F+iy)
    call greenlinearsoc(ep,etap,s,kp,gf,gvg)
    gfuu (:,:,:,:,2) = gf (       1:s%nOrb ,       1: s%nOrb,:,:)
    gfud (:,:,:,:,2) = gf (       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gfdu (:,:,:,:,2) = gf (s%nOrb+1:s%nOrb2,       1: s%nOrb,:,:)
    gfdd (:,:,:,:,2) = gf (s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)
    gvguu(:,:,:,:,2) = gvg(       1:s%nOrb ,       1: s%nOrb,:,:)
    gvgud(:,:,:,:,2) = gvg(       1:s%nOrb ,s%nOrb+1:s%nOrb2,:,:)
    gvgdu(:,:,:,:,2) = gvg(s%nOrb+1:s%nOrb2,       1: s%nOrb,:,:)
    gvgdd(:,:,:,:,2) = gvg(s%nOrb+1:s%nOrb2,s%nOrb+1:s%nOrb2,:,:)

    do j=1,s%nAtoms; do i=1,s%nAtoms; do gamma=1,s%nOrb; do mu=1,s%nOrb; do xi=1,s%nOrb; do nu=1,s%nOrb
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
    end do; end do; end do; end do; end do; end do

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
