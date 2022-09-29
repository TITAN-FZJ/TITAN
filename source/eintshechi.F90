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

  complex(dp), dimension(dimens,dimens) :: Fint

  integer(int32), dimension(4) :: index1, index2
  integer(int32) :: i,j,mu,nu,gama,xi,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j
  real(dp)       :: kp(3),kpq(3)
  real(dp)       :: weight, ep
  integer(int64) :: ix
  integer(int64) :: ix2, nep,nkp
  integer(int64) :: real_points
  integer(int32) :: ncount
  complex(dp), dimension(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms) :: gf,gfq
  complex(dp), dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2) :: gfuu,gfud,gfdu,gfdd

  external :: MPI_Allreduce

  ncount=dimens*dimens

  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  real_points = 0
  if(abs(e) >= 1.e-15_dp) real_points = int(pn2*realBZ%workload,8)

  ! Build local hamiltonian
  if(.not.llineargfsoc) call hamilt_local(s)

  ! Starting to calculate energy integral
  Fint = cZero
  chiorb_hf = cZero

  !$omp parallel default(none) &
  !$omp& private(ix,ix2,i,j,mu,nu,gama,xi,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j,nep,nkp,ep,kp,kpq,weight,gf,gfq,gfuu,gfud,gfdu,gfdd,index1,index2) &
  !$omp& shared(bzs,s,calc_green,realBZ,local_points,q,e,y,wght,x2,p2,pn2,real_points,E_k_imag_mesh,eta,etap,dimens,sigmaimunu2i,Fint,chiorb_hf)

  !$omp do schedule(dynamic) reduction(+:chiorb_hf)
  do ix = 1, local_points
    ep = y( E_k_imag_mesh(1,ix) )
    kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
    weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix))

    ! Green function at (k+q,E_F+E+iy)
    kpq = kp+q
    call calc_green(s%Ef+e,ep+eta,s,kpq,gfq)

    ! Green function at (k,E_F+iy)
    call calc_green(s%Ef,ep+etap,s,kp,gf)

    do j=1,s%nAtoms; do i=1,s%nAtoms
      nOrb_i = s%Types(s%Basis(i)%Material)%nOrb
      nOrb_j = s%Types(s%Basis(j)%Material)%nOrb
      nOrb2_i = 2*nOrb_i
      nOrb2_j = 2*nOrb_j

      gfuu(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)

      gfuu(1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud(1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu(1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd(1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
    end do; end do


    do j=1,s%nAtoms; do i=1,s%nAtoms; do gama=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(j)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      index1(1:4) = sigmaimunu2i(1:4,i,mu,nu)
      index2(1:4) = sigmaimunu2i(1:4,j,gama,xi)
      chiorb_hf(index1(1),index2(1)) = chiorb_hf(index1(1),index2(1)) + weight*(gfdd(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(1),index2(2)) = chiorb_hf(index1(1),index2(2)) + weight*(gfdu(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))
      chiorb_hf(index1(1),index2(3)) = chiorb_hf(index1(1),index2(3)) + weight*(gfdd(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(1),index2(4)) = chiorb_hf(index1(1),index2(4)) + weight*(gfdu(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))

      chiorb_hf(index1(2),index2(1)) = chiorb_hf(index1(2),index2(1)) + weight*(gfud(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(2),index2(2)) = chiorb_hf(index1(2),index2(2)) + weight*(gfuu(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))
      chiorb_hf(index1(2),index2(3)) = chiorb_hf(index1(2),index2(3)) + weight*(gfud(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(2),index2(4)) = chiorb_hf(index1(2),index2(4)) + weight*(gfuu(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))

      chiorb_hf(index1(3),index2(1)) = chiorb_hf(index1(3),index2(1)) + weight*(gfdd(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(3),index2(2)) = chiorb_hf(index1(3),index2(2)) + weight*(gfdu(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))
      chiorb_hf(index1(3),index2(3)) = chiorb_hf(index1(3),index2(3)) + weight*(gfdd(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(3),index2(4)) = chiorb_hf(index1(3),index2(4)) + weight*(gfdu(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))

      chiorb_hf(index1(4),index2(1)) = chiorb_hf(index1(4),index2(1)) + weight*(gfud(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(4),index2(2)) = chiorb_hf(index1(4),index2(2)) + weight*(gfuu(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))
      chiorb_hf(index1(4),index2(3)) = chiorb_hf(index1(4),index2(3)) + weight*(gfud(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(4),index2(4)) = chiorb_hf(index1(4),index2(4)) + weight*(gfuu(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))
    end do; end do; end do; end do; end do; end do
    !call zaxpy(dimens*dimens,weight,df1,1,Fint,1)
    !Fint = Fint + df1*weight
  end do !end ix <= pn1 loop
  !$omp end do nowait

  !$omp do schedule(dynamic) reduction(+:Fint)
  do ix2 = 1, real_points ! Third integration (on the real axis)
    nep = (ix2-1) / int(realBZ % workload,8) + 1
    nkp = mod(ix2-1, int(realBZ % workload,8)) + 1
    ep = x2(nep)
    kp = realBZ % kp(:,nkp)
    weight = p2(nep) * realBZ % w(nkp)

    ! Green function at (k+q,Ep+E+i.eta)
    kpq = kp+q
    call calc_green(ep+e,eta,s,kpq,gfq)

    ! Green function at (k,Ep+i.eta)
    call calc_green(ep,etap,s,kp,gf)

    do j=1,s%nAtoms; do i=1,s%nAtoms
      nOrb_i = s%Types(s%Basis(i)%Material)%nOrb
      nOrb_j = s%Types(s%Basis(j)%Material)%nOrb
      nOrb2_i = 2*nOrb_i
      nOrb2_j = 2*nOrb_j

      gfuu(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd(1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)

      gfuu(1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud(1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu(1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd(1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
    end do; end do

    do j=1,s%nAtoms; do i=1,s%nAtoms; do gama=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(j)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      index1(1:4) = sigmaimunu2i(1:4,i,mu,nu)
      index2(1:4) = sigmaimunu2i(1:4,j,gama,xi)
      Fint(index1(1),index2(1)) = Fint(index1(1),index2(1)) - weight * cI*(gfdd(nu,gama,i,j,1) - conjg(gfdd(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(1),index2(2)) = Fint(index1(1),index2(2)) - weight * cI*(gfdu(nu,gama,i,j,1) - conjg(gfud(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(1),index2(3)) = Fint(index1(1),index2(3)) - weight * cI*(gfdd(nu,gama,i,j,1) - conjg(gfdd(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
      Fint(index1(1),index2(4)) = Fint(index1(1),index2(4)) - weight * cI*(gfdu(nu,gama,i,j,1) - conjg(gfud(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

      Fint(index1(2),index2(1)) = Fint(index1(2),index2(1)) - weight * cI*(gfud(nu,gama,i,j,1) - conjg(gfdu(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(2),index2(2)) = Fint(index1(2),index2(2)) - weight * cI*(gfuu(nu,gama,i,j,1) - conjg(gfuu(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(2),index2(3)) = Fint(index1(2),index2(3)) - weight * cI*(gfud(nu,gama,i,j,1) - conjg(gfdu(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
      Fint(index1(2),index2(4)) = Fint(index1(2),index2(4)) - weight * cI*(gfuu(nu,gama,i,j,1) - conjg(gfuu(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

      Fint(index1(3),index2(1)) = Fint(index1(3),index2(1)) - weight * cI*(gfdd(nu,gama,i,j,1) - conjg(gfdd(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(3),index2(2)) = Fint(index1(3),index2(2)) - weight * cI*(gfdu(nu,gama,i,j,1) - conjg(gfud(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(3),index2(3)) = Fint(index1(3),index2(3)) - weight * cI*(gfdd(nu,gama,i,j,1) - conjg(gfdd(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
      Fint(index1(3),index2(4)) = Fint(index1(3),index2(4)) - weight * cI*(gfdu(nu,gama,i,j,1) - conjg(gfud(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))

      Fint(index1(4),index2(1)) = Fint(index1(4),index2(1)) - weight * cI*(gfud(nu,gama,i,j,1) - conjg(gfdu(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(4),index2(2)) = Fint(index1(4),index2(2)) - weight * cI*(gfuu(nu,gama,i,j,1) - conjg(gfuu(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(4),index2(3)) = Fint(index1(4),index2(3)) - weight * cI*(gfud(nu,gama,i,j,1) - conjg(gfdu(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
      Fint(index1(4),index2(4)) = Fint(index1(4),index2(4)) - weight * cI*(gfuu(nu,gama,i,j,1) - conjg(gfuu(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
    end do; end do; end do; end do; end do; end do

    ! Locally add up df1
    ! call zaxpy(dimens*dimens,(p2(ix2)*s%wkbz(iz)),df1,1,Fint,1)
    !Fint = Fint + df1*(p2(ix2)*s%wkbz(iz))
  end do !end pn1+1, nepoints loop
  !$omp end do nowait
  !$omp end parallel

  chiorb_hf = (chiorb_hf + Fint) / tpi

  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1), ierr)

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

  integer(int32), dimension(4) :: index1, index2
  integer(int32) :: i,j,mu,nu,gama,xi,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j
  integer(int64) :: ix,ix2, nep, nkp
  real(dp) :: kp(3), kpq(3), ep
  real(dp) :: weight
  complex(dp), dimension(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms) :: gf,gfq,gvg,gvgq
  complex(dp), dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2) :: gfuu,gfud,gfdu,gfdd
  complex(dp), dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2) :: gvguu,gvgud,gvgdu,gvgdd
  complex(dp), dimension(dimens,dimens) :: Fint,Fintlsoc
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
  Fint = cZero
  Fintlsoc = cZero
  chiorb_hf = cZero
  chiorb_hflsoc = cZero

  !$omp parallel default(none) &
  !$omp& private(ix,ix2,i,j,mu,nu,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j,nep,nkp,gama,xi,kp,kpq,ep,weight,gf,gfq,gfuu,gfud,gfdu,gfdd,gvg,gvgq,gvguu,gvgud,gvgdu,gvgdd,index1,index2) &
  !$omp& shared(local_points,s,bzs,realBZ,real_points,E_k_imag_mesh,q,e,y,wght,x2,p2,eta,etap,dimens,sigmaimunu2i,Fint,Fintlsoc,chiorb_hf,chiorb_hflsoc)

  ! Starting to calculate energy integral
  !$omp do schedule(dynamic) reduction(+:chiorb_hf,chiorb_hflsoc)
  do ix = 1, local_points
    ep = y( E_k_imag_mesh(1,ix) )
    kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
    weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix))

    ! Green function at (k+q,E_F+E+iy)
    kpq = kp+q
    call greenlinearsoc(s%Ef+e,ep+eta,s,kpq,gfq,gvgq)

    ! Green function at (k,E_F+iy)
    call greenlinearsoc(s%Ef,ep+etap,s,kp,gf,gvg)

    do j=1,s%nAtoms; do i=1,s%nAtoms
      nOrb_i = s%Types(s%Basis(i)%Material)%nOrb
      nOrb_j = s%Types(s%Basis(j)%Material)%nOrb
      nOrb2_i = 2*nOrb_i
      nOrb2_j = 2*nOrb_j

      gfuu (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
      gvguu(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gvgud(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gvgdu(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gvgdd(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)

      gfuu (1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud (1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu (1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd (1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
      gvguu(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gvgud(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gvgdu(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gvgdd(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
    end do; end do

    do j=1,s%nAtoms; do i=1,s%nAtoms; do gama=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(j)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      index1(1:4) = sigmaimunu2i(1:4,i,mu,nu)
      index2(1:4) = sigmaimunu2i(1:4,j,gama,xi)
      chiorb_hf(index1(1),index2(1)) = chiorb_hf(index1(1),index2(1)) + weight*(gfdd(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(1),index2(2)) = chiorb_hf(index1(1),index2(2)) + weight*(gfdu(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))
      chiorb_hf(index1(1),index2(3)) = chiorb_hf(index1(1),index2(3)) + weight*(gfdd(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(1),index2(4)) = chiorb_hf(index1(1),index2(4)) + weight*(gfdu(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))

      chiorb_hf(index1(2),index2(1)) = chiorb_hf(index1(2),index2(1)) + weight*(gfud(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(2),index2(2)) = chiorb_hf(index1(2),index2(2)) + weight*(gfuu(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))
      chiorb_hf(index1(2),index2(3)) = chiorb_hf(index1(2),index2(3)) + weight*(gfud(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(2),index2(4)) = chiorb_hf(index1(2),index2(4)) + weight*(gfuu(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))

      chiorb_hf(index1(3),index2(1)) = chiorb_hf(index1(3),index2(1)) + weight*(gfdd(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(3),index2(2)) = chiorb_hf(index1(3),index2(2)) + weight*(gfdu(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))
      chiorb_hf(index1(3),index2(3)) = chiorb_hf(index1(3),index2(3)) + weight*(gfdd(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1)))
      chiorb_hf(index1(3),index2(4)) = chiorb_hf(index1(3),index2(4)) + weight*(gfdu(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfud(gama,nu,j,i,1)))

      chiorb_hf(index1(4),index2(1)) = chiorb_hf(index1(4),index2(1)) + weight*(gfud(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(4),index2(2)) = chiorb_hf(index1(4),index2(2)) + weight*(gfuu(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))
      chiorb_hf(index1(4),index2(3)) = chiorb_hf(index1(4),index2(3)) + weight*(gfud(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1)))
      chiorb_hf(index1(4),index2(4)) = chiorb_hf(index1(4),index2(4)) + weight*(gfuu(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1)))


      chiorb_hflsoc(index1(1),index2(1)) = chiorb_hflsoc(index1(1),index2(1)) + weight*(gvgdd(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1))  +  gfdd(nu,gama,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvgdd(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(1),index2(2)) = chiorb_hflsoc(index1(1),index2(2)) + weight*(gvgdu(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfud(gama,nu,j,i,1))  +  gfdu(nu,gama,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvgud(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(1),index2(3)) = chiorb_hflsoc(index1(1),index2(3)) + weight*(gvgdd(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1))  +  gfdd(nu,gama,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvgdd(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(1),index2(4)) = chiorb_hflsoc(index1(1),index2(4)) + weight*(gvgdu(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfud(gama,nu,j,i,1))  +  gfdu(nu,gama,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvgud(gama,nu,j,i,1)))

      chiorb_hflsoc(index1(2),index2(1)) = chiorb_hflsoc(index1(2),index2(1)) + weight*(gvgud(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1))  +  gfud(nu,gama,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvgdu(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(2),index2(2)) = chiorb_hflsoc(index1(2),index2(2)) + weight*(gvguu(nu,gama,i,j,1)*gfuu(xi,mu,j,i,2) + conjg(gvguu(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1))  +  gfuu(nu,gama,i,j,1)*gvguu(xi,mu,j,i,2) + conjg(gfuu(mu,xi,i,j,2)*gvguu(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(2),index2(3)) = chiorb_hflsoc(index1(2),index2(3)) + weight*(gvgud(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1))  +  gfud(nu,gama,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvgdu(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(2),index2(4)) = chiorb_hflsoc(index1(2),index2(4)) + weight*(gvguu(nu,gama,i,j,1)*gfdu(xi,mu,j,i,2) + conjg(gvgud(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1))  +  gfuu(nu,gama,i,j,1)*gvgdu(xi,mu,j,i,2) + conjg(gfud(mu,xi,i,j,2)*gvguu(gama,nu,j,i,1)))

      chiorb_hflsoc(index1(3),index2(1)) = chiorb_hflsoc(index1(3),index2(1)) + weight*(gvgdd(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1))  +  gfdd(nu,gama,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvgdd(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(3),index2(2)) = chiorb_hflsoc(index1(3),index2(2)) + weight*(gvgdu(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfud(gama,nu,j,i,1))  +  gfdu(nu,gama,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvgud(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(3),index2(3)) = chiorb_hflsoc(index1(3),index2(3)) + weight*(gvgdd(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfdd(gama,nu,j,i,1))  +  gfdd(nu,gama,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvgdd(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(3),index2(4)) = chiorb_hflsoc(index1(3),index2(4)) + weight*(gvgdu(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfud(gama,nu,j,i,1))  +  gfdu(nu,gama,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvgud(gama,nu,j,i,1)))

      chiorb_hflsoc(index1(4),index2(1)) = chiorb_hflsoc(index1(4),index2(1)) + weight*(gvgud(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1))  +  gfud(nu,gama,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvgdu(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(4),index2(2)) = chiorb_hflsoc(index1(4),index2(2)) + weight*(gvguu(nu,gama,i,j,1)*gfud(xi,mu,j,i,2) + conjg(gvgdu(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1))  +  gfuu(nu,gama,i,j,1)*gvgud(xi,mu,j,i,2) + conjg(gfdu(mu,xi,i,j,2)*gvguu(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(4),index2(3)) = chiorb_hflsoc(index1(4),index2(3)) + weight*(gvgud(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfdu(gama,nu,j,i,1))  +  gfud(nu,gama,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvgdu(gama,nu,j,i,1)))
      chiorb_hflsoc(index1(4),index2(4)) = chiorb_hflsoc(index1(4),index2(4)) + weight*(gvguu(nu,gama,i,j,1)*gfdd(xi,mu,j,i,2) + conjg(gvgdd(mu,xi,i,j,2)*gfuu(gama,nu,j,i,1))  +  gfuu(nu,gama,i,j,1)*gvgdd(xi,mu,j,i,2) + conjg(gfdd(mu,xi,i,j,2)*gvguu(gama,nu,j,i,1)))
    end do; end do; end do; end do; end do; end do

    ! Fint = Fint + df1 * weight
    ! Fintlsoc = Fintlsoc + df1lsoc * weight

  end do
  !$omp end do nowait

  !$omp do schedule(dynamic) reduction(+:Fint,Fintlsoc)
  do ix2 = 1, realBZ%workload ! Third integration (on the real axis)
    nep = (ix2-1) / int(realBZ % workload,8) + 1
    nkp = mod(ix2-1, int(realBZ % workload,8))+1
    ep = x2(nep)
    kp = realBZ % kp(:,nkp)
    weight = p2(nep) * realBZ % w(nkp)

    ! Green function at (k+q,E_F+E+iy)
    kpq = kp+q
    call greenlinearsoc(ep+e,eta,s,kpq,gfq,gvgq)

    ! Green function at (k,E_F+iy)
    call greenlinearsoc(ep,etap,s,kp,gf,gvg)

    do j=1,s%nAtoms; do i=1,s%nAtoms
      nOrb_i = s%Types(s%Basis(i)%Material)%nOrb
      nOrb_j = s%Types(s%Basis(j)%Material)%nOrb
      nOrb2_i = 2*nOrb_i
      nOrb2_j = 2*nOrb_j

      gfuu (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd (1:nOrb_i,1:nOrb_j,i,j,1) = gfq(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
      gvguu(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gvgud(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gvgdu(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gvgdd(1:nOrb_i,1:nOrb_j,i,j,1) = gvgq(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)

      gfuu (1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gfud (1:nOrb_i,1:nOrb_j,i,j,2) = gf(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gfdu (1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gfdd (1:nOrb_i,1:nOrb_j,i,j,2) = gf(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
      gvguu(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(       1:nOrb_i ,       1:nOrb_j ,i,j)
      gvgud(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(       1:nOrb_i ,nOrb_j+1:nOrb2_j,i,j)
      gvgdu(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(nOrb_i+1:nOrb2_i,       1:nOrb_j ,i,j)
      gvgdd(1:nOrb_i,1:nOrb_j,i,j,2) = gvg(nOrb_i+1:nOrb2_i,nOrb_j+1:nOrb2_j,i,j)
    end do; end do

    do j=1,s%nAtoms; do i=1,s%nAtoms; do gama=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(j)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      index1(1:4) = sigmaimunu2i(1:4,i,mu,nu)
      index2(1:4) = sigmaimunu2i(1:4,j,gama,xi)
      Fint(index1(1),index2(1)) = Fint(index1(1),index2(1)) - weight * cI*(gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(1),index2(2)) = Fint(index1(1),index2(2)) - weight * cI*(gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(1),index2(3)) = Fint(index1(1),index2(3)) - weight * cI*(gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
      Fint(index1(1),index2(4)) = Fint(index1(1),index2(4)) - weight * cI*(gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

      Fint(index1(2),index2(1)) = Fint(index1(2),index2(1)) - weight * cI*(gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(2),index2(2)) = Fint(index1(2),index2(2)) - weight * cI*(gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))
      Fint(index1(2),index2(3)) = Fint(index1(2),index2(3)) - weight * cI*(gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))
      Fint(index1(2),index2(4)) = Fint(index1(2),index2(4)) - weight * cI*(gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))

      Fint(index1(3),index2(1)) = Fint(index1(3),index2(1)) - weight * cI*(gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(3),index2(2)) = Fint(index1(3),index2(2)) - weight * cI*(gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(3),index2(3)) = Fint(index1(3),index2(3)) - weight * cI*(gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
      Fint(index1(3),index2(4)) = Fint(index1(3),index2(4)) - weight * cI*(gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))

      Fint(index1(4),index2(1)) = Fint(index1(4),index2(1)) - weight * cI*(gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(4),index2(2)) = Fint(index1(4),index2(2)) - weight * cI*(gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))
      Fint(index1(4),index2(3)) = Fint(index1(4),index2(3)) - weight * cI*(gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))
      Fint(index1(4),index2(4)) = Fint(index1(4),index2(4)) - weight * cI*(gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))


      Fintlsoc(index1(1),index2(1)) = Fintlsoc(index1(1),index2(1)) - weight * cI*((gvgdd(nu,gama,i,j,1)-conjg(gvgdd(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
      Fintlsoc(index1(1),index2(2)) = Fintlsoc(index1(1),index2(2)) - weight * cI*((gvgdu(nu,gama,i,j,1)-conjg(gvgud(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
      Fintlsoc(index1(1),index2(3)) = Fintlsoc(index1(1),index2(3)) - weight * cI*((gvgdd(nu,gama,i,j,1)-conjg(gvgdd(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))
      Fintlsoc(index1(1),index2(4)) = Fintlsoc(index1(1),index2(4)) - weight * cI*((gvgdu(nu,gama,i,j,1)-conjg(gvgud(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))

      Fintlsoc(index1(2),index2(1)) = Fintlsoc(index1(2),index2(1)) - weight * cI*((gvgud(nu,gama,i,j,1)-conjg(gvgdu(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
      Fintlsoc(index1(2),index2(2)) = Fintlsoc(index1(2),index2(2)) - weight * cI*((gvguu(nu,gama,i,j,1)-conjg(gvguu(gama,nu,j,i,1)))*conjg(gfuu(mu,xi,i,j,2))  +  (gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gvguu(mu,xi,i,j,2)))
      Fintlsoc(index1(2),index2(3)) = Fintlsoc(index1(2),index2(3)) - weight * cI*((gvgud(nu,gama,i,j,1)-conjg(gvgdu(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))
      Fintlsoc(index1(2),index2(4)) = Fintlsoc(index1(2),index2(4)) - weight * cI*((gvguu(nu,gama,i,j,1)-conjg(gvguu(gama,nu,j,i,1)))*conjg(gfud(mu,xi,i,j,2))  +  (gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gvgud(mu,xi,i,j,2)))

      Fintlsoc(index1(3),index2(1)) = Fintlsoc(index1(3),index2(1)) - weight * cI*((gvgdd(nu,gama,i,j,1)-conjg(gvgdd(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
      Fintlsoc(index1(3),index2(2)) = Fintlsoc(index1(3),index2(2)) - weight * cI*((gvgdu(nu,gama,i,j,1)-conjg(gvgud(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
      Fintlsoc(index1(3),index2(3)) = Fintlsoc(index1(3),index2(3)) - weight * cI*((gvgdd(nu,gama,i,j,1)-conjg(gvgdd(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfdd(nu,gama,i,j,1)-conjg(gfdd(gama,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
      Fintlsoc(index1(3),index2(4)) = Fintlsoc(index1(3),index2(4)) - weight * cI*((gvgdu(nu,gama,i,j,1)-conjg(gvgud(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfdu(nu,gama,i,j,1)-conjg(gfud(gama,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))

      Fintlsoc(index1(4),index2(1)) = Fintlsoc(index1(4),index2(1)) - weight * cI*((gvgud(nu,gama,i,j,1)-conjg(gvgdu(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
      Fintlsoc(index1(4),index2(2)) = Fintlsoc(index1(4),index2(2)) - weight * cI*((gvguu(nu,gama,i,j,1)-conjg(gvguu(gama,nu,j,i,1)))*conjg(gfdu(mu,xi,i,j,2))  +  (gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,j,2)))
      Fintlsoc(index1(4),index2(3)) = Fintlsoc(index1(4),index2(3)) - weight * cI*((gvgud(nu,gama,i,j,1)-conjg(gvgdu(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfud(nu,gama,i,j,1)-conjg(gfdu(gama,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
      Fintlsoc(index1(4),index2(4)) = Fintlsoc(index1(4),index2(4)) - weight * cI*((gvguu(nu,gama,i,j,1)-conjg(gvguu(gama,nu,j,i,1)))*conjg(gfdd(mu,xi,i,j,2))  +  (gfuu(nu,gama,i,j,1)-conjg(gfuu(gama,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,j,2)))
    end do; end do; end do; end do; end do; end do
  end do
  !$omp end do nowait
  !$omp end parallel

  chiorb_hf = (chiorb_hf + Fint) / tpi
  chiorb_hflsoc = (chiorb_hflsoc + Fintlsoc) / tpi
  
  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hf    , ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1), ierr)
  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hflsoc, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1), ierr)

end subroutine eintshechilinearsoc
