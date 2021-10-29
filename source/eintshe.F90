! ---------- Parallel spin current: Energy integration ---------
subroutine eintshe(q,e)
  use mod_kind,          only: dp,int64
  use mod_constants,     only: cZero,cI,tpi
  use mod_parameters,    only: dimens,sigmaimunu2i,eta,etap,sigmai2i
  use mod_SOC,           only: llineargfsoc
  use EnergyIntegration, only: y, wght, x2, p2, generate_real_epoints, pn2
  use mod_system,        only: s => sys
  use mod_BrillouinZone, only: realBZ
  use mod_prefactors,    only: prefactor !, lxpt, lypt, lzpt, tlxp, tlyp, tlzp
  use mod_disturbances,  only: tchiorbiikl
  use ElectricField,     only: ElectricFieldVector
  use mod_hamiltonian,   only: hamilt_local
  use mod_greenfunction, only: calc_green
  use adaptiveMesh,      only: bzs,E_k_imag_mesh,local_points
  use mod_mpi_pars,      only: rFreq,MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr,abortProgram
  ! use mod_currents,      only: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl !TODO: Re-Include
  ! use mod_system,        only: n0sc1, n0sc2, n0sc, r_nn, npln
  implicit none

  real(dp),intent(in)   :: e,q(3)

  complex(dp), dimension(dimens,4) :: tFintiikl1,tFintiikl2
  !complex(dp), dimension(n0sc1:n0sc2,dimspinAtoms,4) :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl TODO:Re-Include

  integer :: i,j,l,mu,nu,gamma,xi,sigma,sigmap,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j !,neighbor
  real(dp) :: kp(3),kpq(3),ep
  complex(dp) :: wkbzc
  !complex(dp) :: expikr(n0sc1:n0sc2) TODO:Re-Include
  complex(dp),dimension(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms) :: gf,gfq
  complex(dp),dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms)   :: dtdk
  complex(dp),dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2) :: gfuu,gfud,gfdu,gfdd
  complex(dp),dimension(dimens,4) :: df1iikl,pfdf1iikl
  !complex(dp),dimension(n0sc1:n0sc2,s%nAtoms,s%nOrb, s%nOrb)    :: prett,preLxtt,preLytt,preLztt !TODO: Re-Include
  !complex(dp),dimension(n0sc1:n0sc2,dimspinAtoms,4), intent(out) :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl !TODO:Re-Include
  !--------------------- begin MPI vars --------------------
  integer(int64) :: ix, ix2, nep,nkp
  integer(int64) :: real_points
  integer   :: ncountkl !,nncountkl !TODO: Re-Include

  external :: dtdksub,zgemm,MPI_Reduce

  ncountkl = dimens*4
  !nncountkl = n0sc*dimspinAtoms*4 !TODO: Re-Include
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  real_points = 0
  if(abs(e) >= 1.e-15_dp) real_points = int(pn2*realBZ%workload,8)

  ! Starting to calculate energy integral
  ! ttchiorbiikl    = cZero !TODO: Re-Include
  ! Lxttchiorbiikl  = cZero !TODO: Re-Include
  ! Lyttchiorbiikl  = cZero !TODO: Re-Include
  ! Lzttchiorbiikl  = cZero !TODO: Re-Include

  tFintiikl1    = cZero
  tFintiikl2    = cZero
  !ttFintiikl   = cZero !TODO: Re-Include
  !LxttFintiikl = cZero !TODO: Re-Include
  !LyttFintiikl = cZero !TODO: Re-Include
  !LzttFintiikl = cZero !TODO: Re-Include


  ! Build local hamiltonian
  if(.not.llineargfsoc) call hamilt_local(s)

  !$omp parallel default(none) &
  !$omp& private(ix,ix2,wkbzc,ep,kp,kpq,nep,nkp,i,j,l,mu,nu,gamma,xi,sigma,sigmap,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j,dtdk,gf,gfq,gfuu,gfud,gfdu,gfdd,df1iikl,pfdf1iikl) & !,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,expikr,prett,preLxtt,preLytt,preLztt) &
  !$omp& shared(s,calc_green,bzs,realBZ,E_k_imag_mesh,ElectricFieldVector,prefactor,q,e,y,x2,wght,p2,pn2,eta,etap,local_points,real_points,sigmai2i,sigmaimunu2i,dimens,tchiorbiikl,tFintiikl1,tFintiikl2) !,n0sc1,n0sc2,lxpt,lypt,lzpt,tlxp,tlyp,tlzp)

  !$omp do schedule(dynamic) reduction(+:tFintiikl1)
  do ix = 1, local_points
    ep = y( E_k_imag_mesh(1,ix) )
    kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
    wkbzc = cmplx(wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix)),0._dp,dp)

    df1iikl = cZero

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(ElectricFieldVector, kp, dtdk)

    ! TODO: Fix currents
    !do neighbor = n0sc1,n0sc2
    !  expikr(neighbor) = exp(cI * dot_product(kp, r_nn(:,neighbor)))
    !end do

    ! Calculating the prefactor (L).t.exp - t.(L).exp
    ! TODO: Fix Currents
    ! do nu=1,s%nOrb
    !   do mu=1,s%nOrb
    !     do neighbor=n0sc1,n0sc2
    !       do i = 1, s%nAtoms
    !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i)*expikr(neighbor)-(t0i(mu,nu,neighbor,i)*conjg(expikr(neighbor)))
    !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !       end do
    !     end do
    !   end do
    ! end do

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

    do j=1,s%nAtoms; do i=1,s%nAtoms; do l=1,s%nAtoms; do gamma=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(l)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      ! if(abs(j-l) > npln) cycle
      df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1) + (gfdd(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2) + (gfdu(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3) + (gfdd(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4) + (gfdu(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1) + (gfud(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2) + (gfuu(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3) + (gfud(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4) + (gfuu(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1) + (gfdd(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2) + (gfdu(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3) + (gfdd(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4) + (gfdu(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1) + (gfud(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2) + (gfuu(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3) + (gfud(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4) + (gfuu(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
    end do; end do; end do; end do; end do; end do; end do

    ! Multiplying the prefactor by the susceptibility and the k-point weight
    ! pfdf1iikl = prefactor*df1iikl
    call zgemm('n','n',dimens,4,dimens,wkbzc,prefactor,dimens,df1iikl,dimens,cZero,pfdf1iikl,dimens)

    tFintiikl1 = tFintiikl1 + pfdf1iikl

    ! TODO: Fix currents
    ! do sigmap = 1, 4
    !   do nu = 1, s%nOrb
    !     do mu = 1, s%nOrb
    !       do i = 1, s%nAtoms
    !         do sigma = 1, 4
    !           do neighbor = n0sc1, n0sc2
    !             ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) = ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) + (prett  (neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             ! Orbital angular momentum currents
    !             LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLxtt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLytt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLztt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !           end do
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do

  end do
  !$omp end do nowait

  !$omp do schedule(dynamic) reduction(+:tFintiikl2)
  do ix2 = 1, real_points ! Third integration (on the real axis)
    nep = (ix2-1) / realBZ % workload + 1
    nkp = mod(ix2-1, int(realBZ % workload,8))+1
    ep = x2(nep)
    kp = realBZ % kp(:,nkp)
    wkbzc = cmplx(p2(nep) * realBZ % w(nkp), 0._dp,dp)

    df1iikl = cZero

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(ElectricFieldVector, kp, dtdk)

    ! TODO: Fix currents
    !do neighbor = n0sc1,n0sc2
    !  expikr(neighbor) = exp(cI * dot_product(kp, r_nn(:,neighbor)))
    !end do

    ! Calculating the prefactor (L).t.exp - t.(L).exp
    ! TODO: Fix Currents
    ! do nu=1,s%nOrb
    !   do mu=1,s%nOrb
    !     do neighbor=n0sc1,n0sc2
    !       do i = 1, s%nAtoms
    !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i)*expikr(neighbor)-(t0i(mu,nu,neighbor,i)*conjg(expikr(neighbor)))
    !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !       end do
    !     end do
    !   end do
    ! end do

    ! Green function at (k+q,E'+E+i.eta)
    kpq = kp+q
    call calc_green(ep+e,eta,s,kpq,gfq)

    ! Green function at (k,E'+i.eta)
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

    do j=1,s%nAtoms; do i=1,s%nAtoms; do l=1,s%nAtoms; do gamma=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(l)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
    end do; end do; end do; end do; end do; end do; end do

    ! Multiplying the prefactor by the susceptibility and the k-point weight
    ! pfdf1iikl = prefactor*df1iikl
    call zgemm('n','n',dimens,4,dimens,wkbzc,prefactor,dimens,df1iikl,dimens,cZero,pfdf1iikl,dimens)

    ! Integrating matrices
    tFintiikl2 = tFintiikl2 + pfdf1iikl

    ! TODO: Fix currents
    ! do sigmap=1,4
    !   do nu=1,s%nOrb
    !     do mu=1,s%nOrb
    !       do i=1,Npl
    !         do sigma=1,4
    !           do neighbor=n0sc1,n0sc2
    !             ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) = ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) + (prett  (neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             ! Orbital angular momentum currents
    !             LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLxtt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLytt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLztt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !           end do
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do
  end do
  !$omp end do nowait
  !$omp end parallel

  tchiorbiikl = (tFintiikl1+tFintiikl2)/tpi
  !ttchiorbiikl = ttchiorbiikl + ttFintiikl
  !Lxttchiorbiikl = Lxttchiorbiikl + LxttFintiikl
  !Lyttchiorbiikl = Lyttchiorbiikl + LyttFintiikl
  !Lzttchiorbiikl = Lzttchiorbiikl + LzttFintiikl

  if(rFreq(1) == 0) then
    call MPI_Reduce(MPI_IN_PLACE, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr)
    !call MPI_Reduce(MPI_IN_PLACE, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
  else
    call MPI_Reduce(tchiorbiikl, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr)
    !call MPI_Reduce(ttchiorbiikl, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(Lxttchiorbiikl, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(Lyttchiorbiikl, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(Lzttchiorbiikl, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
  end if

end subroutine eintshe


subroutine eintshelinearsoc(q,e)
  use mod_kind,          only: dp,int64
  use mod_constants,     only: cZero,cOne,cI,tpi
  use mod_parameters,    only: dimens,eta,etap,sigmai2i,sigmaimunu2i
  use mod_system,        only: s => sys
  use mod_BrillouinZone, only: realBZ
  use EnergyIntegration, only: y,wght,x2,p2,generate_real_epoints,pn2
  use mod_prefactors,    only: prefactor,prefactorlsoc
  use mod_disturbances,  only: tchiorbiikl
  use ElectricField,     only: ElectricFieldVector
  use mod_greenfunction, only: greenlinearsoc
  use adaptiveMesh,      only: bzs,E_k_imag_mesh,local_points
  use mod_mpi_pars,      only: rFreq,MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr,abortProgram

  !use mod_system, only: r_nn, npln, nkpt, kbz, wkbz, n0sc1, n0sc2
  !use mod_currents,         only: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl !TODO: Re-Include
  !use mod_system,           only: n0sc1, n0sc2, n0sc
  implicit none
  real(dp),intent(in) :: e,q(3)

  complex(dp), dimension(dimens,4) :: tFintiikl1,tFintiikl2
  !complex(dp), dimension(n0sc1:n0sc2,dimspinAtoms,4) :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl TODO:Re-Include

  integer :: i,j,l,mu,nu,gamma,xi,sigma,sigmap,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j !,neighbor
  real(dp) :: kp(3), kpq(3), ep
  integer(int64) :: nep, nkp
  !complex(dp) :: expikr(n0sc1:n0sc2) !TODO: Re-Include
  complex(dp) :: wkbzc

  complex(dp),dimension(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms) :: gf,gfq,gvg,gvgq
  complex(dp),dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms)   :: dtdk
  complex(dp),dimension(s%nOrb,s%nOrb,s%nAtoms,s%nAtoms,2) :: gfuu,gfud,gfdu,gfdd,gvguu,gvgud,gvgdu,gvgdd
  complex(dp),dimension(dimens,4) :: df1iikl,pfdf1iikl,df1lsoc
  !complex(dp),dimension(n0sc1:n0sc2,Npl,s%nOrb,s%nOrb)    :: prett,preLxtt,preLytt,preLztt !TODO: Re-Include

  integer(int64) :: ix,ix2, iz
  integer(int64) :: real_points
  integer :: ncountkl !,nncountkl !TODO: Re-Include

  external :: dtdksub,zgemm,MPI_Reduce

  ncountkl = dimens*4
  !nncountkl = n0sc*dimspinAtoms*4 !TODO: Re-Include

  ! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  real_points = 0
  if(abs(e) >= 1.e-15_dp) real_points = int(pn2*realBZ%workload,8)

  ! Starting to calculate energy integral
  ! ttchiorbiikl    = cZero !TODO: Re-Include
  ! Lxttchiorbiikl  = cZero !TODO: Re-Include
  ! Lyttchiorbiikl  = cZero !TODO: Re-Include
  ! Lzttchiorbiikl  = cZero !TODO: Re-Include

  tFintiikl1    = cZero
  tFintiikl2    = cZero
  ! ttFintiikl   = cZero !TODO: Re-Include
  ! LxttFintiikl = cZero !TODO: Re-Include
  ! LyttFintiikl = cZero !TODO: Re-Include
  ! LzttFintiikl = cZero !TODO: Re-Include


  !$omp parallel default(none) &
  !$omp& private(iz,wkbzc,kp,kpq,ep,nep,nkp,df1iikl,pfdf1iikl,df1lsoc,dtdk,gf,gfq,gfuu,gfud,gfdu,gfdd,gvg,gvgq,gvguu,gvgud,gvgdu,gvgdd,sigma,sigmap,i,j,l,mu,nu,gamma,xi,nOrb_i,nOrb_j,nOrb2_i,nOrb2_j) &    !,expikr,prett,preLxtt,preLytt,preLztt
  !$omp& shared(local_points,real_points,s,bzs,E_k_imag_mesh,ElectricFieldVector,tchiorbiikl,prefactor,prefactorlsoc,realBZ,q,e,y,x2,wght,p2,eta,etap,sigmai2i,sigmaimunu2i,dimens,tFintiikl1,tFintiikl2)                                                                                 !,n0sc1,n0sc2,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,lxpt,lypt,lzpt,tlxp,tlyp,tlzp

  !$omp do schedule(dynamic) reduction(+:tFintiikl1)
  do ix = 1, local_points
    ep = y( E_k_imag_mesh(1,ix) )
    kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
    wkbzc = cmplx(wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix)),0._dp,dp)

    df1iikl = cZero
    df1lsoc = cZero

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(ElectricFieldVector, kp,dtdk)

    ! do neighbor=n0sc1,n0sc2 !TODO: Re-Include
    !   expikr(neighbor) = exp(cI*dot_product(kp, r_nn(:,neighbor)))
    ! end do

    ! Calculating the prefactor (L).t.exp - t.(L).exp
    ! do nu=1,s%nOrb !TODO: Re-Include
    !   do mu=1,s%nOrb
    !     do i = 1, s%nAtoms
    !       do neighbor=n0sc1,n0sc2
    !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i)*expikr(neighbor)-(t0i(mu,nu,neighbor,i)*conjg(expikr(neighbor)))
    !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !       end do
    !     end do
    !   end do
    ! end do

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

    do j=1,s%nAtoms;  do i=1,s%nAtoms; do l=1,s%nAtoms; do gamma=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(l)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      ! if(abs(j-l) > npln) cycle
      df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1) + (gfdd(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2) + (gfdu(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3) + (gfdd(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4) + (gfdu(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1) + (gfud(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2) + (gfuu(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3) + (gfud(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4) + (gfuu(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1) + (gfdd(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2) + (gfdu(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3) + (gfdd(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4) + (gfdu(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1) + (gfud(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2) + (gfuu(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3) + (gfud(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4) + (gfuu(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)


      df1lsoc(sigmaimunu2i(1,i,mu,nu),1) = df1lsoc(sigmaimunu2i(1,i,mu,nu),1) + (gvgdd(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gvguu(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvguu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gvgdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(1,i,mu,nu),2) = df1lsoc(sigmaimunu2i(1,i,mu,nu),2) + (gvgdu(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gvguu(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvguu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gvgud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(1,i,mu,nu),3) = df1lsoc(sigmaimunu2i(1,i,mu,nu),3) + (gvgdd(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gvgud(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvgdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gvgdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(1,i,mu,nu),4) = df1lsoc(sigmaimunu2i(1,i,mu,nu),4) + (gvgdu(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gvgud(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvgdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gvgud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1lsoc(sigmaimunu2i(2,i,mu,nu),1) = df1lsoc(sigmaimunu2i(2,i,mu,nu),1) + (gvgud(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gvguu(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvguu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gvgdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(2,i,mu,nu),2) = df1lsoc(sigmaimunu2i(2,i,mu,nu),2) + (gvguu(nu,gamma,i,j,1)*gfuu(xi,mu,l,i,2) + conjg(gvguu(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvguu(xi,mu,l,i,2) + conjg(gfuu(mu,xi,i,l,2)*gvguu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(2,i,mu,nu),3) = df1lsoc(sigmaimunu2i(2,i,mu,nu),3) + (gvgud(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gvgud(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvgdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gvgdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(2,i,mu,nu),4) = df1lsoc(sigmaimunu2i(2,i,mu,nu),4) + (gvguu(nu,gamma,i,j,1)*gfdu(xi,mu,l,i,2) + conjg(gvgud(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvgdu(xi,mu,l,i,2) + conjg(gfud(mu,xi,i,l,2)*gvguu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1lsoc(sigmaimunu2i(3,i,mu,nu),1) = df1lsoc(sigmaimunu2i(3,i,mu,nu),1) + (gvgdd(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gvgdu(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvgud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gvgdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(3,i,mu,nu),2) = df1lsoc(sigmaimunu2i(3,i,mu,nu),2) + (gvgdu(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gvgdu(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvgud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gvgud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(3,i,mu,nu),3) = df1lsoc(sigmaimunu2i(3,i,mu,nu),3) + (gvgdd(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gvgdd(mu,xi,i,l,2)*gfdd(gamma,nu,j,i,1))  +  gfdd(nu,gamma,i,j,1)*gvgdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gvgdd(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(3,i,mu,nu),4) = df1lsoc(sigmaimunu2i(3,i,mu,nu),4) + (gvgdu(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gvgdd(mu,xi,i,l,2)*gfud(gamma,nu,j,i,1))  +  gfdu(nu,gamma,i,j,1)*gvgdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gvgud(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)

      df1lsoc(sigmaimunu2i(4,i,mu,nu),1) = df1lsoc(sigmaimunu2i(4,i,mu,nu),1) + (gvgud(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gvgdu(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvgud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gvgdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(4,i,mu,nu),2) = df1lsoc(sigmaimunu2i(4,i,mu,nu),2) + (gvguu(nu,gamma,i,j,1)*gfud(xi,mu,l,i,2) + conjg(gvgdu(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvgud(xi,mu,l,i,2) + conjg(gfdu(mu,xi,i,l,2)*gvguu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(4,i,mu,nu),3) = df1lsoc(sigmaimunu2i(4,i,mu,nu),3) + (gvgud(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gvgdd(mu,xi,i,l,2)*gfdu(gamma,nu,j,i,1))  +  gfud(nu,gamma,i,j,1)*gvgdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gvgdu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(4,i,mu,nu),4) = df1lsoc(sigmaimunu2i(4,i,mu,nu),4) + (gvguu(nu,gamma,i,j,1)*gfdd(xi,mu,l,i,2) + conjg(gvgdd(mu,xi,i,l,2)*gfuu(gamma,nu,j,i,1))  +  gfuu(nu,gamma,i,j,1)*gvgdd(xi,mu,l,i,2) + conjg(gfdd(mu,xi,i,l,2)*gvguu(gamma,nu,j,i,1)))*dtdk(gamma,xi,j,l)
    end do; end do; end do; end do; end do; end do; end do


    ! Multiplying the prefactor by the susceptibility and the k-point weight
    ! pfdf1iikl = prefactorlsoc*df1iikl
    call zgemm('n','n',dimens,4,dimens,wkbzc,prefactorlsoc,dimens,df1iikl,dimens,cZero,pfdf1iikl,dimens)
    ! pfdf1lsoc = prefactorlsoc*df1iikl + prefactor*df1lsoc
    call zgemm('n','n',dimens,4,dimens,wkbzc,prefactor,dimens,df1lsoc,dimens,cOne,pfdf1iikl,dimens)


    ! Integrating matrices
    tFintiikl1 = tFintiikl1 + pfdf1iikl

    ! do sigmap=1,4 !TODO: Re-Include
    !   do i = 1, s%nAtoms
    !     do sigma=1,4
    !       do nu=1,s%nOrb
    !         do mu=1,s%nOrb
    !           do neighbor = n0sc1, n0sc2
    !             ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) = ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) + (prett  (neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             ! Orbital angular momentum currents
    !             LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLxtt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLytt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLztt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !           end do
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do

  end do
  !$omp end do nowait

  !$omp do schedule(dynamic) reduction(+:tFintiikl2)
  do ix2 = 1, realBZ%workload ! Third integration (on the real axis)
    nep = (ix2-1) / realBZ % workload + 1
    nkp = mod(ix2-1, int(realBZ % workload,8))+1
    ep = x2(nep)
    kp = realBZ % kp(:,nkp)
    wkbzc = cmplx(p2(nep) * realBZ % w(nkp),0._dp,dp)

    df1iikl = cZero
    df1lsoc = cZero

    ! Calculating derivative of in-plane and n.n. inter-plane hoppings
    call dtdksub(ElectricFieldVector, kp,dtdk)

    ! do neighbor=n0sc1,n0sc2 !TODO: Re-Include
    !   expikr(neighbor) = exp(cI*dot_product(kp, r_nn(:,neighbor)))
    ! end do

    ! Calculating the prefactor (L).t.exp - t.(L).exp
    ! do nu=1,s%nOrb !TODO: Re-Include
    !   do mu=1,s%nOrb
    !     do neighbor=n0sc1,n0sc2
    !       do i = 1, s%nAtoms
    !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i)*expikr(neighbor)-(t0i(mu,nu,neighbor,i)*conjg(expikr(neighbor)))
    !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
    !       end do
    !     end do
    !   end do
    ! end do

    ! Green function at (k+q,E'+E+i.eta)
    kpq = kp+q
    call greenlinearsoc(ep+e,eta,s,kpq,gfq,gvgq)

    ! Green function at (k,E'+i.eta)
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

    do j=1,s%nAtoms;  do i=1,s%nAtoms; do l=1,s%nAtoms; do gamma=1,s%Types(s%Basis(j)%Material)%nOrb; do mu=1,s%Types(s%Basis(i)%Material)%nOrb; do xi=1,s%Types(s%Basis(l)%Material)%nOrb; do nu=1,s%Types(s%Basis(i)%Material)%nOrb
      !if(abs(j-l) > npln) cycle
      df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3)-cI*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4)-cI*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

      df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3)-cI*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
      df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4)-cI*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)


      df1lsoc(sigmaimunu2i(1,i,mu,nu),1) = df1lsoc(sigmaimunu2i(1,i,mu,nu),1)-cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(1,i,mu,nu),2) = df1lsoc(sigmaimunu2i(1,i,mu,nu),2)-cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(1,i,mu,nu),3) = df1lsoc(sigmaimunu2i(1,i,mu,nu),3)-cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(1,i,mu,nu),4) = df1lsoc(sigmaimunu2i(1,i,mu,nu),4)-cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)

      df1lsoc(sigmaimunu2i(2,i,mu,nu),1) = df1lsoc(sigmaimunu2i(2,i,mu,nu),1)-cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(2,i,mu,nu),2) = df1lsoc(sigmaimunu2i(2,i,mu,nu),2)-cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(2,i,mu,nu),3) = df1lsoc(sigmaimunu2i(2,i,mu,nu),3)-cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(2,i,mu,nu),4) = df1lsoc(sigmaimunu2i(2,i,mu,nu),4)-cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)

      df1lsoc(sigmaimunu2i(3,i,mu,nu),1) = df1lsoc(sigmaimunu2i(3,i,mu,nu),1)-cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(3,i,mu,nu),2) = df1lsoc(sigmaimunu2i(3,i,mu,nu),2)-cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(3,i,mu,nu),3) = df1lsoc(sigmaimunu2i(3,i,mu,nu),3)-cI*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(3,i,mu,nu),4) = df1lsoc(sigmaimunu2i(3,i,mu,nu),4)-cI*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)

      df1lsoc(sigmaimunu2i(4,i,mu,nu),1) = df1lsoc(sigmaimunu2i(4,i,mu,nu),1)-cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(4,i,mu,nu),2) = df1lsoc(sigmaimunu2i(4,i,mu,nu),2)-cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(4,i,mu,nu),3) = df1lsoc(sigmaimunu2i(4,i,mu,nu),3)-cI*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
      df1lsoc(sigmaimunu2i(4,i,mu,nu),4) = df1lsoc(sigmaimunu2i(4,i,mu,nu),4)-cI*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
    end do; end do; end do; end do; end do; end do; end do

    ! Multiplying the prefactor by the susceptibility and the k-point weight
    ! pfdf1iikl = prefactorlsoc*df1iikl
    call zgemm('n','n',dimens,4,dimens,wkbzc,prefactorlsoc,dimens,df1iikl,dimens,cZero,pfdf1iikl,dimens)
    ! pfdf1lsoc = prefactorlsoc*df1iikl + prefactor*df1lsoc
    call zgemm('n','n',dimens,4,dimens,wkbzc,prefactor,dimens,df1lsoc,dimens,cOne,pfdf1iikl,dimens)

    ! Integrating matrices
    tFintiikl2 = tFintiikl2 + pfdf1iikl

    ! do sigmap=1,4 !TODO: Re-Include
    !   do i = 1, s%nAtoms
    !     do sigma = 1, 4
    !       do nu=1,s%nOrb
    !         do mu=1,s%nOrb
    !           do neighbor = n0sc1, n0sc2
    !             ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) = ttFintiikl  (neighbor,sigmai2i(sigma,i),sigmap) + (prett  (neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             ! Orbital angular momentum currents
    !             LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LxttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLxtt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LyttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLytt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !             LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) = LzttFintiikl(neighbor,sigmai2i(sigma,i),sigmap) + (preLztt(neighbor,i,mu,nu)*pfdf1iikl(sigmaimunu2i(sigma,i,mu,nu),sigmap))
    !           end do
    !         end do
    !       end do
    !     end do
    !   end do
    ! end do

  end do
  !$omp end do nowait
  !$omp end parallel

  tchiorbiikl = (tFintiikl1+tFintiikl2)/tpi
  ! ttchiorbiikl = ttchiorbiikl + ttFintiikl !TODO: Re-Include
  ! Lxttchiorbiikl = Lxttchiorbiikl + LxttFintiikl !TODO: Re-Include
  ! Lyttchiorbiikl = Lyttchiorbiikl + LyttFintiikl !TODO: Re-Include
  ! Lzttchiorbiikl = Lzttchiorbiikl + LzttFintiikl !TODO: Re-Include

  if(rFreq(1) == 0) then
    call MPI_Reduce(MPI_IN_PLACE, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr)
    !call MPI_Reduce(MPI_IN_PLACE, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
  else
    call MPI_Reduce(tchiorbiikl, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr)
    !call MPI_Reduce(ttchiorbiikl, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(Lxttchiorbiikl, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(Lyttchiorbiikl, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
    !call MPI_Reduce(Lzttchiorbiikl, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, FreqComm(1), ierr) TODO: Re-Include
  end if

end subroutine eintshelinearsoc
