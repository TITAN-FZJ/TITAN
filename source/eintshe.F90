! ---------- Parallel spin current: Energy integration ---------
subroutine eintshe(e)
  use mod_f90_kind, only: double
  use mod_constants,     only: zero, zum, zi, tpi
  use mod_parameters, only: Ef, dim, sigmaimunu2i, eta, sigmai2i, offset
  use TightBinding, only: nOrb,nOrb2
  use mod_SOC, only: llineargfsoc
  use EnergyIntegration, only: y, wght, nepoints, x2, p2, generate_real_epoints, pn1
  use mod_system, only: s => sys !, n0sc1, n0sc2
  use mod_prefactors,    only: prefactor !, lxpt, lypt, lzpt, tlxp, tlyp, tlzp
  use mod_disturbances,     only: tchiorbiikl
  use mod_mpi_pars

  !use mod_currents,         only: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl !TODO: Re-Include
  !use mod_system,           only: n0sc1, n0sc2, n0sc
  !use mod_system,        only: r_nn, npln, nkpt, kbz, wkbz, n0sc1, n0sc2
  implicit none

  real(double),intent(in)   :: e

  integer :: AllocateStatus
  complex(double), dimension(:,:),allocatable         :: tFintiikl
  !complex(double), dimension(:,:,:), allocatable      :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl TODO:Re-Include

  integer :: i,j,l,mu,nu,gamma,xi,sigma,sigmap, neighbor
  real(double) :: kp(3)
  complex(double) :: wkbzc
  !complex(double) :: expikr(n0sc1:n0sc2) TODO:Re-Include
  complex(double),dimension(:,:,:,:),allocatable    :: gf,dtdk
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd
  complex(double),dimension(:,:),allocatable        :: df1iikl,pfdf1iikl
  !complex(double),dimension(n0sc1:n0sc2,s%nAtoms,nOrb, nOrb)    :: prett,preLxtt,preLytt,preLztt !TODO: Re-Include
  !complex(double),dimension(n0sc1:n0sc2,dimsigmaNpl,4), intent(out) :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl !TODO:Re-Include
  !--------------------- begin MPI vars --------------------
  integer :: ix, ix2, iz
  integer :: start, end, work, remainder, start1, end1, start2, end2
  integer :: ncountkl !,nncountkl !TODO: Re-Include
  ncountkl = dim*4
  !nncountkl = n0sc*dimsigmaNpl*4 !TODO: Re-Include
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
  tchiorbiikl     = zero
  ! ttchiorbiikl    = zero !TODO: Re-Include
  ! Lxttchiorbiikl  = zero !TODO: Re-Include
  ! Lyttchiorbiikl  = zero !TODO: Re-Include
  ! Lzttchiorbiikl  = zero !TODO: Re-Include


  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,ix2,iz,wkbzc,kp,i,j,l,mu,nu,gamma,xi,sigma,sigmap,neighbor,dtdk,gf,gfuu,gfud,gfdu,gfdd,df1iikl,pfdf1iikl,tFintiikl) & !,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,expikr,prett,preLxtt,preLytt,preLztt) &
  !$omp& shared(s,llineargfsoc,prefactor,e,y,x2,wght,p2,Ef,eta,start1,end1,start2,end2,sigmai2i,sigmaimunu2i,dim,offset, tchiorbiikl) !,n0sc1,n0sc2,lxpt,lypt,lzpt,tlxp,tlyp,tlzp)

  allocate(df1iikl(dim,4),pfdf1iikl(dim,4), &
           gf(nOrb2,nOrb2, s%nAtoms,s%nAtoms), &
           dtdk(nOrb,nOrb,s%nAtoms,s%nAtoms), &
           gfuu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfud(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
           gfdd(nOrb,nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus)
  if(AllocateStatus /= 0) call abortProgram("[sumk] Not enough memory for: df1iikl,pfdf1iikl,gf,dtdk,gfuu,gfud,gfdu,gfdd")

  allocate( tFintiikl(dim,4), STAT = AllocateStatus ) !, &
            !ttFintiikl   (n0sc1:n0sc2, dimsigmaNpl, 4), &
            !LxttFintiikl (n0sc1:n0sc2, dimsigmaNpl, 4), &
            !LyttFintiikl (n0sc1:n0sc2, dimsigmaNpl, 4), &
            !LzttFintiikl (n0sc1:n0sc2, dimsigmaNpl, 4), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[eintshe] Not enough memory for: tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl")

  tFintiikl    = zero
  !ttFintiikl   = zero !TODO: Re-Include
  !LxttFintiikl = zero !TODO: Re-Include
  !LyttFintiikl = zero !TODO: Re-Include
  !LzttFintiikl = zero !TODO: Re-Include

  !$omp do schedule(static) collapse(2)
  do ix = start1, end1 ! First and second integrations (in the complex plane)
    do iz=1, s%nkpt
      kp = s%kbz(:,iz)
      wkbzc = cmplx(s%wkbz(iz) * wght(ix),0.d0)
      df1iikl = zero

      ! Calculating derivative of in-plane and n.n. inter-plane hoppings
      call dtdksub(kp, dtdk)

      ! TODO: Fix currents
      !do neighbor = n0sc1,n0sc2
      !  expikr(neighbor) = exp(zi * dot_product(kp, r_nn(:,neighbor)))
      !end do

      ! Calculating the prefactor (L).t.exp - t.(L).exp
      ! TODO: Fix Currents
      ! do nu=1,9
      !   do mu=1,9
      !     do neighbor=n0sc1,n0sc2
      !       do i = 1, s%nAtoms
      !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i+offset)*expikr(neighbor)-(t0i(mu,nu,neighbor,i+offset)*conjg(expikr(neighbor)))
      !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !       end do
      !     end do
      !   end do
      ! end do

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

      do nu=1,nOrb
        do mu=1,nOrb
          do i=1, s%nAtoms
            do xi=1,nOrb
              do gamma=1,nOrb
                do l=1, s%nAtoms
                  do j=1, s%nAtoms
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
                  end do
                end do
              end do
            end do
          end do
        end do
      end do

      ! Multiplying the prefactor by the susceptibility and the k-point weight
      ! pfdf1iikl = prefactor*df1iikl
      call zgemm('n','n',dim,4,dim,wkbzc,prefactor,dim,df1iikl,dim,zero,pfdf1iikl,dim)

      tFintiikl = tFintiikl + pfdf1iikl

      ! TODO: Fix currents
      ! do sigmap = 1, 4
      !   do nu = 1, 9
      !     do mu = 1, 9
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

    end do ! iz
  end do
  !$omp end do nowait

  !$omp do schedule(static) collapse(2)
  do ix2 = start2, end2  ! Third integration (on the real axis)
    do iz=1, s%nkpt
      kp = s%kbz(:,iz)
      wkbzc = cmplx(s%wkbz(iz)*p2(ix2),0.d0)

      df1iikl = zero

      ! Calculating derivative of in-plane and n.n. inter-plane hoppings
      call dtdksub(kp, dtdk)

      ! TODO: Fix currents
      !do neighbor = n0sc1,n0sc2
      !  expikr(neighbor) = exp(zi * dot_product(kp, r_nn(:,neighbor)))
      !end do

      ! Calculating the prefactor (L).t.exp - t.(L).exp
      ! TODO: Fix Currents
      ! do nu=1,9
      !   do mu=1,9
      !     do neighbor=n0sc1,n0sc2
      !       do i = 1, s%nAtoms
      !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i+offset)*expikr(neighbor)-(t0i(mu,nu,neighbor,i+offset)*conjg(expikr(neighbor)))
      !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !       end do
      !     end do
      !   end do
      ! end do

        ! Green function at (k+q,E_F+E+iy)
        if(llineargfsoc) then
          call greenlineargfsoc(x2(ix2)+e,eta,kp,gf)
        else
          call green(x2(ix2)+e,eta,kp,gf)
        end if
        gfuu(:,:,:,:,1) = gf(     1:  nOrb,     1:  nOrb, :,:)
        gfud(:,:,:,:,1) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
        gfdu(:,:,:,:,1) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
        gfdd(:,:,:,:,1) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

        ! Green function at (k,E_F+iy)
        if(llineargfsoc) then
          call greenlineargfsoc(x2(ix2),eta,kp,gf)
        else
          call green(x2(ix2),eta,kp,gf)
        end if
        gfuu(:,:,:,:,2) = gf(     1:  nOrb,     1:  nOrb, :,:)
        gfud(:,:,:,:,2) = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
        gfdu(:,:,:,:,2) = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
        gfdd(:,:,:,:,2) = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

        do nu=1,nOrb
          do gamma=1,nOrb
            do mu=1,nOrb
              do i=1,s%nAtoms
                do xi=1,nOrb
                  do l=1,s%nAtoms
                    do j=1,s%nAtoms
                      df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

                      df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

                      df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

                      df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                      df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do


      ! Multiplying the prefactor by the susceptibility and the k-point weight
      ! pfdf1iikl = prefactor*df1iikl
      call zgemm('n','n',dim,4,dim,wkbzc,prefactor,dim,df1iikl,dim,zero,pfdf1iikl,dim)

      ! Integrating matrices
      tFintiikl = tFintiikl + pfdf1iikl

      ! TODO: Fix currents
      ! do sigmap=1,4
      !   do nu=1,9
      !     do mu=1,9
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
    end do ! iz
  end do
  !$omp end do nowait

  tFintiikl    = tFintiikl/tpi
  !ttFintiikl   = ttFintiikl/tpi
  !LxttFintiikl = LxttFintiikl/tpi
  !LyttFintiikl = LyttFintiikl/tpi
  !LzttFintiikl = LzttFintiikl/tpi

  !$omp critical
    tchiorbiikl = tchiorbiikl + tFintiikl
    !ttchiorbiikl = ttchiorbiikl + ttFintiikl
    !Lxttchiorbiikl = Lxttchiorbiikl + LxttFintiikl
    !Lyttchiorbiikl = Lyttchiorbiikl + LyttFintiikl
    !Lzttchiorbiikl = Lzttchiorbiikl + LzttFintiikl
  !$omp end critical

  deallocate(df1iikl,pfdf1iikl)
  deallocate(gf,dtdk,gfuu,gfud,gfdu,gfdd)
  deallocate(tFintiikl) !, ttFintiikl, LxttFintiikl, LyttFintiikl, LzttFintiikl)
  !$omp end parallel


  if(myrank_row == 0) then
    call MPI_Reduce(MPI_IN_PLACE, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
    !call MPI_Reduce(MPI_IN_PLACE, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
  else
    call MPI_Reduce(tchiorbiikl, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
    !call MPI_Reduce(ttchiorbiikl, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(Lxttchiorbiikl, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(Lyttchiorbiikl, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(Lzttchiorbiikl, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
  end if
  deallocate(tFintiikl) !,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl) TODO: Re-Incude

  return
end subroutine eintshe


subroutine eintshelinearsoc(e)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, zum, zi, tpi
  use mod_parameters, only: dim, Ef, eta, sigmai2i, sigmaimunu2i, offset
  use TightBinding, only: nOrb,nOrb2
  use mod_system, only: s => sys
  use EnergyIntegration, only: y, wght, nepoints, x2, p2, generate_real_epoints, pn1
  use mod_prefactors, only: prefactor, prefactorlsoc
  use mod_disturbances, only: tchiorbiikl
  use mod_mpi_pars


  !use mod_system, only: r_nn, npln, nkpt, kbz, wkbz, n0sc1, n0sc2
  !use mod_currents,         only: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl !TODO: Re-Include
  !use mod_system,           only: n0sc1, n0sc2, n0sc
  implicit none
  real(double),intent(in) :: e
  integer :: AllocateStatus
  complex(double), dimension(:,:),allocatable :: tFintiikl
  !complex(double), dimension(:,:,:), allocatable :: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl !TODO: Re-Include


  integer :: i,j,l,mu,nu,gamma,xi,sigma,sigmap,neighbor
  real(double) :: kp(3)
  !complex(double) :: expikr(n0sc1:n0sc2) !TODO: Re-Include
  complex(double) :: wkbzc
  complex(double),dimension(:,:,:,:),allocatable    :: gf,dtdk,gvg
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd,gvguu,gvgud,gvgdu,gvgdd
  complex(double),dimension(:,:),allocatable        :: df1iikl,pfdf1iikl,df1lsoc
  !complex(double),dimension(n0sc1:n0sc2,Npl,9,9)    :: prett,preLxtt,preLytt,preLztt !TODO: Re-Include

!--------------------- begin MPI vars --------------------
  integer :: ix,ix2, iz
  integer :: start, end, work, remainder, start1, end1, start2, end2
  integer :: ncountkl !,nncountkl !TODO: Re-Include
  ncountkl = dim*4
  !nncountkl = n0sc*dimsigmaNpl*4 !TODO: Re-Include
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
  tchiorbiikl     = zero
  ! ttchiorbiikl    = zero !TODO: Re-Include
  ! Lxttchiorbiikl  = zero !TODO: Re-Include
  ! Lyttchiorbiikl  = zero !TODO: Re-Include
  ! Lzttchiorbiikl  = zero !TODO: Re-Include

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,iz,wkbzc,kp,tFintiikl,df1iikl,pfdf1iikl,df1lsoc,dtdk,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,sigma,sigmap,i,j,l,mu,nu,gamma,xi,neighbor) &    !,expikr,prett,preLxtt,preLytt,preLztt
  !$omp& shared(start1,end1,start2,end2,tchiorbiikl,prefactor,prefactorlsoc,s,e,Ef,y,x2,wght,p2,eta,sigmai2i,sigmaimunu2i,dim,offset)                                                                                 !,n0sc1,n0sc2,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,lxpt,lypt,lzpt,tlxp,tlyp,tlzp

  allocate( df1iikl(dim,4),pfdf1iikl(dim,4),df1lsoc(dim,4), &
            dtdk(nOrb,nOrb,s%nAtoms,s%nAtoms), &
            gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
            gfuu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gfud(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gfdu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gfdd(nOrb,nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[sumklinearsoc] Not enough memory for: df1iikl,pfdf1iikl,df1lsoc,gf,dtdk,gfuu,gfud,gfdu,gfdd")

  allocate( gvg(nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
            gvguu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gvgud(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gvgdu(nOrb,nOrb,s%nAtoms,s%nAtoms,2), &
            gvgdd(nOrb,nOrb,s%nAtoms,s%nAtoms,2), STAT = AllocateStatus  )
  if (AllocateStatus/=0) call abortProgram("[sumklinearsoc] Not enough memory for: gvg,gvguu,gvgud,gvgdu,gvgdd")

  allocate( tFintiikl(dim,4), STAT = AllocateStatus ) !, & !TODO: Re-Include
            !ttFintiikl   (n0sc1:n0sc2, dimsigmaNpl, 4), &
            !LxttFintiikl (n0sc1:n0sc2, dimsigmaNpl, 4), &
            !LyttFintiikl (n0sc1:n0sc2, dimsigmaNpl, 4), &
            !LzttFintiikl (n0sc1:n0sc2, dimsigmaNpl, 4), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[eintshe] Not enough memory for: tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl")
  tFintiikl    = zero
  ! ttFintiikl   = zero !TODO: Re-Include
  ! LxttFintiikl = zero !TODO: Re-Include
  ! LyttFintiikl = zero !TODO: Re-Include
  ! LzttFintiikl = zero !TODO: Re-Include

  !$omp do schedule(auto)
  do ix = start1, end1 ! First and second integrations (in the complex plane)
    do iz = 1, s%nkpt
      kp = s%kbz(:,iz)
      wkbzc = cmplx(s%wkbz(iz)*wght(ix),0.d0)

      df1iikl = zero
      df1lsoc = zero

      ! Calculating derivative of in-plane and n.n. inter-plane hoppings
      call dtdksub(kp,dtdk)

      ! do neighbor=n0sc1,n0sc2 !TODO: Re-Include
      !   expikr(neighbor) = exp(zi*dot_product(kp, r_nn(:,neighbor)))
      ! end do

      ! Calculating the prefactor (L).t.exp - t.(L).exp
      ! do nu=1,9 !TODO: Re-Include
      !   do mu=1,9
      !     do i = 1, s%nAtoms
      !       do neighbor=n0sc1,n0sc2
      !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i+offset)*expikr(neighbor)-(t0i(mu,nu,neighbor,i+offset)*conjg(expikr(neighbor)))
      !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !       end do
      !     end do
      !   end do
      ! end do

        ! Green function at (k+q,E_F+E+iy)
        call greenlinearsoc(Ef+e,y(ix),kp,gf,gvg)
        gfuu(:,:,:,:,1)  = gf(     1:  nOrb,     1:  nOrb, :,:)
        gfud(:,:,:,:,1)  = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
        gfdu(:,:,:,:,1)  = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
        gfdd(:,:,:,:,1)  = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)
        gvguu(:,:,:,:,1) = gvg(     1:  nOrb,     1:  nOrb, :,:)
        gvgud(:,:,:,:,1) = gvg(     1:  nOrb,nOrb+1:nOrb2, :,:)
        gvgdu(:,:,:,:,1) = gvg(nOrb+1:nOrb2,     1:  nOrb, :,:)
        gvgdd(:,:,:,:,1) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

        ! Green function at (k,E_F+iy)
        call greenlinearsoc(Ef,y(ix),kp,gf,gvg)
        gfuu(:,:,:,:,2)  = gf(     1:  nOrb,     1:  nOrb, :,:)
        gfud(:,:,:,:,2)  = gf(     1:  nOrb,nOrb+1:nOrb2, :,:)
        gfdu(:,:,:,:,2)  = gf(nOrb+1:nOrb2,     1:  nOrb, :,:)
        gfdd(:,:,:,:,2)  = gf(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)
        gvguu(:,:,:,:,2) = gvg(     1:  nOrb,     1:  nOrb, :,:)
        gvgud(:,:,:,:,2) = gvg(     1:  nOrb,nOrb+1:nOrb2, :,:)
        gvgdu(:,:,:,:,2) = gvg(nOrb+1:nOrb2,     1:  nOrb, :,:)
        gvgdd(:,:,:,:,2) = gvg(nOrb+1:nOrb2,nOrb+1:nOrb2, :,:)

        do nu = 1, 9
          do mu = 1, 9
            do i = 1, s%nAtoms
              do xi = 1, 9
                do gamma = 1, 9
                  do l = 1, s%nAtoms
                    do j = 1, s%nAtoms
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

                    end do
                  end do
                end do
              end do
            end do
          end do
        end do


      ! Multiplying the prefactor by the susceptibility and the k-point weight
      ! pfdf1iikl = prefactorlsoc*df1iikl
      call zgemm('n','n',dim,4,dim,wkbzc,prefactorlsoc,dim,df1iikl,dim,zero,pfdf1iikl,dim)
      ! pfdf1lsoc = prefactorlsoc*df1iikl + prefactor*df1lsoc
      call zgemm('n','n',dim,4,dim,wkbzc,prefactor,dim,df1lsoc,dim,zum,pfdf1iikl,dim)


      ! Integrating matrices

      tFintiikl = tFintiikl + pfdf1iikl

      ! do sigmap=1,4 !TODO: Re-Include
      !   do i = 1, s%nAtoms
      !     do sigma=1,4
      !       do nu=1,9
      !         do mu=1,9
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

    end do ! iz
  end do
  !$omp end do nowait

  !$omp do schedule(static) collapse(2)
  do ix2 = start2, end2 ! Third integration (on the real axis)
    do iz=1,s%nkpt
      kp = s%kbz(:,iz)
      wkbzc = cmplx(s%wkbz(iz)*p2(ix2),0.d0)

      df1iikl = zero
      df1lsoc = zero

      ! Calculating derivative of in-plane and n.n. inter-plane hoppings
      call dtdksub(kp,dtdk)

      ! do neighbor=n0sc1,n0sc2 !TODO: Re-Include
      !   expikr(neighbor) = exp(zi*dot_product(kp, r_nn(:,neighbor)))
      ! end do

      ! Calculating the prefactor (L).t.exp - t.(L).exp
      ! do nu=1,9 !TODO: Re-Include
      !   do mu=1,9
      !     do neighbor=n0sc1,n0sc2
      !       do i = 1, s%nAtoms
      !         prett  (neighbor,i,mu,nu) = t0i(nu,mu,neighbor,i+offset)*expikr(neighbor)-(t0i(mu,nu,neighbor,i+offset)*conjg(expikr(neighbor)))
      !         preLxtt(neighbor,i,mu,nu) = lxpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlxp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLytt(neighbor,i,mu,nu) = lypt(i,neighbor,mu,nu)*expikr(neighbor)-(tlyp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !         preLztt(neighbor,i,mu,nu) = lzpt(i,neighbor,mu,nu)*expikr(neighbor)-(tlzp(i,neighbor,mu,nu)*conjg(expikr(neighbor)))
      !       end do
      !     end do
      !   end do
      ! end do

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

      do nu=1,9
        do mu=1,9
          do i = 1, s%nAtoms
            do gamma=1,9
              do xi = 1, 9
                do l = 1, s%nAtoms
                  do j = 1, s%nAtoms
                    !if(abs(j-l) > npln) cycle
                    df1iikl(sigmaimunu2i(1,i,mu,nu),1) = df1iikl(sigmaimunu2i(1,i,mu,nu),1)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(1,i,mu,nu),2) = df1iikl(sigmaimunu2i(1,i,mu,nu),2)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(1,i,mu,nu),3) = df1iikl(sigmaimunu2i(1,i,mu,nu),3)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(1,i,mu,nu),4) = df1iikl(sigmaimunu2i(1,i,mu,nu),4)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

                    df1iikl(sigmaimunu2i(2,i,mu,nu),1) = df1iikl(sigmaimunu2i(2,i,mu,nu),1)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(2,i,mu,nu),2) = df1iikl(sigmaimunu2i(2,i,mu,nu),2)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(2,i,mu,nu),3) = df1iikl(sigmaimunu2i(2,i,mu,nu),3)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(2,i,mu,nu),4) = df1iikl(sigmaimunu2i(2,i,mu,nu),4)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

                    df1iikl(sigmaimunu2i(3,i,mu,nu),1) = df1iikl(sigmaimunu2i(3,i,mu,nu),1)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(3,i,mu,nu),2) = df1iikl(sigmaimunu2i(3,i,mu,nu),2)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(3,i,mu,nu),3) = df1iikl(sigmaimunu2i(3,i,mu,nu),3)-zi*(gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(3,i,mu,nu),4) = df1iikl(sigmaimunu2i(3,i,mu,nu),4)-zi*(gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)

                    df1iikl(sigmaimunu2i(4,i,mu,nu),1) = df1iikl(sigmaimunu2i(4,i,mu,nu),1)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(4,i,mu,nu),2) = df1iikl(sigmaimunu2i(4,i,mu,nu),2)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(4,i,mu,nu),3) = df1iikl(sigmaimunu2i(4,i,mu,nu),3)-zi*(gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)
                    df1iikl(sigmaimunu2i(4,i,mu,nu),4) = df1iikl(sigmaimunu2i(4,i,mu,nu),4)-zi*(gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))*dtdk(gamma,xi,j,l)


                    df1lsoc(sigmaimunu2i(1,i,mu,nu),1) = df1lsoc(sigmaimunu2i(1,i,mu,nu),1)-zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(1,i,mu,nu),2) = df1lsoc(sigmaimunu2i(1,i,mu,nu),2)-zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(1,i,mu,nu),3) = df1lsoc(sigmaimunu2i(1,i,mu,nu),3)-zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(1,i,mu,nu),4) = df1lsoc(sigmaimunu2i(1,i,mu,nu),4)-zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)

                    df1lsoc(sigmaimunu2i(2,i,mu,nu),1) = df1lsoc(sigmaimunu2i(2,i,mu,nu),1)-zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(2,i,mu,nu),2) = df1lsoc(sigmaimunu2i(2,i,mu,nu),2)-zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfuu(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvguu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(2,i,mu,nu),3) = df1lsoc(sigmaimunu2i(2,i,mu,nu),3)-zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(2,i,mu,nu),4) = df1lsoc(sigmaimunu2i(2,i,mu,nu),4)-zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfud(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgud(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)

                    df1lsoc(sigmaimunu2i(3,i,mu,nu),1) = df1lsoc(sigmaimunu2i(3,i,mu,nu),1)-zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(3,i,mu,nu),2) = df1lsoc(sigmaimunu2i(3,i,mu,nu),2)-zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(3,i,mu,nu),3) = df1lsoc(sigmaimunu2i(3,i,mu,nu),3)-zi*((gvgdd(nu,gamma,i,j,1)-conjg(gvgdd(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfdd(nu,gamma,i,j,1)-conjg(gfdd(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(3,i,mu,nu),4) = df1lsoc(sigmaimunu2i(3,i,mu,nu),4)-zi*((gvgdu(nu,gamma,i,j,1)-conjg(gvgud(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfdu(nu,gamma,i,j,1)-conjg(gfud(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)

                    df1lsoc(sigmaimunu2i(4,i,mu,nu),1) = df1lsoc(sigmaimunu2i(4,i,mu,nu),1)-zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(4,i,mu,nu),2) = df1lsoc(sigmaimunu2i(4,i,mu,nu),2)-zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdu(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdu(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(4,i,mu,nu),3) = df1lsoc(sigmaimunu2i(4,i,mu,nu),3)-zi*((gvgud(nu,gamma,i,j,1)-conjg(gvgdu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfud(nu,gamma,i,j,1)-conjg(gfdu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                    df1lsoc(sigmaimunu2i(4,i,mu,nu),4) = df1lsoc(sigmaimunu2i(4,i,mu,nu),4)-zi*((gvguu(nu,gamma,i,j,1)-conjg(gvguu(gamma,nu,j,i,1)))*conjg(gfdd(mu,xi,i,l,2))  +  (gfuu(nu,gamma,i,j,1)-conjg(gfuu(gamma,nu,j,i,1)))*conjg(gvgdd(mu,xi,i,l,2)))*dtdk(gamma,xi,j,l)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do


      ! Multiplying the prefactor by the susceptibility and the k-point weight
      ! pfdf1iikl = prefactorlsoc*df1iikl
      call zgemm('n','n',dim,4,dim,wkbzc,prefactorlsoc,dim,df1iikl,dim,zero,pfdf1iikl,dim)
      ! pfdf1lsoc = prefactorlsoc*df1iikl + prefactor*df1lsoc
      call zgemm('n','n',dim,4,dim,wkbzc,prefactor,dim,df1lsoc,dim,zum,pfdf1iikl,dim)

      ! Integrating matrices
      tFintiikl = tFintiikl + pfdf1iikl

      ! do sigmap=1,4 !TODO: Re-Include
      !   do i = 1, s%nAtoms
      !     do sigma = 1, 4
      !       do nu=1,9
      !         do mu=1,9
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

    end do ! iz
  end do
  !$omp end do nowait

  tFintiikl    = tFintiikl/tpi
  ! ttFintiikl   = ttFintiikl/tpi !TODO: Re-Include
  ! LxttFintiikl = LxttFintiikl/tpi !TODO: Re-Include
  ! LyttFintiikl = LyttFintiikl/tpi !TODO: Re-Include
  ! LzttFintiikl = LzttFintiikl/tpi !TODO: Re-Include

  !omp critical
    tchiorbiikl = tchiorbiikl + tFintiikl
    ! ttchiorbiikl = ttchiorbiikl + ttFintiikl !TODO: Re-Include
    ! Lxttchiorbiikl = Lxttchiorbiikl + LxttFintiikl !TODO: Re-Include
    ! Lyttchiorbiikl = Lyttchiorbiikl + LyttFintiikl !TODO: Re-Include
    ! Lzttchiorbiikl = Lzttchiorbiikl + LzttFintiikl !TODO: Re-Include
  !omp end critical

  deallocate(df1iikl,pfdf1iikl,df1lsoc)
  deallocate(gvg,gvguu,gvgud,gvgdu,gvgdd)
  deallocate(gf,dtdk,gfuu,gfud,gfdu,gfdd)
!$omp end parallel


  if(myrank_row == 0) then
    call MPI_Reduce(MPI_IN_PLACE, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
    !call MPI_Reduce(MPI_IN_PLACE, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(MPI_IN_PLACE, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
  else
    call MPI_Reduce(tchiorbiikl, tchiorbiikl, ncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr)
    !call MPI_Reduce(ttchiorbiikl, ttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(Lxttchiorbiikl, Lxttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(Lyttchiorbiikl, Lyttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
    !call MPI_Reduce(Lzttchiorbiikl, Lzttchiorbiikl, nncountkl, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_Comm_Row, ierr) TODO: Re-Include
  end if
  deallocate(tFintiikl) !,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl) TODO: Re-Incude

  return
end subroutine eintshelinearsoc
