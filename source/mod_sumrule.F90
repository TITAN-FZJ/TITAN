module mod_sumrule
  use mod_kind, only: dp
  implicit none

  complex(dp), dimension(:,:,:,:), allocatable :: Smunuiivec

contains

  ! This subroutine checks if the sum rules are being satisfied
  subroutine sumrule(chiorb_hf)
    use mod_kind,       only: dp
    use mod_System,     only: s => sys
    use mod_constants,  only: levi_civita,StoC,CtoS,cZero
    use mod_magnet,     only: mxd,myd,mzd
    use mod_parameters, only: nOrb,output,dimens,sigmaimunu2i
    use mod_mpi_pars,   only: abortProgram,rField
    implicit none
    complex(dp), dimension(dimens,dimens), intent(in)    :: chiorb_hf

    integer :: AllocateStatus
    integer :: i,j,r,t,p,q,mu,nu,gamma,xi,m,n,k
    real(dp)   , dimension(3,s%nAtoms)             :: mvec
    complex(dp), dimension(dimens)                 :: Beff
    complex(dp), dimension(:,:,:,:,:),allocatable  :: lhs,rhs
    complex(dp), dimension(:,:)      ,allocatable  :: chiorb_hf_cart

    if(rField == 0) &
    write(output%unit_loop,"('[sumrule] Checking if the sum rule is satisfied... ')", advance='no')

    allocate(chiorb_hf_cart(dimens,dimens),lhs(3,3,nOrb,nOrb,s%nAtoms),rhs(3,3,nOrb,nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[sumrule] Not enough memory for: chiorb_hf_cart,lhs,rhs")

    mvec(1,:) = mxd(:)
    mvec(2,:) = myd(:)
    mvec(3,:) = mzd(:)

    call calcSmunu()

    chiorb_hf_cart = cZero
    ! Transforming the susceptibility to cartesian components
    do i = 1, s%nAtoms
      do j = 1, s%nAtoms
        do r = 1, 4
          do t = 1, 4
            do mu = 1, nOrb
              do nu = 1, nOrb
                do gamma = 1, nOrb
                  do xi = 1, nOrb
                    do p = 1, 4
                      do q = 1, 4
                        if(abs(StoC(p,r)) < 1.e-15_dp .or. abs(CtoS(t,q))  < 1.e-15_dp) cycle
                        chiorb_hf_cart(sigmaimunu2i(p,i,mu,nu),sigmaimunu2i(q,j,gamma,xi)) = chiorb_hf_cart(sigmaimunu2i(p,i,mu,nu),sigmaimunu2i(q,j,gamma,xi)) &
                                                                          + StoC(p,r) * chiorb_hf(sigmaimunu2i(r,i,mu,nu),sigmaimunu2i(t,j,gamma,xi))* CtoS(t,q)
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    call Beffective(mvec,s%nAtoms,nOrb,Beff)

    lhs = cZero
    rhs = cZero
    do i = 1, s%nAtoms
      do mu = 1, nOrb
        do nu = 1, nOrb
          do m = 1, 3
            do n = 1, 3
              do k = 1, 3
                if(abs(levi_civita(m,n,k)) > 1.e-15_dp) &
                  lhs(m,n,mu,nu,i) = lhs(m,n,mu,nu,i) + levi_civita(m,n,k)*Smunuiivec(mu,nu,k,i)
              end do

              do j = 1, s%nAtoms
                do gamma = 1, nOrb
                  do xi = 1, nOrb
                    do q = 1,3
                      do p = 1,3
                        if(abs(levi_civita(n,q,p)) > 1.e-15_dp) &
                          rhs(m,n,mu,nu,i) = rhs(m,n,mu,nu,i) + levi_civita(n,q,p)*chiorb_hf_cart(sigmaimunu2i(m+1,i,mu,nu),sigmaimunu2i(p+1,j,gamma,xi))*Beff(sigmaimunu2i(q+1,j,gamma,xi))
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    rhs = 2._dp*rhs

    if(sum(abs(lhs - rhs)) < 1.e-5_dp) then
      if(rField == 0) write(output%unit_loop,"(' YES! ')")
    else
      if(rField == 0) write(output%unit_loop,"(' NO! Difference: ',es16.9)") sum(abs(lhs - rhs))
    end if

    ! deallocate(lhs,rhs)
    deallocate(Smunuiivec)

  end subroutine sumrule

  ! This subroutine calculates the effective field on the Hamiltonian
  ! H = H_0 + \sigma.B_eff
  subroutine Beffective(mvec,nAtoms,nOrb,Beff)
    use mod_kind,       only: dp
    use mod_constants,  only: cZero
    use mod_parameters, only: dimens, sigmaimunu2i, offset
    use mod_SOC,        only: SOC, socscale
    use mod_parameters, only: Um
    use mod_magnet,     only: lvec, lfield, hhw
    use mod_System,     only: s => sys
    implicit none
    integer :: i,mu,nu,sigma
    integer                         , intent(in)  :: nAtoms, nOrb
    real(dp)   , dimension(3,nAtoms), intent(in)  :: mvec
    complex(dp), dimension(dimens)  , intent(out) :: Beff

    Beff = cZero
    do i=1,nAtoms
      do sigma=1,3
        do mu=1, nOrb
          if(lfield) Beff(sigmaimunu2i(sigma+1,i,mu,mu)) = Beff(sigmaimunu2i(sigma+1,i,mu,mu)) + hhw(sigma,i)
          ! p block
          if(SOC) then
            if((mu>=2).and.(mu<=4)) then
              do nu=2,4
                Beff(sigmaimunu2i(sigma+1,i,mu,nu))  = Beff(sigmaimunu2i(sigma+1,i,mu,nu)) + 0.5_dp*socscale * s%Types(s%Basis(i)%Material)%LambdaP * lvec(mu,nu,sigma)
              end do
            end if
          end if
          ! d block
          if(mu>=5) then
            Beff(sigmaimunu2i(sigma+1,i,mu,mu)) = Beff(sigmaimunu2i(sigma+1,i,mu,mu)) - 0.5_dp*Um(i+offset) * mvec(sigma,i)
            if(SOC) then
              do nu=5,nOrb
                Beff(sigmaimunu2i(sigma+1,i,mu,nu))  = Beff(sigmaimunu2i(sigma+1,i,mu,nu)) + 0.5_dp*socscale * s%Types(s%Basis(i)%Material)%LambdaD * lvec(mu,nu,sigma)
              end do
            end if
          end if
        end do
      end do
    end do

  end subroutine Beffective

  subroutine calcSmunu()
    !! Calculates the expectation value <S^munu_ii> = <c^+_imu c_inu>
    use mod_kind,          only: dp,int64
    use mod_constants,     only: cZero,cI,tpi
    use mod_System,        only: s => sys
    use EnergyIntegration, only: y, wght
    use mod_parameters,    only: nOrb, nOrb2, eta
    use mod_hamiltonian,   only: hamilt_local
    use mod_greenfunction, only: calc_green
    use adaptiveMesh,      only: local_points,E_k_imag_mesh,activeComm,bzs
    use mod_mpi_pars,      only: MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,ierr,abortProgram
    implicit none
    integer         :: AllocateStatus
    integer(int64)  :: ix
    integer         :: i,mu,nu,mup,nup
    real(dp) :: kp(3)
    real(dp) :: weight, ep
    complex(dp), dimension(:,:,:,:), allocatable :: gf
    complex(dp), dimension(:,:,:),   allocatable :: imguu,imgdd,imgud,imgdu
    integer :: ncount

    external :: MPI_Allreduce
    
    ncount=s%nAtoms*nOrb*nOrb

    allocate(imguu(nOrb, nOrb,s%nAtoms),imgdd(nOrb, nOrb,s%nAtoms),imgud(nOrb, nOrb,s%nAtoms),imgdu(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[calcSmunu] Not enough memory for: imguu,imgdd,imgud,imgdu")
    allocate(Smunuiivec(nOrb,nOrb,3,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[calcSmunu] Not enough memory for: Smunuiivec")

    imguu = cZero
    imgdd = cZero
    imgud = cZero
    imgdu = cZero

    ! Build local hamiltonian
    call hamilt_local(s)

    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat = AllocateStatus)
    gf = cZero

    !$omp parallel do schedule(static) &
    !$omp& default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(calc_green,local_points,s,nOrb,nOrb2,E_k_imag_mesh,bzs,eta,y,wght) &
    !$omp& reduction(+:imguu,imgdd,imgud,imgdu)
    do ix = 1, local_points
      ep = y(E_k_imag_mesh(1,ix))
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
      !Green function on energy Ef + iy, and wave vector kp
      call calc_green(s%Ef,ep+eta,s,kp,gf)
      do i=1,s%nAtoms
        do mu=1,nOrb
          mup = mu+nOrb
          do nu=1,nOrb
            nup = nu+nOrb

            imguu(mu,nu,i) = imguu(mu,nu,i) + ( gf(nu ,mu ,i,i) + conjg(gf(mu ,nu ,i,i)) ) * weight
            imgdd(mu,nu,i) = imgdd(mu,nu,i) + ( gf(nup,mup,i,i) + conjg(gf(mup,nup,i,i)) ) * weight
            imgud(mu,nu,i) = imgud(mu,nu,i) + ( gf(nu ,mup,i,i) + conjg(gf(mup,nu ,i,i)) ) * weight
            imgdu(mu,nu,i) = imgdu(mu,nu,i) + ( gf(nup,mu ,i,i) + conjg(gf(mu ,nup,i,i)) ) * weight
          end do
        end do
      end do
    end do
    !$omp end parallel do
    imguu = imguu / tpi
    imgdd = imgdd / tpi
    imgud = imgud / tpi
    imgdu = imgdu / tpi

    call MPI_Allreduce(MPI_IN_PLACE, imguu, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgdd, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgud, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgdu, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    do mu=1,nOrb
      imguu(mu,mu,:) = 0.5_dp + imguu(mu,mu,:)
      imgdd(mu,mu,:) = 0.5_dp + imgdd(mu,mu,:)
    end do

    Smunuiivec(:,:,1,:) =  imgud(:,:,:) + imgdu(:,:,:)
    Smunuiivec(:,:,2,:) = (imgud(:,:,:) - imgdu(:,:,:))*cI
    Smunuiivec(:,:,3,:) =  imguu(:,:,:) - imgdd(:,:,:)

    deallocate(gf)
    deallocate(imguu,imgdd,imgud,imgdu)

  end subroutine calcSmunu

end module mod_sumrule



! To be used in sumrule for test

! complex(dp) :: test1,test2


! if(rField == 0) then
! do i = 1, s%nAtoms ; do k = 1, 3 ; do nu = 1, nOrb ; do mu = 1, nOrb
! if(abs(Smunuiivec(mu,nu,k,i)) > 1.e-15_dp) write(15,"(a,2x,4(i0,2x),2(es16.9,2x))")'Smunuiivec',mu,nu,k,i,real(Smunuiivec(mu,nu,k,i)),dimag(Smunuiivec(mu,nu,k,i))
! end do ; end do ; end do ; end do
! do i = 1, s%nAtoms ; do j = 1, s%nAtoms ; do p = 1,3 ; do m = 1,3 ; do mu = 1, nOrb ; do nu = 1, nOrb ; do gamma = 1, nOrb ; do xi = 1, nOrb
! if(abs(chiorb_hf_cart(sigmaimunu2i(m+1,i,mu,nu),sigmaimunu2i(p+1,j,gamma,xi))) > 1.e-15_dp) write(16,"(a,2x,8(i0,2x),2(es16.9,2x))")'chiorb_hf_cart',i,j,m+1,p+1,mu,nu,gamma,xi,real(chiorb_hf_cart(sigmaimunu2i(m+1,i,mu,nu),sigmaimunu2i(p+1,j,gamma,xi))),dimag(chiorb_hf_cart(sigmaimunu2i(m+1,i,mu,nu),sigmaimunu2i(p+1,j,gamma,xi)))
! end do ; end do ; end do ; end do ; end do ; end do ; end do ; end do
! do j = 1, s%nAtoms ; do q = 1, 3 ; do gamma = 1, nOrb ; do xi = 1, nOrb
! if(abs(Beff(sigmaimunu2i(q+1,j,gamma,xi))) > 1.e-15_dp) write(17,"(a,2x,4(i0,2x),2(es16.9,2x))") 'Beff',q,j,gamma,xi,real(Beff(sigmaimunu2i(q+1,j,gamma,xi))),dimag(Beff(sigmaimunu2i(q+1,j,gamma,xi)))
! end do ; end do ; end do ; end do
! do i = 1, s%nAtoms ; do m = 1, 3 ; do n = 1, 3 ; do nu = 1, nOrb ; do mu = 1, nOrb
! if(abs(lhs(m,n,mu,nu,i) - rhs(m,n,mu,nu,i)) > 1.e-15_dp) write(18,"(a,2x,5(i0,2x),2(es16.9,2x))") 'sum',m,n,mu,nu,i,real(lhs(m,n,mu,nu,i) - rhs(m,n,mu,nu,i)),dimag(lhs(m,n,mu,nu,i) - rhs(m,n,mu,nu,i))
! end do ; end do ; end do ; end do ; end do
! do i = 1, s%nAtoms ; do m = 1, 3 ; do n = 1, 3 ; do nu = 1, nOrb ; do mu = 1, nOrb
! if(abs(lhs(m,n,mu,nu,i)) > 1.e-15_dp) write(19,"(a,2x,5(i0,2x),2(es16.9,2x))") 'lhs',m,n,mu,nu,i,real(lhs(m,n,mu,nu,i)),dimag(lhs(m,n,mu,nu,i))
! end do ; end do ; end do ; end do ; end do

! do i = 1, s%nAtoms ; do m = 1, 3 ; do n = 1, 3 ; do nu = 1, nOrb ; do mu = 1, nOrb
! if(abs(rhs(m,n,mu,nu,i)) > 1.e-15_dp) write(20,"(a,2x,5(i0,2x),2(es16.9,2x))") 'rhs',m,n,mu,nu,i,real(rhs(m,n,mu,nu,i)),dimag(rhs(m,n,mu,nu,i))
! end do ; end do ; end do ; end do ; end do
! test1=cZero
! test2=cZero
! do mu=1,nOrb
! test1 = test1 + lhs(1,2,mu,mu,1)
! test2 = test2 + lhs(2,1,mu,mu,1)
! end do
! write(*,"(a,2x,3(es16.9,2x))") 'mz ',mvec_cartesian(3,1),real(test1),dimag(test1)
! write(*,"(a,2x,3(es16.9,2x))") 'mz ',mvec_cartesian(3,1),real(test2),dimag(test2)
! test1=cZero
! test2=cZero
! do mu=5,nOrb
! test1 = test1 + lhs(1,2,mu,mu,1)
! test2 = test2 + lhs(2,1,mu,mu,1)
! end do
! write(*,"(a,2x,3(es16.9,2x))") 'mzd',mvec(3,1),real(test1),dimag(test1)
! write(*,"(a,2x,3(es16.9,2x))") 'mzd',mvec(3,1),real(test2),dimag(test2)
! write(*,*) sum(abs(lhs(:,:,:,:,:) - rhs(:,:,:,:,:)))
! end if
