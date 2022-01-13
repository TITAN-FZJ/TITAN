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
    use mod_parameters, only: output,dimens,sigmaimunu2i
    use mod_mpi_pars,   only: abortProgram,rField
    implicit none
    complex(dp), dimension(dimens,dimens), intent(in)    :: chiorb_hf

    integer :: AllocateStatus
    integer :: i,j,r,t,p,q,mu,nu,gama,xi,m,n,k
    real(dp)   , dimension(3,s%nAtoms)             :: mvec
    complex(dp), dimension(dimens)                 :: Beff
    complex(dp), dimension(:,:,:,:,:),allocatable  :: lhs,rhs
    complex(dp), dimension(:,:)      ,allocatable  :: chiorb_hf_cart

    if(rField == 0) &
    write(output%unit_loop,"('[sumrule] Checking if the sum rule is satisfied... ')", advance='no')

    allocate(chiorb_hf_cart(dimens,dimens),lhs(3,3,s%nOrb,s%nOrb,s%nAtoms),rhs(3,3,s%nOrb,s%nOrb,s%nAtoms), stat = AllocateStatus)
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
            do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
              do nu = 1,s%Types(s%Basis(i)%Material)%nOrb
                do gama = 1,s%Types(s%Basis(j)%Material)%nOrb
                  do xi = 1,s%Types(s%Basis(j)%Material)%nOrb
                    do p = 1, 4
                      do q = 1, 4
                        if(abs(StoC(p,r)) < 1.e-15_dp .or. abs(CtoS(t,q))  < 1.e-15_dp) cycle
                        chiorb_hf_cart(sigmaimunu2i(p,i,mu,nu),sigmaimunu2i(q,j,gama,xi)) = chiorb_hf_cart(sigmaimunu2i(p,i,mu,nu),sigmaimunu2i(q,j,gama,xi)) &
                                                                          + StoC(p,r) * chiorb_hf(sigmaimunu2i(r,i,mu,nu),sigmaimunu2i(t,j,gama,xi))* CtoS(t,q)
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

    call Beffective(mvec,s,Beff)

    lhs = cZero
    rhs = cZero
    do i = 1,s%nAtoms
      do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
        do nu = 1,s%Types(s%Basis(i)%Material)%nOrb
          do m = 1,3
            do n = 1,3
              do k = 1,3
                if(abs(levi_civita(m,n,k)) > 1.e-15_dp) &
                  lhs(m,n,mu,nu,i) = lhs(m,n,mu,nu,i) + levi_civita(m,n,k)*Smunuiivec(mu,nu,k,i)
              end do

              do j = 1,s%nAtoms
                do gama = 1,s%Types(s%Basis(j)%Material)%nOrb
                  do xi = 1,s%Types(s%Basis(j)%Material)%nOrb
                    do q = 1,3
                      do p = 1,3
                        if(abs(levi_civita(n,q,p)) > 1.e-15_dp) &
                          rhs(m,n,mu,nu,i) = rhs(m,n,mu,nu,i) + levi_civita(n,q,p)*chiorb_hf_cart(sigmaimunu2i(m+1,i,mu,nu),sigmaimunu2i(p+1,j,gama,xi))*Beff(sigmaimunu2i(q+1,j,gama,xi))
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
  subroutine Beffective(mvec,s,Beff)
    use mod_kind,       only: dp
    use mod_constants,  only: cZero
    use mod_parameters, only: dimens,sigmaimunu2i
    use mod_SOC,        only: SOC
    use mod_magnet,     only: lfield, hhw
    use mod_System,     only: System_type
    implicit none
    type(System_type),                  intent(in)  :: s
    real(dp)   , dimension(3,s%nAtoms), intent(in)  :: mvec
    complex(dp), dimension(dimens)    , intent(out) :: Beff
    integer :: i,mu,nu,mup,mud,nup,nud,sigma

    Beff = cZero
    do i=1,s%nAtoms
      do sigma=1,3
        ! Exchange term
        do mud=1,s%Types(s%Basis(i)%Material)%ndOrb
          mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
          Beff(sigmaimunu2i(sigma+1,i,mu,mu)) = Beff(sigmaimunu2i(sigma+1,i,mu,mu)) - 0.5_dp*s%Basis(i)%Um * mvec(sigma,i)
        end do
        ! External field term
        if(lfield) then
          do mu=1,s%Types(s%Basis(i)%Material)%nOrb
            Beff(sigmaimunu2i(sigma+1,i,mu,mu)) = Beff(sigmaimunu2i(sigma+1,i,mu,mu)) + hhw(sigma,i)
          end do
        end if
        ! SOC term
        if(SOC) then
          ! p block
          do nup=1,s%Types(s%Basis(i)%Material)%npOrb
            nu = s%Types(s%Basis(i)%Material)%pOrbs(nup)
            do mup=1,s%Types(s%Basis(i)%Material)%npOrb
              mu = s%Types(s%Basis(i)%Material)%pOrbs(mup)
              Beff(sigmaimunu2i(sigma+1,i,mu,nu))  = Beff(sigmaimunu2i(sigma+1,i,mu,nu)) + 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * s%Types(s%Basis(i)%Material)%lvec(mu,nu,sigma)
            end do
          end do
          ! d block
          do nud=1,s%Types(s%Basis(i)%Material)%ndOrb
            nu = s%Types(s%Basis(i)%Material)%dOrbs(nud)
            do mud=1,s%Types(s%Basis(i)%Material)%ndOrb
              mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
              Beff(sigmaimunu2i(sigma+1,i,mu,nu))  = Beff(sigmaimunu2i(sigma+1,i,mu,nu)) + 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * s%Types(s%Basis(i)%Material)%lvec(mu,nu,sigma)
            end do
          end do
        end if
      end do
    end do

  end subroutine Beffective

  subroutine calcSmunu()
    !! Calculates the expectation value <S^munu_ii> = <c^+_imu c_inu>
    use mod_kind,          only: dp,int64
    use mod_constants,     only: cZero,cI,tpi
    use mod_System,        only: s => sys
    use EnergyIntegration, only: y, wght
    use mod_parameters,    only: eta
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
    complex(dp), dimension(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms) :: gf
    complex(dp), dimension(s%nOrb, s%nOrb,s%nAtoms) :: imguu,imgdd,imgud,imgdu
    integer :: ncount

    external :: MPI_Allreduce
    
    ncount=s%nAtoms*s%nOrb*s%nOrb

    allocate(Smunuiivec(s%nOrb,s%nOrb,3,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[calcSmunu] Not enough memory for: Smunuiivec")

    imguu = cZero
    imgdd = cZero
    imgud = cZero
    imgdu = cZero

    ! Build local hamiltonian
    call hamilt_local(s)

    !$omp parallel do schedule(static) &
    !$omp& default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(calc_green,local_points,s,E_k_imag_mesh,bzs,eta,y,wght) &
    !$omp& reduction(+:imguu,imgdd,imgud,imgdu)
    do ix = 1, local_points
      ep = y(E_k_imag_mesh(1,ix))
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
      !Green function on energy Ef + iy, and wave vector kp
      call calc_green(s%Ef,ep+eta,s,kp,gf)
      do i=1,s%nAtoms
        do mu=1,s%Types(s%Basis(i)%Material)%nOrb
          mup = mu+s%Types(s%Basis(i)%Material)%nOrb
          do nu=1,s%Types(s%Basis(i)%Material)%nOrb
            nup = nu+s%Types(s%Basis(i)%Material)%nOrb

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

    do i=1,s%nAtoms
      do mu=1,s%Types(s%Basis(i)%Material)%nOrb
        imguu(mu,mu,i) = 0.5_dp + imguu(mu,mu,i)
        imgdd(mu,mu,i) = 0.5_dp + imgdd(mu,mu,i)
      end do
    end do

    Smunuiivec(:,:,1,:) =  imgud(:,:,:) + imgdu(:,:,:)
    Smunuiivec(:,:,2,:) = (imgud(:,:,:) - imgdu(:,:,:))*cI
    Smunuiivec(:,:,3,:) =  imguu(:,:,:) - imgdd(:,:,:)

  end subroutine calcSmunu

end module mod_sumrule
