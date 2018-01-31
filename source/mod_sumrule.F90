module mod_sumrule
  use mod_f90_kind, only: double
  implicit none

  complex(double), dimension(:,:,:,:), allocatable :: Smunuiivec

contains

  ! This subroutine checks if the sum rules are being satisfied
  subroutine sumrule(chiorb_hf)
    use mod_f90_kind,   only: double
    use mod_System,     only: s => sys
    use mod_constants,  only: levi_civita, StoC, CtoS, cZero
    use TightBinding,   only: nOrb
    use mod_magnet!,     only: mxd,myd,mzd
    use mod_parameters, only: output, dim, sigmaimunu2i
    use mod_mpi_pars,   only: abortProgram,rField
    implicit none
    integer :: AllocateStatus
    integer :: i,j,r,t,p,q,mu,nu,gamma,xi,m,n,k
    real(double)   , dimension(3,s%nAtoms)           :: mvec
    complex(double), dimension(dim)                  :: Beff
    complex(double), dimension(dim,dim)              :: chiorb_hf_cart
    complex(double), dimension(dim,dim), intent(in)  :: chiorb_hf
    complex(double), dimension(:,:,:,:,:),allocatable  :: lhs,rhs

    if(rField == 0) write(output%unit_loop,"('[sumrule] Checking if the sum rule is satisfied... ')")

    allocate(lhs(3,3,nOrb,nOrb,s%nAtoms),rhs(3,3,nOrb,nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[sumrule] Not enough memory for: lhs,rhs")


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
                        if(StoC(p,r) == cZero .or. CtoS(t,q) == 0) cycle
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
                if(levi_civita(m,n,k) /= 0.d0) lhs(m,n,mu,nu,i) = lhs(m,n,mu,nu,i) + levi_civita(m,n,k)*Smunuiivec(mu,nu,k,i)
              end do

              do j = 1, s%nAtoms
                do gamma = 1, nOrb
                  do xi = 1, nOrb
                    do q = 1,3
                      do p = 1,3
                        if(levi_civita(n,q,p) == 0.d0) cycle
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
    rhs = 2.d0*rhs

    if(sum(abs(lhs - rhs)) < 1.d-10) then
      if(rField == 0) write(output%unit_loop,"('[sumrule] YES! ')")
    else
      if(rField == 0) write(output%unit_loop,"('[sumrule] NO! ')")
    end if

    return
  end subroutine sumrule

  ! This subroutine calculates the effective field on the Hamiltonian
  ! H = H_0 + \sigma.B_eff
  subroutine Beffective(mvec,nAtoms,nOrb,Beff)
    use mod_f90_kind,   only: double
    use mod_constants,  only: cZero
    use mod_parameters, only:  dim, sigmaimunu2i, offset
    use mod_SOC,        only: SOC
    use mod_parameters, only: U
    use mod_magnet,     only: l, lfield, hhw
    use mod_System,     only: s => sys
    implicit none
    integer                             , intent(in)  :: nAtoms, nOrb
    real(double)   , dimension(3,nAtoms), intent(in)  :: mvec
    complex(double), dimension(dim)     , intent(out) :: Beff
    integer         :: i,mu,nu,sigma

    Beff = cZero
    do i=1,nAtoms
      do sigma=1,3
        do mu=1, nOrb
          if(lfield) Beff(sigmaimunu2i(sigma+1,i,mu,mu)) = Beff(sigmaimunu2i(sigma+1,i,mu,mu)) + hhw(sigma,i)
          if(mu>=5)  Beff(sigmaimunu2i(sigma+1,i,mu,mu)) = Beff(sigmaimunu2i(sigma+1,i,mu,mu)) - 0.5d0*U(i+offset) * mvec(sigma,i)
          if(SOC) then
            do nu=1,nOrb
              Beff(sigmaimunu2i(sigma+1,i,mu,nu))  = Beff(sigmaimunu2i(sigma+1,i,mu,mu)) + 0.5d0*s%Types(s%Basis(i)%Material)%Lambda * l(mu,nu,sigma)
            end do
          end if
        end do
      end do
    end do

    return
  end subroutine Beffective

  subroutine calcSmunu()
    !! Calculates the expectation value <S^munu_ii> = <c^+_imu c_inu>
    use mod_f90_kind, only: double
    use mod_constants, only: cZero,tpi
    use mod_System, only: s => sys
    use adaptiveMesh
    use TightBinding, only: nOrb,nOrb2
    use EnergyIntegration, only: y, wght
    use mod_mpi_pars
    implicit none
    integer      :: AllocateStatus
    integer*8    :: ix
    integer      :: i,mu,nu,mup,nup
    real(double) :: kp(3)
    real(double) :: weight, ep
    complex(double), dimension(:,:,:,:), allocatable :: gf
    complex(double), dimension(:,:,:),   allocatable :: imguu,imgdd,imgud
    !--------------------- begin MPI vars --------------------
    integer :: ncount
    ncount=s%nAtoms*nOrb*nOrb
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    allocate(imguu(nOrb, nOrb,s%nAtoms),imgdd(nOrb, nOrb,s%nAtoms),imgud(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[calcSmunu] Not enough memory for: imguu,imgdd,imgud")
    allocate(Smunuiivec(nOrb,nOrb,3,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[calcSmunu] Not enough memory for: Smunuiivec")

    imguu = cZero
    imgdd = cZero
    imgud = cZero

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(local_points,s,E_k_imag_mesh,bzs,y,wght,imguu,imgdd,imgud)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat = AllocateStatus)
    gf = cZero

    !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd) reduction(+:imgud)
    do ix = 1, local_points
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      ep = y(E_k_imag_mesh(1,ix))
      weight = bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix)) * wght(E_k_imag_mesh(1,ix))
      !Green function on energy Ef + iy, and wave vector kp
      call green(s%Ef,ep,kp,gf)

      do i=1,s%nAtoms
        do nu=1,nOrb
          do mu=1,nOrb
            mup = mu+nOrb
            nup = nu+nOrb

            imguu(mu,nu,i) = imguu(mu,nu,i) + ( gf(mu ,nu ,i,i) + conjg(gf(nu ,mu ,i,i)) ) * weight
            imgdd(mu,nu,i) = imgdd(mu,nu,i) + ( gf(mup,nup,i,i) + conjg(gf(nup,mup,i,i)) ) * weight
            imgud(mu,nu,i) = imgud(mu,nu,i) + ( gf(mup,nu ,i,i) + conjg(gf(nu ,mup,i,i)) ) * weight
          end do
        end do
      end do
    end do
    !$omp end do

    deallocate(gf)
    !$omp end parallel

    imguu = imguu / tpi
    imgdd = imgdd / tpi
    imgud = imgud / tpi
    call MPI_Allreduce(MPI_IN_PLACE, imguu, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgdd, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgud, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)


    do mu=1,nOrb
      imguu(mu,mu,:) = 0.5d0 + imguu(mu,mu,:)
      imgdd(mu,mu,:) = 0.5d0 + imgdd(mu,mu,:)
    end do

    Smunuiivec(:,:,1,:) = real (imgud(:,:,:))
    Smunuiivec(:,:,2,:) = aimag(imgud(:,:,:))
    Smunuiivec(:,:,3,:) = imguu(:,:,:) - imgdd(:,:,:)

    deallocate(imguu,imgdd,imgud)

    return
  end subroutine calcSmunu

end module mod_sumrule