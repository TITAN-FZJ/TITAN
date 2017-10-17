module mod_self_consistency
  implicit none
  character(len=300)  :: default_file

contains

  ! This subroutine performs the self-consistency
  subroutine do_self_consistency()
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_parameters, only: outputunit_loop, lslatec, lnojac, mag_tol
    use mod_magnet, only: eps1, mx, my, mz, mabs, mtheta, mphi, mvec_cartesian, mvec_spherical, &
                          hw_count, iter
    use mod_mpi_pars, only: myrank_row_hw,mpitag
    use mod_system, only: s => sys
    use mod_dnsqe
    use mod_susceptibilities, only: lrot
    implicit none
    real(double),allocatable      :: fvec(:),jac(:,:),wa(:),sc_solu(:)
    real(double),allocatable      :: diag(:),qtf(:)
    real(double)                  :: epsfcn,factor
#if !defined(_OSX) && !defined(_JUQUEEN)
    real(double)                  :: ruser(1)
    integer                       :: iuser(1)
#else
    real(double),allocatable :: w(:,:)
#endif
    integer                       :: i,neq,maxfev,ml,mr,mode,nfev,njev,lwa,ifail=0

    neq = 4*s%nAtoms
    allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )

    ! Putting read eps1 existing solutions into esp1_solu (first guess of the subroutine)
    sc_solu(1:s%nAtoms)         = eps1
    sc_solu(s%nAtoms+1:2*s%nAtoms)   = mx
    sc_solu(2*s%nAtoms+1:3*s%nAtoms) = my
    sc_solu(3*s%nAtoms+1:4*s%nAtoms) = mz
    iter  = 1
    !mpitag = (Npl-Npl_i)*total_hw_npt1 + hw_count
    mpitag = hw_count
    if(myrank_row_hw==0) write(outputunit_loop,"('[self_consistency] Starting self-consistency:')")

#if defined(_OSX) || defined(_JUQUEEN)
    if(lslatec) then
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa),w(neq,4) )
      if(lnojac) then
        call dnsqe(sc_eqs_old,sc_jac_old,2,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs_old,sc_jac_old,1,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      end if
      ifail = ifail-1
    else
      lwa=neq*(neq+1)/2
      allocate( wa(lwa),w(neq,4) )
      if(lnojac) then
!         call c05nbf(sc_equations,neq,sc_solu,fvec,mag_tol,wa,lwa,ifail)
        maxfev = 200*(neq+1)
        ml = neq-1
        mr = neq-1
        epsfcn = 1.d-5
        mode = 1
        factor = 100.d0
        call c05ncf(sc_eqs_old,neq,sc_solu,fvec,mag_tol,maxfev,ml,mr,epsfcn,diag,mode,factor,0,nfev,jac,neq,wa,lwa,qtf,w,ifail)
      else
!         call c05pbf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,neq,mag_tol,wa,lwa,ifail)
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
!         diag(Npl+1:4*Npl) = 100.d0
        factor = 100.d0
        call c05pcf(sc_eqs_and_jac_old,neq,sc_solu,fvec,jac,neq,mag_tol,maxfev,diag,mode,factor,0,nfev,njev,wa,lwa,qtf,w,ifail)
      end if
    end if
    deallocate( w )
#else
    if(lslatec) then
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa) )
      if(lnojac) then
        call dnsqe(sc_eqs_old,sc_jac_old,2,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs_old,sc_jac_old,1,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      end if
      ifail = ifail-1
    else
      lwa=neq*(neq+1)/2
      allocate( wa(lwa) )
      if(lnojac) then
        maxfev = 200*(neq+1)
        ml = neq-1
        mr = neq-1
        epsfcn = 1.d-5
        mode = 1
        factor = 100.d0
        call c05qcf(sc_equations,neq,sc_solu,fvec,mag_tol,maxfev,ml,mr,epsfcn,mode,diag,factor,0,nfev,jac,wa,qtf,iuser,ruser,ifail)
      else
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
!         diag(npl+1:4*npl) = 100.d0
        factor = 100.d0
        call c05rcf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,mag_tol,maxfev,mode,diag,factor,0,nfev,njev,wa,qtf,iuser,ruser,ifail)
      end if
    end if
#endif

    deallocate(sc_solu,diag,qtf,fvec,jac,wa)

    ! Calculating new angles of GS magnetization in units of pi and magnetization vector
    do i = 1,s%nAtoms
      mabs(i)   = sqrt((mx(i)**2)+(my(i)**2)+(mz(i)**2))
      mtheta(i) = acos(mz(i)/mabs(i))/pi
      if(abs(mtheta(i))>1.d-8) then
        if(abs(abs(mtheta(i))-1.d0)>1.d-8) then
          mphi(i)   = atan2(my(i),mx(i))/pi
        else
          mphi(i) = 0.d0
        end if
        lrot = .true. ! Susceptibilities need to be rotated
      else
        mphi(i) = 0.d0
      end if
      mvec_cartesian(i,1) = mx(i)
      mvec_cartesian(i,2) = my(i)
      mvec_cartesian(i,3) = mz(i)
      mvec_spherical(i,1) = mabs(i)
      mvec_spherical(i,2) = mtheta(i)
      mvec_spherical(i,3) = mphi(i)
    end do

    return
  end subroutine do_self_consistency

  subroutine calcMagnetization(n, mx, my, mz, mp, size)
    !! Calculates occupation density and magnetization.
    use mod_f90_kind, only: double
    use mod_constants, only: cI, pi, cZero
    use mod_parameters, only: Ef
    use mod_SOC, only: llinearsoc, llineargfsoc
    use EnergyIntegration, only: y, wght
    use mod_system, only: s => sys
    !use mod_BrillouinZone, only: BZ
    use adaptiveMesh
    use TightBinding, only: nOrb,nOrb2
    use mod_mpi_pars
    !$  use omp_lib
    implicit none
    integer, intent(in) :: size
    real(double), dimension(s%nAtoms), intent(inout) :: n, mx, my, mz
    complex(double), dimension(s%nAtoms), intent(inout) :: mp

    integer  :: i,j
    real(double), dimension(3) :: kp
    real(double), dimension(s%nAtoms,nOrb) :: n_orb_u, n_orb_d
    real(double), dimension(s%nAtoms,nOrb) :: gdiaguur,gdiagddr
    complex(double), dimension(s%nAtoms,nOrb) :: gdiagud,gdiagdu
    complex(double), dimension(:,:,:,:), allocatable :: gf, gf_loc
    !--------------------- begin MPI vars --------------------
    integer :: ix
    integer :: ncount,ncount2
    integer :: mu,mup
    real(double) :: weight, ep
    ncount = s%nAtoms * 9
    ncount2 = size * size

    n_orb_u = 0.d0
    n_orb_d = 0.d0
    mp = cZero

    gdiaguur= 0.d0
    gdiagddr= 0.d0
    gdiagud = cZero
    gdiagdu = cZero

    !$omp parallel default(none) &
    !$omp& private(ix,ep,kp,weight,i,mu,mup,gf,gf_loc) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,wght,s,bzs,E_k_imag_mesh,Ef,y,gdiaguur,gdiagddr,gdiagud,gdiagdu)
    allocate(gf    (nOrb2,nOrb2,s%nAtoms,s%nAtoms))
    allocate(gf_loc(nOrb2,nOrb2,s%nAtoms,s%nAtoms))
    gf = cZero
    gf_loc = cZero

    if((llineargfsoc).or.(llinearsoc)) then
      !$omp do schedule(static)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         call greenlineargfsoc(Ef,ep,kp,gf_loc)
         gf = gf + gf_loc * weight
      end do
      !$omp end do nowait
    else
      !$omp do schedule(static)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         call green(Ef,ep,kp,gf_loc)
         gf = gf + gf_loc * weight
      end do
      !$omp end do nowait
    end if
    !$omp critical
    do i=1,s%nAtoms
      do mu=1,nOrb
        mup = mu+nOrb
        gdiaguur(i,mu) = gdiaguur(i,mu) + real(gf(mu,mu,i,i))
        gdiagddr(i,mu) = gdiagddr(i,mu) + real(gf(mup,mup,i,i))
        gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i)
        gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i)
      end do
    end do
    !$omp end critical

    deallocate(gf_loc, gf)
    !$omp end parallel

    do j=1,s%nAtoms
      mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
    end do

    call MPI_Allreduce(gdiaguur, n_orb_u, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
    call MPI_Allreduce(gdiagddr, n_orb_d, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp, s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)

    n_orb_u = 0.5d0 + n_orb_u/pi
    n_orb_d = 0.5d0 + n_orb_d/pi
    mp      = mp/pi
    mx      = real(mp)
    my      = aimag(mp)

    !n_orb_t = n_orb_u + n_orb_d
    !mag_orb = n_orb_u - n_orb_d

    do i = 1, s%nAtoms
      n(i) = sum(n_orb_u(i,:)) + sum(n_orb_d(i,:))
      mz(i) = sum(n_orb_u(i,5:9)) - sum(n_orb_d(i,5:9))
    end do
    return
  end subroutine calcMagnetization

  subroutine calcJacobian(jacobian, N)
    !! Calculated the Jacobian of the spin magnetization
    use mod_f90_kind, only: double
    use mod_constants, only: cI, pi, identorb18, cZero, pauli_dorb, cOne
    use mod_parameters, only: U, Ef, offset
    use mod_SOC, only: llinearsoc, llineargfsoc
    use EnergyIntegration, only: y, wght
    use mod_system, only: s => sys
    !use mod_BrillouinZone, only: BZ
    use adaptiveMesh
    use TightBinding, only: nOrb,nOrb2
    use mod_mpi_pars
    !$  use omp_lib
    implicit none
    integer, intent(in) :: N
    real(double), intent(inout), dimension(N,N) :: jacobian
    complex(double), dimension(nOrb2, nOrb2, s%nAtoms, s%nAtoms) :: gf,gvg

    integer :: i,j
    integer :: AllocateStatus
    integer :: i0,j0,sigma,sigmap

    real(double) :: kp(3), ep

    complex(double) :: mhalfU(4,s%nAtoms), weight
    complex(double), dimension(nOrb2, nOrb2, 4) :: pauli_components1,pauli_components2
    complex(double), dimension(nOrb2, nOrb2, 4) :: temp1,temp2
    complex(double), dimension(nOrb2, nOrb2) :: gij,gji,temp,paulitemp
    complex(double), dimension(:,:,:,:,:,:), allocatable :: gdHdxg,gvgdHdxgvg

    !--------------------- begin MPI vars --------------------
    integer :: ix, iz
    integer :: ncount,ncount2
    integer :: mu

    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount=s%nAtoms*9
    ncount2=N*N

    pauli_components1 = cZero
    pauli_components2 = cZero
!   Includes d orbitals in the charge component
    pauli_components1(:,:,1) = identorb18(:,:)
    pauli_components1(:,:,2) = pauli_dorb(1,:,:)
    pauli_components1(:,:,3) = pauli_dorb(2,:,:)
    pauli_components1(:,:,4) = pauli_dorb(3,:,:)
!   Excludes d orbitals in the charge component
    pauli_components2(5:9, 5:9, 1) = identorb18(5:9, 5:9)
    pauli_components2(14:18, 14:18, 1) = identorb18(14:18, 14:18)
    pauli_components2(:,:,2) = pauli_dorb(1,:,:)
    pauli_components2(:,:,3) = pauli_dorb(2,:,:)
    pauli_components2(:,:,4) = pauli_dorb(3,:,:)

    ! Prefactor -U/2 in dH/dm and 1 in dH/deps1
    do i=1,s%nAtoms
      mhalfU(1,i) = cOne
      mhalfU(2:4,i) = -0.5d0*U(i+offset)
    end do

    jacobian = 0.d0

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,j,i0,j0,mu,sigma,sigmap,ep,kp,weight,gf,gvg,gij,gji,temp,temp1,temp2,paulitemp,gdHdxg,gvgdHdxgvg) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,s,bzs,E_k_imag_mesh,Ef,y,wght,mhalfU,pauli_components1,pauli_components2,jacobian, myrank)
    allocate( gdHdxg(nOrb2,nOrb2,4,4,s%nAtoms,s%nAtoms), gvgdHdxgvg(nOrb2,nOrb2,4,4,s%nAtoms,s%nAtoms) , STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumk_jacobian] Not enough memory for: gdHdxg,gvgdHdxgvg")
    gdHdxg = cZero
    gvgdHdxgvg = cZero

   !$omp do schedule(static)
   do ix = 1, local_points
      ep = y(E_k_imag_mesh(1,ix))
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      weight = cmplx(1.d0,0.d0) * wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
        ! Green function on energy Ef + iy, and wave vector kp
        if((llineargfsoc).or.(llinearsoc)) then
          call greenlinearsoc(Ef,ep,kp,gf,gvg)
          gf = gf + gvg
        else
          call green(Ef,ep,kp,gf)
        end if

        do j=1,s%nAtoms
          do i=1,s%nAtoms
            gij = gf(:,:,i,j)
            gji = gf(:,:,j,i)

            do sigma = 1,4
              ! temp1 =  pauli*g_ij
              paulitemp = pauli_components1(:,:, sigma)
              call zgemm('n','n',18,18,18,cOne,paulitemp,18,gij,18,cZero,temp,18)
              temp1(:,:, sigma) = temp
            end do

            do sigma = 1,4
              ! temp2 = (-U/2) * sigma* g_ji
              paulitemp = pauli_components2(:,:, sigma)
              call zgemm('n','n',18,18,18,mhalfU(sigma,j),paulitemp,18,gji,18,cZero,temp,18)
              temp2(:,:, sigma) = temp
            end do

            do sigmap = 1,4
              do sigma = 1,4
                ! gdHdxg = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji
                gij = temp1(:,:, sigma)
                gji = temp2(:,:, sigmap)
                call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)
                gdHdxg(:,:,sigma,sigmap,i,j) = gdHdxg(:,:,sigma,sigmap,i,j) + temp
              end do
            end do


            if((llineargfsoc).or.(llinearsoc)) then ! non-linear term
              gij = gvg(:,:,i,j)
              gji = gvg(:,:,j,i)

              do sigma = 1,4
                ! temp1 = wkbz* pauli*gvg_ij
                paulitemp = pauli_components1(:,:,sigma)
                call zgemm('n','n',18,18,18,cOne,paulitemp,18,gij,18,cZero,temp,18)
                temp1(:,:,sigma) = temp
              end do

              do sigmap = 1,4
                ! temp2 = (-U/2) * sigma* gvg_ji
                paulitemp = pauli_components2(:,:,sigmap)
                call zgemm('n','n',18,18,18,mhalfU(sigmap,j),paulitemp,18,gji,18,cZero,temp,18)
                temp2(:,:,sigmap) = temp
              end do
              do sigmap = 1,4
                do sigma = 1,4
                  ! gdHdxg = temp1*temp2 = wkbz* pauli*gvg_ij*(-U/2)*sigma* gvg_ji
                  gij = temp1(:,:,sigma)
                  gji = temp2(:,:,sigmap)
                  call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)
                  gvgdHdxgvg(:,:,sigma,sigmap,i,j) = gvgdHdxgvg(:,:,sigma,sigmap,i,j) + temp
                end do
              end do
            end if ! End linear part
          end do ! End nAtoms i loop
        end do ! End nAtoms j loop

    end do ! End pn1 loop
    !$omp end do

    ! removing non-linear SOC term
    if((llineargfsoc).or.(llinearsoc)) gdHdxg = gdHdxg - gvgdHdxgvg

    !$omp critical
    do mu=1, nOrb2
      do j=1,s%nAtoms
        do i=1,s%nAtoms
          do sigmap=1,4
            do sigma=1,4
              i0 = (sigma-1)*s%nAtoms + i
              j0 = (sigmap-1)*s%nAtoms + j
              ! Trace over orbitals and spins of the real part
              jacobian(i0,j0) = jacobian(i0,j0) + real(gdHdxg(mu,mu,sigma,sigmap,i,j))
            end do
          end do
        end do
      end do
    end do
    !$omp end critical

    deallocate(gdHdxg,gvgdHdxgvg)
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, jacobian, ncount2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
    jacobian = jacobian/pi
    do i = s%nAtoms+1, 4*s%nAtoms
      jacobian(i,i) = jacobian(i,i) - 1.d0
    end do
    return
  end subroutine calcJacobian

  subroutine calcLGS()
    !! Calculates the orbital angular momentum ground state
    use mod_f90_kind, only: double
    use mod_constants, only: cZero,pi
    use mod_System, only: s => sys
    !use mod_BrillouinZone, only: BZ
    use adaptiveMesh
    use TightBinding, only: nOrb,nOrb2
    use mod_parameters, only: outputunit_loop, Ef
    use mod_magnet
    use EnergyIntegration, only: y, wght
    use mod_mpi_pars
    !$  use omp_lib
    implicit none

    integer :: AllocateStatus
    integer :: ix
    integer :: i,mu,nu,mup,nup
    real(double) :: kp(3)
    complex(double), dimension(:,:,:,:), allocatable :: gf
    complex(double), dimension(:,:,:), allocatable :: gupgd
    complex(double), dimension(:,:,:), allocatable :: gupgdint
    real(double) :: weight, ep
    !--------------------- begin MPI vars --------------------
    integer :: ncount
    ncount=s%nAtoms*nOrb*nOrb

    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    allocate( lxm(s%nAtoms), &
              lym(s%nAtoms), &
              lzm(s%nAtoms), &
              lxpm(s%nAtoms), &
              lypm(s%nAtoms), &
              lzpm(s%nAtoms), stat = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[L_gs] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")

    allocate(gupgdint(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[L_gs] Not enough meory for: gupgdint")

    if(myrank_row_hw==0) write(outputunit_loop,"('[L_gs] Calculating Orbital Angular Momentum ground state... ')")

    ! Calculating the number of particles for each spin and orbital using a complex integral

    gupgdint  = cZero

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf,gupgd) &
    !$omp& shared(local_points,s,E_k_imag_mesh,bzs,Ef,y,wght,gupgdint)

    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
             gupgd(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    gupgd   = cZero

    !$omp do schedule(static)
    do ix = 1, local_points
        kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
        ep = y(E_k_imag_mesh(1,ix))
        weight = bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix)) * wght(E_k_imag_mesh(1,ix))
        !Green function on energy Ef + iy, and wave vector kp
        call green(Ef,ep,kp,gf)

        do i=1,s%nAtoms
          do mu=1,nOrb
            mup = mu+nOrb
            do nu=1,nOrb
              nup = nu+nOrb
              gupgd(mu,nu,i) = gupgd(mu,nu,i) + (gf(mu,nu,i,i) + gf(mup,nup,i,i)) * weight
            end do
          end do
        end do
      end do
    !$omp end do nowait

    gupgd = gupgd/pi
    !$omp critical
      gupgdint = gupgdint + gupgd
    !$omp end critical

    deallocate(gf, gupgd)
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, gupgdint, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)


    lxpm = 0.d0
    lypm = 0.d0
    lzpm = 0.d0
    lxm = 0.d0
    lym = 0.d0
    lzm = 0.d0
    do nu=5,9
      do mu=5,9
        do i=1,s%nAtoms
          lxpm(i) = lxpm(i) + real(lxp(mu,nu,i)*gupgdint(nu,mu,i))
          lypm(i) = lypm(i) + real(lyp(mu,nu,i)*gupgdint(nu,mu,i))
          lzpm(i) = lzpm(i) + real(lzp(mu,nu,i)*gupgdint(nu,mu,i))
          lxm(i)  = lxm(i)  + real(lx (mu,nu)*gupgdint(nu,mu,i))
          lym(i)  = lym(i)  + real(ly (mu,nu)*gupgdint(nu,mu,i))
          lzm(i)  = lzm(i)  + real(lz (mu,nu)*gupgdint(nu,mu,i))
        end do
      end do
    end do

    ! Calculating angles of GS OAM (in units of pi)
    do i = 1,s%nAtoms
      labs(i)   = sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
      ltheta(i) = acos(lzm(i)/sqrt(lxm(i)**2+lym(i)**2+lzm(i)**2))/pi
      lphi(i)   = atan2(lym(i),lxm(i))/pi
      lpabs(i)  = sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2))
      lptheta(i)= acos(lzpm(i)/sqrt(lxpm(i)**2+lypm(i)**2+lzpm(i)**2))/pi
      lpphi(i)  = atan2(lypm(i),lxpm(i))/pi
    end do

    deallocate(gupgdint)
    return
  end subroutine calcLGS

  ! Tries to read eps1 and m if available - includes hdel, hdelp and hdelm calculations
  subroutine read_previous_results(lsuccess)
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_magnet, only: mabs, mx, my, mz, mp, mtheta, mphi, &
                          mvec_cartesian, mvec_spherical, &
                          eps1, hw_count, hw_list, lfield, &
                          hdel, hdelm, hdelp
    use mod_parameters, only: skipsc, outputunit_loop, lselfcon, U,&
                              magaxis, magaxisvec, offset, layertype
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_mpi_pars, only: myrank_row_hw, myrank, abortProgram
    use mod_susceptibilities, only: lrot
    use mod_Umatrix
    implicit none
    integer             :: i,err
    logical,intent(out) :: lsuccess

    lsuccess = .false.
    call read_sc_results(err,lsuccess)

    if(lsuccess) then
      if(err==0) then ! Same parameters
        if(skipsc) then ! Skip option ON
          if(myrank_row_hw==0) write(outputunit_loop,"('[read_previous_results] Existing results for the same parameters were read. Skipping self-consistency...')")
          lselfcon = .false.
        else ! Skip option OFF
          if(myrank_row_hw==0) write(outputunit_loop,"('[read_previous_results] Existing results for the same parameters were read. Updating values...')")
          lselfcon = .true.
        end if
      else ! Other parameters
        if(myrank_row_hw==0) write(outputunit_loop,"('[read_previous_results] Existing results for other parameters were read. Updating values...')")
        lselfcon = .true.
      end if
      ! Calculating angles of GS magnetization in units of pi and magnetization vector
      do i = 1,s%nAtoms
        mabs(i)   = sqrt((mx(i)**2)+(my(i)**2)+(mz(i)**2))
        mtheta(i) = acos(mz(i)/mabs(i))/pi
        if(abs(mtheta(i))>1.d-8) then
          if(abs(abs(mtheta(i))-1.d0)>1.d-8) then
            mphi(i)   = atan2(my(i),mx(i))/pi
          else
            mphi(i) = 0.d0
          end if
          lrot = .true. ! Susceptibilities need to be rotated
        else
          mphi(i) = 0.d0
        end if
        mvec_cartesian(i,1) = mx(i)
        mvec_cartesian(i,2) = my(i)
        mvec_cartesian(i,3) = mz(i)
        mvec_spherical(i,1) = mabs(i)
        mvec_spherical(i,2) = mtheta(i)
        mvec_spherical(i,3) = mphi(i)
      end do
    else !  If file doesn't exist
      if(myrank_row_hw==0) then
        write(outputunit_loop,"('[read_previous_results] Self-consistency file does not exist:')")
        write(outputunit_loop,"('[read_previous_results] ',a)") trim(default_file)
      end if
      lselfcon = .true.
      ! Parameters: center of band, magnetization, exchange split
      eps1 = 0.d0
      if(magaxis == -1) then
        continue
      else if(magaxis == -2) then
        magaxisvec = magaxisvec(1) * s%a1 + magaxisvec(2) * s%a2 + magaxisvec(3) * s%a3
      else if(magaxis == -3) then
        magaxisvec = [cos(magaxisvec(2)*pi/180)*sin(magaxisvec(1)*pi/180), sin(magaxisvec(2)*pi/180)*sin(magaxisvec(1)*pi/180), cos(magaxisvec(1)*pi/180)]
      else if(magaxis == 0) then
        magaxisvec = [0.d0, 0.d0, sign(1.0d0, hw_list(hw_count,1))]
      else if(magaxis >=1 .and. magaxis <= s%nAtoms) then
        !magaxisvec(1:3) = c_nn(1:3, magaxis)
        stop "Not Implemented"
      else
        if(myrank.eq.0) call abortProgram("[read_previous_results] Unknown magnetization direction!")
      end if
      magaxisvec = magaxisvec / sqrt(dot_product(magaxisvec, magaxisvec))
      magaxisvec = magaxisvec * 0.5d0

      mx = magaxisvec(1)
      my = magaxisvec(2)
      mz = magaxisvec(3)

      do i=1,s%nAtoms
        if(layertype(i+offset)==2) then
          mx(i) = mx(i) * sign(4.d0,hw_list(hw_count,1))
          my(i) = my(i) * sign(4.d0,hw_list(hw_count,1))
          mz(i) = mz(i) * sign(4.d0,hw_list(hw_count,1))
        end if
      end do

      if(lfield .and. magaxis == 0) then
        mx = mz * sin(hw_list(hw_count,2)*pi) * cos(hw_list(hw_count,3)*pi)
        my = mz * sin(hw_list(hw_count,2)*pi) * sin(hw_list(hw_count,3)*pi)
        mz = mz * cos(hw_list(hw_count,2)*pi)
      end if

      mp = cmplx(mx,my,double)

      ! Variables used in the hamiltonian
      do i=1,s%nAtoms
        hdel(i)   = 0.5d0*U(i+offset)*mz(i)
        hdelp(i)  = 0.5d0*U(i+offset)*mp(i)
      end do
      hdelm = conjg(hdelp)
    end if

    call init_Umatrix(eps1,hdel,hdelm,hdelp,s%nAtoms,nOrb)

    return
  end subroutine read_previous_results

  subroutine rotate_magnetization_to_field()
  !! Rotate the magnetization to the direction of the field (useful for SOC=F)
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_magnet, only: hw_count, hw_list, hhwx, hhwy, hhwz, &
                          mx, my, mz, mabs, mp

    use mod_parameters, only: outputunit_loop
    use mod_System, only: s => sys
    use mod_mpi_pars, only: myrank_row_hw
    implicit none
    integer :: i,sign
    real(double) :: mdotb

    if(myrank_row_hw==0) write(outputunit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    do i=1,s%nAtoms
      mdotb   = hhwx(i)*mx(i)+hhwy(i)*my(i)+hhwz(i)*mz(i)
      sign    = dble(mdotb/abs(mdotb))
      mabs(i) = sqrt(abs(mp(i))**2+(mz(i)**2))
      mx(i)   = sign*mabs(i)*sin(hw_list(hw_count,2)*pi)*cos(hw_list(hw_count,3)*pi)
      my(i)   = sign*mabs(i)*sin(hw_list(hw_count,2)*pi)*sin(hw_list(hw_count,3)*pi)
      mz(i)   = sign*mabs(i)*cos(hw_list(hw_count,2)*pi)
      mp(i)   = cmplx(mx(i),my(i),double)
    end do

    ! Writing new eps1 and rotated mag to file (without self-consistency)
    if(myrank_row_hw==0) call write_sc_results()

    ! Writing self-consistency results on screen
    if(myrank_row_hw==0)  call print_sc_results()

    return
  end subroutine rotate_magnetization_to_field

  ! Writes the self-consistency results on the screen
  subroutine print_sc_results()
    use mod_parameters, only: outputunit_loop,lGSL
    use mod_system, only: s => sys
    !use mod_mpi_pars
    use mod_magnet, only: eps1, mx, my, mz, mp, mphi, mtheta, mabs, &
                          lxpm, lypm, lzpm, lpphi, lptheta, lxm, lym, lzm, lpabs, labs
    implicit none
    integer :: i

    write(outputunit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    write(outputunit_loop,"(11x,' *************** Center of d bands: ***************')")
    do i=1,s%nAtoms
      write(outputunit_loop,"(26x,'eps1(',i2.0,')=',f11.8)") i,eps1(i)
    end do
    write(outputunit_loop,"(11x,' *********** Magnetization components: **********')")
    do i=1,s%nAtoms
      write(outputunit_loop,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mx(i),i,my(i),i,mz(i)
      if(abs(mp(i))/=0) write(outputunit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") mtheta(i),mphi(i)
    end do
    if(lGSL) then
      write(outputunit_loop,"(11x,' *** Orbital components in local frame:  ***')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(4x,'Lxp(',i2.0,')=',f11.8,4x,'Lyp(',i2.0,')=',f11.8,4x,'Lzp(',i2.0,')=',f11.8)") i,lxpm(i),i,lypm(i),i,lzpm(i)
        if(sqrt(lxpm(i)**2+lypm(i)**2)/=0) write(outputunit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") lptheta(i),lpphi(i)
      end do
      write(outputunit_loop,"(11x,' *** Orbital components in global frame: ***')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(4x,'Lx (',i2.0,')=',f11.8,4x,'Ly (',i2.0,')=',f11.8,4x,'Lz (',i2.0,')=',f11.8)") i,lxm(i),i,lym(i),i,lzm(i)
      end do
      write(outputunit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(4x,'M (',i2.0,') =',f11.8,4x,'Lp (',i2.0,')=',f11.8,4x,'L (',i2.0,') =',f11.8)") i,mabs(i),i,lpabs(i),i,labs(i)
      end do
    else
      write(outputunit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(27x,'M (',i2.0,') =',f11.8)") i,mabs(i)
      end do
    end if
    write(outputunit_loop,"('|----------=============================================================----------|')")

    return
  end subroutine print_sc_results




  ! This subroutine reads previous band-shifting and magnetization results
  subroutine read_sc_results(err,lsuccess)
    use mod_f90_kind, only: double
    use mod_constants, only: cI
    use mod_parameters, only: offset, fieldpart, eta, U,Utype,scfile, outputunit_loop, strSites, dfttype
    use EnergyIntegration, only: parts
    use mod_magnet, only: eps1, hdel, hdelm, hdelp, mp, mz, hw_count, mx, my, mz
    use mod_SOC, only: SOCc, socpart
    use mod_mpi_pars
    use mod_system, only: s => sys
    use mod_BrillouinZone, only: BZ
    implicit none
    character(len=300)  :: file = ""
    integer,intent(out) :: err
    logical,intent(out) :: lsuccess
    integer             :: i
    real(double)        :: previous_results(s%nAtoms,4)

    if(trim(scfile) /= "") then
      open(unit=99,file=scfile,status="old",iostat=err)
      if(err/=0) then
        if(myrank_row_hw==0) write(outputunit_loop,"('*** WARNING: Self-consistency file given on input file does not exist! Using default... ***')")
        scfile = " "
      end if
      close(99)
    end if

    lsuccess = .false.
    !   Reading previous results (mx, my, mz and eps1) from files (if available)
    if(trim(scfile)=="") then ! If a filename is not given in inputcard (or don't exist), use the default one
      write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc,trim(strSites),dfttype,parts,Utype,trim(fieldpart),BZ%nkpt,eta,trim(socpart)
      open(unit=99,file=file,status="old",iostat=err)
      if((err==0).and.(myrank_row_hw==0)) then
        write(outputunit_loop,"('[read_sc_results] Self-consistency file already exists. Reading it now...')")
        write(outputunit_loop,"(a)") trim(file)
      else
        default_file = trim(file)
      end if
    else ! If filename in inputcard exists or 2nd+ angular iteration
      if(((hw_count)==1)) then !.and.(Npl==Npl_i)) then ! Filename in inputcard (1st iteration on loop)
        open(unit=99,file=scfile,status="old",iostat=err)
        if((err==0).and.(myrank_row_hw==0)) then
          write(outputunit_loop,"('[read_sc_results] Using filename given in input file for self-consistency:')")
          write(outputunit_loop,"(a)") trim(scfile)
        end if
      else ! 2nd+ iteration, cheking if default file exists
        write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc,trim(strSites),dfttype,parts,Utype,trim(fieldpart),BZ%nkpt,eta,trim(socpart)
        open(unit=99,file=file,status="old",iostat=err)
        if(err==0) then ! Reading file for the same parameters
          if(myrank_row_hw==0) then
            write(outputunit_loop,"('[read_sc_results] Self-consistency file already exists. Reading it now...')")
            write(outputunit_loop,"(a)") trim(file)
          end if
        else ! Reading results from previous iteration
          open(unit=99,file=scfile,status="old",iostat=err)
          if((err==0).and.(myrank_row_hw==0)) then
            write(outputunit_loop,"('[read_sc_results] Using results from previous iteration as input for self-consistency:')")
            write(outputunit_loop,"(a)") trim(scfile)
          end if
          lsuccess   = .true. ! something was read
        end if
      end if
    end if
    if(err==0) then
      if(myrank_row_hw==0) then
        do i=1,s%nAtoms
          read(99,fmt=*) previous_results(i,1)
          read(99,fmt=*) previous_results(i,2)
          read(99,fmt=*) previous_results(i,3)
          read(99,fmt=*) previous_results(i,4)
        end do
      end if
      call MPI_Bcast(previous_results,4*s%nAtoms,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
      eps1(:) = previous_results(:,1)
      mx  (:) = previous_results(:,2)
      my  (:) = previous_results(:,3)
      mz  (:) = previous_results(:,4)
      mp  = mx + cI*my
      do i=1,s%nAtoms
        hdel(i)   = 0.5d0*U(i+offset)*mz(i)
        hdelp(i)  = 0.5d0*U(i+offset)*mp(i)
      end do
      hdelm = conjg(hdelp)
      if(lsuccess) then
        err = 1   ! Read different parameters
      else
        lsuccess   = .true. ! Read same parameters (err=0)
      end if
    else
      ! If file does not exist, try to read for parts-1
      close(99)
      write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc,trim(strSites),dfttype,parts-1,Utype,trim(fieldpart),BZ%nkpt,eta,trim(socpart)
      open(unit=99,file=file,status="old",iostat=err)
      if(err==0) then
        if(myrank_row_hw==0) then
          write(outputunit_loop,"('[read_sc_results] Self-consistency file does not exist. Reading results for parts-1 now...')")
          write(outputunit_loop,"('[read_sc_results] Updating values obtained for parts-1...')")
          write(outputunit_loop,"(a)") file
        end if
        do i=1,s%nAtoms
          read(99,*) eps1(i)
          read(99,*) mx(i)
          read(99,*) my(i)
          read(99,*) mz(i)
        end do
        mp  = mx + cI*my
        do i=1,s%nAtoms
          hdel(i)   = 0.5d0*U(i+offset)*mz(i)
          hdelp(i)  = 0.5d0*U(i+offset)*mp(i)
        end do
        hdelm = conjg(hdelp)
        lsuccess = .true. ! Read...
        err = 1           ! ... different parameters
      end if
    end if
    close(99)
    return
  end subroutine read_sc_results

  subroutine write_sc_results()
    !! Writes the self-consistency results into files and broadcasts the scfile for the next iteration.
    use mod_parameters, only: fieldpart, eta, Utype,scfile, outputunit_loop, strSites, dfttype
    use EnergyIntegration, only: parts
    use mod_magnet, only: eps1, mx, my, mz
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use mod_BrillouinZone, only: BZ
    use mod_mpi_pars
    implicit none
    integer :: i

    if(myrank_row_hw == 0) then
      ! Writing new results (mx, my, mz and eps1) and mz to file
      write(outputunit_loop,"('[write_sc_results] Writing new eps1, mx, my and mz to file...')")
      write(scfile,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(strSites),dfttype,parts,Utype,trim(fieldpart),BZ%nkpt,eta,trim(socpart)
      open (unit=99,status='replace',file=scfile)
      do i=1,s%nAtoms
        write(99,"(es21.11,2x,'! eps1')") eps1(i)
        write(99,"(es21.11,2x,'! mx  ')") mx(i)
        write(99,"(es21.11,2x,'! my  ')") my(i)
        write(99,"(es21.11,2x,'! mz  ')") mz(i)
      end do
      close(99)
    end if

    call MPI_Bcast(scfile, len(scfile), MPI_CHARACTER, 0, MPI_Comm_Row_hw,ierr)

    return
  end subroutine write_sc_results

  ! This subroutine calculates the self-consistency equations
  !  n  - n0    = 0
  !  mx - mx_in = 0
  !  my - my_in = 0
  !  mz - mz_in = 0
  ! and the correspondent jacobian
#if !defined(_OSX) && !defined(_JUQUEEN)
  subroutine sc_equations_and_jacobian(N,x,fvec,selfconjac,iuser,ruser,iflag)
    use mod_f90_kind, only: double
    use mod_constants, only: cI
    use mod_parameters, only: offset, U, outputunit, outputunit_loop, lontheflysc
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx, my, mz
    use mod_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag
    integer, intent(inout) :: iuser(*)
    real(double), intent(inout) :: ruser(*)
    real(double),dimension(N) :: x,fvec
    real(double),dimension(N,N) :: selfconjac
    real(double),dimension(s%nAtoms) :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms) :: mp_in

    ! Values used in the hamiltonian
    eps1  = x(           1:  s%nAtoms)
    mx_in = x(  s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in + cI*my_in

    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    call update_Umatrix(eps1, hdel, hdelm, hdelp, s%nAtoms, nOrb)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    select case (iflag)
    case(1)
      call calcMagnetization(n_t, mx, my, mz, mp, N)
      do i = 1, s%nAtoms
        fvec(i) = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation
        fvec(i+1*s%nAtoms) = mx(i) - mx_in(i)
        fvec(i+2*s%nAtoms) = my(i) - my_in(i)
        fvec(i+3*s%nAtoms) = mz(i) - mz_in(i)
      end do
      if(myrank_row_hw==0) then
        do i=1,s%nAtoms
          if(abs(mp(i))>1.d-10) then
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          else
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          end if
        end do
      end if
      if(lontheflysc) call write_sc_results()
    case(2)
      call calcJacobian(selfconjac, N)
    case default
      write(outputunit,"('[sc_equations_and_jacobian] Problem in self-consistency! iflag = ',I0)") iflag
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select

    iter = iter + 1

    return
  end subroutine sc_equations_and_jacobian

  ! For a given value of center of band eps1 it calculates the
  ! occupation number and the magnetic moment
  subroutine sc_equations(N,x,fvec,iuser,ruser,iflag)
    use mod_constants, only: cI
    use mod_parameters, only: offset, U, outputunit_loop, lontheflysc
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx ,my ,mz
    use mod_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag
    integer     , intent(inout)         :: iuser(*)
    real(double), intent(inout)         :: ruser(*)
    real(double),dimension(N)           :: x,fvec
    real(double),dimension(s%nAtoms)         :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)      :: mp_in

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+cI*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    call update_Umatrix(eps1, hdel, hdelm, hdelp, s%nAtoms, nOrb)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    call calcMagnetization(n_t, mx, my, mz, mp, N)
    do i = 1, s%nAtoms
      fvec(i) = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation
      fvec(i+1*s%nAtoms) = mx(i) - mx_in(i)
      fvec(i+2*s%nAtoms) = my(i) - my_in(i)
      fvec(i+3*s%nAtoms) = mz(i) - mz_in(i)
    end do

    if(myrank_row_hw==0) then
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
        end if
      end do
    end if

    if(lontheflysc) call write_sc_results()

    iter = iter + 1

    return
  end subroutine sc_equations

#endif

  ! This subroutine calculates the self-consistency equations
  !  n  - n0    = 0
  !  mx - mx_in = 0
  !  my - my_in = 0
  !  mz - mz_in = 0
  ! and the correspondent jacobian
  subroutine sc_eqs_and_jac_old(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_constants, only: cI
    use mod_parameters, only: offset, U, outputunit_loop, outputunit, lontheflysc
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx, my, mz
    use mod_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag,ldfjac
    real(double),dimension(N)        :: x,fvec
    real(double),dimension(ldfjac,N) :: selfconjac
    real(double),dimension(s%nAtoms)      :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)   :: mp_in

  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+cI*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    call update_Umatrix(eps1, hdel, hdelm, hdelp, s%nAtoms, nOrb)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    flag: select case (iflag)
    case(1)
      ! Calculating the number of particles for each spin and orbital using a complex integral
      call calcMagnetization(n_t, mx, my, mz, mp, N)
      do i = 1, s%nAtoms
        fvec(i) = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation
        fvec(i+1*s%nAtoms) = mx(i) - mx_in(i)
        fvec(i+2*s%nAtoms) = my(i) - my_in(i)
        fvec(i+3*s%nAtoms) = mz(i) - mz_in(i)
      end do

      if(myrank_row_hw==0) then
        do i=1,s%nAtoms
          if(abs(mp(i))>1.d-10) then
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          else
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          end if
        end do
      end if
      if(lontheflysc) call write_sc_results()
    case(2)
      call calcJacobian(selfconjac, N)
    case default
      write(outputunit,"('[sc_eqs_and_jac_old] Problem in self-consistency! iflag = ',I0)") iflag
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select flag

    iter = iter + 1

    return
  end subroutine sc_eqs_and_jac_old

  ! For a given value of center of band eps1 it calculates the
  ! occupation number and the magnetic moment
  subroutine sc_eqs_old(N,x,fvec,iflag)
    use mod_f90_kind, only: double
    use mod_constants, only: cI
    use mod_parameters, only: offset, U, outputunit_loop, lontheflysc
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx, my, mz
    use mod_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag
    real(double),dimension(N)           :: x,fvec
    real(double),dimension(s%nAtoms)         :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)      :: mp_in

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+cI*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    call update_Umatrix(eps1, hdel, hdelm, hdelp, s%nAtoms, nOrb)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    call calcMagnetization(n_t, mx, my, mz, mp, N)
    do i = 1, s%nAtoms
      fvec(i) = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation
      fvec(i+1*s%nAtoms) = mx(i) - mx_in(i)
      fvec(i+2*s%nAtoms) = my(i) - my_in(i)
      fvec(i+3*s%nAtoms) = mz(i) - mz_in(i)
    end do

    if(myrank_row_hw==0) then
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
        end if
      end do
    end if

    if(lontheflysc) call write_sc_results()

    iter = iter + 1

    return
  end subroutine sc_eqs_old

  subroutine sc_jac_old(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_f90_kind, only: double
    use mod_constants, only: cI
    use mod_parameters, only: offset, U
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp
    use mod_Umatrix
    use mod_mpi_pars
    implicit none
    integer       :: N,ldfjac,i,iflag
    real(double)  :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(double),dimension(s%nAtoms)         :: mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)      :: mp_in
    !--------------------- begin MPI vars --------------------

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+cI*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    call update_Umatrix(eps1, hdel, hdelm, hdelp, s%nAtoms, nOrb)

    fvec=fvec

    call calcJacobian(selfconjac, N)

    iter = iter + 1

    return
  end subroutine sc_jac_old

end module mod_self_consistency
