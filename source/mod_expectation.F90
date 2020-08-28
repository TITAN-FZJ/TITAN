!! This module contains the calculation of the ground state expectation values using the GF
!! and using the eigenstates of the hamiltonian
module mod_expectation
  implicit none

  procedure(expectation_values_sub), pointer :: expectation_values => expectation_values_greenfunction
  procedure(calcLGS_sub),            pointer :: calcLGS => calcLGS_greenfunction

  abstract interface
    subroutine expectation_values_sub(s,rho,mp,mx,my,mz,deltas)
      use mod_kind, only: dp
      use mod_parameters, only: nOrb
      use mod_System,     only: System_type
      implicit none
      type(System_type),                          intent(in)  :: s
      real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
      real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: deltas
      complex(dp), dimension(nOrb,s%nAtoms), intent(out) :: mp
    end subroutine
  end interface

  abstract interface
    subroutine calcLGS_sub()
      implicit none
    end subroutine
  end interface

contains
  subroutine expectation_values_greenfunction(s,rho,mp,mx,my,mz,deltas)
    !! Calculates ground state (occupation and magnetization) quantities using the Green functions
    use mod_kind, only: dp
    use mod_constants,     only: pi,cZero
    use mod_SOC,           only: llinearsoc,llineargfsoc
    use EnergyIntegration, only: y,wght
    use mod_system,        only: System_type
    use adaptiveMesh,      only: bzs,E_k_imag_mesh,activeComm,local_points
    use mod_parameters,    only: nOrb,nOrb2,eta
    use mod_hamiltonian,   only: hamilt_local
    use mod_superconductivity, only: lsupercond
    use mod_mpi_pars
    implicit none
    type(System_type),                          intent(in)  :: s
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: deltas
    complex(dp), dimension(nOrb,s%nAtoms), intent(out) :: mp

    integer  :: i,j, AllocateStatus
    real(dp),    dimension(3)                    :: kp
    complex(dp), dimension(:,:),     allocatable :: gdiagud,gdiagdu
    real(dp),    dimension(:,:),     allocatable :: imguu,imgdd
    complex(dp), dimension(:,:,:,:), allocatable :: gf
    !--------------------- begin MPI vars --------------------
    integer(int64) :: ix
    integer :: ncount
    integer :: mu,mup
    real(dp) :: weight, ep
    ncount = s%nAtoms * nOrb

    allocate(imguu(nOrb,s%nAtoms),imgdd(nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) &
      call abortProgram("[expectation_values_greenfunction] Not enough memory for: imguu,imgdd")

    allocate(gdiagud(s%nAtoms,nOrb), gdiagdu(s%nAtoms,nOrb), stat = AllocateStatus)
    if(AllocateStatus /= 0) &
      call abortProgram("[expectation_values_greenfunction] Not enough memory for: gdiagdu, gdiagud")

    if(lsupercond) &
      call abortProgram("[expectation_values_greenfunction] Calculation of superconducting parameter Delta is not yet implemented with Green Functions.")

    imguu   = 0._dp
    imgdd   = 0._dp
    gdiagud = cZero
    gdiagdu = cZero
    deltas = 0._dp

    ! Build local hamiltonian
    if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

    !$omp parallel default(none) &
    !$omp& private(ix,ep,kp,weight,i,mu,mup,gf,AllocateStatus) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,eta,wght,s,nOrb,nOrb2,bzs,E_k_imag_mesh,y,gdiagud,gdiagdu,imguu,imgdd)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat=AllocateStatus)
    if(AllocateStatus /= 0) &
    call AbortProgram("[expectation_values_greenfunction] Not enough memory for: gf")
    gf = cZero

    if(llineargfsoc .or. llinearsoc) then
      !$omp do schedule(dynamic) reduction(+:imguu) reduction(+:imgdd) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         call greenlineargfsoc(s%Ef,ep+eta,s,kp,gf)
         do i=1,s%nAtoms
           do mu=1,nOrb
             mup = mu+nOrb
             gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i) * weight
             gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i) * weight

             imguu(mu,i) = imguu(mu,i) + real(gf(mu ,mu ,i,i)) * weight
             imgdd(mu,i) = imgdd(mu,i) + real(gf(mup,mup,i,i)) * weight
           end do
         end do
      end do
      !$omp end do
    else
      !$omp do schedule(dynamic) reduction(+:imguu) reduction(+:imgdd) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         call green(0._dp,ep+eta,s,kp,gf)
         do i=1,s%nAtoms
           do mu=1,nOrb
             mup = mu+nOrb
             gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i) * weight
             gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i) * weight

             imguu(mu,i) = imguu(mu,i) + real(gf(mu ,mu ,i,i)) * weight
             imgdd(mu,i) = imgdd(mu,i) + real(gf(mup,mup,i,i)) * weight
           end do
         end do
      end do
      !$omp end do
    end if

    deallocate(gf)
    !$omp end parallel
    imguu = imguu / pi
    imgdd = imgdd / pi

    do j=1,s%nAtoms
      mp(:,j)= gdiagdu(j,:) + conjg(gdiagud(j,:))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, imguu, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgdd, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp,    ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    mp      = mp/pi
    mx      = real(mp)
    my      = aimag(mp)

    do i = 1, s%nAtoms
      do mu=1,nOrb
        imguu(mu,i) = 0.5_dp + imguu(mu,i)
        imgdd(mu,i) = 0.5_dp + imgdd(mu,i)
        rho(mu,i) = imguu(mu,i) + imgdd(mu,i)
        mz (mu,i) = imguu(mu,i) - imgdd(mu,i)
      end do
    end do

    deallocate(imguu,imgdd)
    deallocate(gdiagdu, gdiagud)
  end subroutine expectation_values_greenfunction


  !   Calculates ground state quantities from eigenstates
  subroutine expectation_values_eigenstates(s,rho,mp,mx,my,mz,deltas)
    use mod_kind, only: dp,int64
    use mod_BrillouinZone,     only: realBZ
    use mod_parameters,        only: nOrb,nOrb2,output
    use mod_system,            only: System_type
    use mod_tools,             only: itos
    use mod_superconductivity, only: lsuperCond, superCond
    use mod_constants,         only: cZero
    use mod_hamiltonian,       only: hamilt_local,hamiltk
    use mod_mpi_pars,          only: abortProgram,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    type(System_type),                     intent(in)  :: s
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: deltas
    complex(dp), dimension(nOrb,s%nAtoms), intent(out) :: mp

    integer(int64)                           :: iz
    integer                                  :: info, ncount, lwork, dimH
    real(dp),    dimension(nOrb,s%nAtoms)    :: expec_0, expec_z
    complex(dp), dimension(nOrb,s%nAtoms)    :: expec_p
    real(dp),    dimension(nOrb,s%nAtoms)    :: expec_singlet
    real(dp),    dimension(:),  allocatable  :: rwork(:), eval(:)
    complex(dp),                allocatable  :: work(:), hk(:,:)

    dimH  = (s%nAtoms)*nOrb2*superCond
    lwork = 21*dimH
    ncount = nOrb*s%nAtoms

    allocate( hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork) )

    call hamilt_local(s)

    !$omp parallel default(none) &
    !$omp& firstprivate(lwork) &
    !$omp& private(iz,hk,eval,work,rwork,info,expec_0, expec_p, expec_z,expec_singlet) &
    !$omp& shared(s,dimH,output,realBZ,rho,mp,mz,lsuperCond,deltas)
    rho = 0._dp
    mz  = 0._dp
    mp  = cZero
    deltas = 0._dp
    !!$acc kernels 
    !!$acc parallel loop ! firstprivate(lwork) private(iz,hk,eval,work,rwork,info,expec_0, expec_p, expec_z,expec_singlet) shared(s,dimH,output,realBZ,rho,mp,mz,lsuperCond,deltas) reduction(+:rho,mp,mz,deltas)
    !$omp do reduction(+:rho,mp,mz,deltas) schedule(dynamic)
    do iz = 1,realBZ%workload
      ! Calculating the hamiltonian for a given k-point
      call hamiltk(s,realBZ%kp(1:3,iz),hk)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call zheev('V','L', dimH,hk,dimH,eval,work,lwork,rwork,info)

      if(info/=0) &
        call abortProgram("[expectation_values_eigenstates] Problem with diagonalization. info = " // itos(info))

      ! Calculating expectation values for a given k-point
      call expec_val(s, dimH, hk, eval, expec_0, expec_p, expec_z, expec_singlet)

      rho = rho + expec_0*realBZ%w(iz)
      mp  = mp  + expec_p*realBZ%w(iz)
      mz  = mz  + expec_z*realBZ%w(iz)

      ! Superconducting order parameter
      deltas = deltas + expec_singlet*realBZ%w(iz)
    end do
    !$omp end do
    !$omp end parallel
    !!$acc end kernels
    !!$acc end parallel loop

    !Gather and sum all the results from the different processes using an Allreduce clause
    call MPI_Allreduce(MPI_IN_PLACE, rho   , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mz    , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp    , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, deltas, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    mx = real(mp)
    my = aimag(mp)

    deallocate(hk,rwork,eval,work)

  end subroutine expectation_values_eigenstates



  !   Calculates ground state quantities from eigenstates
  subroutine expectation_eigenstates_fullhk(s,rho,mp,mx,my,mz,deltas)
    use mod_kind, only: dp,int64
    use mod_BrillouinZone,     only: realBZ
    use mod_parameters,        only: nOrb,dimH,output
    use mod_system,            only: System_type
    use mod_tools,             only: itos
    use mod_superconductivity, only: lsuperCond, superCond
    use mod_constants,         only: cZero
    use mod_hamiltonian,       only: hamilt_local,h0,fullhk
    use mod_mpi_pars,          only: abortProgram,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    type(System_type),                     intent(in)  :: s
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: deltas
    complex(dp), dimension(nOrb,s%nAtoms), intent(out) :: mp

    integer(int64)                           :: iz
    integer                                  :: info, ncount, lwork, dimHsc
    real(dp),    dimension(nOrb,s%nAtoms)    :: expec_0, expec_z
    complex(dp), dimension(nOrb,s%nAtoms)    :: expec_p
    real(dp),    dimension(nOrb,s%nAtoms)    :: expec_singlet
    real(dp),    dimension(:),  allocatable  :: rwork(:), eval(:)
    complex(dp),                allocatable  :: work(:), hk(:,:)

    dimHsc  = dimH*superCond
    lwork = 21*dimHsc
    ncount = nOrb*s%nAtoms

    allocate( hk(dimHsc,dimHsc),rwork(3*dimHsc-2),eval(dimHsc),work(lwork) )

    call hamilt_local(s)

    !$omp parallel default(none) &
    !$omp& firstprivate(lwork) &
    !$omp& private(iz,hk,eval,work,rwork,info,expec_0, expec_p, expec_z,expec_singlet) &
    !$omp& shared(s,dimHsc,h0,fullhk,output,realBZ,rho,mp,mz,lsuperCond,deltas)
    rho = 0._dp
    mz  = 0._dp
    mp  = cZero
    deltas = 0._dp
    !!$acc kernels 
    !!$acc parallel loop ! firstprivate(lwork) private(iz,hk,eval,work,rwork,info,expec_0, expec_p, expec_z,expec_singlet) shared(s,dimH,output,realBZ,rho,mp,mz,lsuperCond,deltas) reduction(+:rho,mp,mz,deltas)
    !$omp do reduction(+:rho,mp,mz,deltas) schedule(dynamic)
    do iz = 1,realBZ%workload
      ! hamiltonian for a given k-point
      hk = h0 + fullhk(:,:,iz)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call zheev('V','L', dimHsc,hk,dimHsc,eval,work,lwork,rwork,info)

      if(info/=0) &
        call abortProgram("[expectation_values_eigenstates] Problem with diagonalization. info = " // itos(info))

      ! Calculating expectation values for a given k-point
      call expec_val(s, dimHsc, hk, eval, expec_0, expec_p, expec_z, expec_singlet)

      rho = rho + expec_0*realBZ%w(iz)
      mp  = mp  + expec_p*realBZ%w(iz)
      mz  = mz  + expec_z*realBZ%w(iz)

      ! Superconducting order parameter
      deltas = deltas + expec_singlet*realBZ%w(iz)
    end do
    !$omp end do
    !$omp end parallel
    !!$acc end kernels
    !!$acc end parallel loop

    !Gather and sum all the results from the different processes using an Allreduce clause
    call MPI_Allreduce(MPI_IN_PLACE, rho   , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mz    , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp    , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, deltas, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    mx = real(mp)
    my = aimag(mp)

    deallocate(hk,rwork,eval,work)

  end subroutine expectation_eigenstates_fullhk



  ! subroutine expectation value of the operators 1 (occupation), Sp and Sz:
  subroutine expec_val(s, dimens, hk, eval, expec_0, expec_p, expec_z, expec_singlet)
    use mod_kind, only: dp
    use mod_constants,         only: cZero,pi,pauli_mat
    use mod_parameters,        only: nOrb, eta, isigmamu2n
    use mod_distributions,     only: fd_dist
    use mod_system,            only: System_type
    use mod_superconductivity, only: lsuperCond
    implicit none
    integer,                               intent(in)  :: dimens
    type(System_type),                     intent(in)  :: s
    real(dp),    dimension(dimens),        intent(in)  :: eval
    complex(dp), dimension(dimens,dimens), intent(in)  :: hk
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: expec_0, expec_z
    complex(dp), dimension(nOrb,s%nAtoms), intent(out) :: expec_p
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: expec_singlet

    real(dp)    :: fermi, beta
    integer     :: i, n, sigma, sigmap, mu
    real(dp)    :: f_n(dimens),f_n_negative(dimens),tanh_n(dimens)
    complex(dp) :: evec_isigmamu, evec_isigmamu_cong !dimens = 2*nOrb*nAtoms

    beta = 1._dp/(pi*eta)
    expec_0 = 0._dp
    expec_z = 0._dp
    expec_p = cZero
    expec_singlet = 0._dp

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)

    do concurrent (n = 1:dimens)
      f_n(n) = fd_dist(fermi, beta, eval(n))
    end do

    do concurrent (n = 1:dimens, lsuperCond)
      f_n_negative(n) = fd_dist(fermi, beta, -eval(n))
      tanh_n(n) = tanh(eval(n)*beta/2._dp)
    end do


    if(.not.lsupercond) then
      do concurrent(n = 1:dimens, i = 1:s%nAtoms, mu = 1:nOrb, sigma = 1:2)
        evec_isigmamu = hk(isigmamu2n(i,sigma,mu),n)
        evec_isigmamu_cong = conjg( evec_isigmamu )

        ! Charge
        expec_0(mu,i) = expec_0(mu,i) + f_n(n)*real( evec_isigmamu_cong*evec_isigmamu )

        !$OMP SIMD
        do sigmap = 1, 2
          evec_isigmamu = hk(isigmamu2n(i,sigmap,mu),n)
          ! M_p
          expec_p(mu,i) = expec_p(mu,i) + f_n(n)*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)*evec_isigmamu

          ! M_z
          expec_z(mu,i) = expec_z(mu,i) + f_n(n)*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3)*evec_isigmamu )
        end do
      end do
      ! !$OMP SIMD
      ! do i = 1, s%nAtoms
      !   !$OMP SIMD
      !   do mu = 1, nOrb
      !     !$OMP SIMD
      !     do sigma = 1, 2
      !       evec_isigmamu = hk(isigmamu2n(i,sigma,mu),n)
      !       evec_isigmamu_cong = conjg( evec_isigmamu )

      !       ! Charge
      !       expec_0(mu,i) = expec_0(mu,i) + f_n(n)*real( evec_isigmamu_cong*evec_isigmamu )

      !       !$OMP SIMD
      !       do sigmap = 1, 2
      !         evec_isigmamu = hk(isigmamu2n(i,sigmap,mu),n)
      !         ! M_p
      !         expec_p(mu,i) = expec_p(mu,i) + f_n(n)*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)*evec_isigmamu

      !         ! M_z
      !         expec_z(mu,i) = expec_z(mu,i) + f_n(n)*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3)*evec_isigmamu )
      !       end do ! sigmap
      !     end do ! sigma
      !   end do ! mu
      ! end do ! i
    else
      do concurrent(n = 1:dimens, i = 1:s%nAtoms, mu = 1:nOrb)
        ! up spin (using u's) + down spin (using v's)
        expec_0(mu,i) = expec_0(mu,i) + f_n(n)*real( conjg(hk(nOrb*2*(i-1)+mu,n))*hk(nOrb*2*(i-1)+mu,n) ) + f_n_negative(n)*real( conjg(hk(nOrb*s%nAtoms*2+mu+nOrb+(i-1)*nOrb*2,n))*hk(nOrb*s%nAtoms*2+mu+nOrb+(i-1)*nOrb*2,n) )

        expec_singlet(mu,i) = expec_singlet(mu,i) + 0.5_dp*s%Types(s%Basis(i)%Material)%lambda(mu)*tanh_n(n)*real( conjg(hk(isigmamu2n(i,1,mu)+nOrb*2*s%nAtoms,n))*hk(isigmamu2n(i,2,mu),n) )
      end do

      do concurrent(n = 1:dimens, i = 1:s%nAtoms, mu = 1:nOrb, sigma = 1:2, sigmap = 1:2)
        evec_isigmamu_cong = conjg( hk(isigmamu2n(i,sigma,mu),n) )
        evec_isigmamu = hk(isigmamu2n(i,sigmap,mu),n)
        ! M_p
        expec_p(mu,i) = expec_p(mu,i) + f_n(n)*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)*evec_isigmamu

        ! M_z
        expec_z(mu,i) = expec_z(mu,i) + f_n(n)*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3)*evec_isigmamu )
      end do
    end if

  end subroutine expec_val


  subroutine expec_val_n(s, dimens, evec, eval, expec_0, expec_p, expec_z, expec_singlet)
    !! Calculate the expectation value of the operators 1 (occupation), \sigma^+ and \sigma^z
    !! for a given state n (evec) with eigenenergy eval
    use mod_kind, only: dp
    use mod_constants,         only: cZero,pi,pauli_mat
    use mod_parameters,        only: nOrb, eta, isigmamu2n
    use mod_distributions,     only: fd_dist
    use mod_system,            only: System_type
    use mod_superconductivity, only: lsuperCond
    implicit none

    type(System_type),                     intent(in)  :: s
    integer,                               intent(in)  :: dimens
    complex(dp), dimension(dimens),        intent(in)  :: evec
    real(dp),                              intent(in)  :: eval
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: expec_0, expec_z
    complex(dp), dimension(nOrb,s%nAtoms), intent(out) :: expec_p
    real(dp),    dimension(nOrb,s%nAtoms), intent(out) :: expec_singlet

    integer     :: i, sigma, sigmap, mu
    real(dp)    :: f_n, f_n_negative, tanh_n
    real(dp)    :: fermi, beta
    complex(dp) :: evec_isigmamu, evec_isigmamu_cong

    beta = 1._dp/(pi*eta)
    expec_0 = 0._dp
    expec_z = 0._dp
    expec_p = cZero
    expec_singlet = 0._dp

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    ! Fermi-Dirac:
    f_n = fd_dist(fermi, beta, eval)

    if(.not.lsupercond) then
      do concurrent(i = 1:s%nAtoms, mu = 1:nOrb, sigma = 1:2)
        evec_isigmamu = evec(isigmamu2n(i,sigma,mu))
        evec_isigmamu_cong = conjg( evec_isigmamu )

        ! Charge
        expec_0(mu,i) = expec_0(mu,i) + f_n*real( evec_isigmamu_cong*evec_isigmamu )

        !$OMP SIMD
        do sigmap = 1, 2
          evec_isigmamu = evec(isigmamu2n(i,sigmap,mu))
          ! M_p
          expec_p(mu,i) = expec_p(mu,i) + f_n*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)*evec_isigmamu

          ! M_z
          expec_z(mu,i) = expec_z(mu,i) + f_n*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3)*evec_isigmamu )
        end do
      end do


      ! !$OMP SIMD
      ! do i = 1, s%nAtoms
      !   !$OMP SIMD
      !   do mu = 1, nOrb
      !     !$OMP SIMD
      !     do sigma = 1, 2
      !       evec_isigmamu = evec(isigmamu2n(i,sigma,mu))
      !       evec_isigmamu_cong = conjg( evec_isigmamu )

      !       ! Charge
      !       expec_0(mu,i) = expec_0(mu,i) + f_n*real( evec_isigmamu_cong*evec_isigmamu )

      !       !$OMP SIMD
      !       do sigmap = 1, 2
      !         evec_isigmamu = evec(isigmamu2n(i,sigmap,mu))
      !         ! M_p
      !         expec_p(mu,i) = expec_p(mu,i) + f_n*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)*evec_isigmamu

      !         ! M_z
      !         expec_z(mu,i) = expec_z(mu,i) + f_n*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3)*evec_isigmamu )
      !       end do
      !     end do
      !   end do
      ! end do
    else

      f_n_negative = fd_dist(fermi, beta, -eval)
      tanh_n = tanh(eval*beta/2._dp)
      do concurrent( i = 1:s%nAtoms, mu = 1:nOrb)
        ! up spin (using u's) + down spin (using v's)
        expec_0(mu,i) = expec_0(mu,i) + f_n*real( conjg(evec(nOrb*2*(i-1)+mu))*evec(nOrb*2*(i-1)+mu) ) + f_n_negative*real( conjg(evec(nOrb*s%nAtoms*2+mu+nOrb+(i-1)*nOrb*2))*evec(nOrb*s%nAtoms*2+mu+nOrb+(i-1)*nOrb*2) )

        expec_singlet(mu,i) = expec_singlet(mu,i) + 0.5_dp*s%Types(s%Basis(i)%Material)%lambda(mu)*tanh_n*real( conjg(evec(isigmamu2n(i,1,mu)+nOrb*2*s%nAtoms))*evec(isigmamu2n(i,2,mu)) )
      end do

      do concurrent( i = 1:s%nAtoms, mu = 1:nOrb, sigma = 1:2, sigmap = 1:2)
        evec_isigmamu_cong = conjg( evec(isigmamu2n(i,sigma,mu)) )
        evec_isigmamu = evec(isigmamu2n(i,sigmap,mu))
        ! M_p
        expec_p(mu,i) = expec_p(mu,i) + f_n*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)*evec_isigmamu

        ! M_z
        expec_z(mu,i) = expec_z(mu,i) + f_n*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3)*evec_isigmamu )
      end do
      ! !$OMP SIMD
      ! do i = 1, s%nAtoms
      !   !$OMP SIMD
      !   do mu = 1, nOrb
      !     ! up spin (using u's) + down spin (using v's)
      !     expec_0(mu,i) = expec_0(mu,i) + f_n*real( conjg(evec(nOrb*2*(i-1)+mu))*evec(nOrb*2*(i-1)+mu) ) + f_n_negative*real( conjg(evec(nOrb*s%nAtoms*2+mu+nOrb+(i-1)*nOrb*2))*evec(nOrb*s%nAtoms*2+mu+nOrb+(i-1)*nOrb*2) )

      !     expec_singlet(mu,i) = expec_singlet(mu,i) + 0.5_dp*s%Types(s%Basis(i)%Material)%lambda(mu)*tanh_n*real( conjg(evec(isigmamu2n(i,1,mu)+nOrb*2*s%nAtoms))*evec(isigmamu2n(i,2,mu)) )

      !     !$OMP SIMD
      !     do sigma = 1, 2
      !       evec_isigmamu_cong = conjg( evec(isigmamu2n(i,sigma,mu)) )

      !       !$OMP SIMD
      !       do sigmap = 1, 2
      !         evec_isigmamu = evec(isigmamu2n(i,sigmap,mu))
      !         ! M_p
      !         expec_p(mu,i) = expec_p(mu,i) + f_n*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)*evec_isigmamu

      !         ! M_z
      !         expec_z(mu,i) = expec_z(mu,i) + f_n*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3)*evec_isigmamu )
      !       end do
      !     end do
      !   end do
      ! end do

    end if

  end subroutine expec_val_n


  subroutine groundstate_L()
    !! Calculates the expectation value of the orbital angular momentum in the ground state
    use mod_kind, only: dp
    use mod_parameters,    only: output
    use mod_constants,     only: rad2deg
    use mod_System,        only: s => sys
    use mod_mpi_pars,      only: rField,abortProgram
    use mod_magnet,        only: labs,ltheta,lphi,lpabs,lptheta,lpphi,lxm,lym,lzm,lxpm,lypm,lzpm
    implicit none
    integer      :: i,AllocateStatus

    if(rField == 0) &
      write(output%unit_loop,"('[calcLGS] Calculating Orbital Angular Momentum ground state... ')")

    if(allocated(lxm)) deallocate(lxm)
    if(allocated(lym)) deallocate(lym)
    if(allocated(lzm)) deallocate(lzm)
    if(allocated(lxpm)) deallocate(lxpm)
    if(allocated(lypm)) deallocate(lypm)
    if(allocated(lzpm)) deallocate(lzpm)
    allocate( lxm(s%nAtoms), &
              lym(s%nAtoms), &
              lzm(s%nAtoms), &
              lxpm(s%nAtoms), &
              lypm(s%nAtoms), &
              lzpm(s%nAtoms), stat = AllocateStatus )
    if (AllocateStatus/=0) &
      call abortProgram("[calcLGS] Not enough memory for: lxm,lym,lzm,lxpm,lypm,lzpm")

    call calcLGS()

    ! Calculating angles of GS OAM (in degrees)
    do i = 1,s%nAtoms
      labs(i)   = sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
      if(labs(i)>1.e-8_dp) then
        ltheta(i) = acos(lzm(i)/labs(i))*rad2deg
      else
        ltheta(i) = 0._dp
      end if
      if(abs(ltheta(i))>1.e-8_dp) then
        if(abs(abs(ltheta(i))-180._dp)>1.e-8_dp) then
          lphi(i)   = atan2(lym(i),lxm(i))*rad2deg
        else
          lphi(i) = 0._dp
        end if
      else
        lphi(i) = 0._dp
      end if
      lpabs(i)  = sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2))
      if(abs(lpabs(i))>1.e-8_dp) then
        lptheta(i)= acos(lzpm(i)/lpabs(i))*rad2deg
      else
        lptheta(i) = 0._dp
      end if
      if(abs(lptheta(i))>1.e-8_dp) then
        if(abs(abs(lptheta(i))-180._dp)>1.e-8_dp) then
          lpphi(i)   = atan2(lypm(i),lxpm(i))*rad2deg
        else
          lpphi(i) = 0._dp
        end if
      else
        lpphi(i) = 0._dp
      end if
    end do

  end subroutine groundstate_L


  subroutine calcLGS_greenfunction()
    !! Calculates the expectation value of the orbital angular momentum in the ground state using green functions
    use mod_kind, only: dp,int64
    use mod_constants,     only: cZero,pi
    use mod_System,        only: s => sys
    use mod_parameters,    only: nOrb, nOrb2, eta
    use EnergyIntegration, only: y, wght
    use mod_magnet,        only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz
    use mod_hamiltonian,   only: hamilt_local
    use adaptiveMesh,      only: local_points,activeComm,E_k_imag_mesh,bzs
    use mod_mpi_pars,      only: abortProgram,MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,ierr
    implicit none
    integer(int64)    :: ix
    integer      :: AllocateStatus
    integer      :: i,mu,nu,mup,nup
    real(dp) :: kp(3)
    real(dp) :: weight, ep
    complex(dp), dimension(:,:,:,:), allocatable :: gf
    complex(dp), dimension(:,:,:),   allocatable :: gupgd
    !--------------------- begin MPI vars --------------------
    integer :: ncount
    ncount=s%nAtoms*nOrb*nOrb
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    allocate(gupgd(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) &
      call abortProgram("[calcLGS_greenfunction] Not enough memory for: gupgd")

    ! Build local hamiltonian
    call hamilt_local(s)

    ! Calculating the jacobian using a complex integral
    gupgd  = cZero
    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(local_points,s,nOrb,nOrb2,E_k_imag_mesh,bzs,eta,y,wght,gupgd)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat = AllocateStatus)
    if (AllocateStatus/=0) &
      call abortProgram("[calcLGS_greenfunction] Not enough memory for: gf")

    gf = cZero
    !$omp do schedule(dynamic) reduction(+:gupgd)
    do ix = 1, local_points
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      ep = y(E_k_imag_mesh(1,ix))
      weight = bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix)) * wght(E_k_imag_mesh(1,ix))
      !Green function on energy Ef + iy, and wave vector kp
      call green(s%Ef,ep+eta,s,kp,gf)

      site_i: do i=1,s%nAtoms
        orb_mu: do mu=1,nOrb
          mup = mu+nOrb
          orb_nu: do nu=1,nOrb
            nup = nu+nOrb
            gupgd(mu,nu,i) = gupgd(mu,nu,i) + (gf(mu,nu,i,i) + gf(mup,nup,i,i)) * weight
          end do orb_nu
        end do orb_mu
      end do site_i
    end do
    !$omp end do

    deallocate(gf)
    !$omp end parallel

    gupgd = gupgd / pi
    call MPI_Allreduce(MPI_IN_PLACE, gupgd, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    lxpm = 0._dp
    lypm = 0._dp
    lzpm = 0._dp
    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp

    do nu=5,9
      do mu=5,9
        do i=1,s%nAtoms
          lxpm(i) = lxpm(i) + real( lxp(mu,nu,i)*gupgd(nu,mu,i) )
          lypm(i) = lypm(i) + real( lyp(mu,nu,i)*gupgd(nu,mu,i) )
          lzpm(i) = lzpm(i) + real( lzp(mu,nu,i)*gupgd(nu,mu,i) )
          lxm (i) = lxm (i) + real( lx (mu,nu  )*gupgd(nu,mu,i) )
          lym (i) = lym (i) + real( ly (mu,nu  )*gupgd(nu,mu,i) )
          lzm (i) = lzm (i) + real( lz (mu,nu  )*gupgd(nu,mu,i) )
        end do
      end do
    end do

    deallocate(gupgd)
  end subroutine calcLGS_greenfunction


  !! Calculate the expectation value of the orbital momentum in the propagated states:
  subroutine expec_L_n(s, dimens, evec, eval, lxm, lym, lzm)
    use mod_kind, only: dp
    use mod_constants,     only: pi
    use mod_parameters,    only: nOrb,eta,isigmamu2n
    use mod_System,        only: System_type
    use mod_magnet,        only: lx,ly,lz
    use mod_distributions, only: fd_dist
    implicit none
    type(System_type),                intent(in)  :: s
    real(dp),                         intent(in)  :: eval
    integer,                          intent(in)  :: dimens
    real(dp),    dimension(s%nAtoms), intent(out) :: lxm, lym, lzm
    complex(dp), dimension(dimens),   intent(in)  :: evec
    integer                                       :: i, mu, nu, sigma
    real(dp)                                      :: f_n
    complex(dp)                                   :: prod, evec_isigmamu, evec_isigmamu_cong

    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp

    ! Fermi-Dirac:
    f_n = fd_dist(s%Ef, 1._dp/(pi*eta), eval)

    sites_loop: do i = 1, s%nAtoms
      do sigma = 1, 2
        do nu = 1, nOrb
          do mu = 1, nOrb
            evec_isigmamu = evec(isigmamu2n(i,sigma,mu))
            evec_isigmamu_cong = conjg( evec_isigmamu )
            prod = f_n*evec_isigmamu_cong*evec_isigmamu
            lxm (i) = lxm (i) + real( prod*lx (mu,nu) ) !> angular momentum at atomic site (i)
            lym (i) = lym (i) + real( prod*ly (mu,nu) )
            lzm (i) = lzm (i) + real( prod*lz (mu,nu) )
          end do
        end do
      end do
    end do sites_loop

  end subroutine expec_L_n


  !! Calculate the expectation value of the time-dependent Hamiltonian in the propagated states
  subroutine expec_H_n(s, kp, t, dimens, evec, eval, E_0)
    use mod_kind, only: dp
    use mod_constants,     only: pi
    use mod_parameters,    only: eta
    use mod_System,        only: System_type
    use mod_distributions, only: fd_dist
    use mod_imRK4,         only: build_td_hamiltonian
    implicit none

    type(System_type),           intent(in)  :: s
    integer,                     intent(in)  :: dimens
    complex(dp), dimension(dimens), intent(in)  :: evec
    real(dp),                    intent(in)  :: eval, t, kp(3)
    integer                                  :: i, j
    real(dp)                                 :: f_n, expec_H_0, E_0
    complex(dp)                              :: hamilt_t(dimens,dimens), hamilt_0(dimens,dimens)

    ! Fermi-Dirac:
    f_n = fd_dist(s%Ef, 1._dp/(pi*eta), eval)

    call build_td_hamiltonian(s,t,kp,eval,hamilt_t,hamilt_0)

    E_0 = 0._dp

    do i=1, dimens
      do j=1, dimens
        ! expec_H_0 = real( conjg( evec(i) ) * hamilt_t(i,j) * evec(j) )
        expec_H_0 = real( conjg( evec(i) ) * hamilt_0(i,j) * evec(j) )
        E_0       =  E_0 + f_n * expec_H_0
      end do
    end do

  end subroutine expec_H_n


  !   Calculates ground state quantities from eigenstates
  subroutine calcLGS_eigenstates()
    use mod_kind, only: dp,int64
    use mod_BrillouinZone,     only: realBZ
    use mod_constants,         only: pi,cZero
    use mod_parameters,        only: nOrb,dimH,output,eta,isigmamu2n
    use mod_System,            only: s => sys
    use mod_magnet,            only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz
    use mod_distributions,     only: fd_dist
    use mod_tools,             only: itos
    use mod_superconductivity, only: lsuperCond,superCond
    use mod_hamiltonian,       only: hamiltk,hamilt_local
    use mod_mpi_pars,          only: abortProgram,MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    integer(int64)                             :: iz
    integer                                    :: lwork,dimHsc,info,n,i,mu,nu,sigma
    real(dp)                                   :: fermi,beta
    real(dp),    dimension(:),     allocatable :: rwork,eval,f_n,f_n_negative
    complex(dp),                   allocatable :: work(:),hk(:,:),prod(:,:,:)

    dimHsc = dimH*superCond
    lwork  = 21*dimHsc
    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    beta = 1._dp/(pi*eta)

    allocate( hk(dimHsc,dimHsc),rwork(3*dimHsc-2),eval(dimHsc),f_n(dimHsc),f_n_negative(dimHsc),work(lwork),prod(nOrb,nOrb,s%nAtoms) )

    call hamilt_local(s)

    !$omp parallel default(none) &
    !$omp& firstprivate(lwork) &
    !$omp& private(iz,n,f_n,f_n_negative,i,sigma,mu,nu,hk,eval,work,rwork,info) &
    !$omp& shared(s,nOrb,dimH,dimHsc,output,realBZ,fermi,beta,eta,isigmamu2n,prod,lsupercond)

    prod = cZero
    !$omp do reduction(+:prod) schedule(dynamic)
    kloop: do iz = 1,realBZ%workload
      ! Calculating the hamiltonian for a given k-point
      call hamiltk(s,realBZ%kp(1:3,iz),hk)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call zheev('V','L',dimHsc,hk,dimHsc,eval,work,lwork,rwork,info)
      if(info/=0) &
        call abortProgram("[calcLGS_eigenstates] Problem with diagonalization. info = " // itos(info))

      do concurrent (n = 1:dimHsc)
        f_n(n) = fd_dist(fermi, beta, eval(n))
      end do

      if(.not.lsupercond) then

        do n = 1,dimHsc
          do i = 1,s%nAtoms
            do nu = 1,nOrb
              do mu = 1,nOrb
                do sigma=1,2
                  prod(mu,nu,i) = prod(mu,nu,i) + f_n(n)*conjg( hk(isigmamu2n(i,sigma,mu),n) )*hk(isigmamu2n(i,sigma,nu),n)*realBZ%w(iz)
                end do
              end do
            end do
          end do
        end do

      else
        do concurrent (n = 1:dimHsc)
          f_n_negative(n) = fd_dist(fermi, beta, -eval(n))
        end do

        do n = 1,dimHsc
          do i = 1,s%nAtoms
            do nu = 1,nOrb
              do mu = 1,nOrb
                prod(mu,nu,i) = prod(mu,nu,i) + ( f_n(n)*conjg( hk(isigmamu2n(i,1,mu),n) )*hk(isigmamu2n(i,1,nu),n)+f_n_negative(n)*conjg( hk(isigmamu2n(i,2,mu)+dimH,n))*hk(isigmamu2n(i,2,nu)+dimH,n) )*realBZ%w(iz)
              end do
            end do
          end do
        end do      
      end if
    end do kloop
    !$omp end do
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, prod, nOrb*nOrb*s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1) , ierr)

    ! Building different components of the orbital angular momentum
    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp
    lxpm = 0._dp
    lypm = 0._dp
    lzpm = 0._dp
    do i = 1, s%nAtoms
      do nu = 1, nOrb
        do mu = 1, nOrb
          lxm (i) = lxm (i) + real( prod(mu,nu,i)*lx (mu,nu  ) )
          lym (i) = lym (i) + real( prod(mu,nu,i)*ly (mu,nu  ) )
          lzm (i) = lzm (i) + real( prod(mu,nu,i)*lz (mu,nu  ) )
          lxpm(i) = lxpm(i) + real( prod(mu,nu,i)*lxp(mu,nu,i) )
          lypm(i) = lypm(i) + real( prod(mu,nu,i)*lyp(mu,nu,i) )
          lzpm(i) = lzpm(i) + real( prod(mu,nu,i)*lzp(mu,nu,i) )
        end do
      end do
    end do

    deallocate(hk,rwork,eval,work)

  end subroutine calcLGS_eigenstates

end module mod_expectation
