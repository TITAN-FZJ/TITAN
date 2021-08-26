module mod_expectation
!> This module contains the calculation of the ground state expectation values using the GF
!> and using the eigenstates of the hamiltonian
  implicit none

  procedure(expectation_values_sub), pointer :: expectation_values => expectation_values_greenfunction
  procedure(calc_GS_L_and_E_sub),    pointer :: calc_GS_L_and_E => calc_GS_L_and_E_greenfunction

  abstract interface
    subroutine expectation_values_sub(s,rho,mp,mx,my,mz,deltas)
      use mod_kind,       only: dp
      use mod_System,     only: System_type
      implicit none
      type(System_type),                       intent(in)  :: s
      real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
      real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: deltas
      complex(dp), dimension(s%nOrb,s%nAtoms), intent(out) :: mp
    end subroutine expectation_values_sub
  end interface

  abstract interface
    subroutine calc_GS_L_and_E_sub()
      implicit none
    end subroutine calc_GS_L_and_E_sub
  end interface

contains
  subroutine expectation_values_greenfunction(s,rho,mp,mx,my,mz,deltas)
  !> Calculates ground state (occupation and magnetization) quantities using the Green functions
    use mod_kind,              only: dp,int64
    use mod_constants,         only: pi,cZero
    use mod_SOC,               only: llinearsoc,llineargfsoc
    use EnergyIntegration,     only: y,wght
    use mod_system,            only: System_type
    use adaptiveMesh,          only: bzs,E_k_imag_mesh,activeComm,local_points
    use mod_parameters,        only: eta
    use mod_hamiltonian,       only: hamilt_local
    use mod_greenfunction,     only: calc_green
    use mod_superconductivity, only: lsuperCond
    use mod_mpi_pars,          only: abortProgram,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,ierr
    implicit none
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: rho,mx,my,mz
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: deltas
    complex(dp), dimension(s%nOrb,s%nAtoms), intent(out) :: mp

    integer  :: i,j
    real(dp),    dimension(3)                                 :: kp
    complex(dp), dimension(s%nAtoms,s%nOrb)                   :: gdiagud,gdiagdu
    real(dp),    dimension(s%nOrb,s%nAtoms)                   :: imguu,imgdd
    complex(dp), dimension(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms) :: gf

    !--------------------- begin MPI vars --------------------
    integer(int64) :: ix
    integer  :: ncount
    integer  :: mu,mup
    real(dp) :: weight,ep,fermi

    external :: MPI_Allreduce

    ncount = s%nAtoms * s%nOrb

    ! allocate(imguu(s%nOrb,s%nAtoms),imgdd(s%nOrb,s%nAtoms), stat = AllocateStatus)
    ! if(AllocateStatus/=0) &
    !   call abortProgram("[expectation_values_greenfunction] Not enough memory for: imguu,imgdd")
    !
    ! allocate(gdiagud(s%nAtoms,s%nOrb), gdiagdu(s%nAtoms,s%nOrb), stat = AllocateStatus)
    ! if(AllocateStatus /= 0) &
      ! call abortProgram("[expectation_values_greenfunction] Not enough memory for: gdiagdu, gdiagud")
    !
    ! allocate(gf(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms), stat=AllocateStatus)
    ! if(AllocateStatus /= 0) &
    !   call AbortProgram("[expectation_values_greenfunction] Not enough memory for: gf")

    gf = cZero

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)

    imguu   = 0._dp
    imgdd   = 0._dp
    gdiagud = cZero
    gdiagdu = cZero
    deltas  = 0._dp

    ! Build local hamiltonian
    if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

    !$omp parallel  &
    !$omp& default(none) &
    !$omp& firstprivate(gf) &
    !$omp& private(ix,ep,kp,weight,i,mu,mup) &
    !$omp& shared(calc_green,deltas,lsupercond,local_points,fermi,eta,wght,s,bzs,E_k_imag_mesh,y,imguu,imgdd,gdiagud,gdiagdu)
    gf = cZero

    !$omp do schedule(dynamic) reduction(+:imguu,imgdd,gdiagud,gdiagdu,deltas)
    do ix = 1, local_points
       ep = y(E_k_imag_mesh(1,ix))
       kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
       weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
       call calc_green(fermi,ep+eta,s,kp,gf)
       do i=1,s%nAtoms
         do mu=1,s%nOrb
           mup = mu+s%nOrb
           gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i) * weight
           gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i) * weight

           imguu(mu,i) = imguu(mu,i) + real(gf(mu ,mu ,i,i)) * weight
           imgdd(mu,i) = imgdd(mu,i) + real(gf(mup,mup,i,i)) * weight

           if ( lsuperCond ) then
             mup = mu + s%nOrb2 + s%norb
             deltas(mu,i) = deltas(mu,i) + real(gf(mu,mup,i,i)) * weight
           end if
         end do
       end do
    end do
    !$omp end do
    !$omp end parallel
    imguu = imguu / pi
    imgdd = imgdd / pi
    deltas = deltas / pi

    do j=1,s%nAtoms
      mp(:,j)= gdiagdu(j,:) + conjg(gdiagud(j,:))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, imguu, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgdd, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp   , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, activeComm, ierr)
    if ( lsuperCond ) &
      call MPI_Allreduce(MPI_IN_PLACE, deltas, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)

    mp      = mp/pi
    mx      = real(mp)
    my      = aimag(mp)

    do i=1,s%nAtoms
      do mu=1,s%nOrb
        imguu(mu,i) = 0.5_dp + imguu(mu,i)
        imgdd(mu,i) = 0.5_dp + imgdd(mu,i)
        rho(mu,i) = imguu(mu,i) + imgdd(mu,i)
        mz (mu,i) = imguu(mu,i) - imgdd(mu,i)
        if ( lsuperCond ) &
          deltas(mu,i) = s%Types(s%Basis(i)%Material)%lambda(mu)*deltas(mu,i)
      end do
    end do

  end subroutine expectation_values_greenfunction


  subroutine expectation_values_eigenstates(s,rho,mp,mx,my,mz,deltas)
  !>  Calculates ground state quantities from eigenstates
    use mod_kind,          only: dp,int64
    use mod_BrillouinZone, only: realBZ
    use mod_parameters,    only: output,dimHsc
    use mod_system,        only: System_type
    use mod_tools,         only: diagonalize,lwork
    use mod_constants,     only: cZero
    use mod_hamiltonian,   only: hamilt_local,h0,calchk
    use mod_mpi_pars,      only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: deltas
    complex(dp), dimension(s%nOrb,s%nAtoms), intent(out) :: mp

    integer(int64)                             :: iz
    integer                                    :: ncount,ilaenv
    real(dp),    dimension(s%nOrb,s%nAtoms)    :: expec_0,expec_z
    complex(dp), dimension(s%nOrb,s%nAtoms)    :: expec_p
    real(dp),    dimension(s%nOrb,s%nAtoms)    :: expec_d
    real(dp),    dimension(dimHsc)             :: eval
    complex(dp), dimension(dimHsc,dimHsc)      :: hk

    external :: MPI_Allreduce,ilaenv

    ncount = s%nOrb*s%nAtoms

    ! Getting lwork for diagonalization
    lwork = (ilaenv( 1, 'zhetrd', 'VU', dimHsc, -1, -1, -1 )+1)*dimHsc

    call hamilt_local(s)

    rho = 0._dp
    mz  = 0._dp
    mp  = cZero
    deltas = 0._dp

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp& private(iz,hk,eval,expec_0,expec_p,expec_z,expec_d) &
    !$omp& shared(s,h0,dimHsc,output,realBZ) &
    !$omp& reduction(+:rho,mp,mz,deltas)
    do iz = 1,realBZ%workload
      ! Calculating the hamiltonian for a given k-point
      hk = h0 + calchk(s,realBZ%kp(1:3,iz))

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call diagonalize(dimHsc,hk,eval)

      ! Calculating expectation values for a given k-point
      call expec_val(s,dimHsc,hk,eval,expec_0,expec_p,expec_z,expec_d)

      ! Occupation
      rho = rho + expec_0*realBZ%w(iz)
      ! Spin moments
      mp  = mp  + expec_p*realBZ%w(iz)
      mz  = mz  + expec_z*realBZ%w(iz)
      ! Superconducting order parameter
      deltas = deltas + expec_d*realBZ%w(iz)
    end do
    !$omp end parallel do

    !Gather and sum all the results from the different processes using an Allreduce clause
    call MPI_Allreduce(MPI_IN_PLACE, rho   , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mz    , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp    , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, deltas, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    mx = real(mp)
    my = aimag(mp)

  end subroutine expectation_values_eigenstates



  subroutine expectation_eigenstates_fullhk(s,rho,mp,mx,my,mz,deltas)
  !> Calculates ground state quantities from eigenstates using
  !> full hamiltonian matrix
    use mod_kind,              only: dp,int64
    use mod_BrillouinZone,     only: realBZ
    use mod_parameters,        only: dimHsc,output
    use mod_system,            only: System_type
    use mod_tools,             only: diagonalize,lwork
    use mod_constants,         only: cZero
    use mod_hamiltonian,       only: hamilt_local,h0,fullhk
    use mod_mpi_pars,          only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: deltas
    complex(dp), dimension(s%nOrb,s%nAtoms), intent(out) :: mp

    integer(int64)                           :: iz
    integer                                  :: ncount,ilaenv
    real(dp),    dimension(s%nOrb,s%nAtoms)  :: expec_0,expec_z
    complex(dp), dimension(s%nOrb,s%nAtoms)  :: expec_p
    real(dp),    dimension(s%nOrb,s%nAtoms)  :: expec_d
    real(dp),    dimension(dimHsc)           :: eval
    complex(dp), dimension(dimHsc,dimHsc)    :: hk

    external :: MPI_Allreduce,ilaenv

    ncount = s%nOrb*s%nAtoms

    ! Getting lwork for diagonalization
    lwork = (ilaenv( 1, 'zhetrd', 'VU', dimHsc, -1, -1, -1 )+1)*dimHsc

    call hamilt_local(s)

    rho = 0._dp
    mz  = 0._dp
    mp  = cZero
    deltas = 0._dp

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp& private(iz,hk,eval,expec_0,expec_p,expec_z,expec_d) &
    !$omp& shared(s,dimHsc,h0,fullhk,output,realBZ) &
    !$omp& reduction(+:rho,mp,mz,deltas)
    !!$acc kernels
    !$acc parallel loop private(iz,hk,eval,work,rwork,info,expec_0, expec_p, expec_z,expec_singlet) ! firstprivate(lwork) shared(s,dimHsc,output,realBZ,rho,mp,mz,lsuperCond,deltas) reduction(+:rho,mp,mz,deltas)
    do iz = 1,realBZ%workload
      ! hamiltonian for a given k-point
      hk = h0 + fullhk(:,:,iz)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call diagonalize(dimHsc,hk,eval)

      ! Calculating expectation values for a given k-point
      call expec_val(s,dimHsc,hk,eval,expec_0,expec_p,expec_z,expec_d)

      ! Occupation
      rho = rho + expec_0*realBZ%w(iz)
      ! Spin moments
      mp  = mp  + expec_p*realBZ%w(iz)
      mz  = mz  + expec_z*realBZ%w(iz)
      ! Superconducting order parameter
      deltas = deltas + expec_d*realBZ%w(iz)
    end do
    !$omp end parallel do

    !Gather and sum all the results from the different processes using an Allreduce clause
    call MPI_Allreduce(MPI_IN_PLACE, rho   , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mz    , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp    , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, deltas, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    mx = real(mp)
    my = aimag(mp)

  end subroutine expectation_eigenstates_fullhk


#ifdef _GPU
  subroutine expectation_eigenstates_fullhk_gpu(s,rho,mp,mx,my,mz,deltas)
  !>  Calculates ground state quantities - on the GPUs - from eigenstates
  !>  using full hamiltonian matrix
    use mod_kind,          only: dp,int64
    use mod_BrillouinZone, only: realBZ
    use mod_parameters,    only: dimHsc,output
    use mod_system,        only: System_type
    use mod_cuda,          only: h,diagonalize_gpu
    use nvtx,              only: nvtxStartRange,nvtxEndRange
    use mod_constants,     only: cZero
    use mod_hamiltonian,   only: hamilt_local_gpu,h0_d,fullhk_d
    use mod_mpi_pars,      only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: rho,mx,my,mz
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: deltas
    complex(dp), dimension(s%nOrb,s%nAtoms), intent(out) :: mp

    integer(int64)                                  :: iz
    integer                                         :: ncount,i,j,mu

    real(dp),    dimension(s%nOrb,s%nAtoms), device :: expec_0_d,rho_d
    real(dp),    dimension(s%nOrb,s%nAtoms), device :: expec_z_d,mz_d
    complex(dp), dimension(s%nOrb,s%nAtoms), device :: expec_p_d,mp_d
    real(dp),    dimension(s%nOrb,s%nAtoms), device :: expec_d_d,deltas_d
    real(dp),    dimension(s%nOrb,s%nAtoms), device :: mx_d,my_d

    complex(dp), dimension(:,:), allocatable,   device :: hk_d
    real(dp),    dimension(:),   allocatable,   device :: eval_d
    real(dp) :: weight_d

    external :: MPI_Allreduce

    ncount = s%nOrb*s%nAtoms

    allocate(hk_d(dimHsc,dimHsc),eval_d(dimHsc))

    call hamilt_local_gpu(s)

    rho_d = 0._dp
    mz_d  = 0._dp
    mp_d  = cZero
    deltas_d = 0._dp

    do iz = 1,realBZ%workload
      weight_d = realBZ%w(iz)

      !$cuf kernel do(2) <<< *, * >>>
      do j = 1, dimHsc
        do i = 1, dimHsc
          hk_d(i,j) = h0_d(i,j) + fullhk_d(i,j,iz)
        end do
      end do

      ! Starting marker of diagonalization for profiler
      call nvtxStartRange("Diagonalization",2)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call diagonalize_gpu(dimHsc,hk_d,eval_d,h)

      ! End of diagonalization marker
      call nvtxEndRange

      ! Starting marker of expectation value for profiler
      call nvtxStartRange("Expectation",3)

      ! Calculating expectation values for a given k-point
      call expec_val_gpu(s,dimHsc,hk_d,eval_d,expec_0_d,expec_p_d,expec_z_d,expec_d_d)

      ! End of expectation value marker
      call nvtxEndRange

      !$cuf kernel do(2) <<< (1,*), (9,*) >>>
      do j = 1,s%nAtoms
        do mu = 1,s%nOrb
          ! Occupation
          rho_d(mu,j) = rho_d(mu,j) + expec_0_d(mu,j)*weight_d
          ! Spin moments
          mp_d(mu,j)  = mp_d(mu,j)  + expec_p_d(mu,j)*weight_d
          mz_d(mu,j)  = mz_d(mu,j)  + expec_z_d(mu,j)*weight_d
          ! Superconducting order parameter
          deltas_d(mu,j) = deltas_d(mu,j) + expec_d_d(mu,j)*weight_d
        end do
      end do

    end do

    ! Starting marker of MPI reduction for profiler
    call nvtxStartRange("MPIReduction",4)

    !Gather and sum all the results from the different processes using an Allreduce clause
    call MPI_Allreduce(MPI_IN_PLACE, rho_d   , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mz_d    , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp_d    , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, deltas_d, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    ! End of MPI reduction marker
    call nvtxEndRange

    ! Starting marker of DevicetoHost copy for profiler
    call nvtxStartRange("DevicetoHost",5)

    rho = rho_d
    mp  = mp_d
    mz  = mz_d
    deltas = deltas_d
    mx = real(mp)
    my = aimag(mp)

    deallocate(hk_d,eval_d)

    ! End of DevicetoHost copy marker
    call nvtxEndRange

  end subroutine expectation_eigenstates_fullhk_gpu


  subroutine expec_val_gpu(s,dimens,hk_d,eval_d,expec_0_d,expec_p_d,expec_z_d,expec_d_d)
  !> Subroutine to calculate - on the GPUS - the expectation value of the operators
  !> 1 (occupation), Sp, Sz, and superconducting delta
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,pi,pauli_mat_d
    use mod_parameters,        only: eta,isigmamu2n_d
    use mod_system,            only: System_type
    use mod_superconductivity, only: lsuperCond
    use mod_distributions,     only: fd_dist_gpu
    implicit none
    integer,                                 intent(in)  :: dimens
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(dimens),          intent(in),  device :: eval_d
    complex(dp), dimension(dimens,dimens),   intent(in),  device :: hk_d
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out), device :: expec_0_d, expec_z_d
    complex(dp), dimension(s%nOrb,s%nAtoms), intent(out), device :: expec_p_d
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out), device :: expec_d_d

    integer     :: i,n,sigma,sigmap,mu,hdimens
    real(dp)    :: fermi,beta,x
    complex(dp) :: evec_isigmamu, evec_isigmamu_cong !dimens = 2*nOrb*nAtoms
    real(dp), device :: f_n_d(dimens),f_n_negative_d(dimens),tanh_n_d(dimens)

    real(dp), dimension(s%nOrb,s%nAtoms), device :: lambda_d
    complex(dp)  :: sum_c1
    real(dp)     :: sum_r1,sum_r2

    beta = 1._dp/(pi*eta)

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)

    !$acc parallel loop
    do n = 1,dimens
      f_n_d(n) = fd_dist_gpu(fermi, beta, eval_d(n))
    end do
    !$acc end parallel loop

    if(.not.lsupercond) then
      !$cuf kernel do(2) <<< (1,*), (9,*) >>>
      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          sum_r1 = 0._dp
          sum_r2 = 0._dp
          sum_c1 = cZero

          do n = 1,dimens
            do sigma = 1,2
              do sigmap = 1, 2

                evec_isigmamu = hk_d(isigmamu2n_d(i,sigma,mu),n)
                evec_isigmamu_cong = conjg( evec_isigmamu )

                ! Charge
                sum_r1 = sum_r1 + f_n_d(n)*real( evec_isigmamu_cong*evec_isigmamu )

                evec_isigmamu = evec_isigmamu_cong * hk_d(isigmamu2n_d(i,sigmap,mu),n)
                ! M_p
                sum_c1 = sum_c1 + f_n_d(n)*pauli_mat_d(sigma,sigmap,4)*evec_isigmamu

                ! M_z
                sum_r2 = sum_r2 + f_n_d(n)*real( pauli_mat_d(sigma,sigmap,3)*evec_isigmamu )
              end do
            end do
          end do
          expec_0_d(mu,i) = sum_r1 !expec_0_d(mu,i) + f_n_d(n)*real( conjg(hk_d(isigmamu2n_d(i,1,mu),n))*hk_d(isigmamu2n_d(i,1,mu),n) ) + f_n_negative_d(n)*real( conjg(hk_d(isigmamu2n_d(i,2,mu)+hdimens,n))*hk_d(isigmamu2n_d(i,2,mu)+hdimens,n) )
          expec_p_d(mu,i) = sum_c1 !expec_p_d(mu,i) + f_n_d(n)*evec_isigmamu_cong*pauli_mat_d(sigma,sigmap,4)*evec_isigmamu
          expec_z_d(mu,i) = sum_r2 !expec_z_d(mu,i) + f_n_d(n)*real( evec_isigmamu_cong*pauli_mat_d(sigma,sigmap,3)*evec_isigmamu )

        end do
      end do
    else
      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          lambda_d(mu,i) = s%Types(s%Basis(i)%Material)%lambda(mu)
        end do
      end do

      hdimens=dimens/2

      !$acc parallel loop
      do n = 1,dimens
        f_n_negative_d(n) = fd_dist_gpu(fermi, beta, -eval_d(n))
        tanh_n_d(n) = tanh(eval_d(n)*beta*0.5_dp)
      end do
      !$acc end parallel loop

      !$cuf kernel do(2) <<< (1,*), (9,*) >>>
      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          sum_r1 = 0._dp
          sum_r2 = 0._dp
          do n = 1,dimens
            ! up spin (using u's) + down spin (using v's)
            sum_r1 = sum_r1 + f_n_d(n)*real( conjg(hk_d(isigmamu2n_d(i,1,mu),n))*hk_d(isigmamu2n_d(i,1,mu),n) ) + f_n_negative_d(n)*real( conjg(hk_d(isigmamu2n_d(i,2,mu)+hdimens,n))*hk_d(isigmamu2n_d(i,2,mu)+hdimens,n) )

            sum_r2 = sum_r2 + 0.5_dp*lambda_d(mu,i)*tanh_n_d(n)*real( conjg(hk_d(isigmamu2n_d(i,1,mu)+hdimens,n))*hk_d(isigmamu2n_d(i,2,mu),n) )
          end do
          expec_0_d(mu,i) = sum_r1 !expec_0_d(mu,i) + f_n_d(n)*real( conjg(hk_d(isigmamu2n_d(i,1,mu),n))*hk_d(isigmamu2n_d(i,1,mu),n) ) + f_n_negative_d(n)*real( conjg(hk_d(isigmamu2n_d(i,2,mu)+hdimens,n))*hk_d(isigmamu2n_d(i,2,mu)+hdimens,n) )
          expec_d_d(mu,i) = sum_r2 !expec_d_d(mu,i) + 0.5_dp*lambda_d(mu,i)*tanh_n_d(n)*real( conjg(hk_d(isigmamu2n_d(i,1,mu)+hdimens,n))*hk_d(isigmamu2n_d(i,2,mu),n) )
        end do
      end do


      !$cuf kernel do(2) <<< (1,*), (9,*) >>>
      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          sum_c1 = cZero
          sum_r1 = 0._dp
          do n = 1,dimens
            do sigma = 1,2
              do sigmap = 1,2
                evec_isigmamu_cong = conjg( hk_d(isigmamu2n_d(i,sigma,mu),n) ) * hk_d(isigmamu2n_d(i,sigmap,mu),n)

                ! M_p
                sum_c1 = sum_c1 + f_n_d(n)*evec_isigmamu_cong*pauli_mat_d(sigma,sigmap,4)

                ! M_z
                sum_r1 = sum_r1 + f_n_d(n)*real( evec_isigmamu_cong*pauli_mat_d(sigma,sigmap,3) )
              end do
            end do
          end do
          expec_p_d(mu,i) = sum_c1 !expec_p_d(mu,i) + f_n_d(n)*evec_isigmamu_cong*pauli_mat_d(sigma,sigmap,4)*evec_isigmamu
          expec_z_d(mu,i) = sum_r1 !expec_z_d(mu,i) + f_n_d(n)*real( evec_isigmamu_cong*pauli_mat_d(sigma,sigmap,3)*evec_isigmamu )
        end do
      end do
    end if

  end subroutine expec_val_gpu
#endif


  subroutine expec_val(s,dimens,hk,eval,expec_0,expec_p,expec_z,expec_d)
  !> Subroutine to calculate the expectation value of the operators
  !> 1 (occupation), Sp, Sz, and superconducting delta
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,pi,pauli_mat
    use mod_parameters,        only: eta,isigmamu2n
    use mod_distributions,     only: fd_dist
    use mod_system,            only: System_type
    use mod_superconductivity, only: lsuperCond
    implicit none
    integer,                                 intent(in)  :: dimens
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(dimens),          intent(in)  :: eval
    complex(dp), dimension(dimens,dimens),   intent(in)  :: hk
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: expec_0, expec_z
    complex(dp), dimension(s%nOrb,s%nAtoms), intent(out) :: expec_p
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: expec_d

    real(dp)    :: fermi,beta
    integer     :: i,n,sigma,sigmap,mu,hdimens
    real(dp)    :: f_n(dimens),f_n_negative(dimens),tanh_n(dimens)
    complex(dp) :: evec_isigmamu, evec_isigmamu_cong !dimens = 2*nOrb*nAtoms

    beta = 1._dp/(pi*eta)
    expec_0 = 0._dp
    expec_z = 0._dp
    expec_p = cZero
    expec_d = 0._dp

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)

    do n = 1,dimens
      f_n(n) = fd_dist(fermi, beta, eval(n))
    end do

    if(.not.lsupercond) then
      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          do n = 1,dimens
            do sigma = 1,2
              evec_isigmamu = hk(isigmamu2n(i,sigma,mu),n)
              evec_isigmamu_cong = conjg( evec_isigmamu )

              ! Charge
              expec_0(mu,i) = expec_0(mu,i) + f_n(n)*real( evec_isigmamu_cong*evec_isigmamu )

              do sigmap = 1, 2
                evec_isigmamu = evec_isigmamu_cong*hk(isigmamu2n(i,sigmap,mu),n)
                ! M_p
                expec_p(mu,i) = expec_p(mu,i) + f_n(n)*pauli_mat(sigma,sigmap,4)*evec_isigmamu

                ! M_z
                expec_z(mu,i) = expec_z(mu,i) + f_n(n)*real( pauli_mat(sigma,sigmap,3)*evec_isigmamu )
              end do
            end do
          end do
        end do
      end do
    else
      hdimens=dimens/2

      do n = 1,dimens
        f_n_negative(n) = fd_dist(fermi, beta, -eval(n))
        tanh_n(n) = tanh(eval(n)*beta*0.5_dp)
      end do

      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          do n = 1,dimens
            ! up spin (using u's) + down spin (using v's)
            expec_0(mu,i) = expec_0(mu,i) + f_n(n)*real( conjg(hk(isigmamu2n(i,1,mu),n))*hk(isigmamu2n(i,1,mu),n) ) + f_n_negative(n)*real( conjg(hk(isigmamu2n(i,2,mu)+hdimens,n))*hk(isigmamu2n(i,2,mu)+hdimens,n) )

            expec_d(mu,i) = expec_d(mu,i) + 0.5_dp*s%Types(s%Basis(i)%Material)%lambda(mu)*tanh_n(n)*real( conjg(hk(isigmamu2n(i,1,mu)+hdimens,n))*hk(isigmamu2n(i,2,mu),n) )
          end do
        end do
      end do

      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          do n = 1,dimens
            do sigma = 1,2
              do sigmap = 1,2
                evec_isigmamu_cong = conjg( hk(isigmamu2n(i,sigma,mu),n) )*hk(isigmamu2n(i,sigmap,mu),n)

                ! M_p
                expec_p(mu,i) = expec_p(mu,i) + f_n(n)*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)

                ! M_z
                expec_z(mu,i) = expec_z(mu,i) + f_n(n)*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3) )
              end do
            end do
          end do
        end do
      end do
    end if

  end subroutine expec_val


  subroutine expec_val_n(s, dimens, evec, eval, expec_0, expec_p, expec_z, expec_d)
  !> Calculate the expectation value of the operators
  !> 1 (occupation), Sp, Sz, and superconducting delta
  !> for a given state [n] ([evec]) with eigen-energy [eval]
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,pi,pauli_mat
    use mod_parameters,        only: eta,isigmamu2n
    use mod_distributions,     only: fd_dist
    use mod_system,            only: System_type
    use mod_superconductivity, only: lsuperCond
    implicit none

    type(System_type),                       intent(in)  :: s
    integer,                                 intent(in)  :: dimens
    complex(dp), dimension(dimens),          intent(in)  :: evec
    real(dp),                                intent(in)  :: eval
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: expec_0, expec_z
    complex(dp), dimension(s%nOrb,s%nAtoms), intent(out) :: expec_p
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(out) :: expec_d

    integer     :: i,sigma,sigmap,mu,hdimens
    real(dp)    :: f_n,f_n_negative,tanh_n
    real(dp)    :: fermi,beta
    complex(dp) :: evec_isigmamu,evec_isigmamu_cong

    beta = 1._dp/(pi*eta)
    expec_0 = 0._dp
    expec_z = 0._dp
    expec_p = cZero
    expec_d = 0._dp

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    ! Fermi-Dirac:
    f_n = fd_dist(fermi, beta, eval)

    if(.not.lsupercond) then
      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          do sigma = 1,2

            evec_isigmamu = evec(isigmamu2n(i,sigma,mu))
            evec_isigmamu_cong = conjg( evec_isigmamu )

            ! Charge
            expec_0(mu,i) = expec_0(mu,i) + f_n*real( evec_isigmamu_cong*evec_isigmamu )

            do sigmap = 1, 2
              evec_isigmamu = evec_isigmamu_cong*evec(isigmamu2n(i,sigmap,mu))
              ! M_p
              expec_p(mu,i) = expec_p(mu,i) + f_n*pauli_mat(sigma,sigmap,4)*evec_isigmamu

              ! M_z
              expec_z(mu,i) = expec_z(mu,i) + f_n*real( pauli_mat(sigma,sigmap,3)*evec_isigmamu )
            end do
          end do
        end do
      end do
    else
      hdimens=dimens/2

      f_n_negative = fd_dist(fermi, beta, -eval)
      tanh_n = tanh(eval*beta*0.5_dp)

      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          ! up spin (using u's) + down spin (using v's)
          expec_0(mu,i) = expec_0(mu,i) + f_n*real( conjg(evec(isigmamu2n(i,1,mu)))*evec(isigmamu2n(i,1,mu)) ) + f_n_negative*real( conjg(evec(isigmamu2n(i,2,mu)+hdimens))*evec(isigmamu2n(i,2,mu)+hdimens) )

          expec_d(mu,i) = expec_d(mu,i) + 0.5_dp*s%Types(s%Basis(i)%Material)%lambda(mu)*tanh_n*real( conjg(evec(isigmamu2n(i,1,mu)+hdimens))*evec(isigmamu2n(i,2,mu)) )
        end do
      end do

      do i = 1,s%nAtoms
        do mu = 1,s%nOrb
          do sigma = 1,2
            do sigmap = 1,2
              evec_isigmamu_cong = conjg( evec(isigmamu2n(i,sigma,mu)) )*evec(isigmamu2n(i,sigmap,mu))

              ! M_p
              expec_p(mu,i) = expec_p(mu,i) + f_n*evec_isigmamu_cong*pauli_mat(sigma,sigmap,4)

              ! M_z
              expec_z(mu,i) = expec_z(mu,i) + f_n*real( evec_isigmamu_cong*pauli_mat(sigma,sigmap,3) )
            end do
          end do
        end do
      end do
    end if

  end subroutine expec_val_n


  subroutine expec_torque_n(s,dimens,evec,eval,tso_op,txc_op,expec_tso,expec_txc)
  !> Calculate the expectation value of the operators
  !> 1 (occupation), Sp, Sz, and superconducting delta
  !> for a given state [n] ([evec]) with eigen-energy [eval]
    use mod_kind,              only: dp
    use mod_constants,         only: pi
    use mod_parameters,        only: eta,isigmamu2n
    use mod_distributions,     only: fd_dist
    use mod_system,            only: System_type
    use mod_superconductivity, only: lsuperCond
    implicit none

    type(System_type),                                  intent(in)  :: s
    integer,                                            intent(in)  :: dimens
    complex(dp), dimension(dimens),                     intent(in)  :: evec
    real(dp),                                           intent(in)  :: eval
    complex(dp), dimension(s%nOrb2,s%nOrb2,3,s%nAtoms), intent(in)  :: tso_op,txc_op
    real(dp),    dimension(3,2,s%nAtoms),               intent(out) :: expec_tso
    real(dp),    dimension(3,s%nOrb,s%nAtoms),          intent(out) :: expec_txc

    integer     :: i,sigma,sigmap,mu,nu,mup,nup,mud,nud,mus,nus,xyz
    real(dp)    :: f_n
    real(dp)    :: fermi,beta
    complex(dp) :: evec_isigmamu,evec_isigmapnu_cong

    beta = 1._dp/(pi*eta)
    expec_tso = 0._dp
    expec_txc = 0._dp

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    ! Fermi-Dirac:
    f_n = fd_dist(fermi, beta, eval)

    sites_loop: do i = 1,s%nAtoms
      ! <Tso> for s orbitals is zero L(l=0)=0

      ! <Tso> for p orbitals:
      do mup = 1,s%npOrb
        mu = s%pOrbs(mup)

        do sigma = 1,2
          mus = (sigma-1)*s%nOrb+mu
          evec_isigmamu = conjg( evec(isigmamu2n(i,sigma,mu)) )

          do nup = 1,s%npOrb
            nu = s%pOrbs(nup)
            do sigmap = 1,2
              nus = (sigmap-1)*s%nOrb+nu
              evec_isigmapnu_cong = evec_isigmamu*evec(isigmamu2n(i,sigmap,nu))
              do xyz = 1,3
                expec_tso(xyz,1,i) = expec_tso(xyz,1,i) + f_n*real( tso_op(mus,nus,xyz,i)*evec_isigmapnu_cong )
              end do
            end do
          end do
        end do
      end do

      ! <Tso> for d orbitals:
      do mud = 1,s%ndOrb
        mu = s%dOrbs(mud)

        do sigma = 1,2
          mus = (sigma-1)*s%nOrb+mu
          evec_isigmamu = conjg( evec(isigmamu2n(i,sigma,mu)) )

          do nud = 1,s%ndOrb
            nu = s%dOrbs(nud)
            do sigmap = 1,2
              nus = (sigmap-1)*s%nOrb+nu
              evec_isigmapnu_cong = evec_isigmamu*evec(isigmamu2n(i,sigmap,nu))
              do xyz = 1,3
                expec_tso(xyz,2,i) = expec_tso(xyz,2,i) + f_n*real( tso_op(mus,nus,xyz,i)*evec_isigmapnu_cong )
              end do
            end do
          end do
        end do
      end do

      ! <Txc> for each orbital:
      do mu = 1,s%nOrb
        do sigma = 1,2
          mus = (sigma-1)*s%nOrb+mu
          evec_isigmamu = conjg( evec(isigmamu2n(i,sigma,mu)))
          do sigmap = 1, 2
            nus = (sigmap-1)*s%nOrb+mu
            evec_isigmapnu_cong = evec_isigmamu*evec(isigmamu2n(i,sigmap,mu))

            do xyz = 1,3
              ! T_xc
              expec_txc(xyz,mu,i) = expec_txc(xyz,mu,i) + f_n*real( txc_op(mus,nus,xyz,i)*evec_isigmapnu_cong )
            end do
          end do
        end do
      end do
    end do sites_loop
  end subroutine expec_torque_n


  subroutine groundstate_L_and_E()
  !> Wrapper for the calculation of the expectation value
  !> of the orbital angular momentum and the band energy in the ground state
    use mod_kind,          only: dp
    use mod_parameters,    only: output
    use mod_constants,     only: rad2deg
    use mod_System,        only: s => sys
    use mod_mpi_pars,      only: rField,abortProgram
    use mod_magnet,        only: labs,ltheta,lphi,lpabs,lptheta,lpphi,lxm,lym,lzm,lxpm,lypm,lzpm
    implicit none
    integer      :: i,AllocateStatus

    if(rField == 0) &
      write(output%unit_loop,"('[groundstate_L_and_E] Calculating Orbital Angular Momentum and Band energy of ground state... ')")

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
      call abortProgram("[groundstate_L_and_E] Not enough memory for: lxm,lym,lzm,lxpm,lypm,lzpm")

    call calc_GS_L_and_E()

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

  end subroutine groundstate_L_and_E


  subroutine calc_GS_L_and_E_greenfunction()
  !> Calculates the expectation value of the orbital angular momentum
  !> in the ground state using green functions
    use mod_kind,          only: dp,int64
    use mod_constants,     only: cZero,pi
    use mod_System,        only: s => sys
    use mod_parameters,    only: eta,output
    use EnergyIntegration, only: y, wght
    use mod_magnet,        only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz
    use mod_hamiltonian,   only: hamilt_local,energy
    use mod_greenfunction, only: green
    use adaptiveMesh,      only: local_points,activeComm,E_k_imag_mesh,bzs
    use mod_superconductivity, only: lsuperCond
    use mod_SOC,               only: llinearsoc,llineargfsoc
    use mod_mpi_pars,      only: abortProgram,MPI_IN_PLACE,MPI_DOUBLE_COMPLEX,MPI_SUM,ierr,myrank
    implicit none
    integer(int64)    :: ix

    integer      :: i,mu,nu,mup,nup
    real(dp) :: kp(3)
    real(dp) :: weight, ep, fermi
    complex(dp), dimension(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms) :: gf
    complex(dp), dimension(s%nOrb,s%nOrb,s%nAtoms)            :: gupgd

    integer :: ncount

    external :: MPI_Allreduce

    ncount=s%nAtoms*s%nOrb*s%nOrb

    ! Build local hamiltonian
    if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)

    if(myrank == 0) &
      write(output%unit, "('[Warning] [calc_GS_L_and_E_greenfunction] Band energy not implemented with greenfunctions.')")
    energy = 0._dp

    gupgd  = cZero
    ! Calculating the jacobian using a complex integral
    !$omp parallel default(none) &
    !$omp& private(ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(local_points,fermi,s,lsupercond,E_k_imag_mesh,bzs,eta,y,wght,gupgd)


    gf = cZero
    !$omp do schedule(dynamic) reduction(+:gupgd)
    do ix = 1, local_points
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      ep = y(E_k_imag_mesh(1,ix))
      weight = bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix)) * wght(E_k_imag_mesh(1,ix))
      !Green function on energy Ef + iy, and wave vector kp
      call green(fermi,ep+eta,s,kp,gf)

      site_i: do i=1,s%nAtoms
        orb_mu: do mu=1,s%nOrb
          mup = mu+s%nOrb
          orb_nu: do nu=1,s%nOrb
            nup = nu+s%nOrb
            gupgd(mu,nu,i) = gupgd(mu,nu,i) + (gf(mu,nu,i,i) + gf(mup,nup,i,i)) * weight
          end do orb_nu
        end do orb_mu
      end do site_i
    end do
    !$omp end do
    !$omp end parallel

    gupgd = gupgd / pi
    call MPI_Allreduce(MPI_IN_PLACE, gupgd, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    lxpm = 0._dp
    lypm = 0._dp
    lzpm = 0._dp
    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp

    do nu=1,s%nOrb
      do mu=1,s%nOrb
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

  end subroutine calc_GS_L_and_E_greenfunction


  subroutine expec_L_n(s,dimens,evec,eval,lxm,lym,lzm)
  !> Calculate the expectation value of the orbital momentum
  !> for a given state [n] ([evec]) with eigen-energy [eval]
    use mod_kind,              only: dp
    use mod_constants,         only: pi
    use mod_parameters,        only: eta,isigmamu2n
    use mod_System,            only: System_type
    use mod_magnet,            only: lx,ly,lz
    use mod_distributions,     only: fd_dist
    use mod_superconductivity, only: lsuperCond
    implicit none
    type(System_type),                  intent(in)  :: s
    real(dp),                           intent(in)  :: eval
    integer,                            intent(in)  :: dimens
    real(dp),    dimension(2,s%nAtoms), intent(out) :: lxm,lym,lzm
    complex(dp), dimension(dimens),     intent(in)  :: evec
    integer                                         :: i,mu,nu,mup,nup,mud,nud,sigma
    real(dp)                                        :: f_n,fermi,beta
    complex(dp)                                     :: prod

    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    beta = 1._dp/(pi*eta)

    ! Fermi-Dirac:
    f_n = fd_dist(fermi, beta, eval)

    sites_loop: do i = 1,s%nAtoms
      ! <L> for s orbitals is zero (l=0)

      ! <L> for p orbitals:
      do nup = 1,s%npOrb
        nu = s%pOrbs(nup)
        do mup = 1,s%npOrb
          mu = s%pOrbs(mup)
          do sigma = 1,2
            prod = f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*evec(isigmamu2n(i,sigma,nu))
            lxm (1,i) = lxm (1,i) + real( prod*lx (mu,nu) ) !> angular momentum at atomic site (i)
            lym (1,i) = lym (1,i) + real( prod*ly (mu,nu) )
            lzm (1,i) = lzm (1,i) + real( prod*lz (mu,nu) )
          end do
        end do
      end do

      ! <L> for d orbitals:
      do nud = 1,s%ndOrb
        nu = s%dOrbs(nud)
        do mud = 1,s%ndOrb
          mu = s%dOrbs(mud)
          do sigma = 1,2
            prod = f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*evec(isigmamu2n(i,sigma,nu))
            lxm (2,i) = lxm (2,i) + real( prod*lx (mu,nu) ) !> angular momentum at atomic site (i)
            lym (2,i) = lym (2,i) + real( prod*ly (mu,nu) )
            lzm (2,i) = lzm (2,i) + real( prod*lz (mu,nu) )
          end do
        end do
      end do
    end do sites_loop

  end subroutine expec_L_n


  subroutine expec_H_n(s,b_field,A_t,hk,kp,evec,eval,E_0)
  !> Calculate the expectation value of the time-dependent Hamiltonian
  !> in the propagated states
    use mod_kind,          only: dp
    use mod_constants,     only: pi
    use mod_parameters,    only: eta,dimH
    use mod_System,        only: System_type
    use mod_distributions, only: fd_dist
    use mod_hamiltonian,   only: build_hext
    implicit none
    type(System_type),                 intent(in) :: s
    complex(dp), dimension(dimH),      intent(in) :: evec
    complex(dp), dimension(dimH,dimH), intent(in) :: hk
    real(dp),                          intent(in) :: eval
    real(dp),    dimension(3),         intent(in) :: b_field,A_t,kp
    integer     :: i, j
    real(dp)    :: f_n, expec_H_0, E_0
    complex(dp) :: hext_t(dimH,dimH),hamilt_0(dimH,dimH)

    ! Fermi-Dirac:
    f_n = fd_dist(s%Ef, 1._dp/(pi*eta), eval)

    ! Building time dependent hamiltonian
    call build_hext(kp,b_field,A_t,hext_t)
    hamilt_0 = hk + hext_t

    E_0 = 0._dp

    do i=1, dimH
      do j=1, dimH
        expec_H_0 = real( conjg( evec(i) ) * hamilt_0(i,j) * evec(j) )
        E_0       =  E_0 + f_n * expec_H_0
      end do
    end do

  end subroutine expec_H_n


  subroutine calc_GS_L_and_E_eigenstates()
  !> Calculates the expectation value of the orbital angular momentum
  !> and the band energy in the ground state
    use mod_kind,              only: dp,int64
    use mod_BrillouinZone,     only: realBZ
    use mod_constants,         only: pi,cZero
    use mod_parameters,        only: dimHsc,eta,isigmamu2n
    use mod_System,            only: s => sys
    use mod_magnet,            only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz
    use mod_distributions,     only: fd_dist
    use mod_tools,             only: diagonalize,lwork
    use mod_superconductivity, only: lsuperCond
    use mod_hamiltonian,       only: calchk,h0,hamilt_local,energy
    use mod_mpi_pars,          only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    integer(int64)                 :: iz
    integer                        :: n,i,mu,nu,sigma,ilaenv
    real(dp)                       :: fermi,beta
    real(dp),    dimension(dimHsc) :: eval,f_n
    complex(dp)                    :: hk(dimHsc,dimHsc),prod(s%nOrb,s%nOrb,s%nAtoms)

    external :: MPI_Allreduce,ilaenv

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    beta = 1._dp/(pi*eta)

    ! Getting lwork for diagonalization
    lwork = (ilaenv( 1, 'zhetrd', 'VU', dimHsc, -1, -1, -1 )+1)*dimHsc

    call hamilt_local(s)

    prod = cZero
    energy = 0._dp

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp& private(iz,n,f_n,i,sigma,mu,nu,hk,eval) &
    !$omp& shared(s,dimHsc,realBZ,h0,fermi,beta,isigmamu2n) &
    !$omp& reduction(+:prod,energy)
    kloop: do iz = 1,realBZ%workload
      ! Calculating the hamiltonian for a given k-point
      hk = h0 + calchk(s,realBZ%kp(1:3,iz))

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call diagonalize(dimHsc,hk,eval)

      do n = 1,dimHsc
        ! Fermi-Dirac distrubution
        f_n(n) = fd_dist(fermi, beta, eval(n))
        energy = energy + f_n(n)*eval(n)*realBZ%w(iz)
      end do

      do i = 1,s%nAtoms
        do nu = 1,s%nOrb
          do mu = 1,s%nOrb
            do sigma=1,2
              do n = 1,dimHsc
                prod(mu,nu,i) = prod(mu,nu,i) + f_n(n)*conjg( hk(isigmamu2n(i,sigma,mu),n) )*hk(isigmamu2n(i,sigma,nu),n)*realBZ%w(iz)
              end do
            end do
          end do
        end do
      end do

    end do kloop
    !$omp end parallel do

    call MPI_Allreduce(MPI_IN_PLACE, prod   , s%nOrb*s%nOrb*s%nAtoms, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, energy , 1                     , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    ! Building different components of the orbital angular momentum
    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp
    lxpm = 0._dp
    lypm = 0._dp
    lzpm = 0._dp
    do i = 1,s%nAtoms
      do nu = 1,s%nOrb
        do mu = 1,s%nOrb
          lxm (i) = lxm (i) + real( prod(mu,nu,i)*lx (mu,nu  ) )
          lym (i) = lym (i) + real( prod(mu,nu,i)*ly (mu,nu  ) )
          lzm (i) = lzm (i) + real( prod(mu,nu,i)*lz (mu,nu  ) )
          lxpm(i) = lxpm(i) + real( prod(mu,nu,i)*lxp(mu,nu,i) )
          lypm(i) = lypm(i) + real( prod(mu,nu,i)*lyp(mu,nu,i) )
          lzpm(i) = lzpm(i) + real( prod(mu,nu,i)*lzp(mu,nu,i) )
        end do
      end do
    end do

  end subroutine calc_GS_L_and_E_eigenstates

#ifdef _GPU
  subroutine calc_GS_L_and_E_fullhk_gpu()
  !> Calculates the expectation value of the orbital angular momentum
  !> and the band energy in the ground state using the full hamiltonian matrix
    use mod_kind,              only: dp,int64
    use mod_BrillouinZone,     only: realBZ
    use mod_constants,         only: pi,cZero
    use mod_parameters,        only: dimHsc,eta,isigmamu2n_d
    use mod_System,            only: s => sys
    use mod_magnet,            only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz
    use mod_distributions,     only: fd_dist_gpu
    use mod_cuda,              only: h,diagonalize_gpu
    use mod_superconductivity, only: lsuperCond
    use mod_hamiltonian,       only: hamilt_local_gpu,h0_d,fullhk_d,energy
    use mod_mpi_pars,          only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    integer(int64)                                     :: iz
    integer                                            :: n,i,j,mu,nu,sigma
    real(dp)                                           :: fermi,beta
    complex(dp), dimension(:,:,:), allocatable         :: prod
    complex(dp), dimension(:,:,:), allocatable, device :: prod_d
    real(dp),    dimension(dimHsc),             device :: f_n_d
    complex(dp), dimension(:,:),   allocatable, device :: hk_d
    real(dp),    dimension(:),     allocatable, device :: eval_d
    complex(dp) :: sum_c
    real(dp)    :: weight_d

    external :: MPI_Allreduce

    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    beta = 1._dp/(pi*eta)

    allocate(hk_d(dimHsc,dimHsc),eval_d(dimHsc))
    allocate(prod(s%nOrb,s%nOrb,s%nAtoms))
    allocate(prod_d(s%nOrb,s%nOrb,s%nAtoms))

    call hamilt_local_gpu(s)

    prod_d = cZero
    energy = 0._dp

    kloop: do iz = 1,realBZ%workload
      weight_d = realBZ%w(iz)

      !$cuf kernel do(2) <<< *, * >>>
      do j = 1, dimHsc
        do i = 1, dimHsc
          hk_d(i,j) = h0_d(i,j) + fullhk_d(i,j,iz)
        end do
      end do
      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call diagonalize_gpu(dimHsc,hk_d,eval_d,h)

      !$acc parallel loop
      do n = 1,dimHsc
        f_n_d(n) = fd_dist_gpu(fermi, beta, eval_d(n))
        energy = energy + f_n_d(n) * eval_d(n) * weight_d
      end do
      !$acc end parallel loop

      !$cuf kernel do(3) <<< (1,1,*), (9,9,*) >>>
      do i = 1,s%nAtoms
        do nu = 1,s%nOrb
          do mu = 1,s%nOrb
            sum_c = cZero
            do sigma=1,2
              do n = 1,dimHsc
                sum_c = sum_c + f_n_d(n)*conjg( hk_d(isigmamu2n_d(i,sigma,mu),n) )*hk_d(isigmamu2n_d(i,sigma,nu),n)*weight_d
              end do
            end do
            prod_d(mu,nu,i) = sum_c
          end do
        end do
      end do
    end do kloop

    call MPI_Allreduce(MPI_IN_PLACE, prod_d , s%nOrb*s%nOrb*s%nAtoms, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, energy , 1                     , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    prod = prod_d

    ! Building different components of the orbital angular momentum
    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp
    lxpm = 0._dp
    lypm = 0._dp
    lzpm = 0._dp
    do i = 1,s%nAtoms
      do nu = 1,s%nOrb
        do mu = 1,s%nOrb
          lxm (i) = lxm (i) + real( prod(mu,nu,i)*lx (mu,nu  ) )
          lym (i) = lym (i) + real( prod(mu,nu,i)*ly (mu,nu  ) )
          lzm (i) = lzm (i) + real( prod(mu,nu,i)*lz (mu,nu  ) )
          lxpm(i) = lxpm(i) + real( prod(mu,nu,i)*lxp(mu,nu,i) )
          lypm(i) = lypm(i) + real( prod(mu,nu,i)*lyp(mu,nu,i) )
          lzpm(i) = lzpm(i) + real( prod(mu,nu,i)*lzp(mu,nu,i) )
        end do
      end do
    end do

  deallocate(prod_d)
  deallocate(prod)
  deallocate(hk_d,eval_d)

  end subroutine calc_GS_L_and_E_fullhk_gpu
#else
  subroutine calc_GS_L_and_E_fullhk()
  !> Calculates the expectation value of the orbital angular momentum
  !> and the band energy in the ground state using the full hamiltonian matrix
    use mod_kind,              only: dp,int64
    use mod_BrillouinZone,     only: realBZ
    use mod_constants,         only: pi,cZero
    use mod_parameters,        only: dimHsc,eta,isigmamu2n
    use mod_System,            only: s => sys
    use mod_magnet,            only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz
    use mod_distributions,     only: fd_dist
    use mod_tools,             only: diagonalize,lwork
    use mod_superconductivity, only: lsuperCond
    use mod_hamiltonian,       only: hamilt_local,h0,fullhk,energy
    use mod_mpi_pars,          only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,ierr
    implicit none
    integer(int64)                 :: iz
    integer                        :: n,i,mu,nu,sigma,ilaenv
    real(dp)                       :: fermi,beta
    real(dp),    dimension(dimHsc) :: eval,f_n
    complex(dp)                    :: hk(dimHsc,dimHsc),prod(s%nOrb,s%nOrb,s%nAtoms)

    external :: MPI_Allreduce,ilaenv


    !If lsupercond is true then fermi is 0.0 otherwise is s%Ef
    fermi = merge(0._dp,s%Ef,lsuperCond)
    beta = 1._dp/(pi*eta)

    ! Getting lwork for diagonalization
    lwork = (ilaenv( 1, 'zhetrd', 'VU', dimHsc, -1, -1, -1 )+1)*dimHsc

    call hamilt_local(s)

    prod = cZero
    energy = 0._dp

    !$omp parallel do default(none) schedule(dynamic) &
    !$omp& private(iz,n,f_n,i,sigma,mu,nu,hk,eval) &
    !$omp& shared(s,dimHsc,h0,fullhk,realBZ,fermi,beta,isigmamu2n) &
    !$omp& reduction(+:prod,energy)
    kloop: do iz = 1,realBZ%workload
      ! hamiltonian for a given k-point
      hk = h0 + fullhk(:,:,iz)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call diagonalize(dimHsc,hk,eval)

      do n = 1,dimHsc
        f_n(n) = fd_dist(fermi, beta, eval(n))
        energy = energy + f_n(n)*eval(n)*realBZ%w(iz)
      end do

      do i = 1,s%nAtoms
        do nu = 1,s%nOrb
          do mu = 1,s%nOrb
            do sigma=1,2
              do n = 1,dimHsc
                prod(mu,nu,i) = prod(mu,nu,i) + f_n(n)*conjg( hk(isigmamu2n(i,sigma,mu),n) )*hk(isigmamu2n(i,sigma,nu),n)*realBZ%w(iz)
              end do
            end do
          end do
        end do
      end do

    end do kloop
    !$omp end parallel do

    call MPI_Allreduce(MPI_IN_PLACE, prod   , s%nOrb*s%nOrb*s%nAtoms, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, energy , 1                     , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

    ! Building different components of the orbital angular momentum
    lxm  = 0._dp
    lym  = 0._dp
    lzm  = 0._dp
    lxpm = 0._dp
    lypm = 0._dp
    lzpm = 0._dp
    do i = 1, s%nAtoms
      do nu = 1,s%nOrb
        do mu = 1,s%nOrb
          lxm (i) = lxm (i) + real( prod(mu,nu,i)*lx (mu,nu  ) )
          lym (i) = lym (i) + real( prod(mu,nu,i)*ly (mu,nu  ) )
          lzm (i) = lzm (i) + real( prod(mu,nu,i)*lz (mu,nu  ) )
          lxpm(i) = lxpm(i) + real( prod(mu,nu,i)*lxp(mu,nu,i) )
          lypm(i) = lypm(i) + real( prod(mu,nu,i)*lyp(mu,nu,i) )
          lzpm(i) = lzpm(i) + real( prod(mu,nu,i)*lzp(mu,nu,i) )
        end do
      end do
    end do

  end subroutine calc_GS_L_and_E_fullhk

#endif

  subroutine calc_E_dc()
    use mod_kind,        only: dp
    use mod_parameters,  only: output
    use mod_hamiltonian, only: energy_dc,energy_dc_n,energy_constr
    use mod_system,      only: s => sys
    use mod_magnet,      only: lconstraining_field,mzd0,mpd0,bc,mabsd,m_fix,mx,my,mz,rhod,rhod0,rho,rho0
    use mod_mpi_pars,    only: rField
    implicit none
    real(dp), dimension(3,s%nAtoms) :: m_tot
    integer                         :: i,mu,mud

    if(rField==0) &
      write(output%unit_loop,"('[calc_E_dc] Calculating double-counting energies... ')")

    energy_dc     = 0.0_dp
    energy_dc_n   = 0.0_dp
    energy_constr = 0.0_dp
    do i=1,s%nAtoms
      energy_dc   = energy_dc   + 0.25_dp*s%Basis(i)%Um*(mabsd(i)**2-mzd0(i)**2-abs(mpd0(i))**2)
      energy_dc_n = energy_dc_n + 0.25_dp*s%Basis(i)%Un*(rhod(i)**2-rhod0(i)**2)
      do mud=1,s%ndOrb
        mu = s%dOrbs(mud)
        energy_dc_n = energy_dc_n - 0.5_dp*s%Basis(i)%Un*(rho(mu,i)**2-rho0(mu,i)**2)
      end do
    end do
    if(lconstraining_field) then
      do i=1,s%nAtoms
        m_tot(1,i)= sum(mx(:,i))
        m_tot(2,i)= sum(my(:,i))
        m_tot(3,i)= sum(mz(:,i))
        energy_constr = energy_constr - dot_product(m_tot(:,i),bc(:,i)) + dot_product(m_tot(:,i),m_fix(:,i))*dot_product(m_fix(:,i),bc(:,i))
      end do
      if((energy_constr > 1.e-7_dp).and.(rField==0)) &
        write(output%unit_loop,"('[calc_E_dc] Constraninig energy not 0! energy_constr = ',es16.9)") energy_constr

    end if

  end subroutine calc_E_dc

end module mod_expectation
