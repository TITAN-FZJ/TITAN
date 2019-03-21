!! This module contains the calculation of the ground state expectation values using the GF
!! and using the eigenstates of the hamiltonian
module mod_expectation
  implicit none

contains
  subroutine expectation_values_greenfunction(s,rho,mp,mx,my,mz)
    !! Calculates ground state (occupation and magnetization) quantities using the Green functions
    use mod_f90_kind,      only: double
    use mod_constants,     only: cI,pi,cZero
    use mod_SOC,           only: llinearsoc,llineargfsoc
    use EnergyIntegration, only: y,wght
    use mod_system,        only: System
    use adaptiveMesh,      only: bzs,E_k_imag_mesh,activeComm,local_points
    use TightBinding,      only: nOrb,nOrb2
    use mod_parameters,    only: eta
    use ElectricField,     only: EshiftBZ,ElectricFieldVector
    use mod_mpi_pars
    implicit none
    type(System),                              intent(in)  :: s
    real(double),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    complex(double), dimension(nOrb,s%nAtoms), intent(out) :: mp

    integer  :: i,j, AllocateStatus
    real(double),    dimension(3)                    :: kp
    complex(double), dimension(:,:),     allocatable :: gdiagud,gdiagdu
    real(double),    dimension(:,:),     allocatable :: imguu,imgdd
    complex(double), dimension(:,:,:,:), allocatable :: gf
    !--------------------- begin MPI vars --------------------
    integer*8 :: ix
    integer :: ncount
    integer :: mu,mup
    real(double) :: weight, ep
    ncount = s%nAtoms * nOrb

    allocate(imguu(nOrb,s%nAtoms),imgdd(nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) &
    call abortProgram("[expectation_values_greenfunction] Not enough memory for: imguu,imgdd")

    allocate(gdiagud(s%nAtoms,nOrb), gdiagdu(s%nAtoms,nOrb), stat = AllocateStatus)
    if(AllocateStatus /= 0) &
    call abortProgram("[expectation_values_greenfunction] Not enough memory for: gdiagdu, gdiagud")

    imguu   = 0.d0
    imgdd   = 0.d0
    gdiagud = cZero
    gdiagdu = cZero

    !$omp parallel default(none) &
    !$omp& private(ix,ep,kp,weight,i,mu,mup,gf,AllocateStatus) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,eta,wght,s,nOrb,nOrb2,bzs,E_k_imag_mesh,y,gdiagud,gdiagdu,imguu,imgdd,EshiftBZ,ElectricFieldVector)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat=AllocateStatus)
    if(AllocateStatus /= 0) &
    call AbortProgram("[expectation_values_greenfunction] Not enough memory for: gf")
    gf = cZero

    if(llineargfsoc .or. llinearsoc) then
      !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix)) + EshiftBZ*ElectricFieldVector
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
      !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix)) + EshiftBZ*ElectricFieldVector
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         call green(s%Ef,ep+eta,s,kp,gf)
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
        imguu(mu,i) = 0.5d0 + imguu(mu,i)
        imgdd(mu,i) = 0.5d0 + imgdd(mu,i)
        rho(mu,i) = imguu(mu,i) + imgdd(mu,i)
        mz (mu,i) = imguu(mu,i) - imgdd(mu,i)
      end do
    end do

    deallocate(imguu,imgdd)
    deallocate(gdiagdu, gdiagud)
  end subroutine expectation_values_greenfunction


  !   Calculates ground state quantities from eigenstates
  subroutine expectation_values_eigenstates(s,rho,mp,mx,my,mz)
    use mod_f90_kind,       only: double
    use mod_BrillouinZone,  only: realBZ
    use mod_parameters,     only: output
    use mod_system,         only: System
    use TightBinding,       only: nOrb,nOrb2
    use ElectricField,      only: EshiftBZ,ElectricFieldVector
    implicit none
    type(System),                              intent(in)  :: s
    real(double),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    complex(double), dimension(nOrb,s%nAtoms), intent(out) :: mp

    integer                                      :: iz, info !, mu,i
    integer                                      :: lwork,dimH
    real(double)                                 :: weight, kp(3)
    real(double),    dimension(nOrb,s%nAtoms)    :: expec_0, expec_z
    complex(double), dimension(nOrb,s%nAtoms)    :: expec_p
    real(double),    dimension(:),  allocatable  :: rwork(:), eval(:)
    complex(double),                allocatable  :: work(:), hk(:,:)

    dimH  = (s%nAtoms)*nOrb2
    lwork = 21*dimH

    allocate( hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork) )

    !$omp parallel default(none) &
    !$omp& firstprivate(lwork) &
    !$omp& private(iz,kp,weight,hk,eval,work,rwork,info,expec_0, expec_p, expec_z) &
    !$omp& shared(s,dimH,output,nOrb2,realBZ,rho,mp,mz,EshiftBZ,ElectricFieldVector)

    rho = 0.d0
    mp  = 0.d0
    mz  = 0.d0

    !$omp do reduction(+:rho,mp,mz)
    do iz = 1,realBZ%workload
      kp = realBZ%kp(1:3,iz) + EshiftBZ*ElectricFieldVector
      weight = realBZ%w(iz)

      ! Calculating the hamiltonian for a given k-point
      call hamiltk(s,kp,hk)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)

      if(info/=0) then
        write(output%unit_loop,"('[expectation_values_eigenstates] Problem with diagonalization. info = ',i0)") info
        stop
      end if

      ! Calculating expectation values for a given k-point
      call expec_val(s, dimH, hk, eval, expec_0, expec_p, expec_z)
      rho = rho + expec_0*weight
      mp  = mp  + expec_p*weight
      mz  = mz  + expec_z*weight
    end do
    !$omp end do
    !$omp end parallel
    mx = real(mp)
    my = aimag(mp)

    deallocate(hk,rwork,eval,work)

  end subroutine expectation_values_eigenstates


  ! subroutine expectation value of the operators 1 (occupation), Sp and Sz:
  subroutine expec_val(s, dim, hk, eval, expec_0, expec_p, expec_z)
    use mod_f90_kind,      only: double 
    use mod_constants,     only: cOne,cZero,pi,pauli_mat
    use mod_parameters,    only: eta, isigmamu2n
    use mod_distributions, only: fd_dist
    use TightBinding,      only: nOrb
    use mod_system,        only: System
    implicit none
    integer,                                   intent(in)  :: dim
    type(System),                              intent(in)  :: s
    real(double),    dimension(dim),           intent(in)  :: eval
    complex(double), dimension(dim,dim),       intent(in)  :: hk
    real(double),    dimension(nOrb,s%nAtoms), intent(out) :: expec_0, expec_z
    complex(double), dimension(nOrb,s%nAtoms), intent(out) :: expec_p

    integer                                     :: i, n, sigma, sigmap, mu
    real(double)                                :: f_n
    complex(double)                             :: evec(dim)

    expec_0 = 0.d0
    expec_z = 0.d0
    expec_p = cZero
    do n = 1, dim
      ! Fermi-Dirac:
      f_n = fd_dist(s%Ef, 1.d0/(pi*eta), eval(n))

      ! Getting eigenvector and its transpose conjugate
      evec(:) = hk(:,n)

      do i = 1, s%nAtoms
        do mu = 1, nOrb
          do sigma = 1, 2
            ! Charge
            expec_0(mu,i) = expec_0(mu,i) + f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*evec(isigmamu2n(i,sigma,mu))

            do sigmap = 1, 2
              ! M_p
              expec_p(mu,i) = expec_p(mu,i) + f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*pauli_mat(sigma,sigmap,4)*evec(isigmamu2n(i,sigmap,mu))

              ! M_z
              expec_z(mu,i) = expec_z(mu,i) + f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*pauli_mat(sigma,sigmap,3)*evec(isigmamu2n(i,sigmap,mu))
            end do
          end do
        end do
      end do
    end do


  end subroutine expec_val


  ! subroutine expectation value of the operators 1 (occupation), Sp and Sz for a given state n:
  subroutine expec_val_n(s, dim, evec, eval, expec_0, expec_p, expec_z)
    use mod_f90_kind,      only: double 
    use mod_constants,     only: cOne,cZero,pi,pauli_mat
    use mod_parameters,    only: eta, isigmamu2n
    use mod_distributions, only: fd_dist
    use TightBinding,      only: nOrb
    use mod_system,        only: System
    implicit none
    
    type(System),                              intent(in)  :: s
    integer,                                   intent(in)  :: dim
    complex(double), dimension(dim),           intent(in)  :: evec
    real(double),                              intent(in)  :: eval
    real(double),    dimension(nOrb,s%nAtoms), intent(out) :: expec_0, expec_z
    complex(double), dimension(nOrb,s%nAtoms), intent(out) :: expec_p

    integer :: i, sigma, sigmap, mu
    real(double) :: f_n

    expec_0 = 0.d0
    expec_z = 0.d0
    expec_p = cZero
  
    ! Fermi-Dirac:
    f_n = fd_dist(s%Ef, 1.d0/(pi*eta), eval)

    do i = 1, s%nAtoms
      do mu = 1, nOrb
        do sigma = 1, 2
          ! Charge
          expec_0(mu,i) = expec_0(mu,i) + f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*evec(isigmamu2n(i,sigma,mu))

          do sigmap = 1, 2
            ! M_p
            expec_p(mu,i) = expec_p(mu,i) + f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*pauli_mat(sigma,sigmap,4)*evec(isigmamu2n(i,sigmap,mu))

            ! M_z
            expec_z(mu,i) = expec_z(mu,i) + f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*pauli_mat(sigma,sigmap,3)*evec(isigmamu2n(i,sigmap,mu))
          end do
        end do
      end do
    end do

  end subroutine expec_val_n

  
end module mod_expectation
