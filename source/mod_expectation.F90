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
    use mod_parameters,    only: nOrb,nOrb2,eta
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
    use mod_f90_kind,      only: double
    use mod_BrillouinZone, only: realBZ
    use mod_parameters,    only: nOrb,nOrb2,output
    use mod_system,        only: System
    use ElectricField,     only: EshiftBZ,ElectricFieldVector
    use mod_tools,         only: itos
    use mod_mpi_pars
    implicit none
    type(System),                              intent(in)  :: s
    real(double),    dimension(nOrb,s%nAtoms), intent(out) :: rho, mx, my, mz
    complex(double), dimension(nOrb,s%nAtoms), intent(out) :: mp

    integer                                      :: iz, info, ncount !, mu,i
    integer                                      :: lwork,dimH
    real(double)                                 :: weight, kp(3)
    real(double),    dimension(nOrb,s%nAtoms)    :: expec_0, expec_z
    complex(double), dimension(nOrb,s%nAtoms)    :: expec_p
    real(double),    dimension(:),  allocatable  :: rwork(:), eval(:)
    complex(double),                allocatable  :: work(:), hk(:,:)

    dimH  = (s%nAtoms)*nOrb2
    lwork = 21*dimH
    ncount = nOrb*s%nAtoms

    allocate( hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork) )

    !$omp parallel default(none) &
    !$omp& firstprivate(lwork) &
    !$omp& private(iz,kp,weight,hk,eval,work,rwork,info,expec_0, expec_p, expec_z) &
    !$omp& shared(s,dimH,output,realBZ,rho,mp,mz,EshiftBZ,ElectricFieldVector)
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

      if(info/=0) &
        call abortProgram("[expectation_values_eigenstates] Problem with diagonalization. info = " // itos(info))

      ! Calculating expectation values for a given k-point
      call expec_val(s, dimH, hk, eval, expec_0, expec_p, expec_z)

      rho = rho + expec_0*weight
      mp  = mp  + expec_p*weight
      mz  = mz  + expec_z*weight
    end do
    !$omp end do
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, rho, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mz , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)

    mx = real(mp)
    my = aimag(mp)

    deallocate(hk,rwork,eval,work)

  end subroutine expectation_values_eigenstates


  ! subroutine expectation value of the operators 1 (occupation), Sp and Sz:
  subroutine expec_val(s, dim, hk, eval, expec_0, expec_p, expec_z)
    use mod_f90_kind,      only: double 
    use mod_constants,     only: cOne,cZero,pi,pauli_mat
    use mod_parameters,    only: nOrb, eta, isigmamu2n
    use mod_distributions, only: fd_dist
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


  subroutine expec_val_n(s, dim, evec, eval, expec_0, expec_p, expec_z)
    !! Calculate the expectation value of the operators 1 (occupation), \sigma^+ and \sigma^z
    !! for a given state n (evec) with eigenenergy eval
    use mod_f90_kind,      only: double 
    use mod_constants,     only: cOne,cZero,pi,pauli_mat
    use mod_parameters,    only: nOrb, eta, isigmamu2n
    use mod_distributions, only: fd_dist
    use mod_system,        only: System
    implicit none
    
    type(System),                              intent(in)  :: s
    integer,                                   intent(in)  :: dim
    complex(double), dimension(dim),           intent(in)  :: evec
    real(double),                              intent(in)  :: eval
    real(double),    dimension(nOrb,s%nAtoms), intent(out) :: expec_0, expec_z
    complex(double), dimension(nOrb,s%nAtoms), intent(out) :: expec_p

    integer      :: i, sigma, sigmap, mu
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


  subroutine calcLGS()
    !! Calculates the expectation value of the orbital angular momentum in the ground state
    use mod_parameters,    only: output, leigenstates
    use mod_constants,     only: cZero,rad2deg
    use mod_System,        only: s => sys
    use mod_mpi_pars,      only: rField,abortProgram
    use mod_magnet
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

    if(leigenstates) then
      call calcLGS_eigenstates()
    else
      call calcLGS_greenfunction()
    end if

    ! Calculating angles of GS OAM (in degrees)
    do i = 1,s%nAtoms
      labs(i)   = sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
      if(labs(i)>1.d-8) then
        ltheta(i) = acos(lzm(i)/labs(i))*rad2deg
      else
        ltheta(i) = 0.d0
      end if
      if(abs(ltheta(i))>1.d-8) then
        if(abs(abs(ltheta(i))-180.d0)>1.d-8) then
          lphi(i)   = atan2(lym(i),lxm(i))*rad2deg
        else
          lphi(i) = 0.d0
        end if
      else
        lphi(i) = 0.d0
      end if
      lpabs(i)  = sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2))
      if(abs(lpabs(i))>1.d-8) then
        lptheta(i)= acos(lzpm(i)/lpabs(i))*rad2deg
      else
        lptheta(i) = 0.d0
      end if
      if(abs(lptheta(i))>1.d-8) then
        if(abs(abs(lptheta(i))-180.d0)>1.d-8) then
          lpphi(i)   = atan2(lypm(i),lxpm(i))*rad2deg
        else
          lpphi(i) = 0.d0
        end if
      else
        lpphi(i) = 0.d0
      end if
    end do

  end subroutine calcLGS


  subroutine calcLGS_greenfunction()
    !! Calculates the expectation value of the orbital angular momentum in the ground state using green functions
    use mod_f90_kind,      only: double
    use mod_constants,     only: cZero,pi
    use mod_System,        only: s => sys
    use mod_parameters,    only: nOrb, nOrb2, eta
    use EnergyIntegration, only: y, wght
    use ElectricField,     only: EshiftBZ,ElectricFieldVector
    use mod_magnet,        only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz 
    use adaptiveMesh
    use mod_mpi_pars
    implicit none
    integer*8    :: ix
    integer      :: AllocateStatus
    integer      :: i,mu,nu,mup,nup
    real(double) :: kp(3)
    real(double) :: weight, ep
    complex(double), dimension(:,:,:,:), allocatable :: gf
    complex(double), dimension(:,:,:),   allocatable :: gupgd
    !--------------------- begin MPI vars --------------------
    integer :: ncount
    ncount=s%nAtoms*nOrb*nOrb
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    allocate(gupgd(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) &
    call abortProgram("[calcLGS_greenfunction] Not enough memory for: gupgd")

    ! Calculating the jacobian using a complex integral
    gupgd  = cZero
    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(local_points,s,nOrb,nOrb2,E_k_imag_mesh,bzs,eta,y,wght,gupgd,EshiftBZ,ElectricFieldVector)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat = AllocateStatus)
    if (AllocateStatus/=0) &
    call abortProgram("[calcLGS_greenfunction] Not enough memory for: gf")

    gf = cZero
    !$omp do schedule(static) reduction(+:gupgd)
    do ix = 1, local_points
        kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix)) + EshiftBZ*ElectricFieldVector
        ep = y(E_k_imag_mesh(1,ix))
        weight = bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix)) * wght(E_k_imag_mesh(1,ix))
        !Green function on energy Ef + iy, and wave vector kp
        call green(s%Ef,ep+eta,s,kp,gf)

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
    !$omp end do

    deallocate(gf)
    !$omp end parallel

    gupgd = gupgd / pi
    call MPI_Allreduce(MPI_IN_PLACE, gupgd, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    lxpm = 0.d0
    lypm = 0.d0
    lzpm = 0.d0
    lxm  = 0.d0
    lym  = 0.d0
    lzm  = 0.d0

    do nu=5,9
      do mu=5,9
        do i=1,s%nAtoms
          lxpm(i) = lxpm(i) + real(lxp(mu,nu,i)*gupgd(nu,mu,i))
          lypm(i) = lypm(i) + real(lyp(mu,nu,i)*gupgd(nu,mu,i))
          lzpm(i) = lzpm(i) + real(lzp(mu,nu,i)*gupgd(nu,mu,i))
          lxm (i) = lxm (i) + real(lx (mu,nu  )*gupgd(nu,mu,i))
          lym (i) = lym (i) + real(ly (mu,nu  )*gupgd(nu,mu,i))
          lzm (i) = lzm (i) + real(lz (mu,nu  )*gupgd(nu,mu,i))
        end do
      end do
    end do

    deallocate(gupgd)
  end subroutine calcLGS_greenfunction


!! Calculate the expectation value of the orbital momentum in the propagated states:
subroutine expec_L_n(s, dim, evec, eval, lxm, lxpm, lym, lypm, lzm, lzpm)
  use mod_f90_kind,      only: double
  use mod_constants,     only: pi
  use mod_parameters,    only: nOrb,eta,isigmamu2n
  use mod_System,        only: system
  use mod_magnet,        only: lxp,lyp,lzp,lx,ly,lz
  use mod_distributions, only: fd_dist
  implicit none 
  type(System),                         intent(in)  :: s
  real(double),                         intent(in)  :: eval
  integer,                              intent(in)  :: dim
  real(double),    dimension(s%nAtoms), intent(out) :: lxm, lxpm, lym, lypm, lzm, lzpm
  complex(double), dimension(dim),      intent(in)  :: evec
  integer                                           :: i, mu, nu, sigma
  real(double)                                      :: f_n
  complex(double)                                   :: prod

  lxm  = 0.d0
  lym  = 0.d0
  lzm  = 0.d0
  lxpm = 0.d0
  lypm = 0.d0
  lzpm = 0.d0

  ! Fermi-Dirac:
  f_n = fd_dist(s%Ef, 1.d0/(pi*eta), eval)

  sites_loop: do i = 1, s%nAtoms
    do sigma = 1, 2
      do nu = 1, nOrb
        do mu = 1, nOrb
          prod = f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*evec(isigmamu2n(i,sigma,nu))
          lxm (i) = lxm (i) + prod*lx (mu,nu  ) !> angular momentum at atomic site (i)
          lym (i) = lym (i) + prod*ly (mu,nu  )
          lzm (i) = lzm (i) + prod*lz (mu,nu  )
          lxpm(i) = lxpm(i) + prod*lxp(mu,nu,i)
          lypm(i) = lypm(i) + prod*lyp(mu,nu,i)
          lzpm(i) = lzpm(i) + prod*lzp(mu,nu,i)
        end do
      end do
    end do
  end do sites_loop

end subroutine expec_L_n


!! Calculate the expectation value of the time-dependent Hamiltonian in the propagated states
subroutine expec_H_n(s, kp, t, dim, evec, eval, E_0)
  use mod_f90_kind,      only: double
  use mod_constants,     only: pi
  use mod_parameters,    only: eta
  use mod_System,        only: system
  use mod_distributions, only: fd_dist
  use mod_imRK4,         only: build_td_hamiltonian
  implicit none

  type(System),                    intent(in)  :: s
  integer,                         intent(in)  :: dim
  complex(double), dimension(dim), intent(in)  :: evec
  real(double),                    intent(in)  :: eval, t, kp(3)
  integer                                      :: i, j
  real(double)                                 :: f_n, expec_H_0, E_0
  complex(double)                              :: hamilt_t(dim,dim), hamilt_0(dim,dim)

  ! Fermi-Dirac:
  f_n = fd_dist(s%Ef, 1.d0/(pi*eta), eval)

  call build_td_hamiltonian(s,t,kp,eval,hamilt_t,hamilt_0)

  E_0 = 0.d0

  do i=0, dim
    do j=0, dim
      ! expec_H_0 = real( conjg( evec(i) ) * hamilt_t(i,j) * evec(j) )
      expec_H_0 = real( conjg( evec(i) ) * hamilt_0(i,j) * evec(j) )
      E_0       =  E_0 + f_n * expec_H_0
    end do 
  end do
  
end subroutine expec_H_n


  !   Calculates ground state quantities from eigenstates
  subroutine calcLGS_eigenstates()
    use mod_f90_kind,      only: double
    use mod_BrillouinZone, only: realBZ
    use mod_constants,     only: pi,cZero
    use mod_parameters,    only: nOrb,nOrb2,output,eta,isigmamu2n
    use mod_System,        only: s => sys
    use ElectricField,     only: EshiftBZ,ElectricFieldVector
    use mod_magnet,        only: lxm,lym,lzm,lxpm,lypm,lzpm,lxp,lyp,lzp,lx,ly,lz
    use mod_distributions, only: fd_dist
    use mod_tools,         only: itos
    use mod_mpi_pars
    implicit none
    integer                                        :: iz, info , n, i, mu, nu, sigma
    integer                                        :: lwork,dimH
    real(double)                                   :: weight, kp(3), f_n
    complex(double), dimension(:,:,:), allocatable :: prod
    real(double),    dimension(:),     allocatable :: rwork(:), eval(:)
    complex(double),                   allocatable :: work(:), hk(:,:), evec(:)

    dimH  = (s%nAtoms)*nOrb2
    lwork = 21*dimH

    allocate( hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),evec(dimH),work(lwork),prod(nOrb,nOrb,s%nAtoms) )

    !$omp parallel default(none) &
    !$omp& firstprivate(lwork) &
    !$omp& private(iz,n,i,sigma,mu,nu,kp,weight,hk,eval,f_n,evec,work,rwork,info) &
    !$omp& shared(s,nOrb,dimH,output,realBZ,eta,isigmamu2n,prod,EshiftBZ,ElectricFieldVector)

    prod = cZero
    !$omp do reduction(+:prod) schedule(static)
    kloop: do iz = 1,realBZ%workload
      kp = realBZ%kp(1:3,iz) + EshiftBZ*ElectricFieldVector
      weight = realBZ%w(iz)
      ! Calculating the hamiltonian for a given k-point
      call hamiltk(s,kp,hk)

      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)
      if(info/=0) &
        call abortProgram("[calcLGS_eigenstates] Problem with diagonalization. info = " // itos(info))

      eval_loop: do n = 1, dimH
        ! Fermi-Dirac:
        f_n = fd_dist(s%Ef, 1.d0/(pi*eta), eval(n))

        ! Getting eigenvector and its transpose conjugate
        evec(:) = hk(:,n)

        sites_loop: do i = 1, s%nAtoms
          do nu = 1, nOrb
            do mu = 1, nOrb
              do sigma = 1, 2
                prod(mu,nu,i) = prod(mu,nu,i) + f_n*conjg( evec(isigmamu2n(i,sigma,mu)) )*evec(isigmamu2n(i,sigma,nu))*weight
              end do    
            end do
          end do
        end do sites_loop

      end do eval_loop
    end do kloop
    !$omp end do
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, prod, nOrb*nOrb*s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, FreqComm(1) , ierr)

    ! Building different components of the orbital angular momentum
    lxm  = 0.d0
    lym  = 0.d0
    lzm  = 0.d0
    lxpm = 0.d0
    lypm = 0.d0
    lzpm = 0.d0
    do i = 1, s%nAtoms
      do nu = 1, nOrb
        do mu = 1, nOrb
          lxm (i) = lxm (i) + prod(mu,nu,i)*lx (mu,nu  )
          lym (i) = lym (i) + prod(mu,nu,i)*ly (mu,nu  )
          lzm (i) = lzm (i) + prod(mu,nu,i)*lz (mu,nu  )
          lxpm(i) = lxpm(i) + prod(mu,nu,i)*lxp(mu,nu,i)
          lypm(i) = lypm(i) + prod(mu,nu,i)*lyp(mu,nu,i)
          lzpm(i) = lzpm(i) + prod(mu,nu,i)*lzp(mu,nu,i)
        end do
      end do
    end do

    deallocate(hk,rwork,eval,work)

  end subroutine calcLGS_eigenstates
  
end module mod_expectation
