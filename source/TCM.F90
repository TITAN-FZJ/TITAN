subroutine calculate_gilbert_damping()
  use mod_f90_kind, only: double
  use mod_System, only: s => sys
  use mod_mpi_pars, only: myrank
  implicit none
  complex(double), dimension(s%nAtoms,s%nAtoms,3,3) :: alpha
  integer :: m,n,i,j
  call TCM(alpha)
  if(myrank == 0) then
    write(*,*) "# m    n    i    j"
    do m = 1, 3
      do n = 1, 3
        do i = 1, s%nAtoms
          do j = 1, s%nAtoms
            print *, m,n,i,j,alpha(j,i,n,m)
          end do
        end do
      end do
    end do
  end if
  return
end subroutine calculate_gilbert_damping

subroutine TCM(alpha)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, pi, zum
  use mod_System, only: s => sys
  use TightBinding, only: nOrb
  use mod_parameters, only: Ef, eta
  use mod_magnet, only: mabs
  implicit none
  complex(double), dimension(s%nAtoms,s%nAtoms,3,3), intent(out) :: alpha
  !! Contains the Gilbert Damping as matrix in sites and cartesian coordinates (nAtoms,nAtoms,3,3)
  integer :: i,j,m,n, iz, mu
  real(double), dimension(3) :: kp
  real(double) :: wght
  complex(double), dimension(2*nOrb,2*nOrb,3,s%nAtoms) :: torque
  complex(double),dimension(:,:,:,:),allocatable :: gf
  complex(double),dimension(:,:),allocatable :: temp1, temp2, temp3

  call local_SO_torque(torque)

  !! Alpha = g/(m*pi)sum_k wkbz * ( Tr(T^nu Im(G(Ef)) T^mu Im(G(Ef))))
  !! Alpha = g/(m*pi)sum_k wkbz * ( sum_ij (T_ii^nu Im(G_ij(Ef)) T_jj^mu Im(G_ji(Ef)))) maybe
  !! g = 2
  !! m = magnetization amplitude

  ! Calculate full greens function at Ef+i*eta

  alpha = zero

  !$omp parallel default(none) &
  !$omp& private(m,n,i,j,mu,iz,kp,wght,gf,temp1,temp2,temp3) &
  !$omp& shared(s,Ef,eta,torque,alpha)
  allocate(gf(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb), &
           temp1(2*nOrb,2*nOrb), temp2(2*nOrb,2*nOrb), temp3(2*nOrb,2*nOrb))
  gf = zero

  !$omp do schedule(static) reduction(+:alpha)
  do iz = 1, s%nkpt
    kp = s%kbz(:,iz)
    wght = s%wkbz(iz)
    gf = zero
    call green(Ef, eta, kp, gf)

    do m = 1, 3
      do n = 1, 3
        do i = 1, s%nAtoms
          do j = 1, s%nAtoms
            temp1 = zero
            temp2 = zero
            temp3 = zero
            ! alpha^{mn}_{ij} = Tr ( Torque^m_i Im(G_ij(Ef)) * Torque^n_j * Im(G_ji(Ef)) )
            temp2 = aimag(gf(i,j,:,:)) ! Im(G_ji(EF))
            temp3 = torque(:,:,n,j)
            call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,zum,temp3,2*nOrb,temp2,2*nOrb,zero,temp1,2*nOrb) ! Torque^n_j * Im(G_ji(Ef))
            temp2 = aimag(gf(j,i,:,:)) ! Im(G_ji(EF))
            call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,zum,temp2,2*nOrb,temp1,2*nOrb,zero,temp3,2*nOrb) ! Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))
            temp2 = torque(:,:,m,i)
            call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,zum,temp2,2*nOrb,temp3,2*nOrb,zero,temp1,2*nOrb) ! Torque^m_i * Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))

            do mu = 1, 2*nOrb
              alpha(j,i,n,m) = alpha(j,i,n,m) + temp1(mu,mu) * wght
            end do
          end do
        end do
      end do
    end do
  end do
  !$omp end do nowait

  deallocate(gf)
  !$omp end parallel

  do i = 1, s%nAtoms
    do j = 1, s%nAtoms
      alpha(j,i,:,:) = 2.d0 * alpha(j,i,:,:) / (sqrt(mabs(i)*mabs(j))*pi)
    end do
  end do

  return
end subroutine TCM


subroutine local_SO_torque(torque)
  use mod_f90_kind, only: double
  use mod_constants, only: zero, zum, zi, levi_civita
  use mod_System, only: s => sys
  use mod_magnet, only: Lxp, Lyp, Lzp
  use TightBinding, only: nOrb
  use mod_mpi_pars, only: myrank
  implicit none
  integer :: i,m,n,k
  complex(double), dimension(2*nOrb,2*nOrb,3, s%nAtoms), intent(out) :: torque
  complex(double), dimension(nOrb, nOrb, 3, s%nAtoms) :: L
  complex(double), dimension(2,2,3) :: sigma

  do i = 1, s%nAtoms
    L(:,:,1, i) = Lxp(:,:,i)
    L(:,:,2, i) = Lyp(:,:,i)
    L(:,:,3, i) = Lzp(:,:,i)
  end do

  sigma(:,1,1) = [zero,zum]
  sigma(:,2,1) = [zum, zero]
  sigma(:,1,2) = [zero,-zi]
  sigma(:,2,2) = [zi,zero]
  sigma(:,1,3) = [zum,zero]
  sigma(:,2,3) = [zero,-zum]

  torque = zero

  !! [sigma, H_SO]
  !! Lambda * i * ( L_i^y sigma^z *e_x - L_i^x*sigma^z*e_y)

  !! t_m = i*lambda_i * sum_alpha_beta * sum_nk_gamma,zeta eps_mnk L^n_gamma,zeta * sigma^k_alpha_beta
  do i = 1, s%nAtoms
    do m = 1, 3
      do n = 1, 3
        do k = 1, 3
          torque(     1:  nOrb,     1:  nOrb,m,i) = torque(     1:  nOrb,     1:  nOrb,m,i) + L(:,:,n,i) * sigma(1,1,k) * levi_civita(m,n,k)
          torque(nOrb+1:2*nOrb,     1:  nOrb,m,i) = torque(nOrb+1:2*nOrb,     1:  nOrb,m,i) + L(:,:,n,i) * sigma(2,1,k) * levi_civita(m,n,k)
          torque(     1:  nOrb,nOrb+1:2*nOrb,m,i) = torque(     1:  nOrb,nOrb+1:2*nOrb,m,i) + L(:,:,n,i) * sigma(1,2,k) * levi_civita(m,n,k)
          torque(nOrb+1:2*nOrb,nOrb+1:2*nOrb,m,i) = torque(nOrb+1:2*nOrb,nOrb+1:2*nOrb,m,i) + L(:,:,n,i) * sigma(2,2,k) * levi_civita(m,n,k)
        end do
      end do
    end do
    torque(:,:,:,i) = torque(:,:,:,i) * zi * s%Types(s%Basis(i)%Material)%Lambda
  end do

  if (myrank == 0) then
    do i=1,2*nOrb
      print *, (torque(n,i,3,1), n = 1, 2*nOrb)
    end do
  end if
  return
end subroutine local_SO_torque
