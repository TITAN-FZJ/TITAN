module mod_gilbert_damping

contains
subroutine calculate_gilbert_damping()
  use mod_f90_kind, only: double
  use mod_System, only: s => sys
  use mod_parameters, only: outputunit
  use mod_mpi_pars, only: myrank
  implicit none
  complex(double), dimension(s%nAtoms,s%nAtoms,3,3) :: alpha
  integer :: m,n,i,j

  write(outputunit, *) "[calculate_gilbert_damping] Start..."

  if(myrank == 0) print *, "SO Torque"
  call TCM(alpha, local_SO_torque)

  if(myrank == 0) then
    write(*,*) "# m          i"
    do m = 1, 3
      do i = 1, s%nAtoms
        print *, m,i,alpha(i,i,m,m)
      end do
    end do
  end if

  if(myrank == 0) print *, "xc Torque"
  call TCM(alpha, local_xc_torque)

  if(myrank == 0) then
    write(*,*) "# m          i"
    do m = 1, 3
      do i = 1, s%nAtoms
        print *, m,i,alpha(i,i,m,m)
      end do
    end do
  end if


  return
end subroutine calculate_gilbert_damping

subroutine TCM(alpha, torque_fct)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, pi, cOne
  use mod_System, only: s => sys
  use TightBinding, only: nOrb
  use mod_parameters, only: Ef, eta
  use mod_magnet, only: mabs
  use mod_mpi_pars
  implicit none
  interface
    subroutine torque_fct(torque)
      use mod_f90_kind, only: double
      use TightBinding, only: nOrb
      use mod_System, only: sys
      implicit none
      complex(double), dimension(2*nOrb,2*nOrb,3,sys%nAtoms), intent(out) :: torque
    end subroutine torque_fct
  end interface

  complex(double), dimension(s%nAtoms,s%nAtoms,3,3), intent(out) :: alpha
  !! Contains the Gilbert Damping as matrix in sites and cartesian coordinates (nAtoms,nAtoms,3,3)
  integer :: i,j,m,n, iz, mu
  real(double), dimension(3) :: kp
  real(double) :: wght
  complex(double), dimension(2*nOrb,2*nOrb,3,s%nAtoms) :: torque
  complex(double),dimension(:,:,:,:),allocatable :: gf
  complex(double),dimension(:,:),allocatable :: temp1, temp2, temp3
  complex(double), dimension(:,:,:,:), allocatable :: alpha_loc
  integer :: start, end, work, remainder

  ! Calculate workload for each MPI process
  remainder = mod(s%nkpt,numprocs)
  if(myrank < remainder) then
    work = ceiling(dble(s%nkpt) / dble(numprocs))
    start = myrank*work + 1
    end = (myrank+1) * work
  else
    work = floor(dble(s%nkpt) / dble(numprocs))
    start = myrank*work + 1 + remainder
    end = (myrank+1) * work + remainder
  end if


  call torque_fct(torque)

  !! Alpha = g/(m*pi)sum_k wkbz * ( Tr(T^nu Im(G(Ef)) T^mu Im(G(Ef))))
  !! Alpha = g/(m*pi)sum_k wkbz * ( sum_ij (T_ii^nu Im(G_ij(Ef)) T_jj^mu Im(G_ji(Ef)))) maybe
  !! g = 2
  !! m = magnetization amplitude

  ! Calculate full greens function at Ef+i*eta

  alpha = cZero

  !$omp parallel default(none) &
  !$omp& private(m,n,i,j,mu,iz,kp,wght,gf,temp1,temp2,temp3, alpha_loc) &
  !$omp& shared(s,start, end,Ef,eta,torque,alpha)
  allocate(gf(2*nOrb,2*nOrb,s%nAtoms,s%nAtoms), &
           temp1(2*nOrb,2*nOrb), temp2(2*nOrb,2*nOrb), temp3(2*nOrb,2*nOrb), &
           alpha_loc(s%nAtoms,s%nAtoms,3,3))
  gf = cZero
  alpha_loc = cZero
  !$omp do schedule(static)
  do iz = start, end
    kp = s%kbz(:,iz)
    wght = s%wkbz(iz)
    gf = cZero
    call green(Ef, eta, kp, gf)

    do m = 1, 3
      do n = 1, 3
        do i = 1, s%nAtoms
          do j = 1, s%nAtoms
            temp1 = cZero
            temp2 = cZero
            temp3 = cZero
            ! alpha^{mn}_{ij} = Tr ( Torque^m_i Im(G_ij(Ef)) * Torque^n_j * Im(G_ji(Ef)) )
            temp2 = aimag(gf(:,:,j,i)) ! Im(G_ji(EF))
            temp3 = torque(:,:,n,j)
            call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,cOne,temp3,2*nOrb,temp2,2*nOrb,cZero,temp1,2*nOrb) ! Torque^n_j * Im(G_ji(Ef))
            temp2 = aimag(gf(:,:,i,j)) ! Im(G_ij(EF))
            call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,cOne,temp2,2*nOrb,temp1,2*nOrb,cZero,temp3,2*nOrb) ! Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))
            temp2 = torque(:,:,m,i)
            call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,cOne,temp2,2*nOrb,temp3,2*nOrb,cZero,temp1,2*nOrb) ! Torque^m_i * Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))

            do mu = 1, 2*nOrb
              alpha_loc(j,i,n,m) = alpha_loc(j,i,n,m) + temp1(mu,mu) * wght
            end do
          end do
        end do
      end do
    end do
  end do
  !$omp end do nowait

  !$omp critical
    alpha = alpha + alpha_loc
  !$omp end critical

  deallocate(gf)
  deallocate(alpha_loc)
  !$omp end parallel

  do i = 1, s%nAtoms
    do j = 1, s%nAtoms
      alpha(j,i,:,:) = 2.d0 * alpha(j,i,:,:) / (sqrt(mabs(i)*mabs(j))*pi)
    end do
  end do

  call MPI_Allreduce(MPI_IN_PLACE, alpha, s%nAtoms*s%nAtoms*3*3, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

  return
end subroutine TCM


subroutine local_xc_torque(torque)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cOne, cI, levi_civita
  use TightBinding, only: nOrb
  use mod_System, only: s => sys
  use mod_parameters, only: U
  use mod_magnet, only: mx,my,mz

  complex(double), dimension(2*nOrb,2*nOrb,3,s%nAtoms), intent(out) :: torque
  complex(double), dimension(2,2,3) :: sigma
  complex(double), dimension(nOrb,nOrb) :: ident
  real(double), dimension(3,s%nAtoms) :: mag
  integer :: i,mu,m,n,k

  ident = cZero
  do i = 5, nOrb
    ident(i,i) = cOne
  end do

  sigma(:,1,1) = [cZero,cOne]
  sigma(:,2,1) = [cOne, cZero]
  sigma(:,1,2) = [cZero,-cI]
  sigma(:,2,2) = [cI,cZero]
  sigma(:,1,3) = [cOne,cZero]
  sigma(:,2,3) = [cZero,-cOne]

  do i = 1, s%nAtoms
    mag(:,i) = [mx(i),my(i),mz(i)]
  end do

  torque = cZero
  do i = 1, s%nAtoms
    do m = 1, 3
      do n = 1, 3
        do k = 1, 3
          torque(     1:  nOrb,     1:  nOrb,m,i) = torque(     1:  nOrb,     1:  nOrb,m,i) + mag(n,i) * ident(:,:) * sigma(1,1,k) * levi_civita(m,n,k)
          torque(nOrb+1:2*nOrb,     1:  nOrb,m,i) = torque(nOrb+1:2*nOrb,     1:  nOrb,m,i) + mag(n,i) * ident(:,:) * sigma(2,1,k) * levi_civita(m,n,k)
          torque(     1:  nOrb,nOrb+1:2*nOrb,m,i) = torque(     1:  nOrb,nOrb+1:2*nOrb,m,i) + mag(n,i) * ident(:,:) * sigma(1,2,k) * levi_civita(m,n,k)
          torque(nOrb+1:2*nOrb,nOrb+1:2*nOrb,m,i) = torque(nOrb+1:2*nOrb,nOrb+1:2*nOrb,m,i) + mag(n,i) * ident(:,:) * sigma(2,2,k) * levi_civita(m,n,k)
        end do
      end do
    end do
    torque(:,:,:,i) = - U(i) * torque(:,:,:,i)
  end do

  return
end subroutine local_xc_torque

subroutine local_SO_torque(torque)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cOne, cI, levi_civita
  use mod_System, only: s => sys
  use mod_magnet, only: Lxp, Lyp, Lzp
  use TightBinding, only: nOrb
  use mod_mpi_pars, only: myrank
  use mod_SOC, only: SOC
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

  sigma(:,1,1) = [cZero,cOne]
  sigma(:,2,1) = [cOne, cZero]
  sigma(:,1,2) = [cZero,-cI]
  sigma(:,2,2) = [cI,cZero]
  sigma(:,1,3) = [cOne,cZero]
  sigma(:,2,3) = [cZero,-cOne]

  torque = cZero

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
    torque(:,:,:,i) = torque(:,:,:,i) * cI * s%Types(s%Basis(i)%Material)%Lambda
  end do
  return
end subroutine local_SO_torque

end module mod_gilbert_damping
