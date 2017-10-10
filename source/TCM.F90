module mod_gilbert_damping

contains
subroutine calculate_gilbert_damping()
  use mod_f90_kind, only: double
  use mod_System, only: s => sys
  use mod_BrillouinZone, only: BZ
  use mod_parameters, only: outputunit_loop, strSites, eta, Utype, fieldpart, suffix
  use mod_SOC, only: socpart, SOCc
  use ElectricField, only: strElectricField
  use mod_mpi_pars, only: myrank, myrank_row
  implicit none
  complex(double), dimension(s%nAtoms,s%nAtoms,3,3) :: alpha
  integer :: m,n,i
  character(len=500)  :: varm


  if(myrank_row==0) write(outputunit_loop, *) "[calculate_gilbert_damping] Start..."

  call TCM(alpha, local_SO_torque)

  if(myrank == 0) then
    write(varm,"('./results/',a1,'SOC/',a,'/A/TCM/TCM_SO_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(strSites),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
    open (unit=634893, file=varm, status='replace', form='formatted')
    write(unit=634893, fmt="('# i , m , Re(A_mx), Im(A_mx), Re(A_my), Im(A_my), Re(A_mz), Im(A_mz) ')")
    do i = 1, s%nAtoms
      do m = 1, 3
        write(634893,*) i, m, (real(alpha(i,i,m,n)), aimag(alpha(i,i,m,n)), n = 1, 3)
      end do
    end do
    close(unit=634893)
  end if

  call TCM(alpha, local_xc_torque)

  if(myrank == 0) then
    write(varm,"('./results/',a1,'SOC/',a,'/A/TCM/TCM_XC_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(strSites),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
    open (unit=634893, file=varm, status='replace', form='formatted')
    write(unit=634893, fmt="('# i , m , Re(A_mx), Im(A_mx), Re(A_my), Im(A_my), Re(A_mz), Im(A_mz) ')")
    do i = 1, s%nAtoms
      do m = 1, 3
        write(634893,*) i, m, (real(alpha(i,i,m,n)), aimag(alpha(i,i,m,n)), n = 1,3)
      end do
    end do
    close(unit=634893)
  end if


  return
end subroutine calculate_gilbert_damping

subroutine TCM(alpha, torque_fct)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, pi, cOne, cI
  use mod_System, only: s => sys
  use mod_BrillouinZone, only: BZ

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
  integer :: i,j,k,l,m,n, iz, mu
  real(double), dimension(3) :: kp
  real(double) :: wght
  complex(double), dimension(2*nOrb,2*nOrb,3,s%nAtoms) :: torque
  complex(double), dimension(2*nOrb,2*nOrb,3) :: tsq
  complex(double),dimension(:,:,:,:),allocatable :: gf
  complex(double),dimension(:,:),allocatable :: temp1, temp2, temp3
  complex(double), dimension(:,:,:,:), allocatable :: alpha_loc
  integer :: start, end, work, remainder
  ! Calculate workload for each MPI process
  remainder = mod(BZ%nkpt,numprocs)
  if(myrank < remainder) then
    work = ceiling(dble(BZ%nkpt) / dble(numprocs))
    start = myrank*work + 1
    end = (myrank+1) * work
  else
    work = floor(dble(BZ%nkpt) / dble(numprocs))
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
  !$omp& shared(s,BZ,start, end,Ef,eta,torque,alpha)
  allocate(gf(2*nOrb,2*nOrb,s%nAtoms,s%nAtoms), &
           temp1(2*nOrb,2*nOrb), temp2(2*nOrb,2*nOrb), temp3(2*nOrb,2*nOrb), &
           alpha_loc(s%nAtoms,s%nAtoms,3,3))
  gf = cZero
  alpha_loc = cZero
  !$omp do schedule(static)
  do iz = start, end
    kp = BZ%kp(1:3,iz)
    wght = BZ%w(iz)
    gf  = cZero
    call green(Ef, eta, kp, gf)
    do m = 1, 3
      do n = 1, 3
        do i = 1, s%nAtoms
          do j = 1, s%nAtoms
            temp1 = cZero
            temp2 = cZero
            temp3 = cZero
            ! alpha^{mn}_{ij} = Tr ( Torque^m_i Im(G_ij(Ef)) * Torque^n_j * Im(G_ji(Ef)) )
            temp2 = cI * (gf(:,:,j,i) - transpose(conjg(gf(:,:,i,j))))
            temp3 = torque(:,:,n,j)
            call zgemm('n','n',2*nOrb,2*nOrb,2*nOrb,cOne,temp3,2*nOrb,temp2,2*nOrb,cZero,temp1,2*nOrb) ! Torque^n_j * Im(G_ji(Ef))
            temp2 = cI * (gf(:,:,i,j) - transpose(conjg(gf(:,:,i,j))))
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
  use mod_constants, only: cZero, cOne, cI, levi_civita, sigma => pauli_mat
  use TightBinding, only: nOrb
  use mod_System, only: s => sys
  use mod_parameters, only: U
  use mod_magnet, only: mx,my,mz

  complex(double), dimension(2*nOrb,2*nOrb,3,s%nAtoms), intent(out) :: torque
  complex(double), dimension(nOrb,nOrb) :: ident
  real(double), dimension(3,s%nAtoms) :: mag
  integer :: i,m,n,k

  ident = cZero
  do i = 5, nOrb
    ident(i,i) = cOne
  end do

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
  use mod_constants, only: cZero, levi_civita, sigma => pauli_mat, cI
  use mod_System, only: s => sys
  use mod_magnet, only: Lx, Ly, Lz
  use TightBinding, only: nOrb, nOrb2
  use mod_mpi_pars
  implicit none
  integer :: i,m,n,k, j
  complex(double), dimension(nOrb2,nOrb2,3, s%nAtoms), intent(out) :: torque
  complex(double), dimension(nOrb, nOrb, 3, s%nAtoms) :: L
  complex(double), dimension(nOrb, nOrb) :: test1
  complex(double), dimension(2,2) :: test2

  do i = 1, s%nAtoms
    L(1:nOrb,1:nOrb,1, i) = Lx(1:nOrb,1:nOrb)
    L(1:nOrb,1:nOrb,2, i) = Ly(1:nOrb,1:nOrb)
    L(1:nOrb,1:nOrb,3, i) = Lz(1:nOrb,1:nOrb)
  end do

  torque = cZero

  !! [sigma, H_SO]
  !! t_i^m = Lambda_i * sum_nk eps_mnk L_imunu^n * S_imunu^k
  !! t_m = lambda_i * sum_alpha_beta * sum_nk_gamma,zeta eps_mnk L^n_gamma,zeta * sigma^k_alpha_beta
  do i = 1, s%nAtoms
    do m = 1, 3
      do n = 1, 3
        do k = 1, 3
          torque(     1:nOrb ,     1:nOrb ,m,i) = torque(     1:nOrb ,     1:nOrb ,m,i) + L(1:nOrb,1:nOrb,n,i) * sigma(1,1,k) * levi_civita(m,n,k)
          torque(nOrb+1:nOrb2,     1:nOrb ,m,i) = torque(nOrb+1:nOrb2,     1:nOrb ,m,i) + L(1:nOrb,1:nOrb,n,i) * sigma(2,1,k) * levi_civita(m,n,k)
          torque(     1:nOrb ,nOrb+1:nOrb2,m,i) = torque(     1:nOrb ,nOrb+1:nOrb2,m,i) + L(1:nOrb,1:nOrb,n,i) * sigma(1,2,k) * levi_civita(m,n,k)
          torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) = torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) + L(1:nOrb,1:nOrb,n,i) * sigma(2,2,k) * levi_civita(m,n,k)
        end do
      end do
    end do
    torque(:,:,:,i) = torque(:,:,:,i) * s%Types(s%Basis(i)%Material)%Lambda
  end do
  return
end subroutine local_SO_torque

end module mod_gilbert_damping
