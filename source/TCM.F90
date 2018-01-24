module mod_gilbert_damping
  use mod_f90_kind, only: double
  implicit none
  character(len=5), private :: folder = "A/TCM"
  character(len=2), dimension(2), private :: filename = ["SO", "XC"]

contains

  subroutine openTCMFiles()
    use mod_parameters, only: output
    implicit none
    integer :: n,iw
    character(len=500)  :: varm

    do n = 1, size(filename)
      iw = 634893 + n

      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(n)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
      open (unit=iw, file=varm, status='replace', form='formatted')
      write(unit=iw, fmt="('# i , m , Re(A_mx), Im(A_mx), Re(A_my), Im(A_my), Re(A_mz), Im(A_mz) ')")

    end do
    return
  end subroutine openTCMFiles

  subroutine closeTCMFiles()
    implicit none
    integer :: i, iw
    do i = 1, size(filename)
      iw = 634893 + i
      close(iw)
    end do
    return
  end subroutine closeTCMFiles

  subroutine writeTCM()

  end subroutine writeTCM

subroutine calculate_gilbert_damping()
  use mod_f90_kind, only: double
  use mod_parameters, only: output
  use mod_System, only: s => sys
  use mod_mpi_pars, only: rField
  implicit none
  complex(double), dimension(s%nAtoms,s%nAtoms,3,3) :: alphaSO, alphaXC
  integer :: i,j,k

  if(rField == 0) write(output%unit_loop, *) "[calculate_gilbert_damping] Start..."

  call TCM(alphaSO, local_SO_torque)
  call TCM(alphaXC, local_xc_torque)

  if(rField == 0) then
    call openTCMFiles()
  end if

  if(rField == 0) then
    do i = 1, s%nAtoms
      do j = 1, 3
        write(634894,*) i, j, (real(alphaSO(i,i,j,k)), aimag(alphaSO(i,i,j,k)), k = 1, 3)
        write(634895,*) i, j, (real(alphaXC(i,i,j,k)), aimag(alphaXC(i,i,j,k)), k = 1, 3)
      end do
    end do
  end if

  if(rField == 0) call closeTCMFiles()

  return
end subroutine calculate_gilbert_damping

subroutine TCM(alpha, torque_fct)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, pi, cOne, cI
  use mod_System, only: s => sys
  use mod_BrillouinZone, only: realBZ

  use TightBinding, only: nOrb
  use mod_parameters, only: eta
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
  integer :: i,j,m,n, mu
  real(double), dimension(3) :: kp
  real(double) :: wght
  complex(double), dimension(2*nOrb,2*nOrb,3,s%nAtoms) :: torque
  complex(double),dimension(:,:,:,:),allocatable :: gf
  complex(double),dimension(:,:),allocatable :: temp1, temp2, temp3
  complex(double), dimension(:,:,:,:), allocatable :: alpha_loc
  integer*8 :: firstPoint, lastPoint, iz
  ! Calculate workload for each MPI process

  call realBZ % setup_fraction(rField, sField, FieldComm)

  call torque_fct(torque)

  !! Alpha = g/(m*pi)sum_k wkbz * ( Tr(T^nu Im(G(Ef)) T^mu Im(G(Ef))))
  !! Alpha = g/(m*pi)sum_k wkbz * ( sum_ij (T_ii^nu Im(G_ij(Ef)) T_jj^mu Im(G_ji(Ef)))) maybe
  !! g = 2
  !! m = magnetization amplitude

  ! Calculate full greens function at Ef+i*eta

  alpha = cZero

  !$omp parallel default(none) &
  !$omp& private(m,n,i,j,mu,iz,kp,wght,gf,temp1,temp2,temp3, alpha_loc) &
  !$omp& shared(s,realBZ,firstPoint,lastPoint,eta,torque,alpha)
  allocate(gf(2*nOrb,2*nOrb,s%nAtoms,s%nAtoms), &
           temp1(2*nOrb,2*nOrb), temp2(2*nOrb,2*nOrb), temp3(2*nOrb,2*nOrb), &
           alpha_loc(s%nAtoms,s%nAtoms,3,3))
  gf = cZero
  alpha_loc = cZero
  !$omp do schedule(static)
  do iz = 1, realBZ%workload
    kp = realBZ%kp(1:3,iz)
    wght = realBZ%w(iz)
    gf  = cZero
    call green(s%Ef, eta, kp, gf)
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

  call MPI_Allreduce(MPI_IN_PLACE, alpha, s%nAtoms*s%nAtoms*3*3, MPI_DOUBLE_COMPLEX, MPI_SUM, FieldComm, ierr)

  return
end subroutine TCM


subroutine local_xc_torque(torque)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cOne, cI, levi_civita, sigma => pauli_mat
  use TightBinding, only: nOrb
  use mod_System, only: s => sys
  use mod_parameters, only: U
  use mod_magnet, only: mvec_cartesian

  complex(double), dimension(2*nOrb,2*nOrb,3,s%nAtoms), intent(out) :: torque
  complex(double), dimension(nOrb,nOrb) :: ident
  integer :: i,m,n,k

  ident = cZero
  do i = 5, nOrb
    ident(i,i) = cOne
  end do

  torque = cZero
  do i = 1, s%nAtoms
    do m = 1, 3
      do n = 1, 3
        do k = 1, 3
          torque(     1:  nOrb,     1:  nOrb,m,i) = torque(     1:  nOrb,     1:  nOrb,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(1,1,k) * levi_civita(m,n,k)
          torque(nOrb+1:2*nOrb,     1:  nOrb,m,i) = torque(nOrb+1:2*nOrb,     1:  nOrb,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(2,1,k) * levi_civita(m,n,k)
          torque(     1:  nOrb,nOrb+1:2*nOrb,m,i) = torque(     1:  nOrb,nOrb+1:2*nOrb,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(1,2,k) * levi_civita(m,n,k)
          torque(nOrb+1:2*nOrb,nOrb+1:2*nOrb,m,i) = torque(nOrb+1:2*nOrb,nOrb+1:2*nOrb,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(2,2,k) * levi_civita(m,n,k)
        end do
      end do
    end do
    torque(:,:,:,i) = - 0.5d0 * U(i) * torque(:,:,:,i)
  end do

  return
end subroutine local_xc_torque

subroutine local_SO_torque(torque)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, levi_civita, sigma => pauli_mat, cI
  use mod_System, only: s => sys
  use mod_magnet, only: l
  use TightBinding, only: nOrb, nOrb2
  use mod_mpi_pars
  implicit none
  integer :: i,m,n,k
  complex(double), dimension(nOrb2,nOrb2,3, s%nAtoms), intent(out) :: torque

  torque = cZero

  !! [sigma, H_SO]
  !! t_i^m = Lambda_i * sum_nk eps_mnk L_imunu^n * S_imunu^k
  !! t_m = lambda_i * sum_alpha_beta * sum_nk_gamma,zeta eps_mnk L^n_gamma,zeta * sigma^k_alpha_beta
  do i = 1, s%nAtoms
    do m = 1, 3
      do n = 1, 3
        do k = 1, 3
          torque(     1:nOrb ,     1:nOrb ,m,i) = torque(     1:nOrb ,     1:nOrb ,m,i) + l(1:nOrb,1:nOrb,n) * sigma(1,1,k) * levi_civita(m,n,k)
          torque(nOrb+1:nOrb2,     1:nOrb ,m,i) = torque(nOrb+1:nOrb2,     1:nOrb ,m,i) + l(1:nOrb,1:nOrb,n) * sigma(2,1,k) * levi_civita(m,n,k)
          torque(     1:nOrb ,nOrb+1:nOrb2,m,i) = torque(     1:nOrb ,nOrb+1:nOrb2,m,i) + l(1:nOrb,1:nOrb,n) * sigma(1,2,k) * levi_civita(m,n,k)
          torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) = torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) + l(1:nOrb,1:nOrb,n) * sigma(2,2,k) * levi_civita(m,n,k)
        end do
      end do
    end do
    torque(:,:,:,i) = 0.25d0 * torque(:,:,:,i) * s%Types(s%Basis(i)%Material)%Lambda
  end do
  return
end subroutine local_SO_torque

end module mod_gilbert_damping
