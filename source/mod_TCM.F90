!> This subroutine includes the calculation of the Gilbert damping
!> and the related subroutines.
module mod_TCM
  use mod_kind, only: dp
  implicit none
  character(len=5),               private :: folder = "A/TCM"
  character(len=3), dimension(4), private :: filename = ["SO ", "SOi", "XC ", "XCi"]

contains

  !> Open files for the results obtained with the torque-correlation method
  subroutine openTCMFiles()
    use mod_parameters, only: output
    use mod_System,     only: s => sys
    implicit none
    integer :: i,iw
    character(len=500)  :: varm

    do i = 1, size(filename)
      iw = 634893 + i
      ! Skip files for integrated alpha if not in 3D
      if((s%isysdim/=3).and.(scan(filename(i),"i")>0)) cycle
      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(i)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=iw, file=varm, status='replace', form='formatted')
      write(unit=iw, fmt="('# (k), i , m , Re(A_mx), Im(A_mx), Re(A_my), Im(A_my), Re(A_mz), Im(A_mz) ')")
    end do

  end subroutine openTCMFiles


  !> Write the results for the damping calculated with the TCM
  subroutine writeTCM(unit,alpha,ndiffk,diff_k,ialpha,iwght)
    use mod_System,     only: s => sys
    implicit none
    integer, intent(in) :: unit
    !! File unit
    integer, intent(in) :: ndiffk
    !! Number of different kzs
    real(dp), dimension(ndiffk), intent(in) :: diff_k
    !! Different kzs
    real(dp), dimension(ndiffk), intent(in) :: iwght
    !! Weights of different kzs
    real(dp), dimension(s%nAtoms,3,ndiffk), intent(in) :: ialpha
    !! Contributions of each kz to alpha
    integer :: i,j,k
    complex(dp), dimension(s%nAtoms,s%nAtoms,3,3),intent(in) :: alpha
    !! Total Gilbert damping matrix alpha

    !! Writing Gilbert damping
    do i = 1, s%nAtoms
      do j = 1, 3
        write(unit=unit  ,fmt=*) i, j, (real(alpha(i,i,j,k)), aimag(alpha(i,i,j,k)), k = 1, 3)
      end do
    end do

    !! Writing k-dependent integrated Gilbert damping alpha (3D only)
    if(s%isysdim==3) then
      do k=1,size(diff_k)
        write(unit=unit+1,fmt=*) diff_k(k), ( (ialpha(i,j,k), j=1,3), i=1,s%nAtoms ), iwght(k)
      end do
    end if
  end subroutine writeTCM

  !> Close files for the results obtained with the torque-correlation method
  subroutine closeTCMFiles()
    use mod_System,     only: s => sys
    implicit none
    integer :: i, iw

    do i = 1, size(filename)
      iw = 634893 + i
      ! Skip files for integrated alpha if not in 3D
      if((s%isysdim/=3).and.(scan(filename(i),"i")>0)) cycle
      close(iw)
    end do
  end subroutine closeTCMFiles

  !> Driver routine to call the different torque-correlation methods and their I/O
  subroutine calculate_TCM()
    use mod_kind, only: dp
    use mod_parameters, only: output
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: rField
    use mod_progress,   only: write_time
    implicit none
    integer :: ndiffk
    real(dp),    dimension(:),     allocatable    :: diff_k, iwght
    real(dp),    dimension(:,:,:), allocatable    :: ialphaSO, ialphaXC
    complex(dp), dimension(s%nAtoms,s%nAtoms,3,3) :: alphaSO, alphaXC

    if(rField == 0) &
      write(output%unit_loop, "('[calculate_gilbert_damping] Starting to calculate Gilbert damping within the TCM:')")

    if(rField == 0) &
      call openTCMFiles()

    if(rField == 0) &
      write(output%unit_loop, "('[calculate_gilbert_damping] SO-TCM...')")

    ! SO-TCM
    call TCM(local_SO_torque, alphaSO, ndiffk, diff_k, ialphaSO, iwght)

    if(rField == 0) then
      write(output%unit_loop, "('[calculate_gilbert_damping] Writing results of SO-TCM to file...')")
      call writeTCM(634894, alphaSO, ndiffk, diff_k, ialphaSO, iwght)
      call write_time(output%unit_loop,'[main] Time after SO-TCM: ')
      write(output%unit_loop, "('[calculate_gilbert_damping] XC-TCM...')")
    end if

    ! XC-TCM
    call TCM(local_xc_torque, alphaXC, ndiffk, diff_k, ialphaXC, iwght)

    if(rField == 0) then
      write(output%unit_loop, "('[calculate_gilbert_damping] Writing results of XC-TCM to file...')")
      call writeTCM(634896, alphaXC, ndiffk, diff_k, ialphaXC, iwght)
      call write_time(output%unit_loop,'[main] Time after XC-TCM: ')
      call closeTCMFiles()
    end if

  end subroutine calculate_TCM


  !> This subroutine calculates the Gilbert damping using the torque-correlation method
  subroutine TCM(torque_fct, alpha, ndiffk, diff_k, ialpha, iwght)
    use mod_kind, only: dp
    use mod_constants,     only: cZero, pi, cOne, cI
    use mod_System,        only: s => sys
    use mod_BrillouinZone, only: realBZ, store_diff
    use mod_parameters,    only: nOrb2, eta, output, kptotal_in
    use mod_magnet,        only: mabs
    use mod_tools,         only: sort, itos, rtos
    use mod_hamiltonian,   only: hamilt_local,h0
    use mod_mpi_pars
    implicit none
    interface
      subroutine torque_fct(torque)
        use mod_kind, only: dp
        use mod_parameters, only: nOrb2
        use mod_System,     only: sys
        implicit none
        complex(dp), dimension(nOrb2,nOrb2,3,sys%nAtoms), intent(out) :: torque
      end subroutine torque_fct
    end interface

    complex(dp), dimension(s%nAtoms,s%nAtoms,3,3), intent(out) :: alpha
    !! Contains the Gilbert Damping as matrix in sites and cartesian coordinates (nAtoms,nAtoms,3,3)
    integer(int64) :: iz
    integer   :: i, j, m, n, k, mu, ios

    integer   :: ndiffk
    !! Number of different abs(k_z) in the k points list
    integer(int64) :: nk
    !! Total number of k points
    real(dp), dimension(:), allocatable :: diff_k_unsrt
    !! Different values of abs(k_z) - unsorted
    real(dp), dimension(:), allocatable :: diff_k
    !! Different values of abs(k_z) - sorted
    integer,      dimension(:), allocatable :: order
    !! Order of increasing abs(k_z) in the different k points list

    real(dp), dimension(3) :: kp
    !! k-point
    real(dp) :: wght
    !! k-point weight
    real(dp), dimension(:), allocatable :: iwght
    !! integrated weight (up to k_z^max)
    real(dp), dimension(:,:,:), allocatable, intent(out) :: ialpha
    !! Integrated alpha (as a function of the maximum kz summed) $ \alpha = \sum_{k_z}^{k_z^\text{max}} \alpha(k_z)$
    complex(dp), dimension(nOrb2,nOrb2,3,s%nAtoms) :: torque
    !! Torque operator (SO or XC)
    complex(dp),dimension(:,:,:,:),allocatable :: gf
    !! Green function
    complex(dp),dimension(:,:),allocatable :: temp1, temp2, temp3
    complex(dp) :: alphatemp
    !! Temporary k-dependent alpha
    complex(dp), dimension(:,:,:,:), allocatable :: alpha_loc
    !! Local (in-processor) alpha


    ! Calculate workload for each MPI process
    if(rField == 0) &
      write(output%unit_loop, "('[TCM] Generating local k-points...')")
    call realBZ % setup_fraction(s,rField, sField, FieldComm)

    ! Calculating integrated k-dependent alpha (3D only)
    if(s%isysdim==3) then
      ! Obtaining the different kz
      if(rField == 0) &
        write(output%unit_loop, "('[TCM] Obtaining different kz: ')", advance='no')
      open (unit=99999, file='diffkz_'//trim(itos(kptotal_in)),status='old', iostat = ios)
      if(ios /= 0) then
        if(rField == 0) &
          write(output%unit_loop, "('Generating...')")
        nk = int(realBZ%nkpt_x*realBZ%nkpt_y*realBZ%nkpt_z,kind(nk))
        call store_diff(nk, s%b1, s%b2, s%b3, 3, ndiffk, diff_k_unsrt)
      else
        if(rField == 0) &
          write(output%unit_loop, "('Reading from file...')")
        read(unit=99999, fmt=*) ndiffk
        allocate( diff_k_unsrt(ndiffk) )
        do i=1,ndiffk
          read(unit=99999, fmt=*) diff_k_unsrt(i)
        end do
        close(unit=99999)
      end if

      ! Sorting the k-points for increasing abs(kz)
      if(rField == 0) &
        write(output%unit_loop, "('[TCM] Sorting different kz...')")
      allocate( order(ndiffk) )
      if(.not.allocated(diff_k)) allocate( diff_k(ndiffk) )
      if(.not.allocated(iwght))  allocate( iwght(ndiffk) )
      call sort( diff_k_unsrt(:), ndiffk, order(:) )
      sort_diffk: do k = 1, ndiffk
        diff_k(k) = diff_k_unsrt(order(k))
      end do sort_diffk

      iwght = 0._dp

      !! k-dependent alpha(k_z)
      !! $\alpha(k_z) = \frac{1}{N}\sum_{k_x,k_y}\alpha(k_x,k_y,k_z)$

      !! Allocating integrated alpha
      allocate( ialpha(s%nAtoms,3,ndiffk) )
      ialpha = 0._dp
    end if

    call torque_fct(torque)

    !! XC-TCM:
    !! $\alpha= \frac{1}{\gamma M \pi N} \sum_{\mathbf{k}} \operatorname{Tr}\{\operatorname{Im} G(\mathbf{k},\epsilon_\text{F})\hat{T}^-_\text{xc} \operatorname{Im} G(\mathbf{k},\epsilon_F) \hat{T}^+_\text{xc}\}$
    !! SO-TCM:
    !! $\alpha - \alpha_{\text{noSOI}} = \frac{1}{\gamma M \pi N} \sum_{\mathbf{k}}\operatorname{Tr}\{\operatorname{Im} G(\mathbf{k},\epsilon_F)\hat{T}^-_{\text{SOI}} \operatorname{Im} G(\mathbf{k},\epsilon_F) \hat{T}^+_{\text{SOI}}\} $
    !! g = 2
    !! M = magnetic moment
    alpha(:,:,:,:) = cZero

    ! Build local hamiltonian
    call hamilt_local(s)

    if(rField == 0) &
      write(output%unit_loop, "('[TCM] Calculating alpha...')")

    !$omp parallel default(none) reduction(+:ialpha,iwght) &
    !$omp& private(m,n,i,j,k,mu,iz,kp,wght,gf,temp1,temp2,temp3, alpha_loc,alphatemp) &
    !$omp& shared(s,realBZ,nOrb2,eta,torque,alpha,rField,order,ndiffk,diff_k)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), &
             temp1(nOrb2,nOrb2), temp2(nOrb2,nOrb2), temp3(nOrb2,nOrb2), &
             alpha_loc(s%nAtoms,s%nAtoms,3,3))
    gf = cZero
    alpha_loc  = cZero

    !$omp do schedule(static)
    kpoints: do iz = 1, realBZ%workload
      kp   = realBZ%kp(1:3,iz)
      wght = realBZ%w(iz)
      gf   = cZero
      call green(s%Ef,eta,s,kp,gf)
      dir_m: do m = 1, 3
        dir_n: do n = 1, 3
          site_i: do i = 1, s%nAtoms
            site_j: do j = 1, s%nAtoms
              alphatemp = cZero
              temp1 = cZero
              temp2 = cZero
              temp3 = cZero
              ! alpha^{mn}_{ij} = Tr ( Torque^m_i Im(G_ij(Ef)) * Torque^n_j * Im(G_ji(Ef)) )
              temp2 = cI * (gf(:,:,j,i) - transpose(conjg(gf(:,:,i,j))))
              temp3 = torque(:,:,n,j)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp3,nOrb2,temp2,nOrb2,cZero,temp1,nOrb2) ! Torque^n_j * Im(G_ji(Ef))
              temp2 = cI * (gf(:,:,i,j) - transpose(conjg(gf(:,:,i,j))))
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp2,nOrb2,temp1,nOrb2,cZero,temp3,nOrb2) ! Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))
              temp2 = torque(:,:,m,i)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp2,nOrb2,temp3,nOrb2,cZero,temp1,nOrb2) ! Torque^m_i * Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))

              do mu = 1, nOrb2
                alphatemp = alphatemp + temp1(mu,mu) * wght
              end do
              alpha_loc(j,i,n,m) = alpha_loc(j,i,n,m) + alphatemp

              if((s%isysdim==3).and.(i == j).and.(m == n)) then
                ! Testing if kz is not on the list of different kz calculated before
                diffk: do k = 1, ndiffk
                  if ( abs(abs(kp(3)) - diff_k(k)) < 1.e-15_dp ) then
                    ialpha(i,m,k) = ialpha(i,m,k) + real(alphatemp)
                    if( (i==1) .and. (m==1) ) iwght(k) = iwght(k) + wght
                    exit diffk
                  end if
                end do diffk
              end if

            end do site_j
          end do site_i
        end do dir_n
      end do dir_m
    end do kpoints
    !$omp end do nowait

    !$omp critical
      alpha = alpha + alpha_loc
    !$omp end critical

    deallocate(gf,temp1,temp2,temp3)
    deallocate(alpha_loc)
    !$omp end parallel

    do i = 1, s%nAtoms
      do j = 1, s%nAtoms
        alpha(j,i,:,:) = 2._dp * alpha(j,i,:,:) / (sqrt(mabs(i)*mabs(j))*pi)
      end do
      if(s%isysdim==3) ialpha(i,:,:) = 2._dp * ialpha(i,:,:) / (mabs(i)*pi)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, alpha , s%nAtoms*s%nAtoms*3*3, MPI_DOUBLE_COMPLEX  , MPI_SUM, FieldComm, ierr)
    if(s%isysdim==3) then
      call MPI_Allreduce(MPI_IN_PLACE, ialpha, s%nAtoms*3*ndiffk    , MPI_DOUBLE_PRECISION, MPI_SUM, FieldComm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, iwght , ndiffk               , MPI_DOUBLE_PRECISION, MPI_SUM, FieldComm, ierr)
      deallocate(order,diff_k_unsrt)
    end if

    ! Deallocate local hamiltonian
    deallocate(h0)

  end subroutine TCM


  !> This subroutine defines the exchange-correlation torque operator/matrix
  subroutine local_xc_torque(torque)
    use mod_kind, only: dp
    use mod_constants,  only: cZero, cOne, levi_civita, sigma => pauli_mat
    use mod_System,     only: s => sys
    use mod_parameters, only: nOrb, nOrb2, Um
    use mod_magnet,     only: mvec_cartesian

    complex(dp), dimension(nOrb2,nOrb2,3,s%nAtoms), intent(out) :: torque
    complex(dp), dimension(nOrb,nOrb) :: ident
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
            torque(     1: nOrb,     1: nOrb,m,i) = torque(     1: nOrb,     1: nOrb,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(1,1,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,     1: nOrb,m,i) = torque(nOrb+1:nOrb2,     1: nOrb,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(2,1,k) * levi_civita(m,n,k)
            torque(     1: nOrb,nOrb+1:nOrb2,m,i) = torque(     1: nOrb,nOrb+1:nOrb2,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(1,2,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) = torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) + mvec_cartesian(n,i) * ident(:,:) * sigma(2,2,k) * levi_civita(m,n,k)
          end do
        end do
      end do
      torque(:,:,:,i) = - 0.5_dp * Um(i) * torque(:,:,:,i)
    end do
  end subroutine local_xc_torque


  !> This subroutine defines the spin-orbit torque operator/matrix
  subroutine local_SO_torque(torque)
    use mod_kind, only: dp
    use mod_constants,  only: cZero, levi_civita, sigma => pauli_mat
    use mod_System,     only: s => sys
    use mod_magnet,     only: lvec
    use mod_parameters, only: nOrb, nOrb2
    use mod_mpi_pars
    implicit none
    integer :: i,m,n,k
    complex(dp), dimension(nOrb2,nOrb2,3, s%nAtoms), intent(out) :: torque

    torque = cZero

    !! [sigma, H_SO]
    !! t_i^m = Lambda_i * sum_nk eps_mnk L_imunu^n * S_imunu^k
    !! t_m = lambda_i * sum_alpha_beta * sum_nk_gamma,zeta eps_mnk L^n_gamma,zeta * sigma^k_alpha_beta
    do i = 1, s%nAtoms
      do m = 1, 3
        do n = 1, 3
          do k = 1, 3
            torque(     1:nOrb ,     1:nOrb ,m,i) = torque(     1:nOrb ,     1:nOrb ,m,i) + lvec(1:nOrb,1:nOrb,n) * sigma(1,1,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,     1:nOrb ,m,i) = torque(nOrb+1:nOrb2,     1:nOrb ,m,i) + lvec(1:nOrb,1:nOrb,n) * sigma(2,1,k) * levi_civita(m,n,k)
            torque(     1:nOrb ,nOrb+1:nOrb2,m,i) = torque(     1:nOrb ,nOrb+1:nOrb2,m,i) + lvec(1:nOrb,1:nOrb,n) * sigma(1,2,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) = torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) + lvec(1:nOrb,1:nOrb,n) * sigma(2,2,k) * levi_civita(m,n,k)
          end do
        end do
      end do

      ! p-block
      torque( 2: 4, 2: 4,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque( 2: 4, 2: 4,:,i)
      torque( 2: 4,11:13,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque( 2: 4,11:13,:,i)
      torque(11:13, 2: 4,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque(11:13, 2: 4,:,i)
      torque(11:13,11:13,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque(11:13,11:13,:,i)

      ! d-block
      torque( 5: 9, 5: 9,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque( 5: 9, 5: 9,:,i)
      torque( 5: 9,14:18,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque( 5: 9,14:18,:,i)
      torque(14:18, 5: 9,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque(14:18, 5: 9,:,i)
      torque(14:18,14:18,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque(14:18,14:18,:,i)
    end do
  end subroutine local_SO_torque

end module mod_TCM
