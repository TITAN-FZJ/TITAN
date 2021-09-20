module mod_TCM
!> This subroutine includes the calculation of the Gilbert damping
!> and the related subroutines.
  use mod_kind, only: dp
  implicit none
  character(len=5),               private :: folder = "A/TCM"
  character(len=3), dimension(4), private :: filename = ["SO ", "SOi", "XC ", "XCi"]

contains

  subroutine openTCMFiles()
  !> Open files for the results obtained with the torque-correlation method
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


  subroutine writeTCM(unit,alpha,ndiffk,diff_k,ialpha,iwght)
  !> Write the results for the damping calculated with the TCM
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

  subroutine closeTCMFiles()
  !> Close files for the results obtained with the torque-correlation method
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

  subroutine calculate_TCM()
  !> Driver routine to call the different torque-correlation methods and their I/O
    use mod_kind,       only: dp
    use mod_parameters, only: output
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: rField
    use mod_progress,   only: write_time
    use mod_torques,    only: SO_torque_operator,xc_torque_operator
    implicit none
    integer :: ndiffk
    real(dp),    dimension(:),     allocatable    :: diff_k, iwght
    real(dp),    dimension(:,:,:), allocatable    :: ialphaSO, ialphaXC
    complex(dp), dimension(s%nAtoms,s%nAtoms,3,3) :: alphaSO, alphaXC

    if(rField == 0) &
      write(output%unit_loop, "('[calculate_TCM] Starting to calculate Gilbert damping within the TCM:')")

    if(rField == 0) &
      call openTCMFiles()

    !! What is calculated in TCM() is:
    !! $\alpha= \frac{1}{M \pi N} \sum_{\mathbf{k}} \operatorname{Tr}\{\operatorname{Im} G(\mathbf{k},\epsilon_\text{F})\hat{T}^-_\text{xc} \operatorname{Im} G(\mathbf{k},\epsilon_F) \hat{T}^+_\text{xc}\}$
    if(rField == 0) &
      write(output%unit_loop, "('[calculate_TCM] SO-TCM...')")

    !! SO-TCM: The missing pre-factor is 2*\gamma, which, for g = 2 results in 4.0

    call TCM(SO_torque_operator, alphaSO, ndiffk, diff_k, ialphaSO, iwght)

    if(rField == 0) then
      alphaSO = 4._dp * alphaSO
      ialphaSO = 4._dp * ialphaSO
      write(output%unit_loop, "('[calculate_TCM] Writing results of SO-TCM to file...')")
      call writeTCM(634894, alphaSO, ndiffk, diff_k, ialphaSO, iwght)
      call write_time('[calculate_TCM] Time after SO-TCM: ',output%unit_loop)
      write(output%unit_loop, "('[calculate_TCM] XC-TCM...')")
    end if

    !! XC-TCM: The missing pre-factor is \gamma/2, which, for g = 2 results in 1.0
    call TCM(xc_torque_operator, alphaXC, ndiffk, diff_k, ialphaXC, iwght)

    if(rField == 0) then
      write(output%unit_loop, "('[calculate_TCM] Writing results of XC-TCM to file...')")
      call writeTCM(634896, alphaXC, ndiffk, diff_k, ialphaXC, iwght)
      call write_time('[calculate_TCM] Time after XC-TCM: ',output%unit_loop)
      call closeTCMFiles()
    end if

  end subroutine calculate_TCM


  subroutine TCM(torque_operator,alpha,ndiffk,diff_k,ialpha,iwght)
  !> This subroutine calculates the Gilbert damping using the torque-correlation method
    use mod_kind,          only: dp,int64
    use mod_constants,     only: cZero,pi,cOne,cI
    use mod_System,        only: s => sys
    use mod_BrillouinZone, only: realBZ,store_diff
    use mod_parameters,    only: eta,output,kptotal_in
    use mod_magnet,        only: mabs
    use mod_tools,         only: sort,itos,rtos
    use mod_hamiltonian,   only: hamilt_local
    use mod_greenfunction, only: green
    use mod_mpi_pars,      only: rField,sField,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FieldComm,ierr
    implicit none
    interface
      subroutine torque_operator(torque)
        use mod_kind,       only: dp
        use mod_System,     only: s => sys
        implicit none
        complex(dp), dimension(s%nOrb2,s%nOrb2,3,s%nAtoms), intent(out) :: torque
      end subroutine torque_operator
    end interface

    complex(dp), dimension(s%nAtoms,s%nAtoms,3,3), intent(out) :: alpha
    !! Contains the Gilbert Damping as matrix in sites and cartesian coordinates (nAtoms,nAtoms,3,3)
    integer(int64) :: iz
    integer   :: i, j, m, n, k, mu, ios, nOrb2_i, nOrb2_j

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
    complex(dp), dimension(s%nOrb2,s%nOrb2,3,s%nAtoms) :: torque
    !! Torque operator (SO or XC)
    complex(dp),dimension(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms) :: gf
    !! Green function
    complex(dp),dimension(s%nOrb2,s%nOrb2) :: temp1, temp2, temp3, temp4
    complex(dp) :: alphatemp
    !! Temporary k-dependent alpha
    complex(dp), dimension(s%nAtoms,s%nAtoms,3,3) :: alpha_loc
    !! Local (in-processor) alpha

    external :: zgemm,MPI_Allreduce

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

    call torque_operator(torque)

    !! XC-TCM:
    !! $\alpha= \frac{\gamma}{2 M \pi N} \sum_{\mathbf{k}} \operatorname{Tr}\{\operatorname{Im} G(\mathbf{k},\epsilon_\text{F})\hat{T}^-_\text{xc} \operatorname{Im} G(\mathbf{k},\epsilon_F) \hat{T}^+_\text{xc}\}$
    !! SO-TCM:
    !! $\alpha - \alpha_{\text{noSOI}} = 2 \frac{\gamma}{M \pi N} \sum_{\mathbf{k}}\operatorname{Tr}\{\operatorname{Im} G(\mathbf{k},\epsilon_F)\hat{T}^-_{\text{SOI}} \operatorname{Im} G(\mathbf{k},\epsilon_F) \hat{T}^+_{\text{SOI}}\} $
    !! g = 2
    !! M = magnetic moment
    alpha(:,:,:,:) = cZero

    ! Build local hamiltonian
    call hamilt_local(s)

    if(rField == 0) &
      write(output%unit_loop, "('[TCM] Calculating alpha...')")

    !$omp parallel default(none) reduction(+:ialpha,iwght) &
    !$omp& private(m,n,i,j,k,mu,nOrb2_i,nOrb2_j,iz,kp,wght,gf,temp1,temp2,temp3,temp4,alpha_loc,alphatemp) &
    !$omp& shared(s,realBZ,eta,torque,alpha,rField,order,ndiffk,diff_k)
    alpha_loc  = cZero

    !$omp do schedule(static)
    kpoints: do iz = 1, realBZ%workload
      kp   = realBZ%kp(1:3,iz)
      wght = realBZ%w(iz)
      gf   = cZero
      call green(s%Ef,eta,s,kp,gf)
      site_i: do i = 1, s%nAtoms
        site_j: do j = 1, s%nAtoms
          nOrb2_i = 2*s%Types(s%Basis(i)%Material)%nOrb
          nOrb2_j = 2*s%Types(s%Basis(j)%Material)%nOrb
          temp2 = 0.5_dp * cI * ( transpose(conjg(gf(1:nOrb2_i,1:nOrb2_j,i,j))) - gf(1:nOrb2_j,1:nOrb2_i,j,i) ) ! Im(G_ji(Ef)) )
          temp4 = 0.5_dp * cI * ( transpose(conjg(gf(1:nOrb2_j,1:nOrb2_i,j,i))) - gf(1:nOrb2_i,1:nOrb2_j,i,j) ) ! Im(G_ij(Ef))

          dir_m: do m = 1, 3
            dir_n: do n = 1, 3

              alphatemp = cZero
              ! alpha^{mn}_{ij} = Tr ( Torque^m_i * Im(G_ij(Ef)) * Torque^n_j * Im(G_ji(Ef)) )
              !                      temp2(below) *     temp4    *    temp3   * temp2(above)
              temp3 = torque(1:nOrb2_j,1:nOrb2_j,n,j)
              call zgemm('n','n',nOrb2_j,nOrb2_i,nOrb2_j,cOne,temp3,s%nOrb2,temp2,s%nOrb2,cZero,temp1,s%nOrb2) ! Torque^n_j * Im(G_ji(Ef))
              call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,cOne,temp4,s%nOrb2,temp1,s%nOrb2,cZero,temp3,s%nOrb2) ! Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))
              temp2 = torque(1:nOrb2_i,1:nOrb2_i,m,i)
              call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_i,cOne,temp2,s%nOrb2,temp3,s%nOrb2,cZero,temp1,s%nOrb2) ! Torque^m_i * Im(G_ij(Ef) * Torque^n_j * Im(G_ji(Ef))

              do mu = 1, nOrb2_i
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

            end do dir_n
          end do dir_m
        end do site_j
      end do site_i
    end do kpoints
    !$omp end do nowait

    !$omp critical
    alpha = alpha + alpha_loc
    !$omp end critical

    !$omp end parallel

    !! XC-TCM:
    !! $\alpha= \frac{\gamma}{2 M \pi N} \sum_{\mathbf{k}} \operatorname{Tr}\{\operatorname{Im} G(\mathbf{k},\epsilon_\text{F})\hat{T}^-_\text{xc} \operatorname{Im} G(\mathbf{k},\epsilon_F) \hat{T}^+_\text{xc}\}$
    !! SO-TCM:
    !! $\alpha - \alpha_{\text{noSOI}} = 2 \frac{\gamma}{M \pi N} \sum_{\mathbf{k}}\operatorname{Tr}\{\operatorname{Im} G(\mathbf{k},\epsilon_F)\hat{T}^-_{\text{SOI}} \operatorname{Im} G(\mathbf{k},\epsilon_F) \hat{T}^+_{\text{SOI}}\} $
    !! g = 2
    !! M = magnetic moment
    do i = 1, s%nAtoms
      do j = 1, s%nAtoms
        alpha(j,i,:,:) = alpha(j,i,:,:) / (sqrt(mabs(i)*mabs(j))*pi) ! The sqrt shows up because the i=j term multiplies twice
      end do
      if(s%isysdim==3) ialpha(i,:,:) = ialpha(i,:,:) / (mabs(i)*pi)
    end do

    call MPI_Allreduce(MPI_IN_PLACE, alpha , s%nAtoms*s%nAtoms*3*3, MPI_DOUBLE_COMPLEX  , MPI_SUM, FieldComm, ierr)
    if(s%isysdim==3) then
      call MPI_Allreduce(MPI_IN_PLACE, ialpha, s%nAtoms*3*ndiffk    , MPI_DOUBLE_PRECISION, MPI_SUM, FieldComm, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, iwght , ndiffk               , MPI_DOUBLE_PRECISION, MPI_SUM, FieldComm, ierr)
      deallocate(order,diff_k_unsrt)
    end if

  end subroutine TCM

end module mod_TCM
