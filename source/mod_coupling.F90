module mod_Coupling
  use mod_kind, only: dp
  implicit none
  character(len=3), private :: folder = "Jij"
  character(len=3), dimension(1), private :: filename = ["Jij"]

  ! Exchange interaction
  complex(dp), dimension(:,:,:,:), allocatable :: Jij
  complex(dp), dimension(:,:,:,:,:), allocatable :: Jij_q
  real(dp), dimension(:,:), allocatable :: trJij
  real(dp), dimension(:,:,:,:), allocatable :: Jijs,Jija

contains

  subroutine allocateCoupling()
    use mod_System,     only: s => sys
    implicit none

    if(allocated(trJij)) deallocate(trJij)
    if(allocated(Jij)) deallocate(Jij)
    if(allocated(Jijs)) deallocate(Jijs)
    if(allocated(Jija)) deallocate(Jija)
    allocate(trJij(s%nAtoms,s%nAtoms))
    allocate(Jij(s%nAtoms,s%nAtoms,3,3))
    allocate(Jijs(s%nAtoms,s%nAtoms,3,3))
    allocate(Jija(s%nAtoms,s%nAtoms,3,3))

  end subroutine allocateCoupling

  subroutine deallocateCoupling()
    implicit none
    if(allocated(trJij)) deallocate(trJij)
    if(allocated(Jij)) deallocate(Jij)
    if(allocated(Jijs)) deallocate(Jijs)
    if(allocated(Jija)) deallocate(Jija)
  end subroutine deallocateCoupling

  subroutine allocateRealCoupling()
    use mod_System,     only: s => sys
    use mod_BrillouinZone, only: realBZ
    implicit none

    if(allocated(Jij)) deallocate(Jij)
    if(allocated(Jij_q)) deallocate(Jij_q)
    allocate(Jij(s%nAtoms,s%nAtoms,3,3))
    allocate(Jij_q(realBZ%workload,s%nAtoms,s%nAtoms,3,3))

  end subroutine allocateRealCoupling

  subroutine deallocateRealCoupling()
    implicit none
    if(allocated(Jij)) deallocate(Jij)
    if(allocated(Jij_q)) deallocate(Jij_q)
  end subroutine deallocateRealCoupling

  subroutine openCouplingFiles()
    use mod_System,     only: s => sys
    use mod_parameters, only: output,kdirection
    use mod_io,         only: write_header
    implicit none
    character(len=400) :: varm
    integer :: i, j, iw

    do j=1,s%nAtoms
      do i=1,s%nAtoms
        iw = 2000 + (j-1) * s%nAtoms * 2 + (i-1) * 2 + 1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,'_',i0,'_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(1)),i,j,trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open (unit=iw, file=varm,status='replace')
        call write_header(iw,"#   energy      ,   q-vector    ,  ((real(Jij(i,j,m,n)),aimag(Jij(i,j,m,n)), m=1,3), n=1,3) ")
        ! Anisotropy energy is given by K^a = 2*J_ii^aa
        ! omega_res ~ gamma*m_i*J_ii (*2?) ,
        ! where J_ii is the one calculated here
      end do
    end do
  end subroutine openCouplingFiles

  subroutine closeCouplingFiles()
    use mod_System,     only: s => sys
    implicit none
    integer :: i, j, iw
    do j=1,s%nAtoms
      do i=1,s%nAtoms
        iw = 2000 + (j-1) * s%nAtoms * 2 + (i-1)*2 + 1
        close (iw)
        if(i/=j) close(iw+1)
      end do
    end do
  end subroutine closeCouplingFiles

  subroutine openRealCouplingFiles()
    use mod_System,     only: s => sys
    use mod_parameters, only: output
    use mod_io,         only: write_header
    implicit none
    character(len=400) :: varm
    integer :: i, j, iw

    do j=1,s%nAtoms
      do i=1,s%nAtoms
        iw = 2000 + (j-1) * s%nAtoms * 2 + (i-1) * 2 + 1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,'_',i0,'_kdir=',a,a,a,a,'_real.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(1)),i,j,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open (unit=iw, file=varm,status='replace')
        call write_header(iw,"#     Rji_x     ,      Rji_y      ,      Rji_z      ,      |Rji|      ,  Jii_xx         ," // &
        &"  Jii_xy         ,  Jii_xz         ,  Jii_yz         ,  Jii_yy         ,  Jii_yz         ," // &
        &"  Jii_zx         ,  Jii_zy         ,  Jii_zz         ")
        ! Anisotropy energy is given by K^a = 2*J_ii^aa
        ! omega_res ~ gamma*m_i*J_ii (*2?) ,
        ! where J_ii is the one calculated here
      end do
    end do
  end subroutine openRealCouplingFiles

  subroutine writeCoupling(e,q)
    use mod_kind, only: dp
    use mod_System,     only: s => sys
    implicit none
    real(dp), intent(in) :: e,q
    integer :: i, j, m, n, iw

    do j=1,s%nAtoms
      do i=1,s%nAtoms
        iw = 2000 + (j-1) * s%nAtoms * 2 + (i-1) * 2 + 1
        write(unit=iw,fmt="(20(es16.9,2x))") e,q, ((real(Jij(i,j,m,n)),aimag(Jij(i,j,m,n)), m=1,3), n=1,3)
      end do
    end do
  end subroutine writeCoupling

! This subroutine calculates Js and Ks from the inverse static susceptibility
! It only works for FM or AFM alignments
  subroutine get_J_K_from_chi()
    use mod_kind,             only: dp
    use mod_susceptibilities, only: schi
    use mod_parameters,       only: output,sigmai2i
    use mod_system,           only: s => sys
    use mod_magnet,           only: mabs,mvec_cartesian
    use mod_mpi_pars,         only: abortProgram
    use mod_tools,            only: invers
    implicit none
    integer :: AllocateStatus
    integer,     dimension(:),     allocatable :: fmalign
    real(dp),    dimension(:),     allocatable :: Kx_minus_Ky, Kx_plus_Ky, Kx, Ky, Jtotal
    real(dp),    dimension(:,:),   allocatable :: Jij_chi
    complex(dp), dimension(:,:),   allocatable :: chi_inv
    integer :: i, j

    write(output%unit_loop, "('**********************************************************************')")
    write(output%unit_loop, "('[get_J_K_from_chi] Calculating J and K from the static susceptibility:')")

    allocate(chi_inv(4*s%nAtoms,4*s%nAtoms), &
            Jij_chi(s%nAtoms,s%nAtoms), &
            fmalign(s%nAtoms), &
            Jtotal(s%nAtoms), &
            Kx_minus_Ky(s%nAtoms), &
            Kx_plus_Ky(s%nAtoms), &
            Kx(s%nAtoms), &
            Ky(s%nAtoms), &
            stat = AllocateStatus)
    if(AllocateStatus /= 0) &
      call abortProgram("[get_J_K_from_chi] Not enough memory for: chi_inv, Jij_chi, fmalign, Jtotal, Kx_minus_Ky, Kx_plus_Ky, Kx, Ky")

    ! Inverting the spin-spin susceptibility to obtain $(\chi_{SS})^{-1}$
    chi_inv = schi
    call invers(chi_inv,4*s%nAtoms)

    ! Obtaning the couplings J between atoms in the unit cell (in meV)
    ! Definition: $E = -(J/2) \sum_{ij} m_i.m_j = -(J/2) \sum_{ij} M_i.M_j/(|M_i| |M_j|)$
    ! i.e.: J > 0 - FM coupling
    !       J < 0 - AFM coupling
    Jtotal  = 0._dp
    Jij_chi = 0._dp
    do j = 1,s%nAtoms
      do i = 1,s%nAtoms
        if( (i == j) .or. (mabs(i) < 1.e-8_dp) .or. (mabs(j) < 1.e-8_dp) ) cycle

        if(dot_product(mvec_cartesian(:,i),mvec_cartesian(:,j))>0._dp) then
          ! FM alignment:
          ! $J_{ij} = [\chi_{SS}^{-1}]^{+-}_{ij} |M_i| |M_j| / 2$
          Jij_chi(i,j) = 0.5_dp*real(chi_inv(sigmai2i(1,i),sigmai2i(1,j)))*mabs(i)*mabs(j)*13606._dp
          write(output%unit_loop,"('J(',i0,',',i0,') = ',es16.9,' meV (FM alignment)')") i,j,Jij_chi(i,j)
          Jtotal(i) = Jtotal(i) + Jij_chi(i,j)

          write(output%unit_loop,"('Results indicate that sites have ')", advance="no")
          if(Jij_chi(i,j)>0) then
            write(output%unit_loop,"('FM alignment with FM ground state')")
          else
            write(output%unit_loop,"('FM alignment with AFM ground state')")
          end if
        else
          ! AFM alignment:
          ! $J_{ij} = -[\chi_{SS}^{-1}]^{++}_{ij} |M_i| |M_j| / 2$
          Jij_chi(i,j) = -0.5_dp*real(chi_inv(sigmai2i(1,i),sigmai2i(4,j)))*mabs(i)*mabs(j)*13606._dp
          write(output%unit_loop,"('J(',i0,',',i0,') = ',es16.9,' meV (AFM alignment)')") i,j,Jij_chi(i,j)
          Jtotal(i) = Jtotal(i) - Jij_chi(i,j)

          write(output%unit_loop,"('Results indicate that sites have ')", advance="no")
          if(Jij_chi(i,j)>0) then
            write(output%unit_loop,"('AFM alignment with FM ground state')")
          else
            write(output%unit_loop,"('AFM alignment with AFM ground state')")
          end if
        end if

      end do
    end do

    ! Obtaning the anisotropies Kx and Ky (Kz=0) for each atom in the unit cell
    ! Definition: $E = - \sum_{i\mu} K_{i,\mu} (m_i.\hat{e}_\mu)^2 = - \sum_{i\mu} K_{i,\mu}/M_i^2 (M_i.\hat{e}_\mu)^2$
    ! i.e.: Kx,Ky > 0 - z is hard axis
    !       Kx,Ky < 0 - z is easy axis
    do i = 1,s%nAtoms
      if(mabs(i) < 1.e-8_dp) cycle

      Kx_plus_Ky(i)  = 0.5_dp*real(chi_inv(sigmai2i(1,i),sigmai2i(1,i)))*mabs(i)*mabs(i)*13606._dp + Jtotal(i)
      Kx_minus_Ky(i) = 0.5_dp*real(chi_inv(sigmai2i(1,i),sigmai2i(4,i)))*mabs(i)*mabs(i)*13606._dp

      Kx(i) = 0.5_dp*( Kx_plus_Ky(i) + Kx_minus_Ky(i) )
      Ky(i) = 0.5_dp*( Kx_plus_Ky(i) - Kx_minus_Ky(i) )

      write(output%unit_loop,"('Kx(',i0,') = ',es16.9,' meV , Ky(',i0,') = ',es16.9,' meV ')") i, Kx(i), i, Ky(i)

      if(Kx_minus_Ky(i)<1.e-8_dp) &
        write(output%unit_loop,"('-> Uniaxial anisotropy: Kx(',i0,') = Ky(',i0,')')") i,i
      if((Kx(i)<0._dp) .and. (Ky(i)<0._dp)) then
        write(output%unit_loop,"('-> z is easy axis: Kx(',i0,'), Ky(',i0,') < 0')") i,i
      else if((Kx(i)>0._dp) .and. (Ky(i)>0._dp)) then
        write(output%unit_loop,"('-> z is hard axis: Kx(',i0,'), Ky(',i0,') > 0')") i,i
      end if

    end do

! TODO: Calculate Kx and Ky (for Kz=0)
! TODO: Add ground-state direction (easy/hard axis)
! TODO: Add Dzyaloshinskii-Moriya

!-dot_product(mvec_cartesian(:,i),mvec_cartesian(:,j))


    deallocate(chi_inv, Jij_chi, Kx_minus_Ky, Kx_plus_Ky, Kx, Ky)
    write(output%unit_loop, "('**********************************************************************')")

  end subroutine get_J_K_from_chi



end module mod_Coupling
