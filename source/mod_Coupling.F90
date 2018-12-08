module mod_Coupling
  use mod_f90_kind, only: double
  implicit none
  character(len=3), private :: folder = "Jij"
  character(len=3), dimension(3), private :: filename = ["Jij", "J  ", "Dz "]

  ! Exchange interaction
  real(double), dimension(:,:), allocatable :: trJij
  real(double), dimension(:,:,:,:), allocatable :: Jij,Jijs,Jija

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

  subroutine openCouplingFiles()
    use mod_System,     only: s => sys
    use mod_parameters, only: output
    implicit none
    character(len=400) :: varm
    integer :: i, j, iw

    do j=1,s%nAtoms
      do i=1,s%nAtoms
        iw = 2000 + (j-1) * s%nAtoms * 2 + (i-1) * 2
        if(i==j) then
          iw = iw + 1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(1)),i,trim(output%info),trim(output%BField),trim(output%SOC)
          open (unit=iw, file=varm,status='replace')
          write(unit=iw, fmt="('#   energy      ,  Jii_xx           ,   Jii_yy  ')")
          ! Anisotropy energy is given by K^a = 2*J_ii^aa
          ! omega_res ~ gamma*m_i*J_ii (*2?) ,
          ! where J_ii is the one calculated here
        else
          iw = iw + 1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,'_',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(2)),i,j,trim(output%info),trim(output%BField),trim(output%SOC)
          open (unit=iw, file=varm,status='replace')
          write(unit=iw, fmt="('#   energy      ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
          iw = iw + 1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,'_',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(3)),i,j,trim(output%info),trim(output%BField),trim(output%SOC)
          open (unit=iw, file=varm,status='replace')
          write(unit=iw, fmt="('#   energy      , Dz = (Jxy - Jyx)/2       ')")
        end if
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

  subroutine writeCoupling(e)
    use mod_f90_kind, only: double
    use mod_System,     only: s => sys
    implicit none
    real(double), intent(in) :: e
    integer :: i, j, iw

    do j=1,s%nAtoms
      do i=1,s%nAtoms
        iw = 2000 + (j-1) * s%nAtoms * 2 + (i-1) * 2
        if(i==j) then
          iw = iw + 1
          write(unit=iw,fmt="(3(es16.9,2x))") e,Jij(i,j,1,1),Jij(i,j,2,2)
        else
          iw = iw + 1
          write(unit=iw,fmt="(4(es16.9,2x))") e,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
          iw = iw + 1
          write(unit=iw,fmt="(2(es16.9,2x))") e,Jija(i,j,1,2)
        end if
      end do
    end do
  end subroutine writeCoupling

! This subroutine calculates Js and Ks from the inverse static susceptibility
! It only works for FM or AFM alignments
  subroutine get_J_K_from_chi()
    use mod_f90_kind,         only: double
    use mod_susceptibilities, only: schi
    use mod_parameters,       only: output,sigmai2i
    use mod_system,           only: s => sys
    use mod_magnet,           only: mabs,mvec_cartesian
    use mod_mpi_pars,         only: abortProgram
    implicit none
    integer :: AllocateStatus
    integer,         dimension(:),     allocatable :: fmalign
    real(double),    dimension(:),     allocatable :: Kx_minus_Ky, Kx_plus_Ky, Kx, Ky, Jtotal
    real(double),    dimension(:,:),   allocatable :: Jij_chi
    complex(double), dimension(:,:),   allocatable :: chi_inv
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
    Jtotal  = 0.d0
    Jij_chi = 0.d0
    do j = 1,s%nAtoms
      do i = 1,s%nAtoms
        if( (i == j) .or. (mabs(i) < 1.d-8) .or. (mabs(j) < 1.d-8) ) cycle

        if(dot_product(mvec_cartesian(:,i),mvec_cartesian(:,j))>0.d0) then
          ! FM alignment:
          ! $J_{ij} = [\chi_{SS}^{-1}]^{+-}_{ij} |M_i| |M_j| / 2$
          Jij_chi(i,j) = 0.5d0*real(chi_inv(sigmai2i(1,i),sigmai2i(1,j)))*mabs(i)*mabs(j)*13606.d0
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
          Jij_chi(i,j) = -0.5d0*real(chi_inv(sigmai2i(1,i),sigmai2i(4,j)))*mabs(i)*mabs(j)*13606.d0
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
      if(mabs(i) < 1.d-8) cycle

      Kx_plus_Ky(i)  = 0.5d0*real(chi_inv(sigmai2i(1,i),sigmai2i(1,i)))*mabs(i)*mabs(i)*13606.d0 + Jtotal(i)
      Kx_minus_Ky(i) = 0.5d0*real(chi_inv(sigmai2i(1,i),sigmai2i(4,i)))*mabs(i)*mabs(i)*13606.d0

      Kx(i) = 0.5d0*( Kx_plus_Ky(i) + Kx_minus_Ky(i) )
      Ky(i) = 0.5d0*( Kx_plus_Ky(i) - Kx_minus_Ky(i) )

      write(output%unit_loop,"('Kx(',i0,') = ',es16.9,' meV , Ky(',i0,') = ',es16.9,' meV ')") i, Kx(i), i, Ky(i)

      if(Kx_minus_Ky(i)<1.d-8) &
        write(output%unit_loop,"('-> Uniaxial anisotropy: Kx(',i0,') = Ky(',i0,')')") i,i
      if((Kx(i)<0.d0) .and. (Ky(i)<0.d0)) then
        write(output%unit_loop,"('-> z is easy axis: Kx(',i0,'), Ky(',i0,') < 0')") i,i
      else if((Kx(i)>0.d0) .and. (Ky(i)>0.d0)) then
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
