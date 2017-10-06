module mod_alpha
  use mod_f90_kind, only: double
  implicit none

  complex(double), dimension(:,:),   allocatable :: m_chi, m_chi_hf, m_chi_inv, m_chi_hf_inv, m_chiorb_inv, m_chiorbhf_inv
  character(len=7), parameter, private :: folder = "A/Slope"
  character(len=8), dimension(5), parameter, private :: filename = ["chi     ", "chihf   ", "chiinv  ", "chihfinv", "sumrule "]

  complex(double) :: chi0_0
contains

  subroutine allocate_alpha()
    use mod_System, only: s => sys
    use mod_parameters, only: dim
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer :: AllocateStatus
    allocate(m_chi       (4*s%nAtoms,4*s%nAtoms), &
             m_chi_hf    (4*s%nAtoms,4*s%nAtoms), &
             m_chi_inv   (4*s%nAtoms,4*s%nAtoms), &
             m_chi_hf_inv(4*s%nAtoms,4*s%nAtoms), &
             m_chiorb_inv(dim,dim), &
             m_chiorbhf_inv(dim,dim), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocate_alpha] Not enough memory for: m_chi, m_chi_hf, m_chi_inv, m_chi_hf_inv")

    return
  end subroutine allocate_alpha

  subroutine deallocate_alpha()
    implicit none

    if(allocated(m_chi)) deallocate(m_chi)
    if(allocated(m_chi_hf)) deallocate(m_chi_hf)
    if(allocated(m_chi_inv)) deallocate(m_chi_inv)
    if(allocated(m_chi_hf_inv)) deallocate(m_chi_hf_inv)

    return
  end subroutine deallocate_alpha

  subroutine create_alpha_files()
    use mod_System, only: s => sys
    use mod_BrillouinZone, only: BZ
    use mod_parameters, only: strSites, eta, Utype, suffix, fieldpart
    use mod_SOC, only: SOCc, socpart
    use EnergyIntegration, only: strEnergyParts
    implicit none
    character(len=500)  :: varm
    integer :: i,j

    do i=1, s%nAtoms
      do j = 1, size(filename)-1
         write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder),trim(filename(j)),i,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
         open (unit=55+(j-1)*s%nAtoms+i, file=varm, status='replace', form='formatted')
         write(unit=55+(j-1)*s%nAtoms+i, fmt="('#     energy    ,  alpha   ,  gamma  ,  alpha/gamma')")
      end do
      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder),trim(filename(5)),i,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
      open (unit=55+4*s%nAtoms+i, file=varm, status='replace', form='formatted')
      write(unit=55+4*s%nAtoms+i, fmt="('#     energy    ,  gammaM   ,  (ReX(w))^(-2),     U^2,   v1  ,  v2  ,  v3,    v4')")
    end do

    return

  end subroutine create_alpha_files

  subroutine open_alpha_files()
    use mod_parameters, only: fieldpart, suffix, eta, Utype, strSites
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use mod_BrillouinZone, only: BZ
    use EnergyIntegration, only: strEnergyParts
    use mod_mpi_pars, only: abortProgram
    implicit none

    character(len=500) :: varm
    character(len=500) :: missing_files
    integer :: i, j,err, errt
    errt = 0

    do i=1, s%nAtoms
      do j = 1, size(filename)
         write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder),trim(filename(j)),i,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
         open (unit=55+(j-1)*s%nAtoms+i, file=varm, status='old', position='append', form='formatted', iostat=err)
         errt = errt + err
         if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
      end do
    end do
    if(errt/=0) call abortProgram("[open_alpha_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))

    return
  end subroutine open_alpha_files

  subroutine close_alpha_files()
    use mod_system, only: s => sys
    implicit none
    integer :: i

    do i = 1, s%nAtoms*size(filename)
      close(unit=55+i)
    end do

    return
  end subroutine close_alpha_files

  subroutine write_alpha(e)
    use mod_f90_kind, only: double
    use mod_constants, only: cI, cZero
    use mod_susceptibilities, only: schi, schihf, chiorb, chiorb_hf
    use mod_parameters, only: sigmai2i, U, dim, sigmaimunu2i
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: mabs
    implicit none
    real(double) :: e
    real(double) :: gammaM, alpha_v1, alpha_v2, alpha_v3, alpha_v4
    complex(double) :: axx, axy, axx_hf, axy_hf, axx_inv, axy_inv, axx_hf_inv, axy_hf_inv
    integer :: i,j,sigma, sigmap, mu, nu

    call open_alpha_files()

    m_chiorb_inv = chiorb
    m_chiorbhf_inv = chiorb_hf

    do i = 1, s%nAtoms
      do j = 1, s%nAtoms
        do sigma = 1, 4
          do sigmap = 1, 4
            m_chi   (sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schi  (sigma,sigmap,i,j)
            m_chi_hf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schihf(sigma,sigmap,i,j)
          end do
        end do
      end do
    end do
   !  m_chi_inv    = m_chi
   !  m_chi_hf_inv = m_chi_hf
   !  call invers(m_chi_inv,    4*s%nAtoms)
   !  call invers(m_chi_hf_inv, 4*s%nAtoms)

    m_chi_inv = cZero
    m_chi_hf_inv = cZero

    do i = 1, s%nAtoms
      do j = 1, s%nAtoms
         do sigma = 1, 4
            do sigmap = 1, 4
               do mu = 5, nOrb
                  do nu = 5, nOrb
                     m_chi_inv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = m_chi_inv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                     m_chi_hf_inv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = m_chi_hf_inv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                  end do
               end do
            end do
         end do
      end do
    end do

    call invers(m_chi_inv, 4*s%nAtoms)
    call invers(m_chi_hf_inv, 4*s%nAtoms)


    do i = 1, s%nAtoms
      axx_inv    = 0.5d0   *(m_chi_inv(sigmai2i(1,i),sigmai2i(4,i)) + m_chi_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_inv(sigmai2i(4,i),sigmai2i(4,i)) + m_chi_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axy_inv    = 0.5d0*cI*(m_chi_inv(sigmai2i(1,i),sigmai2i(4,i)) - m_chi_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_inv(sigmai2i(4,i),sigmai2i(4,i)) - m_chi_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axx_hf_inv = 0.5d0   *(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(4,i)) + m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf_inv(sigmai2i(4,i),sigmai2i(4,i)) + m_chi_hf_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axy_hf_inv = 0.5d0*cI*(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(4,i)) - m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf_inv(sigmai2i(4,i),sigmai2i(4,i)) - m_chi_hf_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axx        = 0.5d0   *(m_chi(sigmai2i(1,i),sigmai2i(4,i)) + m_chi(sigmai2i(1,i),sigmai2i(1,i)) + m_chi(sigmai2i(4,i),sigmai2i(4,i)) + m_chi(sigmai2i(4,i),sigmai2i(1,i)))
      axy        = 0.5d0*cI*(m_chi(sigmai2i(1,i),sigmai2i(4,i)) - m_chi(sigmai2i(1,i),sigmai2i(1,i)) + m_chi(sigmai2i(4,i),sigmai2i(4,i)) - m_chi(sigmai2i(4,i),sigmai2i(1,i)))
      axx_hf     = 0.5d0   *(m_chi_hf(sigmai2i(1,i),sigmai2i(4,i)) + m_chi_hf(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf(sigmai2i(4,i),sigmai2i(4,i)) + m_chi_hf(sigmai2i(4,i),sigmai2i(1,i)))
      axy_hf     = 0.5d0*cI*(m_chi_hf(sigmai2i(1,i),sigmai2i(4,i)) - m_chi_hf(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf(sigmai2i(4,i),sigmai2i(4,i)) - m_chi_hf(sigmai2i(4,i),sigmai2i(1,i)))

      gammaM = e / aimag(axy_inv)

      if(e == 0.d0) then
         chi0_0 = 1.d0/real(m_chi_hf(sigmai2i(1,i),sigmai2i(1,i)))**2
      end if
      alpha_v1 = - gammaM * aimag(m_chi_inv(sigmai2i(1,i),sigmai2i(1,i))) / e
      alpha_v2 = - gammaM * aimag(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i))) / e
      alpha_v3 =   gammaM * aimag(m_chi_hf(sigmai2i(1,i),sigmai2i(1,i))) / (e * chi0_0)
      alpha_v4 =   gammaM * aimag(m_chi_hf(sigmai2i(1,i),sigmai2i(1,i))) / e * U(i)**2

      !print *, real(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i)))**2
      !print *, 1.d0/real(m_chi_hf(sigmai2i(1,i),sigmai2i(1,i)))**2
      write(55 +            i,"(4(es16.9,2x))") e, -1.d0 * aimag(axx       )/aimag(axy       ), e / (mabs(i) * aimag(axy       )), -mabs(i)*aimag(axx       ) / e
      write(55 +   s%nAtoms+i,"(4(es16.9,2x))") e, -1.d0 * aimag(axx_hf    )/aimag(axy_hf    ), e / (mabs(i) * aimag(axy_hf    )), -mabs(i)*aimag(axx_hf    ) / e
      write(55 + 2*s%nAtoms+i,"(4(es16.9,2x))") e, -1.d0 * aimag(axx_inv   )/aimag(axy_inv   ), e / (mabs(i) * aimag(axy_inv   )), -mabs(i)*aimag(axx_inv   ) / e
      write(55 + 3*s%nAtoms+i,"(4(es16.9,2x))") e, -1.d0 * aimag(axx_hf_inv)/aimag(axy_hf_inv), e / (mabs(i) * aimag(axy_hf_inv)), -mabs(i)*aimag(axx_hf_inv) / e
      write(55 + 4*s%nAtoms+i,"(8(es16.9,2x))") e, gammaM, 1.d0/real(m_chi_hf(sigmai2i(1,i),sigmai2i(1,i)))**2, U(i)**2, alpha_v1, alpha_v2, alpha_v3, alpha_v4
    end do

    call close_alpha_files()

    return
  end subroutine write_alpha




end module mod_alpha
