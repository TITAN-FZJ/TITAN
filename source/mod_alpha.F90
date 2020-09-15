module mod_alpha
  use mod_kind, only: dp
  implicit none

  complex(dp), dimension(:,:),   allocatable :: m_chi, m_chi_hf, m_chi_inv, m_chi_hf_inv
  character(len=7), parameter, private :: folder = "A/Slope"
  character(len=8), dimension(10), parameter, private :: filename = ["chi     ", "chihf   ", "chiinv  ", "chihfinv", "sumrule ", "chiLS   ", "chiSL   ", "chiLL   ", "chiSt   ", "totalchi"]

contains

  subroutine allocate_alpha()
    use mod_System, only: s => sys
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer :: AllocateStatus
    allocate(m_chi       (4*s%nAtoms,4*s%nAtoms), &
             m_chi_hf    (4*s%nAtoms,4*s%nAtoms), &
             m_chi_inv   (4*s%nAtoms,4*s%nAtoms), &
             m_chi_hf_inv(4*s%nAtoms,4*s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocate_alpha] Not enough memory for: m_chi, m_chi_hf, m_chi_inv, m_chi_hf_inv")

  end subroutine allocate_alpha

  subroutine deallocate_alpha()
    implicit none

    if(allocated(m_chi)) deallocate(m_chi)
    if(allocated(m_chi_hf)) deallocate(m_chi_hf)
    if(allocated(m_chi_inv)) deallocate(m_chi_inv)
    if(allocated(m_chi_hf_inv)) deallocate(m_chi_hf_inv)

  end subroutine deallocate_alpha

  subroutine create_alpha_files()
    use mod_System, only: s => sys
    use mod_parameters, only: output
    implicit none
    character(len=500)  :: varm
    integer :: i,j

    do i=1, s%nAtoms
      do j = 1, size(filename)-6
         write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(j)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
         open (unit=55+(j-1)*s%nAtoms+i, file=varm, status='replace', form='formatted')
         write(unit=55+(j-1)*s%nAtoms+i, fmt="('#     energy    ,  alpha   ,  gamma  ,  alpha/gamma  ,  ((real[chi(j,i)], imag[chi(j,i)], j=1,4),i=1,4)  ')")
      end do
      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(5)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=55+4*s%nAtoms+i, file=varm, status='replace', form='formatted')
      write(unit=55+4*s%nAtoms+i, fmt="('#     energy    ,  gammaM   ,  (ReX(w))^(-2),     Um^2,   v1  ,  v2  ,  v3,    v4')")

      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(6)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=55+5*s%nAtoms+i, file=varm, status='replace', form='formatted')
      write(unit=55+5*s%nAtoms+i, fmt="('#     energy    ,  gammaM   ,  alpha_LS')")

      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(7)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=55+6*s%nAtoms+i, file=varm, status='replace', form='formatted')
      write(unit=55+6*s%nAtoms+i, fmt="('#     energy    ,  gammaM   ,  alpha_SL')")

      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(8)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=55+7*s%nAtoms+i, file=varm, status='replace', form='formatted')
      write(unit=55+7*s%nAtoms+i, fmt="('#     energy    ,  gammaM   ,  alpha_LL')")

      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(9)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=55+8*s%nAtoms+i, file=varm, status='replace', form='formatted')
      write(unit=55+8*s%nAtoms+i, fmt="('#     energy    ,  gammaM   ,  alpha_St')")

      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(10)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=55+9*s%nAtoms+i, file=varm, status='replace', form='formatted')
      write(unit=55+9*s%nAtoms+i, fmt="('#     energy    ,  gammaM   ,  alpha_t')")

    end do
  end subroutine create_alpha_files

  subroutine open_alpha_files()
    use mod_parameters, only: output
    use mod_system, only: s => sys
    use mod_mpi_pars, only: abortProgram
    implicit none

    character(len=500) :: varm
    character(len=500) :: missing_files
    integer :: i, j,err, errt
    errt = 0

    do i=1, s%nAtoms
      do j = 1, size(filename)
         write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(j)),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
         open (unit=55+(j-1)*s%nAtoms+i, file=varm, status='old', position='append', form='formatted', iostat=err)
         errt = errt + err
         if(err /= 0) missing_files = trim(missing_files) // " " // trim(varm)
      end do
    end do
    if(errt/=0) call abortProgram("[open_alpha_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))

  end subroutine open_alpha_files

  subroutine close_alpha_files()
    use mod_system, only: s => sys
    implicit none
    integer :: i

    do i = 1, s%nAtoms*size(filename)
      close(unit=55+i)
    end do

  end subroutine close_alpha_files

  subroutine write_alpha(e)
    use mod_kind,             only: dp
    use mod_constants,        only: StoC, CtoS
    use mod_susceptibilities, only: schi, schihf, schiLS, schiSL, schiLL
    use mod_parameters,       only: sigmai2i, Um
    use mod_system,           only: s => sys
    use mod_magnet,           only: mabs
    use mod_SOC,              only: SOC
    use mod_tools,            only: invers
    implicit none
    real(dp) :: e
    real(dp) :: gammaM, alpha_v1, alpha_v2, alpha_v3, alpha_v4
    complex(dp), dimension(4,4) :: acart, acarthf, acarthfinv, acartinv
    complex(dp), dimension(3,3) :: acart_t_inv,acart_St_inv,acart_LS_inv,acart_SL_inv,acart_LL_inv
    integer :: i, p,q,r,t

    call allocate_alpha()

    call open_alpha_files()

    m_chi = schi
    m_chi_hf = schihf
    m_chi_inv = schi
    m_chi_hf_inv = schihf
    call invers(m_chi_inv,    4*s%nAtoms)
    call invers(m_chi_hf_inv, 4*s%nAtoms)

    acart = 0._dp
    acarthf = 0._dp
    acartinv = 0._dp
    acarthfinv = 0._dp
    do i = 1, s%nAtoms

      do p = 1, 4
        do q = 1, 4
          do r = 1, 4
            do t = 1, 4
              acart(p,q) = acart(p,q) + StoC(p,r) * m_chi(sigmai2i(r,i), sigmai2i(t,i)) * CtoS(t,q)
              acarthf(p,q) = acarthf(p,q) + StoC(p,r) * m_chi_hf(sigmai2i(r,i), sigmai2i(t,i)) * CtoS(t,q)
              acartinv(p,q) = acartinv(p,q) + StoC(p,r) * m_chi_inv(sigmai2i(t,i),sigmai2i(r,i)) * CtoS(t,q)
              acarthfinv(p,q) = acarthfinv(p,q) + StoC(p,r) * m_chi_hf_inv(sigmai2i(r,i), sigmai2i(t,i)) * CtoS(t,q)
            end do
          end do
        end do
      end do


      gammaM = e / aimag(acartinv(2,3))
      alpha_v1 = - gammaM * aimag(m_chi_inv(sigmai2i(1,i),sigmai2i(1,i))) / e
      alpha_v2 = - gammaM * aimag(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i))) / e
      alpha_v3 =   gammaM * aimag(m_chi_hf(sigmai2i(1,i),sigmai2i(1,i))) / e * real(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i)))**2
      alpha_v4 =   gammaM * aimag(m_chi_hf(sigmai2i(1,i),sigmai2i(1,i))) / e * Um(i)**2

      write(55 +            i,"(36(es16.9,2x))") e, -1._dp * aimag(acart(2,2))/aimag(acart(2,3)), e / (mabs(i) * aimag(acart(2,3))), -mabs(i)*aimag(acart(2,2)) / e, &
                                                ((real(acart(q,p)), aimag(acart(q,p)), q = 1, 4), p = 1, 4)
      write(55 +   s%nAtoms+i,"(36(es16.9,2x))") e, -1._dp * aimag(acarthf(2,2))/aimag(acarthf(2,3)), e / (mabs(i) * aimag(acarthf(2,3))), -mabs(i)*aimag(acarthf(2,2)) / e, &
                                                ((real(acarthf(q,p)), aimag(acarthf(q,p)), q = 1, 4), p = 1, 4)
      write(55 + 2*s%nAtoms+i,"(36(es16.9,2x))") e, -1._dp * aimag(acartinv(2,2))/aimag(acartinv(2,3)), e / (mabs(i) * aimag(acartinv(2,3))), -mabs(i)*aimag(acartinv(2,2)) / e, &
                                                ((real(acartinv(q,p)), aimag(acartinv(q,p)), q = 1, 4), p = 1, 4)
      write(55 + 3*s%nAtoms+i,"(36(es16.9,2x))") e, -1._dp * aimag(acarthfinv(2,2))/aimag(acarthfinv(2,3)), e / (mabs(i) * aimag(acarthfinv(2,3))), -mabs(i)*aimag(acarthfinv(2,2)) / e, &
                                                ((real(acarthfinv(q,p)), aimag(acarthfinv(q,p)), q = 1, 4), p = 1, 4)
      write(55 + 4*s%nAtoms+i,"(8(es16.9,2x))") e, gammaM, real(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i)))**2, Um(i)**2, alpha_v1, alpha_v2, alpha_v3, alpha_v4

      if(SOC) then
        do p = 1, 3
          do q = 1, 3
            acart_LS_inv(p,q) = schiLS(sigmai2i(p,i),sigmai2i(q,i))
            acart_SL_inv(p,q) = schiSL(sigmai2i(p,i),sigmai2i(q,i))
            acart_LL_inv(p,q) = schiLL(sigmai2i(p,i),sigmai2i(q,i))
            acart_St_inv(p,q) = 4._dp*acart(p+1,q+1) + 2._dp*schiSL(sigmai2i(p,i),sigmai2i(q,i))
            acart_t_inv(p,q)  = 4._dp*acart(p+1,q+1) + 2._dp*schiLS(sigmai2i(p,i),sigmai2i(q,i)) + 2._dp*schiSL(sigmai2i(p,i),sigmai2i(q,i)) + schiLL(sigmai2i(p,i),sigmai2i(q,i))
          end do
        end do

        call invers(acart_LS_inv,3)
        call invers(acart_SL_inv,3)
        call invers(acart_LL_inv,3)
        call invers(acart_St_inv,3)
        call invers(acart_t_inv,3)

        write(55 + 5*s%nAtoms+i,"(3(es16.9,2x))") e, e / aimag(acart_LS_inv(2,3)),  -1._dp * aimag(acart_LS_inv(2,2))/aimag(acart_LS_inv(2,3))
        write(55 + 6*s%nAtoms+i,"(3(es16.9,2x))") e, e / aimag(acart_SL_inv(2,3)),  -1._dp * aimag(acart_SL_inv(2,2))/aimag(acart_SL_inv(2,3))
        write(55 + 7*s%nAtoms+i,"(3(es16.9,2x))") e, e / aimag(acart_LL_inv(2,3)),  -1._dp * aimag(acart_LL_inv(2,2))/aimag(acart_LL_inv(2,3))
        write(55 + 8*s%nAtoms+i,"(3(es16.9,2x))") e, e / aimag(acart_St_inv(2,3)),  -1._dp * aimag(acart_St_inv(2,2))/aimag(acart_St_inv(2,3))
        write(55 + 9*s%nAtoms+i,"(3(es16.9,2x))") e, e / aimag(acart_t_inv(2,3)),  -1._dp * aimag(acart_t_inv(2,2))/aimag(acart_t_inv(2,3))
      end if
    end do

    call close_alpha_files()

    call deallocate_alpha()

  end subroutine write_alpha

end module mod_alpha
