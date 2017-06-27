module mod_alpha
  use mod_f90_kind, only: double
  implicit none

  complex(double), dimension(:,:),   allocatable :: m_chi, m_chi_hf, m_chi_inv, m_chi_hf_inv

contains
  subroutine open_alpha_files()
    use mod_parameters
    use mod_system, only:nkpt
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=1)    :: SOCc
    integer :: i,j,sigma,iw,iflag,err,errt=0

    fieldpart = ""
    socpart   = ""
    if(SOC) then
      if(llinearsoc) then
        SOCc = "L"
      else
        SOCc = "T"
      end if
      if(abs(socscale-1.d0)>1.d-6) write(socpart,"('_socscale=',f5.2)") socscale
    else
      SOCc = "F"
    end if
    if(lfield) then
      write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
      if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
      if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
      if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
      if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"
    end if

    do i=1,Npl
      write(varm,"('./results/',a1,'SOC/',a,'/A/chi_',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
      open (unit=55+i, file=varm, status='replace', form='formatted')
      write(unit=55+i, fmt="('#     energy    ,  alpha   ,  gamma  ,  alpha/gamma  ')")
      write(varm,"('./results/',a1,'SOC/',a,'/A/chi_hf',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
      open (unit=55+Npl+i, file=varm, status='replace', form='formatted')
      write(unit=55+Npl+i, fmt="('#     energy    ,  alpha   ,  gamma  ,  alpha/gamma  ')")
      write(varm,"('./results/',a1,'SOC/',a,'/A/chi_inv',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
      open (unit=55+2*Npl+i, file=varm, status='replace', form='formatted')
      write(unit=55+2*Npl+i, fmt="('#     energy    ,  alpha   ,  gamma  ,  alpha/gamma  ')")
      write(varm,"('./results/',a1,'SOC/',a,'/A/chi_hf_inv',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
      open (unit=55+3*Npl+i, file=varm, status='replace', form='formatted')
      write(unit=55+3*Npl+i, fmt="('#     energy    ,  alpha   ,  gamma  ,  alpha/gamma  ')")
    end do

    allocate(m_chi(4*Npl,4*Npl), m_chi_hf(4*Npl,4*Npl), m_chi_inv(4*Npl,4*Npl), m_chi_hf_inv(4*Npl,4*Npl))

    return
  end subroutine open_alpha_files

  subroutine close_alpha_files()
    use mod_parameters, only: Npl
    implicit none
    integer :: i

    do i = 1, Npl
      close(unit=55+i)
    end do

    deallocate(m_chi, m_chi_hf, m_chi_inv, m_chi_hf_inv)

    return
  end subroutine close_alpha_files

  subroutine write_alpha(e)
    use mod_f90_kind, only: double
    use mod_constants, only: zi
    use mod_susceptibilities, only: schi, schihf
    use mod_parameters, only: sigmai2i, Npl
    use mod_magnet, only: mabs
    implicit none
    real(double) :: e
    complex(double) :: axx, axy, axx_hf, axy_hf, axx_inv, axy_inv, axx_hf_inv, axy_hf_inv
    integer :: i,j,sigma, sigmap
    do i = 1, Npl
      do j = 1, Npl
        do sigma = 1, 4
          do sigmap = 1, 4
            m_chi   (sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schi  (sigma,sigmap,i,j)
            m_chi_hf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schihf(sigma,sigmap,i,j)
          end do
        end do
      end do
    end do
    m_chi_inv    = m_chi
    m_chi_hf_inv = m_chi_hf
    call invers(m_chi_inv, 4*Npl)
    call invers(m_chi_hf_inv, 4*Npl)
    do i = 1, Npl
      axx_inv    = 0.5d0   *(m_chi_inv(sigmai2i(1,i),sigmai2i(4,i)) + m_chi_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_inv(sigmai2i(4,i),sigmai2i(4,i)) + m_chi_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axy_inv    = 0.5d0*zi*(m_chi_inv(sigmai2i(1,i),sigmai2i(4,i)) - m_chi_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_inv(sigmai2i(4,i),sigmai2i(4,i)) - m_chi_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axx_hf_inv = 0.5d0   *(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(4,i)) + m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf_inv(sigmai2i(4,i),sigmai2i(4,i)) + m_chi_hf_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axy_hf_inv = 0.5d0*zi*(m_chi_hf_inv(sigmai2i(1,i),sigmai2i(4,i)) - m_chi_hf_inv(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf_inv(sigmai2i(4,i),sigmai2i(4,i)) - m_chi_hf_inv(sigmai2i(4,i),sigmai2i(1,i)))
      axx        = 0.5d0   *(m_chi(sigmai2i(1,i),sigmai2i(4,i)) + m_chi(sigmai2i(1,i),sigmai2i(1,i)) + m_chi(sigmai2i(4,i),sigmai2i(4,i)) + m_chi(sigmai2i(4,i),sigmai2i(1,i)))
      axy        = 0.5d0*zi*(m_chi(sigmai2i(1,i),sigmai2i(4,i)) - m_chi(sigmai2i(1,i),sigmai2i(1,i)) + m_chi(sigmai2i(4,i),sigmai2i(4,i)) - m_chi(sigmai2i(4,i),sigmai2i(1,i)))
      axx_hf     = 0.5d0   *(m_chi_hf(sigmai2i(1,i),sigmai2i(4,i)) + m_chi_hf(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf(sigmai2i(4,i),sigmai2i(4,i)) + m_chi_hf(sigmai2i(4,i),sigmai2i(1,i)))
      axy_hf     = 0.5d0*zi*(m_chi_hf(sigmai2i(1,i),sigmai2i(4,i)) - m_chi_hf(sigmai2i(1,i),sigmai2i(1,i)) + m_chi_hf(sigmai2i(4,i),sigmai2i(4,i)) - m_chi_hf(sigmai2i(4,i),sigmai2i(1,i)))

      write(55 +       i,*) e, -1.d0 * aimag(axx       )/aimag(axy       ), e / (mabs(i) * aimag(axy       )), -e*mabs(i)*aimag(axx       )
      write(55 +   Npl+i,*) e, -1.d0 * aimag(axx_hf    )/aimag(axy_hf    ), e / (mabs(i) * aimag(axy_hf    )), -e*mabs(i)*aimag(axx_hf    )
      write(55 + 2*Npl+i,*) e, -1.d0 * aimag(axx_inv   )/aimag(axy_inv   ), e / (mabs(i) * aimag(axy_inv   )), -e*mabs(i)*aimag(axx_inv   )
      write(55 + 3*Npl+i,*) e, -1.d0 * aimag(axx_hf_inv)/aimag(axy_hf_inv), e / (mabs(i) * aimag(axy_hf_inv)), -e*mabs(i)*aimag(axx_hf_inv)

    end do

  end subroutine write_alpha




end module mod_alpha
