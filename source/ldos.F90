! This subroutine calculates LDOS
subroutine ldos()
  use mod_f90_kind
  use mod_parameters
  use mod_mpi_pars, only: mpitag
  implicit none
  character(len=400) :: varm
  character(len=50)  :: fieldpart,socpart
  character(len=1)   :: SOCc
  integer            :: i,iw,j,mu
  real(double)       :: e
  real(double),dimension(:,:),allocatable :: ldosu,ldosd

  allocate(ldosu(Npl,9),ldosd(Npl,9))

  write(outputunit_loop,"('CALCULATING LDOS')")

  fieldpart = ""
  socpart   = ""
  if(SOC) then
    if(llinearsoc) then
      SOCc = "L"
    else
      SOCc = "T"
    end if
    write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
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

  ! Opening files
  ! LDOS
  do i=1,Npl
    iw = 1000+(mpitag-1)*Npl*2 + (i-1)*2 + 1
    write(varm,"('./results/',a1,'SOC/',a,'/LDOS/ldosu_layer',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),i,ncp,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=iw, file=varm,status='unknown')
    iw = iw+1
    write(varm,"('./results/',a1,'SOC/',a,'/LDOS/ldosd_layer',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),i,ncp,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=iw, file=varm,status='unknown')
  end do

  ldos_energy_loop: do count=1,npt1
    e = emin + (count-1)*deltae
    write(outputunit_loop,"('[ldos] ',i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e

    call ldos_jij_energy(e,ldosu,ldosd)

    ! Writing into files
    ! LDOS
    ldos_writing_plane_loop: do i=1,Npl
        iw = 1000+(mpitag-1)*Npl*2 + (i-1)*2 + 1
        write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
        iw = iw+1
        write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
    end do ldos_writing_plane_loop

  end do ldos_energy_loop

  ! Closing files
  do i=1,Npl
    iw = 1000+(mpitag-1)*Npl*2 + (i-1)*2 + 1
    close (iw)
    iw = iw+1
    close (iw)
  end do

  deallocate(ldosu,ldosd)

  return
end subroutine ldos