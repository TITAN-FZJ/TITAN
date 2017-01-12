module mod_sha
  use mod_f90_kind
  implicit none
  ! Spin Hall Angle
  real(double), allocatable :: sha_re(:,:)
  real(double) :: sha_re_total(4)
  complex(double), allocatable :: sha_complex(:,:)
  complex(double) :: sha_complex_total(4)
contains

  ! This subroutine allocates variables related to the sha calculation
  subroutine allocate_sha()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters, only: Npl,outputunit
    implicit none
    integer           :: AllocateStatus

    if(myrank_row.eq.0) then
      allocate( sha_re(4,Npl),sha_complex(4,Npl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
         write(outputunit,"('[allocate_currents] Not enough memory for: sha_re,sha_complex')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if

    return
  end subroutine allocate_sha

  ! This subroutine allocates variables related to the sha calculation
  subroutine deallocate_sha()
    use mod_f90_kind
    use mod_mpi_pars
    implicit none

    if(myrank_row.eq.0) deallocate(sha_re,sha_complex)

    return
  end subroutine deallocate_sha

  ! This subroutine opens and closes all the files needed for the sha
  subroutine openclose_sha_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=1)    :: SOCc
    character(len=6)    :: folder(8)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,iflag,err,errt=0

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

    folder(1) = "CC"
    folder(2) = "SC"
    folder(3) = "SC"
    folder(4) = "SC"
    folder(5) = "LC"
    folder(6) = "LC"
    folder(7) = "LC"
    folder(8) = "SHA"
    if(lhfresponses) then
      do i=1,8
        folder(i) = trim(folder(i)) // "_HF"
      end do
    end if

    filename(1) = "Ich"
    filename(2) = "Isx"
    filename(3) = "Isy"
    filename(4) = "Isz"
    filename(5) = "Ilx"
    filename(6) = "Ily"
    filename(7) = "Ilz"

    if(iflag.eq.0) then
      ! Header for SHA
      do j=2,4
        do i=1,Npl
          iw = 8200+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),filename(j),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#     energy     ,  SHA reIs/reIc  ,  SHA abs(Is/Ic) , SHA atan(Is/Ic) ,  SHA re(Is/Ic)  ,  SHA im(Is/Ic)  ')")
          close(unit=iw)
        end do
        iw = 8200+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),filename(j),parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#     energy     ,  SHA reIs/reIc  ,  SHA abs(Is/Ic) , SHA atan(Is/Ic) ,  SHA re(Is/Ic)  ,  SHA im(Is/Ic)  ')")
        close(unit=iw)
      end do

    else if(iflag.eq.1) then
      ! SHA
      do j=2,4
        do i=1,Npl
          iw = 8200+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),filename(j),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end do
        iw = 8200+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),filename(j),parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do

      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(outputunit,"(a,i0,a)") "[openclose_sha_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      ! Closing SHA
      do j=2,4
        do i=1,Npl
          iw = 8200+(i-1)*3+j-1
          close(unit=iw)
        end do
        iw = 8200+Npl*3+j-1
        close(unit=iw)
      end do
    end if

    return
  end subroutine openclose_sha_files

  ! This subroutine write all the SHA into files
  ! (already opened with openclose_sha_files(1))
  subroutine write_sha(e)
    use mod_f90_kind
    use mod_parameters, only: Npl
    implicit none
    integer  :: i,j,iw
    real(double),intent(in) :: e

    do j=2,4
      do i=1,Npl
        iw = 8200+(i-1)*3+j-1
        write(unit=iw,fmt="(6(es16.9,2x))") e , sha_re(j,i) , abs(sha_complex(j,i)) , atan2(aimag(sha_complex(j,i)),real(sha_complex(j,i))) , real(sha_complex(j,i)) , aimag(sha_complex(j,i))
      end do
      iw = 8200+Npl*3+j-1
      write(unit=iw,fmt="(6(es16.9,2x))") e , sha_re_total(j) , abs(sha_complex_total(j)) , atan2(aimag(sha_complex_total(j)),real(sha_complex_total(j))) , real(sha_complex_total(j)) , aimag(sha_complex_total(j))
    end do

    return
  end subroutine write_sha

  ! This subroutine opens and closes all the files needed for the sha
  subroutine openclose_dc_sha_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=1)    :: SOCc
    character(len=6)    :: folder(8)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,iflag,err,errt=0

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

    if(dcfield_dependence.ne.7) then
      if((dcfield_dependence.ne.1).and.(dcfield_dependence.ne.4).and.(dcfield_dependence.ne.5)) write(fieldpart,"(a,'_hwa=',es9.2)") trim(fieldpart),hwa
      if((dcfield_dependence.ne.2).and.(dcfield_dependence.ne.4).and.(dcfield_dependence.ne.6)) write(fieldpart,"(a,'_hwt=',f5.2)") trim(fieldpart),hwt
      if((dcfield_dependence.ne.3).and.(dcfield_dependence.ne.5).and.(dcfield_dependence.ne.6)) write(fieldpart,"(a,'_hwp=',f5.2)") trim(fieldpart),hwp
    end if
    if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
    if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
    if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
    if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"

    folder(1) = "CC"
    folder(2) = "SC"
    folder(3) = "SC"
    folder(4) = "SC"
    folder(5) = "LC"
    folder(6) = "LC"
    folder(7) = "LC"
    folder(8) = "SHA"
    if(lhfresponses) then
      do i=1,8
        folder(i) = trim(folder(i)) // "_HF"
      end do
    end if

    filename(1) = "Ich"
    filename(2) = "Isx"
    filename(3) = "Isy"
    filename(4) = "Isz"
    filename(5) = "Ilx"
    filename(6) = "Ily"
    filename(7) = "Ilz"

    if(iflag.eq.0) then
      ! Header for SHA
      do j=2,4
        do i=1,Npl
          iw = 82000+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#',a,'    SHA reIs/reIc   ,    SHA abs(Is/Ic)   ,    SHA atan(Is/Ic)   ,    SHA re(Is/Ic)   ,    SHA im(Is/Ic)   ')") trim(dc_header)
          close(unit=iw)
        end do
        iw = 82000+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#',a,'    SHA reIs/reIc   ,    SHA abs(Is/Ic)   ,    SHA atan(Is/Ic)   ,    SHA re(Is/Ic)   ,    SHA im(Is/Ic)   ')") trim(dc_header)
        close(unit=iw)
      end do

    else if(iflag.eq.1) then
      ! SHA
      do j=2,4
        do i=1,Npl
          iw = 82000+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          if(.not.lsha) errt = errt + err
        end do
        iw = 82000+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(8)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        if(.not.lsha) errt = errt + err
      end do

      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(outputunit,"(a,i0,a)") "[openclose_dc_sha_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      ! SHA
      do j=2,4
        do i=1,Npl
          iw = 82000+(i-1)*3+j-1
          close(unit=iw)
        end do
        iw = 82000+Npl*3+j-1
        close(unit=iw)
      end do

    end if

    return
  end subroutine openclose_dc_sha_files

  ! This subroutine write all the sha in the dc limit into files
  ! (already opened with openclose_dc_sha_files(1))
  ! Some information is also written on the screen
  subroutine write_dc_sha()
    use mod_f90_kind
    use mod_parameters
    use mod_magnet, only: mvec_spherical
    implicit none
    integer  :: i,j,iw

    do j=2,4
      do i=1,Npl
        iw = 82000+(i-1)*3+j-1
        write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , sha_re(j,i) , abs(sha_complex(j,i)) , atan2(aimag(sha_complex(j,i)),real(sha_complex(j,i))) , real(sha_complex(j,i)) , aimag(sha_complex(j,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      end do
      iw = 82000+Npl*3+j-1
      write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , sha_re_total(j) , abs(sha_complex_total(j)) , atan2(aimag(sha_complex_total(j)),real(sha_complex_total(j))) , real(sha_complex_total(j)) , aimag(sha_complex_total(j)) , mvec_spherical(i,2) , mvec_spherical(i,3)
    end do

    return
  end subroutine write_dc_sha

  ! This subroutine sorts SHA files
  subroutine sort_sha()
    use mod_f90_kind
    use mod_parameters, only: Npl,itype
    use mod_tools, only: sort_file
    implicit none
    integer  :: i,j,iw,idc=1

    ! Opening SHA files
    if(itype.eq.9) then
      idc=10
      call openclose_dc_sha_files(1)
    else
      call openclose_sha_files(1)
    end if

    ! SHA
    do j=2,4
      do i=1,Npl
        iw = 8200*idc+(i-1)*3+j-1
        call sort_file(iw,.true.)
      end do
      iw = 8200*idc+Npl*3+j-1
      call sort_file(iw,.true.)
    end do

    ! Closing SHA files
    if(itype.eq.9) then
      call openclose_dc_sha_files(2)
    else
      call openclose_sha_files(2)
    end if

    return
  end subroutine sort_sha

end module mod_sha