module mod_torques
  use mod_f90_kind
  implicit none
  ! Torques (sot,xc-torque,external torque ; x,y,z ; layers)
  integer       :: ntypetorque=3 ! Number of types of torques implemented
  complex(double),allocatable   :: torques(:,:,:),rtorques(:,:,:)
contains

  ! This subroutine allocates variables related to the torques calculation
  subroutine allocate_torques()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters
    implicit none
    integer           :: AllocateStatus

    ! Turning off external torque if static field is zero
    if((lfield).or.(total_hw_npt1/=1)) then
      ntypetorque=3
    else
      ntypetorque=2
    end if

    if(myrank_row==0) then
      allocate( torques(ntypetorque,3,Npl), STAT = AllocateStatus )
      if (AllocateStatus/=0) then
        write(outputunit,"('[allocate_torques] Not enough memory for: torques')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      if(renorm) then
        allocate( rtorques(ntypetorque,3,Npl), STAT = AllocateStatus )
        if (AllocateStatus/=0) then
          write(outputunit,"('[allocate_torques] Not enough memory for: rtorques')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if

    return
  end subroutine allocate_torques

  ! This subroutine deallocates variables related to the torques calculation
  subroutine deallocate_torques()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters, only: renorm
    implicit none

    if(myrank_row==0) then
      deallocate(torques)
      if(renorm) deallocate(rtorques)
    end if

    return
  end subroutine deallocate_torques

  ! This subroutine opens and closes all the files needed for the disturbances
  subroutine openclose_torque_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    use mod_system, only: nkpt
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=6)    :: folder
    character(len=3)    :: filename(ntypetorque)
    character(len=1)    :: direction(3),SOCc
    integer :: i,sigma,typetorque,iw,iflag,err,errt=0

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

    folder = "SOT"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename(1) = "SOT"
    filename(2) = "XCT"
    if((lfield).or.(total_hw_npt1/=1)) filename(3) = "EXT"

    direction(1) = "x"
    direction(2) = "y"
    direction(3) = "z"

    if(iflag==0) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 9000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
        close(unit=iw)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/r',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
          close(unit=iw)
        end if
      end do ; end do ; end do
    else if (iflag==1) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 9000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/r',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do ; end do
      ! Stop if some file does not exist
      if(errt/=0) then
        write(outputunit,"(a,i0,a)") "[openclose_torque_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 9000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do ; end do
    end if

    return
  end subroutine openclose_torque_files

  ! This subroutine write all the torques into files
  ! (already opened with openclose_torque_files(1))
  ! Some information may also be written on the screen
  subroutine write_torques(e)
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm
    use mod_magnet, only: mvec_spherical
    implicit none
    integer  :: i,iw,sigma,typetorque
    real(double) :: phase,sine,cosine
    real(double),intent(in) :: e

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
      iw = 9000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i

      if(abs(torques(typetorque,sigma,i))>=1.d-10) then
        phase  = atan2(aimag(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i)))
        sine   = real(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
        cosine = aimag(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
      else
        phase  = 0.d0
        sine   = 0.d0
        cosine = 0.d0
      end if

      write(unit=iw,fmt="(7(es16.9,2x))") e , abs(torques(typetorque,sigma,i)) , real(torques(typetorque,sigma,i)) , aimag(torques(typetorque,sigma,i)) , phase , sine , cosine

      ! Writing renormalized torques
      if(renorm) then
        iw = iw+1000

        if(abs(rtorques(typetorque,sigma,i))>=1.d-10) then
          phase  = atan2(aimag(rtorques(typetorque,sigma,i)),real(rtorques(typetorque,sigma,i)))
          sine   = real(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
          cosine = aimag(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
        else
          phase  = 0.d0
          sine   = 0.d0
          cosine = 0.d0
        end if

        write(unit=iw,fmt="(7(es16.9,2x))") e , abs(rtorques(typetorque,sigma,i)) , real(rtorques(typetorque,sigma,i)) , aimag(rtorques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)
      end if
    end do ; end do ; end do

    return
  end subroutine write_torques

  ! This subroutine opens and closes all the files needed for the disturbances
  subroutine openclose_dc_torque_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    use mod_system, only:nkpt
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=6)    :: folder
    character(len=3)    :: filename(ntypetorque)
    character(len=1)    :: direction(3),SOCc
    integer :: i,sigma,typetorque,iw,iflag,err,errt=0

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

    if(dcfield_dependence/=7) then
      if((dcfield_dependence/=1).and.(dcfield_dependence/=4).and.(dcfield_dependence/=5)) write(fieldpart,"(a,'_hwa=',es9.2)") trim(fieldpart),hwa
      if((dcfield_dependence/=2).and.(dcfield_dependence/=4).and.(dcfield_dependence/=6)) write(fieldpart,"(a,'_hwt=',f5.2)") trim(fieldpart),hwt
      if((dcfield_dependence/=3).and.(dcfield_dependence/=5).and.(dcfield_dependence/=6)) write(fieldpart,"(a,'_hwp=',f5.2)") trim(fieldpart),hwp
    end if
    if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
    if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
    if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
    if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"

    folder = "SOT"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename(1) = "SOT"
    filename(2) = "XCT"
    if((lfield).or.(total_hw_npt1/=1)) filename(3) = "EXT"

    direction(1) = "x"
    direction(2) = "y"
    direction(3) = "z"

    if(iflag==0) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 90000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
        close(unit=iw)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'r',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
          close(unit=iw)
        end if
      end do ; end do ; end do
    else if (iflag==1) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 90000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'r',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do ; end do
      ! Stop if some file does not exist
      if(errt/=0) then
        write(outputunit,"(a,i0,a)") "[openclose_dc_torque_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 90000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do ; end do
    end if

    return
  end subroutine openclose_dc_torque_files

  ! This subroutine write all the torques into files
  ! (already opened with openclose_torque_files(1))
  ! Some information may also be written on the screen
  subroutine write_dc_torques()
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm,dc_fields,hw_count
    use mod_magnet, only: mvec_spherical
    implicit none
    integer      :: i,iw,sigma,typetorque
    real(double) :: phase,sine,cosine

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
      iw = 90000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i

      if(abs(torques(typetorque,sigma,i))>=1.d-10) then
        phase  = atan2(aimag(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i)))
        sine   = real(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
        cosine = aimag(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
      else
        phase  = 0.d0
        sine   = 0.d0
        cosine = 0.d0
      end if
      write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(torques(typetorque,sigma,i)) , real(torques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)

      ! Writing renormalized torques
      if(renorm) then
        iw = iw+1000

        if(abs(rtorques(typetorque,sigma,i))>=1.d-10) then
          phase  = atan2(aimag(rtorques(typetorque,sigma,i)),real(rtorques(typetorque,sigma,i)))
          sine   = real(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
          cosine = aimag(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
        else
          phase  = 0.d0
          sine   = 0.d0
          cosine = 0.d0
        end if

        write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rtorques(typetorque,sigma,i)) , real(rtorques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)
      end if
    end do ; end do ; end do

    return
  end subroutine write_dc_torques

  ! This subroutine sorts torque files
  subroutine sort_torques()
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm,itype
    use mod_tools, only: sort_file
    implicit none
    integer  :: i,iw,sigma,typetorque,idc=1

    ! Opening torque files
    if(itype==9) then
      idc=10
      call openclose_dc_torque_files(1)
    else
      call openclose_torque_files(1)
    end if

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
      iw = 9000*idc+(typetorque-1)*Npl*3+(sigma-1)*Npl+i

      ! Sorting torque files
      call sort_file(iw,.true.)

      ! Sorting renormalized torque files
      if(renorm) then
        iw = iw+1000
        call sort_file(iw,.true.)
      end if
    end do ; end do ; end do

    ! Closing torque files
    if(itype==9) then
      call openclose_dc_torque_files(2)
    else
      call openclose_torque_files(2)
    end if

    return
  end subroutine sort_torques

end module mod_torques
