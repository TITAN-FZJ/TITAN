module mod_torques
  use mod_f90_kind
  implicit none
  ! Torques (sot,xc-torque,external torque ; x,y,z ; layers)
  integer       :: ntypetorque=3 ! Number of types of torques implemented
  complex(double),allocatable   :: torques(:,:,:),rtorques(:,:,:)
contains

  ! This subroutine allocates variables related to the disturbance calculation
  subroutine allocate_torques()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters, only: Npl,renorm,outputunit,lfield
    implicit none
    integer           :: AllocateStatus

    ! Turning off external torque if static field is zero
    if(.not.lfield) then
      ntypetorque=2
    else
      ntypetorque=3
    end if

    if(myrank_row.eq.0) then
      allocate( torques(ntypetorque,3,Npl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(outputunit,"('[allocate_torques] Not enough memory for: torques')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      if(renorm) then
        allocate( rtorques(ntypetorque,3,Npl), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          write(outputunit,"('[allocate_torques] Not enough memory for: rtorques')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if

    return
  end subroutine allocate_torques

  ! This subroutine allocates variables related to the disturbance calculation
  subroutine deallocate_torques()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters, only: renorm
    implicit none

    if(myrank_row.eq.0) then
      deallocate(torques)
      if(renorm) deallocate(rtorques)
    end if

    return
  end subroutine deallocate_torques

  ! This subroutine opens and closes all the files needed for the disturbances
  subroutine openclose_torque_files(iflag)
    use mod_parameters
    use mod_mpi_pars
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
      write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
    else
      SOCc = "F"
    end if
    if(lfield) then
      write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
      if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
    end if

    folder = "SOT"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename(1) = "SOT"
    filename(2) = "XCT"
    if(lfield) filename(3) = "EXT"

    direction(1) = "x"
    direction(2) = "y"
    direction(3) = "z"

    if(iflag.eq.0) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 9000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   energy  , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
        close(unit=iw)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/r',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#   energy  , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
          close(unit=iw)
        end if
      end do ; end do ; end do
    else if (iflag.eq.1) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 9000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/r',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(filename(typetorque)),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        end if
      end do ; end do ; end do
      ! Stop if some file does not exist
      if(errt.ne.0) then
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
    implicit none
    integer  :: i,iw,sigma,typetorque
    real(double),intent(in) :: e

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
      iw = 9000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i

      write(unit=iw,fmt="(7(es16.9,2x))") e,abs(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i)),aimag(torques(typetorque,sigma,i)),atan2(aimag(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i))),real(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i)),aimag(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))

      ! Writing renormalized torques
      if(renorm) then
        iw = iw+1000
        write(unit=iw,fmt="(7(es16.9,2x))") e,abs(rtorques(typetorque,sigma,i)),real(rtorques(typetorque,sigma,i)),aimag(rtorques(typetorque,sigma,i)),atan2(aimag(rtorques(typetorque,sigma,i)),real(rtorques(typetorque,sigma,i))),real(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i)),aimag(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
      end if
    end do ; end do ; end do

    return
  end subroutine write_torques

  ! This subroutine opens and closes all the files needed for the disturbances
  subroutine openclose_dc_torque_files(iflag)
    use mod_parameters
    use mod_mpi_pars
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
      write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
    else
      SOCc = "F"
    end if

    if(dcfield_dependence.ne.7) then
      if((dcfield_dependence.ne.1).and.(dcfield_dependence.ne.4).and.(dcfield_dependence.ne.5)) write(fieldpart,"(a,'_hwa=',es9.2)") trim(fieldpart),hwa
      if((dcfield_dependence.ne.2).and.(dcfield_dependence.ne.4).and.(dcfield_dependence.ne.6)) write(fieldpart,"(a,'_hwt=',f5.2)") trim(fieldpart),hwt
      if((dcfield_dependence.ne.3).and.(dcfield_dependence.ne.5).and.(dcfield_dependence.ne.6)) write(fieldpart,"(a,'_hwp=',f5.2)") trim(fieldpart),hwp
    end if
    if(ltesla) fieldpart = trim(fieldpart) // "_tesla"

    folder = "SOT"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename(1) = "SOT"
    filename(2) = "XCT"
    if(lfield) filename(3) = "EXT"

    direction(1) = "x"
    direction(2) = "y"
    direction(3) = "z"

    if(iflag.eq.0) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 90000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/',a,a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(dcprefix),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
        close(unit=iw)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/',a,'r',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(dcprefix),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
          close(unit=iw)
        end if
      end do ; end do ; end do
    else if (iflag.eq.1) then
      do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
        iw = 90000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/',a,a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(dcprefix),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/',a,'/',a,'r',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),trim(dcprefix),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        end if
      end do ; end do ; end do
      ! Stop if some file does not exist
      if(errt.ne.0) then
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
    use mod_magnet, only: mtheta,mphi
    implicit none
    integer  :: i,iw,sigma,typetorque

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,Npl
      iw = 90000+(typetorque-1)*Npl*3+(sigma-1)*Npl+i

      write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(torques(typetorque,sigma,i)) , real(torques(typetorque,sigma,i)) , atan2(aimag(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i))) , real(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i)) , aimag(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i)) , mtheta(i) , mphi(i)

      ! Writing renormalized torques
      if(renorm) then
        iw = iw+1000
        write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rtorques(typetorque,sigma,i)) , real(rtorques(typetorque,sigma,i)) , atan2(aimag(rtorques(typetorque,sigma,i)),real(rtorques(typetorque,sigma,i))) , real(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i)) , aimag(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i)) , mtheta(i) , mphi(i)
      end if
    end do ; end do ; end do

    return
  end subroutine write_dc_torques

end module mod_torques