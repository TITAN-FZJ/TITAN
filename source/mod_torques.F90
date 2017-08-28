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
    use mod_magnet, only: lfield, total_hw_npt1
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

  subroutine create_torque_files()
  !! This subroutine creates all the files needed for the disturbances
  use mod_parameters, only: fieldpart, lhfresponses, Npl_folder, eta, Utype, suffix, renorm, missing_files
  use mod_SOC, only: SOCc, socpart
  use mod_magnet, only: lfield, total_hw_npt1
  use mod_mpi_pars
  use mod_system, only: s => sys
  use electricfield, only: strElectricField
  use EnergyIntegration, only: strEnergyParts
  implicit none

  character(len=500)  :: varm
  character(len=6)    :: folder
  character(len=3)    :: filename(ntypetorque)
  character(len=1)    :: direction(3)
  integer :: i,sigma,typetorque,iw,iflag,err,errt=0

  folder = "SOT"
  if(lhfresponses) folder = trim(folder) // "_HF"
  filename(1) = "SOT"
  filename(2) = "XCT"
  if((lfield).or.(total_hw_npt1/=1)) filename(3) = "EXT"

  direction(1) = "x"
  direction(2) = "y"
  direction(3) = "z"

  do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
    iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
    write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
    open (unit=iw, file=varm, status='replace', form='formatted')
    write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
    close(unit=iw)
    if(renorm) then
      iw = iw+1000
      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/r',a,a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
      open (unit=iw, file=varm, status='replace', form='formatted')
      write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
      close(unit=iw)
    end if
  end do ; end do ; end do

  return

  end subroutine create_torque_files

  subroutine open_torque_files()
  !! This subroutine opens all the files needed for the disturbances
  use mod_parameters, only: fieldpart, lhfresponses, Npl_folder, eta, Utype, suffix, renorm, missing_files
  use mod_SOC, only: SOCc, socpart
  use mod_magnet, only: lfield, total_hw_npt1
  use mod_mpi_pars
  use mod_system, only: s => sys
  use electricfield, only: strElectricField
  use EnergyIntegration, only: strEnergyParts
  implicit none

  character(len=500)  :: varm
  character(len=6)    :: folder
  character(len=3)    :: filename(ntypetorque)
  character(len=1)    :: direction(3)
  integer :: i,sigma,typetorque,iw,iflag,err,errt=0

  folder = "SOT"
  if(lhfresponses) folder = trim(folder) // "_HF"
  filename(1) = "SOT"
  filename(2) = "XCT"
  if((lfield).or.(total_hw_npt1/=1)) filename(3) = "EXT"

  direction(1) = "x"
  direction(2) = "y"
  direction(3) = "z"

  do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
    iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
    write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
    open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
    errt = errt + err
    if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
    if(renorm) then
      iw = iw+1000
      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/r',a,a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename(typetorque)),direction(sigma),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
      open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
    end if
  end do ; end do ; end do
  ! Stop if some file does not exist
  if(errt/=0) call abortProgram("[openclose_torque_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))

  return

  end subroutine open_torque_files

  subroutine close_torque_files()
  !! This subroutine closes all the files needed for the disturbances
  use mod_parameters, only: fieldpart, lhfresponses, Npl_folder, eta, Utype, suffix, renorm, missing_files
  use mod_SOC, only: SOCc, socpart
  use mod_magnet, only: lfield, total_hw_npt1
  use mod_mpi_pars
  use mod_system, only: s => sys
  use electricfield, only: strElectricField
  use EnergyIntegration, only: strEnergyParts
  implicit none

  character(len=500)  :: varm
  character(len=6)    :: folder
  character(len=3)    :: filename(ntypetorque)
  character(len=1)    :: direction(3)
  integer :: i,sigma,typetorque,iw,iflag,err,errt=0

  do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
    iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
    close(unit=iw)

    if(renorm) then
      iw = iw+1000
      close(unit=iw)
    end if
  end do ; end do ; end do

  return

  end subroutine close_torque_files


  ! This subroutine write all the torques into files
  ! (already opened with openclose_torque_files(1))
  ! Some information may also be written on the screen
  subroutine write_torques(e)
    use mod_f90_kind
    use mod_parameters, only: renorm
    use mod_magnet, only: mvec_spherical
    use mod_System, only: s => sys
    implicit none
    integer  :: i,iw,sigma,typetorque
    real(double) :: phase,sine,cosine
    real(double),intent(in) :: e

    call open_torque_files()

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
      iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i

      if(abs(torques(typetorque,sigma,i))>=1.d-10) then
        phase  = atan2(aimag(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i)))
        sine   = real(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
        cosine = aimag(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
      else
        phase  = 0.d0
        sine   = 0.d0
        cosine = 0.d0
      end if

      write(unit=iw,fmt="(9(es16.9,2x))") e , abs(torques(typetorque,sigma,i)) , real(torques(typetorque,sigma,i)) , aimag(torques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)

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

        write(unit=iw,fmt="(9(es16.9,2x))") e , abs(rtorques(typetorque,sigma,i)) , real(rtorques(typetorque,sigma,i)) , aimag(rtorques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)
      end if
    end do ; end do ; end do

    call close_torque_files()

    return
  end subroutine write_torques


  subroutine create_dc_torque_files()
  !! This subroutine creates all the files needed for the disturbances
  use mod_parameters, only: dcfieldpart, lhfresponses, Npl_folder, count, missing_files, eta, Utype, suffix, renorm
  use mod_magnet, only: lfield, total_hw_npt1, dcprefix, dcfield_dependence, dcfield, dc_header
  use mod_mpi_pars
  use mod_SOC, only: SOCc, socpart
  use mod_system, only: s => sys
  use electricfield, only: strElectricField
  use EnergyIntegration, only: strEnergyParts
  implicit none

  character(len=500)  :: varm
  character(len=6)    :: folder
  character(len=3)    :: filename(ntypetorque)
  character(len=1)    :: direction(3)
  integer :: i,sigma,typetorque,iw,iflag,err,errt=0

  folder = "SOT"
  if(lhfresponses) folder = trim(folder) // "_HF"
  filename(1) = "SOT"
  filename(2) = "XCT"
  if((lfield).or.(total_hw_npt1/=1)) filename(3) = "EXT"

  direction(1) = "x"
  direction(2) = "y"
  direction(3) = "z"

  do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
    iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
    write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
    open (unit=iw, file=varm, status='replace', form='formatted')
    write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
    close(unit=iw)
    if(renorm) then
      iw = iw+1000
      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'r',a,a,'_',a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
      open (unit=iw, file=varm, status='replace', form='formatted')
      write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
      close(unit=iw)
    end if
  end do ; end do ; end do

  return

  end subroutine create_dc_torque_files

  subroutine open_dc_torque_files()
  !! This subroutine opens all the files needed for the disturbances
  use mod_parameters, only: dcfieldpart, lhfresponses, Npl_folder, count, missing_files, eta, Utype, suffix, renorm
  use mod_magnet, only: lfield, total_hw_npt1, dcprefix, dcfield_dependence, dcfield, dc_header
  use mod_mpi_pars
  use mod_SOC, only: SOCc, socpart
  use mod_system, only: s => sys
  use electricfield, only: strElectricField
  use EnergyIntegration, only: strEnergyParts
  implicit none

  character(len=500)  :: varm
  character(len=6)    :: folder
  character(len=3)    :: filename(ntypetorque)
  character(len=1)    :: direction(3)
  integer :: i,sigma,typetorque,iw,iflag,err,errt=0

  folder = "SOT"
  if(lhfresponses) folder = trim(folder) // "_HF"
  filename(1) = "SOT"
  filename(2) = "XCT"
  if((lfield).or.(total_hw_npt1/=1)) filename(3) = "EXT"

  direction(1) = "x"
  direction(2) = "y"
  direction(3) = "z"

  do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
    iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
    write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
    open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
    errt = errt + err
    if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
    if(renorm) then
      iw = iw+1000
      write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'r',a,a,'_',a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
      open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
    end if
  end do ; end do ; end do
  ! Stop if some file does not exist
  if(errt/=0) call abortProgram("[openclose_dc_torque_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))

  return

  end subroutine open_dc_torque_files

  subroutine close_dc_torque_files()
  !! This subroutine closes all the files needed for the disturbances
  use mod_parameters, only: dcfieldpart, lhfresponses, Npl_folder, count, missing_files, eta, Utype, suffix, renorm
  use mod_magnet, only: lfield, total_hw_npt1, dcprefix, dcfield_dependence, dcfield, dc_header
  use mod_mpi_pars
  use mod_SOC, only: SOCc, socpart
  use mod_system, only: s => sys
  use electricfield, only: strElectricField
  use EnergyIntegration, only: strEnergyParts
  implicit none

  character(len=500)  :: varm
  character(len=6)    :: folder
  character(len=3)    :: filename(ntypetorque)
  character(len=1)    :: direction(3)
  integer :: i,sigma,typetorque,iw,iflag,err,errt=0

  do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
    iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
    close(unit=iw)

    if(renorm) then
      iw = iw+1000
      close(unit=iw)
    end if
  end do ; end do ; end do

  return

  end subroutine close_dc_torque_files

  ! This subroutine write all the torques into files
  ! (already opened with openclose_torque_files(1))
  ! Some information may also be written on the screen
  subroutine write_dc_torques()
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm
    use mod_magnet, only: mvec_spherical,dc_fields,hw_count
    use mod_System, only: s => sys
    implicit none
    integer      :: i,iw,sigma,typetorque
    real(double) :: phase,sine,cosine

    call open_dc_torque_files()

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
      iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i

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

    call close_dc_torque_files()

    return
  end subroutine write_dc_torques

  ! This subroutine sorts torque files
  subroutine sort_torques()
    use mod_f90_kind
    use mod_parameters, only: renorm,itype
    use mod_tools, only: sort_file
    use mod_System, only: s => sys
    implicit none
    integer  :: i,iw,sigma,typetorque,idc=1

    ! Opening torque files
    if(itype==9) then
      idc=10
      call open_dc_torque_files()
    else
      call open_torque_files()
    end if

    do typetorque=1,ntypetorque ; do sigma=1,3 ; do i=1,s%nAtoms
      iw = 9000*idc+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i

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
      call close_dc_torque_files()
    else
      call close_torque_files()
    end if

    return
  end subroutine sort_torques

end module mod_torques
