module mod_torques
  use mod_kind, only: dp
  implicit none
  integer       :: ntypetorque=2 ! Number of types of torques implemented
  ! Torques (sot,xc-torque,external torque ; x,y,z ; layers)
  complex(dp),allocatable   :: torques(:,:,:),total_torques(:,:), rtorques(:,:,:)

  character(len=6), parameter, private :: folder = "SOT"
  character(len=3), dimension(3), parameter, private :: filename = ["SOT", "XCT", "EXT"]
  character(len=1), dimension(3), parameter, private    :: direction = ["x", "y", "z"]

contains

  subroutine allocate_torques()
  !> This subroutine allocates variables related to the torques calculation
    use mod_parameters, only: renorm
    use mod_System,     only: s => sys
    use mod_magnet,     only: lfield, total_hw_npt1
    use mod_mpi_pars,   only: rFreq,abortProgram
    implicit none
    integer :: AllocateStatus

    ! Turning off external torque if static field is cZero
    if((lfield).or.(total_hw_npt1/=1)) then
      ntypetorque=3
    else
      ntypetorque=2
    end if

    if(rFreq(1)==0) then
      allocate( torques(ntypetorque,3, s%nAtoms), total_torques(ntypetorque,3), STAT = AllocateStatus )
      if (AllocateStatus/=0) call abortProgram("[allocate_torques] Not enough memory for: torques,total_torques")
      if(renorm) then
        allocate( rtorques(ntypetorque,3,s%nAtoms), STAT = AllocateStatus )
        if (AllocateStatus/=0) call abortProgram("[allocate_torques] Not enough memory for: rtorques")
      end if
    end if

  end subroutine allocate_torques


  subroutine deallocate_torques()
  !> This subroutine deallocates variables related to the torques calculation
    implicit none

    if(allocated(torques)) deallocate(torques)
    if(allocated(total_torques)) deallocate(total_torques)
    if(allocated(rtorques)) deallocate(rtorques)

  end subroutine deallocate_torques


  subroutine create_torque_files()
  !> This subroutine creates all the files needed for the torques
    use mod_parameters, only: output, renorm
    use mod_system,     only: s => sys
    implicit none

    character(len=500)  :: varm
    integer :: i,sigma,typetorque,iw

    do typetorque = 1, ntypetorque
      do sigma = 1, 3
        do i = 1, s%nAtoms
          iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
          write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(filename(typetorque)),direction(sigma),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
          close(unit=iw)
          if(renorm) then
            iw = iw+1000
            write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/r',a,a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(filename(typetorque)),direction(sigma),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
            open (unit=iw, file=varm, status='replace', form='formatted')
            write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
            close(unit=iw)
          end if
        end do
        ! Total torque files
        iw = 9500+(typetorque-1)*3+sigma
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(filename(typetorque)),direction(sigma),trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
        close(unit=iw)
      end do
    end do
  end subroutine create_torque_files


  subroutine open_torque_files()
  !> This subroutine opens all the files needed for the torques
    use mod_parameters, only: output, renorm, missing_files
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none

    character(len=500)  :: varm
    integer :: i,sigma,typetorque,iw,err,errt=0

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
          iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
          write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(filename(typetorque)),direction(sigma),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
          if(renorm) then
            iw = iw+1000
            write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/r',a,a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(filename(typetorque)),direction(sigma),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
            open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
            errt = errt + err
            if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
          end if
        end do
        ! Total torque files
        iw = 9500+(typetorque-1)*3+sigma
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(filename(typetorque)),direction(sigma),trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
      end do
    end do
    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openclose_torque_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
  end subroutine open_torque_files


  subroutine close_torque_files()
  !> This subroutine closes all the files needed for the torques
    use mod_parameters, only: renorm
    use mod_System, only: s => sys
    implicit none

    integer :: i,sigma,typetorque,iw

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
          iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
          close(unit=iw)

          if(renorm) then
            iw = iw+1000
            close(unit=iw)
          end if
        end do
        ! Total torque files
        iw = 9500+(typetorque-1)*3+sigma
        close(unit=iw)
      end do
    end do

  end subroutine close_torque_files


  subroutine write_torques(e)
  !> This subroutine write all the torques into files
  !> (already opened with openclose_torque_files(1))
  !> Some information may also be written on the screen
    use mod_kind, only: dp
    use mod_parameters, only: renorm
    use mod_magnet,     only: mvec_spherical,mtotal_spherical
    use mod_System,     only: s => sys
    implicit none
    integer  :: i,iw,sigma,typetorque
    real(dp) :: phase,sine,cosine
    real(dp),intent(in) :: e

    call open_torque_files()

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
          iw = 9000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i

          if(abs(torques(typetorque,sigma,i))>=1.e-15_dp) then
            phase  = atan2(aimag(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i)))
            sine   = real(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
            cosine = aimag(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
          else
            phase  = 0._dp
            sine   = 0._dp
            cosine = 0._dp
          end if

          write(unit=iw,fmt="(9(es16.9,2x))") e , abs(torques(typetorque,sigma,i)) , real(torques(typetorque,sigma,i)) , aimag(torques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)

          ! Writing renormalized torques
          if(renorm) then
            iw = iw+1000

            if(abs(rtorques(typetorque,sigma,i))>=1.e-15_dp) then
              phase  = atan2(aimag(rtorques(typetorque,sigma,i)),real(rtorques(typetorque,sigma,i)))
              sine   = real(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
              cosine = aimag(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
            else
              phase  = 0._dp
              sine   = 0._dp
              cosine = 0._dp
            end if

            write(unit=iw,fmt="(9(es16.9,2x))") e , abs(rtorques(typetorque,sigma,i)) , real(rtorques(typetorque,sigma,i)) , aimag(rtorques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)
          end if
        end do
        ! Writing total torques
        iw = 9500+(typetorque-1)*3+sigma

        if(abs(total_torques(typetorque,sigma))>=1.e-15_dp) then
          phase  = atan2(aimag(total_torques(typetorque,sigma)),real(total_torques(typetorque,sigma)))
          sine   = real(total_torques(typetorque,sigma))/abs(total_torques(typetorque,sigma))
          cosine = aimag(total_torques(typetorque,sigma))/abs(total_torques(typetorque,sigma))
        else
          phase  = 0._dp
          sine   = 0._dp
          cosine = 0._dp
        end if

        write(unit=iw,fmt="(9(es16.9,2x))") e , abs(total_torques(typetorque,sigma)) , real(total_torques(typetorque,sigma)) , aimag(total_torques(typetorque,sigma)) , phase , sine , cosine , mtotal_spherical(2) , mtotal_spherical(3)

      end do
    end do

    call close_torque_files()

  end subroutine write_torques


  subroutine create_dc_torque_files()
  !> This subroutine creates all the files needed for the torques
    use mod_parameters, only: output, kount, renorm
    use mod_magnet, only: dcprefix, dcfield_dependence, dcfield, dc_header
    use mod_system, only: s => sys
    implicit none

    character(len=500)  :: varm
    integer :: i,sigma,typetorque,iw

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
          iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
          write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(dcprefix(kount)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
          close(unit=iw)
          if(renorm) then
            iw = iw+1000
            write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'r',a,a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(dcprefix(kount)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
            open (unit=iw, file=varm, status='replace', form='formatted')
            write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
            close(unit=iw)
          end if
        end do
        ! Total torque files
        iw = 95000+(typetorque-1)*3+sigma
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,a,'_',a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(dcprefix(kount)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,'  imaginary part of ',a,a,' ,  real part of ',a,a,'  , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma),trim(filename(typetorque)),direction(sigma)
        close(unit=iw)
      end do
    end do

  end subroutine create_dc_torque_files


  subroutine open_dc_torque_files()
  !> This subroutine opens all the files needed for the torques
    use mod_parameters, only: output, kount, missing_files, renorm
    use mod_magnet,     only: dcprefix, dcfield_dependence, dcfield
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none

    character(len=500)  :: varm
    integer :: i,sigma,typetorque,iw,err,errt=0

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
          iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
          write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(dcprefix(kount)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
          if(renorm) then
            iw = iw+1000
            write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'r',a,a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(dcprefix(kount)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
            open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
            errt = errt + err
            if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
          end if
        end do
        ! Total torque files
        iw = 95000+(typetorque-1)*3+sigma
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,a,'_',a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(output%hfr),trim(dcprefix(kount)),trim(filename(typetorque)),direction(sigma),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
      end do
    end do
    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openclose_dc_torque_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))

  end subroutine open_dc_torque_files


  subroutine close_dc_torque_files()
  !> This subroutine closes all the files needed for the torques
    use mod_parameters, only: renorm
    use mod_system, only: s => sys
    implicit none

    integer :: i,sigma,typetorque,iw

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
          iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i
          close(unit=iw)

          if(renorm) then
            iw = iw+1000
            close(unit=iw)
          end if
        end do
        ! Total torque files
        iw = 95000+(typetorque-1)*3+sigma
        close(unit=iw)
      end do
    end do

  end subroutine close_dc_torque_files


  subroutine write_dc_torques()
  !> This subroutine write all the torques into files
  !> (already opened with openclose_torque_files(1))
  !> Some information may also be written on the screen
    use mod_kind, only: dp
    use mod_parameters, only: renorm
    use mod_magnet,     only: mvec_spherical,mtotal_spherical,dc_fields,hw_count
    use mod_System,     only: s => sys
    implicit none
    integer      :: i,iw,sigma,typetorque
    real(dp) :: phase,sine,cosine

    call open_dc_torque_files()

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
          iw = 90000+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i

          if(abs(torques(typetorque,sigma,i))>=1.e-15_dp) then
            phase  = atan2(aimag(torques(typetorque,sigma,i)),real(torques(typetorque,sigma,i)))
            sine   = real(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
            cosine = aimag(torques(typetorque,sigma,i))/abs(torques(typetorque,sigma,i))
          else
            phase  = 0._dp
            sine   = 0._dp
            cosine = 0._dp
          end if
          write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(torques(typetorque,sigma,i)) , real(torques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)

          ! Writing renormalized torques
          if(renorm) then
            iw = iw+1000

            if(abs(rtorques(typetorque,sigma,i))>=1.e-15_dp) then
              phase  = atan2(aimag(rtorques(typetorque,sigma,i)),real(rtorques(typetorque,sigma,i)))
              sine   = real(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
              cosine = aimag(rtorques(typetorque,sigma,i))/abs(rtorques(typetorque,sigma,i))
            else
              phase  = 0._dp
              sine   = 0._dp
              cosine = 0._dp
            end if

            write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rtorques(typetorque,sigma,i)) , real(rtorques(typetorque,sigma,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)
          end if
        end do
        ! Writing total torques
      iw = 95000+(typetorque-1)*3+sigma

      if(abs(total_torques(typetorque,sigma))>=1.e-15_dp) then
        phase  = atan2(aimag(total_torques(typetorque,sigma)),real(total_torques(typetorque,sigma)))
        sine   = real(total_torques(typetorque,sigma))/abs(total_torques(typetorque,sigma))
        cosine = aimag(total_torques(typetorque,sigma))/abs(total_torques(typetorque,sigma))
      else
        phase  = 0._dp
        sine   = 0._dp
        cosine = 0._dp
      end if
      write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(total_torques(typetorque,sigma)) , real(total_torques(typetorque,sigma)) , phase , sine , cosine , mtotal_spherical(2) , mtotal_spherical(3)

      end do
    end do

    call close_dc_torque_files()

  end subroutine write_dc_torques


  subroutine sort_torques()
  !> This subroutine sorts torque files
    use mod_parameters, only: renorm,itype
    use mod_tools,      only: sort_file
    use mod_System,     only: s => sys
    implicit none
    integer  :: i,iw,sigma,typetorque,idc=1

    ! Opening torque files
    if(itype==9) then
      idc=10
      call open_dc_torque_files()
    else
      call open_torque_files()
    end if

    do typetorque=1,ntypetorque
      do sigma=1,3
        do i=1,s%nAtoms
            iw = 9000*idc+(typetorque-1)*s%nAtoms*3+(sigma-1)*s%nAtoms+i

            ! Sorting torque files
            call sort_file(iw)

            ! Sorting renormalized torque files
            if(renorm) then
              iw = iw+1000
              call sort_file(iw)
            end if
          end do
          iw = 9500*idc+(typetorque-1)*3+sigma
        ! Sorting total torque files
        call sort_file(iw)
      end do
    end do

    ! Closing torque files
    if(itype==9) then
      call close_dc_torque_files()
    else
      call close_torque_files()
    end if

  end subroutine sort_torques

  subroutine xc_torque_operator(torque)
  !> This subroutine defines the exchange-correlation torque operator/matrix at site i as
  !> $\vec{T}^xc_m = U \vec{M}\times\vec{S} $ , 
  !> with M being the magnetic moment and S = \sigma/2 being the spin operator.
    use mod_kind,       only: dp
    use mod_constants,  only: cZero, cOne, levi_civita, sigma => pauli_mat
    use mod_System,     only: s => sys
    use mod_magnet,     only: mdvec_cartesian
    implicit none
    complex(dp), dimension(s%nOrb2,s%nOrb2,3,s%nAtoms), intent(out) :: torque
    complex(dp), dimension(s%nOrb,s%nOrb) :: ident
    integer :: i,m,n,k,nOrb,nOrb2

    torque = cZero
    do i = 1, s%nAtoms
      nOrb =s%Types(s%Basis(i)%Material)%nOrb
      nOrb2=s%Types(s%Basis(i)%Material)%nOrb2

      ident = cZero
      do m=1,s%Types(s%Basis(i)%Material)%ndOrb
        ident(s%Types(s%Basis(i)%Material)%dOrbs(m),s%Types(s%Basis(i)%Material)%dOrbs(m)) = cOne
      end do

      do m = 1, 3
        do n = 1, 3
          do k = 1, 3
            torque(     1:nOrb ,     1:nOrb ,m,i) = torque(     1:nOrb ,     1:nOrb ,m,i) + mdvec_cartesian(n,i) * ident(1:nOrb,1:nOrb) * sigma(1,1,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,     1:nOrb ,m,i) = torque(nOrb+1:nOrb2,     1:nOrb ,m,i) + mdvec_cartesian(n,i) * ident(1:nOrb,1:nOrb) * sigma(2,1,k) * levi_civita(m,n,k)
            torque(     1:nOrb ,nOrb+1:nOrb2,m,i) = torque(     1:nOrb ,nOrb+1:nOrb2,m,i) + mdvec_cartesian(n,i) * ident(1:nOrb,1:nOrb) * sigma(1,2,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) = torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) + mdvec_cartesian(n,i) * ident(1:nOrb,1:nOrb) * sigma(2,2,k) * levi_civita(m,n,k)
          end do
        end do
      end do
      torque(:,:,:,i) = - 0.5_dp * s%Basis(i)%Um * torque(:,:,:,i)
    end do
  end subroutine xc_torque_operator


  subroutine SO_torque_operator(torque)
  !> This subroutine defines the spin-orbit torque operator/matrix at site i as
  !> $\vec{T}^so_i = lambda_i \vec{L}_i\times\vec{S} $ , 
  !> with L being the orbital operator and S = \sigma/2 being the spin operator.
    use mod_kind,       only: dp
    use mod_constants,  only: cZero, levi_civita, sigma => pauli_mat
    use mod_System,     only: s => sys
    implicit none
    integer :: i,m,n,k,mu,nu,mup,nup,mupp,nupp,nOrb,nOrb2
    complex(dp), dimension(s%nOrb2,s%nOrb2,3,s%nAtoms), intent(out) :: torque

    torque = cZero

    !! [sigma, H_SO]
    !! t_i^m = Lambda_i * sum_nk eps_mnk L_imunu^n * S_imunu^k
    !! t_m = lambda_i * sum_alpha_beta * sum_nk_gamma,zeta eps_mnk L^n_gamma,zeta * (sigma^k_alpha_beta/2)
    do i = 1, s%nAtoms
      nOrb =s%Types(s%Basis(i)%Material)%nOrb
      nOrb2=s%Types(s%Basis(i)%Material)%nOrb2
      do m = 1, 3
        do n = 1, 3
          do k = 1, 3
            torque(     1:nOrb ,     1:nOrb ,m,i) = torque(     1:nOrb ,     1:nOrb ,m,i) + s%Types(s%Basis(i)%Material)%lvec(:,:,n) * sigma(1,1,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,     1:nOrb ,m,i) = torque(nOrb+1:nOrb2,     1:nOrb ,m,i) + s%Types(s%Basis(i)%Material)%lvec(:,:,n) * sigma(2,1,k) * levi_civita(m,n,k)
            torque(     1:nOrb ,nOrb+1:nOrb2,m,i) = torque(     1:nOrb ,nOrb+1:nOrb2,m,i) + s%Types(s%Basis(i)%Material)%lvec(:,:,n) * sigma(1,2,k) * levi_civita(m,n,k)
            torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) = torque(nOrb+1:nOrb2,nOrb+1:nOrb2,m,i) + s%Types(s%Basis(i)%Material)%lvec(:,:,n) * sigma(2,2,k) * levi_civita(m,n,k)
          end do
        end do
      end do

      ! p block
      do nupp=1,s%Types(s%Basis(i)%Material)%npOrb
        nu = s%Types(s%Basis(i)%Material)%pOrbs(nupp)
        nup = nu + nOrb 
        do mupp=1,s%Types(s%Basis(i)%Material)%npOrb
          mu = s%Types(s%Basis(i)%Material)%pOrbs(mupp)
          mup = mu + nOrb 
          torque(mu ,nu ,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque(mu ,nu ,:,i)
          torque(mu ,nup,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque(mu ,nup,:,i)
          torque(mup,nu ,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque(mup,nu ,:,i)
          torque(mup,nup,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaP * torque(mup,nup,:,i)
        end do
      end do

      ! d block
      do nupp=1,s%Types(s%Basis(i)%Material)%ndOrb
        nu = s%Types(s%Basis(i)%Material)%dOrbs(nupp)
        nup = nu + nOrb 
        do mupp=1,s%Types(s%Basis(i)%Material)%ndOrb
          mu = s%Types(s%Basis(i)%Material)%dOrbs(mupp)
          mup = mu + nOrb 
          torque(mu ,nu ,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque(mu ,nu ,:,i)
          torque(mu ,nup,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque(mu ,nup,:,i)
          torque(mup,nu ,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque(mup,nu ,:,i)
          torque(mup,nup,:,i) = 0.5_dp * s%Types(s%Basis(i)%Material)%LambdaD * torque(mup,nup,:,i)
        end do
      end do
    end do
  end subroutine SO_torque_operator

end module mod_torques
