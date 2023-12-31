module mod_disturbances
  use mod_kind, only: dp
  implicit none
  ! Disturbances and renormalized disturbances: chd,sdx,sdy,sdz,ldx,ldy,ldz
  complex(dp),allocatable   :: disturbances(:,:)
  !! Disturbances
  complex(dp),allocatable   :: total_disturbances(:)
  !! Total Disturbances
  complex(dp),allocatable   :: rdisturbances(:,:)
  !! Renormalized disturbances
  complex(dp), dimension(:,:),   allocatable :: tchiorbiikl
  !! Full disturbance response function
  complex(dp), dimension(:,:,:), allocatable :: ldmat
  !! Angular Momentum Disturbance Matrix
  complex(dp),dimension(:),allocatable       :: sdmat
  !! Spin Disturbance Matrix

  character(len=5), dimension(7), parameter, private :: folder = ["CD", "SD", "SD", "SD", "LD", "LD", "LD"]
  character(len=2), dimension(7), parameter, private :: filename = ["Cd", "Sx", "Sy", "Sz", "Lx", "Ly", "Lz"]

contains

  subroutine allocate_disturbances()
  !! This subroutine allocates variables related to the disturbance calculation
    use mod_parameters, only: renorm,dimens,dimspinAtoms
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: rFreq, abortProgram
    implicit none
    integer           :: AllocateStatus

    if(rFreq(1) == 0) then
      allocate( disturbances(7,s%nAtoms),total_disturbances(7),sdmat(dimspinAtoms),ldmat(s%nAtoms,s%nOrb,s%nOrb), STAT = AllocateStatus )
      if (AllocateStatus/=0) call abortProgram("[allocate_disturbances] Not enough memory for: disturbances,total_disturbances,dmat,ldmat")

      if(renorm) then
        allocate( rdisturbances(7,s%nAtoms), STAT = AllocateStatus )
        if (AllocateStatus/=0) call abortProgram("[allocate_disturbances] Not enough memory for: rdisturbances")
      end if
    end if

    allocate( tchiorbiikl(dimens,4), STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[allocate_disturbances] Not enough memory for: tchiorbiikl")

  end subroutine allocate_disturbances

  subroutine deallocate_disturbances()
  !! This subroutine deallocates variables related to the disturbance calculation
    implicit none

    if(allocated(disturbances)) deallocate(disturbances)
    if(allocated(sdmat)) deallocate(sdmat)
    if(allocated(ldmat)) deallocate(ldmat)
    if(allocated(rdisturbances)) deallocate(rdisturbances)
    if(allocated(tchiorbiikl)) deallocate(tchiorbiikl)

  end subroutine deallocate_disturbances

  subroutine create_disturbance_files()
  !! This subroutine creates all the files needed for the disturbances
    use mod_parameters, only: output, renorm, renormnb
    use mod_system,     only: s => sys
    implicit none

    character(len=500)  :: varm
    integer :: i,j,iw

    do j=1,7
      do i=1,s%nAtoms
        iw = 3000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),filename(j),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)
        if(renorm) then
        iw = iw+s%nAtoms*7+1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/r',a,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),filename(j),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),renormnb,trim(output%suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)
        end if
      end do
      ! Total disturbances files
      iw = j+2*s%nAtoms*7+1
      write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),filename(j),trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
      open (unit=iw, file=varm, status='replace', form='formatted')
      write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
      close(unit=iw)
    end do

  end subroutine create_disturbance_files

  subroutine open_disturbance_files()
  !! This subroutine opens all the files needed for the disturbances
    use mod_parameters, only: output, renorm, renormnb, missing_files
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none

    character(len=500)  :: varm
    integer :: i,j,iw,err,errt=0

    do j=1,7
      do i=1,s%nAtoms
        iw = 3000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),filename(j),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        if(renorm) then
          iw = iw+s%nAtoms*7+1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/r',a,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),filename(j),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),renormnb,trim(output%suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end if
      end do
      ! Total disturbances files
      iw = j+2*s%nAtoms*7+1
      write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),filename(j),trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%EField),trim(output%suffix)
      open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
    end do
    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openclose_disturbance_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))

  end subroutine open_disturbance_files

  subroutine close_disturbance_files()
  !! This subroutine closes all the files needed for the disturbances
    use mod_parameters, only: renorm
    use mod_system, only: s => sys
    implicit none

    integer :: i,j,iw

    do j=1,7
      do i=1,s%nAtoms
        iw = 3000+(i-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+s%nAtoms*7+1
          close(unit=iw)
        end if
      end do
      ! Total disturbances files
      iw = j+2*s%nAtoms*7+1
      close(unit=iw)

    end do
  end subroutine close_disturbance_files

  ! This subroutine write all the disturbances into files
  ! (already opened with openclose_disturbance_files(1))
  ! Some information may also be written on the screen
  subroutine write_disturbances(e)
    use mod_kind,       only: dp
    use mod_parameters, only: output,lwriteonscreen
    use mod_magnet,     only: mvec_spherical,mtotal_spherical
    use mod_System,     only: s => sys
    implicit none
    integer  :: i,j,iw
    real(dp),intent(in) :: e
    real(dp) :: phase,sine,cosine

    call open_disturbance_files()

    if(lwriteonscreen) write(output%unit_loop,"(' ################# Disturbances: #################')")
    ! Writing Spin, Charge and Orbital disturbances
    do i=1,s%nAtoms
      if(lwriteonscreen) then
        write(output%unit_loop,"('|--------------- Energy = ',es11.4,' , Plane: ',i0,' ---------------|')") e,i

        write(output%unit_loop,"('     Cd  = (',es16.9,') + i(',es16.9,')')") real(disturbances(1,i)),aimag(disturbances(1,i))
        write(output%unit_loop,"(' abs(Cd) = ',es16.9)") abs(disturbances(1,i))
        write(output%unit_loop,"('atan(Cd) = ',es16.9)") atan2(aimag(disturbances(1,i)),real(disturbances(1,i)))

        write(output%unit_loop,"('     Sdx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(2,i)),aimag(disturbances(2,i))
        write(output%unit_loop,"(' abs(Sdx) = ',es16.9)") abs(disturbances(2,i))
        write(output%unit_loop,"('atan(Sdx) = ',es16.9)") atan2(aimag(disturbances(2,i)),real(disturbances(2,i)))

        write(output%unit_loop,"('     Sdy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(3,i)),aimag(disturbances(3,i))
        write(output%unit_loop,"(' abs(Sdy) = ',es16.9)") abs(disturbances(3,i))
        write(output%unit_loop,"('atan(Sdy) = ',es16.9)") atan2(aimag(disturbances(3,i)),real(disturbances(3,i)))

        write(output%unit_loop,"('     Sdz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(4,i)),aimag(disturbances(4,i))
        write(output%unit_loop,"(' abs(Sdz) = ',es16.9)") abs(disturbances(4,i))
        write(output%unit_loop,"('atan(Sdz) = ',es16.9)") atan2(aimag(disturbances(4,i)),real(disturbances(4,i)))

        write(output%unit_loop,"('     Ldx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(5,i)),aimag(disturbances(5,i))
        write(output%unit_loop,"(' abs(Ldx) = ',es16.9)") abs(disturbances(5,i))
        write(output%unit_loop,"('atan(Ldx) = ',es16.9)") atan2(aimag(disturbances(5,i)),real(disturbances(5,i)))

        write(output%unit_loop,"('     Ldy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(6,i)),aimag(disturbances(6,i))
        write(output%unit_loop,"(' abs(Ldy) = ',es16.9)") abs(disturbances(6,i))
        write(output%unit_loop,"('atan(Ldy) = ',es16.9)") atan2(aimag(disturbances(6,i)),real(disturbances(6,i)))

        write(output%unit_loop,"('     Ldz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(7,i)),aimag(disturbances(7,i))
        write(output%unit_loop,"(' abs(Ldz) = ',es16.9)") abs(disturbances(7,i))
        write(output%unit_loop,"('atan(Ldz) = ',es16.9)") atan2(aimag(disturbances(7,i)),real(disturbances(7,i)))
      end if


      ! Writing disturbances
      iw = 3000+(i-1)*size(filename)
      do j=1,size(filename)
        if(abs(disturbances(j,i))>=1.e-15_dp) then
          phase  = atan2(aimag(disturbances(j,i)),real(disturbances(j,i)))
          if (phase <=1.e-15_dp) phase = 0._dp
          sine   = real(disturbances(j,i))/abs(disturbances(j,i))
          if (sine <=1.e-15_dp) sine = 0._dp
          cosine = aimag(disturbances(j,i))/abs(disturbances(j,i))
          if (cosine <=1.e-15_dp) cosine = 0._dp
        else
          phase  = 0._dp
          sine   = 0._dp
          cosine = 0._dp
        end if

        write(unit=iw+j,fmt="(9(es16.9,2x))") e, abs(disturbances(j,i)) , real(disturbances(j,i)) , aimag(disturbances(j,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing renormalized disturbances
        ! if(renorm) then
        !   if(abs(rdisturbances(j,i))>=1.e-15_dp) then
        !     phase  = atan2(aimag(rdisturbances(j,i)),real(rdisturbances(j,i)))
        !     if (phase <=1.e-15_dp) phase = 0._dp
        !     sine   = real(rdisturbances(j,i))/abs(rdisturbances(j,i))
        !     if (sine <=1.e-15_dp) sine = 0._dp
        !     cosine = aimag(rdisturbances(j,i))/abs(rdisturbances(j,i))
        !     if (cosine <=1.e-15_dp) cosine = 0._dp
        !   else
        !     phase  = 0._dp
        !     sine   = 0._dp
        !     cosine = 0._dp
        !   end if

        !   write(unit=iw+1000+j,fmt="(9(es16.9,2x))") e, abs(rdisturbances(j,i)) , real(rdisturbances(j,i)) , aimag(rdisturbances(j,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! end if

      end do
    end do

    ! Writing total charge disturbance
    iw = 2*s%nAtoms*7+1
    do j=1,size(filename)
      if(abs(total_disturbances(j))>=1.e-15_dp) then
        phase  = atan2(aimag(total_disturbances(j)),real(total_disturbances(j)))
        if (phase <=1.e-15_dp) phase = 0._dp
        sine   = real(total_disturbances(j))/abs(total_disturbances(j))
        if (sine <=1.e-15_dp) sine = 0._dp
        cosine = aimag(total_disturbances(j))/abs(total_disturbances(j))
        if (cosine <=1.e-15_dp) cosine = 0._dp
      else
        phase  = 0._dp
        sine   = 0._dp
        cosine = 0._dp
      end if

      write(unit=iw+j,fmt="(9(es16.9,2x))") e, abs(total_disturbances(j)) , real(total_disturbances(j)) , aimag(total_disturbances(j)) , phase , sine , cosine , mtotal_spherical(2) , mtotal_spherical(3)
    end do

    call close_disturbance_files()
  end subroutine write_disturbances

  subroutine create_dc_disturbance_files
  !! This subroutine creates all the files needed for the dc-limit disturbances
    use mod_parameters, only: kount, output, renorm, renormnb
    use mod_magnet,     only: dcprefix, dcfield_dependence, dcfield, dc_header
    use mod_system,     only: s => sys
    implicit none
    character(len=500)  :: varm
    integer :: i,j,iw

    do i=1,s%nAtoms
      do j=1,7
        iw = 30000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,', amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)
        if(renorm) then
          iw = iw+s%nAtoms*7+1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'r',a,'_',a,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),renormnb,trim(output%suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#',a,', amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do
      ! Total disturbances files
      iw = j+10*2*s%nAtoms*7+1
      write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_',a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
      open (unit=iw, file=varm, status='replace', form='formatted')
      write(unit=iw, fmt="('#',a,', amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
      close(unit=iw)
    end do
  end subroutine create_dc_disturbance_files

  subroutine open_dc_disturbance_files
  ! This subroutine opens all the files needed for the dc-limit disturbances
    use mod_parameters, only: output, kount, renorm, renormnb, missing_files
    use mod_magnet,     only: dcprefix, dcfield_dependence, dcfield
    use mod_mpi_pars,   only: abortProgram
    use mod_system,     only: s => sys
    implicit none

    character(len=500)  :: varm
    integer :: i,j,iw,err,errt=0

    do j=1,7
      do i=1,s%nAtoms
        iw = 30000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        if(renorm) then
          iw = iw+s%nAtoms*7+1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,'r',a,'_',a,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),renormnb,trim(output%suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end if
      end do
      ! Total disturbances files
      iw = j+10*2*s%nAtoms*7+1
      write(varm,"('./results/',a1,'SOC/',a,'/',a,a,'/',a,a,'_',a,'_total',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder(j)),trim(output%hfr),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%EField),trim(output%suffix)
      open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
    end do

    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openclose_dc_disturbance_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
  end subroutine open_dc_disturbance_files

  subroutine close_dc_disturbance_files()
  !! This subroutine closes all the files needed for the dc-limit disturbances
    use mod_parameters, only: renorm
    use mod_system,     only: s => sys
    implicit none
    integer :: i,j,iw

    do j=1,7
      do i=1,s%nAtoms
        iw = 30000+(i-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+s%nAtoms*7+1
          close(unit=iw)
        end if
      end do
      ! Total disturbances files
      iw = j+10*2*s%nAtoms*7+1
      close(unit=iw)
    end do
  end subroutine close_dc_disturbance_files

  subroutine write_dc_disturbances()
    !! This subroutine write all the dc-limit disturbances into files
    !! (already opened with openclose_dc_disturbance_files(1))
    !! Some information may also be written on the screen
    use mod_parameters, only: output,lwriteonscreen
    use mod_magnet,     only: mvec_spherical, mtotal_spherical, dcfield, dcfield_dependence, dc_fields, hw_count
    use mod_System,     only: s => sys
    implicit none
    integer  :: i,j,iw
    real(dp) :: phase,sine,cosine

    call open_dc_disturbance_files()

    if(lwriteonscreen) write(output%unit_loop,"(' ################# Disturbances: #################')")
    ! Writing Spin, Charge and Orbital disturbances
    do i=1, s%nAtoms
      if(lwriteonscreen) then
        write(output%unit_loop,"('|--------------- ',a,' = ',a,' , Plane: ',i0,' ---------------|')") trim(dcfield(dcfield_dependence)),trim(dc_fields(hw_count)),i

        write(output%unit_loop,"('     Cd  = (',es16.9,') + i(',es16.9,')')") real(disturbances(1,i)),aimag(disturbances(1,i))
        write(output%unit_loop,"(' abs(Cd) = ',es16.9)") abs(disturbances(1,i))
        write(output%unit_loop,"('atan(Cd) = ',es16.9)") atan2(aimag(disturbances(1,i)),real(disturbances(1,i)))

        write(output%unit_loop,"('     Sdx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(2,i)),aimag(disturbances(2,i))
        write(output%unit_loop,"(' abs(Sdx) = ',es16.9)") abs(disturbances(2,i))
        write(output%unit_loop,"('atan(Sdx) = ',es16.9)") atan2(aimag(disturbances(2,i)),real(disturbances(2,i)))

        write(output%unit_loop,"('     Sdy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(3,i)),aimag(disturbances(3,i))
        write(output%unit_loop,"(' abs(Sdy) = ',es16.9)") abs(disturbances(3,i))
        write(output%unit_loop,"('atan(Sdy) = ',es16.9)") atan2(aimag(disturbances(3,i)),real(disturbances(3,i)))

        write(output%unit_loop,"('     Sdz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(4,i)),aimag(disturbances(4,i))
        write(output%unit_loop,"(' abs(Sdz) = ',es16.9)") abs(disturbances(4,i))
        write(output%unit_loop,"('atan(Sdz) = ',es16.9)") atan2(aimag(disturbances(4,i)),real(disturbances(4,i)))

        write(output%unit_loop,"('     Ldx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(5,i)),aimag(disturbances(5,i))
        write(output%unit_loop,"(' abs(Ldx) = ',es16.9)") abs(disturbances(5,i))
        write(output%unit_loop,"('atan(Ldx) = ',es16.9)") atan2(aimag(disturbances(5,i)),real(disturbances(5,i)))

        write(output%unit_loop,"('     Ldy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(6,i)),aimag(disturbances(6,i))
        write(output%unit_loop,"(' abs(Ldy) = ',es16.9)") abs(disturbances(6,i))
        write(output%unit_loop,"('atan(Ldy) = ',es16.9)") atan2(aimag(disturbances(6,i)),real(disturbances(6,i)))

        write(output%unit_loop,"('     Ldz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(7,i)),aimag(disturbances(7,i))
        write(output%unit_loop,"(' abs(Ldz) = ',es16.9)") abs(disturbances(7,i))
        write(output%unit_loop,"('atan(Ldz) = ',es16.9)") atan2(aimag(disturbances(7,i)),real(disturbances(7,i)))
      end if

      ! Writing disturbances
      iw = 30000+(i-1)*size(filename)
      do j=1,size(filename)
        if(abs(disturbances(j,i))>=1.e-15_dp) then
          phase  = atan2(aimag(disturbances(j,i)),real(disturbances(j,i)))
          if (phase <=1.e-15_dp) phase = 0._dp
          sine   = real(disturbances(j,i))/abs(disturbances(j,i))
          if (sine <=1.e-15_dp) sine = 0._dp
          cosine = aimag(disturbances(j,i))/abs(disturbances(j,i))
          if (cosine <=1.e-15_dp) cosine = 0._dp
        else
          phase  = 0._dp
          sine   = 0._dp
          cosine = 0._dp
        end if

        write(unit=iw+j,fmt="(a,2x,8(es16.9,2x))") trim(dc_fields(hw_count)) , abs(disturbances(j,i)) , real(disturbances(j,i)) , aimag(disturbances(j,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing renormalized disturbances
        ! if(renorm) then
        !   if(abs(rdisturbances(j,i))>=1.e-15_dp) then
        !     phase  = atan2(aimag(rdisturbances(j,i)),real(rdisturbances(j,i)))
        !     if (phase <=1.e-15_dp) phase = 0._dp
        !     sine   = real(rdisturbances(j,i))/abs(rdisturbances(j,i))
        !     if (sine <=1.e-15_dp) sine = 0._dp
        !     cosine = aimag(rdisturbances(j,i))/abs(rdisturbances(j,i))
        !     if (cosine <=1.e-15_dp) cosine = 0._dp
        !   else
        !     phase  = 0._dp
        !     sine   = 0._dp
        !     cosine = 0._dp
        !   end if

        !   write(unit=iw+1000+j,fmt="(a,2x,8(es16.9,2x))") trim(dc_fields(hw_count)) , abs(rdisturbances(j,i)) , real(rdisturbances(j,i)) , aimag(rdisturbances(j,i)) , phase , sine , cosine , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! end if
      end do
    end do

    ! Writing total charge disturbance
    iw = 10*2*s%nAtoms*7+1
    do j=1,size(filename)
      if(abs(total_disturbances(j))>=1.e-15_dp) then
        phase  = atan2(aimag(total_disturbances(j)),real(total_disturbances(j)))
        if (phase <=1.e-15_dp) phase = 0._dp
        sine   = real(total_disturbances(j))/abs(total_disturbances(j))
        if (sine <=1.e-15_dp) sine = 0._dp
        cosine = aimag(total_disturbances(j))/abs(total_disturbances(j))
        if (cosine <=1.e-15_dp) cosine = 0._dp
      else
        phase  = 0._dp
        sine   = 0._dp
        cosine = 0._dp
      end if

      write(unit=iw+j,fmt="(a,2x,8(es16.9,2x))") trim(dc_fields(hw_count)) , abs(total_disturbances(j)) , real(total_disturbances(j)) , aimag(total_disturbances(j)) , phase , sine , cosine , mtotal_spherical(2) , mtotal_spherical(3)
    end do

    call close_dc_disturbance_files()

  end subroutine write_dc_disturbances

  ! This subroutine sorts disturbance files
  subroutine sort_disturbances()
    use mod_parameters, only: renorm,itype
    use mod_tools, only: sort_file
    use mod_System, only: s => sys
    implicit none
    integer :: i,j,iw,idc=1

    ! Opening disturbance files
    if(itype==9) then
      idc=10
      call open_dc_disturbance_files()
    else
      call open_disturbance_files()
    end if

    ! Sorting total disturbance
    iw = idc*2*s%nAtoms*7+1
    do j=1,7
      call sort_file(iw+j)
    end do


    do i=1, s%nAtoms
      ! Sorting disturbance files
      iw = 3000*idc+(i-1)*7
      do j=1,7
        call sort_file(iw+j)
      end do

      ! Sorting renormalized disturbances
      if(renorm) then
        do j=1+2*s%nAtoms*7+1,7+2*s%nAtoms*7+1
          call sort_file(iw+j)
        end do
      end if
    end do

    ! Closing disturbance files
    if(itype==9) then
      call close_dc_disturbance_files()
    else
      call close_disturbance_files()
    end if

  end subroutine sort_disturbances

end module mod_disturbances
