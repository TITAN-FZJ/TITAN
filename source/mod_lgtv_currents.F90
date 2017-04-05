module mod_lgtv_currents
  use mod_f90_kind
  implicit none
  ! Longitudinal and transverse currents
  complex(double),allocatable   :: long_currents(:,:),total_long_currents(:),transv_currents(:,:),total_transv_currents(:)

contains

  ! This subroutine allocates variables related to the longitudinal and transverse currents calculation
  subroutine allocate_lgtv_currents()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters, only: Npl,outputunit
    implicit none
    integer           :: AllocateStatus

    if(myrank_row==0) then
      allocate( long_currents(7,Npl),total_long_currents(7),transv_currents(7,Npl),total_transv_currents(7), STAT = AllocateStatus )
      if (AllocateStatus/=0) then
         write(outputunit,"('[allocate_lgtv_currents] Not enough memory for: long_currents,total_long_currents,transv_currents,total_transv_currents')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if

    return
  end subroutine allocate_lgtv_currents

  ! This subroutine allocates variables related to the longitudinal and transverse currents calculation
  subroutine deallocate_lgtv_currents()
    use mod_f90_kind
    use mod_mpi_pars
    implicit none

    if(myrank_row==0) deallocate(long_currents,total_long_currents,transv_currents,total_transv_currents)

    return
  end subroutine deallocate_lgtv_currents

  ! This subroutine opens and closes all the files needed for the longitudinal and transverse currents
  subroutine openclose_lgtv_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=1)    :: SOCc
    character(len=6)    :: folder(7),typec(2)
    character(len=3)    :: filename(7)
    integer :: i,j,k,iw,iflag,err,errt=0

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

    folder(1) = "CC"
    folder(2) = "SC"
    folder(3) = "SC"
    folder(4) = "SC"
    folder(5) = "LC"
    folder(6) = "LC"
    folder(7) = "LC"
    if(lhfresponses) then
      do i=1,7
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

    typec(1) = "longit"
    typec(2) = "transv"

    if(iflag==0) then
      ! Header for longitudinal and transverse currents per plane
      do i=1,Npl ; do j=1,7 ; do k=1,2
        iw = 8300+(i-1)*7*2+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),typec(k),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#     energy     , real part of ',a,', imag part of ',a,',   phase of ',a,'  ')") filename(j),filename(j),filename(j)
        close(unit=iw)
      end do ; end do ; end do
      ! Header for total longitudinal and transverse currents
      do j=1,7 ; do k=1,2
        iw = 8500+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,a,'_total_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),typec(k),parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#     energy     , real part of ',a,', imag part of ',a,',   phase of ',a,'  ')") filename(j),filename(j),filename(j)
        close(unit=iw)
      end do ; end do
    else if(iflag==1) then
      ! Longitudinal and transverse currents per plane
      do i=1,Npl ; do j=1,7 ; do k=1,2
        iw = 8300+(i-1)*7*2+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),typec(k),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do ; end do
      ! Total longitudinal and transverse currents
      do j=1,7 ; do k=1,2
        iw = 8500+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,a,'_total_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),typec(k),parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do
      ! Stop if some file does not exist
      if(errt/=0) then
        write(outputunit,"(a,i0,a)") "[openclose_lgtv_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      ! Closing longitudinal and transverse currents
      do i=1,Npl ; do j=1,7 ; do k=1,2
        iw = 8300+(i-1)*7*2+(j-1)*2+k
        close(unit=iw)
      end do ; end do ; end do
      do j=1,7 ; do k=1,2
        iw = 8500+(j-1)*2+k
        close(unit=iw)
      end do ; end do
    end if

    return
  end subroutine openclose_lgtv_files

  ! This subroutine opens and closes all the files needed for the field dependent longitudinal and transverse currents
  subroutine openclose_dc_lgtv_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=1)    :: SOCc
    character(len=6)    :: folder(7),typec(2)
    character(len=3)    :: filename(7)
    integer :: i,j,k,iw,iflag,err,errt=0

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

    folder(1) = "CC"
    folder(2) = "SC"
    folder(3) = "SC"
    folder(4) = "SC"
    folder(5) = "LC"
    folder(6) = "LC"
    folder(7) = "LC"
    if(lhfresponses) then
      do i=1,7
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

    typec(1) = "longit"
    typec(2) = "transv"

    if(iflag==0) then
      ! Header for longitudinal and transverse currents per plane
      do i=1,Npl ; do j=1,7 ; do k=1,2
        iw = 83000+(i-1)*7*2+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),typec(k),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#',a,' real part of ',a,', imag part of ',a,',   phase of ',a,'  , mag angle theta ,  mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j)
        close(unit=iw)
      end do ; end do ; end do
      ! Header for total longitudinal and transverse currents
      do j=1,7 ; do k=1,2
        iw = 85000+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,a,'_',a,'_total_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),typec(k),trim(dcfield(dcfield_dependence)),parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#',a,' real part of ',a,', imag part of ',a,',   phase of ',a,'  , mag angle theta ,  mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j)
        close(unit=iw)
      end do ; end do
    else if(iflag==1) then
      ! Header for longitudinal and transverse currents per plane
      do i=1,Npl ; do j=1,7 ; do k=1,2
        iw = 83000+(i-1)*7*2+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),typec(k),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do ; end do
      ! Header for total longitudinal and transverse currents
      do j=1,7 ; do k=1,2
        iw = 85000+(j-1)*2+k
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,a,'_',a,'_total_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),typec(k),trim(dcfield(dcfield_dependence)),parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do

      ! Stop if some file does not exist
      if(errt/=0) then
        write(outputunit,"(a,i0,a)") "[openclose_dc_sha_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      ! SHA
      do i=1,Npl ; do j=1,7 ; do k=1,2
        iw = 83000+(i-1)*7*2+(j-1)*2+k
        close(unit=iw)
      end do ; end do ; end do
      do j=1,7 ; do k=1,2
        iw = 85000+(j-1)*2+k
        close(unit=iw)
      end do ; end do

    end if

    return
  end subroutine openclose_dc_lgtv_files


  ! This subroutine write all the longitudinal and transverse currents into files
  ! (already opened with openclose_lgtv_files(1))
  subroutine write_lgtv_currents(e)
    use mod_f90_kind
    use mod_parameters, only: Npl, mmlayermag
    use mod_magnet, only: mvec_spherical
    implicit none
    integer  :: i,j,iw
    real(double),intent(in) :: e

    do i=1,Npl ; do j=1,7
      iw = 8300+(i-1)*7*2+(j-1)*2+1
      write(unit=iw,fmt="(6(es16.9,2x))") e , real(long_currents(j,i)) , aimag(long_currents(j,i)) , atan2(aimag(long_currents(j,i)),real(long_currents(j,i))) , mvec_spherical(i,2) , mvec_spherical(i,3)
      iw = iw+1
      write(unit=iw,fmt="(6(es16.9,2x))") e , real(transv_currents(j,i)) , aimag(transv_currents(j,i)) , atan2(aimag(transv_currents(j,i)),real(transv_currents(j,i))) , mvec_spherical(i,2) , mvec_spherical(i,3)
    end do ; end do
    do j=1,7
      iw = 8500+(j-1)*2+1
      write(unit=iw,fmt="(6(es16.9,2x))") e , real(total_long_currents(j)) , aimag(total_long_currents(j)) , atan2(aimag(total_long_currents(j)),real(total_long_currents(j))) , mvec_spherical(mmlayermag(1)-1,2) , mvec_spherical(mmlayermag(1)-1,3)
      iw = iw+1
      write(unit=iw,fmt="(6(es16.9,2x))") e , real(total_transv_currents(j)) , aimag(total_transv_currents(j)) , atan2(aimag(total_transv_currents(j)),real(total_transv_currents(j))) , mvec_spherical(mmlayermag(1)-1,2) , mvec_spherical(mmlayermag(1)-1,3)
    end do

    return
  end subroutine write_lgtv_currents


  ! This subroutine write all the longitudinal and transverse currents for fixed frequency into files
  ! (already opened with openclose_dc_lgtv_files(1))
  subroutine write_dc_lgtv_currents()
    use mod_f90_kind
    use mod_parameters
    use mod_magnet, only: mvec_spherical
    implicit none
    integer  :: i,j,iw

    do i=1,Npl ; do j=1,7
      iw = 83000+(i-1)*7*2+(j-1)*2+1
      write(unit=iw,fmt="(a,2x,5(es16.9,2x))") trim(dc_fields(hw_count)) , real(long_currents(j,i)) , aimag(long_currents(j,i)) , atan2(aimag(long_currents(j,i)),real(long_currents(j,i))) , mvec_spherical(i,2) , mvec_spherical(i,3)
      iw = iw+1
      write(unit=iw,fmt="(a,2x,5(es16.9,2x))") trim(dc_fields(hw_count)) , real(transv_currents(j,i)) , aimag(transv_currents(j,i)) , atan2(aimag(transv_currents(j,i)),real(transv_currents(j,i))) , mvec_spherical(i,2) , mvec_spherical(i,3)
    end do ; end do
    do j=1,7
      iw = 85000+(j-1)*2+1
      write(unit=iw,fmt="(a,2x,5(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_long_currents(j)) , aimag(total_long_currents(j)) , atan2(aimag(total_long_currents(j)),real(total_long_currents(j))) , mvec_spherical(mmlayermag(1)-1,2) , mvec_spherical(mmlayermag(1)-1,3)
      iw = iw+1
      write(unit=iw,fmt="(a,2x,5(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_transv_currents(j)) , aimag(total_transv_currents(j)) , atan2(aimag(total_transv_currents(j)),real(total_transv_currents(j))) , mvec_spherical(mmlayermag(1)-1,2) , mvec_spherical(mmlayermag(1)-1,3)
    end do

    return
  end subroutine write_dc_lgtv_currents

  ! This subroutine sorts longitudinal and transverse currents files
  subroutine sort_lgtv_currents()
    use mod_f90_kind
    use mod_parameters, only: Npl,itype
    use mod_tools, only: sort_file
    implicit none
    integer  :: i,j,k,iw,idc=1

    ! Opening longitudinal and transverse currents files
    if(itype==9) then
      idc=10
      call openclose_dc_lgtv_files(1)
    else
      call openclose_lgtv_files(1)
    end if

    ! Sorting longitudinal and transverse currents files
    do i=1,Npl ; do j=1,7 ; do k=1,2
      iw = 8300*idc+(i-1)*7*2+(j-1)*2+k
      call sort_file(iw,.true.)
    end do ; end do ; end do
    do j=1,7 ; do k=1,2
      iw = 8500*idc+(j-1)*2+k
      call sort_file(iw,.true.)
    end do ; end do

    ! Closing longitudinal and transverse currents files
    if(itype==9) then
      call openclose_dc_lgtv_files(2)
    else
      call openclose_lgtv_files(2)
    end if

    return
  end subroutine sort_lgtv_currents

  ! This subroutine reads currents from existing files
  ! and calculates longitudinal and transverse currents
  subroutine read_calculate_lgtv_currents()
    use mod_parameters
    use mod_constants
    use mod_magnet
    use mod_progress
    use mod_currents
    use mod_mpi_pars
    use mod_lattice
    use mod_tools
    implicit none
    character(len=50) :: formatvar
    integer           :: i,j,k,iw,neighbor,rows,cols,idc,ie,iemin
    real(double)      :: idia
    real(double)   , allocatable :: data(:,:),x(:,:),e(:),mangles(:,:,:)
    complex(double), allocatable :: currents_from_file(:,:,:,:),total_currents_from_file(:,:,:)

    write(outputunit,"('[read_calculate_lgtv_currents] Reading current files and calculating longitudinal and transverse currents...')")
    ! Opening disturbance files
    select case(itype)
    case(8)
      ! Allocating current and longitudinal and transverse currents variables
      call allocate_currents()
      call allocate_lgtv_currents()

      ! Creating longitudinal and transverse currents files with headers
      call openclose_lgtv_files(0)

      ! Opening current and longitudinal and transverse currents files
      call openclose_lgtv_files(1)
      call openclose_currents_files(1)

      ! Reading currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j

        if(.not.allocated(data)) then ! Getting information from first file
          ! Obtaining number of rows and cols in the file
          call number_of_rows_cols(iw,rows,cols)
          ! Allocating variables to read and store data
          allocate(data(rows,cols),e(rows),currents_from_file(7,n0sc1:n0sc2,Npl,rows),total_currents_from_file(7,n0sc1:n0sc2,rows))

          ! Reading data and storing to variable 'data'
          call read_data(iw,rows,cols,data)
          ! Getting energy list
          e(:) = data(:,1)

          ! Minimum value of energy of the file
          iemin = minloc(abs(e),1)
        else
          ! Reading data and storing to variable 'data'
          call read_data(iw,rows,cols,data)

          ! Checking if energy list is the same
          if(sum(abs(e(:)-data(:,1)))>1.d-8) then
            write(outputunit,"('[read_calculate_lgtv_currents] Different energies on current files!')")
            call MPI_Finalize(ierr)
            stop
          end if
        end if
        ! Obtaining diamagnetic current (if the minimum energy is small enough)
        if(e(iemin)<=1.d-5) then
          idia = data(iemin,4)
          data(:,4) = data(:,4) - idia
        else
          write(outputunit,"('[read_calculate_lgtv_currents] WARNING: Diamagnetic current was not removed!')")
        end if

        currents_from_file(j,neighbor,i,:) = cmplx(data(:,3),data(:,4))
      end do ; end do ; end do

      ! Reading total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 7000+(neighbor-1)*7+j

        ! Reading data and storing to variable 'data'
        call read_data(iw,rows,cols,data)

        ! Checking if energy list is the same
        if(sum(abs(e(:)-data(:,1)))>1.d-8) then
          write(outputunit,"('[read_calculate_lgtv_currents] Different energies on current files!')")
          call MPI_Finalize(ierr)
          stop
        end if
        ! Obtaining diamagnetic current (if the minimum energy is small enough)
        if(e(iemin)<=1.d-5) then
          idia = data(iemin,4)
          data(:,4) = data(:,4) - idia
        else
          write(outputunit,"('[read_calculate_lgtv_currents] WARNING: Diamagnetic current was not removed!')")
        end if

        total_currents_from_file(j,neighbor,:) = cmplx(data(:,3),data(:,4))
      end do ; end do

      ! Calculating longitudinal and transverse currents for all energies
      do ie=1,rows
        currents(:,:,:) = currents_from_file(:,:,:,ie)
        total_currents(:,:) = total_currents_from_file(:,:,ie)
! write(*,*) real(currents(1,1,1)),real(currents(1,5,1)),real(currents(1,6,1))

        call calculate_lgtv_currents()
        call write_lgtv_currents(e(ie))
! write(*,*) real(currents(1,1,1)),real(currents(1,4,1))
! write(*,*) e(ie),real(long_currents(1,1))
      end do

      ! Closing files
      call openclose_currents_files(2)
      call openclose_lgtv_files(2)

      ! Sorting files
      ! call sort_lgtv_currents()

      ! Deallocating variables
      deallocate(data,e,currents_from_file,total_currents_from_file)
      call deallocate_currents()
      call deallocate_lgtv_currents()
    case(9)
      if(hw_count==1) then

        ! Allocating current and longitudinal and transverse currents variables
        call allocate_currents()
        call allocate_lgtv_currents()

        energy_loop: do count=1,npt1

          ! Creating longitudinal and transverse currents files with headers
          call openclose_dc_lgtv_files(0)

          ! Opening current and longitudinal and transverse currents files
          call openclose_dc_lgtv_files(1)
          call openclose_dc_currents_files(1)

          idc = ceiling((dble(dcfield_dependence)-0.01d0)/3.d0)
          ! Currents per plane per neighbor
          do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
            iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7+j

            if(.not.allocated(data)) then ! Getting information from first file
              ! Obtaining number of rows and cols in the file
              call number_of_rows_cols(iw,rows,cols)

              ! If number of rows is different than total_hw_npt1, correct dimensions
              if(rows/=total_hw_npt1) then
                write(outputunit,"('[read_calculate_lgtv_currents] WARNING: Number of rows different than field values!')")
                deallocate(dc_fields)
                allocate(dc_fields(rows))
                total_hw_npt1 = rows
              end if

              ! Allocating variables to read and sort data
              allocate(data(rows,cols),x(rows,idc),currents_from_file(7,n0sc1:n0sc2,Npl,rows),total_currents_from_file(7,n0sc1:n0sc2,rows),mangles(rows,Npl,2))
              mangles = 999.d0
              ! Reading data and storing to variable 'data'
              call read_data(iw,rows,cols,data)
              ! Getting fields list
              do k=1,idc
                x(:,k) = data(:,k)
              end do
              ! Getting magnetization angles
              mangles(:,i,:) = data(:,cols-1:cols)
            else
              ! Reading data and storing to variable 'data'
              call read_data(iw,rows,cols,data)

              ! Checking if fields list is the same
              if(sum(abs(x(:,1:idc)-data(:,1:idc)))>1.d-8) then
                write(outputunit,"('[read_calculate_lgtv_currents] Different fields on current files!')")
                call MPI_Finalize(ierr)
                stop
              end if

              ! Checking if angles are the same or getting next plane values
              if(sum(abs(mangles(:,i,:)))<10.d0) then
                if(sum(abs(mangles(:,i,:)-data(:,cols-1:cols)))>1.d-8) then
                  write(outputunit,"('[read_calculate_lgtv_currents] Different angles on current files!')")
                  call MPI_Finalize(ierr)
                  stop
                end if
              else
                ! Getting magnetization angles
                mangles(:,i,:) = data(:,cols-1:cols)
              end if
            end if

            currents_from_file(j,neighbor,i,:) = cmplx(data(:,1+idc),data(:,2+idc))
          end do ; end do ; end do

          ! Total currents for each neighbor direction
          do neighbor=n0sc1,n0sc2 ; do j=1,7
            iw = 70000+(neighbor-1)*7+j

            ! Reading data and storing to variable 'data'
            call read_data(iw,rows,cols,data)

            ! Checking if fields list is the same
            if(sum(abs(x(:,1:idc)-data(:,1:idc)))>1.d-8) then
              write(outputunit,"('[read_calculate_lgtv_currents] Different fields on current files!')")
              call MPI_Finalize(ierr)
              stop
            end if

            ! Checking if angles are the same
            if(sum(abs(mangles(:,mmlayermag(1)-1,:)-data(:,cols-1:cols)))>1.d-8) then
              write(outputunit,"('[read_calculate_lgtv_currents] Different angles on current files!')")
              call MPI_Finalize(ierr)
              stop
            end if

            total_currents_from_file(j,neighbor,:) = cmplx(data(:,1+idc),data(:,2+idc))
          end do ; end do

          ! Calculating longitudinal and transverse currents for all fields
          write(formatvar,fmt="(a,i0,a)") '(',idc,'(es16.9,2x))'
          do ie=1,rows
            ! Filling variables used in writing routing
            hw_count = ie
            mvec_spherical(:,2) = mangles(ie,:,1)
            mvec_spherical(:,3) = mangles(ie,:,2)
            write(dc_fields(hw_count),fmt=formatvar) (x(ie,k),k=1,idc)

            currents(:,:,:) = currents_from_file(:,:,:,ie)
            total_currents(:,:) = total_currents_from_file(:,:,ie)

            call calculate_lgtv_currents()
            call write_dc_lgtv_currents()
          end do

          hw_count = 1
          ! Closing files
          call openclose_dc_currents_files(2)
          call openclose_dc_lgtv_files(2)

          ! Sorting files
          call sort_lgtv_currents()

          ! Deallocating variables
          deallocate(data,x,currents_from_file,total_currents_from_file,mangles)
        end do energy_loop
      call deallocate_currents()
      call deallocate_lgtv_currents()
      end if
    case default
      write(outputunit,"('[read_calculate_lgtv_currents] Not possible to calculate longitudinal and transverse currents! itype = ',i0)") itype
      call MPI_Finalize(ierr)
      stop
    end select
    write(outputunit,"('[read_calculate_lgtv_currents] Finished calculating longitudinal and transverse currents!')")

    return
  end subroutine read_calculate_lgtv_currents

  subroutine calculate_lgtv_currents()
    use mod_parameters
    use mod_constants, only: zero
    use mod_currents
    use mod_lattice
    implicit none
    integer      :: i,j,neighbor

    long_currents = zero
    total_long_currents = zero
    transv_currents = zero
    total_transv_currents = zero
    do neighbor=n0sc1,n0sc2
      ! Summing currents flowing on longitudinal direction
      if(any(neighbor==sha_longitudinal(1:longitudinal_neighbors))) then
        do j=1,7
          do i=1,Npl
            long_currents(j,i) = long_currents(j,i) + currents(j,neighbor,i)*long_cos(neighbor)
          end do
          total_long_currents(j) = total_long_currents(j) + total_currents(j,neighbor)*long_cos(neighbor)
        end do
      end if

      ! Summing currents flowing on transverse direction
      if(any(neighbor==sha_transverse(1:transverse_neighbors))) then
        do j=1,7
          do i=1,Npl
            transv_currents(j,i) = transv_currents(j,i) + currents(j,neighbor,i)*transv_cos(neighbor)
          end do
          total_transv_currents(j) = total_transv_currents(j) + total_currents(j,neighbor)*transv_cos(neighbor)
        end do
      end if
    end do

    return
  end subroutine calculate_lgtv_currents
end module mod_lgtv_currents
