module mod_sha
  !! Spin Hall Angle variables and procedures
  use, intrinsic :: iso_fortran_env
  implicit none
  ! Spin Hall Angle
  real(dp), allocatable :: sha_re(:,:)
  real(dp) :: sha_re_total(4)
  complex(dp), allocatable :: sha_complex(:,:)
  complex(dp) :: sha_complex_total(4)

contains

subroutine allocate_sha()
  !! This subroutine allocates variables related to the sha calculation
    use, intrinsic :: iso_fortran_env
    use mod_mpi_pars
    use mod_parameters, only: Npl,outputunit
    implicit none
    integer           :: AllocateStatus

    if(myrank_row==0) then
      allocate( sha_re(4,Npl),sha_complex(4,Npl), STAT = AllocateStatus )
      if (AllocateStatus/=0) then
        write(outputunit,"('[allocate_sha] Not enough memory for: sha_re,sha_complex')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if

  end subroutine allocate_sha

  subroutine deallocate_sha()
    !! This subroutine allocates variables related to the sha calculation
    use, intrinsic :: iso_fortran_env
    use mod_mpi_pars
    implicit none

    if(myrank_row==0) deallocate(sha_re,sha_complex)

  end subroutine deallocate_sha

  subroutine openclose_sha_files(iflag)
    !! This subroutine opens and closes all the files needed for the sha
    use mod_parameters, only: fieldpart,output
    use mod_SOC, only: SOCc, socpart
    use mod_mpi_pars
    use mod_system, only: nkpt
    use electricfield, only: strElectricField
    implicit none

    character(len=500)  :: varm
    character(len=6)    :: folder(8)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,iflag,err,errt=0

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

    if(iflag==0) then
      ! Header for SHA
      do j=2,4
        do i=1,Npl
          iw = 8200+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),filename(j),i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#     energy     ,  SHA reIs/reIc  ,  SHA abs(Is/Ic) , SHA atan(Is/Ic) ,  SHA re(Is/Ic)  ,  SHA im(Is/Ic)  ')")
          close(unit=iw)
        end do
        iw = 8200+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),filename(j),trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy     ,  SHA reIs/reIc  ,  SHA abs(Is/Ic) , SHA atan(Is/Ic) ,  SHA re(Is/Ic)  ,  SHA im(Is/Ic)  ')")
        close(unit=iw)
      end do

    else if(iflag==1) then
      ! SHA
      do j=2,4
        do i=1,Npl
          iw = 8200+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),filename(j),i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end do
        iw = 8200+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/SHA_',a,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),filename(j),trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
      end do

      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_sha_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
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

  end subroutine openclose_sha_files

  subroutine write_sha(e)
  !! This subroutine write all the SHA into files
  !! (already opened with openclose_sha_files(1))
    use, intrinsic :: iso_fortran_env
    use mod_parameters, only: Npl
    use mod_magnet, only: mvec_spherical,mtotal_spherical
    implicit none
    integer  :: i,j,iw
    real(dp),intent(in) :: e

    do j=2,4
      do i=1,Npl
        iw = 8200+(i-1)*3+j-1
        write(unit=iw,fmt="(8(es16.9,2x))") e , sha_re(j,i) , abs(sha_complex(j,i)) , atan2(aimag(sha_complex(j,i)),real(sha_complex(j,i))) , real(sha_complex(j,i)) , aimag(sha_complex(j,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
      end do
      iw = 8200+Npl*3+j-1
      write(unit=iw,fmt="(8(es16.9,2x))") e , sha_re_total(j) , abs(sha_complex_total(j)) , atan2(aimag(sha_complex_total(j)),real(sha_complex_total(j))) , real(sha_complex_total(j)) , aimag(sha_complex_total(j)) , mtotal_spherical(2) , mtotal_spherical(3)
    end do

  end subroutine write_sha

  subroutine openclose_dc_sha_files(iflag)
    !! This subroutine opens and closes all the files needed for the sha
    use mod_parameters, only: dcfieldpart,output
    use mod_SOC, only: SOCc, socpart
    use mod_mpi_pars, only: errorcode, ierr, MPI_COMM_WORLD
    use mod_system, only: nkpt
    use electricfield, only: strElectricField
    implicit none

    character(len=500)  :: varm
    character(len=6)    :: folder(8)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,iflag,err,errt=0

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

    if(iflag==0) then
      ! Header for SHA
      do j=2,4
        do i=1,Npl
          iw = 82000+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#',a,'    SHA reIs/reIc   ,    SHA abs(Is/Ic)   ,    SHA atan(Is/Ic)   ,    SHA re(Is/Ic)   ,    SHA im(Is/Ic)   ')") trim(dc_header)
          close(unit=iw)
        end do
        iw = 82000+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,'    SHA reIs/reIc   ,    SHA abs(Is/Ic)   ,    SHA atan(Is/Ic)   ,    SHA re(Is/Ic)   ,    SHA im(Is/Ic)   ')") trim(dc_header)
        close(unit=iw)
      end do

    else if(iflag==1) then
      ! SHA
      do j=2,4
        do i=1,Npl
          iw = 82000+(i-1)*3+j-1
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          if(.not.lsha) then
            errt = errt + err
            if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
          end if
        end do
        iw = 82000+Npl*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'SHA_',a,'_',a,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(8)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        if(.not.lsha) then
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end if
      end do

      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_dc_sha_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
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

  end subroutine openclose_dc_sha_files

  subroutine write_dc_sha()
  !! This subroutine write all the sha in the dc limit into files
  !! (already opened with openclose_dc_sha_files(1))
    use, intrinsic :: iso_fortran_env
    use mod_parameters
    use mod_magnet, only: mvec_spherical,mtotal_spherical
    implicit none
    integer  :: i,j,iw

    do j=2,4
      do i=1,Npl
        iw = 82000+(i-1)*3+j-1
        write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , sha_re(j,i) , abs(sha_complex(j,i)) , atan2(aimag(sha_complex(j,i)),real(sha_complex(j,i))) , real(sha_complex(j,i)) , aimag(sha_complex(j,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
      end do
      iw = 82000+Npl*3+j-1
      write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , sha_re_total(j) , abs(sha_complex_total(j)) , atan2(aimag(sha_complex_total(j)),real(sha_complex_total(j))) , real(sha_complex_total(j)) , aimag(sha_complex_total(j)) , mtotal_spherical(2) , mtotal_spherical(3)
    end do

  end subroutine write_dc_sha

  subroutine sort_sha()
    !! This subroutine sorts SHA files
    use, intrinsic :: iso_fortran_env
    use mod_parameters, only: Npl,itype
    use mod_tools, only: sort_file
    implicit none
    integer  :: i,j,iw,idc=1

    ! Opening SHA files
    if(itype==9) then
      idc=10
      call openclose_dc_sha_files(1)
    else
      call openclose_sha_files(1)
    end if

    ! SHA
    do j=2,4
      do i=1,Npl
        iw = 8200*idc+(i-1)*3+j-1
        call sort_file(iw)
      end do
      iw = 8200*idc+Npl*3+j-1
      call sort_file(iw)
    end do

    ! Closing SHA files
    if(itype==9) then
      call openclose_dc_sha_files(2)
    else
      call openclose_sha_files(2)
    end if

  end subroutine sort_sha

  subroutine read_currents_and_calculate_sha()
    !! This subroutine reads currents from existing files
    !! and calculates the spin hall angle
    use mod_parameters
    use mod_constants
    use mod_magnet
    use mod_progress
    use mod_currents
    use mod_mpi_pars
    use mod_system, only: n0sc1, n0sc2
    use mod_tools
    implicit none
    character(len=50) :: formatvar
    integer           :: i,j,k,iw,neighbor,rows,cols,idc,ie,iemin
    real(dp)      :: idia
    real(dp)   , allocatable :: data(:,:),x(:,:),e(:),mangles(:,:,:)
    complex(dp), allocatable :: currents_from_file(:,:,:,:),total_currents_from_file(:,:,:)

    write(outputunit,"('[read_currents_and_calculate_sha] Reading current files and calculating SHA...')")
    ! Opening disturbance files
    select case(itype)
    case(8)
      ! Allocating current and SHA variables
      call allocate_currents()
      call allocate_sha()

      ! Creating SHA files with headers
      call openclose_sha_files(0)

      ! Opening current and SHA files
      call openclose_sha_files(1)
      call openclose_currents_files(1)

      ! Reading currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j

        if(.not.allocated(data)) then ! Getting information from first file
          ! Obtaining number of rows and cols in the file
          call number_of_rows_cols(iw,rows,cols)
          ! Allocating variables to read and store data
          allocate( data (rows,cols), &
                    e (rows), &
                    currents_from_file       (7,n0sc1:n0sc2,Npl,rows), &
                    total_currents_from_file (7,n0sc1:n0sc2,rows))

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
          if(sum(abs(e(:)-data(:,1)))>1.e-8_dp) then
            write(outputunit,"('[read_currents_and_calculate_sha] Different energies on current files!')")
            call MPI_Finalize(ierr)
            stop
          end if
        end if
        ! Obtaining diamagnetic current (if the minimum energy is small enough)
        if(e(iemin)<=1.e-5_dp) then
          idia = data(iemin,4)
          data(:,4) = data(:,4) - idia
        else
          write(outputunit,"('[read_currents_and_calculate_sha] WARNING: Diamagnetic current was not removed!')")
        end if

        currents_from_file(j,neighbor,i,:) = cmplx(data(:,3),data(:,4))
      end do ; end do ; end do

      ! Reading total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 7000+(neighbor-1)*7+j

        ! Reading data and storing to variable 'data'
        call read_data(iw,rows,cols,data)

        ! Checking if energy list is the same
        if(sum(abs(e(:)-data(:,1)))>1.e-8_dp) then
          write(outputunit,"('[read_currents_and_calculate_sha] Different energies on current files!')")
          call MPI_Finalize(ierr)
          stop
        end if
        ! Obtaining diamagnetic current (if the minimum energy is small enough)
        if(e(iemin)<=1.e-5_dp) then
          idia = data(iemin,4)
          data(:,4) = data(:,4) - idia
        else
          write(outputunit,"('[read_currents_and_calculate_sha] WARNING: Diamagnetic current was not removed!')")
        end if

        total_currents_from_file(j,neighbor,:) = cmplx(data(:,3),data(:,4))
      end do ; end do

      ! Calculating SHA for all energies
      do ie=1,rows
        currents(:,:,:) = currents_from_file(:,:,:,ie)
        total_currents(:,:) = total_currents_from_file(:,:,ie)

        call calculate_sha()
        call write_sha(e(ie))
      end do

      ! Closing files
      call openclose_currents_files(2)
      call openclose_sha_files(2)
      ! Deallocating variables
      deallocate(data,e,currents_from_file,total_currents_from_file)
      call deallocate_currents()
      call deallocate_sha()
    case(9)
      if(hw_count==1) then

        ! Allocating current and SHA variables
        call allocate_currents()
        call allocate_sha()
        energy_loop: do kount=1,nEner1
          ! Creating SHA files with headers
          call openclose_dc_sha_files(0)

          ! Opening current and SHA files
          call openclose_dc_sha_files(1)
          call openclose_dc_currents_files(1)

          idc = ceiling((dble(dcfield_dependence)-0.01_dp)/3._dp)
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
              allocate( data (rows,cols), &
                        x (rows,idc), &
                        currents_from_file       (7,n0sc1:n0sc2,Npl,rows), &
                        total_currents_from_file (7,n0sc1:n0sc2,rows), &
                        mangles (rows,Npl,2))
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
              if(sum(abs(x(:,1:idc)-data(:,1:idc)))>1.e-8_dp) then
                write(outputunit,"('[read_currents_and_calculate_sha] Different fields on current files!')")
                call MPI_Finalize(ierr)
                stop
              end if

              ! Checking if angles are the same or getting next plane values
              if(sum(abs(mangles(:,i,:)))<10._dp) then
                if(sum(abs(mangles(:,i,:)-data(:,cols-1:cols)))>1.e-8_dp) then
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
            if(sum(abs(x(:,1:idc)-data(:,1:idc)))>1.e-8_dp) then
              write(outputunit,"('[read_currents_and_calculate_sha] Different fields on current files!')")
              call MPI_Finalize(ierr)
              stop
            end if

            ! Checking if angles are the same
            if(sum(abs(mangles(:,1,:)-data(:,cols-1:cols)))>1.e-8_dp) then
              write(outputunit,"('[read_calculate_lgtv_currents] Different angles on current files!')")
              call MPI_Finalize(ierr)
              stop
            end if

            total_currents_from_file(j,neighbor,:) = cmplx(data(:,1+idc),data(:,2+idc))
          end do ; end do

          ! Calculating SHA for all fields
          write(formatvar,fmt="(a,i0,a)") '(',idc,'(es16.9,2x))'
          do ie=1,rows
            ! Filling variables used in writing routing
            hw_count = ie
            mvec_spherical(2,:) = mangles(ie,:,1)
            mvec_spherical(3,:) = mangles(ie,:,2)
            write(dc_fields(hw_count),fmt=formatvar) (x(ie,k),k=1,idc)
            currents(:,:,:) = currents_from_file(:,:,:,ie)
            total_currents(:,:) = total_currents_from_file(:,:,ie)

            call calculate_sha()
            call write_dc_sha()
          end do

          hw_count = 1
          ! Closing files
          call openclose_dc_currents_files(2)
          call openclose_dc_sha_files(2)

          ! Deallocating variables
          deallocate(data,x,currents_from_file,total_currents_from_file,mangles)
        end do energy_loop

        call deallocate_currents()
        call deallocate_sha()
      end if
    case default
      write(outputunit,"('[read_currents_and_calculate_sha] Not possible to calculate spin Hall angle! itype = ',i0)") itype
      call MPI_Finalize(ierr)
      stop
    end select
    write(outputunit,"('[read_currents_and_calculate_sha] Finished calculating SHA!')")

  end subroutine read_currents_and_calculate_sha


  subroutine calculate_sha()
    !! This subroutine calculates the spin hall angle given by
    !! SHA = (2._dp*e/\hbar) I_S/I_C
    !! NOTE: THIS IS DEFINED AS A CURRENT RATIO, NOT CURRENT DENSITY!
    use mod_parameters
    use mod_constants, only: cZero
    use mod_currents
    use mod_system, only: n0sc1, n0sc2
    implicit none
    integer      :: i,j,neighbor

    sha_re = 0._dp
    sha_complex = cZero
    sha_re_total = 0._dp
    sha_complex_total = cZero
    do neighbor=n0sc1,n0sc2
      ! Summing charge currents flowing on longitudinal direction
      if(any(neighbor==sha_longitudinal(1:longitudinal_neighbors))) then
        do i=1,Npl
          sha_re(1,i) = sha_re(1,i) + real(currents(1,neighbor,i))*long_cos(neighbor)
          sha_complex(1,i) = sha_complex(1,i) + currents(1,neighbor,i)*long_cos(neighbor)
        end do
        sha_re_total(1) = sha_re_total(1) + real(total_currents(1,neighbor))*long_cos(neighbor)
        sha_complex_total(1) = sha_complex_total(1) + total_currents(1,neighbor)*long_cos(neighbor)
      end if

      ! Summing spin currents flowing on transverse direction
      if(any(neighbor==sha_transverse(1:transverse_neighbors))) then
        do j=2,4
          do i=1,Npl
            sha_re(j,i) = sha_re(j,i) + real(currents(j,neighbor,i))*transv_cos(neighbor)
            sha_complex(j,i) = sha_complex(j,i) + currents(j,neighbor,i)*transv_cos(neighbor)
          end do
          sha_re_total(j) = sha_re_total(j) + real(total_currents(j,neighbor))*transv_cos(neighbor)
          sha_complex_total(j) = sha_complex_total(j) + total_currents(j,neighbor)*transv_cos(neighbor)
        end do
      end if
    end do
    ! Calculating real and imaginary part of spin Hall angles
    do i=1,Npl
      sha_re(2:4,i) = 2._dp*sha_re(2:4,i)/sha_re(1,i)
      sha_complex(2:4,i) = 2._dp*sha_complex(2:4,i)/sha_complex(1,i)
    end do
    sha_re_total(2:4) = 2._dp*sha_re_total(2:4)/sha_re_total(1)
    sha_complex_total(2:4) = 2._dp*sha_complex_total(2:4)/sha_complex_total(1)

  end subroutine calculate_sha
end module mod_sha
