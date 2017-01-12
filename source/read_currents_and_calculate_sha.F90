! This subroutine reads currents from existing files
! and calculates the spin hall angle
subroutine read_currents_and_calculate_sha()
  use mod_parameters
  use mod_constants
  use mod_magnet
  use mod_progress
  use mod_currents
  use mod_sha
  use mod_mpi_pars
  use mod_lattice
  use mod_tools
  implicit none
  character(len=50) :: formatvar
  integer           :: i,j,k,iw,neighbor,rows,cols,idc,ie,iemin
  real(double)      :: idia
  real(double)   , allocatable :: data(:,:),x(:,:),e(:)
  complex(double), allocatable :: currents_from_file(:,:,:,:),total_currents_from_file(:,:,:)

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
        if(sum(abs(e(:)-data(:,1))).gt.1.d-8) then
          write(outputunit,"('[read_currents_and_calculate_sha] Different energies on current files!')")
          call MPI_Finalize(ierr)
          stop
        end if
      end if
      ! Obtaining diamagnetic current (if the minimum energy is small enough)
      if(e(iemin).le.1.d-5) then
        idia = data(iemin,4)
        data(:,4) = data(:,4) - idia
      end if

      currents_from_file(j,neighbor,i,:) = cmplx(data(:,3),data(:,4))
    end do ; end do ; end do

    ! Reading total currents for each neighbor direction
    do neighbor=n0sc1,n0sc2 ; do j=1,7
      iw = 7000+(neighbor-1)*7+j

      ! Reading data and storing to variable 'data'
      call read_data(iw,rows,cols,data)

      ! Checking if energy list is the same
      if(sum(abs(e(:)-data(:,1))).gt.1.d-8) then
        write(outputunit,"('[read_currents_and_calculate_sha] Different energies on current files!')")
        call MPI_Finalize(ierr)
        stop
      end if
      ! Obtaining diamagnetic current (if the minimum energy is small enough)
      if(e(iemin).le.1.d-5) then
        idia = data(iemin,4)
        data(:,4) = data(:,4) - idia
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
    if(hw_count.eq.1) then
      energy_loop: do count=1,npt1

        ! Allocating current and SHA variables
        call allocate_currents()
        call allocate_sha()

        ! Creating SHA files with headers
        call openclose_dc_sha_files(0)

        ! Opening current and SHA files
        call openclose_dc_sha_files(1)
        call openclose_dc_currents_files(1)

        idc = ceiling((dble(dcfield_dependence)-0.01d0)/3.d0)
        write(formatvar,fmt="(a,i0,a)") '(',idc+4,'(es16.9,2x))'
        ! Currents per plane per neighbor
        do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
          iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7+j

          if(.not.allocated(data)) then ! Gettingn information from first file
            ! Obtaining number of rows and cols in the file
            call number_of_rows_cols(iw,rows,cols)
            ! Allocating variables to read and sort data
            allocate(data(rows,cols),x(rows,idc),currents_from_file(7,n0sc1:n0sc2,Npl,rows),total_currents_from_file(7,n0sc1:n0sc2,rows))
            ! Reading data and storing to variable 'data'
            call read_data(iw,rows,cols,data)
            ! Getting fields list
            do k=1,idc
              x(:,k) = data(:,k)
            end do

          else
            ! Reading data and storing to variable 'data'
            call read_data(iw,rows,cols,data)

            ! Checking if fields list is the same
            if(sum(abs(x(:,1:idc)-data(:,1:idc))).gt.1.d-8) then
              write(outputunit,"('[read_currents_and_calculate_sha] Different fields on current files!')")
              call MPI_Finalize(ierr)
              stop
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
          if(sum(abs(x(:,1:idc)-data(:,1:idc))).gt.1.d-8) then
            write(outputunit,"('[read_currents_and_calculate_sha] Different fields on current files!')")
            call MPI_Finalize(ierr)
            stop
          end if

          total_currents_from_file(j,neighbor,:) = cmplx(data(:,1+idc),data(:,2+idc))
        end do ; end do


        ! Calculating SHA for all fields
        do ie=1,rows
          currents(:,:,:) = currents_from_file(:,:,:,ie)
          total_currents(:,:) = total_currents_from_file(:,:,ie)

          call calculate_sha()

          do j=2,4
            do i=1,Npl
              iw = 82000+(i-1)*3+j-1
              write(unit=iw,fmt=formatvar) (x(ie,k),k=1,idc) , sha_re(j,i) , abs(sha_complex(j,i)) , real(sha_complex(j,i)) , aimag(sha_complex(j,i))
            end do
            iw = 82000+Npl*3+j-1
            write(unit=iw,fmt=formatvar) (x(ie,k),k=1,idc) , sha_re_total(j) , abs(sha_complex_total(j)) , real(sha_complex_total(j)) , aimag(sha_complex_total(j))
          end do
        end do

        ! Closing files
        call openclose_dc_currents_files(2)
        call openclose_dc_sha_files(2)
      end do energy_loop
    end if
    ! Deallocating variables
    deallocate(data,x,currents_from_file,total_currents_from_file)
    call deallocate_currents()
    call deallocate_sha()
  case default
    write(outputunit,"('[read_currents_and_calculate_sha] Not possible to calculate spin Hall angle! itype = ',i0)") itype
    call MPI_Finalize(ierr)
    stop
  end select
  write(outputunit,"('[read_currents_and_calculate_sha] Finished calculating SHA!')")

  return
end subroutine read_currents_and_calculate_sha
