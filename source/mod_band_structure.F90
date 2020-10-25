module mod_band_structure
  implicit none

contains

  ! Read reciprocal points from file kbands
  ! IN UNITS OF 2pi/a0
  subroutine read_band_points(kbands, a0, b1, b2, b3)
    use mod_kind,       only: dp
    use mod_parameters, only: bands,band_cnt,qbasis
    use mod_mpi_pars,   only: abortProgram
    use mod_tools,      only: itos,number_of_lines
    use mod_constants,  only: tpi
    implicit none
    real(dp), dimension(:,:), allocatable, intent(out) :: kbands
    real(dp), dimension(3),                intent(in)  :: b1, b2, b3
    real(dp),                              intent(in)  :: a0
    character(len=40)  :: band_file = "kbands"
    integer :: file_unit = 666999

    character(len=200) :: line
    integer :: ios, line_count, i, j
    logical :: found
    type :: band_point
      character(len=5) :: name
      real(dp), dimension(3) :: kp
    end type
    type(band_point), dimension(:), allocatable :: kband

    ! Opening file kbands
    open(unit = file_unit, file = trim(band_file), status='old', iostat=ios)
    if(ios /= 0) call abortProgram("[read_band_points] File 'kbands' not found!")

    ! Count non commented lines
    call number_of_lines(file_unit,i,line_count)

    ! Reading point name and position from file kbands
    allocate(kband(line_count))
    i = 0
    ios=0
    do while (ios == 0)
      line = ""
      read(unit=file_unit, fmt='(A)', iostat=ios) line
      if(len_trim(line) <= 0 .or. len_trim(line) >= 200 .or. index(adjustl(line),"#") == 1) cycle
      if(index(line, "#") > 0) line = line(1:index(line, "#")-1)
      i = i + 1
      read(unit = line, fmt = *, iostat=ios) kband(i)%name, (kband(i)%kp(j), j = 1,3)
      ! read(unit=file_unit, fmt='(A)', iostat=ios) line
    end do

    ! Storing path points defined in input
    allocate(kbands(3,band_cnt))
    do i = 1, band_cnt
      found = .false.
      do j = 1, line_count
        if(trim(kband(j)%name) == trim(bands(i))) then
          found = .true.
          if(qbasis(1:1) == "b") then
            kbands(1:3,i) = (kband(j)%kp(1) * b1 + kband(j)%kp(2) * b2 + kband(j)%kp(3) * b3)!*tpi/a0
          else if (qbasis(1:1) == "c") then
            kbands(1:3,i) = kband(j)%kp(:)*tpi/a0
          else
            call abortProgram("[read_band_points] Basis " // trim(qbasis) //" not defined!")
          end if
          exit
        endif
      end do
      if(.not. found) then
        call abortProgram("[read_band_points] No point " // trim(bands(i)) //" found!")
      endif
    end do
    close(file_unit)
  end subroutine read_band_points


  !   Calculates band structure
  subroutine band_structure(s)
    use mod_kind,              only: dp
    use mod_parameters,        only: output,dimHsc,kdirection,nQvec,nQvec1,kpoints,bsfile,wsfile,deltak
    use mod_system,            only: System_type
    use mod_tools,             only: cross,diagonalize,lwork
    use mod_mpi_pars,          only: rField
    use mod_io,                only: write_header
    use mod_hamiltonian,       only: hamilt_local,hamiltk
    implicit none
    type(System_type), intent(in) :: s
    integer :: i,kount,f_unit=666,w_unit=667,n
    integer :: ilaenv
    real(dp), dimension(:), allocatable :: eval
    complex(dp), allocatable :: hk(:,:)
    character(len=30) :: formatvar1,formatvar2

    external :: zheev,ilaenv
    
    if(rField == 0) write(output%unit_loop,"('CALCULATING THE BAND STRUCTURE')")

    allocate( hk(dimHsc,dimHsc),eval(dimHsc) )

    ! Getting lwork for diagonalization
    lwork = (ilaenv( 1, 'zhetrd', 'VU', dimHsc, -1, -1, -1 )+1)*dimHsc

    call hamilt_local(s)

    ! Band structure file
    write(bsfile,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=f_unit, file=bsfile, status='replace')
    call write_header(f_unit,"# dble((kount-1._dp)*deltak), (eval(i),i=1,dimHsc)",s%Ef)

    write(wsfile,"('./results/',a1,'SOC/',a,'/BS/weights_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=w_unit, file=wsfile, status='replace')
    call write_header(w_unit,"# ((abs(hk(i,n))**2,i=1,dimHsc),n=1,dimHsc)",s%Ef)

    write(unit=f_unit, fmt="(a,2x,i3)") "# dimHsc ",dimHsc
    write(formatvar1,fmt="(a,i0,a)") '(',1+dimHsc,'(es16.8e3,2x))'
    write(formatvar2,fmt="(a,i0,a)") '(',dimHsc*dimHsc,'(es10.3e3,2x))'

    do kount=1,nQvec1
      write(output%unit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") kount,nQvec1,dble((kount-1._dp)/nQvec)

      ! Calculating the hamiltonian for a given k-point
      call hamiltk(s,kpoints(:,kount),hk)
      ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
      call diagonalize(dimHsc,hk,eval)

      write(unit=f_unit,fmt=formatvar1) dble((kount-1._dp)*deltak), (eval(i),i=1,dimHsc)
      write(unit=w_unit,fmt=formatvar2) ((abs(hk(i,n))**2,i=1,dimHsc),n=1,dimHsc)
    end do

    close(f_unit)
    close(w_unit)
    deallocate(hk,eval)

  end subroutine band_structure
end module mod_band_structure
