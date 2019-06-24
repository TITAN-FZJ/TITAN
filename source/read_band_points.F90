! Read reciprocal points from file kbands
! IN UNITS OF 2pi/a0
subroutine read_band_points(kbands, a0, b1, b2, b3)
  use mod_f90_kind,   only: double
  use mod_parameters, only: bands, band_cnt, qbasis
  use mod_mpi_pars,   only: abortProgram
  use mod_tools,      only: itos, number_of_lines
  use mod_constants,  only: pi
  implicit none
  real(double), dimension(:,:), allocatable, intent(out) :: kbands
  real(double), dimension(3),                intent(in)  :: b1, b2, b3
  real(double),                              intent(in)  :: a0
  character(len=40)  :: band_file = "kbands"
  integer :: file_unit = 666999

  character(len=200) :: line
  integer :: ios, line_count, i, j
  logical :: found
  type :: band_point
    character(len=5) :: name
    real(double), dimension(3) :: kp
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
          kbands(1:3,i) = (kband(j)%kp(1) * b1 + kband(j)%kp(2) * b2 + kband(j)%kp(3) * b3)*pi/a0
        else if (qbasis(1:1) == "c") then
          kbands(1:3,i) = kband(j)%kp(:)*pi/a0
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