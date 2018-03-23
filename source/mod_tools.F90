module mod_tools
  implicit none

contains

  ! --------------------------------------------------------------------
  ! double precision function cross():
  !    This subroutine calculates the cross product of arrays a and b
  ! --------------------------------------------------------------------
  function cross(a, b)
    use mod_f90_kind, only: double
    implicit none
    real(double), dimension(3) :: cross
    real(double), dimension(3), intent(in) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
      end function cross

  function vecDist(a,b)
  !! Function calculating the distance of two 3D points
     use mod_f90_kind, only: double
     implicit none
     real(double) :: vecDist
     real(double), dimension(3), intent(in) :: a, b

     vecDist = sqrt(dot_product(a - b, a - b))
       end function vecDist

  ! --------------------------------------------------------------------
  ! logical function is_perpendicular(a,b):
  !  This subroutines return whether a and b are perpendicular or not
  ! --------------------------------------------------------------------

  function is_parallel(a,b)
    use mod_f90_kind, only: double
    implicit none
    logical :: is_parallel
    real(double), dimension(3) :: crs
    real(double), dimension(3), intent(in) :: a, b

    is_parallel = .false.
    crs = cross(a,b)
    if( 1d-9 > dot_product(crs(1:3),crs(1:3)) ) then
       is_parallel = .true.
    end if
      end function is_parallel

  ! --------------------------------------------------------------------
  ! logical function is_parallel(a,b):
  !  This subroutines return whether a and b are parallel or not
  ! --------------------------------------------------------------------

  function is_perpendicular(a,b)
    use mod_f90_kind, only: double
    implicit none
    logical :: is_perpendicular
    real(double), dimension(3), intent(in) :: a, b

    is_perpendicular = .false.
    if( 1.d-9 > abs(dot_product(a,b)) ) then
       is_perpendicular = .true.
    end if
      end function is_perpendicular

  ! --------------------------------------------------------------------
  ! double precision function cross_unit():
  !    This subroutine calculates the unit vector in the direction
  ! of the cross product of arrays a and b
  ! --------------------------------------------------------------------
  function cross_unit(a, b)
    use mod_f90_kind, only: double
    implicit none
    real(double), dimension(3) :: cross_unit
    real(double), dimension(3), intent(in) :: a, b

    cross_unit(1) = a(2) * b(3) - a(3) * b(2)
    cross_unit(2) = a(3) * b(1) - a(1) * b(3)
    cross_unit(3) = a(1) * b(2) - a(2) * b(1)

    cross_unit = cross_unit/sqrt(dot_product(cross_unit,cross_unit))
      end function cross_unit

  ! --------------------------------------------------------------------
  ! subroutine Sort():
  !    This subroutine receives an array x() and returns an integer array
  !  'order' with the positions of ascending numbers of x.
  ! --------------------------------------------------------------------
  subroutine sort(x, size, order)
    use mod_f90_kind, only: double
    implicit  none
    integer, intent(in)                        :: size
    real(double), dimension(size), intent(in)  :: x
    integer     , dimension(size), intent(out) :: order
    integer                                    :: i
    logical, dimension(size) :: mask

    mask = .true.
    do i = 1,size
       order(i) = minloc(x,1,mask)
       mask(order(i)) = .false.
    end do

      end subroutine sort

  ! --------------------------------------------------------------------
  ! subroutine number_of_lines():
  !    this subroutine returns the number of lines (commented and
  ! not commented) on a file (given by "unit", already opened).
  ! Blank lines are ignored.
  ! --------------------------------------------------------------------
  subroutine number_of_lines(unit,total_lines,non_commented)
    implicit  none
    integer, intent(out) :: total_lines, non_commented
    integer, intent(in)  :: unit
    character(len=20)    :: stringtemp
    integer :: ios

    rewind unit

    ! Counting the number of lines
    non_commented  = 0
    total_lines    = 0
    do
       read (unit=unit,fmt=*,iostat=ios) stringtemp
       if (ios/=0) exit
       if (stringtemp=="") cycle ! If the line is blank, ignore
       ! Total number of non-empty lines
       total_lines = total_lines + 1

       ! Getting the number of non-commented lines
       if ((stringtemp(1:1)=="#").or.(stringtemp(1:1)=="!")) cycle
       non_commented = non_commented + 1
    end do

      end subroutine number_of_lines


  ! --------------------------------------------------------------------
  ! subroutine number_of_rows_cols():
  !    this subroutine returns the number of rows and cols of a data file
  ! (given by "unit", already opened).
  ! --------------------------------------------------------------------
  subroutine number_of_rows_cols(unit,rows,cols)
    implicit  none
    integer, intent(out) :: rows, cols
    integer, intent(in)  :: unit
    character(len=900)   :: stringtemp
    integer :: ios,i

    rewind unit
    ! Counting the number of lines
    rows  = 0
    do
       read (unit=unit,fmt='(A)',iostat=ios) stringtemp
       if (ios/=0) exit
       ! Getting the number of rows
       if ((stringtemp(1:1)=="#").or.(stringtemp(1:1)=="!").or.(stringtemp=="")) cycle
       rows = rows + 1
    end do
    cols = count([( stringtemp(i:i), i=1,len(stringtemp) )] == "E")

      end subroutine number_of_rows_cols


  ! --------------------------------------------------------------------
  ! subroutine number_of_lines():
  !    this subroutine returns the number of lines (commented and
  ! not commented) on a file (given by "unit", already opened).
  ! Blank lines are ignored.
  ! --------------------------------------------------------------------
  subroutine read_data(unit,rows,cols,data)
    use mod_f90_kind, only: double
    implicit  none
    integer     , intent(in)  :: unit,rows,cols
    real(double), intent(out) :: data(rows,cols)
    character(len=900)        :: stringtemp
    integer :: ios,i,j

    rewind unit
    i = 0
    do
       read(unit=unit,fmt='(A)',iostat=ios) stringtemp
       if (ios/=0) exit
       if ((stringtemp(1:1)=="#").or.(stringtemp(1:1)=="!").or.(stringtemp=="")) cycle
       i=i+1
       read(unit=stringtemp,fmt=*,iostat=ios) (data(i,j),j=1,cols)
    end do
    ! Writing data
    ! do i=1,rows
    !   write(*,"(10(es16.9,2x))") (data(i,j),j=1,cols)
    ! end do

      end subroutine read_data

  ! --------------------------------------------------------------------
  ! subroutine sort_file():
  !    This subroutine sorts the lines of file 'unit'
  ! option to skip first line: header = .true.
  ! --------------------------------------------------------------------
  subroutine sort_file(unit,header)
    use mod_f90_kind, only: double
    implicit  none
    integer, intent(in) :: unit
    logical, intent(in) :: header
    character(len=50)   :: colformat
    integer             :: i,j,rows,cols,ios
    integer     , allocatable :: order(:)
    real(double), allocatable :: data(:,:),x(:)

    ! Obtaining number of rows and cols in the file
    call number_of_rows_cols(unit,rows,cols)

    ! Allocating variables to read and sort data
    allocate(data(rows,cols),x(rows),order(rows))

    ! Reading data and storing to variable 'data'
    call read_data(unit,rows,cols,data)

    ! Sorting by the first column
    x(:) = data(:,1)
    call sort(x,rows,order)

    ! Restarting the file and optionally skipping first line (header)
    rewind unit
    if(header) read(unit=unit,fmt='(A)',iostat=ios) colformat

    ! Rewriting data, now sorted
    write(colformat,fmt="(a,i0,a)") '(',cols,'(es16.9,2x))'
    do i=1,rows
       write(unit=unit,fmt=trim(colformat)) (data(order(i),j),j=1,cols)
    end do

    deallocate(data,x,order)

      end subroutine sort_file

  character(len=900) function ItoS(i)
     implicit none
     integer :: i
     write(Itos, "(i0)") i
       end function ItoS

end module mod_tools
