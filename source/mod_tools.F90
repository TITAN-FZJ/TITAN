module mod_tools
  implicit none

  interface itos
    module procedure i4tos, &
                      i8tos
  end interface itos

  interface sort
    module procedure sort_int, &
                      sort_double_int
  end interface sort

  interface vec_norm
    module procedure vec_norm_real, &
                      vec_norm_complex
  end interface vec_norm

  interface normalize
    module procedure normalize_real, &
                      normalize_complex
  end interface normalize
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

  ! --------------------------------------------------------------------
  ! double precision function vecDist():
  !    Calculating the distance of two 3D points a and b
  ! --------------------------------------------------------------------
  function vecDist(a,b)
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
  ! subroutine sort_int():
  !    This subroutine receives an array x() and returns an integer array
  !  'order' with the positions of ascending numbers of x.
  ! --------------------------------------------------------------------
  subroutine sort_int(x, size, order)
    use mod_f90_kind, only: double
    implicit none
    integer,                       intent(in)  :: size
    real(double), dimension(size), intent(in)  :: x
    integer     , dimension(size), intent(out) :: order
    integer                                    :: i
    logical,      dimension(:), allocatable    :: mask

    allocate( mask(size) )

    mask = .true.
    do i = 1,size
      order(i) = minloc( x, 1, mask(:) )
      mask(order(i)) = .false.
    end do

    deallocate( mask )

  end subroutine sort_int

  ! --------------------------------------------------------------------
  ! subroutine sort_double_int():
  !    This subroutine receives an array x() and returns an double integer array
  !  'order' with the positions of ascending numbers of x.
  ! --------------------------------------------------------------------
  subroutine sort_double_int(x, size, order)
    use mod_f90_kind, only: double
    implicit none
    integer*8,                     intent(in)  :: size
    real(double), dimension(size), intent(in)  :: x
    integer*8   , dimension(size), intent(out) :: order
    integer*8                                  :: i
    logical,      dimension(:), allocatable    :: mask

    allocate( mask(size) )

    mask(:) = .true.
    do i = 1,size
      order(i) = minloc( x, 1, mask(:) )
      mask(order(i)) = .false.
    end do

    deallocate( mask )

  end subroutine sort_double_int


  ! --------------------------------------------------------------------
  ! subroutine number_of_lines():
  !    this subroutine returns the number of lines (commented and
  ! not commented) on a file (given by "unit", already opened).
  ! Notes: - Blank lines are ignored.
  !        - automatically rewinds the file
  ! --------------------------------------------------------------------
  subroutine number_of_lines(unit,total_lines,non_commented)
    implicit none
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

    rewind unit
  end subroutine number_of_lines


  ! --------------------------------------------------------------------
  ! subroutine number_of_rows_cols():
  !    this subroutine returns the number of rows and cols of a data file
  ! (given by "unit", already opened), as well as the number of commented 
  ! rows
  ! --------------------------------------------------------------------
  subroutine number_of_rows_cols(unit,rows,cols,commented_rows)
    implicit none
    integer, intent(out) :: rows, cols, commented_rows
    integer, intent(in)  :: unit
    character(len=900)   :: stringtemp
    integer :: ios,i

    rewind unit
    ! Counting the number of lines
    commented_rows = 0
    rows  = 0
    do
      read (unit=unit,fmt='(A)',iostat=ios) stringtemp
      if (ios/=0) exit
      ! Getting the number of rows
      if ((stringtemp(1:1)=="#").or.(stringtemp(1:1)=="!").or.(stringtemp=="")) then
        commented_rows = commented_rows + 1
        cycle
      end if
      rows = rows + 1
    end do
    cols = count([( stringtemp(i:i), i=1,len(stringtemp) )] == "E")

  end subroutine number_of_rows_cols


  ! --------------------------------------------------------------------
  ! subroutine read_data():
  !    this subroutine reads a table of data (size rows,cols) from a file
  ! (given by "unit", already opened) and returns it on "data".
  ! Blank lines are ignored.
  ! --------------------------------------------------------------------
  subroutine read_data(unit,rows,cols,data)
    use mod_f90_kind, only: double
    use mod_mpi_pars, only: abortProgram
    implicit none
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
      if (ios/=0)  call abortProgram("[read_data] Incorrect number of cols: " // trim(itos(j)) // " when expecting " // trim(itos(cols)))
    end do

    if(i/=rows) call abortProgram("[read_data] Incorrect number of rows: " // trim(itos(i)) // " when expecting " // trim(itos(rows)))
    ! Writing data
    ! do i=1,rows
    !   write(*,"(10(es16.9,2x))") (data(i,j),j=1,cols)
    ! end do

  end subroutine read_data


  ! --------------------------------------------------------------------
  ! subroutine read_data_with_comments():
  !    this subroutine reads a table of data (size rows,cols) from a file
  ! (given by "unit", already opened) and returns it on "data".
  ! The comments are returned in the variable "comments".
  ! "Mask" is an array of size total number of rows with .true. where data
  ! is, and .false. for the comments
  ! --------------------------------------------------------------------
  subroutine read_data_with_comments(unit,rows,cols,commented_rows,data,comments,mask)
    use mod_f90_kind, only: double
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer     , intent(in)  :: unit,rows,cols,commented_rows
    real(double), intent(out) :: data(rows,cols)
    character(len=900), dimension(commented_rows),      intent(out) :: comments
    logical,            dimension(rows+commented_rows), intent(out) :: mask

    character(len=900)        :: stringtemp
    integer :: ios,i,j,k,l

    rewind unit
    i = 0
    k = 0
    l = 0
    do
      l = l+1
      read(unit=unit,fmt='(A)',iostat=ios) stringtemp
      if (ios/=0) exit
      if ((stringtemp(1:1)=="#").or.(stringtemp(1:1)=="!").or.(stringtemp=="")) then
        k = k + 1
        comments(k) = trim(stringtemp)
        mask(l) = .false.
        cycle
      end if
      i=i+1
      read(unit=stringtemp,fmt=*,iostat=ios) (data(i,j),j=1,cols)
      if (ios/=0)  call abortProgram("[read_data] Incorrect number of cols: " // trim(itos(j)) // " when expecting " // trim(itos(cols)))
      mask(l) = .true.
    end do

    if(i/=rows) call abortProgram("[read_data] Incorrect number of data rows: " // trim(itos(i)) // " when expecting " // trim(itos(rows)))
    if(k/=commented_rows) call abortProgram("[read_data] Incorrect number of commented rows: " // trim(itos(k)) // " when expecting " // trim(itos(commented_rows)))
    ! Writing data
    ! do i=1,rows
    !   write(*,"(10(es16.9,2x))") (data(i,j),j=1,cols)
    ! end do

  end subroutine read_data_with_comments

  ! --------------------------------------------------------------------
  ! subroutine sort_file():
  !    This subroutine sorts the lines containing data of file 'unit'
  ! --------------------------------------------------------------------
  subroutine sort_file(unit)
    use mod_f90_kind, only: double
    use mod_mpi_pars
    ! use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: unit
    character(len=50)   :: colformat
    integer             :: i,j,k,l,rows,cols,commented_rows
    character(len=900), allocatable :: comments(:)
    integer     ,       allocatable :: order(:)
    real(double),       allocatable :: data(:,:),x(:)
    logical,            allocatable :: mask(:)

    ! Obtaining number of rows and cols in the file
    call number_of_rows_cols(unit,rows,cols,commented_rows)

    ! Allocating variables to read and sort data
    allocate(data(rows,cols),x(rows),order(rows),mask(rows+commented_rows))

    ! Reading data and storing to variable 'data'
    if(commented_rows==0) then
      call read_data(unit,rows,cols,data)
      mask = .false.
    else
      allocate( comments(commented_rows) )
      call read_data_with_comments(unit,rows,cols,commented_rows,data,comments,mask)
    end if

    ! Sorting by the first column
    x(:) = data(:,1)
    call sort(x,rows,order)

    ! Restarting the file
    rewind unit
    ! Rewriting data, now sorted
    write(colformat,fmt="(a,i0,a)") '(',cols,'(es16.9,2x))'
    i = 0
    k = 0
    do l=1,rows+commented_rows
      if(mask(l)) then
        i = i + 1
        write(unit=unit,fmt=trim(colformat)) (data(order(i),j),j=1,cols)
      else
        k = k + 1
        write(unit=unit,fmt=*) trim(comments(k))
      end if
    end do
    if(i/=rows) call abortProgram("[sort_file] Incorrect number of data rows: " // trim(itos(i)) // " when expecting " // trim(itos(rows)))
    if(k/=commented_rows) call abortProgram("[sort_file] Incorrect number of commented rows: " // trim(itos(k)) // " when expecting " // trim(itos(commented_rows)))

    deallocate(data,x,order,mask)
    if(commented_rows/=0) deallocate(comments)

  end subroutine sort_file

  ! --------------------------------------------------------------------
  ! function ItoS():
  !    This function transforms an integer i into a character variable ItoS
  ! --------------------------------------------------------------------
  character(len=100) function I4toS(i)
    implicit none
    integer :: i
    write(I4toS, "(i0)") i
  end function I4toS

  character(len=100) function I8toS(i)
    implicit none
    integer*8 :: i
    write(I8toS, "(i0)") i
  end function I8toS

  ! --------------------------------------------------------------------
  ! function RtoS():
  !    This function transforms a real r into a character variable RtoS.
  ! It also cuts leading spaces (on the left)
  ! --------------------------------------------------------------------
  character(len=900) function RtoS(r,format)
    use mod_f90_kind, only: double
    implicit none
    real(double)     :: r
    character(len=*) :: format
    write(Rtos, fmt=format) r
    RtoS = adjustl(RtoS)
  end function RtoS


  ! --------------------------------------------------------------------
  ! subroutine KronProd():
  !    This subroutine calculates the outer product (Kronecker product)
  ! of two matrices A(nax,nay) and B(nbx,nby)
  ! --------------------------------------------------------------------
  subroutine KronProd(nax,nay,nbx,nby,A,B,AB)
    use mod_f90_kind, only: double
    implicit none
    integer,         intent(in)  :: nax, nay, nbx, nby
    complex(double), intent(in)  :: A(nax,nay), B(nbx,nby) 
    complex(double), intent(out) :: AB(nax*nbx,nay*nby) 
    integer                      :: i, j, p, q, l, m                                              
    do i = 1,nax
      do j = 1,nay
        l=(i-1)*nbx + 1
        m=l+nbx-1
        p=(j-1)*nby + 1
        q=p+nby-1
        AB(l:m,p:q) = A(i,j)*B
      end do
    end do
  end subroutine KronProd

  ! --------------------------------------------------------------------
  ! function vec_norm_complex():
  !    This function calculates the norm of a complex vector v of dimension
  ! dim_v.
  ! --------------------------------------------------------------------
  function vec_norm_complex(v, dim_v)
    use mod_f90_kind, only: double
    implicit none
    integer,                           intent(in) :: dim_v ! vector dimension
    complex(double), dimension(dim_v), intent(in) :: v ! vector v
    real(double)                                  :: vec_norm_complex
    real(double)                                  :: sum
    integer                                       :: i

    sum= 0.d0
    do i = 1, dim_v
      sum = sum + v(i)*conjg(v(i))
    end do
    vec_norm_complex= sqrt(sum)
  end function vec_norm_complex


  ! --------------------------------------------------------------------
  ! function vec_norm_real():
  !    This function calculates the norm of a real vector v of dimension
  ! dim_v.
  ! --------------------------------------------------------------------
  function vec_norm_real(v, dim_v)
    use mod_f90_kind, only: double
    implicit none
    integer,                        intent(in) :: dim_v ! vector dimension
    real(double), dimension(dim_v), intent(in) :: v ! vector v
    real(double)                               :: vec_norm_real

    vec_norm_real = sqrt(dot_product(v,v))
  end function vec_norm_real


  ! --------------------------------------------------------------------
  ! function normalize_complex():
  !    This function normalizes the complex vector v of dimension
  ! dim_v.
  ! --------------------------------------------------------------------
  function normalize_complex(v, dim_v)
    use mod_f90_kind, only: double
    implicit none
    integer,                           intent(in) :: dim_v ! vector dimension
    complex(double), dimension(dim_v), intent(in) :: v ! vector v
    complex(double), dimension(dim_v)             :: normalize_complex

    normalize_complex = v/vec_norm(v,dim_v)
  end function normalize_complex


  ! --------------------------------------------------------------------
  ! function normalize_real():
  !    This function normalizes the real vector v of dimension
  ! dim_v.
  ! --------------------------------------------------------------------
  function normalize_real(v, dim_v)
    use mod_f90_kind, only: double
    implicit none
    integer,                        intent(in) :: dim_v ! vector dimension
    real(double), dimension(dim_v), intent(in) :: v ! vector v
    real(double), dimension(dim_v)             :: normalize_real

    normalize_real = v/vec_norm(v,dim_v)
  end function normalize_real


  ! --------------------------------------------------------------------
  ! subroutine LS_solver():
  !    This subruotine is a simple interface for the LAPACK linear
  ! system solver A*X = B.
  ! --------------------------------------------------------------------
  subroutine LS_solver(n,A,b)
    use mod_f90_kind, only: double
    implicit none
    
    integer,         intent(in)      :: n
    complex(double), intent(in)      :: A(n,n)
    complex(double), intent(out)     :: b(n)
    ! Workspace variables
    integer                          :: nrhs, lda, ldb, ldx, info
    integer,          allocatable    :: ipiv(:)
    complex,          allocatable    :: swork(:)
    complex(double),  allocatable    :: work(:), X(:)
    double precision, allocatable    :: rwork(:)

    lda= n
    ldb=n
    ldx= n
    nrhs=1

    allocate(rwork(n),swork(n*(n+1)),work(n), ipiv(n), X(n) )
    ! call zcgesv( n, nrhs, a, lda, ipiv, b, ldb, X, ldx, work, swork, rwork, iter, info )
    ! b = X
    call zgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
    
    if (info /= 0) write(*,'("eigensolver: failure in zheev")')
    deallocate(work,swork,rwork,ipiv,X)
    
    end subroutine LS_solver

end module mod_tools


