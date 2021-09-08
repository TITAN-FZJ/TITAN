module mod_tools
  implicit none

  integer :: lwork
  
  interface itos
    module procedure i4tos, &
                      i8tos
  end interface itos

  interface StoI
    module procedure StoI_scalar, &
                      StoI_array
  end interface StoI

  interface StoR
    module procedure StoR_scalar, &
                      StoR_array
  end interface StoR

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

  interface print_matrix
    module procedure print_matrix_real, &
                      print_matrix_complex
  end interface print_matrix

contains

  function cross(a, b)
  !> --------------------------------------------------------------------
  !> double precision function cross():
  !>    This subroutine calculates the cross product of arrays a and b
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    real(dp), dimension(3) :: cross
    real(dp), dimension(3), intent(in) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross

  function vecDist(a,b)
  !> --------------------------------------------------------------------
  !> double precision function vecDist():
  !>    Calculating the distance of two 3D points a and b
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    real(dp) :: vecDist
    real(dp), dimension(3), intent(in) :: a, b

    vecDist = sqrt(dot_product(a - b, a - b))
  end function vecDist

  function is_numeric(string)
  !> --------------------------------------------------------------------
  !> logical function is_numeric(string):
  !>  This subroutines return whether a string is numeric (.true.) or not (.false.)
  !> --------------------------------------------------------------------
    implicit none
    character(len=*), intent(in) :: string
    logical :: is_numeric
    real    :: x
    integer :: err
    read(string,*,iostat=err) x
    is_numeric = err == 0
  end function is_numeric


  function is_parallel(a,b)
  !> --------------------------------------------------------------------
  !> logical function is_perpendicular(a,b):
  !>  This subroutines return whether a and b are perpendicular or not
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    logical :: is_parallel
    real(dp), dimension(3) :: crs
    real(dp), dimension(3), intent(in) :: a, b

    is_parallel = .false.
    crs = cross(a,b)
    if( 1e-9_dp > dot_product(crs(1:3),crs(1:3)) ) then
      is_parallel = .true.
    end if
  end function is_parallel


  function is_perpendicular(a,b)
  !> --------------------------------------------------------------------
  !> logical function is_parallel(a,b):
  !>  This subroutines return whether a and b are parallel or not
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    logical :: is_perpendicular
    real(dp), dimension(3), intent(in) :: a, b

    is_perpendicular = .false.
    if( 1.e-9_dp > abs(dot_product(a,b)) ) then
      is_perpendicular = .true.
    end if
  end function is_perpendicular

  function cross_unit(a, b)
  !> --------------------------------------------------------------------
  !> double precision function cross_unit():
  !>    This subroutine calculates the unit vector in the direction
  !> of the cross product of arrays a and b
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    real(dp), dimension(3) :: cross_unit
    real(dp), dimension(3), intent(in) :: a, b

    cross_unit(1) = a(2) * b(3) - a(3) * b(2)
    cross_unit(2) = a(3) * b(1) - a(1) * b(3)
    cross_unit(3) = a(1) * b(2) - a(2) * b(1)

    cross_unit = cross_unit/sqrt(dot_product(cross_unit,cross_unit))
  end function cross_unit


  subroutine sort_int(x, isize, order)
  !> --------------------------------------------------------------------
  !> subroutine sort_int():
  !>    This subroutine receives an array x() and returns an integer array
  !>  'order' with the positions of ascending numbers of x.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,                    intent(in)  :: isize
    real(dp), dimension(isize), intent(in)  :: x
    integer,  dimension(isize), intent(out) :: order
    integer                                 :: i
    logical,  dimension(isize)              :: mask

    mask = .true.
    do i = 1,isize
      order(i) = minloc( x, 1, mask(:) )
      mask(order(i)) = .false.
    end do

  end subroutine sort_int

  subroutine sort_double_int(x, isize, order)
  !> --------------------------------------------------------------------
  !> subroutine sort_double_int():
  !>    This subroutine receives an array x() and returns an double integer array
  !>  'order' with the positions of ascending numbers of x.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp,int64
    implicit none
    integer(int64),                   intent(in)  :: isize
    real(dp),       dimension(isize), intent(in)  :: x
    integer(int64), dimension(isize), intent(out) :: order
    integer(int64)                                :: i
    logical,        dimension(isize)              :: mask

    mask(:) = .true.
    do i = 1,isize
      order(i) = minloc( x, 1, mask(:) )
      mask(order(i)) = .false.
    end do

  end subroutine sort_double_int


  subroutine number_of_lines(unit,total_lines,non_commented)
  !> --------------------------------------------------------------------
  !> subroutine number_of_lines():
  !>    this subroutine returns the number of lines (commented and
  !> not commented) on a file (given by "unit", already opened).
  !> Notes: - Blank lines are ignored.
  !>        - automatically rewinds the file
  !> --------------------------------------------------------------------
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


  subroutine number_of_rows_cols(unit,rows,cols,commented_rows)
  !> --------------------------------------------------------------------
  !> subroutine number_of_rows_cols():
  !>    this subroutine returns the number of rows and cols of a data file
  !> (given by "unit", already opened), as well as the number of commented 
  !> rows
  !> --------------------------------------------------------------------
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


  subroutine read_data(unit,rows,cols,data)
  !> --------------------------------------------------------------------
  !> subroutine read_data():
  !>    this subroutine reads a table of data (size rows,cols) from a file
  !> (given by "unit", already opened) and returns it on "data".
  !> Blank lines are ignored.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer     , intent(in)  :: unit,rows,cols
    real(dp), intent(out)     :: data(rows,cols)
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


  subroutine read_data_with_comments(unit,rows,cols,commented_rows,data,comments,mask)
  !> --------------------------------------------------------------------
  !> subroutine read_data_with_comments():
  !>    this subroutine reads a table of data (size rows,cols) from a file
  !> (given by "unit", already opened) and returns it on "data".
  !> The comments are returned in the variable "comments".
  !> "Mask" is an array of size total number of rows with .true. where data
  !> is, and .false. for the comments
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer     , intent(in)  :: unit,rows,cols,commented_rows
    real(dp), intent(out) :: data(rows,cols)
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


  subroutine sort_file(unit)
  !> --------------------------------------------------------------------
  !> subroutine sort_file():
  !>    This subroutine sorts the lines containing data of file 'unit'
  !> --------------------------------------------------------------------
    use mod_kind,     only: dp
    use mod_mpi_pars, only: abortProgram
    ! use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: unit
    character(len=50)   :: colformat
    integer             :: i,j,k,l,rows,cols,commented_rows
    character(len=900), allocatable :: comments(:)
    integer     ,       allocatable :: order(:)
    real(dp),           allocatable :: data(:,:),x(:)
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


  subroutine sort_command(filename,ncols,cols)
  !> --------------------------------------------------------------------
  !> subroutine sort_command():
  !>    This subroutine uses bash "sort" command to sort a file "filename"
  !> using "cols" columns
  !> --------------------------------------------------------------------
    implicit none
    character(len=400),        intent(in) :: filename
    integer,                   intent(in) :: ncols
    integer, dimension(ncols), intent(in) :: cols
    integer :: i
    character(len=500) :: command

    command = "sort"
    do i=1,ncols
      command = trim(command) // " -k" // trim(itos(cols(i))) // "g"
    end do
    command = trim(command) // " -o " // trim(filename) // " " // trim(filename)

    call execute_command_line(trim(command))

  end subroutine sort_command


  function next_line(procedure,f_unit,item)
  !> --------------------------------------------------------------------
  !> function next_line():
  !>    This function reads the next non-comment line in file f_unit
  !> --------------------------------------------------------------------
    use mod_mpi_pars, only: abortProgram
    implicit none
    character(len=*), optional, intent(in) :: procedure
    integer,                    intent(in) :: f_unit
    character(len=*), optional, intent(in) :: item
    integer             :: ios
    character(len=200)  :: next_line

    do
      read(f_unit, fmt='(A)', iostat = ios) next_line
      next_line = adjustl(next_line)
      if(.not.(len_trim(next_line) == 0 .or. (next_line(1:1) == "!"))) exit
    end do

    if(ios/=0) then
      if(present(item)) then
        call abortProgram("[" // trim(procedure) // "] Something went wrong when reading " // trim(item) // " in file unit " // trim(itos(f_unit)) )
      else
        call abortProgram("[" // trim(procedure) // "] Something went wrong when reading the file " // trim(itos(f_unit)) )
      end if
    end if
  end function next_line


  pure recursive function replaceStr(string,search,substitute) result(modifiedString)
  !> This function replaces a a given part "search" of a string "string"
  !> by another string "substitute".
  !> Function found in https://stackoverflow.com/questions/58938347/how-do-i-replace-a-character-in-the-string-with-another-charater-in-fortran
    implicit none
    character(len=*), intent(in)  :: string, search, substitute
    character(len=:), allocatable :: modifiedString
    integer                       :: i, stringLen, searchLen
    stringLen = len(string)
    searchLen = len(search)
    if (stringLen==0 .or. searchLen==0) then
      modifiedString = ""
      return
    elseif (stringLen<searchLen) then
      modifiedString = string
      return
    end if
    i = 1
    do
      if (string(i:i+searchLen-1)==search) then
        modifiedString = string(1:i-1) // substitute // replaceStr(string(i+searchLen:stringLen),search,substitute)
        exit
      end if
      if (i+searchLen>stringLen) then
        modifiedString = string
        exit
      end if
      i = i + 1
      cycle
    end do
  end function replaceStr


  character(len=20) function I4toS(i)
  !> --------------------------------------------------------------------
  !> function I4toS():
  !>    This function transforms a 4-bit integer i 
  !> into a character variable ItoS
  !> --------------------------------------------------------------------
    implicit none
    integer, intent(in) :: i
    write(I4toS, "(i0)") i
  end function I4toS

  character(len=20) function I8toS(i)
  !> --------------------------------------------------------------------
  !> function I4toS():
  !>    This function transforms an 8-bit integer i 
  !> into a character variable ItoS
  !> --------------------------------------------------------------------
    use mod_kind, only: int64
    implicit none
    integer(int64), intent(in) :: i
    write(I8toS, "(i0)") i
  end function I8toS


  character(len=900) function RtoS(r,format)
  !> --------------------------------------------------------------------
  !> function RtoS():
  !>    This function transforms a real r into a character variable RtoS.
  !> It also cuts leading spaces (on the left)
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    real(dp),         intent(in) :: r
    character(len=*), intent(in) :: format
    write(Rtos, fmt=format) r
    RtoS = adjustl(RtoS)
  end function RtoS


  function StoI_scalar(string)
  !> --------------------------------------------------------------------
  !> function StoI():
  !>    This function transforms a character variable into an integer StoI
  !> --------------------------------------------------------------------
    use mod_mpi_pars, only: abortProgram
    implicit none
    character(len=*) , intent(in):: string
    integer  :: ios
    integer  :: StoI_scalar

    read(unit= string, fmt=*, iostat=ios) StoI_scalar
    if(ios/=0) call abortProgram("[StoI_scalar] Something went wrong when reading the string: " // trim(string))
  end function StoI_scalar


  function StoI_array(string,dim_v)
  !> --------------------------------------------------------------------
  !> function StoI_array():
  !>    This function transforms a character variable into an integer array StoI
  !> --------------------------------------------------------------------
    use mod_mpi_pars, only: abortProgram
    implicit none
    character(len=*), intent(in) :: string
    integer,          intent(in) :: dim_v
    integer  :: ios
    integer  :: StoI_array(dim_v)

    read(unit= string, fmt=*, iostat=ios) StoI_array(1:dim_v)
    if(ios/=0) call abortProgram("[StoI_array] Something went wrong when reading the string: " // trim(string))
  end function StoI_array


  function StoR_scalar(string)
  !> --------------------------------------------------------------------
  !> function StoR_scalar():
  !>    This function transforms a character variable into a real StoR
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    use mod_mpi_pars, only: abortProgram
    implicit none
    character(len=*), intent(in) :: string
    integer  :: ios
    real(dp) :: StoR_scalar

    read(unit= string, fmt=*, iostat=ios) StoR_scalar
    if(ios/=0) call abortProgram("[StoR_scalar] Something went wrong when reading the string: " // trim(string))
  end function StoR_scalar


  function StoR_array(string,dim_v)
  !> --------------------------------------------------------------------
  !> function StoR_array():
  !>    This function transforms a character variable into a real array StoR
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    use mod_mpi_pars, only: abortProgram
    implicit none
    character(len=*), intent(in) :: string
    integer,          intent(in) :: dim_v
    integer  :: ios
    real(dp) :: StoR_array(dim_v)

    read(unit= string, fmt=*, iostat=ios) StoR_array(1:dim_v)
    if(ios/=0) call abortProgram("[StoR_array] Something went wrong when reading the string: " // trim(string))
  end function StoR_array


  function StoArray(string,dim_v)
  !> --------------------------------------------------------------------
  !> function StoArray():
  !>    This function transforms a character variable into a string array StoArray
  !> --------------------------------------------------------------------
    use mod_mpi_pars, only: abortProgram
    implicit none
    character(len=*), intent(in) :: string
    integer,          intent(in) :: dim_v
    integer           :: ios,i
    character(len=20) :: StoArray(dim_v)

    StoArray = ""
    read(unit=string, fmt=*, iostat=ios) (StoArray(i), i=1,dim_v)
  end function StoArray


  function KronProd(nax,nay,nbx,nby,A,B) result(AB)
  !> --------------------------------------------------------------------
  !> subroutine KronProd():
  !>    This subroutine calculates the outer product (Kronecker product)
  !> of two matrices A(nax,nay) and B(nbx,nby)
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,     intent(in)  :: nax, nay, nbx, nby
    complex(dp), intent(in)  :: A(nax,nay), B(nbx,nby) 
    complex(dp) :: AB(nax*nbx,nay*nby) 
    integer     :: i, j, p, q, l, m                                              
    do i = 1,nax
      do j = 1,nay
        l=(i-1)*nbx + 1
        m=l+nbx-1
        p=(j-1)*nby + 1
        q=p+nby-1
        AB(l:m,p:q) = A(i,j)*B(:,:)
      end do
    end do
  end function KronProd

  function vec_norm_complex(v, dim_v)
  !> --------------------------------------------------------------------
  !> function vec_norm_complex():
  !>    This function calculates the norm of a complex vector v of dimension
  !> dim_v.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,                       intent(in) :: dim_v ! vector dimension
    complex(dp), dimension(dim_v), intent(in) :: v ! vector v
    real(dp)                                  :: vec_norm_complex
    real(dp)                                  :: rsum
    integer                                   :: i

    rsum= 0._dp
    !$omp simd reduction(+:rsum)
    do i = 1, dim_v
      rsum = rsum + real(v(i)*conjg(v(i)))
    end do
    vec_norm_complex= sqrt(rsum)
  end function vec_norm_complex


  function vec_norm_real(v, dim_v)
  !> --------------------------------------------------------------------
  !> function vec_norm_real():
  !>    This function calculates the norm of a real vector v of dimension
  !> dim_v.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,                        intent(in) :: dim_v ! vector dimension
    real(dp), dimension(dim_v), intent(in) :: v ! vector v
    real(dp)                               :: vec_norm_real

    vec_norm_real = sqrt(dot_product(v,v))
  end function vec_norm_real


  function normalize_complex(v, dim_v)
  !> --------------------------------------------------------------------
  !> function normalize_complex():
  !>    This function normalizes the complex vector v of dimension
  !> dim_v.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,                           intent(in) :: dim_v ! vector dimension
    complex(dp), dimension(dim_v), intent(in) :: v ! vector v
    complex(dp), dimension(dim_v)             :: normalize_complex

    normalize_complex = v/vec_norm(v,dim_v)
  end function normalize_complex


  function normalize_real(v, dim_v)
  !> --------------------------------------------------------------------
  !> function normalize_real():
  !>    This function normalizes the real vector v of dimension
  !> dim_v.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,                    intent(in) :: dim_v ! vector dimension
    real(dp), dimension(dim_v), intent(in) :: v ! vector v
    real(dp), dimension(dim_v)             :: normalize_real

    normalize_real = v/vec_norm(v,dim_v)
  end function normalize_real


  subroutine transvComponent(vec,dir0)
  !> --------------------------------------------------------------------
  !> subroutine transvComponent(vec,dir0):
  !>    This subroutine calculates the transverse components of vec
  !> with respect to the direction dir0 (and returns it in vec)
  !> --------------------------------------------------------------------
    use mod_kind,     only: dp
    implicit none
    real(dp), dimension(3), intent(inout) :: vec
    real(dp), dimension(3), intent(in)    :: dir0
    real(dp), dimension(3)                :: dir, vecp, vect

    ! Normalize the direction
    dir = dir0/vec_norm(dir0,3)
    ! Get the component of vec in the direction dir (longitudinal)
    vecp(:) = dot_product(vec,dir)*dir(:)
    ! Subtract the longitudinal component to obtain the transverse ones
    vect(:) = vec(:)-vecp(:)
    ! Return in vec
    vec(:) = vect(:)
  end subroutine transvComponent


  subroutine LS_solver(n,a,b)
  !> --------------------------------------------------------------------
  !> subroutine LS_solver(n,a,b):
  !>    This subroutine is a simple interface for the LAPACK linear
  !> system solver A*X = B.
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,     intent(in)      :: n
    complex(dp), intent(inout)   :: a(n,n)
    complex(dp), intent(inout)   :: b(n)
    ! Workspace variables
    integer           :: info
    integer           :: ipiv(n)

    external :: zgesv

    call zgesv( n, 1, a, n, ipiv, b, n, info )
        
  end subroutine LS_solver


  subroutine diagonalize(n,a,eval)
  !> --------------------------------------------------------------------
  !> subroutine diagonalize():
  !>    This subruotine is a simple interface for the LAPACK zheev subroutine
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,     intent(in)    :: n
    complex(dp), intent(inout) :: a(n,n)
    real(dp),    intent(out)   :: eval(n)
    ! Workspace variables
    real(dp)    :: rwork(3*n-2)
    integer     :: info
    complex(dp) :: work(2*n-1) ! dimension of lwork below

    external :: zheev    

    call zheev('V','U',n,a,n,eval,work,lwork,rwork,info)
        
  end subroutine diagonalize


  subroutine invers(matriz,nn)
  !> --------------------------------------------------------------------
  !> function invers(matriz,nn):
  !>    This subroutine calculates the inverse of nn x nn matrix 'matriz'
  !> --------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env
    use mod_mpi_pars
    use mod_parameters, only: output
    implicit none
    integer :: nn,info
    integer,     dimension(nn)    :: ipiv
    complex(dp), dimension(nn,nn) :: matriz
    complex(dp), dimension(nn*4)  :: work

    external :: zgetrf,zgetri
    lwork = 4*nn
    info = 0
    call zgetrf(nn,nn,matriz,nn,ipiv,info)
    if (info>0) then
      write(output%unit,"('[invers] Singular matrix! info = ',i0,' on rank ',i0)") info,myrank
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    if (info<0) then
      write(output%unit,"('[invers] Illegal value of argument ',i0,' on rank ',i0)") -info,myrank
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    call zgetri(nn,matriz,nn,ipiv,work,lwork,info)

    if (info/=0) then
      write(output%unit,"('[invers] Singular matrix! info = ',i0,' on rank ',i0)") info,myrank
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if

  end subroutine invers


  subroutine print_matrix_real(unit,n,mat)
  !> --------------------------------------------------------------------
  !> function print_matrix_real():
  !>      This subroutine prints the real elements of matrix 'mat'
  !> --------------------------------------------------------------------
    use mod_kind, only: dp
    implicit none
    integer,                  intent(in) :: unit,n
    real(dp), dimension(n,n), intent(in) :: mat
    integer :: i,j

    do j = 1, n
      do i = 1, n
        write(unit=unit,fmt=*) mat(i,j)
      end do
    end do

  end subroutine print_matrix_real


  subroutine print_matrix_complex(unit,n,mat)
  !> --------------------------------------------------------------------
  !> function print_matrix_real():
  !>      This subroutine prints the real and imaginary parts of 
  !>  the elements of matrix 'mat'
  !> --------------------------------------------------------------------
    use mod_kind,   only: dp
    implicit none
    integer,                     intent(in) :: unit,n
    complex(dp), dimension(n,n), intent(in) :: mat
    integer :: i,j

    do j = 1, n
      do i = 1, n
        write(unit=unit,fmt=*) mat(i,j)%re, mat(i,j)%im
      end do
    end do

  end subroutine print_matrix_complex


  function get_memory(units,mem) result(success)
  !> --------------------------------------------------------------------
  !> function get_memory():
  !>    This function returns the free memory in the cpu or gpu where
  !> the job is running
  !> --------------------------------------------------------------------
    use mod_kind,       only: int32
    use mod_parameters, only: output
#ifdef _GPU
    use mod_mpi_pars,   only: myrank,rField
#else
    use mod_mpi_pars,   only: rField
#endif
    ! use mod_io, only: log_warning
    implicit none
    character(len=1), intent(in)  :: units
    !> In which units the result will be returned: g = GB, m = MB, k = KB
    integer(int32),   intent(out) :: mem
    !> Amount of free memory
    logical :: success

#ifdef _GPU
    character(len=20)  :: filename = "./gpu_memory", keyword="MiB"
    integer            :: isize =3 , memunit=1
#else
    character(len=13)  :: filename = "/proc/meminfo", keyword="MemFree:"
    integer            :: isize =0 , memunit=1024
#endif
    integer, parameter :: max_elements=50
    integer            :: i,ios,pos,file_unit = 931
    character(200)     :: line
    character(len=80), dimension(max_elements) :: str_arr
    character(len=80)  :: str_tmp

#ifdef _GPU
    filename = trim(filename) // itos(myrank)
    call execute_command_line( "nvidia-smi | sed 's|/||g' > " // trim(filename) )
#endif

    open (unit=file_unit, file=trim(filename), status='old', action='read', iostat=ios)
    if(ios/=0) then
      if(rField == 0) write(unit=output%unit, fmt="('[Warning] [get_memory] Could not get memory from file ',a,'.')") trim(filename)
      success = .false.
      return
    end if

    lines: do
      str_arr = ""
      read(unit=file_unit, fmt='(A)', iostat=ios) line
      if(ios/=0) exit
      read(unit=line, fmt=*, iostat=ios) (str_arr(i), i = 1, max_elements)
      do i = 1, max_elements-1
        if(len_trim(str_arr(i)) == 0) cycle
        pos = index(str_arr(i),trim(keyword))
        if(pos/=0) then
          str_tmp = str_arr(i+1)(1:len_trim(str_arr(i+1))-isize)
          select case(units)
          case("m")
            mem = int(dble(StoI(str_tmp))/memunit)
          case("g")
            mem = int(dble(StoI(str_tmp))/memunit/1024)
          end select
#ifdef _GPU
          call execute_command_line('rm ' // trim(filename))
#endif
          close(file_unit)
          success = .true.
          return
        end if
      end do
    end do lines
    if(rField == 0) write(unit=output%unit, fmt="('[Warning] [get_memory] Could not find ',a,' on file ',a,'.')") trim(keyword),trim(filename)

#ifdef _GPU
    call execute_command_line('rm '// trim(filename) )
#endif
    success = .false.
    close(file_unit)
  end function get_memory

end module mod_tools


