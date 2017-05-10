module mod_input
  implicit none
  integer :: nlines
  integer, parameter :: word_length = 50
  integer, parameter :: line_length = 300
  integer, parameter :: max_lines   = 100
  integer, parameter :: max_elements = 20
  character(line_length) :: out_file = ""
  integer :: out_unit = -1
  character(len=line_length), dimension(max_lines) :: key, val

  interface get_parameter
     module procedure get_int, get_int_array, &
          get_real, get_real_array, &
          get_string, get_string_array, &
          get_complex, get_complex_array, &
          get_logical
  end interface get_parameter
contains

  function enable_input_logging(filename) result(success)
    implicit none
    character(len=*), intent(in) :: filename
    logical :: success
    integer :: eof
    success = .false.
    open(unit = 123456788, file=trim(filename), status='replace', iostat=eof)
    if(eof /= 0) return
    out_unit = 123456788
    success = .true.

    return
  end function enable_input_logging

  function disable_input_logging() result(success)
    implicit none
    logical :: success
    integer :: eof
    success = .true.
    if(out_unit < 0) return
    success = .false.
    close(unit=out_unit, iostat=eof)
    if(eof /= 0) return
    success = .true.
    return
  end function disable_input_logging

  subroutine log_parameter(key_string, val_string)
    use mod_mpi_pars
    implicit none
    character(len=*), intent(in) :: key_string, val_string

    if(myrank /= 0) return
    if(out_unit < 0) return

    write(out_unit, "(' ',a,' = ',a,' ')") trim(adjustl(key_string)), trim(adjustl(val_string))

    return
  end subroutine log_parameter

  function read_file(filename) result(success)
    implicit none
    character(len=200), intent(in) :: filename
    character(len=line_length) :: line
    integer :: eof, iunit
    logical :: success
    integer :: eq_pos, com_pos

    eof = 0
    iunit = 92412
    success = .false.

    nlines = 0

    val = " "
    open(unit=iunit, file=trim(filename), iostat=eof)
    if(eof /= 0) return
    do while(eof == 0)
       read(unit=iunit, fmt='(A)', iostat=eof) line
       if(line(1:2) /= "->") cycle
       nlines = nlines + 1
       eq_pos = index(line, "=")
       com_pos = index(line, "!")
       key(nlines) = adjustl(line(3:eq_pos-1))
       if(com_pos > 0) then
          val(nlines) = adjustl(line(eq_pos+1:com_pos-1))
       else
          val(nlines) = adjustl(line(eq_pos+1:))
       end if

    end do
    close(unit=iunit, iostat = eof)
    if(eof /= 0) return

    success = .true.
    return

  end function read_file

  function find_val(key_val, ind) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    integer, intent(out) :: ind
    logical success
    integer :: i
    success = .false.

    do i = 1,nlines
       if(trim(key_val) == trim(key(i))) then
          ind = i
          success = .true.
          call log_parameter(key(i), val(i))
          exit
       else
       end if
    end do

  end function find_val

  function get_int(key_val, ret_val) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    integer, intent(out) :: ret_val
    logical :: success
    integer :: ind , ios

    success = find_val(key_val, ind)
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) ret_val
       if(ios /= 0) success = .false.
    end if
  end function get_int

  function get_int_array(key_val, ret_val, ret_cnt) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    integer, allocatable, intent(out) :: ret_val(:)
    integer, intent(out) :: ret_cnt
    logical :: success
    integer :: ind , ios, i
    integer :: tmp_arr(max_elements)
    character(len=word_length) :: str_arr(max_elements)

    success = find_val(key_val, ind)
    ret_cnt = 0
    do i=1,max_elements
       str_arr(i) = ""
    end do
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) (str_arr(i), i=1,max_elements)
       do i = 1, max_elements
          if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
          ret_cnt = ret_cnt + 1
          read(unit=str_arr(i), fmt=*,iostat=ios ) tmp_arr(ret_cnt)
       end do
       allocate(ret_val(ret_cnt))
       ret_val = tmp_arr(:ret_cnt)
    end if
  end function get_int_array

  function get_real(key_val, ret_val) result(success)
    use mod_f90_kind, only:double
    implicit none
    character(len=*), intent(in) :: key_val
    real(double), intent(out) :: ret_val
    logical :: success
    integer :: ind , ios

    success = find_val(key_val, ind)
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) ret_val
       if(ios /= 0) success = .false.
    end if
  end function get_real

  function get_real_array(key_val, ret_val, ret_cnt) result(success)
    use mod_f90_kind, only: double
    implicit none
    character(len=*), intent(in) :: key_val
    real(double), allocatable, intent(out) :: ret_val(:)
    integer, intent(out) :: ret_cnt
    logical :: success
    integer :: ind , ios, i
    real(double) :: tmp_arr(max_elements)
    character(len=word_length) :: str_arr(max_elements)

    success = find_val(key_val, ind)
    ret_cnt = 0
    do i=1,max_elements
       str_arr(i) = ""
    end do
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) (str_arr(i), i=1,max_elements)
       do i = 1, max_elements
          if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
          ret_cnt = ret_cnt + 1
          read(unit=str_arr(i), fmt=*,iostat=ios ) tmp_arr(ret_cnt)
       end do
       allocate(ret_val(ret_cnt))
       ret_val = tmp_arr(:ret_cnt)
    end if
  end function get_real_array

  function get_complex(key_val, ret_val) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    double complex, intent(out) :: ret_val
    logical :: success
    integer :: ind , ios

    success = find_val(key_val, ind)
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) ret_val
       if(ios /= 0) success = .false.
    end if
  end function get_complex

  function get_complex_array(key_val, ret_val, ret_cnt) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    double complex, allocatable, intent(out) :: ret_val(:)
    integer, intent(out) :: ret_cnt
    logical :: success
    integer :: ind , ios, i
    double complex :: tmp_arr(max_elements)
    character(len=word_length) :: str_arr(max_elements)

    success = find_val(key_val, ind)
    ret_cnt = 0
    do i=1,max_elements
       str_arr(i) = ""
    end do
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) (str_arr(i), i=1,max_elements)
       do i = 1, max_elements
          if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
          ret_cnt = ret_cnt + 1
          read(unit=str_arr(i), fmt=*,iostat=ios ) tmp_arr(ret_cnt)
       end do
       allocate(ret_val(ret_cnt))
       ret_val = tmp_arr(:ret_cnt)
    end if
  end function get_complex_array

  function get_string(key_val, ret_val) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    character(len=*), intent(out) :: ret_val
    logical :: success
    integer :: ind , ios

    success = find_val(key_val, ind)
    if(success) then
       read(unit=val(ind), fmt='(A)', iostat=ios) ret_val
       if(ios /= 0) success = .false.
    end if
  end function get_string

  function get_string_array(key_val, ret_val, ret_cnt) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    character(len=*), allocatable, intent(out) :: ret_val(:)
    integer, intent(out) :: ret_cnt
    logical :: success
    integer :: ind , ios, i
    character(len=word_length) :: tmp_arr(max_elements)
    character(len=word_length) :: str_arr(max_elements)
    do i=1,max_elements
       str_arr(i) = ""
       tmp_arr(i) = ""
    end do
    success = find_val(key_val, ind)
    ret_cnt = 0

    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) (str_arr(i), i=1,max_elements)
       do i = 1, max_elements
          if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
          ret_cnt = ret_cnt + 1
          read(unit=str_arr(i), fmt=*,iostat=ios ) tmp_arr(ret_cnt)
       end do
       allocate(ret_val(ret_cnt))
       ret_val(1:ret_cnt) = tmp_arr(1:ret_cnt)
    end if
  end function get_string_array

  function get_logical(key_val, ret_val) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    logical, intent(out) :: ret_val
    logical :: success
    integer :: ind , ios

    success = find_val(key_val, ind)
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) ret_val
       if(ios /= 0) success = .false.
    end if
  end function get_logical

end module mod_input
