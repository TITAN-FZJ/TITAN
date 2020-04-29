module mod_input
  implicit none
  integer :: nlines
  integer, parameter :: word_length = 50
  integer, parameter :: line_length = 300
  integer, parameter :: max_lines   = 100
  integer, parameter :: max_elements = 300
  character(line_length) :: out_file = ""
  integer :: out_unit = -1
  character(len=line_length), dimension(max_lines) :: key, val

  interface get_parameter
     module procedure get_int, &
                      get_int_default, &
                      get_int4_array, &
                      get_int8_array, &
                      get_real, &
                      get_real_default, &
                      get_real_array, &
                      get_string, &
                      get_string_default, &
                      get_string_array, &
                      get_complex, &
                      get_complex_default, &
                      get_complex_array, &
                      get_logical, &
                      get_logical_default
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
  end function disable_input_logging

  subroutine log_parameter(key_string, val_string)
    use mod_mpi_pars
    implicit none
    character(len=*), intent(in) :: key_string, val_string

    if(myrank /= 0) return
    if(out_unit < 0) return

    write(out_unit, "(' ',a,' = ',a,' ')") trim(adjustl(key_string)), trim(adjustl(val_string))
  end subroutine log_parameter

  function read_file(filename) result(success)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=line_length) :: line
    logical   :: success
    integer*4 :: eof=0, iunit=92412
    integer*4 :: eq_pos, com_pos

    nlines = 0

    success = .false.

    val = " "
    open(unit=iunit, file=trim(filename), status='old', iostat=eof)
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
    close(unit=iunit)

    success = .true.
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

  function get_int_default(key_val, ret_val, default) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    integer, intent(in) :: default
    integer, intent(out) :: ret_val
    logical :: success
    character(len=50) :: default_string

    success = get_int(key_val, ret_val)
    if(.not. success) then
      ret_val = default
      write(default_string, *) default
      call log_parameter(key_val, default_string)
    end if
  end function get_int_default

  function get_int4_array(key_val, ret_val, ret_cnt) result(success)
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
  end function get_int4_array

  function get_int8_array(key_val, ret_val, ret_cnt) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    integer*8, allocatable, intent(out) :: ret_val(:)
    integer, intent(out) :: ret_cnt
    logical :: success
    integer :: ind , ios, i
    integer*8 :: tmp_arr(max_elements)
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
  end function get_int8_array

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

  function get_real_default(key_val, ret_val, default) result(success)
    use mod_f90_kind, only:double
    implicit none
    character(len=*), intent(in) :: key_val
    real(double), intent(out) :: ret_val
    real(double), intent(in) :: default
    logical :: success
    character(len=50) :: default_string

    success = get_real(key_val, ret_val)
    if(.not. success) then
      ret_val = default
      write(default_string, *) default
      call log_parameter(key_val, default_string)
    end if
  end function get_real_default

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
    use mod_f90_kind, only: double
    implicit none
    character(len=*), intent(in) :: key_val
    complex(double), intent(out) :: ret_val
    logical :: success
    integer :: ind , ios

    success = find_val(key_val, ind)
    if(success) then
       read(unit=val(ind), fmt=*, iostat=ios) ret_val
       if(ios /= 0) success = .false.
    end if
  end function get_complex

  function get_complex_default(key_val, ret_val, default) result(success)
    use mod_f90_kind, only: double
    implicit none
    character(len=*), intent(in) :: key_val
    complex(double), intent(out) :: ret_val
    complex(double), intent(in) :: default
    logical :: success
    character(len=50) :: default_string

    success = get_complex(key_val, ret_val)
    if(.not. success) then
      ret_val = default
      write(default_string, *) default
      call log_parameter(key_val, default_string)
    end if
  end function get_complex_default

  function get_complex_array(key_val, ret_val, ret_cnt) result(success)
    use mod_f90_kind, only: double
    implicit none
    character(len=*), intent(in) :: key_val
    complex(double), allocatable, intent(out) :: ret_val(:)
    integer, intent(out) :: ret_cnt
    logical :: success
    integer :: ind , ios, i
    complex(double) :: tmp_arr(max_elements)
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

  function get_string_default(key_val, ret_val, default) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    character(len=*), intent(out) :: ret_val
    character(len=*), intent(in) :: default
    logical :: success
    character(len=50) :: default_string

    success = get_string(key_val, ret_val)
    if(.not. success) then
      ret_val = default
      write(default_string, *) default
      call log_parameter(key_val, default_string)
    end if
  end function get_string_default

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

  function get_logical_default(key_val, ret_val, default) result(success)
    implicit none
    character(len=*), intent(in) :: key_val
    logical, intent(out) :: ret_val
    logical, intent(in) :: default
    logical :: success
    character(len=50) :: default_string

    success = get_logical(key_val, ret_val)
    if(.not. success) then
      ret_val = default
      write(default_string, *) default
      call log_parameter(key_val, default_string)
    end if
  end function get_logical_default
end module mod_input
