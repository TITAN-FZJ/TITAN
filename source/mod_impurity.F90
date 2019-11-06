module mod_impurity
implicit none
character(len=200) :: impurity_file = ""
logical :: limpurity = .false.

contains

    subroutine read_impurity(filename, host, s)
      use mod_f90_kind,  only: double
      use mod_system,    only: System
      use mod_mpi_pars,  only: myrank,abortProgram
      use mod_constants, only: tpi
      use mod_tools,     only: cross, vec_norm
      use mod_io,        only: log_warning
      implicit none
  !
      character(len=*), intent(in)    :: filename
      type(System),     intent(inout) :: s
      type(System),     intent(in)    :: host
  !
      integer, parameter :: line_length = 300, word_length = 50, max_elements = 50
      integer :: i, j, k, l
      integer :: f_unit = 99, ios
      character(len=word_length), dimension(max_elements) :: str_arr
      character(len=line_length) :: line
      character(len=1)           :: coord_type
      integer, dimension(50)     :: type_count
      real(double), dimension(3) :: zdir, ydir

      ! write(*,*) "Hi again"

      open(unit = f_unit, file=trim(filename), status='old', iostat = ios)
      if((ios /= 0).and.(myrank==0)) call abortProgram("[read_impurity] Error reading " // trim(filename))

      ! Counting number elements
      str_arr = ""
      read(f_unit, fmt='(A)', iostat=ios) line
      read(unit=line, fmt=*, iostat=ios) (str_arr(i), i = 1, max_elements)

      ! write(*,*) line
      do i = 1, max_elements
        if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
        s%nTypes = s%nTypes + 1
      end do

      ! Reading name of elements
      allocate(s%Types(s%nTypes))
      j = 1
      do i = 1, max_elements
        if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
        s%Types(j)%Name = str_arr(i)
        ! write(*,*) str_arr(i)
        j = j + 1
      end do

      ! Reading quantity of each element
      read(f_unit, fmt='(A)', iostat=ios) line
      read(unit=line, fmt=*, iostat=ios) (type_count(i), i = 1, s%nTypes)

      ! write(*,*) line
      ! write(*,*) type_count(:)

      ! Reading units of coordinates
      read(f_unit, fmt='(A)', iostat=ios) line
      read(unit=line, fmt=*, iostat=ios) coord_type

      ! write(*,*) coord_type

      ! Testing if total number of atoms is correct
      s%nAtoms = sum(type_count(1:s%nTypes))
      if(s%nAtoms <= 0) call abortProgram("[read_basis] No basis atoms given in '" // trim(filename) // "'!")
      allocate(s%Basis(s%nAtoms))

      ! write(*,*) type_count(1)

      ! Reading list of atoms in the unit cell and multiplying by the correct units
      k = 0
      do i = 1, s%nTypes
        do j = 1, type_count(i)
          line = ""
          do while(len_trim(line) == 0 .or. len_trim(line) == 200)
            read(f_unit, fmt='(A)', iostat=ios) line
            if((ios /= 0).and.(myrank==0)) &
              call abortProgram("[read_basis] Not enough basis atoms given in '" // trim(filename) // "'!")
          end do
          k = k + 1
          s%Basis(k)%Material = i
          read(unit=line, fmt=*, iostat=ios) (s%Basis(k)%Position(l), l=1,3)
          if(coord_type == 'C' .or. coord_type == 'c' .or. coord_type == 'K' .or. coord_type == 'k') then
            ! Positions of atoms given in Cartesian coordinates
            s%Basis(k)%Position = s%Basis(k)%Position * host%a0
          else
            ! Position of atoms given in Bravais (or Direct, Internal, Lattice) coordinates
            if((vec_norm(host%a2,3) <= 1.d-9).and.(myrank==0)) &
              call log_warning("read_basis", "a2 not given in '" // trim(filename) // "', while Bravais used for atom positions!")
            if((vec_norm(host%a3,3) <= 1.d-9).and.(myrank==0)) &
              call log_warning("read_basis", "a3 not given in '" // trim(filename) // "', while Bravais used for atom positions!")
            s%Basis(k)%Position = s%Basis(k)%Position(1) * host%a1 + s%Basis(k)%Position(2) * host%a2 + s%Basis(k)%Position(3) * host%a3
            ! write(*,*) s%Basis(k)%Position
          end if
        end do
      end do

      close(f_unit)

  end subroutine read_impurity

end module mod_impurity
