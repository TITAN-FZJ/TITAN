module mod_io
  use mod_f90_kind
  use mod_parameters
  implicit none
  integer           :: ifile=666666,nlines
#ifdef _JUQUEEN
  character(len=20) :: filename="inputdhe_juqueen"
#else
  character(len=20) :: filename="inputdhe"
#endif
  character(len=200),allocatable :: string1(:),stringt(:)

contains
  subroutine read_input_file()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_tools, only: number_of_lines
    implicit none
    integer                 :: ios,i,j,nlinest

    open(unit=ifile, file=trim(filename), status='old', iostat=ios)
    if(ios.ne.0) then
      if(myrank.eq.0) write(outputunit,"('[read_input_file] File ""',a,'"" not found!')") trim(filename)
      call MPI_Finalize(ierr)
      stop
    end if

    call number_of_lines(ifile,nlinest,nlines)
!     write(outputunit,"('[get_parameters] ""',a,'"" file has ',i0,' non-commented lines and ',i0,' in total (non-blank only)')") trim(filename),nlines,nlinest
    allocate(string1(nlines),stringt(nlinest))

    ! Writing all lines in array stringt(:) and non-commented lines in string1(:)
    rewind ifile
    j = 0
    i = 1
    do
      read(unit=ifile,fmt='(A)',iostat=ios) stringt(i)
      if (ios.ne.0) exit
      if ((stringt(i)(1:1).eq."#").or.(stringt(i)(1:1).eq."!").or.(stringt(i).eq."")) cycle  ! If the line is a comment or blank, ignore
      j = j+1
      string1(j) = stringt(i)
      i = i+1
      if(i.gt.nlinest) exit
!       write(outputunit,*) string1(j)
    end do
    close(ifile)

    return
  end subroutine read_input_file


  subroutine read_output_filename()
    use mod_f90_kind
    use mod_mpi_pars
    implicit none
    integer, parameter      :: iomax=20           ! Maximum number of elements in one line
    integer                 :: ios,i,j,filenameini
    character(len=20)       :: istring1(iomax),istring2(iomax)

    ! Getting variables from non-commented lines
    istring1 = " "
    istring2 = " "
    read(unit=string1(1),fmt=*,iostat=ios) (istring1(i),i=1,iomax) ! Reads first line
    lines_loop: do j=1,nlines
      ! Reads next line for variables defined in two lines
      if(j.ne.nlines) read(unit=string1(j+1),fmt=*,iostat=ios) (istring2(i),i=1,iomax)
      read_line: do i=1,iomax
        ioparams: select case (istring1(i))
!===============================================================================
        case("output","output:","output=")
          if ((index(string1(j),"output:").gt.0).or.(index(string1(j),"output=").gt.0)) then
            filenameini = index(string1(j),"output")+8
          else if (index(string1(j),"output =").gt.0) then
            filenameini = index(string1(j),"output")+9
          else
            filenameini = index(string1(j),"output")+7
          end if
          outputdhe = string1(j)(filenameini:)
          exit lines_loop
!===============================================================================
        case("!")
          exit
        end select ioparams
      end do read_line
      istring1 = istring2
      istring2 = " "
    end do lines_loop

  if((trim(outputdhe).eq."").and.(myrank.eq.0)) then
    write(*,"('[read_output_filename] Output filename not given!')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if

!   if(myrank.eq.0) then
!     if(trim(outputdhe).eq."") then
!       outputunit = 6
!       write(outputunit,"('[read_output_filename] Output file not given. Writing on screen...')")
!     end if
!   end if

    return
  end subroutine read_output_filename


  subroutine get_parameters()
    use mod_f90_kind
    use mod_mpi_pars
    implicit none
    integer, parameter      :: iomax=20           ! Maximum number of elements in one line
    integer                 :: nparams=24         ! Number of required parameters to be read
    integer                 :: nparams2=0         ! Number of secondary parameters required
    integer                 :: ios,i,j,n,filenameini
    character(len=20)       :: istring1(iomax),istring2(iomax)

    ! Getting variables from non-commented lines of input file
    istring1 = " "
    istring2 = " "
    scfile = " "
    runoptions = " "
    if(myrank.eq.0) write(outputunit,"('[get_parameters] Reading parameters from ""',a,'"" file...')") trim(filename)
    read(unit=string1(1),fmt=*,iostat=ios) (istring1(i),i=1,iomax) ! Reads first line
    lines_loop: do j=1,nlines
      ! Reads next line for variables defined in two lines
      if(j.ne.nlines) read(unit=string1(j+1),fmt=*,iostat=ios) (istring2(i),i=1,iomax)
      read_line: do i=1,iomax
!         write(*,*) i,j
        ioparams: select case (istring1(i))
!===============================================================================
        case("itype")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) itype
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) itype
          end if
!           write(outputunit,"('itype = ',i0)") itype
          nparams = nparams-1
        case("itype=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) itype
!           write(outputunit,"('itype = ',i0)") itype
          nparams = nparams-1
!===============================================================================
        case("Options","Options:")
          do n=i+1,iomax
            if(istring1(n).eq."") cycle
            options: select case (istring1(n))
            case ("ry2ev")
              ry2ev = 13.6d0
            case ("tesla")
              tesla = 5.7883817555d-5/13.6d0 ! Ry/T
              ltesla = .true.
            case ("verbose")
              lverbose = .true.
            case ("debug")
              ldebug = .true.
            case ("addresults")
              laddresults = .true.
            case ("createfiles")
              lcreatefiles = .true.
            case ("noUonall")
              if(Utype.eq.1) then
                if(myrank.eq.0) write(outputunit,"('[get_parameters] Runoption ""noUonNM"" is already active')")
              else
                Utype = 0
              end if
            case ("noUonNM")
              if(Utype.eq.0) then
                if(myrank.eq.0) write(outputunit,"('[get_parameters] Runoption ""noUonall"" is already active')")
              else
                Utype = 1
              end if
            case ("slatec")
              lslatec = .true.
            case ("GSL")
              lGSL = .true.
            case ("kpoints")
              lkpoints = .true.
            case ("lineargfsoc")
              llineargfsoc = .true.
            case ("linearsoc")
              llinearsoc = .true.
            case ("nojac")
              lnojac = .true.
            case ("hfresponses")
              lhfresponses  = .true.
            case ("ontheflysc")
              lontheflysc = .true.
            case ("rotatemag")
              lrotatemag = .true.
            case ("nolb")
              lnolb = .true.
            case ("nodiag")
              lnodiag = .true.
            case ("sha")
              lsha = .true.
            case ("writeonscreen")
              lwriteonscreen = .true.
            case ("sortfiles")
              lsortfiles = .true.
            case("!")
              exit
            case default
              if(myrank.eq.0) write(outputunit,"('[get_parameters] Runoption ""',a,'"" not found!')") trim(istring1(n))
              cycle
            end select options
            runoptions  = trim(runoptions) // " " // trim(istring1(n))
          end do
          exit
!===============================================================================
        case("eta")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) eta
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) eta
          end if
!           write(outputunit,"('eta = ',f8.5)") eta
          nparams = nparams-1
        case("eta=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) eta
!           write(outputunit,"('eta = ',f8.5)") eta
          nparams = nparams-1
!===============================================================================
        case("lattice")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            lattice = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            lattice = istring2(i)
          end if
!           write(outputunit,"('lattice = ',a2)") lattice
          nparams = nparams-1
        case("lattice=")
          lattice = istring1(i+1)
!           write(outputunit,"('lattice = ',a2)") lattice
          nparams = nparams-1
!===============================================================================
        case("Npl")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) Npl_i
            if(index(istring1(i+3),"!").eq.0) read(unit=istring1(i+3),fmt=*,iostat=ios) Npl_f
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) Npl_i
            if(index(istring2(i+1),"!").eq.0) read(unit=istring2(i+1),fmt=*,iostat=ios) Npl_f
          end if
!           write(outputunit,"('Npl_i = ',i0)") Npl_i
!           write(outputunit,"('Npl_f = ',i0)") Npl_f
          nparams = nparams-1
        case("Npl=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) Npl_i
          read(unit=istring1(i+2),fmt=*,iostat=ios) Npl_f
!           write(outputunit,"('Npl_i = ',i0)") Npl_i
!           write(outputunit,"('Npl_f = ',i0)") Npl_f
          nparams = nparams-1
!===============================================================================
        case("a0")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) a0
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) a0
          end if
!           write(outputunit,"('a0 = ',f8.5)") a0
          nparams = nparams-1
        case("a0=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) a0
!           write(outputunit,"('a0 = ',f8.5)") a0
          nparams = nparams-1
!===============================================================================
        case("SOC")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) SOC
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) SOC
          end if
          if(SOC) then
!             write(outputunit,"('Spin Orbit Coupling: ACTIVATED')")
            nparams2 = nparams2+2 ! magaxis and socscale
          else
!             write(outputunit,"('Spin Orbit Coupling: DEACTIVATED')")
          end if
          nparams = nparams-1
        case("SOC=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) SOC
          if(SOC) then
!             write(outputunit,"('Spin Orbit Coupling: ACTIVATED')")
            nparams2 = nparams2+2 ! magaxis and socscale
          else
!             write(outputunit,"('Spin Orbit Coupling: DEACTIVATED')")
          end if
          nparams = nparams-1
!-------------------------------------------------------------------------------
        case("socscale")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) socscale
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) socscale
          end if
!           write(outputunit,"('socscale = ',f8.5)") socscale
          nparams2 = nparams2-1
        case("socscale=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) socscale
!           write(outputunit,"('socscale = ',f8.5)") socscale
          nparams2 = nparams2-1
!-------------------------------------------------------------------------------
        case("magaxis")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            magaxis = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            magaxis = istring2(i)
          end if
!           write(outputunit,"('magaxis = ',a2)") magaxis
          nparams2 = nparams2-1
        case("magaxis=")
          magaxis = istring1(i+1)
!           write(outputunit,"('magaxis = ',a2)") magaxis
          nparams2 = nparams2-1
!===============================================================================
        case("FIELD")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) lfield
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) lfield
          end if
          if(lfield) then
!             write(outputunit,"('Static Magnetic field: ACTIVATED')")
            nparams2 = nparams2+3 ! magaxis and socscale
          else
!             write(outputunit,"('Static Magnetic field: DEACTIVATED')")
          end if
          nparams = nparams-1
        case("FIELD=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) lfield
          if(lfield) then
!             write(outputunit,"('Static Magnetic field: ACTIVATED')")
            nparams2 = nparams2+3 ! magaxis and socscale
          else
!             write(outputunit,"('Static Magnetic field: DEACTIVATED')")
          end if
          nparams = nparams-1
!-------------------------------------------------------------------------------
        case("hwa")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwa_i
            if(index(istring1(i+3),"!").eq.0) then
              read(unit=istring1(i+3),fmt=*,iostat=ios) hwa_f
              read(unit=istring1(i+4),fmt=*,iostat=ios) hwa_npts
            end if
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwa_i
            if(index(istring1(i+1),"!").eq.0) then
              read(unit=istring2(i+1),fmt=*,iostat=ios) hwa_f
              read(unit=istring2(i+2),fmt=*,iostat=ios) hwa_npts
            end if
          end if
!           write(outputunit,"('hwa_i = ',f8.5)") hwa_i
!           write(outputunit,"('hwa_f = ',f8.5)") hwa_f
!           write(outputunit,"('hwa_npts = ',f8.5)") hwa_npts
          hwa_npt1 = hwa_npts + 1
          nparams2 = nparams2-1
        case("hwa=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwa_i
          if(index(istring1(i+2),"!").eq.0) then
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwa_f
            read(unit=istring1(i+3),fmt=*,iostat=ios) hwa_npts
          end if
!           write(outputunit,"('hwa_i = ',f8.5)") hwa_i
!           write(outputunit,"('hwa_f = ',f8.5)") hwa_f
!           write(outputunit,"('hwa_npts = ',f8.5)") hwa_npts
          hwa_npt1 = hwa_npts + 1
          nparams2 = nparams2-1
!-------------------------------------------------------------------------------
        case("hwt")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwt_i
            if(index(istring1(i+3),"!").eq.0) then
              read(unit=istring1(i+3),fmt=*,iostat=ios) hwt_f
              read(unit=istring1(i+4),fmt=*,iostat=ios) hwt_npts
            end if
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwt_i
            if(index(istring1(i+1),"!").eq.0) then
              read(unit=istring2(i+1),fmt=*,iostat=ios) hwt_f
              read(unit=istring2(i+2),fmt=*,iostat=ios) hwt_npts
            end if
          end if
!           write(outputunit,"('hwt_i = ',f8.5)") hwt_i
!           write(outputunit,"('hwt_f = ',f8.5)") hwt_f
!           write(outputunit,"('hwt_npts = ',f8.5)") hwt_npts
          hwt_npt1 = hwt_npts + 1
          nparams2 = nparams2-1
        case("hwt=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwt_i
          if(index(istring1(i+2),"!").eq.0) then
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwt_f
            read(unit=istring1(i+3),fmt=*,iostat=ios) hwt_npts
          end if
!           write(outputunit,"('hwt_i = ',f8.5)") hwt_i
!           write(outputunit,"('hwt_f = ',f8.5)") hwt_f
!           write(outputunit,"('hwt_npts = ',f8.5)") hwt_npts
          hwt_npt1 = hwt_npts + 1
          nparams2 = nparams2-1
!-------------------------------------------------------------------------------
        case("hwp")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwp_i
            if(index(istring1(i+3),"!").eq.0) then
              read(unit=istring1(i+3),fmt=*,iostat=ios) hwp_f
              read(unit=istring1(i+4),fmt=*,iostat=ios) hwp_npts
            end if
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwp_i
            if(index(istring1(i+1),"!").eq.0) then
              read(unit=istring2(i+1),fmt=*,iostat=ios) hwp_f
              read(unit=istring2(i+2),fmt=*,iostat=ios) hwp_npts
            end if
          end if
!           write(outputunit,"('hwp_i = ',f8.5)") hwp_i
!           write(outputunit,"('hwp_f = ',f8.5)") hwp_f
!           write(outputunit,"('hwp_s = ',f8.5)") hwp_npts
          hwp_npt1 = hwp_npts + 1
          nparams2 = nparams2-1
        case("hwp=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwp_i
          if(index(istring1(i+2),"!").eq.0) then
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwp_f
            read(unit=istring1(i+3),fmt=*,iostat=ios) hwp_npts
          end if
!           write(outputunit,"('hwp_i = ',f8.5)") hwp_i
!           write(outputunit,"('hwp_f = ',f8.5)") hwp_f
!           write(outputunit,"('hwp_s = ',f8.5)") hwp_npts
          hwp_npt1 = hwp_npts + 1
          nparams2 = nparams2-1
!-------------------------------------------------------------------------------
        case("hwx")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwx
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwx
          end if
!           write(outputunit,"('hwx = ',f8.5)") hwx
          nparams2 = nparams2-1
        case("hwx=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwx
!           write(outputunit,"('hwx = ',f8.5)") hwx
          nparams2 = nparams2-1
!-------------------------------------------------------------------------------
        case("hwy")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwy
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwy
          end if
!           write(outputunit,"('hwy = ',f8.5)") hwy
          nparams2 = nparams2-1
        case("hwy=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwy
!           write(outputunit,"('hwy = ',f8.5)") hwy
          nparams2 = nparams2-1
!-------------------------------------------------------------------------------
        case("hwz")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwz
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwz
          end if
!           write(outputunit,"('hwz = ',f8.5)") hwz
          nparams2 = nparams2-1
        case("hwz=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwz
!           write(outputunit,"('hwz = ',f8.5)") hwz
          nparams2 = nparams2-1
!===============================================================================
        case("hwscale","hwscale:","hwscale=")
          do n=i,iomax
            read(unit=istring1(n+1),fmt=*,iostat=ios) hwscale(n-i+1)
            if(ios.ne.0) exit
          end do
!           write(outputunit,"('hwscale = ',20(es9.2,2x))") hwscale
          exit
!===============================================================================
        case("hwtrotate","hwtrotate:","hwtrotate=")
          do n=i,iomax
            read(unit=istring1(n+1),fmt=*,iostat=ios) hwtrotate(n-i+1)
            if(ios.ne.0) exit
          end do
!           write(outputunit,"('hwtrotate = ',20(es9.2,2x))") hwtrotate
          exit
!===============================================================================
        case("hwprotate","hwprotate:","hwprotate=")
          do n=i,iomax
            read(unit=istring1(n+1),fmt=*,iostat=ios) hwprotate(n-i+1)
            if(ios.ne.0) exit
          end do
!           write(outputunit,"('hwtrotate = ',20(es9.2,2x))") hwprotate
          exit
!===============================================================================
        case("dirEfield")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            dirEfield = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            dirEfield = istring2(i)
          end if
!           write(outputunit,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
        case("dirEfield=")
          dirEfield = istring1(i+1)
!           write(outputunit,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
!===============================================================================
        case("ncp")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) ncp
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) ncp
          end if
!           write(outputunit,"('ncp = ',i0,' x ',i0)") ncp
          nparams = nparams-1
        case("ncp=")
            read(unit=istring1(i+1),fmt=*,iostat=ios) ncp
!           write(outputunit,"('ncp = ',i0,' x ',i0)") ncp
          nparams = nparams-1
!===============================================================================
        case("parts")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) parts
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) parts
          end if
!           write(outputunit,"('parts = ',i0)") parts
          nparams = nparams-1
        case("parts=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) parts
!           write(outputunit,"('parts = ',i0)") parts
          nparams = nparams-1
!===============================================================================
        case("n1gl")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n1gl
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n1gl
          end if
!           write(outputunit,"('n1gl = ',i0)") n1gl
          nparams = nparams-1
        case("n1gl=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n1gl
!           write(outputunit,"('n1gl = ',i0)") n1gl
          nparams = nparams-1
!===============================================================================
        case("parts3")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) parts3
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) parts3
          end if
!           write(outputunit,"('parts3 = ',i0)") parts3
          nparams = nparams-1
        case("parts3=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) parts3
!           write(outputunit,"('parts3 = ',i0)") parts3
          nparams = nparams-1
!===============================================================================
        case("n3gl")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n3gl
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n3gl
          end if
!           write(outputunit,"('n3gl = ',i0)") n3gl
          nparams = nparams-1
        case("n3gl=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n3gl
!           write(outputunit,"('n3gl = ',i0)") n3gl
          nparams = nparams-1
!===============================================================================
        case("emin")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) emin
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) emin
          end if
!           write(outputunit,"('emin = ',f8.5)") emin
          nparams = nparams-1
        case("emin=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) emin
!           write(outputunit,"('emin = ',f8.5)") emin
          nparams = nparams-1
!===============================================================================
        case("emax")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) emax
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) emax
          end if
!           write(outputunit,"('emax = ',f8.5)") emax
          nparams = nparams-1
        case("emax=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) emax
!           write(outputunit,"('emax = ',f8.5)") emax
          nparams = nparams-1
!===============================================================================
        case("skip_steps")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) skip_steps
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) skip_steps
          end if
!           write(outputunit,"('skip_steps = ',i0)") skip_steps
        case("skip_steps=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) skip_steps
!           write(outputunit,"('skip_steps = ',i0)") skip_steps
!===============================================================================
        case("skip_steps_hw")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) skip_steps_hw
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) skip_steps_hw
          end if
!           write(outputunit,"('skip_steps_hw = ',i0)") skip_steps_hw
        case("skip_steps_hw=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) skip_steps_hw
!           write(outputunit,"('skip_steps_hw = ',i0)") skip_steps_hw
!===============================================================================
        case("npts")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) npts
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) npts
          end if
!           write(outputunit,"('npts = ',i0)") npts
          npt1 = npts + 1
          nparams = nparams-1
        case("npts=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) npts
!           write(outputunit,"('npts = ',i0)") npts
          npt1 = npts + 1
          nparams = nparams-1
!===============================================================================
        case("n0sc1")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n0sc1
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n0sc1
          end if
!           write(outputunit,"('n0sc1 = ',i0)") n0sc1
          nparams = nparams-1
        case("n0sc1=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n0sc1
!           write(outputunit,"('n0sc1 = ',i0)") n0sc1
          nparams = nparams-1
!===============================================================================
        case("n0sc2")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n0sc2
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n0sc2
          end if
!           write(outputunit,"('n0sc2 = ',i0)") n0sc2
          nparams = nparams-1
        case("n0sc2=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n0sc2
!           write(outputunit,"('n0sc2 = ',i0)") n0sc2
          nparams = nparams-1
!===============================================================================
        case("renorm")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) renorm
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) renorm
          end if
          if(renorm) then
!             write(outputunit,"('Current renormalization: ACTIVATED')")
            nparams2 = nparams2+1 ! renormnb
!           else
!             write(outputunit,"('Current renormalization: DEACTIVATED')")
          end if
          nparams = nparams-1
        case("renorm=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) renorm
          if(renorm) then
!             write(outputunit,"('Current renormalization: ACTIVATED')")
            nparams2 = nparams2+1 ! renormnb
          else
!             write(outputunit,"('Current renormalization: DEACTIVATED')")
          end if
          nparams = nparams-1
!===============================================================================
        case("renormnb")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) renormnb
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) renormnb
          end if
!           write(outputunit,"('renormnb = ',i0)") renormnb
          nparams2 = nparams2-1
        case("renormnb=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) renormnb
!           write(outputunit,"('renormnb = ',i0)") renormnb
          nparams2 = nparams2-1
!===============================================================================
        case("kdirection")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            kdirection = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            kdirection = istring2(i)
          end if
!           write(outputunit,"('kdirection = ',a2)") kdirection
          nparams = nparams-1
        case("kdirection=")
          kdirection = istring1(i+1)
!           write(outputunit,"('kdirection = ',a2)") kdirection
          nparams = nparams-1
!===============================================================================
        case("skipsc")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) skipsc
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) skipsc
          end if
!           if(skipsc) then
!             write(outputunit,"('Skip self-consistency: ACTIVATED')")
!           else
!             write(outputunit,"('Skip self-consistency: DEACTIVATED')")
!           end if
          nparams = nparams-1
        case("skipsc=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) skipsc
!           if(skipsc) then
!             write(outputunit,"('Skip self-consistency: ACTIVATED')")
!           else
!             write(outputunit,"('Skip self-consistency: DEACTIVATED')")
!           end if
          nparams = nparams-1
!===============================================================================
        case("scfile","scfile:","scfile=")
          if ((index(string1(j),"scfile:").gt.0).or.(index(string1(j),"scfile=").gt.0)) then
            filenameini = index(string1(j),"scfile")+8
          else if (index(string1(j),"scfile =").gt.0) then
            filenameini = index(string1(j),"scfile")+9
          else
            filenameini = index(string1(j),"scfile")+7
          end if
          scfile = string1(j)(filenameini:)
!===============================================================================
        case("dfttype")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            dfttype = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            dfttype = istring2(i)
          end if
!           write(outputunit,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
        case("dfttype=")
          dfttype = istring1(i+1)
!           write(outputunit,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
!===============================================================================
        case("set1")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) set1
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) set1
          end if
!           write(outputunit,"('set1 = ',i0)") set1
          nparams = nparams-1
        case("set1=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) set1
!           write(outputunit,"('set1 = ',i0)") set1
          nparams = nparams-1
!===============================================================================
        case("set2")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) set2
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) set2
          end if
!           write(outputunit,"('set2 = ',i0)") set2
          nparams = nparams-1
        case("set2=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) set2
!           write(outputunit,"('set2 = ',i0)") set2
          nparams = nparams-1
!===============================================================================
        case("addlayers","addlayers:")
          do n=i+1,iomax
            read(unit=istring1(n),fmt=*,iostat=ios) addlayers(naddlayers+1)
            if(ios.ne.0) exit
            naddlayers = naddlayers+1
          end do
!           write(outputunit,"('naddlayers = ',i0,' addlayers = ',10(i0,2x))") naddlayers,addlayers
          exit
!===============================================================================
        case("suffix","suffix:","suffix=")
          if ((index(string1(j),"suffix:").gt.0).or.(index(string1(j),"suffix=").gt.0)) then
            filenameini = index(string1(j),"suffix")+8
          else if (index(string1(j),"suffix =").gt.0) then
            filenameini = index(string1(j),"suffix")+9
          else
            filenameini = index(string1(j),"suffix")+7
          end if
          suffix = string1(j)(filenameini:)
!===============================================================================
        case("!")
          exit
        end select ioparams
      end do read_line
      istring1 = istring2
      istring2 = " "
    end do lines_loop
    deallocate(string1,stringt)
    if(nparams.ne.0) then
      if(myrank.eq.0) write(outputunit,"('[get_parameters] Missing ',i0,' parameter(s)!')") nparams
      call MPI_Finalize(ierr)
      stop
    end if
    if(nparams2.gt.0) then
      if(myrank.eq.0) write(outputunit,"('[get_parameters] Missing ',i0,' secondary parameter(s)!')") nparams2
      call MPI_Finalize(ierr)
      stop
    end if
    if(myrank.eq.0) write(outputunit,"('[get_parameters] Finished reading from ""',a,'"" file')") trim(filename)

!-------------------------------------------------------------------------------
!*********** User manual additions / modifications in the input file **********!
!     Npl_i  = 4
!     Npl_f = 4
!     ncp = 6
!     SOC = .true.
!     magaxis = "5"
!     runoptions = trim(runoptions) // " noUonNM"
!     scfile = "results/selfconsistency/selfconsistency_Npl=4_dfttype=T_parts=2_U= 0.7E-01_hwa= 0.00E+00_hwt= 0.00E+00_hwp= 0.00E+00_ncp=6_eta= 0.5E-03.dat"
!-------------------------------------------------------------------------------
    ! Some consistency checks
    if((renorm).and.((renormnb.lt.n0sc1).or.(renormnb.gt.n0sc2))) then
      if(myrank.eq.0) then
        write(outputunit,"('[get_parameters] Invalid neighbor for renormalization: ',i0,'!')") renormnb
        write(outputunit,"('[get_parameters] Choose a value between ',i0,' and ',i0,'.')") n0sc1,n0sc2
      end if
      call MPI_Finalize(ierr)
      stop
    end if
    if(skip_steps.lt.0) then
      if(myrank.eq.0) write(outputunit,"('[get_parameters] Invalid number of energy steps to skip: ',i0)") skip_steps
      call MPI_Finalize(ierr)
      stop
    end if
    if(skip_steps_hw.lt.0) then
      if(myrank.eq.0) write(outputunit,"('[get_parameters] Invalid number of field steps to skip: ',i0)") skip_steps_hw
      call MPI_Finalize(ierr)
      stop
    end if
    if((lhfresponses).and.(itype.eq.7).and.(myrank.eq.0)) write(outputunit,"('[get_parameters] Susceptibility calculations already include HF responses. Ignoring ""hfresponses"" runoption')")
    ! Adjusting zeeman energy to Ry or eV
    tesla = tesla*ry2ev

    ! Turning off renormalization for non-current calculations
    if(itype.ne.8) renorm = .false.
    n0sc=n0sc2-n0sc1+1 ! Total number of neighbors

    ! Setting up external field variables and loops
    call prepare_field()

    ! Energy loop step
    deltae = (emax - emin)/npts
    if(deltae.le.1.d-14) npt1 = 1

    ! Preparing dc-limit calculation
    if(itype.eq.9) call prepare_dclimit()

    ! Check number of planes
    if(Npl_f.lt.Npl_i) then
      Npl_f = Npl_i
    end if
    ! Add 'naddlayers' to Npl
    if((naddlayers.eq.1).and.(myrank.eq.0)) write(outputunit,"('[get_parameters] WARNING: Added layers must include empty spheres! Only including one layer: naddlayers = ',i0)") naddlayers
    if((set1.eq.9).or.(set2.eq.9)) then
      naddlayers = 0
    end if
    if(naddlayers.ne.0) then
      Npl_i = Npl_i+naddlayers-1
      Npl_f = Npl_f+naddlayers-1
    end if

    tol   = 1.d-10
    pn1=parts*n1gl
    pn2=parts3*n3gl
    pnt=pn1+pn2
    return
  end subroutine get_parameters

  subroutine iowrite()
    use mod_mpi_pars
    use mod_constants, only: pi
    implicit none
    integer :: j,err

#ifdef _OPENMP
    write(outputunit_loop,"(10x,'Running on ',i0,' MPI process(es) WITH openMP')") numprocs
#else
    write(outputunit_loop,"(10x,'Running on ',i0,' MPI process(es) WITHOUT openMP')") numprocs
#endif
    write(outputunit_loop,"('|------------------------------- PARAMETERS: -------------------------------|')")
    write(outputunit_loop,"(10x,'Npl = ',i0)") Npl
    write(outputunit_loop,"(1x,'DFT parameters: ',$)")
    dft_type: select case (dfttype)
    case ("T")
      write(outputunit_loop,"('Tight-binding basis')")
    case ("O")
      write(outputunit_loop,"('Orthogonal basis')")
    end select dft_type
    if(SOC) then
      write(outputunit_loop,"(1x,'Spin Orbit Coupling: ACTIVATED')")
      write(outputunit_loop,"(5x,'socscale =',es9.2)") socscale
    else
      write(outputunit_loop,"(1x,'Spin Orbit Coupling: DEACTIVATED')")
    end if
    write(outputunit_loop,"(8x,'Utype = ',i0)") Utype
    select case (lattice)
    case("bcc110")
      write(outputunit_loop,"(1x,'Spin quantization direction: ',$)")
      write_magnetization_axis_bcc110: select case (magaxis)
      case ("L")
        write(outputunit_loop,"('Long axis')")
      case ("S")
        write(outputunit_loop,"('Short axis')")
      case ("P")
        write(outputunit_loop,"('Perpendicular out-of-plane')")
      case default
        write(outputunit_loop,"(11x,'Other')")
      end select write_magnetization_axis_bcc110
      write(outputunit_loop,"(10x,'phi =',es9.2,'pi')") phi/pi
      write(outputunit_loop,"(8x,'theta =',es9.2,'pi')") theta/pi
      write(outputunit_loop,"(1x,'Electric field direction: ',$)")
      read(dirEfield,fmt=*,iostat=err) j
      if(err.eq.0) then
        write_direction_E_field_bcc110: select case (j)
        case (1:6)
          write(outputunit_loop,"('Neighbor ',i0)") j
        case default
          write(outputunit_loop,"('Other',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)
        end select write_direction_E_field_bcc110
      else
        write_direction_E_field_bcc110_axis: select case (dirEfield)
        case ("L")
          write(outputunit_loop,"('Long axis')")
        case ("S")
          write(outputunit_loop,"('Short axis')")
        case ("O")
          write(outputunit_loop,"('Other',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)
        end select write_direction_E_field_bcc110_axis
      end if
    case("fcc100")
      write(outputunit_loop,"(1x,'Spin quantization direction: ',$)")
      read(unit=magaxis,fmt=*) j
      write_magnetization_axis_fcc100: select case (j)
      case (1:8)   !    In plane neighbors:
        write(outputunit_loop,"('Neighbor ',i0)") j
      case (9)
        write(outputunit_loop,"('Out-of-plane')")
      end select write_magnetization_axis_fcc100
      write(outputunit_loop,"(10x,'phi =',es9.2,'pi')") phi/pi
      write(outputunit_loop,"(8x,'theta =',es9.2,'pi')") theta/pi
      write(outputunit_loop,"(1x,'Electric field direction: ',$)")
      read(unit=dirEfield,fmt=*) j
      write_direction_E_field_fcc100: select case (j)
      case (1:8)   !    In plane neighbors:
        write(outputunit_loop,"('neighbor ',i0)") j
      case default
        write(outputunit_loop,"('Other',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)
      end select write_direction_E_field_fcc100
    end select
    if(renorm) then
      write(outputunit_loop,"(1x,'Current renormalization: ACTIVATED')")
      write(outputunit_loop,"(5x,'renormnb = ',i0)") renormnb
    else
      write(outputunit_loop,"(1x,'Current renormalization: DEACTIVATED')")
    end if
    write(outputunit_loop,"(10x,'ncp = ',i0)") ncp
    write(outputunit_loop,"(8x,'parts = ',i0,'x',i0)") parts,n1gl
    write(outputunit_loop,"(7x,'parts3 = ',i0,'x',i0)") parts3,n3gl
    write(outputunit_loop,"(10x,'eta =',es9.2)") eta
    if(lfield) then
      write(outputunit_loop,"(1x,'Static magnetic field: ACTIVATED')")
      write(outputunit_loop,"(10x,'hwx =',es9.2,5x,'|',5x,'hwa =',es9.2)") hwx,hw_list(hw_count,1)
      write(outputunit_loop,"(10x,'hwy =',es9.2,5x,'|',5x,'hwt =',f6.3)") hwy,hw_list(hw_count,2)
      write(outputunit_loop,"(10x,'hwz =',es9.2,5x,'|',5x,'hwp =',f6.3)") hwz,hw_list(hw_count,3)
    else
      write(outputunit_loop,"(1x,'Static magnetic field: DEACTIVATED')")
    end if
    if(runoptions.ne."") write(outputunit_loop,"(6x,'Activated options:',/,4x,a)") trim(runoptions)

    write(outputunit_loop,"('|------------------------------ TO CALCULATE: ------------------------------|')")
    write_itype: select case (itype)
    case (0)
      write(outputunit_loop,"(1x,'Test before SC')")
      write(outputunit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
      write(outputunit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
      write(outputunit_loop,"(9x,'emin =',es9.2)") emin
      write(outputunit_loop,"(9x,'emax =',es9.2)") emax
      write(outputunit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (1)
      write(outputunit_loop,"(1x,'Self-consistency only')")
    case (2)
      write(outputunit_loop,"(1x,'Test after SC')")
      write(outputunit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
      write(outputunit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
      write(outputunit_loop,"(9x,'emin =',es9.2)") emin
      write(outputunit_loop,"(9x,'emax =',es9.2)") emax
      write(outputunit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (3)
      write(outputunit_loop,"(1x,'LDOS and exchange interactions as a function of energy')")
      write(outputunit_loop,"(9x,'emin =',es9.2)") emin
      write(outputunit_loop,"(9x,'emax =',es9.2)") emax
      write(outputunit_loop,"(1x,'Number of points to calculate: ',i0)") npt1
    case (4)
      write(outputunit_loop,"(1x,'Band structure')")
      write(outputunit_loop,"(3x,'kdirection = ',a)") kdirection
      write(outputunit_loop,"(9x,'Number of points to calculate: ',i0)") npt1
    case (5)
      write(outputunit_loop,"(1x,'Charge and spin density at Fermi surface')")
    case (6)
      write(outputunit_loop,"(1x,'Exhange interactions and anisotropies (full tensor)')")
      if(nmaglayers.eq.1) write(outputunit_loop,"(1x,'Only 1 magnetic layer: calculating only anisotropies')")
      write(outputunit_loop,"(8x,'from Npl = ',i0,' to ',i0)") Npl_i,Npl_f
    case (7)
      write(outputunit_loop,"(1x,'Local susceptibility as a function of energy')")
      write(outputunit_loop,"(9x,'emin =',es9.2)") emin
      write(outputunit_loop,"(9x,'emax =',es9.2)") emax
      write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps of size',es10.3,' each calculating ',i0,' points')") npt1,MPIsteps,MPIdelta,MPIpts
    case (8)
      write(outputunit_loop,"(1x,'Parallel currents, disturbances and local susc. as a function of energy')")
      write(outputunit_loop,"(8x,'n0sc1 = ',i0)") n0sc1
      write(outputunit_loop,"(8x,'n0sc2 = ',i0)") n0sc2
      write(outputunit_loop,"(9x,'emin =',es9.2)") emin
      write(outputunit_loop,"(9x,'emax =',es9.2)") emax
      if(lvdc) then
        write(outputunit_loop,"(7x,'DC voltage calculation: ')")
        if(vdcneighbor(1).ne.0) write(outputunit_loop,"(4x,'Longitudinal neighbor = ',i0)") vdcneighbor(1)
        if(vdcneighbor(2).ne.0) write(outputunit_loop,"(4x,'  Transverse neighbor = ',i0)") vdcneighbor(2)
      end if
      write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps of energy size',es10.3,' each calculating ',i0,' points')") total_hw_npt1*npt1,MPIsteps*MPIsteps_hw,MPIdelta,MPIpts_hw*MPIpts
    case (9)
      write(outputunit_loop,"(1x,'dc limit calculations as a function of ',a)") trim(dcfield(dcfield_dependence))
      write(outputunit_loop,"(1x,'e =',es9.2)") emin
      write(outputunit_loop,"(1x,'hwa_min =',es9.2)") hw_list(1,1)
      write(outputunit_loop,"(1x,'hwa_max =',es9.2)") hw_list(total_hw_npt1,1)
      write(outputunit_loop,"(1x,'hwt_min =',f6.3)") hw_list(1,2)
      write(outputunit_loop,"(1x,'hwt_max =',f6.3)") hw_list(total_hw_npt1,2)
      write(outputunit_loop,"(1x,'hwp_min =',f6.3)") hw_list(1,3)
      write(outputunit_loop,"(1x,'hwp_max =',f6.3)") hw_list(total_hw_npt1,3)
      if(lvdc) then
        write(outputunit_loop,"(7x,'DC voltage calculation: ')")
        if(vdcneighbor(1).ne.0) write(outputunit_loop,"(8x,'Longitudinal neighbor = ',i0)") vdcneighbor(1)
        if(vdcneighbor(2).ne.0) write(outputunit_loop,"(8x,'  Transverse neighbor = ',i0)") vdcneighbor(2)
      end if
      write(outputunit_loop,"(1x,i0,' points divided into ',i0,' steps, each calculating ',i0,' points')") total_hw_npt1*npt1,MPIsteps*MPIsteps_hw,MPIpts_hw*MPIpts
    end select write_itype
    write(outputunit_loop,"('|---------------------------------------------------------------------------|')")
    return
  end subroutine iowrite

end module mod_io
