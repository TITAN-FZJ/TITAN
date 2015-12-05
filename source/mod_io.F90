module mod_io
  use mod_f90_kind
  use mod_parameters
  implicit none
contains
  subroutine ioread()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_magnet, only: hhwx,hhwy,hhwz
    implicit none
    integer, parameter      :: iomax=20           ! Maximum number of elements in one line
    integer                 :: nparams=26         ! Number of required parameters to be read
    integer                 :: nparams2=0         ! Number of secondary parameters required
    integer                 :: ifile,ios,i,j,n,nlines,nlinest
    character(len=20)       :: istring1(iomax),istring2(iomax),stringtemp,filename
    character(len=200),allocatable      :: string1(:),stringt(:)

    ifile=666666
    filename="inputdhe"
    open(unit=ifile, file=filename, status='old', iostat=ios)
    if(ios.ne.0) then
      write(*,"('[ioread] File ',a,' not found!')") filename
      stop
    end if

    ! Counting the number of lines
    nlines  = 0
    nlinest = 0
    do
      read (ifile,*,iostat=ios) stringtemp
      if (ios.ne.0) exit
      if (stringtemp.eq."") cycle ! If the line is blank, ignore
      ! Total number of non-empty lines
      nlinest = nlinest + 1

      ! Getting the number of non-commented lines
      if ((stringtemp(1:1).eq."#").or.(stringtemp(1:1).eq."!")) cycle
      nlines = nlines + 1
    end do
!     write(*,"('[ioread] ""',a,'"" file has ',i0,' non-commented lines and ',i0,' in total (non-blank only)')") trim(filename),nlines,nlinest
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
!       write(*,*) string1(j)
    end do

    ! Getting variables from non-commented lines
    istring1 = " "
    istring2 = " "
    scfile = " "
    runoptions = " "
    if(myrank.eq.0) write(*,"('[ioread] Reading parameters from ""',a,'"" file...')") trim(filename)
    read(unit=string1(1),fmt=*,iostat=ios) (istring1(i),i=1,iomax) ! Reads first line
    lines_loop: do j=1,nlines
      ! Reads next line for variables defined in two lines
      if(j.ne.nlines) read(unit=string1(j+1),fmt=*,iostat=ios) (istring2(i),i=1,iomax)
      read_line: do i=1,iomax
        ioparams: select case (istring1(i))
!===============================================================================
        case("itype")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) itype
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) itype
          end if
!           write(*,"('itype = ',i0)") itype
          nparams = nparams-1
        case("itype=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) itype
!           write(*,"('itype = ',i0)") itype
          nparams = nparams-1
!===============================================================================
        case("Options","Options:")
          do n=i+1,iomax
            if(istring1(n).eq."") cycle
            options: select case (istring1(n))
            case ("ry2ev","verbose","idebug","addresults","createfiles","noUonall","noUonNM","slatec","GSL","kpoints")
              runoptions = trim(runoptions) // " " // trim(istring1(n))
!               write(*,"('Option ""',a,'"" activated.')") trim(istring1(n))
!             case default
!               write(*,"('Option ""',a,'"" not recognized. Ignoring.')") trim(istring1(n))
            case("!")
              exit
            end select options
          end do
          if(index(runoptions,"ry2ev").gt.0) then
            ry2ev = 13.6d0
          else
            ry2ev = 1.d0
          end if
!           write(*,*) (istring1(n),n=i+1,iomax)
          exit
!===============================================================================
        case("eta")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) eta
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) eta
          end if
!           write(*,"('eta = ',f8.5)") eta
          nparams = nparams-1
        case("eta=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) eta
!           write(*,"('eta = ',f8.5)") eta
          nparams = nparams-1
!===============================================================================
        case("lattice")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            lattice = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            lattice = istring2(i)
          end if
!           write(*,"('lattice = ',a2)") lattice
          nparams = nparams-1
        case("lattice=")
          lattice = istring1(i+1)
!           write(*,"('lattice = ',a2)") lattice
          nparams = nparams-1
!===============================================================================
        case("Npl")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) Nplini
            read(unit=istring1(i+3),fmt=*,iostat=ios) Nplfinal
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) Nplini
            read(unit=istring2(i+1),fmt=*,iostat=ios) Nplfinal
          end if
!           write(*,"('Nplini = ',i0)") Nplini
!           write(*,"('Nplfinal = ',i0)") Nplfinal
          nparams = nparams-1
        case("Npl=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) Nplini
          read(unit=istring1(i+2),fmt=*,iostat=ios) Nplfinal
!           write(*,"('Nplini = ',i0)") Nplini
!           write(*,"('Nplfinal = ',i0)") Nplfinal
          nparams = nparams-1
!===============================================================================
        case("a0")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) a0
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) a0
          end if
!           write(*,"('a0 = ',f8.5)") a0
          nparams = nparams-1
        case("a0=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) a0
!           write(*,"('a0 = ',f8.5)") a0
          nparams = nparams-1
!===============================================================================
        case("SOC")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) SOC
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) SOC
          end if
          if(SOC) then
!             write(*,"('Spin Orbit Coupling: ACTIVATED')")
            nparams2 = nparams2+2 ! magaxis and socscale
          else
!             write(*,"('Spin Orbit Coupling: DEACTIVATED')")
          end if
          nparams = nparams-1
        case("SOC=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) SOC
          if(SOC) then
!             write(*,"('Spin Orbit Coupling: ACTIVATED')")
            nparams2 = nparams2+2 ! magaxis and socscale
          else
!             write(*,"('Spin Orbit Coupling: DEACTIVATED')")
          end if
          nparams = nparams-1
!-------------------------------------------------------------------------------
        case("socscale")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) socscale
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) socscale
          end if
!           write(*,"('socscale = ',f8.5)") socscale
          nparams2 = nparams2-1
        case("socscale=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) socscale
!           write(*,"('socscale = ',f8.5)") socscale
          nparams2 = nparams2-1
!-------------------------------------------------------------------------------
        case("magaxis")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            magaxis = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            magaxis = istring2(i)
          end if
!           write(*,"('magaxis = ',a2)") magaxis
          nparams2 = nparams2-1
        case("magaxis=")
          magaxis = istring1(i+1)
!           write(*,"('magaxis = ',a2)") magaxis
          nparams2 = nparams2-1
!===============================================================================
        case("hwx")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwx
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwx
          end if
!           write(*,"('hwx = ',f8.5)") hwx
          nparams = nparams-1
        case("hwx=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwx
!           write(*,"('hwx = ',f8.5)") hwx
          nparams = nparams-1
!-------------------------------------------------------------------------------
        case("hwy")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwy
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwy
          end if
!           write(*,"('hwy = ',f8.5)") hwy
          nparams = nparams-1
        case("hwy=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwy
!           write(*,"('hwy = ',f8.5)") hwy
          nparams = nparams-1
!-------------------------------------------------------------------------------
        case("hwz")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) hwz
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) hwz
          end if
!           write(*,"('hwz = ',f8.5)") hwz
          nparams = nparams-1
        case("hwz=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) hwz
!           write(*,"('hwz = ',f8.5)") hwz
          nparams = nparams-1
!===============================================================================
        case("dirEfield")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            dirEfield = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            dirEfield = istring2(i)
          end if
!           write(*,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
        case("dirEfield=")
          dirEfield = istring1(i+1)
!           write(*,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
!===============================================================================
        case("ncp")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) ncp
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) ncp
          end if
!           write(*,"('ncp = ',i0,' x ',i0)") ncp
          nparams = nparams-1
        case("ncp=")
            read(unit=istring1(i+1),fmt=*,iostat=ios) ncp
!           write(*,"('ncp = ',i0,' x ',i0)") ncp
          nparams = nparams-1
!===============================================================================
        case("parts")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) parts
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) parts
          end if
!           write(*,"('parts = ',i0)") parts
          nparams = nparams-1
        case("parts=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) parts
!           write(*,"('parts = ',i0)") parts
          nparams = nparams-1
!===============================================================================
        case("n1gl")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n1gl
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n1gl
          end if
!           write(*,"('n1gl = ',i0)") parts
          nparams = nparams-1
        case("n1gl=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n1gl
!           write(*,"('n1gl = ',i0)") parts
          nparams = nparams-1
!===============================================================================
        case("parts3")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) parts3
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) parts3
          end if
!           write(*,"('parts3 = ',i0)") parts3
          nparams = nparams-1
        case("parts3=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) parts3
!           write(*,"('parts3 = ',i0)") parts3
          nparams = nparams-1
!===============================================================================
        case("n3gl")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n3gl
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n3gl
          end if
!           write(*,"('n3gl = ',i0)") parts
          nparams = nparams-1
        case("n3gl=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n3gl
!           write(*,"('n3gl = ',i0)") parts
          nparams = nparams-1
!===============================================================================
        case("emin")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) emin
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) emin
          end if
!           write(*,"('emin = ',f8.5)") emin
          nparams = nparams-1
        case("emin=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) emin
!           write(*,"('emin = ',f8.5)") emin
          nparams = nparams-1
!===============================================================================
        case("emax")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) emax
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) emax
          end if
!           write(*,"('emax = ',f8.5)") emax
          nparams = nparams-1
        case("emax=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) emax
!           write(*,"('emax = ',f8.5)") emax
          nparams = nparams-1
!===============================================================================
        case("npts")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) npts
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) npts
          end if
!           write(*,"('npts = ',i0)") npts
          npt1 = npts + 1
          nparams = nparams-1
        case("npts=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) npts
!           write(*,"('npts = ',i0)") npts
          npt1 = npts + 1
          nparams = nparams-1
!===============================================================================
        case("n0sc1")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n0sc1
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n0sc1
          end if
!           write(*,"('n0sc1 = ',i0)") n0sc1
          nparams = nparams-1
        case("n0sc1=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n0sc1
!           write(*,"('n0sc1 = ',i0)") n0sc1
          nparams = nparams-1
!===============================================================================
        case("n0sc2")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) n0sc2
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) n0sc2
          end if
!           write(*,"('n0sc2 = ',i0)") n0sc2
          nparams = nparams-1
        case("n0sc2=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) n0sc2
!           write(*,"('n0sc2 = ',i0)") n0sc2
          nparams = nparams-1
!===============================================================================
        case("renorm")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) renorm
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) renorm
          end if
          if(renorm) then
!             write(*,"('Current renormalization: ACTIVATED')")
            nparams2 = nparams2+1 ! renormnb
!           else
!             write(*,"('Current renormalization: DEACTIVATED')")
          end if
          nparams = nparams-1
        case("renorm=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) renorm
          if(renorm) then
!             write(*,"('Current renormalization: ACTIVATED')")
            nparams2 = nparams2+1 ! renormnb
          else
!             write(*,"('Current renormalization: DEACTIVATED')")
          end if
          nparams = nparams-1
!===============================================================================
        case("renormnb")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) renormnb
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) renormnb
          end if
!           write(*,"('renormnb = ',i0)") renormnb
          nparams2 = nparams2-1
        case("renormnb=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) renormnb
!           write(*,"('renormnb = ',i0)") renormnb
          nparams2 = nparams2-1
!===============================================================================
        case("kdirection")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            kdirection = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            kdirection = istring2(i)
          end if
!           write(*,"('kdirection = ',a2)") kdirection
          nparams = nparams-1
        case("kdirection=")
          kdirection = istring1(i+1)
!           write(*,"('kdirection = ',a2)") kdirection
          nparams = nparams-1
!===============================================================================
        case("skipsc")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) skipsc
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) skipsc
          end if
!           if(skipsc) then
!             write(*,"('Skip self-consistency: ACTIVATED')")
!           else
!             write(*,"('Skip self-consistency: DEACTIVATED')")
!           end if
          nparams = nparams-1
        case("skipsc=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) skipsc
!           if(skipsc) then
!             write(*,"('Skip self-consistency: ACTIVATED')")
!           else
!             write(*,"('Skip self-consistency: DEACTIVATED')")
!           end if
          nparams = nparams-1
!===============================================================================
        case("scfile","scfile:","scfile=")
          if ((index(string1(j),"scfile:").gt.0).or.(index(string1(j),"scfile=").gt.0)) then
            scfileini = index(string1(j),"scfile")+8
          else if (index(string1(j),"scfile =").gt.0) then
            scfileini = index(string1(j),"scfile")+9
          else
            scfileini = index(string1(j),"scfile")+7
          end if
          scfile = string1(j)(scfileini:)
          scfile = trim(scfile)
          open(unit=99,file=scfile,status="old",iostat=ios)
          if(ios.ne.0) then
            if(myrank.eq.0) write(*,"('*** WARNING: Self-consistency file given on input file does not exist! Using default... ***')")
            scfile = " "
          end if
          close(99)
!===============================================================================
        case("dfttype")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            dfttype = istring1(i+2)
          else ! If there's no '=' after keyword, get from next line
            dfttype = istring2(i)
          end if
!           write(*,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
        case("dfttype=")
          dfttype = istring1(i+1)
!           write(*,"('dirEfield = ',a2)") dirEfield
          nparams = nparams-1
!===============================================================================
        case("set1")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) set1
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) set1
          end if
!           write(*,"('set1 = ',i0)") set1
          nparams = nparams-1
        case("set1=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) set1
!           write(*,"('set1 = ',i0)") set1
          nparams = nparams-1
!===============================================================================
        case("set2")
          if(istring1(i+1).eq."=") then ! If after keyword there's an '='
            read(unit=istring1(i+2),fmt=*,iostat=ios) set2
          else ! If there's no '=' after keyword, get from next line
            read(unit=istring2(i),fmt=*,iostat=ios) set2
          end if
!           write(*,"('set2 = ',i0)") set2
          nparams = nparams-1
        case("set2=")
          read(unit=istring1(i+1),fmt=*,iostat=ios) set2
!           write(*,"('set2 = ',i0)") set2
          nparams = nparams-1
!===============================================================================
        case("!")
          exit
        end select ioparams
      end do read_line
      istring1 = istring2
      istring2 = " "
    end do lines_loop
    close(ifile)
    if(nparams.ne.0) then
      if(myrank.eq.0) write(*,"('[ioread] Missing ',i0,' parameter(s)!')") nparams
      stop
    end if
    if(nparams2.gt.0) then
      if(myrank.eq.0) write(*,"('[ioread] Missing ',i0,' secondary parameter(s)!')") nparams2
      stop
    end if
    if(myrank.eq.0) write(*,"('[ioread] Finished reading from ""',a,'"" file')") trim(filename)
  deallocate(string1,stringt)
  if((renorm).and.((renormnb.lt.n0sc1).or.(renormnb.gt.n0sc2))) then
    if(myrank.eq.0) then
      write(*,"('[ioread] Invalid neighbor for renormalization: ',i0,'!')") renormnb
      write(*,"('[ioread] Choose a value between ',i0,' and ',i0,'.')") n0sc1,n0sc2
    end if
    stop
  end if
!-------------------------------------------------------------------------------
!*********** User manual additions / modifications in the input file **********!
!   Nplini  = 4
!   Nplfinal = 4
!   ncp = 6
!   SOC = .true.
!   magaxis = "5"
!   runoptions = trim(runoptions) // " noUonNM"
!   scfile = "results/selfconsistency/selfconsistency_Npl=4_dfttype=T_parts=2_U= 0.7E-01_hwx= 0.0E+00_hwy= 0.0E+00_hwz= 0.0E+00_ncp=6_eta= 0.5E-03.dat"
!-------------------------------------------------------------------------------
  Utype = 2
  if(index(runoptions,"noUonall").ne.0) then
    Utype = 0
  else if(index(runoptions,"noUonNM").ne.0) then
    Utype = 1
  end if
  if((itype.ne.7).and.(itype.ne.8)) renorm = .false.
  n0sc=n0sc2-n0sc1+1
  deltae = (emax - emin)/npts
  if(deltae.le.1.d-14) npt1 = 1
!   if((itype.ne.1).and.(itype.ne.10)) then
!     Nplfinal = Nplini
!   else
    if(Nplfinal.lt.Nplini) then
!       if(myrank.eq.0) then
!         write(*,"('[ioread] Invalid number of planes for exchange calculations: ')")
!         write(*,"('[ioread] Nplini = ',i0,' , Nplfinal = ',i0)") Nplini,Nplfinal
!         write(*,"('[ioread] Using Npl = ',i0)") Nplini
!       end if
      Nplfinal = Nplini
    end if
!   end if
  hhwx  = 0.5d0*hwx
  hhwy  = 0.5d0*hwy
  hhwz  = 0.5d0*hwz
  tol   = 1.d-8
  pn1=parts*n1gl
  pn2=parts3*n3gl
  pnt=pn1+pn2
  end subroutine ioread

  subroutine iowrite()
    use mod_mpi_pars
    use mod_constants, only: pi
    implicit none
    integer :: j

#ifdef _OPENMP
    write(*,"(10x,'Running on ',i0,' MPI process(es) WITH openMP')") numprocs
#else
    write(*,"(10x,'Running on ',i0,' MPI process(es) WITHOUT openMP')") numprocs
#endif
    write(*,"('|------------------------------- PARAMETERS: -------------------------------|')")
    write(*,"(10x,'Npl = ',i0)") Npl
    write(*,"(1x,'DFT parameters: ',$)")
    dft_type: select case (dfttype)
    case ("T")
      write(*,"('Tight-binding basis')")
    case ("O")
      write(*,"('Orthogonal basis')")
    end select dft_type
    if(SOC) then
      write(*,"(1x,'Spin Orbit Coupling: ACTIVATED')")
      write(*,"(5x,'socscale =',e9.2)") socscale
    else
      write(*,"(1x,'Spin Orbit Coupling: DEACTIVATED')")
    end if
    write(*,"(8x,'Utype = ',i0)") Utype
    select case (lattice)
    case("bcc110")
      write(*,"(1x,'Spin quantization direction: ',$)")
      write_magnetization_axis_bcc110: select case (magaxis)
      case ("L")
        write(*,"('Long axis')")
      case ("S")
        write(*,"('Short axis')")
      case ("P")
        write(*,"('Perpendicular out-of-plane')")
      case default
        write(*,"(11x,'Other')")
      end select write_magnetization_axis_bcc110
      write(*,"(10x,'phi =',e9.2,'pi')") phi/pi
      write(*,"(8x,'theta =',e9.2,'pi')") theta/pi
      write(*,"(1x,'Electric field direction: ',$)")
      write_direction_E_field_bcc110: select case (dirEfield)
      case ("L")
        write(*,"('Long axis')")
      case ("S")
        write(*,"('Short axis')")
      case ("O")
        write(*,"('Other',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)
      end select write_direction_E_field_bcc110
    case("fcc100")
      write(*,"(1x,'Spin quantization direction: ',$)")
      read(unit=magaxis,fmt=*) j
      write_magnetization_axis_fcc100: select case (j)
      case (1:8)   !    In plane neighbors:
        write(*,"('Neighbor ',i0)") j
      case (9)
        write(*,"('Out-of-plane')")
      end select write_magnetization_axis_fcc100
      write(*,"(10x,'phi =',e9.2,'pi')") phi/pi
      write(*,"(8x,'theta =',e9.2,'pi')") theta/pi
      write(*,"(1x,'Electric field direction: ',$)")
      read(unit=dirEfield,fmt=*) j
      write_direction_E_field_fcc100: select case (j)
      case (1:8)   !    In plane neighbors:
        write(*,"('neighbor ',i0)") j
      case default
        write(*,"('Other',/,1x,' E = (',f6.3,',',f6.3,',',f6.3,')')") dirEfieldvec(1),dirEfieldvec(2),dirEfieldvec(3)
      end select write_direction_E_field_fcc100
    end select
    if(renorm) then
      write(*,"(1x,'Current renormalization: ACTIVATED')")
      write(*,"(5x,'renormnb = ',i0)") renormnb
    else
      write(*,"(1x,'Current renormalization: DEACTIVATED')")
    end if
    write(*,"(10x,'ncp = ',i0)") ncp
    write(*,"(8x,'parts = ',i0,'x',i0)") parts,n1gl
    write(*,"(7x,'parts3 = ',i0,'x',i0)") parts3,n3gl
    write(*,"(10x,'eta =',e9.2)") eta
    write(*,"(10x,'hwx =',e9.2)") hwx
    write(*,"(10x,'hwy =',e9.2)") hwy
    write(*,"(10x,'hwz =',e9.2)") hwz
    if(runoptions.ne."") write(*,"(6x,'Activated options:',/,4x,a)") trim(runoptions)

    write(*,"('|------------------------------ TO CALCULATE: ------------------------------|')")
    write_itype: select case (itype)
    case (0)
      write(*,"(1x,'Test before SC')")
      write(*,"(8x,'n0sc1 = ',i0)") n0sc1
      write(*,"(8x,'n0sc2 = ',i0)") n0sc2
      write(*,"(9x,'emin =',e9.2)") emin
      write(*,"(9x,'emax =',e9.2)") emax
      write(*,"(1x,'Number of points to calculate: ',i0)") npt1
    case (1)
      write(*,"(1x,'Self-consistency only')")
    case (2)
      write(*,"(1x,'Test after SC')")
      write(*,"(8x,'n0sc1 = ',i0)") n0sc1
      write(*,"(8x,'n0sc2 = ',i0)") n0sc2
      write(*,"(9x,'emin =',e9.2)") emin
      write(*,"(9x,'emax =',e9.2)") emax
      write(*,"(1x,'Number of points to calculate: ',i0)") npt1
    case (3)
      write(*,"(1x,'LDOS and exchange interactions as a function of energy')")
      write(*,"(9x,'emin =',e9.2)") emin
      write(*,"(9x,'emax =',e9.2)") emax
      write(*,"(1x,'Number of points to calculate: ',i0)") npt1
    case (4)
      write(*,"(1x,'Band structure')")
      write(*,"(3x,'kdirection = ',a)") kdirection
      write(*,"(9x,'Number of points to calculate: ',i0)") npt1
    case (5)
      write(*,"(1x,'Local susceptibility as a function of energy')")
      write(*,"(9x,'emin =',e9.2)") emin
      write(*,"(9x,'emax =',e9.2)") emax
      write(*,"(1x,i0,' points divided into ',i0,' steps of size',e10.3,' each calculating ',i0,' points')") npt1,MPIsteps,MPIdelta,MPIpts
    case (6)
      write(*,"(1x,'Disturbances and local susceptibility as a function of energy')")
      write(*,"(9x,'emin =',e9.2)") emin
      write(*,"(9x,'emax =',e9.2)") emax
      write(*,"(1x,i0,' points divided into ',i0,' steps of size',e10.3,' each calculating ',i0,' points')") npt1,MPIsteps,MPIdelta,MPIpts
    case (7)
      write(*,"(1x,'Parallel currents and local susceptibility as a function of energy')")
      write(*,"(8x,'n0sc1 = ',i0)") n0sc1
      write(*,"(8x,'n0sc2 = ',i0)") n0sc2
      write(*,"(9x,'emin =',e9.2)") emin
      write(*,"(9x,'emax =',e9.2)") emax
      write(*,"(1x,i0,' points divided into ',i0,' steps of size',e10.3,' each calculating ',i0,' points')") npt1,MPIsteps,MPIdelta,MPIpts
    case (8)
      write(*,"(1x,'Parallel currents, disturbances and local susc. as a function of energy')")
      write(*,"(8x,'n0sc1 = ',i0)") n0sc1
      write(*,"(8x,'n0sc2 = ',i0)") n0sc2
      write(*,"(9x,'emin =',e9.2)") emin
      write(*,"(9x,'emax =',e9.2)") emax
      write(*,"(1x,i0,' points divided into ',i0,' steps of size',e10.3,' each calculating ',i0,' points')") npt1,MPIsteps,MPIdelta,MPIpts
    case (9)
      write(*,"(1x,'Fermi surface')")
    case (10)
      write(*,"(1x,'Exhange interactions and anisotropies (full tensor)')")
      if(nmaglayers.eq.0) then
        write(*,"(1x,'[main] No magnetic layers!')")
        stop
      end if
      if(nmaglayers.eq.1) write(*,"(1x,'Only 1 magnetic layer: calculating only anisotropies')")
      write(*,"(8x,'from Npl = ',i0,' to ',i0)") Nplini,Nplfinal
    case (11)
      write(*,"(1x,'Probability of spin-flip as a function of position')")
    end select write_itype
    write(*,"('|---------------------------------------------------------------------------|')")
  end subroutine iowrite

end module mod_io