module mod_parameters
	use mod_f90_kind
	implicit none
	integer 			:: ncp
	real(double)	:: eta,hwx,hwy,hwz
	real(double)	:: Ef,q(2)
	integer 			:: dim,dimNpl
!========================================================================================!
! Effective intra-site electron electron interaction
	real(double),allocatable	:: U(:)
	integer				:: Utype
!========================================================================================!
! Lattice and surface direction
	character(len=6) :: lattice
	integer 			:: Npl
!	Lattice parameter (define the units of distance in the program)
	real(double) :: a0
!========================================================================================!
!	Equilibrium magnetization and direction of in-plane applied electric field:
	character(len=1)	:: magaxis,dirEfield
	real(double)	:: theta,phi 					! Euler Angles for the magnetization direction
	real(double)	:: dirEfieldvec(3)		! Direction vector of the electric field
!========================================================================================!
!	type of calculation - defined in input file 'inputshe'
	integer	:: itype
!========================================================================================!
!	Turn on/off SOC
	logical	:: SOC
!========================================================================================!
!	Number of parts to divide energy integral I1+I2 and I3
	integer :: parts,parts3
	integer :: n1gl,n3gl
!========================================================================================!
!	Band structure
	character(len=2)  :: kdirection
!========================================================================================!
!	Number of points and interval of energy/wave vector/position calculations
	integer	:: npts,npt1
	real(double)		:: emin,emax,deltae
	real(double)		:: qxmin,qxmax,qzmin,qzmax
!========================================================================================!
!	Conversion arrays
	integer,allocatable	:: sigmaimunu2i(:,:,:,:),sigmaijmunu2i(:,:,:,:,:)
!========================================================================================!
! Run options are stored in this string
	character(len=100)      		:: runoptions
	real(double)						 		:: ry2ev 							! Optional conversion of ry to eV
!========================================================================================!
!	Activate debug options
	logical	:: verbose
	logical	:: idebug
!========================================================================================!
!	n0sc1 - first neighbor to calculate the in-plane spin and charge current
!	n0sc2 - last neighbor to calculate the in-plane spin and charge current
!	n0sc  - Number of neighbors to calculate currents
	integer	:: n0sc1,n0sc2,n0sc
!========================================================================================!
!	Current renormalization
	logical	:: renorm
	integer	:: renormnb
!========================================================================================!
!	Skip self-consistency
	logical	:: skipsc
!	Give a file to start self-consistency
	character(len=200)      		:: scfile
	integer :: scfileini
!========================================================================================!
! Choose between tight-binding (T) or orthogonal (O) DFT parameters
  character(len=1)	:: dfttype
!========================================================================================!
! Hostname of rank
	character(len=50)	:: host
!========================================================================================!
! Set of tight-binding parameters to be used
! in the first half (set1) and second half (set2) of the slab
! NOTE: the Fermi energy is obtained from the first half.
  integer	:: set1,set2
!========================================================================================!


contains
	subroutine ioread()
		use mod_f90_kind
		use mod_mpi_pars
		use MPI
		implicit none
		integer, parameter 			:: iomax=20 					! Maximum number of elements in one line
		integer						 			:: nparams=26					! Number of required parameters to be read
		integer						 			:: nparams2=0					! Number of secondary parameters required
		integer									:: ifile,ios,i,j,n,nlines,nlinest
		character(len=20)       :: istring1(iomax),istring2(iomax),stringtemp,filename
		character(len=200),allocatable      :: string1(:),stringt(:)

		ifile=666666
		filename="inputshe"
		open(unit=ifile, file=filename, status='old', iostat=ios)

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
! 		write(*,"('[ioread] ""',a,'"" file has ',i0,' non-commented lines and ',i0,' in total (non-blank only)')") trim(filename),nlines,nlinest
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
!  	    write(*,*) string1(j)
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
! 					write(*,"('itype = ',i0)") itype
					nparams = nparams-1
				case("itype=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) itype
! 					write(*,"('itype = ',i0)") itype
					nparams = nparams-1
!===============================================================================
				case("Options","Options:")
					do n=i+1,iomax
						if(istring1(n).eq."") cycle
						options: select case (istring1(n))
						case ("ry2ev","verbose","idebug","addresults","createfiles","noUonall","noUonNM","slatec","GSL")
							runoptions = trim(runoptions) // " " // trim(istring1(n))
! 							write(*,"('Option ""',a,'"" activated.')") trim(istring1(n))
! 						case default
! 							write(*,"('Option ""',a,'"" not recognized. Ignoring.')") trim(istring1(n))
						case("!")
							exit
						end select options
					end do
					if(index(runoptions,"ry2ev").gt.0) then
						ry2ev = 13.6d0
					else
						ry2ev = 1.d0
					end if
! 					write(*,*) (istring1(n),n=i+1,iomax)
					exit
!===============================================================================
				case("eta")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) eta
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) eta
					end if
! 					write(*,"('eta = ',f8.5)") eta
					nparams = nparams-1
				case("eta=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) eta
! 					write(*,"('eta = ',f8.5)") eta
					nparams = nparams-1
!===============================================================================
				case("lattice")
		    	if(istring1(i+1).eq."=") then ! If after keyword there's an '='
		    		lattice = istring1(i+2)
		    	else ! If there's no '=' after keyword, get from next line
		    		lattice = istring2(i)
		    	end if
! 					write(*,"('lattice = ',a2)") lattice
					nparams = nparams-1
				case("lattice=")
					lattice = istring1(i+1)
! 					write(*,"('lattice = ',a2)") lattice
					nparams = nparams-1
!===============================================================================
				case("Npl")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) Npl
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) Npl
					end if
! 					write(*,"('Npl = ',i0)") Npl
					nparams = nparams-1
				case("Npl=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) Npl
! 					write(*,"('Npl = ',i0)") Npl
					nparams = nparams-1
!===============================================================================
				case("a0")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) a0
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) a0
					end if
! 					write(*,"('a0 = ',f8.5)") a0
					nparams = nparams-1
				case("a0=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) a0
! 					write(*,"('a0 = ',f8.5)") a0
					nparams = nparams-1
!===============================================================================
				case("SOC")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) SOC
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) SOC
					end if
					if(SOC) then
! 						write(*,"('Spin Orbit Coupling: ACTIVATED')")
						nparams2 = nparams2+1
					else
! 						write(*,"('Spin Orbit Coupling: DEACTIVATED')")
					end if
					nparams = nparams-1
				case("SOC=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) SOC
					if(SOC) then
! 						write(*,"('Spin Orbit Coupling: ACTIVATED')")
						nparams2 = nparams2+1
					else
! 						write(*,"('Spin Orbit Coupling: DEACTIVATED')")
					end if
					nparams = nparams-1
!-------------------------------------------------------------------------------
				case("magaxis")
		    	if(istring1(i+1).eq."=") then ! If after keyword there's an '='
		    		magaxis = istring1(i+2)
		    	else ! If there's no '=' after keyword, get from next line
		    		magaxis = istring2(i)
		    	end if
! 					write(*,"('magaxis = ',a2)") magaxis
					nparams2 = nparams2-1
				case("magaxis=")
					magaxis = istring1(i+1)
! 					write(*,"('magaxis = ',a2)") magaxis
					nparams2 = nparams2-1
!===============================================================================
				case("hwx")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) hwx
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) hwx
					end if
! 					write(*,"('hwx = ',f8.5)") hwx
					nparams = nparams-1
				case("hwx=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) hwx
! 					write(*,"('hwx = ',f8.5)") hwx
					nparams = nparams-1
!-------------------------------------------------------------------------------
				case("hwy")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) hwy
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) hwy
					end if
! 					write(*,"('hwy = ',f8.5)") hwy
					nparams = nparams-1
				case("hwy=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) hwy
! 					write(*,"('hwy = ',f8.5)") hwy
					nparams = nparams-1
!-------------------------------------------------------------------------------
				case("hwz")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) hwz
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) hwz
					end if
! 					write(*,"('hwz = ',f8.5)") hwz
					nparams = nparams-1
				case("hwz=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) hwz
! 					write(*,"('hwz = ',f8.5)") hwz
					nparams = nparams-1
!===============================================================================
				case("dirEfield")
		    	if(istring1(i+1).eq."=") then ! If after keyword there's an '='
		    		dirEfield = istring1(i+2)
		    	else ! If there's no '=' after keyword, get from next line
		    		dirEfield = istring2(i)
		    	end if
! 					write(*,"('dirEfield = ',a2)") dirEfield
					nparams = nparams-1
				case("dirEfield=")
					dirEfield = istring1(i+1)
! 					write(*,"('dirEfield = ',a2)") dirEfield
					nparams = nparams-1
!===============================================================================
				case("ncp")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) ncp
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) ncp
					end if
! 					write(*,"('ncp = ',i0,' x ',i0)") ncp
					nparams = nparams-1
				case("ncp=")
						read(unit=istring1(i+1),fmt=*,iostat=ios) ncp
! 					write(*,"('ncp = ',i0,' x ',i0)") ncp
					nparams = nparams-1
!===============================================================================
				case("parts")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) parts
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) parts
					end if
! 					write(*,"('parts = ',i0)") parts
					nparams = nparams-1
				case("parts=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) parts
! 					write(*,"('parts = ',i0)") parts
					nparams = nparams-1
!===============================================================================
				case("n1gl")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) n1gl
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) n1gl
					end if
! 					write(*,"('n1gl = ',i0)") parts
					nparams = nparams-1
				case("n1gl=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) n1gl
! 					write(*,"('n1gl = ',i0)") parts
					nparams = nparams-1
!===============================================================================
				case("parts3")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) parts3
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) parts3
					end if
! 					write(*,"('parts3 = ',i0)") parts3
					nparams = nparams-1
				case("parts3=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) parts3
! 					write(*,"('parts3 = ',i0)") parts3
					nparams = nparams-1
!===============================================================================
				case("n3gl")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) n3gl
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) n3gl
					end if
! 					write(*,"('n3gl = ',i0)") parts
					nparams = nparams-1
				case("n3gl=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) n3gl
! 					write(*,"('n3gl = ',i0)") parts
					nparams = nparams-1
!===============================================================================
				case("emin")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) emin
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) emin
					end if
! 					write(*,"('emin = ',f8.5)") emin
					nparams = nparams-1
				case("emin=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) emin
! 					write(*,"('emin = ',f8.5)") emin
					nparams = nparams-1
!===============================================================================
				case("emax")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) emax
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) emax
					end if
! 					write(*,"('emax = ',f8.5)") emax
					nparams = nparams-1
				case("emax=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) emax
! 					write(*,"('emax = ',f8.5)") emax
					nparams = nparams-1
!===============================================================================
				case("npts")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) npts
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) npts
					end if
! 					write(*,"('npts = ',i0)") npts
					npt1 = npts + 1
					nparams = nparams-1
				case("npts=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) npts
! 					write(*,"('npts = ',i0)") npts
					npt1 = npts + 1
					nparams = nparams-1
!===============================================================================
				case("n0sc1")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) n0sc1
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) n0sc1
					end if
! 					write(*,"('n0sc1 = ',i0)") n0sc1
					nparams = nparams-1
				case("n0sc1=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) n0sc1
! 					write(*,"('n0sc1 = ',i0)") n0sc1
					nparams = nparams-1
!===============================================================================
				case("n0sc2")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) n0sc2
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) n0sc2
					end if
! 					write(*,"('n0sc2 = ',i0)") n0sc2
					nparams = nparams-1
				case("n0sc2=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) n0sc2
! 					write(*,"('n0sc2 = ',i0)") n0sc2
					nparams = nparams-1
!===============================================================================
				case("renorm")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) renorm
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) renorm
					end if
					if(renorm) then
! 						write(*,"('Current renormalization: ACTIVATED')")
						nparams2 = nparams2+1
! 					else
! 						write(*,"('Current renormalization: DEACTIVATED')")
					end if
					nparams = nparams-1
				case("renorm=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) renorm
					if(renorm) then
! 						write(*,"('Current renormalization: ACTIVATED')")
						nparams2 = nparams2+1
					else
! 						write(*,"('Current renormalization: DEACTIVATED')")
					end if
					nparams = nparams-1
!===============================================================================
				case("renormnb")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) renormnb
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) renormnb
					end if
! 					write(*,"('renormnb = ',i0)") renormnb
					nparams2 = nparams2-1
				case("renormnb=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) renormnb
! 					write(*,"('renormnb = ',i0)") renormnb
					nparams2 = nparams2-1
!===============================================================================
				case("kdirection")
		    	if(istring1(i+1).eq."=") then ! If after keyword there's an '='
		    		kdirection = istring1(i+2)
		    	else ! If there's no '=' after keyword, get from next line
		    		kdirection = istring2(i)
		    	end if
! 					write(*,"('kdirection = ',a2)") kdirection
					nparams = nparams-1
				case("kdirection=")
					kdirection = istring1(i+1)
! 					write(*,"('kdirection = ',a2)") kdirection
					nparams = nparams-1
!===============================================================================
				case("skipsc")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) skipsc
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) skipsc
					end if
! 					if(skipsc) then
! 						write(*,"('Skip self-consistency: ACTIVATED')")
! 					else
! 						write(*,"('Skip self-consistency: DEACTIVATED')")
! 					end if
					nparams = nparams-1
				case("skipsc=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) skipsc
! 					if(skipsc) then
! 						write(*,"('Skip self-consistency: ACTIVATED')")
! 					else
! 						write(*,"('Skip self-consistency: DEACTIVATED')")
! 					end if
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
! 					write(*,"('dirEfield = ',a2)") dirEfield
					nparams = nparams-1
				case("dfttype=")
					dfttype = istring1(i+1)
! 					write(*,"('dirEfield = ',a2)") dirEfield
					nparams = nparams-1
!===============================================================================
				case("set1")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) set1
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) set1
					end if
! 					write(*,"('set1 = ',i0)") set1
					nparams = nparams-1
				case("set1=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) set1
! 					write(*,"('set1 = ',i0)") set1
					nparams = nparams-1
!===============================================================================
				case("set2")
					if(istring1(i+1).eq."=") then ! If after keyword there's an '='
						read(unit=istring1(i+2),fmt=*,iostat=ios) set2
					else ! If there's no '=' after keyword, get from next line
						read(unit=istring2(i),fmt=*,iostat=ios) set2
					end if
! 					write(*,"('set2 = ',i0)") set2
					nparams = nparams-1
				case("set2=")
					read(unit=istring1(i+1),fmt=*,iostat=ios) set2
! 					write(*,"('set2 = ',i0)") set2
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
	if((renormnb.lt.n0sc1).or.(renormnb.gt.n0sc2)) then
		if(myrank.eq.0) then
			write(*,"('[ioread] Invalid neighbor for renormalization: ',i0,'!')") renormnb
			write(*,"('[ioread] Choose a value between ',i0,' and ',i0,'.')") n0sc1,n0sc2
		end if
		stop
	end if
!-------------------------------------------------------------------------------
!*********** User manual additions / modifications in the input file **********!
! 	Npl = 4
! 	ncp = 6
! 	SOC = .true.
! 	magaxis = "5"
! 	runoptions = trim(runoptions) // " noUonNM"
! 	scfile = "results/selfconsistency/selfconsistency_Npl=4_dfttype=T_parts=2_U= 0.7E-01_hwx= 0.0E+00_hwy= 0.0E+00_hwz= 0.0E+00_ncp=6_eta= 0.5E-03.dat"
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
	end subroutine ioread
end module mod_parameters