subroutine openclose_chi_files(iflag)
	use mod_parameters
	implicit none

	character(len=500)	:: varm
	character(len=50)   :: fieldpart,socpart
	character(len=2)	  :: spin(4)
	character(len=1)	  :: SOCc
	integer :: i,j,sigma,iw,iflag,err,errt=0

	spin(1) = "pm"
	spin(2) = "um"
	spin(3) = "dm"
	spin(4) = "mm"

  fieldpart = ""
  socpart   = ""
	if(SOC) then
		if(llinearsoc) then
			SOCc = "L"
		else
			SOCc = "T"
		end if
		write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
	else
		SOCc = "F"
	end if
  if(FIELD) then
    write(fieldpart,"('_hwx=',e8.1,'_hwy=',e8.1,'_hwz=',e8.1)") hwx,hwy,hwz
    if(index(runoptions,"tesla").gt.0) fieldpart = trim(fieldpart) // "_tesla"
  end if

	if(iflag.eq.0) then
		do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
			! RPA SUSCEPTIBILITIES
			iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/RPA/',a2,'/chi_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  real part of chi ',a,'  ,  imaginary part of chi ',a,'  ,  amplitude of chi ',a,'  ')") spin(sigma),spin(sigma),spin(sigma)
			close(unit=iw)
			! HF SUSCEPTIBILITIES
			iw = iw+1000
 			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/HF/',a2,'/chihf_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
 			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  real part of chi ',a,' HF ,  imaginary part of chi ',a,' HF  ,  amplitude of chi ',a,' HF  ')") spin(sigma),spin(sigma),spin(sigma)
			close(unit=iw)
		end do ; end do ; end do
		! RPA DIAGONALIZATION
		if(nmaglayers.gt.1) then
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/RPA/pm/chi_eval_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
			open (unit=1990, file=varm,status='unknown', form='formatted')
			write(unit=1990,fmt="('#   energy  ,  real part of 1st eigenvalue  ,  imaginary part of 1st eigenvalue  ,  real part of 2nd eigenvalue  ,  imaginary part of 2nd eigenvalue  , ... ')")
			close (unit=1990)
			do i=1,nmaglayers
				write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/RPA/pm/chi_evec',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
				open (unit=1990+i, file=varm,status='unknown', form='formatted')
				write(unit=1990+i,fmt="('#   energy  ,  real part of 1st component  ,  imaginary part of 1st component  ,  real part of 2nd component  ,  imaginary part of 2nd component  , ...   ')")
				close (unit=1990+i)
			end do
		end if

	else if(iflag.eq.1) then
		do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
			! RPA SUSCEPTIBILITIES
			iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/RPA/',a2,'/chi_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
			open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
			errt = errt + err
 			! HF SUSCEPTIBILITIES
			iw = iw+1000
 			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/HF/',a2,'/chihf_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
 			open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
			errt = errt + err
		end do ; end do ; end do
		! RPA DIAGONALIZATION
		if(nmaglayers.gt.1) then
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/RPA/pm/chi_eval_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
			open (unit=1990, file=varm, status='old', position='append', form='formatted', iostat=err)
			errt = errt + err
			do i=1,nmaglayers
				write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/RPA/pm/chi_evec',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
				open (unit=1990+i, file=varm, status='old', position='append', form='formatted', iostat=err)
				errt = errt + err
			end do
		end if
		! Stop if some file does not exist
    if(errt.ne.0) then
      write(*,"(a,i0,a)") "[openclose_chi_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
      stop
    end if

	else
		do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
			! RPA SUSCEPTIBILITIES
			iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
			close(unit=iw)
			! HF SUSCEPTIBILITIES
			iw = iw+1000
			close(unit=iw)
		end do ; end do ; end do
		! RPA DIAGONALIZATION
		if(nmaglayers.gt.1) then
			close (unit=1990)
			do i=1,nmaglayers
				close (unit=1990+i)
			end do
		end if
	end if

	return
end subroutine openclose_chi_files

subroutine openclose_sd_files(iflag)
	use mod_parameters
	implicit none

	character(len=500)	:: varm
	character(len=50)   :: fieldpart,socpart
	character(len=2)		:: folder(7),filename(7)
	character(len=1)	  :: SOCc
	integer :: i,j,iw,iflag,err,errt=0

  fieldpart = ""
  socpart   = ""
	if(SOC) then
		if(llinearsoc) then
			SOCc = "L"
		else
			SOCc = "T"
		end if
		write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
	else
		SOCc = "F"
	end if
  if(FIELD) then
    write(fieldpart,"('_hwx=',e8.1,'_hwy=',e8.1,'_hwz=',e8.1)") hwx,hwy,hwz
    if(index(runoptions,"tesla").gt.0) fieldpart = trim(fieldpart) // "_tesla"
  end if

	folder(1) = "CD"
	folder(2) = "SD"
	folder(3) = "SD"
	folder(4) = "SD"
	folder(5) = "LD"
	folder(6) = "LD"
	folder(7) = "LD"

	filename(1) = "Cd"
	filename(2) = "Sx"
	filename(3) = "Sy"
	filename(4) = "Sz"
	filename(5) = "Lx"
	filename(6) = "Ly"
	filename(7) = "Lz"

	if(iflag.eq.0) then
		do i=1,Npl ; do j=1,7
			iw = 3000+(i-1)*7+j
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
			close(unit=iw)
			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/r',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOCc,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
				open (unit=iw, file=varm, status='unknown', form='formatted')
				write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
				close(unit=iw)
			end if
		end do ; end do
	else if (iflag.eq.1) then
		do i=1,Npl ; do j=1,7
			iw = 3000+(i-1)*7+j
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
			open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
			errt = errt + err
			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/r',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOCc,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
				open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
				errt = errt + err
			end if
		end do ; end do
		! Stop if some file does not exist
    if(errt.ne.0) then
      write(*,"(a,i0,a)") "[openclose_sd_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
      stop
    end if
	else
		do i=1,Npl ; do j=1,7
			iw = 3000+(i-1)*7+j
			close(unit=iw)

			if(renorm) then
				iw = iw+1000
				close(unit=iw)
			end if
		end do ; end do
	end if

	return
end subroutine openclose_sd_files

subroutine openclose_sc_files(iflag)
	use mod_parameters
	implicit none

	character(len=500)	:: varm
	character(len=50)   :: fieldpart,socpart
	character(len=1)	  :: SOCc
	character(len=2)	  :: folder(7)
	character(len=3)	  :: filename(7)
	integer :: i,j,iw,neighbor,iflag,err,errt=0

  fieldpart = ""
  socpart   = ""
	if(SOC) then
		if(llinearsoc) then
			SOCc = "L"
		else
			SOCc = "T"
		end if
		write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
	else
		SOCc = "F"
	end if
  if(FIELD) then
    write(fieldpart,"('_hwx=',e8.1,'_hwy=',e8.1,'_hwz=',e8.1)") hwx,hwy,hwz
    if(index(runoptions,"tesla").gt.0) fieldpart = trim(fieldpart) // "_tesla"
  end if

	folder(1) = "CC"
	folder(2) = "SC"
	folder(3) = "SC"
	folder(4) = "SC"
	folder(5) = "LC"
	folder(6) = "LC"
	folder(7) = "LC"

	filename(1) = "Ich"
	filename(2) = "Isx"
	filename(3) = "Isy"
	filename(4) = "Isz"
	filename(5) = "Ilx"
	filename(6) = "Ily"
	filename(7) = "Ilz"

	if(iflag.eq.0) then
		 do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
			iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/prll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
			close(unit=iw)

			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/rprll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOCc,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
				open (unit=iw, file=varm, status='unknown', form='formatted')
				write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
				close(unit=iw)
			end if
		end do ; end do ; end do
	else if(iflag.eq.1) then
		 do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
			iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/prll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
			open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
			errt = errt + err
			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/rprll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOCc,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
				open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
				errt = errt + err
			end if
		end do ; end do ; end do
		! Stop if some file does not exist
    if(errt.ne.0) then
      write(*,"(a,i0,a)") "[openclose_sc_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
      stop
    end if
	else
		 do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
			iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
			close(unit=iw)

			if(renorm) then
				iw = iw+1000
				close(unit=iw)
			end if
		end do ; end do ; end do
	end if

	return
end subroutine openclose_sc_files


subroutine openclose_beff_files(iflag)
	use mod_parameters
	implicit none

	character(len=500)	:: varm
	character(len=50)   :: fieldpart,socpart
	character(len=1)	  :: direction(4),SOCc
	integer :: i,sigma,iw,iflag,err,errt=0

  fieldpart = ""
  socpart   = ""
	if(SOC) then
		if(llinearsoc) then
			SOCc = "L"
		else
			SOCc = "T"
		end if
		write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
	else
		SOCc = "F"
	end if
  if(FIELD) then
    write(fieldpart,"('_hwx=',e8.1,'_hwy=',e8.1,'_hwz=',e8.1)") hwx,hwy,hwz
    if(index(runoptions,"tesla").gt.0) fieldpart = trim(fieldpart) // "_tesla"
  end if

	direction(1) = "0"
	direction(2) = "x"
	direction(3) = "y"
	direction(4) = "z"

	if(iflag.eq.0) then
		do sigma=1,4 ; do i=1,Npl
			iw = 7000+(sigma-1)*Npl+i
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/Beff/Beff',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  , amplitude of Beff_',a,' , real part of Beff_',a,' , imaginary part of Beff_',a,' , phase of Beff_',a,' , cosine of Beff_',a,'  ,  sine of ',a,'  ')") direction(sigma),direction(sigma),direction(sigma),direction(sigma),direction(sigma),direction(sigma)
			close(unit=iw)
		end do ; end do
	else if (iflag.eq.1) then
		do sigma=1,4 ; do i=1,Npl
			iw = 7000+(sigma-1)*Npl+i
			write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/Beff/Beff',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
			open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
			errt = errt + err
		end do ; end do
		! Stop if some file does not exist
    if(errt.ne.0) then
      write(*,"(a,i0,a)") "[openclose_beff_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
      stop
    end if
	else
		do sigma=1,4 ; do i=1,Npl
			iw = 7000+(sigma-1)*Npl+i
			close(unit=iw)
		end do ; end do
	end if

	return
end subroutine openclose_beff_files
