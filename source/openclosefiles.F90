subroutine openclosechifiles(iflag)
	use mod_parameters
	use mod_magnet
	use mod_lattice
	use mod_generate_epoints
	use mod_tight_binding, only: nmaglayers
	implicit none

	character(len=500)	:: varm
	character(len=2)	  :: spin(4)
	integer :: i,j,sigma,iw,iflag

	spin(1) = "pm"
	spin(2) = "um"
	spin(3) = "dm"
	spin(4) = "mm"

	if(iflag.eq.0) then
		do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
			! RPA SUSCEPTIBILITIES
			iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/RPA/',a2,'/chi_',i0,'_',i0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  real part of chi ',a,'  ,  imaginary part of chi ',a,'  ,  amplitude of chi ',a,'  ')") spin(sigma),spin(sigma),spin(sigma)
			close(unit=iw)
			! HF SUSCEPTIBILITIES
			iw = iw+1000
 			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/HF/',a2,'/chihf_',i0,'_',i0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
 			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  real part of chi ',a,' HF ,  imaginary part of chi ',a,' HF  ,  amplitude of chi ',a,' HF  ')") spin(sigma),spin(sigma),spin(sigma)
			close(unit=iw)
		end do ; end do ; end do
		! RPA DIAGONALIZATION
		if(nmaglayers.gt.1) then
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/RPA/pm/chi_eval_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
			open (unit=1990, file=varm,status='unknown', form='formatted')
			write(unit=1990,fmt="('#   energy  ,  real part of 1st eigenvalue  ,  imaginary part of 1st eigenvalue  ,  real part of 2nd eigenvalue  ,  imaginary part of 2nd eigenvalue  , ... ')")
			close (unit=1990)
			do i=1,nmaglayers
				write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/RPA/pm/chi_evec',i0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
				open (unit=1990+i, file=varm,status='unknown', form='formatted')
				write(unit=1990+i,fmt="('#   energy  ,  real part of 1st component  ,  imaginary part of 1st component  ,  real part of 2nd component  ,  imaginary part of 2nd component  , ...   ')")
				close (unit=1990+i)
			end do
		end if

	else if(iflag.eq.1) then
		do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
			! RPA SUSCEPTIBILITIES
			iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/RPA/',a2,'/chi_',i0,'_',i0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
			open (unit=iw, file=varm, status='old', position='append', form='formatted')
 			! HF SUSCEPTIBILITIES
			iw = iw+1000
 			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/HF/',a2,'/chihf_',i0,'_',i0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
 			open (unit=iw, file=varm, status='old', position='append', form='formatted')
		end do ; end do ; end do
		! RPA DIAGONALIZATION
		if(nmaglayers.gt.1) then
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/RPA/pm/chi_eval_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
			open (unit=1990, file=varm, status='old', position='append', form='formatted')
			do i=1,nmaglayers
				write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/RPA/pm/chi_evec',i0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_qx=',E10.3,'_qz=',E10.3,'.dat')") SOC,Npl,i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,q(1),q(2)
				open (unit=1990+i, file=varm, status='old', position='append', form='formatted')
			end do
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
end subroutine openclosechifiles

subroutine openclosesdfiles(iflag)
	use mod_parameters
	use mod_magnet
	use mod_lattice
	use mod_generate_epoints
	implicit none

	character(len=500)	:: varm
	character(len=2)		:: folder(7),filename(7)
	integer :: i,j,iw,iflag


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
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/',a,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'.dat')") SOC,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield
			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
			close(unit=iw)
			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/r',a,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOC,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield,renormnb
				open (unit=iw, file=varm, status='unknown', form='formatted')
				write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
				close(unit=iw)
			end if
		end do ; end do
	else if (iflag.eq.1) then
		do i=1,Npl ; do j=1,7
			iw = 3000+(i-1)*7+j
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/',a,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'.dat')") SOC,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield
			open (unit=iw, file=varm, status='old', position='append', form='formatted')
			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/r',a,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOC,Npl,folder(j),filename(j),i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield,renormnb
				open (unit=iw, file=varm, status='old', position='append', form='formatted')
			end if
		end do ; end do
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
end subroutine openclosesdfiles

subroutine openclosescfiles(iflag)
	use mod_parameters
	use mod_magnet
	use mod_generate_epoints
	implicit none

	character(len=500)	:: varm
	character(len=2)	:: folder(7)
	character(len=3)	:: filename(7)
	integer :: i,j,iw,neighbor,iflag

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
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/prll',a,'_neighbor=',I0,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'.dat')") SOC,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield
			open (unit=iw, file=varm, status='unknown', form='formatted')
			write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
			close(unit=iw)

			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/rprll',a,'_neighbor=',I0,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOC,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield,renormnb
				open (unit=iw, file=varm, status='unknown', form='formatted')
				write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
				close(unit=iw)
			end if
		end do ; end do ; end do
	else if(iflag.eq.1) then
		 do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
			iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
			write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/prll',a,'_neighbor=',I0,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'.dat')") SOC,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield
			open (unit=iw, file=varm, status='old', position='append', form='formatted')

			if(renorm) then
				iw = iw+1000
				write(varm,"('./results/SOC=',L1,'/Npl=',I0,'/',a,'/rprll',a,'_neighbor=',I0,'_pos=',I0,'_parts=',I0,'_parts3=',I0,'_ncp=',I0,'_eta=',E8.1,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_magaxis=',A,'_dirEfield=',A,'_renormnb=',i0,'.dat')") SOC,Npl,folder(j),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,hwx,hwy,hwz,magaxis,dirEfield,renormnb
				open (unit=iw, file=varm, status='old', position='append', form='formatted')
			end if
		end do ; end do ; end do
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
end subroutine openclosescfiles