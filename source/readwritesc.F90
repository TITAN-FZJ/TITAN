! This subroutine reads previous band-shifting results and writes new ones in a file
subroutine readwritesc(iflag,err)
	use mod_constants
	use mod_parameters
	use mod_magnet
	use mod_mpi_pars, only: myrank
	implicit none
	character(len=200)  :: file
	integer 						:: i,err,iflag
	real(double),dimension(Npl)		:: Sx,Sy

!   Reading previous results (Sx, Sy, mag and eps1) from files (if available)
	if(iflag.eq.0) then
		if(trim(scfile).eq."") then ! If a filename is not given in inputcard, use the default one
			if(SOC) then
				write(file,"('./results/selfconsistency/selfconsistency_Npl=',I0,'_dfttype=',A,'_parts=',I0,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_ncp=',I0,'_eta=',E8.1,'_magaxis=',A,'_socscale=',f5.2,'.dat')") Npl,dfttype,parts,Utype,hwx,hwy,hwz,ncp,eta,magaxis,socscale
			else
				write(file,"('./results/selfconsistency/selfconsistency_Npl=',I0,'_dfttype=',A,'_parts=',I0,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_ncp=',I0,'_eta=',E8.1,'.dat')") Npl,dfttype,parts,Utype,hwx,hwy,hwz,ncp,eta
			end if
			open(unit=99,file=file,status="old",iostat=err)
		else
			open(unit=99,file=scfile,status="old",iostat=err)
		end if
		if(err.eq.0) then
			if(myrank.eq.0) write(*,"('[readwritesc] Self-consistency file already exists. Reading it now...')")
			do i=1,Npl
				read(99,"(e21.11)") eps1(i)
				read(99,"(e21.11)") Sx(i)
				read(99,"(e21.11)") Sy(i)
				read(99,"(e21.11)") mag(i)
			end do
			splus	 = Sx + zi*Sy
	    do i=1,Npl
	      hdel(i)    = 0.5d0*U(i+1)*mag(i)
	      usplus(i)  = U(i+1)*splus(i)
	    end do
      usminus = conjg(usplus)
			iflag   = 1 ! something was read
		else
			! If file does not exist, try to read for parts-1
			close(99)
			if(SOC) then
				write(file,"('./results/selfconsistency/selfconsistency_Npl=',I0,'_dfttype=',A,'_parts=',I0,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_ncp=',I0,'_eta=',E8.1,'_magaxis=',A,'_socscale=',f5.2,'.dat')") Npl,dfttype,parts-1,Utype,hwx,hwy,hwz,ncp,eta,magaxis,socscale
			else
				write(file,"('./results/selfconsistency/selfconsistency_Npl=',I0,'_dfttype=',A,'_parts=',I0,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_ncp=',I0,'_eta=',E8.1,'.dat')") Npl,dfttype,parts-1,Utype,hwx,hwy,hwz,ncp,eta
			end if
			open(unit=99,file=file,status="old",iostat=err)
			if(err.eq.0) then
				if(myrank.eq.0) then
					write(*,"('[readwritesc] Self-consistency file does not exist. Reading results for parts-1 now...')")
					write(*,"('[readwritesc] Updating values obtained for parts-1...')")
				end if
				do i=1,Npl
					read(99,"(e21.11)") eps1(i)
					read(99,"(e21.11)") Sx(i)
					read(99,"(e21.11)") Sy(i)
					read(99,"(e21.11)") mag(i)
				end do
				splus	 = Sx + zi*Sy
		    do i=1,Npl
		      hdel(i)    = 0.5d0*U(i+1)*mag(i)
		      usplus(i)  = U(i+1)*splus(i)
		    end do
	      usminus = conjg(usplus)
				iflag = 1 ! something was read
				err = 1   ! but not the same parameters
			end if
		end if
		close(99)
!   Writing new results (Sx, Sy, mag and eps1) and mag to file
	else
		if(myrank.eq.0) write(*,"('[readwritesc] Writing new eps1, Sx, Sy and mag to file...')")
		if(SOC) then
			write(file,"('./results/selfconsistency/selfconsistency_Npl=',I0,'_dfttype=',A,'_parts=',I0,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_ncp=',I0,'_eta=',E8.1,'_magaxis=',A,'_socscale=',f5.2,'.dat')") Npl,dfttype,parts,Utype,hwx,hwy,hwz,ncp,eta,magaxis,socscale
		else
			write(file,"('./results/selfconsistency/selfconsistency_Npl=',I0,'_dfttype=',A,'_parts=',I0,'_Utype=',i0,'_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1,'_ncp=',I0,'_eta=',E8.1,'.dat')") Npl,dfttype,parts,Utype,hwx,hwy,hwz,ncp,eta
		end if
		open (unit=99,status="unknown",file=file)
		do i=1,Npl
			write(99,"(e21.11,2x,'! eps1')") eps1(i)
			write(99,"(e21.11,2x,'! Sx')") real(splus(i))
			write(99,"(e21.11,2x,'! Sy')") aimag(splus(i))
			write(99,"(e21.11,2x,'! mag')") mag(i)
		end do
		close(99)
	end if
	return
end subroutine readwritesc