! For a given value of center of band eps1 it calculates the
! occupation number and the magnetic moment
subroutine nparticles(N,eps,npart,iflag)
	use mod_constants
	use mod_parameters
	use mod_f90_kind
	use mod_generate_epoints
	use mod_magnet
	use mod_tight_binding, only: npart0
	use mod_progress
	use mod_mpi_pars
	use MPI
	implicit none
	integer         :: N,i,j,iflag
	real(double)    :: eps(N),npart(N),n_t(Npl)
	real(double),dimension(Npl,9)       :: n_orb_u,n_orb_d,n_orb_t,mag_orb
	real(double),dimension(Npl,9)       :: gdiaguur,gdiagddr
	complex(double),dimension(Npl,9)    :: gdiagud,gdiagdu
	!--------------------- begin MPI vars --------------------
	integer :: ix,itask
	integer :: ncount
	ncount=Npl*9
	!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

#ifndef _JUQUEEN
	open(6,carriagecontrol ='fortran')
#endif

	eps1 = eps

	n_orb_u = 0.d0
	n_orb_d = 0.d0

	ix = myrank+1
	itask = numprocs ! Number of tasks done initially

	! Calculating the number of particles for each spin and orbital using a complex integral
	if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
		call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
		gdiaguur = wght(ix)*gdiaguur
		gdiagddr = wght(ix)*gdiagddr
		gdiagud = wght(ix)*gdiagud
		gdiagdu = wght(ix)*gdiagdu

		n_orb_u = gdiaguur
		n_orb_d = gdiagddr

		do j=1,Npl
			splus(j) = (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
		end do

 		if(index(runoptions,"verbose").gt.0) write(*,"('[nparticles] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
		do i=2,pn1
			! Progress bar
      prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
      energy_progress_bar: select case (mod(i,4))
      case(0)
        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '|',prog,char(13)
      case(1)
        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '/',prog,char(13)
      case(2)
        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '-',prog,char(13)
      case(3)
        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '\',prog,char(13)
      end select energy_progress_bar
#else
			elapsed_time = MPI_Wtime() - start_time
			write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),100*(i+1)/pn1
#endif

			call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999,MPI_COMM_WORLD,stat,ierr)
			call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),9998,MPI_COMM_WORLD,stat,ierr)
			call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9997,MPI_COMM_WORLD,stat,ierr)
			call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9996,MPI_COMM_WORLD,stat,ierr)

			n_orb_u = n_orb_u + gdiaguur
			n_orb_d = n_orb_d + gdiagddr

			do j=1,Npl
				splus(j) = splus(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
			end do

			! If the number of processors is less than the total number of points, sends
			! the rest of the points to the ones that finish first
			if (itask.lt.pn1) then
				itask = itask + 1
				call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)
			else
				call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
			end if
		end do
	else
		! Other processors calculate each point of the integral and waits for new points
		do
			if(ix.gt.pn1) exit

			! First and second integrations (in the complex plane)
			call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
			gdiaguur = wght(ix)*gdiaguur
			gdiagddr = wght(ix)*gdiagddr
			gdiagud = wght(ix)*gdiagud
			gdiagdu = wght(ix)*gdiagdu

 			if(index(runoptions,"verbose").gt.0) write(*,"('[nparticles] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
			! Sending results to process 0
			call MPI_Send(gdiaguur,ncount,MPI_DOUBLE_PRECISION,0,9999,MPI_COMM_WORLD,ierr)
			call MPI_Send(gdiagddr,ncount,MPI_DOUBLE_PRECISION,0,9998,MPI_COMM_WORLD,ierr)
			call MPI_Send(gdiagud,ncount,MPI_DOUBLE_COMPLEX,0,9997,MPI_COMM_WORLD,ierr)
			call MPI_Send(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,0,9996,MPI_COMM_WORLD,ierr)
			! Receiving new point or signal to exit
			call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
			if(ix.eq.0) exit
		end do
	end if

	! Send results to all processors
	call MPI_Bcast(n_orb_u,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(n_orb_d,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_Bcast(splus,Npl,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

	n_orb_u = 0.5d0 + n_orb_u/pi
	n_orb_d = 0.5d0 + n_orb_d/pi
	n_orb_t = n_orb_u + n_orb_d
	mag_orb = n_orb_u - n_orb_d
	splus  = splus/tpi
	sminus = conjg(splus)

	!	Number of particles and magnetization:
	do i=1,Npl
		n_t(i)    = sum(n_orb_t(i,:))
		mag(i)    = sum(mag_orb(i,5:9))
		npart(i)  = n_t(i) - npart0(i+1)
	end do

	if(myrank.eq.0) then
		if(index(runoptions,"verbose").gt.0)  then
			do i=1,Npl
				write(*,"('plane ',I0,', eps1(',I0,')=',e16.9,4x,'mag(',I0,')=',e16.9)") i,i,eps1(i),i,mag(i)
				write(*,"('n_t = ',e16.9,', npart0 = ',e16.9,', npart = ',e16.9)") n_t(i),npart0(i+1),npart(i)
				if(abs(splus(i)).gt.1.d-10) then
					write(*,"(30x,'Sx(',I0,')=',e16.9,4x,'Sy(',I0,')=',e16.9)") i,real(splus(i)),i,aimag(splus(i))
				end if
			end do
			write(*,"(24x,'----')")
		end if
	end if

	return
end subroutine nparticles

subroutine nparticlesjac(N,eps,npart,npartjac,ldfjac,iflag)
	use mod_constants
	use mod_parameters
	use mod_f90_kind
	use mod_generate_epoints
	use mod_magnet
	use mod_progress
	use mod_mpi_pars
	use MPI
	implicit none
	integer       :: N,ldfjac,i,j,iflag
	real(double)  :: eps(N),npart(N),npartjac(ldfjac,N)
	real(double)  :: ggr(Npl,Npl)
	!--------------------- begin MPI vars --------------------
	integer :: ix,itask
	integer :: ncount
	!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
	ncount=Npl*Npl

#ifndef _JUQUEEN
	open(6,carriagecontrol ='fortran')
#endif

 	eps1 = eps

	ix = myrank+1
	itask = numprocs ! Number of tasks done initially

	npartjac = 0.d0
	! Calculating the jacobian
	if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
		call sumk_npartjac(Ef,y(ix),ggr)
		npartjac = wght(ix)*ggr

 		if(index(runoptions,"verbose").gt.0) write(*,"('[nparticlesjac] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
		do i=2,pn1
			! Progress bar
      prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
      energy_progress_bar: select case (mod(i,4))
      case(0)
        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '|',prog,char(13)
      case(1)
        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '/',prog,char(13)
      case(2)
        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '-',prog,char(13)
      case(3)
        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '\',prog,char(13)
      end select energy_progress_bar
#else
			elapsed_time = MPI_Wtime() - start_time
			write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),100*(i+1)/pn1
#endif

			call MPI_Recv(ggr,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3333,MPI_COMM_WORLD,stat,ierr)

			npartjac = npartjac + ggr

			! If the number of processors is less than the total number of points, sends
			! the rest of the points to the ones that finish first
			if (itask.lt.pn1) then
				itask = itask + 1
				call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)
			else
				call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
			end if
		end do
	else
		! Other processors calculate each point of the integral and waits for new points
		do
			if(ix.gt.pn1) exit

			! First and second integrations (in the complex plane)
			call sumk_npartjac(Ef,y(ix),ggr)
			ggr = wght(ix)*ggr

 			if(index(runoptions,"verbose").gt.0) write(*,"('[nparticlesjac] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
			! Sending results to process 0
			call MPI_Send(ggr,ncount,MPI_DOUBLE_PRECISION,0,3333,MPI_COMM_WORLD,ierr)
			! Receiving new point or signal to exit
			call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
			if(ix.eq.0) exit
		end do
	end if

	! Send results to all processors
	call MPI_Bcast(npartjac,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	npartjac = npartjac/pi

	return
end subroutine nparticlesjac

! For a given value of center of band eps1 it calculates the
! occupation number and the magnetic moment
subroutine nparticlesjacnag(N,eps,npart,npartjac,ldfjac,iflag)
	use mod_constants
	use mod_parameters
	use mod_f90_kind
	use mod_generate_epoints
	use mod_magnet
	use mod_tight_binding, only: npart0
	use mod_progress
	use mod_mpi_pars
	use MPI
	implicit none
	integer         :: N,ldfjac,i,j,iflag
	real(double)    :: eps(N),npart(N),n_t(Npl),npartjac(ldfjac,N)
	real(double)    :: ggr(Npl,Npl)
	real(double),dimension(Npl,9)       :: n_orb_u,n_orb_d,n_orb_t,mag_orb
	real(double),dimension(Npl,9)       :: gdiaguur,gdiagddr
	complex(double),dimension(Npl,9)    :: gdiagud,gdiagdu
	!--------------------- begin MPI vars --------------------
	integer :: ix,itask
	integer :: ncount,ncount2
	!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
	ncount=Npl*9
	ncount2=Npl*Npl

#ifndef _JUQUEEN
	open(6,carriagecontrol ='fortran')
#endif

	eps1 = eps

	ix = myrank+1
	itask = numprocs ! Number of tasks done initially

	flag: select case (iflag)
	case(1)
		n_orb_u = 0.d0
		n_orb_d = 0.d0
		! Calculating the number of particles for each spin and orbital using a complex integral
		if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
			call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
			gdiaguur = wght(ix)*gdiaguur
			gdiagddr = wght(ix)*gdiagddr
			gdiagud = wght(ix)*gdiagud
			gdiagdu = wght(ix)*gdiagdu

			n_orb_u = gdiaguur
			n_orb_d = gdiagddr

			do j=1,Npl
				splus(j) = (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
			end do

	 		if(index(runoptions,"verbose").gt.0) write(*,"('[nparticlesjacnag] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
			do i=2,pn1
				! Progress bar
	      prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
	      npart_energy_progress_bar: select case (mod(i,4))
	      case(0)
	        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '|',prog,char(13)
	      case(1)
	        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '/',prog,char(13)
	      case(2)
	        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '-',prog,char(13)
	      case(3)
	        write(*,"(a1,2x,i3,'% of nparticles energy sum',a1,$)") '\',prog,char(13)
	      end select npart_energy_progress_bar
#else
				elapsed_time = MPI_Wtime() - start_time
				write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
	      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),100*(i+1)/pn1
#endif

				call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999,MPI_COMM_WORLD,stat,ierr)
				call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),9998,MPI_COMM_WORLD,stat,ierr)
				call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9997,MPI_COMM_WORLD,stat,ierr)
				call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9996,MPI_COMM_WORLD,stat,ierr)

				n_orb_u = n_orb_u + gdiaguur
				n_orb_d = n_orb_d + gdiagddr

				do j=1,Npl
					splus(j) = splus(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
				end do

				! If the number of processors is less than the total number of points, sends
				! the rest of the points to the ones that finish first
				if (itask.lt.pn1) then
					itask = itask + 1
					call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)
				else
					call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
				end if
			end do
		else
			! Other processors calculate each point of the integral and waits for new points
			do
				if(ix.gt.pn1) exit

				! First and second integrations (in the complex plane)
				call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
				gdiaguur = wght(ix)*gdiaguur
				gdiagddr = wght(ix)*gdiagddr
				gdiagud = wght(ix)*gdiagud
				gdiagdu = wght(ix)*gdiagdu

	 			if(index(runoptions,"verbose").gt.0) write(*,"('[nparticlesjacnag] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
				! Sending results to process 0
				call MPI_Send(gdiaguur,ncount,MPI_DOUBLE_PRECISION,0,9999,MPI_COMM_WORLD,ierr)
				call MPI_Send(gdiagddr,ncount,MPI_DOUBLE_PRECISION,0,9998,MPI_COMM_WORLD,ierr)
				call MPI_Send(gdiagud,ncount,MPI_DOUBLE_COMPLEX,0,9997,MPI_COMM_WORLD,ierr)
				call MPI_Send(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,0,9996,MPI_COMM_WORLD,ierr)
				! Receiving new point or signal to exit
				call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
				if(ix.eq.0) exit
			end do
		end if

		! Send results to all processors
		call MPI_Bcast(n_orb_u,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(n_orb_d,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		call MPI_Bcast(splus,Npl,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

		n_orb_u = 0.5d0 + n_orb_u/pi
		n_orb_d = 0.5d0 + n_orb_d/pi
		n_orb_t = n_orb_u + n_orb_d
		mag_orb = n_orb_u - n_orb_d
		splus  = splus/tpi
		sminus = conjg(splus)

		!	Number of particles and magnetization:
		do i=1,Npl
			n_t(i)    = sum(n_orb_t(i,:))
			mag(i)    = sum(mag_orb(i,5:9))
			npart(i)  = n_t(i) - npart0(i+1)
		end do

		if(myrank.eq.0) then
			if(index(runoptions,"verbose").gt.0)  then
				do i=1,Npl
					write(*,"('plane ',I0,', eps1(',I0,')=',e16.9,4x,'mag(',I0,')=',e16.9)") i,i,eps1(i),i,mag(i)
					write(*,"('n_t = ',e16.9,', npart0 = ',e16.9,', npart = ',e16.9)") n_t(i),npart0(i+1),npart(i)
					if(abs(splus(i)).gt.1.d-10) then
						write(*,"(30x,'Sx(',I0,')=',e16.9,4x,'Sy(',I0,')=',e16.9)") i,real(splus(i)),i,aimag(splus(i))
					end if
				end do
				write(*,"(24x,'----')")
			end if
		end if
	case(2)
		npartjac = 0.d0
		! Calculating the jacobian
		if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
			call sumk_npartjac(Ef,y(ix),ggr)
			npartjac = wght(ix)*ggr

	 		if(index(runoptions,"verbose").gt.0) write(*,"('[nparticlesjacnag] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
			do i=2,pn1
				! Progress bar
	      prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
	      jac_energy_progress_bar: select case (mod(i,4))
	      case(0)
	        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '|',prog,char(13)
	      case(1)
	        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '/',prog,char(13)
	      case(2)
	        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '-',prog,char(13)
	      case(3)
	        write(*,"(a1,2x,i3,'% of jacobian energy sum',a1,$)") '\',prog,char(13)
	      end select jac_energy_progress_bar
#else
				elapsed_time = MPI_Wtime() - start_time
				write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
	      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),100*(i+1)/pn1
#endif

				call MPI_Recv(ggr,ncount2,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3333,MPI_COMM_WORLD,stat,ierr)

				npartjac = npartjac + ggr

				! If the number of processors is less than the total number of points, sends
				! the rest of the points to the ones that finish first
				if (itask.lt.pn1) then
					itask = itask + 1
					call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_COMM_WORLD,ierr)
				else
					call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_COMM_WORLD,ierr)
				end if
			end do
		else
			! Other processors calculate each point of the integral and waits for new points
			do
				if(ix.gt.pn1) exit

				! First and second integrations (in the complex plane)
				call sumk_npartjac(Ef,y(ix),ggr)
				ggr = wght(ix)*ggr

	 			if(index(runoptions,"verbose").gt.0) write(*,"('[nparticlesjacnag] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
				! Sending results to process 0
				call MPI_Send(ggr,ncount2,MPI_DOUBLE_PRECISION,0,3333,MPI_COMM_WORLD,ierr)
				! Receiving new point or signal to exit
				call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
				if(ix.eq.0) exit
			end do
		end if

		! Send results to all processors
		call MPI_Bcast(npartjac,ncount2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

		npartjac = npartjac/pi
	case default
		write(*,"('[nparticles] Problem in self-consistency! iflag = ',I0)") iflag
		stop
	end select flag

	return
end subroutine nparticlesjacnag