! Calculates the orbital angular momentum ground state
subroutine L_gs()
	use mod_f90_kind
	use mod_constants, only: zero
	use mod_parameters
	use mod_magnet
	use mod_generate_epoints
	use mod_mpi_pars
	use MPI
	implicit none
	integer			:: i,mu,nu
	complex(double), dimension(9,9)		:: lx,ly,lz
	complex(double),dimension(Npl,9,9)	:: gupgd,gupgdint
!--------------------- begin MPI vars --------------------
	integer :: ix,itask
	integer :: ncount
	ncount=Npl*9*9
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

	gupgdint	= zero

	ix = myrank+1
	itask = numprocs ! Number of tasks done initially

!   Calculating the number of particles for each spin and orbital using a complex integral
	if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
		call sumk_L_gs(Ef,y(ix),gupgd)
		gupgdint = gupgd*wght(ix)
		if(index(runoptions,"verbose").gt.0) write(*,"('[L_gs] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)

		do i=2,pn1
			call MPI_Recv(gupgd,ncount,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,8999,MPI_COMM_WORLD,stat,ierr)
			write(*,"('[L_gs] Point ',i0,' received from ',i0)") i,stat(MPI_SOURCE)

			gupgdint = gupgdint + gupgd

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
			call sumk_L_gs(Ef,y(ix),gupgd)
			gupgd = gupgd*wght(ix)

			if(index(runoptions,"verbose").gt.0) write(*,"('[L_gs] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
			call MPI_Send(gupgd,ncount,MPI_DOUBLE_COMPLEX,0,8999,MPI_COMM_WORLD,ierr)
			call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
			if(ix.eq.0) exit
		end do
	end if

	if(myrank.ne.0) then
		return
	end if

	call l_matrix(lx,ly,lz)

	lxpm = 0.d0
	lypm = 0.d0
	lzpm = 0.d0
	lxm = 0.d0
	lym = 0.d0
	lzm = 0.d0
	do nu=5,9 ; do mu=5,9 ; do i=1,Npl
		lxpm(i) = lxpm(i) + real(lxp(mu,nu)*gupgdint(i,mu,nu))
		lypm(i) = lypm(i) + real(lyp(mu,nu)*gupgdint(i,mu,nu))
		lzpm(i) = lzpm(i) + real(lzp(mu,nu)*gupgdint(i,mu,nu))
		lxm(i)  = lxm(i)  + real(lx(mu,nu)*gupgdint(i,mu,nu))
		lym(i)  = lym(i)  + real(ly(mu,nu)*gupgdint(i,mu,nu))
		lzm(i)  = lzm(i)  + real(lz(mu,nu)*gupgdint(i,mu,nu))
	end do ; end do ; end do

	return
end subroutine L_gs



subroutine sumk_L_gs(e,ep,gupgd)
	use mod_f90_kind
	use mod_constants, only: pi, zero
	use mod_parameters
	use mod_generate_kpoints
!$  use mod_mpi_pars, only: myrank
!$ 	use omp_lib
	implicit none
!$	integer				:: nthreads,mythread
	integer				:: i,mu,nu,mup,nup,iz
	real(double)	:: kp(3),e,ep
	complex(double),dimension(Npl,9,9)	:: gupgd
	complex(double) :: gf(Npl,Npl,18,18)

	gupgd		= zero

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,mu,nu,mup,nup) &
!$omp& shared(kbz,nkpoints,wkbz,e,ep,Npl,gupgd,myrank,nthreads)
!$  mythread = omp_get_thread_num()
!$  if((mythread.eq.0).and.(myrank.eq.0)) then
!$    nthreads = omp_get_num_threads()
!$    write(*,"('Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:gupgd)
	kpoints: do iz=1,nkpoints
		kp = kbz(iz,:)

		!Green function on energy Ef + iy, and wave vector kp
		call green(e,ep,kp,gf)

		do i=1,Npl
			do mu=1,9
				mup = mu+9
				do nu=1,9
					nup = nu+9
					gupgd(i,mu,nu) = gupgd(i,mu,nu) + ((gf(i,i,mu,nu)+gf(i,i,mup,nup))*wkbz(iz))
				end do
			end do
		end do
	end do kpoints
!$omp end do
!$omp end parallel

	gupgd = gupgd/pi

	return
end subroutine sumk_L_gs
