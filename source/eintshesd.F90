! ---------- Spin disturbance: Energy integration ---------
subroutine eintshesd(e)
	use mod_f90_kind
	use mod_constants
	use mod_parameters
	use mod_generate_epoints
  use mod_disturbances, only: tchiorbiikl
	use mod_mpi_pars
	use MPI
	implicit none
	integer						:: AllocateStatus
	integer						:: i
	real(double)							:: start_time,elapsed_time,sizemat,speed
	real(double), intent(in)	:: e
	complex(double), dimension(:,:),allocatable					:: tFintiikl
!--------------------- begin MPI vars --------------------
	integer :: ix,ix2,itask,masterrank=0
	integer :: ncountkl
	ncountkl=dim*4
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

	allocate( tFintiikl(dim,4), STAT = AllocateStatus )
	if (AllocateStatus.ne.0) then
		write(*,"('[eintshesd] Not enough memory for: tFintiikl ')")
		call MPI_Abort(MPIComm_Row,errorcode,ierr)
	end if

! Generating energy points in the real axis for third integration
	call generate_real_epoints(e)

	ix = myrank_row+1
	itask = numprocs ! Number of tasks done initially

! Starting to calculate energy integral
	if (myrank_row.eq.masterrank) then ! Process 0 receives all results and send new tasks if necessary
		call sumkshesd(e,y(ix),tFintiikl,0)
		tchiorbiikl  = tFintiikl*wght(ix)
		if(lverbose) write(*,"('[eintshesd] Finished point ',i0,' in myrank_row ',i0,' myrank_col ',i0,' (',a,')')") ix,myrank_row,myrank_col,trim(host)

		do i=2,nepoints
			if(lverbose) start_time = MPI_Wtime()
			call MPI_Recv(tFintiikl,ncountkl,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,32323232+mpitag,MPIComm_Row,stat,ierr)
			if(lverbose) then
				elapsed_time = MPI_Wtime() - start_time
				write(*,"('Point ',i0,' received from ',i0,'. Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes ')") i,stat(MPI_SOURCE),elapsed_time,elapsed_time/60.d0
				sizemat = ncountkl*16.d0/(1024.d0**2)
				speed 	= sizemat/elapsed_time
				write(*,"('Average speed: ',f11.4,' Mbps / ',f9.4,' Gbps ')") speed,speed/1024.d0
			end if

			tchiorbiikl = tchiorbiikl + tFintiikl

			! If the number of processors is less than the total number of points, sends
			! the rest of the points to the ones that finish first
			if (itask.lt.nepoints) then
				itask = itask + 1
				call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPIComm_Row,ierr)
			else
				call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPIComm_Row,ierr)
			end if
		end do
	else
		! Other processors calculate each point of the integral and waits for new points
		do
			if (ix.le.pn1) then ! First and second integrations (in the complex plane)
				call sumkshesd(e,y(ix),tFintiikl,0)
				tFintiikl = tFintiikl*wght(ix)
			else if ((ix.gt.pn1).and.(ix.le.nepoints)) then ! Third integration (on the real axis)
				ix2 = ix-pn1
				call sumkshesd(e,x2(ix2),tFintiikl,1)
				tFintiikl = tFintiikl*p2(ix2)
			else
				exit
			end if

			if(lverbose) write(*,"('[eintshesd] Finished point ',i0,' in myrank_row ',i0,' myrank_col ',i0,' (',a,')')") ix,myrank_row,myrank_col,trim(host)
			! Sending results to process 0
			call MPI_Send(tFintiikl,ncountkl,MPI_DOUBLE_COMPLEX,masterrank,32323232+mpitag,MPIComm_Row,ierr)
			! Receiving new point or signal to exit
			call MPI_Recv(ix,1,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPIComm_Row,stat,ierr)
			if(ix.eq.0) exit
		end do
	end if

	deallocate(tFintiikl)

	return
end subroutine eintshesd