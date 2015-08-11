! ---------- Parallel spin current: Energy integration ---------
subroutine eintshe(e,tchiorbiikl,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)
	use mod_constants
	use mod_parameters
	use mod_lattice
	use mod_generate_epoints
	use mod_f90_kind
	use mod_mpi_pars
	use MPI
	implicit none
	integer 					:: AllocateStatus
	integer						:: i
	real(double)							:: start_time,elapsed_time,sizemat,speed
	real(double),intent(in)		:: e
	complex(double), dimension(dim,dimNpl),							intent(out) :: tchiorbiikl
	complex(double), dimension(n0sc1:n0sc2,dim,dimNpl),	intent(out)	:: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl
	complex(double), dimension(:,:),allocatable					:: tFintiikl
	complex(double), dimension(:,:,:), allocatable 			:: ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl
!--------------------- begin MPI vars --------------------
	integer :: ix,ix2,itask,masterrank=0
	integer :: ncountkl,nncountkl
	ncountkl=dim*dim*Npl
	nncountkl=n0sc*ncountkl
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

	allocate( tFintiikl(dim,dimNpl),ttFintiikl(n0sc1:n0sc2,dim,dimNpl),LxttFintiikl(n0sc1:n0sc2,dim,dimNpl),LyttFintiikl(n0sc1:n0sc2,dim,dimNpl),LzttFintiikl(n0sc1:n0sc2,dim,dimNpl), STAT = AllocateStatus )
	if (AllocateStatus.ne.0) then
		write(*,"('[eintshe] Not enough memory for: tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl')")
		call MPI_Abort(MPIComm_Row,errorcode,ierr)
	end if

! Generating energy points in the real axis for third integration
	call generate_real_epoints(e)

	ix = myrank_row+1
	itask = numprocs ! Number of tasks done initially

! Starting to calculate energy integral
	if (myrank_row.eq.masterrank) then ! Process 0 calculates one point, receives results from others and send new tasks if necessary
		call sumk(e,y(ix),tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,0)
		tchiorbiikl     = tFintiikl*wght(ix)
		ttchiorbiikl    = ttFintiikl*wght(ix)
		Lxttchiorbiikl  = LxttFintiikl*wght(ix)
		Lyttchiorbiikl  = LyttFintiikl*wght(ix)
		Lzttchiorbiikl  = LzttFintiikl*wght(ix)
		if(index(runoptions,"verbose").gt.0) write(*,"('[eintshe] Finished point ',i0,' in myrank_row ',i0,' myrank_col ',i0,' (',a,')')") ix,myrank_row,myrank_col,trim(host)

		do i=2,nepoints
			if(index(runoptions,"verbose").gt.0) start_time = MPI_Wtime()
			call MPI_Recv(tFintiikl,ncountkl,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,54545454+count,MPIComm_Row,stat,ierr)
			call MPI_Recv(ttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),65656565+count,MPIComm_Row,stat,ierr)
			call MPI_Recv(LxttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),76767676+count,MPIComm_Row,stat,ierr)
			call MPI_Recv(LyttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),87878787+count,MPIComm_Row,stat,ierr)
			call MPI_Recv(LzttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),98989898+count,MPIComm_Row,stat,ierr)
			if(index(runoptions,"verbose").gt.0) then
				elapsed_time = MPI_Wtime() - start_time
				write(*,"('Point ',i0,' received from ',i0,'. Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes ')") i,stat(MPI_SOURCE),elapsed_time,elapsed_time/60.d0
				sizemat = (ncountkl+4.d0*nncountkl)*16.d0/(1024.d0**2)
				speed 	= sizemat/elapsed_time
				write(*,"('Average speed: ',f11.4,' Mbps / ',f9.4,' Gbps ')") speed,speed/1024
			end if

			tchiorbiikl     = tchiorbiikl + tFintiikl
			ttchiorbiikl    = ttchiorbiikl + ttFintiikl
			Lxttchiorbiikl  = Lxttchiorbiikl + LxttFintiikl
			Lyttchiorbiikl  = Lyttchiorbiikl + LyttFintiikl
			Lzttchiorbiikl  = Lzttchiorbiikl + LzttFintiikl

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
				call sumk(e,y(ix),tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,0)
				tFintiikl  = tFintiikl*wght(ix)
				ttFintiikl = ttFintiikl*wght(ix)
				LxttFintiikl = LxttFintiikl*wght(ix)
				LyttFintiikl = LyttFintiikl*wght(ix)
				LzttFintiikl = LzttFintiikl*wght(ix)
			else if ((ix.gt.pn1).and.(ix.le.nepoints)) then ! Third integration (on the real axis)
				ix2 = ix-pn1
				call sumk(e,x2(ix2),tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,1)
				tFintiikl  = tFintiikl*p2(ix2)
				ttFintiikl = ttFintiikl*p2(ix2)
				LxttFintiikl = LxttFintiikl*p2(ix2)
				LyttFintiikl = LyttFintiikl*p2(ix2)
				LzttFintiikl = LzttFintiikl*p2(ix2)
			else
				exit
			end if

			if(index(runoptions,"verbose").gt.0) write(*,"('[eintshe] Finished point ',i0,' in myrank_row ',i0,' myrank_col ',i0,' (',a,')')") ix,myrank_row,myrank_col,trim(host)
			! Sending results to process 0
			call MPI_Send(tFintiikl,ncountkl,MPI_DOUBLE_COMPLEX,masterrank,54545454+count,MPIComm_Row,ierr)
			call MPI_Send(ttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,65656565+count,MPIComm_Row,ierr)
			call MPI_Send(LxttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,76767676+count,MPIComm_Row,ierr)
			call MPI_Send(LyttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,87878787+count,MPIComm_Row,ierr)
			call MPI_Send(LzttFintiikl,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,98989898+count,MPIComm_Row,ierr)
			! Receiving new point or signal to exit
			call MPI_Recv(ix,1,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPIComm_Row,stat,ierr)
			if(ix.eq.0) exit
		end do
	end if

	deallocate(ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl)

	return
end subroutine eintshe


! ---------- Parallel spin current: Energy integration ---------
subroutine eintshe_uff(e,tchiorbiikl,ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)
	use mod_constants
	use mod_parameters
	use mod_lattice
	use mod_generate_epoints
	use mod_f90_kind
	use mod_mpi_pars
	use MPI
	implicit none
	integer 					:: AllocateStatus
	integer						:: i,neighbor
	real(double)							:: start_time,elapsed_time,sizemat,speed
	real(double),intent(in)		:: e
	complex(double), dimension(dim,dimNpl),							intent(out) :: tchiorbiikl
	complex(double), dimension(n0sc1:n0sc2,dim,dimNpl),	intent(out)	:: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl
	complex(double), dimension(:,:),allocatable					:: tFintiikl
	complex(double), dimension(:,:,:), allocatable 			:: ttFintiikl,temp,LxttFintiikl,LyttFintiikl,LzttFintiikl
!--------------------- begin MPI vars --------------------
	integer :: ix,ix2,itask,masterrank=0
	integer :: ncountkl,nncountkl
	ncountkl=dim*dim*Npl
	nncountkl=2*ncountkl
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

	allocate( temp(2,dim,dimNpl),tFintiikl(dim,dimNpl),ttFintiikl(n0sc1:n0sc2,dim,dimNpl),LxttFintiikl(n0sc1:n0sc2,dim,dimNpl),LyttFintiikl(n0sc1:n0sc2,dim,dimNpl),LzttFintiikl(n0sc1:n0sc2,dim,dimNpl), STAT = AllocateStatus )
	if (AllocateStatus.ne.0) then
		write(*,"('[eintshe_uff] Not enough memory for: temp,tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl')")
		call MPI_Abort(MPIComm_Row,errorcode,ierr)
	end if

! Generating energy points in the real axis for third integration
	call generate_real_epoints(e)

	ix = myrank_row+1
	itask = numprocs ! Number of tasks done initially

! Starting to calculate energy integral
	if (myrank_row.eq.masterrank) then ! Process 0 calculates one point, receives results from others and send new tasks if necessary
		call sumk(e,y(ix),tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,0)
		tchiorbiikl     = tFintiikl*wght(ix)
		ttchiorbiikl    = ttFintiikl*wght(ix)
		Lxttchiorbiikl  = LxttFintiikl*wght(ix)
		Lyttchiorbiikl  = LyttFintiikl*wght(ix)
		Lzttchiorbiikl  = LzttFintiikl*wght(ix)
		if(index(runoptions,"verbose").gt.0) write(*,"('[eintshe] Finished point ',i0,' in myrank_row ',i0,' myrank_col ',i0,' (',a,')')") ix,myrank_row,myrank_col,trim(host)

		do i=2,nepoints
			if(index(runoptions,"verbose").gt.0) start_time = MPI_Wtime()
			call MPI_Recv(tFintiikl,ncountkl,MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,54545454+count,MPIComm_Row,stat,ierr)
			neighbor=n0sc1
			do
				call MPI_Recv(temp,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),65656565+neighbor+count*n0sc2,MPIComm_Row,stat,ierr)
				ttFintiikl(neighbor:neighbor+1,:,:)=temp
				call MPI_Recv(temp,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),76767676+neighbor+count*n0sc2,MPIComm_Row,stat,ierr)
				LxttFintiikl(neighbor:neighbor+1,:,:)=temp
				call MPI_Recv(temp,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),87878787+neighbor+count*n0sc2,MPIComm_Row,stat,ierr)
				LyttFintiikl(neighbor:neighbor+1,:,:)=temp
				call MPI_Recv(temp,nncountkl,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),98989898+neighbor+count*n0sc2,MPIComm_Row,stat,ierr)
				LzttFintiikl(neighbor:neighbor+1,:,:)=temp
				neighbor=neighbor+2
				if(neighbor.gt.n0sc2) exit
			end do
			if(index(runoptions,"verbose").gt.0) then
				elapsed_time = MPI_Wtime() - start_time
				write(*,"('Point ',i0,' received from ',i0,'. Elapsed time: ',f11.4,' seconds / ',f9.4,' minutes ')") i,stat(MPI_SOURCE),elapsed_time,elapsed_time/60.d0
				sizemat = (ncountkl+4.d0*nncountkl)*16.d0/(1024.d0**2)
				speed 	= sizemat/elapsed_time
				write(*,"('Average speed: ',f11.4,' Mbps / ',f9.4,' Gbps ')") speed,speed/1024
			end if

			tchiorbiikl     = tchiorbiikl + tFintiikl
			ttchiorbiikl    = ttchiorbiikl + ttFintiikl
			Lxttchiorbiikl  = Lxttchiorbiikl + LxttFintiikl
			Lyttchiorbiikl  = Lyttchiorbiikl + LyttFintiikl
			Lzttchiorbiikl  = Lzttchiorbiikl + LzttFintiikl

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
				call sumk(e,y(ix),tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,0)
				tFintiikl  = tFintiikl*wght(ix)
				ttFintiikl = ttFintiikl*wght(ix)
				LxttFintiikl = LxttFintiikl*wght(ix)
				LyttFintiikl = LyttFintiikl*wght(ix)
				LzttFintiikl = LzttFintiikl*wght(ix)
			else if ((ix.gt.pn1).and.(ix.le.nepoints)) then ! Third integration (on the real axis)
				ix2 = ix-pn1
				call sumk(e,x2(ix2),tFintiikl,ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,1)
				tFintiikl  = tFintiikl*p2(ix2)
				ttFintiikl = ttFintiikl*p2(ix2)
				LxttFintiikl = LxttFintiikl*p2(ix2)
				LyttFintiikl = LyttFintiikl*p2(ix2)
				LzttFintiikl = LzttFintiikl*p2(ix2)
			else
				exit
			end if

			if(index(runoptions,"verbose").gt.0) write(*,"('[eintshe] Finished point ',i0,' in myrank_row ',i0,' myrank_col ',i0,' (',a,')')") ix,myrank_row,myrank_col,trim(host)
			! Sending results to process 0
			call MPI_Send(tFintiikl,ncountkl,MPI_DOUBLE_COMPLEX,masterrank,54545454+count,MPIComm_Row,ierr)
			neighbor=n0sc1
			do
				temp=ttFintiikl(neighbor:neighbor+1,:,:)
				call MPI_Send(temp,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,65656565+neighbor+count*n0sc2,MPIComm_Row,ierr)
				temp=LxttFintiikl(neighbor:neighbor+1,:,:)
				call MPI_Send(temp,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,76767676+neighbor+count*n0sc2,MPIComm_Row,ierr)
				temp=LyttFintiikl(neighbor:neighbor+1,:,:)
				call MPI_Send(temp,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,87878787+neighbor+count*n0sc2,MPIComm_Row,ierr)
				temp=LzttFintiikl(neighbor:neighbor+1,:,:)
				call MPI_Send(temp,nncountkl,MPI_DOUBLE_COMPLEX,masterrank,98989898+neighbor+count*n0sc2,MPIComm_Row,ierr)
				neighbor=neighbor+2
				if(neighbor.gt.n0sc2) exit
			end do
			! Receiving new point or signal to exit
			call MPI_Recv(ix,1,MPI_INTEGER,masterrank,MPI_ANY_TAG,MPIComm_Row,stat,ierr)
			if(ix.eq.0) exit
		end do
	end if

	deallocate(ttFintiikl,LxttFintiikl,LyttFintiikl,LzttFintiikl,temp)

	return
end subroutine eintshe_uff
