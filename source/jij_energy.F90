! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
  subroutine jij_energy(Jij)
    use mod_f90_kind
    use mod_constants, only: pi
    use mod_parameters
    use mod_magnet
    use mod_generate_epoints
    use mod_mpi_pars
    use mod_progress
    use MPI
  integer         :: i
  real(double),dimension(nmaglayers,nmaglayers,3,3) :: Jijint
  real(double),dimension(nmaglayers,nmaglayers,3,3),intent(out) :: Jij
!--------------------- begin MPI vars --------------------
  integer :: ix,itask,ncount
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ncount=nmaglayers*nmaglayers*3*3

  Jij = 0.d0

  ix = myrank+1
  itask = numprocs ! Number of tasks done initially

  ! Calculating the number of particles for each spin and orbital using a complex integral
  if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
    call sumk_jij(Ef,y(ix),Jijint)
    Jij = Jijint*wght(ix)/pi

    if(lverbose) write(outputunit,"('[coupling] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
    do i=2,pn1
      if(lverbose) call progress_bar(outputunit,"energy points",i,pn1)

      call MPI_Recv(Jijint,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,11235,MPI_COMM_WORLD,stat,ierr)

      Jij = Jij + Jijint

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
      call sumk_jij(Ef,y(ix),Jijint)
      Jijint = Jijint*wght(ix)/pi

!       if(lverbose) write(outputunit,"('[coupling] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
      ! Sending results to process 0
      call MPI_Send(Jijint,ncount,MPI_DOUBLE_PRECISION,0,11235,MPI_COMM_WORLD,ierr)
      ! Receiving new point or signal to exit
      call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
      if(ix.eq.0) exit
    end do
  end if

  return
  end subroutine jij_energy