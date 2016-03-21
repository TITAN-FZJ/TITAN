! For a given value of center of band eps1 it calculates the
! occupation number and the magnetic moment
subroutine selfconsistency(N,x,fvec,iflag)
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
  integer         :: N,i,j,iflag,err,iflagrw
  real(double),dimension(N)           :: x,fvec
  real(double),dimension(Npl,9)       :: n_orb_u,n_orb_d,n_orb_t,mag_orb
  real(double),dimension(Npl)         :: n_t,mx_in,my_in,mz_in
  complex(double),dimension(Npl)      :: mp_in
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

  iflag=0
! Values used in the hamiltonian
  eps1  = x(1:Npl)
  mx_in = x(Npl+1:2*Npl)
  my_in = x(2*Npl+1:3*Npl)
  mz_in = x(3*Npl+1:4*Npl)
  mp_in = mx_in+zi*my_in
  do i=1,Npl
    hdel(i)   = 0.5d0*U(i+1)*mz_in(i)
    hdelp(i)  = 0.5d0*U(i+1)*mp_in(i)
  end do
  hdelm = conjg(hdelp)

  n_orb_u = 0.d0
  n_orb_d = 0.d0

  ix = myrank+1
  itask = numprocs ! Number of tasks done initially

  ! Calculating the number of particles for each spin and orbital using a complex integral
  if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
    write(*,"('|------------------- Iteration ',i4,' (densities) ------------------|')") iter
    call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
    gdiaguur = wght(ix)*gdiaguur
    gdiagddr = wght(ix)*gdiagddr
    gdiagud = wght(ix)*gdiagud
    gdiagdu = wght(ix)*gdiagdu

    n_orb_u = gdiaguur
    n_orb_d = gdiagddr

    do j=1,Npl
      mp(j) = (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
    end do

    if(lverbose) write(*,"('[selfconsistency] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
    do i=2,pn1
      ! Progress bar
      prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
      write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of nparticles e-sum on rank ',i0,a1,$)") spiner(mod(i,4)+1),prog,i,pn1,myrank,char(13)
#else
      elapsed_time = MPI_Wtime() - start_program
      write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),prog
#endif

      call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999,MPI_COMM_WORLD,stat,ierr)
      call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),9998,MPI_COMM_WORLD,stat,ierr)
      call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9997,MPI_COMM_WORLD,stat,ierr)
      call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9996,MPI_COMM_WORLD,stat,ierr)

      n_orb_u = n_orb_u + gdiaguur
      n_orb_d = n_orb_d + gdiagddr

      do j=1,Npl
        mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
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

      if(lverbose) write(*,"('[selfconsistency] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
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
  call MPI_Bcast(mp,Npl,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

  n_orb_u = 0.5d0 + n_orb_u/pi
  n_orb_d = 0.5d0 + n_orb_d/pi
  n_orb_t = n_orb_u + n_orb_d
  mag_orb = n_orb_u - n_orb_d
  mp      = mp/pi

  do i=1,Npl
    ! Number of particles
    n_t(i) = sum(n_orb_t(i,:))
    fvec(i)   = n_t(i) - npart0(i+1)
    ! x-component of magnetization
    j = i+Npl
    fvec(j)  = real(mp(i)) - mx_in(i)
    ! y-component of magnetization
    j = j+Npl
    fvec(j)  = aimag(mp(i)) - my_in(i)
    ! z-component of magnetization
    j = j+Npl
    mz(i)    = sum(mag_orb(i,5:9))
    fvec(j)  = mz(i) - mz_in(i)
  end do

  if(myrank.eq.0) then
    do i=1,Npl
      if(abs(mp(i)).gt.1.d-10) then
        write(*,"('Plane ',I2,': eps1(',I2,')=',e16.9,4x,'Mx(',I2,')=',e16.9,4x,'My(',I2,')=',e16.9,4x,'Mz(',I2,')=',e16.9)") i,i,eps1(i),i,real(mp(i)),i,aimag(mp(i)),i,mz(i)
        write(*,"(10x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9)") i,fvec(i),i+Npl,fvec(i+Npl),i+2*Npl,fvec(i+2*Npl),i+3*Npl,fvec(i+3*Npl)
      else
        write(*,"('Plane ',I2,': eps1(',I2,')=',e16.9,4x,'Mz(',I2,')=',e16.9)") i,i,eps1(i),i,mz(i)
        write(*,"(10x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9)") i,fvec(i),i+3*Npl,fvec(i+3*Npl)
      end if
    end do
  end if

  ! Writing new eps1 and mz to file while performing self-consistency
  if(myrank.eq.0) then
    iflagrw = 1
    call readwritesc(iflagrw,err)
  end if

  iter = iter + 1

  return
end subroutine selfconsistency

subroutine selfconsistencyjac(N,x,fvec,selfconjac,ldfjac,iflag)
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
  real(double)  :: x(N),fvec(N),selfconjac(ldfjac,N)
  real(double),dimension(Npl)         :: mx_in,my_in,mz_in
  complex(double),dimension(Npl)      :: mp_in
  real(double),dimension(4*Npl,4*Npl) :: ggr
  !--------------------- begin MPI vars --------------------
  integer :: ix,itask
  integer :: ncount
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
  ncount=16*Npl*Npl

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

  iflag=0
! Values used in the hamiltonian
  eps1  = x(1:Npl)
  mx_in = x(Npl+1:2*Npl)
  my_in = x(2*Npl+1:3*Npl)
  mz_in = x(3*Npl+1:4*Npl)
  mp_in = mx_in+zi*my_in
  do i=1,Npl
    hdel(i)   = 0.5d0*U(i+1)*mz_in(i)
    hdelp(i)  = 0.5d0*U(i+1)*mp_in(i)
  end do
  hdelm = conjg(hdelp)

  fvec=fvec

  ix = myrank+1
  itask = numprocs ! Number of tasks done initially

  selfconjac = 0.d0
  ! Calculating the jacobian
  if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
    write(*,"('|------------------- Iteration ',i4,' (jacobian) -------------------|')") iter
    call sumk_selfconjac(Ef,y(ix),ggr)
    selfconjac = wght(ix)*ggr

    if(lverbose) write(*,"('[selfconsistencyjac] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
    do i=2,pn1
      ! Progress bar
      prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
       write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of jacobian e-sum on rank ',i0,a1,$)") spiner(mod(i,4)+1),prog,i,pn1,myrank,char(13)
#else
      elapsed_time = MPI_Wtime() - start_program
      write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
      write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),prog
#endif

      call MPI_Recv(ggr,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3333,MPI_COMM_WORLD,stat,ierr)

      selfconjac = selfconjac + ggr

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
      call sumk_selfconjac(Ef,y(ix),ggr)
      ggr = wght(ix)*ggr

      if(lverbose) write(*,"('[selfconsistencyjac] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
      ! Sending results to process 0
      call MPI_Send(ggr,ncount,MPI_DOUBLE_PRECISION,0,3333,MPI_COMM_WORLD,ierr)
      ! Receiving new point or signal to exit
      call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
      if(ix.eq.0) exit
    end do
  end if

  ! Send results to all processors
  call MPI_Bcast(selfconjac,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  selfconjac = selfconjac/pi

  iter = iter + 1

  return
end subroutine selfconsistencyjac

! For a given value of center of band eps1 it calculates the
! occupation number and the magnetic moment
subroutine selfconsistencyjacnag(N,x,fvec,selfconjac,ldfjac,iflag)
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
  integer         :: N,ldfjac,i,j,iflag,err,iflagrw
  real(double),dimension(N)           :: x,fvec
  real(double),dimension(ldfjac,N)    :: selfconjac
  real(double),dimension(Npl)         :: n_t,mx_in,my_in,mz_in
  complex(double),dimension(Npl)      :: mp_in
  real(double),dimension(N,N)         :: ggr
  real(double),dimension(Npl,9)       :: n_orb_u,n_orb_d,n_orb_t,mag_orb
  real(double),dimension(Npl,9)       :: gdiaguur,gdiagddr
  complex(double),dimension(Npl,9)    :: gdiagud,gdiagdu
  !--------------------- begin MPI vars --------------------
  integer :: ix,itask
  integer :: ncount,ncount2
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
  ncount=Npl*9
  ncount2=16*Npl*Npl

#ifndef _JUQUEEN
  open(6,carriagecontrol ='fortran')
#endif

! Values used in the hamiltonian
  eps1  = x(1:Npl)
  mx_in = x(Npl+1:2*Npl)
  my_in = x(2*Npl+1:3*Npl)
  mz_in = x(3*Npl+1:4*Npl)
  mp_in = mx_in+zi*my_in
  do i=1,Npl
    hdel(i)   = 0.5d0*U(i+1)*mz_in(i)
    hdelp(i)  = 0.5d0*U(i+1)*mp_in(i)
  end do
  hdelm = conjg(hdelp)

  ix = myrank+1
  itask = numprocs ! Number of tasks done initially

  flag: select case (iflag)
  case(1)
    n_orb_u = 0.d0
    n_orb_d = 0.d0
    ! Calculating the number of particles for each spin and orbital using a complex integral
    if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
      write(*,"('|------------------- Iteration ',i4,' (densities) ------------------|')") iter
      call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
      gdiaguur = wght(ix)*gdiaguur
      gdiagddr = wght(ix)*gdiagddr
      gdiagud = wght(ix)*gdiagud
      gdiagdu = wght(ix)*gdiagdu

      n_orb_u = gdiaguur
      n_orb_d = gdiagddr

      do j=1,Npl
        mp(j) = (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
      end do

      if(lverbose) write(*,"('[selfconsistencyjacnag1] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
      do i=2,pn1
        ! Progress bar
        prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
        write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of nparticles e-sum on rank ',i0,a1,$)") spiner(mod(i,4)+1),prog,i,pn1,myrank,char(13)
#else
        elapsed_time = MPI_Wtime() - start_program
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
        write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),prog
#endif

        call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999+iter,MPI_COMM_WORLD,stat,ierr)
        call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),8998+iter,MPI_COMM_WORLD,stat,ierr)
        call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),7997+iter,MPI_COMM_WORLD,stat,ierr)
        call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),6996+iter,MPI_COMM_WORLD,stat,ierr)

        n_orb_u = n_orb_u + gdiaguur
        n_orb_d = n_orb_d + gdiagddr

        do j=1,Npl
          mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
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

        if(lverbose) write(*,"('[selfconsistencyjacnag1] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
        ! Sending results to process 0
        call MPI_Send(gdiaguur,ncount,MPI_DOUBLE_PRECISION,0,9999+iter,MPI_COMM_WORLD,ierr)
        call MPI_Send(gdiagddr,ncount,MPI_DOUBLE_PRECISION,0,8998+iter,MPI_COMM_WORLD,ierr)
        call MPI_Send(gdiagud,ncount,MPI_DOUBLE_COMPLEX,0,7997+iter,MPI_COMM_WORLD,ierr)
        call MPI_Send(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,0,6996+iter,MPI_COMM_WORLD,ierr)
        ! Receiving new point or signal to exit
        call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
        if(ix.eq.0) exit
      end do
    end if

    ! Send results to all processors
    call MPI_Bcast(n_orb_u,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(n_orb_d,ncount,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(mp,Npl,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

    n_orb_u = 0.5d0 + n_orb_u/pi
    n_orb_d = 0.5d0 + n_orb_d/pi
    n_orb_t = n_orb_u + n_orb_d
    mag_orb = n_orb_u - n_orb_d
    mp      = mp/pi

    do i=1,Npl
      ! Number of particles
      n_t(i) = sum(n_orb_t(i,:))
      fvec(i)   = n_t(i) - npart0(i+1)
      ! x-component of magnetization
      j = i+Npl
      fvec(j)  = real(mp(i)) - mx_in(i)
      ! y-component of magnetization
      j = j+Npl
      fvec(j)  = aimag(mp(i)) - my_in(i)
      ! z-component of magnetization
      j = j+Npl
      mz(i)    = sum(mag_orb(i,5:9))
      fvec(j)  = mz(i) - mz_in(i)
    end do

    if(myrank.eq.0) then
      do i=1,Npl
        if(abs(mp(i)).gt.1.d-10) then
          write(*,"('Plane ',I2,': eps1(',I2,')=',e16.9,4x,'Mx(',I2,')=',e16.9,4x,'My(',I2,')=',e16.9,4x,'Mz(',I2,')=',e16.9)") i,i,eps1(i),i,real(mp(i)),i,aimag(mp(i)),i,mz(i)
          write(*,"(10x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9)") i,fvec(i),i+Npl,fvec(i+Npl),i+2*Npl,fvec(i+2*Npl),i+3*Npl,fvec(i+3*Npl)
        else
          write(*,"('Plane ',I2,': eps1(',I2,')=',e16.9,4x,'Mz(',I2,')=',e16.9)") i,i,eps1(i),i,mz(i)
          write(*,"(10x,'fvec(',I2,')=',e16.9,2x,'fvec(',I2,')=',e16.9)") i,fvec(i),i+3*Npl,fvec(i+3*Npl)
        end if
      end do
    end if

    ! Writing new eps1 and mz to file while performing self-consistency
    if(myrank.eq.0) then
      iflagrw = 1
      call readwritesc(iflagrw,err)
    end if

  case(2)
    selfconjac = 0.d0
    ! Calculating the jacobian
    if (myrank.eq.0) then ! Process 0 receives all results and send new tasks if necessary
      write(*,"('|------------------- Iteration ',i4,' (jacobian) -------------------|')") iter
      call sumk_selfconjac(Ef,y(ix),ggr)
      selfconjac = wght(ix)*ggr

      if(lverbose) write(*,"('[selfconsistencyjacnag2] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
      do i=2,pn1
        ! Progress bar
        prog = floor(i*100.d0/pn1)
#ifdef _JUQUEEN
         write(*,"(a1,2x,i3,'% (',i0,'/',i0,') of jacobian e-sum on rank ',i0,a1,$)") spiner(mod(i,4)+1),prog,i,pn1,myrank,char(13)
#else
        elapsed_time = MPI_Wtime() - start_program
        write(progbar,fmt="( a,i0,a )") "(1h+' ','Total time=',i2,'h:',i2,'m:',i2,'s  ',",1+(i+1)*20/pn1, "a,' ',i0,'%')"
        write(6,fmt=progbar) int(elapsed_time/3600.d0),int(mod(elapsed_time,3600.d0)/60.d0),int(mod(mod(elapsed_time,3600.d0),60.d0)),("|",j=1,1+(i+1)*20/pn1),prog
#endif

        call MPI_Recv(ggr,ncount2,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3333+iter,MPI_COMM_WORLD,stat,ierr)

        selfconjac = selfconjac + ggr

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
        call sumk_selfconjac(Ef,y(ix),ggr)
        ggr = wght(ix)*ggr

        if(lverbose) write(*,"('[selfconsistencyjacnag2] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank,trim(host)
        ! Sending results to process 0
        call MPI_Send(ggr,ncount2,MPI_DOUBLE_PRECISION,0,3333+iter,MPI_COMM_WORLD,ierr)
        ! Receiving new point or signal to exit
        call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_COMM_WORLD,stat,ierr)
        if(ix.eq.0) exit
      end do
    end if

    ! Send results to all processors
    call MPI_Bcast(selfconjac,ncount2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    selfconjac = selfconjac/pi
    do i = Npl+1,4*Npl
      selfconjac(i,i) = selfconjac(i,i) - 1.d0
    end do
  case default
    write(*,"('[selfconsistencyjacnag] Problem in self-consistency! iflag = ',I0)") iflag
    stop
  end select flag

  iter = iter + 1

  return
end subroutine selfconsistencyjacnag