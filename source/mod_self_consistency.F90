module mod_self_consistency
   use mod_f90_kind, only: double
   implicit none
   integer            :: neq
   character(len=300) :: default_file
   real(double)       :: mag_tol = 1.d-10
   character(len=200) :: scfile = ""
   !! Give a file to start self-consistency
   logical :: skipsc
   !! Skip self-consistency
   logical :: lselfcon    = .false.
   logical :: lGSL        = .false.
   logical :: lontheflysc = .false.
   logical :: lslatec     = .false.
   logical :: lnojac      = .false.
   logical :: lrotatemag  = .false.

contains

  subroutine doSelfConsistency()
    use mod_magnet,     only: lp_matrix, mtheta, mphi, lb_matrix, sb_matrix
    use adaptiveMesh,   only: genLocalEKMesh, freeLocalEKMesh
    use mod_mpi_pars,   only: rField, sField, FieldComm
    use mod_SOC,        only: SOC
    use mod_System,     only: s => sys
    use TightBinding,   only: nOrb
    implicit none
    logical :: lsuccess = .false.

    ! Distribute Energy Integration across all points available
    call genLocalEKMesh(s,rField,sField, FieldComm)

    !--------------------------- Self-consistency --------------------------
    ! Trying to read previous densities and Ef from files
    call read_previous_results(lsuccess)
    ! Rotate the magnetization to the direction of the field
    ! (useful for SOC=F)
    ! Only when some previous result was read from file (lsuccess=.true.)
    if(lrotatemag .and. lsuccess) then
      call rotate_magnetization_to_field()
    end if

    if(lselfcon) call calcMagneticSelfConsistency()

    ! Writing new n and mz to file after self-consistency is done
    if(.not. lontheflysc) call write_sc_results()

    ! L matrix in local frame for given quantization direction
    call lp_matrix(mtheta, mphi)

    ! Calculating ground state Orbital Angular Momentum
    if((lGSL).or.(SOC)) call calcLGS()

    ! Writing self-consistency results on screen
    if(rField == 0)  call print_sc_results()

    call freeLocalEKMesh()
  end subroutine doSelfConsistency

  ! Tries to read n and m if available
  subroutine read_previous_results(lsuccess)
    use mod_f90_kind,   only: double
    use mod_constants,  only: deg2rad
    use mod_parameters, only: output, magaxis, magaxisvec, offset, layertype
    use mod_system,     only: s => sys
    use TightBinding,   only: nOrb
    use mod_mpi_pars,   only: abortProgram, rField
    use mod_magnet,     only: mx,my,mz,mxd,myd,mzd,mpd,hw_count,hw_list, &
                              lfield,rho,rhod,rhod0,rho0
    use mod_Umatrix
    implicit none
    integer             :: i,err
    logical,intent(out) :: lsuccess

    lsuccess = .false.
    call read_sc_results(err,lsuccess)

    if(lsuccess) then
      if(err==0) then ! Same parameters
        if(skipsc) then ! Skip option ON
          if(rField == 0) write(output%unit_loop,"('[read_previous_results] Existing results for the same parameters were read. Skipping self-consistency...')")
          lselfcon = .false.
        else ! Skip option OFF
          if(rField == 0) write(output%unit_loop,"('[read_previous_results] Existing results for the same parameters were read. Updating values...')")
          lselfcon = .true.
        end if
      else ! Other parameters
        if(rField == 0) write(output%unit_loop,"('[read_previous_results] Existing results for other parameters were read. Updating values...')")
        lselfcon = .true.
      end if
    else !  If file doesn't exist
      if(rField == 0) then
        write(output%unit_loop,"('[read_previous_results] Self-consistency file does not exist:')")
        write(output%unit_loop,"('[read_previous_results] ',a)") trim(default_file)
      end if
      lselfcon = .true.
      ! Parameters: magnetization, exchange split
      if(magaxis == -1) then
        continue
      else if(magaxis == -2) then
        magaxisvec = magaxisvec(1) * s%a1 + magaxisvec(2) * s%a2 + magaxisvec(3) * s%a3
      else if(magaxis == -3) then
        magaxisvec = [cos(magaxisvec(2)*deg2rad)*sin(magaxisvec(1)*deg2rad), sin(magaxisvec(2)*deg2rad)*sin(magaxisvec(1)*deg2rad), cos(magaxisvec(1)*deg2rad)]
      else if(magaxis == 0) then
        magaxisvec = [0.d0, 0.d0, sign(1.0d0, hw_list(hw_count,1))]
      else if(magaxis >=1 .and. magaxis <= s%nAtoms) then
        !magaxisvec(1:3) = c_nn(1:3, magaxis)
        if(rField == 0) call abortProgram("[read_previous_results] Magaxis along neighbor not implemented!")
      else
        if(rField == 0) call abortProgram("[read_previous_results] Unknown magnetization direction!")
      end if
      magaxisvec = magaxisvec / sqrt(dot_product(magaxisvec, magaxisvec))
      magaxisvec = magaxisvec * 0.5d0

      mxd = magaxisvec(1)
      myd = magaxisvec(2)
      mzd = magaxisvec(3)

      do i=1,s%nAtoms
        if(layertype(i+offset)==2) then
          mxd(i) = mxd(i) * sign(4.d0,hw_list(hw_count,1))
          myd(i) = myd(i) * sign(4.d0,hw_list(hw_count,1))
          mzd(i) = mzd(i) * sign(4.d0,hw_list(hw_count,1))
        end if
      end do

      if(lfield .and. magaxis == 0) then
        mxd = mzd * sin(hw_list(hw_count,2)*deg2rad) * cos(hw_list(hw_count,3)*deg2rad)
        myd = mzd * sin(hw_list(hw_count,2)*deg2rad) * sin(hw_list(hw_count,3)*deg2rad)
        mzd = mzd * cos(hw_list(hw_count,2)*deg2rad)
      end if

      mpd = cmplx(mxd,myd)

      mx = 0.d0
      my = 0.d0
      mz = 0.d0
      ! Initialize U matrix using occupations from elemental files and rho=rho0
      rho = rho0
      do i = 1, s%nAtoms
        mx(5:9,i) = 0.2d0*mxd(i)
        my(5:9,i) = 0.2d0*myd(i)
        mz(5:9,i) = 0.2d0*mzd(i)
        rhod(i)   = s%Types(s%Basis(i)%Material)%OccupationD
      end do

    end if

    call init_Umatrix(mzd,mpd,rhod,rhod0,rho,rho0,s%nAtoms,nOrb)
  end subroutine read_previous_results

  ! This subroutine reads previous band-shifting and magnetization results
  subroutine read_sc_results(err,lsuccess)
    use mod_f90_kind,      only: double
    use mod_parameters,    only: output, dfttype
    use EnergyIntegration, only: parts
    use TightBinding,      only: nOrb
    use mod_system,        only: s => sys
    use mod_magnet,        only: rho, mp, mx, my, mz, rhod, &
                                 mpd, mxd, myd, mzd, hw_count
    use mod_mpi_pars
    implicit none
    character(len=300)  :: file = ""
    integer,intent(out) :: err
    logical,intent(out) :: lsuccess
    integer             :: i,j
    real(double)        :: previous_results(4*nOrb,s%nAtoms), previous_Ef

    if(trim(scfile) /= "") then
      open(unit=99,file=scfile,status="old",iostat=err)
      if(err/=0) then
        if(rField == 0) write(output%unit_loop,"('*** WARNING: Self-consistency file given on input file does not exist! Using default... ***')")
        scfile = " "
      end if
      close(99)
    end if

    lsuccess = .false.
    !   Reading previous results (mx, my, mz and n) from files (if available)
    if(trim(scfile)=="") then ! If a filename is not given in inputcard (or don't exist), use the default one
      write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),dfttype,parts,trim(output%BField),trim(output%info),trim(output%SOC)
      open(unit=99,file=file,status="old",iostat=err)
      if((err==0).and.(rField==0)) then
        write(output%unit_loop,"('[read_sc_results] Self-consistency file already exists. Reading it now...')")
        write(output%unit_loop,"(a)") trim(file)
      else
        default_file = trim(file)
      end if
    else ! If filename in inputcard exists or 2nd+ angular iteration
      if(((hw_count)==1)) then ! Filename in inputcard (1st iteration on loop)
        open(unit = 99,file = scfile, status = "old", iostat = err)
        if(err==0 .and. rField==0) then
          write(output%unit_loop,"('[read_sc_results] Using filename given in input file for self-consistency:')")
          write(output%unit_loop,"(a)") trim(scfile)
        end if
      else ! 2nd+ iteration, cheking if default file exists
        write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),dfttype,parts,trim(output%BField),trim(output%info),trim(output%SOC)
        open(unit=99,file=file,status="old",iostat=err)
        if(err == 0) then ! Reading file for the same parameters
          if(rField == 0) then
            write(output%unit_loop,"('[read_sc_results] Self-consistency file already exists. Reading it now...')")
            write(output%unit_loop,"(a)") trim(file)
          end if
        else ! Reading results from previous iteration
          open(unit=99,file=scfile,status="old",iostat=err)
          if((err==0).and.(rField==0)) then
            write(output%unit_loop,"('[read_sc_results] Using results from previous iteration as input for self-consistency:')")
            write(output%unit_loop,"(a)") trim(scfile)
          end if
          lsuccess   = .true. ! something was read
        end if
      end if
    end if
    if(err==0) then
      if(rField==0) then
        do i=1,s%nAtoms
          read(99,fmt=*) (previous_results(j,i), j=1,4*nOrb)
        end do
        read(99,fmt=*) previous_Ef
      end if

      call MPI_Bcast(previous_results,4*nOrb*s%nAtoms,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
      call MPI_Bcast(previous_Ef,1,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

      rho(:,:) = previous_results(       1:  nOrb,:)
      mx (:,:) = previous_results(  nOrb+1:2*nOrb,:)
      my (:,:) = previous_results(2*nOrb+1:3*nOrb,:)
      mz (:,:) = previous_results(3*nOrb+1:4*nOrb,:)
      mp       = cmplx(mx,my)

      rhod(:) = sum(rho(5:9,:),dim=1)
      mxd(:)  = sum(mx (5:9,:),dim=1)
      myd(:)  = sum(my (5:9,:),dim=1)
      mzd(:)  = sum(mz (5:9,:),dim=1)
      mpd     = cmplx(mxd,myd)

      call calcMagAngle()

      s%Ef = previous_Ef
      if(lsuccess) then
        err = 1   ! Read different parameters
      else
        lsuccess   = .true. ! Read same parameters (err=0)
      end if
    else
      ! If file does not exist, try to read for parts-1
      close(99)
      write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),dfttype,parts-1,trim(output%BField),trim(output%info),trim(output%SOC)
      open(unit=99,file=file,status="old",iostat=err)
      if(err==0) then
        if(rField==0) then
          write(output%unit_loop,"('[read_sc_results] Self-consistency file does not exist. Reading results for parts-1 now...')")
          write(output%unit_loop,"('[read_sc_results] Updating values obtained for parts-1...')")
          write(output%unit_loop,"(a)") file
          do i=1,s%nAtoms
            read(99,fmt=*) (previous_results(j,i), j=1,4*nOrb)
          end do
          read(99,fmt=*) previous_Ef
        end if

        call MPI_Bcast(previous_results,4*nOrb*s%nAtoms,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
        call MPI_Bcast(previous_Ef,1,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

        rho(:,:) = previous_results(       1:  nOrb,:)
        mx (:,:) = previous_results(  nOrb+1:2*nOrb,:)
        my (:,:) = previous_results(2*nOrb+1:3*nOrb,:)
        mz (:,:) = previous_results(3*nOrb+1:4*nOrb,:)
        mp       = cmplx(mx,my)

        rhod(:) = sum(rho(5:9,:),dim=1)
        mxd(:)  = sum(mx (5:9,:),dim=1)
        myd(:)  = sum(my (5:9,:),dim=1)
        mzd(:)  = sum(mz (5:9,:),dim=1)
        mpd     = cmplx(mxd,myd)

        call calcMagAngle()

        s%Ef = previous_Ef
        lsuccess = .true. ! Read...
        err = 1           ! ... different parameters
      end if
    end if
    close(99)
  end subroutine read_sc_results

  subroutine calcMagAngle()
    use mod_constants,        only: pi,rad2deg
    use mod_system,           only: s => sys
    use mod_susceptibilities, only: lrot
    use mod_magnet,           only: mx, my, mz, mabs, &
                                    mtheta, mphi, mvec_cartesian, &
                                    mvec_spherical
    implicit none
    integer :: i

    ! Calculating new angles of GS magnetization in units of pi and magnetization vector
    do i = 1,s%nAtoms
      mabs(i)   = sqrt((sum(mx(:,i))**2)+(sum(my(:,i))**2)+(sum(mz(:,i))**2))
      if(abs(mabs(i))>1.d-8) then
        mtheta(i) = acos(sum(mz(:,i))/mabs(i))*rad2deg
      else
        mtheta(i) = 0.d0
      end if
      if(abs(mtheta(i))>1.d-8) then
        if(abs(abs(mtheta(i))-180.d0)>1.d-8) then
          mphi(i)   = atan2(sum(my(:,i)),sum(mx(:,i)))*rad2deg
        else
          mphi(i) = 0.d0
        end if
        lrot = .true. ! Susceptibilities need to be rotated
      else
        mphi(i) = 0.d0
      end if
      mvec_cartesian(1,i) = sum(mx(:,i))
      mvec_cartesian(2,i) = sum(my(:,i))
      mvec_cartesian(3,i) = sum(mz(:,i))
      mvec_spherical(1,i) = mabs(i)
      mvec_spherical(2,i) = mtheta(i)
      mvec_spherical(3,i) = mphi(i)
    end do
  end  subroutine calcMagAngle

  subroutine calcMagneticSelfConsistency()
  !! This subroutine performs the self-consistency
    use mod_f90_kind,   only: double
    use mod_constants,  only: pi
    use mod_parameters, only: output
    use mod_magnet,     only: iter, rho, mxd, myd, mzd
    use mod_mpi_pars,   only: rField
    use TightBinding,   only: nOrb
    use mod_system,     only: s => sys
    use adaptiveMesh
    use mod_dnsqe
    implicit none
    real(double),allocatable      :: fvec(:),jac(:,:),wa(:),sc_solu(:)
    real(double),allocatable      :: diag(:),qtf(:)
    real(double)                  :: epsfcn,factor
#if !defined(_OSX) && !defined(_JUQUEEN)
    real(double)                  :: ruser(1)
    integer                       :: iuser(1)
#else
    real(double),allocatable :: w(:,:)
#endif
    integer                       :: i,mu,maxfev,ml,mr,mode,nfev,njev,lwa,ifail=0

    neq = 8*s%nAtoms+1
    allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )

    ! Putting read n and m existing solutions into sc_solu (first guess of the subroutine)
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        sc_solu((i-1)*8+(mu-4)) = rho(mu,i)
      end do
      sc_solu((i-1)*8+6) = mxd(i)
      sc_solu((i-1)*8+7) = myd(i)
      sc_solu((i-1)*8+8) = mzd(i)
    end do
    sc_solu(8*s%nAtoms+1) = s%Ef
    iter  = 1

    if(rField == 0) &
    write(output%unit_loop,"('[self_consistency] Starting self-consistency:')")

#if defined(_OSX) || defined(_JUQUEEN)
    if(lslatec) then
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa),w(neq,4) )
      if(lnojac) then
        call dnsqe(sc_eqs_old,sc_jac_old,2,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs_old,sc_jac_old,1,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      end if
      ifail = ifail-1
    else
      lwa=neq*(neq+1)/2
      allocate( wa(lwa),w(neq,4) )
      if(lnojac) then
!         call c05nbf(sc_equations,neq,sc_solu,fvec,mag_tol,wa,lwa,ifail)
        maxfev = 200*(neq+1)
        ml = neq-1
        mr = neq-1
        epsfcn = 1.d-5
        mode = 1
        factor = 0.1d0
        call c05ncf(sc_eqs_old,neq,sc_solu,fvec,mag_tol,maxfev,ml,mr,epsfcn,diag,mode,factor,0,nfev,jac,neq,wa,lwa,qtf,w,ifail)
      else
!         call c05pbf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,neq,mag_tol,wa,lwa,ifail)
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
        factor = 0.1d0
        call c05pcf(sc_eqs_and_jac_old,neq,sc_solu,fvec,jac,neq,mag_tol,maxfev,diag,mode,factor,0,nfev,njev,wa,lwa,qtf,w,ifail)
      end if
    end if
    deallocate( w )
#else
    if(lslatec) then
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa) )
      if(lnojac) then
        call dnsqe(sc_eqs_old,sc_jac_old,2,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs_old,sc_jac_old,1,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      end if
      ifail = ifail-1
    else
      lwa=neq*(neq+1)/2
      allocate( wa(lwa) )
      if(lnojac) then
        maxfev = 200*(neq+1)
        ml = neq-1
        mr = neq-1
        epsfcn = 1.d-5
        mode = 1
        factor = 0.1d0
        call c05qcf(sc_equations,neq,sc_solu,fvec,mag_tol,maxfev,ml,mr,epsfcn,mode,diag,factor,0,nfev,jac,wa,qtf,iuser,ruser,ifail)
      else
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
        factor = 0.1d0
        call c05rcf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,mag_tol,maxfev,mode,diag,factor,0,nfev,njev,wa,qtf,iuser,ruser,ifail)
      end if
    end if
#endif

    deallocate(sc_solu,diag,qtf,fvec,jac,wa)

    ! Calculating the magnetization in cartesian and spherical coordinates
    call calcMagAngle()
  end subroutine calcMagneticSelfConsistency


  subroutine check_jacobian(neq,x)
    use mod_parameters, only: output,lcheckjac
    use mod_mpi_pars,   only: rField,ierr
    use mod_chkder
    implicit none
    integer :: liw,lw
    integer :: i,ifail
    integer     , allocatable :: iw(:)
    real(double), allocatable :: w(:)
    integer     , intent(in)  :: neq
    real(double), intent(in)  :: x(neq)
    real(double) :: fvec(neq),fvecp(neq),jac(neq,neq),xp(neq),err(neq)
#if !defined(_OSX) && !defined(_JUQUEEN)
    real(double) :: ruser(1)
    integer      :: iuser(1)
#endif

    liw = 1
    lw  = (4+neq)*neq
    allocate(iw(liw),w(lw))

    lcheckjac = .false.

    if(rField == 0) &
    write(output%unit_loop,"('[check_jacobian] Checking Jacobian if Jacobian is correct...')", advance='no')

    call e04yaf(neq,neq,lsqfun,x,fvec,jac,neq,iw,liw,w,lw,ifail)

    ! if(rField == 0) write(*,*) ifail

!     call chkder(neq,neq,x,fvec,jac,neq,xp,fvecp,1,err)

! #if defined(_OSX) || defined(_JUQUEEN)
!     call sc_eqs_and_jac_old(neq,x ,fvec ,jac,neq,1)
!     call sc_eqs_and_jac_old(neq,x ,fvec ,jac,neq,2)
!     call sc_eqs_and_jac_old(neq,xp,fvecp,jac,neq,1)
! #else
!     call sc_equations_and_jacobian(neq,x ,fvec ,jac,iuser,ruser,1)
!     call sc_equations_and_jacobian(neq,x ,fvec ,jac,iuser,ruser,2)
!     call sc_equations_and_jacobian(neq,xp,fvecp,jac,iuser,ruser,1)
! #endif

!     call chkder(neq,neq,x,fvec,jac,neq,xp,fvecp,2,err)

    ! do i = 1,neq
    !   fvecp(i) = fvecp(i) - fvec(i)
    ! end do

    ! if(rField == 0) then
    !   ! write(*,*) "fvec"
    !   ! write(*,*) (fvec (i),i=1,neq)
    !   ! write(*,*) "fvecp"
    !   ! write(*,*) (fvecp(i),i=1,neq)
    !   do i=1,neq
    !     if(abs(err(i)-1.d0)>1.d-8) write(*,*) i,err(i)
    !   end do
    ! end if

    if(ifail == 0) then
      if(rField == 0) write(output%unit_loop,"(' YES! ')")
    else
      if(rField == 0) write(output%unit_loop,"(' NO! ifail = ',i0)") ifail
    end if

    lcheckjac = .true.

  end subroutine check_jacobian

  subroutine lsqfun(iflag,M,N,x,fvec,selfconjac,ljc,iw,liw,w,lw)
    use mod_f90_kind,   only: double
    use mod_system,     only: s => sys
    use TightBinding,   only: nOrb
    use mod_magnet,     only: iter,rho,rhod,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix,    only: update_Umatrix
    use mod_tools,      only: itos
    use mod_mpi_pars
    implicit none
    integer  :: M,N,ljc,i,mu,iflag,lw,liw,iw(liw)
    real(double) :: w(lw)
    real(double),   dimension(N)             :: x,fvec
    real(double),   dimension(N,N)           :: selfconjac
    real(double),   dimension(nOrb,s%nAtoms) :: rho_in
    real(double),   dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms)      :: mpd_in

    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*8+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*8+6)
      myd_in(i) = x((i-1)*8+7)
      mzd_in(i) = x((i-1)*8+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i))
    end do
    s%Ef    = x(8*s%nAtoms+1)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)

    ! call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    call calcMagnetization()
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        fvec((i-1)*8+(mu-4)) = rho(mu,i) - rho_in(mu,i)
      end do
      fvec((i-1)*8+6) =  mxd(i) -  mxd_in(i)
      fvec((i-1)*8+7) =  myd(i) -  myd_in(i)
      fvec((i-1)*8+8) =  mzd(i) -  mzd_in(i)
    end do
    fvec(8*s%nAtoms+1) = sum(rho) - s%totalOccupation

    call calcJacobian(selfconjac, N)

  end subroutine lsqfun

  subroutine calcMagnetization()
    !! Calculates occupation density and magnetization.
    use mod_f90_kind,      only: double
    use mod_constants,     only: cI,pi,cZero
    use mod_SOC,           only: llinearsoc,llineargfsoc
    use EnergyIntegration, only: y,wght
    use mod_system,        only: s => sys
    use mod_magnet,        only: mx,my,mz,mp,rho,mxd,myd,mzd,mpd,rhod
    use adaptiveMesh,      only: bzs,E_k_imag_mesh,activeComm,local_points
    use TightBinding,      only: nOrb,nOrb2
    use mod_parameters,    only: eta
    use mod_mpi_pars
    implicit none
    integer  :: i,j, AllocateStatus
    real(double),    dimension(3)                    :: kp
    complex(double), dimension(:,:),     allocatable :: gdiagud,gdiagdu
    real(double),    dimension(:,:),     allocatable :: imguu,imgdd
    complex(double), dimension(:,:,:,:), allocatable :: gf
    !--------------------- begin MPI vars --------------------
    integer*8 :: ix
    integer :: ncount
    integer :: mu,mup
    real(double) :: weight, ep
    ncount = s%nAtoms * nOrb

    allocate(imguu(nOrb,s%nAtoms),imgdd(nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) &
    call abortProgram("[calcMagnetization] Not enough memory for: imguu,imgdd")

    allocate(gdiagud(s%nAtoms,nOrb), gdiagdu(s%nAtoms,nOrb), stat = AllocateStatus)
    if(AllocateStatus /= 0) &
    call abortProgram("[calcMagnetization] Not enough memory for: gdiagdu, gdiagud")

    imguu   = 0.d0
    imgdd   = 0.d0
    gdiagud = cZero
    gdiagdu = cZero

    !$omp parallel default(none) &
    !$omp& private(ix,ep,kp,weight,i,mu,mup,gf,AllocateStatus) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,eta,wght,s,bzs,E_k_imag_mesh,y,gdiagud,gdiagdu,imguu,imgdd)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat=AllocateStatus)
    if(AllocateStatus /= 0) &
    call AbortProgram("[calcMagnetization] Not enough memory for: gf")
    gf = cZero

    if(llineargfsoc .or. llinearsoc) then
      !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         call greenlineargfsoc(s%Ef,ep+eta,s,kp,gf)
         do i=1,s%nAtoms
           do mu=1,nOrb
             mup = mu+nOrb
             gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i) * weight
             gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i) * weight

             imguu(mu,i) = imguu(mu,i) + real(gf(mu ,mu ,i,i)) * weight
             imgdd(mu,i) = imgdd(mu,i) + real(gf(mup,mup,i,i)) * weight
           end do
         end do
      end do
      !$omp end do
    else
      !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         call green(s%Ef,ep+eta,s,kp,gf)
         do i=1,s%nAtoms
           do mu=1,nOrb
             mup = mu+nOrb
             gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i) * weight
             gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i) * weight

             imguu(mu,i) = imguu(mu,i) + real(gf(mu ,mu ,i,i)) * weight
             imgdd(mu,i) = imgdd(mu,i) + real(gf(mup,mup,i,i)) * weight
           end do
         end do
      end do
      !$omp end do
    end if

    deallocate(gf)
    !$omp end parallel
    imguu = imguu / pi
    imgdd = imgdd / pi

    do j=1,s%nAtoms
      mp(:,j)= gdiagdu(j,:) + conjg(gdiagud(j,:))
      mpd(j) = sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9)))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, imguu, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgdd, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp,    ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mpd, s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    mp      = mp/pi
    mpd     = mpd/pi
    mx      = real(mp)
    my      = aimag(mp)
    mxd     = real(mpd)
    myd     = aimag(mpd)

    do i = 1, s%nAtoms
      do mu=1,nOrb
        imguu(mu,i) = 0.5d0 + imguu(mu,i)
        imgdd(mu,i) = 0.5d0 + imgdd(mu,i)
        rho(mu,i) = imguu(mu,i) + imgdd(mu,i)
        mz (mu,i) = imguu(mu,i) - imgdd(mu,i)
      end do
      rhod(i)   = sum(imguu(5:9,i)) + sum(imgdd(5:9,i))
      mzd(i)    = sum(imguu(5:9,i)) - sum(imgdd(5:9,i))
    end do

    deallocate(imguu,imgdd)
    deallocate(gdiagdu, gdiagud)
  end subroutine calcMagnetization

  subroutine calcJacobian(jacobian, N)
    !! Calculated the Jacobian of the spin magnetization
    use mod_f90_kind,      only: double
    use mod_constants,     only: pi, ident_norb2, cZero, pauli_dorb, ident_dorb, cOne
    use mod_parameters,    only: U, offset, eta
    use mod_SOC,           only: llinearsoc, llineargfsoc
    use EnergyIntegration, only: y, wght
    use mod_system,        only: s => sys
    use adaptiveMesh,      only: local_points, E_k_imag_mesh, bzs, activeComm
    use mod_BrillouinZone, only: realBZ
    use TightBinding,      only: nOrb,nOrb2
    use mod_mpi_pars
    implicit none
    integer,                           intent(in)    :: N
    real(double),    dimension(N,N),   intent(inout) :: jacobian
    complex(double), dimension(:,:),     allocatable :: gij,gji,temp,paulitemp
    complex(double), dimension(:,:,:),   allocatable :: temp1,temp2
    complex(double), dimension(:,:,:,:), allocatable :: gf,gvg
    integer :: ix,i,j,mu,nu
    integer :: AllocateStatus
    integer :: sigma,sigmap
    real(double)    :: kp(3), ep
    complex(double) :: weight
    complex(double) :: halfU(s%nAtoms)
    complex(double), dimension(nOrb2, nOrb2, 4) :: pauli_a,pauli_b

    !--------------------- begin MPI vars --------------------
    integer :: ncount2
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount2=N*N

!   Identity and Pauli matrices (excluding d orbitals)
    pauli_a(:,:,1) = ident_norb2(:,:)   ! sigma_0 (includes sp for the last line - total charge neutrality)
    pauli_a(:,:,2) = pauli_dorb(1,:,:)  ! sigma_x
    pauli_a(:,:,3) = pauli_dorb(2,:,:)  ! sigma_y
    pauli_a(:,:,4) = pauli_dorb(3,:,:)  ! sigma_z

    pauli_b(:,:,  1) = ident_dorb  ! excluding d orbitals from charge part
    pauli_b(:,:,2:4) = pauli_a(:,:,2:4)

    ! Prefactor -U/2 in dH/dn and dH/dm
    do i=1,s%nAtoms
      halfU(i) = -0.5d0*U(i+offset)
    end do

    jacobian = 0.d0

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,j,mu,nu,sigma,sigmap,ep,kp,weight,gf,gvg,gij,gji,temp,temp1,temp2,paulitemp) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,s,realBZ,bzs,E_k_imag_mesh,y,eta,wght,halfU,pauli_a,pauli_b,jacobian)
    allocate( temp1(nOrb2, nOrb2, 4), &
              temp2(nOrb2, nOrb2, 4), &
              gij(nOrb2,nOrb2), gji(nOrb2,nOrb2), &
              gf(nOrb2, nOrb2, s%nAtoms, s%nAtoms), &
              temp(nOrb2, nOrb2), paulitemp(nOrb2, nOrb2), stat = AllocateStatus)
    if (AllocateStatus/=0) &
    call abortProgram("[calcJacobian] Not enough memory for: temp1, temp2, gij, gji, gf, temp")
    gf        = cZero
    temp      = cZero

    if(llineargfsoc .or. llinearsoc) then
      allocate(gvg(nOrb2, nOrb2, s%nAtoms, s%nAtoms), STAT = AllocateStatus  )
      if (AllocateStatus/=0) &
      call abortProgram("[calcJacobian] Not enough memory for: gvg")
      gvg = cZero
    end if

   !$omp do schedule(static) reduction(+:jacobian)
   do ix = 1, local_points
      ep = y(E_k_imag_mesh(1,ix))
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
      ! Green function on energy Ef + iy, and wave vector kp
      if(llineargfsoc .or. llinearsoc) then
        call greenlinearsoc(s%Ef,ep+eta,s,kp,gf,gvg)
        gf = gf + gvg
      else
        call green(s%Ef,ep+eta,s,kp,gf)
      end if

      do j=1,s%nAtoms
        do i=1,s%nAtoms

          ! First product: temp1 =  pauli*g_ij
          gij = gf(:,:,i,j)
          do sigma = 1,4
            paulitemp = pauli_a(:,:, sigma)
            call zgemm('n','n',18,18,18,cOne,paulitemp,18,gij,18,cZero,temp,18)
            temp1(:,:,sigma) = temp
          end do

          do nu=5,nOrb
            ! Second product: temp2 = (U_j\delta_{mu,nu} -U_j/2) * sigma * g_ji
            gji = gf(:,:,j,i)
            paulitemp = pauli_b(:,:, 1)
            paulitemp(nu  ,nu  ) = paulitemp(nu  ,nu  ) - 2.d0
            paulitemp(nu+9,nu+9) = paulitemp(nu+9,nu+9) - 2.d0
            call zgemm('n','n',18,18,18,halfU(j),paulitemp,18,gji,18,cZero,temp,18)
            temp2(:,:,1) = temp

            ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

            ! Charge density-charge density part
            gij = temp1(:,:,1)
            gji = temp2(:,:,1)
            call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

            do mu = 1,nOrb
              ! Last line (Total charge neutrality)
              jacobian(  8*s%nAtoms+1,(j-1)*8+(nu-4)) = jacobian(  8*s%nAtoms+1,(j-1)*8+(nu-4)) + real(temp(mu,mu) + temp(mu+9,mu+9))
              if(mu<5) cycle
              jacobian((i-1)*8+(mu-4),(j-1)*8+(nu-4)) = jacobian((i-1)*8+(mu-4),(j-1)*8+(nu-4)) + real(temp(mu,mu) + temp(mu+9,mu+9))
            end do

            ! Magnetic density-charge density part
            do sigma = 2,4
              gij = temp1(:,:,sigma)
              gji = temp2(:,:,1)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

              do mu = 5,nOrb
                jacobian((i-1)*8+4+sigma,(j-1)*8+(nu-4)) = jacobian((i-1)*8+4+sigma,(j-1)*8+(nu-4)) + real(temp(mu,mu) + temp(mu+9,mu+9))
              end do
            end do

          end do

          gji = gf(:,:,j,i)

          do sigmap = 2,4
            ! Second product: temp2 = (-U/2) * sigma* g_ji
            paulitemp = pauli_b(:,:,sigmap)
            call zgemm('n','n',18,18,18,halfU(j),paulitemp,18,gji,18,cZero,temp,18)
            temp2(:,:,sigmap) = temp
          end do

          ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

          ! Charge density-magnetic density part
          do sigmap = 2,4
            gij = temp1(:,:,1)
            gji = temp2(:,:,sigmap)
            call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

            do mu = 1,nOrb
              ! Last line (Total charge neutrality)
              jacobian(  8*s%nAtoms+1,(j-1)*8+4+sigmap) = jacobian(  8*s%nAtoms+1,(j-1)*8+4+sigmap) + real(temp(mu,mu) + temp(mu+9,mu+9))
              if(mu<5) cycle
              jacobian((i-1)*8+(mu-4),(j-1)*8+4+sigmap) = jacobian((i-1)*8+(mu-4),(j-1)*8+4+sigmap) + real(temp(mu,mu) + temp(mu+9,mu+9))
            end do
          end do

          ! Magnetic density-magnetic density part
          do sigma = 2,4
            do sigmap = 2,4
              gij = temp1(:,:,sigma)
              gji = temp2(:,:,sigmap)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

              do mu = 5,nOrb
                jacobian((i-1)*8+4+sigma,(j-1)*8+4+sigmap) = jacobian((i-1)*8+4+sigma,(j-1)*8+4+sigmap) + real(temp(mu,mu) + temp(mu+9,mu+9))
              end do
            end do
          end do

          ! **** Removing non-linear (quadratic) terms: ****
          if(llineargfsoc .or. llinearsoc) then
            gij = gvg(:,:,i,j)

            do sigma = 1,4
              ! First product: temp1 =  pauli*g_ij
              paulitemp = pauli_a(:,:, sigma)
              call zgemm('n','n',18,18,18,cOne,paulitemp,18,gij,18,cZero,temp,18)
              temp1(:,:, sigma) = temp
            end do

            do nu=5,nOrb
              gji = gvg(:,:,j,i)

              ! Second product: temp2 = (-U/2) * sigma* g_ji
              paulitemp = pauli_b(:,:, 1)
              paulitemp(nu  ,nu  ) = paulitemp(nu  ,nu  ) - 2.d0
              paulitemp(nu+9,nu+9) = paulitemp(nu+9,nu+9) - 2.d0
              call zgemm('n','n',18,18,18,halfU(j),paulitemp,18,gji,18,cZero,temp,18)
              temp2(:,:, 1) = temp

              ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

              ! Charge density-charge density part
              gij = temp1(:,:, 1)
              gji = temp2(:,:, 1)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

              do mu = 1,nOrb
                ! Last line (Total charge neutrality)
                jacobian(  8*s%nAtoms+1,(j-1)*8+(nu-4)) = jacobian(  8*s%nAtoms+1,(j-1)*8+(nu-4)) - real(temp(mu,mu) + temp(mu+9,mu+9))
                if(mu<5) cycle
                jacobian((i-1)*8+(mu-4),(j-1)*8+(nu-4)) = jacobian((i-1)*8+(mu-4),(j-1)*8+(nu-4)) - real(temp(mu,mu) + temp(mu+9,mu+9))
              end do

              ! Magnetic density-charge density part
              do sigma = 2,4
                gij = temp1(:,:, sigma)
                gji = temp2(:,:, 1)
                call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

                do mu = 5,nOrb
                  jacobian((i-1)*8+4+sigma,(j-1)*8+(nu-4)) = jacobian((i-1)*8+4+sigma,(j-1)*8+(nu-4)) - real(temp(mu,mu) + temp(mu+9,mu+9))
                end do
              end do

            end do

            gji = gvg(:,:,j,i)

            do sigmap = 2,4
              ! Second product: temp2 = (-U/2) * sigma* g_ji
              paulitemp = pauli_b(:,:, sigmap)
              call zgemm('n','n',18,18,18,halfU(j),paulitemp,18,gji,18,cZero,temp,18)
              temp2(:,:, sigmap) = temp
            end do

            ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

            ! Charge density-magnetic density part
            do sigmap = 2,4
              gij = temp1(:,:, 1)
              gji = temp2(:,:, sigmap)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

              do mu = 1,nOrb
                ! Last line (Total charge neutrality)
                jacobian(  8*s%nAtoms+1,(j-1)*8+4+sigmap) = jacobian(  8*s%nAtoms+1,(j-1)*8+4+sigmap) - real(temp(mu,mu) + temp(mu+9,mu+9))
                if(mu<5) cycle
                jacobian((i-1)*8+(mu-4),(j-1)*8+4+sigmap) = jacobian((i-1)*8+(mu-4),(j-1)*8+4+sigmap) - real(temp(mu,mu) + temp(mu+9,mu+9))
              end do
            end do

            ! Magnetic density-magnetic density part
            do sigma = 2,4
              do sigmap = 2,4
                gij = temp1(:,:, sigma)
                gji = temp2(:,:, sigmap)
                call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

                do mu = 5,nOrb
                  jacobian((i-1)*8+4+sigma,(j-1)*8+4+sigmap) = jacobian((i-1)*8+4+sigma,(j-1)*8+4+sigmap) - real(temp(mu,mu) + temp(mu+9,mu+9))
                end do
              end do
            end do
          end if ! End linear part


        end do ! End nAtoms i loop
      end do ! End nAtoms j loop
    end do ! End Energy+nkpt loop
    !$omp end do nowait

    ! Last column: Derivatives with respect to the Fermi level
    ! No energy integral - only BZ integration of G(Ef)
    !$omp do schedule(static) reduction(+:jacobian)
    do ix = 1, realBZ%workload
      kp = realBZ%kp(1:3,ix)
      weight = cmplx(1.d0,0.d0) * realBZ%w(ix)

      ! Green function at Ef + ieta, and wave vector kp
      if(llineargfsoc .or. llinearsoc) then
        call greenlinearsoc(s%Ef,eta,s,kp,gf,gvg)
        gf = gf + gvg
      else
        call green(s%Ef,eta,s,kp,gf)
      end if

      do i=1,s%nAtoms
        gij = gf(:,:,i,i)

        ! Charge density per orbital lines
        ! temp1 =  pauli*g_ii
        paulitemp = pauli_a(:,:, 1)
        call zgemm('n','n',18,18,18,weight,paulitemp,18,gij,18,cZero,temp,18)

        do mu = 1,nOrb
          ! Last line (Total charge neutrality)
          jacobian(  8*s%nAtoms+1,8*s%nAtoms+1) = jacobian(  8*s%nAtoms+1,8*s%nAtoms+1) - aimag(temp(mu,mu) + temp(mu+9,mu+9))
          if(mu<5) cycle
          jacobian((i-1)*8+(mu-4),8*s%nAtoms+1) = jacobian((i-1)*8+(mu-4),8*s%nAtoms+1) - aimag(temp(mu,mu) + temp(mu+9,mu+9))
        end do

        do sigma = 2,4
          ! temp1 =  pauli*g_ii
          paulitemp = pauli_a(:,:, sigma)
          call zgemm('n','n',18,18,18,weight,paulitemp,18,gij,18,cZero,temp,18)

          do mu = 5,nOrb
            jacobian((i-1)*8+4+sigma,8*s%nAtoms+1) = jacobian((i-1)*8+4+sigma,8*s%nAtoms+1) - aimag(temp(mu,mu) + temp(mu+9,mu+9))
          end do
        end do

        ! No linear correction is needed since it's a single Green function

      end do ! End nAtoms i loop
    end do ! End nkpt loop
    !$omp end do nowait

    deallocate(paulitemp)
    deallocate(temp1, temp2, temp)
    deallocate(gf, gij, gji)
    if(allocated(gvg)) deallocate(gvg)
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, jacobian, ncount2, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)

    jacobian = jacobian/pi
    do i = 1, 8*s%nAtoms
      jacobian(i,i) = jacobian(i,i) - 1.d0
    end do
  end subroutine calcJacobian

  subroutine calcLGS()
    !! Calculates the ground state charge, magnetization and orbital angular momentum ground state
    use mod_f90_kind,      only: double
    use mod_constants,     only: cZero,pi,rad2deg
    use mod_System,        only: s => sys
    use TightBinding,      only: nOrb,nOrb2
    use mod_parameters,    only: output, eta
    use EnergyIntegration, only: y, wght
    use adaptiveMesh
    use mod_magnet
    use mod_mpi_pars
    implicit none
    integer*8    :: ix
    integer      :: AllocateStatus
    integer      :: i,mu,nu,mup,nup
    real(double) :: kp(3)
    real(double) :: weight, ep
    complex(double), dimension(:,:,:,:), allocatable :: gf
    complex(double), dimension(:,:,:),   allocatable :: gupgd
    !--------------------- begin MPI vars --------------------
    integer :: ncount
    ncount=s%nAtoms*nOrb*nOrb
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    if(allocated(lxm)) deallocate(lxm)
    if(allocated(lym)) deallocate(lym)
    if(allocated(lzm)) deallocate(lzm)
    if(allocated(lxpm)) deallocate(lxpm)
    if(allocated(lypm)) deallocate(lypm)
    if(allocated(lzpm)) deallocate(lzpm)
    allocate( lxm(s%nAtoms), &
              lym(s%nAtoms), &
              lzm(s%nAtoms), &
              lxpm(s%nAtoms), &
              lypm(s%nAtoms), &
              lzpm(s%nAtoms), stat = AllocateStatus )
    if (AllocateStatus/=0) &
    call abortProgram("[calcLGS] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")

    allocate(gupgd(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) &
    call abortProgram("[calcLGS] Not enough memory for: gupgd")

    if(rField == 0) &
    write(output%unit_loop,"('[calcLGS] Calculating Orbital Angular Momentum ground state... ')")

    ! Calculating the jacobian using a complex integral
    gupgd  = cZero
    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(local_points,s,E_k_imag_mesh,bzs,eta,y,wght,gupgd)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat = AllocateStatus)
    if (AllocateStatus/=0) &
    call abortProgram("[calcLGS] Not enough memory for: gf")

    gf = cZero
    !$omp do schedule(static) reduction(+:gupgd)
    do ix = 1, local_points
        kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
        ep = y(E_k_imag_mesh(1,ix))
        weight = bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix)) * wght(E_k_imag_mesh(1,ix))
        !Green function on energy Ef + iy, and wave vector kp
        call green(s%Ef,ep+eta,s,kp,gf)

        do i=1,s%nAtoms
          do mu=1,nOrb
            mup = mu+nOrb
            do nu=1,nOrb
              nup = nu+nOrb
              gupgd(mu,nu,i) = gupgd(mu,nu,i) + (gf(mu,nu,i,i) + gf(mup,nup,i,i)) * weight
            end do
          end do
        end do
      end do
    !$omp end do

    deallocate(gf)
    !$omp end parallel

    gupgd = gupgd / pi
    call MPI_Allreduce(MPI_IN_PLACE, gupgd, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    lxpm = 0.d0
    lypm = 0.d0
    lzpm = 0.d0
    lxm  = 0.d0
    lym  = 0.d0
    lzm  = 0.d0

    do nu=5,9
      do mu=5,9
        do i=1,s%nAtoms
          lxpm(i) = lxpm(i) + real(lxp(mu,nu,i)*gupgd(nu,mu,i))
          lypm(i) = lypm(i) + real(lyp(mu,nu,i)*gupgd(nu,mu,i))
          lzpm(i) = lzpm(i) + real(lzp(mu,nu,i)*gupgd(nu,mu,i))
          lxm (i) = lxm (i) + real(lx (mu,nu  )*gupgd(nu,mu,i))
          lym (i) = lym (i) + real(ly (mu,nu  )*gupgd(nu,mu,i))
          lzm (i) = lzm (i) + real(lz (mu,nu  )*gupgd(nu,mu,i))
        end do
      end do
    end do

    ! Calculating angles of GS OAM (in units of pi)
    do i = 1,s%nAtoms
      labs(i)   = sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
      if(abs(labs(i))>1.d-8) then
        ltheta(i) = acos(lzm(i)/labs(i))*rad2deg
      else
        ltheta(i) = 0.d0
      end if
      if(abs(ltheta(i))>1.d-8) then
        if(abs(abs(ltheta(i))-180.d0)>1.d-8) then
          lphi(i)   = atan2(lym(i),lxm(i))*rad2deg
        else
          lphi(i) = 0.d0
        end if
      else
        lphi(i) = 0.d0
      end if
      lpabs(i)  = sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2))
      if(abs(lpabs(i))>1.d-8) then
        lptheta(i)= acos(lzpm(i)/lpabs(i))*rad2deg
      else
        lptheta(i) = 0.d0
      end if
      if(abs(lptheta(i))>1.d-8) then
        if(abs(abs(lptheta(i))-180.d0)>1.d-8) then
          lpphi(i)   = atan2(lypm(i),lxpm(i))*rad2deg
        else
          lpphi(i) = 0.d0
        end if
      else
        lpphi(i) = 0.d0
      end if
    end do

    deallocate(gupgd)
  end subroutine calcLGS

  subroutine rotate_magnetization_to_field()
  !! Rotate the magnetization to the direction of the field (useful for SOC=F)
    use mod_f90_kind,   only: double
    use mod_constants,  only: deg2rad
    use mod_magnet,     only: hw_count,hw_list,hhw,mx,my,mz,mp
    use mod_parameters, only: output
    use mod_System,     only: s => sys
    use TightBinding,   only: nOrb
    use mod_mpi_pars,   only: rField
    implicit none
    integer      :: i,j,sign
    real(double) :: mdotb,mabs(nOrb,s%nAtoms)

    if(rField == 0) &
    write(output%unit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    do i = 1, s%nAtoms
      do j=1,nOrb
        mdotb   = hhw(1,i)*mx(j,i)+hhw(2,i)*my(j,i)+hhw(3,i)*mz(j,i)
        sign    = dble(mdotb/abs(mdotb))
        mabs(j,i) = sqrt((mx(j,i)**2)+(my(j,i)**2)+(mz(j,i)**2))
        mx(j,i)   = sign*mabs(j,i)*sin(hw_list(hw_count,2)*deg2rad)*cos(hw_list(hw_count,3)*deg2rad)
        my(j,i)   = sign*mabs(j,i)*sin(hw_list(hw_count,2)*deg2rad)*sin(hw_list(hw_count,3)*deg2rad)
        mz(j,i)   = sign*mabs(j,i)*cos(hw_list(hw_count,2)*deg2rad)
        mp(j,i)   = cmplx(mx(j,i),my(j,i),double)
      end do
    end do

    ! Writing new n and rotated mag to file (without self-consistency)
    if(rField == 0) call write_sc_results()

    ! Writing self-consistency results on screen
    if(rField == 0) call print_sc_results()
  end subroutine rotate_magnetization_to_field

  ! Writes the self-consistency results on the screen
  subroutine print_sc_results()
    use mod_parameters, only: output
    use mod_system,     only: s => sys
    use mod_SOC,        only: SOC
    use mod_magnet,     only: rho, mvec_cartesian, mp, mvec_spherical, &
                              lxm, lym, lzm, ltheta, lphi, labs
    implicit none
    integer :: i

    write(output%unit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    write(output%unit_loop,"(28x,'Ef=',f11.8)") s%Ef
    write(output%unit_loop,"(11x,' *************** Charge density: ****************')")
    do i=1,s%nAtoms
      write(output%unit_loop,"(4x,'Ns(',i2.0,')=',f11.8,4x,'Np(',i2.0,')=',f11.8,4x,'Nd(',i2.0,')=',f11.8)") i, rho(1,i),i, sum(rho(2:4,i)),i, sum(rho(5:9,i))
    end do
    write(output%unit_loop,"(11x,' *********** Magnetization components: **********')")
    do i=1,s%nAtoms
      write(output%unit_loop,"(4x,'Mx(',i2.0,')=',f11.8,4x,'My(',i2.0,')=',f11.8,4x,'Mz(',i2.0,')=',f11.8)") i,mvec_cartesian(1,i),i,mvec_cartesian(2,i),i,mvec_cartesian(3,i)
      if(abs(sum(mp(:,i)))/=0) &
      write(output%unit_loop,"(12x,'theta = ',f11.6,'  ',4x,'phi = ',f11.6)") mvec_spherical(2,i),mvec_spherical(3,i)
    end do
    if((lGSL).or.(SOC)) then
      write(output%unit_loop,"(11x,' ****** Orbital components in global frame: *****')")
      do i=1,s%nAtoms
        write(output%unit_loop,"(4x,'Lx(',i2.0,')=',f11.8,4x,'Ly(',i2.0,')=',f11.8,4x,'Lz(',i2.0,')=',f11.8)") i,lxm(i),i,lym(i),i,lzm(i)
        if(sqrt(lxm(i)**2+lym(i)**2)/=0) &
        write(output%unit_loop,"(12x,'theta = ',f11.6,'  ',4x,'phi = ',f11.6)") ltheta(i),lphi(i)
      end do
      ! write(output%unit_loop,"(11x,' *** Orbital components in local frame:  ***')")
      ! do i=1,s%nAtoms
      !   write(output%unit_loop,"(4x,'Lxp(',i2.0,')=',f11.8,4x,'Lyp(',i2.0,')=',f11.8,4x,'Lzp(',i2.0,')=',f11.8)") i,lxpm(i),i,lypm(i),i,lzpm(i)
      !   if(sqrt(lxpm(i)**2+lypm(i)**2)/=0) &
      !   write(output%unit_loop,"(12x,'theta = ',f11.6,'  ',4x,'phi = ',f11.6)") lptheta(i),lpphi(i)
      ! end do
      write(output%unit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(output%unit_loop,"(4x,' N(',i2.0,')=',f11.8,4x,' M(',i2.0,')=',f11.8,4x,' L(',i2.0,')=',f11.8)") i,sum(rho(:,i)),i,mvec_spherical(1,i),i,labs(i)
      end do
    else
      write(output%unit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(output%unit_loop,"(9x,' N(',i2.0,') =',f11.8,4x,' M(',i2.0,') =',f11.8)") i,sum(rho(:,i)),i,mvec_spherical(1,i)
      end do
    end if
    write(output%unit_loop,"('|----------=============================================================----------|')")
  end subroutine print_sc_results


  subroutine write_sc_results()
    !! Writes the self-consistency results into files and broadcasts the scfile for the next iteration.
    use mod_parameters,    only: output, dfttype
    use EnergyIntegration, only: parts
    use mod_magnet,        only: rho, mx, my, mz
    use mod_system,        only: s => sys
    use TightBinding,      only: nOrb
    use mod_mpi_pars
    implicit none
    character(len=30) :: formatvar
    integer           :: i,j

    if(rField == 0) then
      ! Writing new results (mx, my, mz and n) to file
      write(output%unit_loop,"('[write_sc_results] Writing new n, mx, my, mz and Ef to file...')")
      write(scfile,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,'.dat')") output%SOCchar, trim(output%Sites),dfttype,parts,trim(output%BField),trim(output%info),trim(output%SOC)
      open (unit=99,status='replace',file=scfile)

      write(formatvar,fmt="(a,i0,a)") '(',4*nOrb,'(es21.11,2x))'

      do i=1,s%nAtoms
        write(99,fmt=formatvar) (rho(j,i), j=1,nOrb),(mx(j,i), j=1,nOrb),(my(j,i), j=1,nOrb),(mz(j,i), j=1,nOrb)
      end do
      write(99,"(es21.11,2x,'! Ef  ')") s%Ef
      write(99,"('! n(1:nOrb), mx(1:nOrb), my(1:nOrb), mz(1:nOrb) per site ')")
      write(99,"('! Ef ')")

      close(99)
    end if

    call MPI_Bcast(scfile, len(scfile), MPI_CHARACTER, 0, FieldComm,ierr)
  end subroutine write_sc_results


  ! Writes the initial values for the self-consistency
  subroutine print_sc_step(n,mx,my,mz,Ef,fvec)
    use mod_f90_kind,   only: double
    use mod_parameters, only: output
    use mod_system,     only: s => sys
    use mod_magnet,     only: iter
    use TightBinding,   only: nOrb
    use mod_mpi_pars
    implicit none
    real(double),dimension(neq), intent(in), optional :: fvec
    real(double),dimension(s%nAtoms),     intent(in) :: n,mx,my,mz
    real(double)                    ,     intent(in) :: Ef
    integer :: i,mu
    real(double) :: fvecsum

    if(rField==0) then
      if(present(fvec)) then
        do i=1,s%nAtoms
          fvecsum = 0.0
          do mu = 5,nOrb
            fvecsum = fvecsum + fvec((i-1)*8+(mu-4))
          end do

          if(abs(cmplx(mx(i),my(i)))>1.d-10) then
            write(output%unit_loop,"('Plane ',I2,': N(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,n(i),i,mx(i),i,my(i),i,mz(i)
            write(output%unit_loop,"(15x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvecsum,(i-1)*8+6,fvec((i-1)*8+6),(i-1)*8+7,fvec((i-1)*8+7),(i-1)*8+8,fvec((i-1)*8+8)
          else
            write(output%unit_loop,"('Plane ',I2,': N(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,n(i),i,mz(i)
            write(output%unit_loop,"(15x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvecsum,(i-1)*8+8,fvec((i-1)*8+8)
          end if
        end do
        write(output%unit_loop,"(13x,'Ef=',es16.9)") Ef
        write(output%unit_loop,"(15x,'fvec(',I2,')=',es16.9)") 8*s%nAtoms+1,fvec(8*s%nAtoms+1)
      else if(iter == 1) then
        write(output%unit_loop,"('|---------------- Starting charge density, magnetization and Ef ----------------|')")
        do i=1,s%nAtoms
          if(abs(cmplx(mx(i),my(i)))>1.d-10) then
            write(output%unit_loop,"('Plane ',I2,': N(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,n(i),i,mx(i),i,my(i),i,mz(i)
          else
            write(output%unit_loop,"('Plane ',I2,': N(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,n(i),i,mz(i)
          end if
        end do
        write(output%unit_loop,"(13x,'Ef=',es16.9)") Ef
      end if
    end if
  end subroutine print_sc_step


  ! This subroutine calculates the self-consistency equations
  !     n  - rho_in    = 0
  !     mx - mx_in   = 0
  !     my - my_in   = 0
  !     mz - mz_in   = 0
  !  sum n - n_total = 0
  ! and the correspondent jacobian
#if !defined(_OSX) && !defined(_JUQUEEN)
  subroutine sc_equations_and_jacobian(N,x,fvec,selfconjac,iuser,ruser,iflag)
    use mod_f90_kind,   only: double
    use mod_system,     only: s => sys
    use TightBinding,   only: nOrb
    use mod_parameters, only: lcheckjac
    use mod_magnet,     only: iter,rho,rhod,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix,    only: update_Umatrix
    use mod_tools,      only: itos
    use mod_mpi_pars
    implicit none
    integer  :: N,i,mu,iflag
    integer     ,  intent(inout)             :: iuser(*)
    real(double),  intent(inout)             :: ruser(*)
    real(double),   dimension(N)             :: x,fvec
    real(double),   dimension(N,N)           :: selfconjac
    real(double),   dimension(nOrb,s%nAtoms) :: rho_in
    real(double),   dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms)      :: mpd_in

    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*8+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*8+6)
      myd_in(i) = x((i-1)*8+7)
      mzd_in(i) = x((i-1)*8+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i))
    end do
    s%Ef    = x(8*s%nAtoms+1)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    select case (iflag)
    case(1)
      call calcMagnetization()
      do i = 1, s%nAtoms
        do mu = 5,nOrb
          fvec((i-1)*8+(mu-4)) = rho(mu,i) - rho_in(mu,i)
        end do
        fvec((i-1)*8+6) =  mxd(i) -  mxd_in(i)
        fvec((i-1)*8+7) =  myd(i) -  myd_in(i)
        fvec((i-1)*8+8) =  mzd(i) -  mzd_in(i)
      end do
      fvec(8*s%nAtoms+1) = sum(rho) - s%totalOccupation

      call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

      if(lontheflysc) call write_sc_results()
    case(2)
      if(lcheckjac) call check_jacobian(neq,x)

      call calcJacobian(selfconjac, N)
    case default
      call abortProgram("[sc_equations_and_jacobian] Problem in self-consistency! iflag = " // trim(itos(iflag)))
    end select

    iter = iter + 1
  end subroutine sc_equations_and_jacobian

  ! This subroutine calculates the self-consistency equations
  !     n  - rho_in    = 0
  !     mx - mx_in   = 0
  !     my - my_in   = 0
  !     mz - mz_in   = 0
  !  sum n - n_total = 0
  subroutine sc_equations(N,x,fvec,iuser,ruser,iflag)
    use mod_f90_kind, only: double
    use mod_constants, only: cI
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter,rho,rhod,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix, only: update_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,mu,iflag
    integer     ,   intent(inout)            :: iuser(*)
    real(double),   intent(inout)            :: ruser(*)
    real(double),   dimension(N)             :: x,fvec
    real(double),   dimension(nOrb,s%nAtoms) :: rho_in
    real(double),   dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms)      :: mpd_in

    iflag=0
    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*8+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*8+6)
      myd_in(i) = x((i-1)*8+7)
      mzd_in(i) = x((i-1)*8+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i))
    end do
    s%Ef    = x(8*s%nAtoms+1)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    call calcMagnetization()
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        fvec((i-1)*8+(mu-4)) = rho(mu,i) - rho_in(mu,i)
      end do
      fvec((i-1)*8+6) =  mxd(i) -  mxd_in(i)
      fvec((i-1)*8+7) =  myd(i) -  myd_in(i)
      fvec((i-1)*8+8) =  mzd(i) -  mzd_in(i)
    end do
    fvec(8*s%nAtoms+1) = sum(rho) - s%totalOccupation

    call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

    if(lontheflysc) call write_sc_results()

    iter = iter + 1
  end subroutine sc_equations

#endif

  ! This subroutine calculates the self-consistency equations
  !     n  - rho_in    = 0
  !     mx - mx_in   = 0
  !     my - my_in   = 0
  !     mz - mz_in   = 0
  !  sum n - n_total = 0
  ! and the correspondent jacobian
  subroutine sc_eqs_and_jac_old(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_f90_kind,   only: double
    use mod_system,     only: s => sys
    use TightBinding,   only: nOrb
    use mod_parameters, only: lcheckjac
    use mod_magnet,     only: iter,rho,rhod,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix,    only: update_Umatrix
    use mod_tools,      only: itos
    use mod_mpi_pars
    implicit none
    integer  :: N,i,mu,iflag,ldfjac
    real(double),   dimension(N)             :: x,fvec
    real(double),   dimension(ldfjac,N)      :: selfconjac
    real(double),   dimension(nOrb,s%nAtoms) :: rho_in
    real(double),   dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms)      :: mpd_in

    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*8+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*8+6)
      myd_in(i) = x((i-1)*8+7)
      mzd_in(i) = x((i-1)*8+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i))
    end do
    s%Ef    = x(8*s%nAtoms+1)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    flag: select case (iflag)
    case(1)
      call calcMagnetization()
      do i = 1, s%nAtoms
        do mu = 5,nOrb
          fvec((i-1)*8+(mu-4)) = rho(mu,i) - rho_in(mu,i)
        end do
        fvec((i-1)*8+6) =  mxd(i) -  mxd_in(i)
        fvec((i-1)*8+7) =  myd(i) -  myd_in(i)
        fvec((i-1)*8+8) =  mzd(i) -  mzd_in(i)
      end do
      fvec(8*s%nAtoms+1) = sum(rho) - s%totalOccupation

      call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

      if(lontheflysc) call write_sc_results()
    case(2)
      if(lcheckjac) call check_jacobian(neq,x)

      call calcJacobian(selfconjac, N)
    case default
      call abortProgram("[sc_eqs_and_jac_old] Problem in self-consistency! iflag = " // trim(itos(iflag)))
    end select flag

    iter = iter + 1
  end subroutine sc_eqs_and_jac_old

  ! This subroutine calculates the self-consistency equations
  !     n  - rho_in    = 0
  !     mx - mx_in   = 0
  !     my - my_in   = 0
  !     mz - mz_in   = 0
  !  sum n - n_total = 0
  subroutine sc_eqs_old(N,x,fvec,iflag)
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter,rho,rhod,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix, only: update_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,mu,iflag
    real(double),   dimension(N)             :: x,fvec
    real(double),   dimension(nOrb,s%nAtoms) :: rho_in
    real(double),   dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms)      :: mpd_in

    iflag=0
    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*8+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*8+6)
      myd_in(i) = x((i-1)*8+7)
      mzd_in(i) = x((i-1)*8+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i))
    end do
    s%Ef    = x(8*s%nAtoms+1)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    call calcMagnetization()
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        fvec((i-1)*8+(mu-4)) = rho(mu,i) - rho_in(mu,i)
      end do
      fvec((i-1)*8+6) =  mxd(i) -  mxd_in(i)
      fvec((i-1)*8+7) =  myd(i) -  myd_in(i)
      fvec((i-1)*8+8) =  mzd(i) -  mzd_in(i)
    end do
    fvec(8*s%nAtoms+1) = sum(rho) - s%totalOccupation

    call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

    if(lontheflysc) call write_sc_results()

    iter = iter + 1
  end subroutine sc_eqs_old

  ! This subroutine calculates the jacobian of the system of equations
  !     n  - rho_in    = 0
  !     mx - mx_in   = 0
  !     my - my_in   = 0
  !     mz - mz_in   = 0
  !  sum n - n_total = 0
  subroutine sc_jac_old(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter,rhod0,rho0,rho
    use mod_Umatrix, only: update_Umatrix
    implicit none
    integer       :: N,i,mu,ldfjac,iflag
    real(double)  :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(double),   dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    real(double),   dimension(nOrb,s%nAtoms) :: rho_in
    complex(double),dimension(s%nAtoms)      :: mpd_in
    !--------------------- begin MPI vars --------------------

    iflag=0
    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*8+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*8+6)
      myd_in(i) = x((i-1)*8+7)
      mzd_in(i) = x((i-1)*8+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i))
    end do
    s%Ef    = x(8*s%nAtoms+1)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)

    fvec=fvec

    call calcJacobian(selfconjac, N)

    iter = iter + 1
  end subroutine sc_jac_old

end module mod_self_consistency
