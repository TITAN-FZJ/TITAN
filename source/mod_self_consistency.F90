module mod_self_consistency
   use mod_f90_kind, only: double
   implicit none
   character(len=300)  :: default_file
   real(double) :: mag_tol = 1.d-10

   character(len=200) :: scfile = ""
   !! Give a file to start self-consistency
   logical :: skipsc
   !! Skip self-consistency
   logical :: lselfcon = .false.
   logical :: lGSL = .false.
   logical :: lontheflysc = .false.
   logical :: lslatec = .false.
   logical :: lnojac = .false.
   logical :: lrotatemag = .false.

contains

  subroutine doSelfConsistency()
    use mod_magnet, only: lp_matrix, mtheta, mphi
    use adaptiveMesh, only: genLocalEKMesh, freeLocalEKMesh
    use mod_mpi_pars, only: rField, sField, FieldComm
    use mod_SOC, only: SOC
    implicit none
    logical :: lsuccess = .false.

    ! Distribute Energy Integration across all points available
    call genLocalEKMesh(rField,sField, FieldComm)

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

    return
  end subroutine doSelfConsistency

  ! Tries to read n and m if available
  subroutine read_previous_results(lsuccess)
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_magnet, only: mx, my, mz, mxd, myd, mzd, mpd, &
                          hw_count, hw_list, lfield, rho, rhod
    use mod_parameters, only: output, magaxis, magaxisvec, offset, layertype
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_mpi_pars, only: abortProgram, rField
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
        magaxisvec = [cos(magaxisvec(2)*pi/180)*sin(magaxisvec(1)*pi/180), sin(magaxisvec(2)*pi/180)*sin(magaxisvec(1)*pi/180), cos(magaxisvec(1)*pi/180)]
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
        mxd = mzd * sin(hw_list(hw_count,2)*pi) * cos(hw_list(hw_count,3)*pi)
        myd = mzd * sin(hw_list(hw_count,2)*pi) * sin(hw_list(hw_count,3)*pi)
        mzd = mzd * cos(hw_list(hw_count,2)*pi)
      end if

      mpd = cmplx(mxd,myd)

      mx = 0.d0
      my = 0.d0
      mz = 0.d0
      do i = 1, s%nAtoms
        mx(5:9,i) = 0.2d0*mxd(i)
        my(5:9,i) = 0.2d0*myd(i)
        mz(5:9,i) = 0.2d0*mzd(i)

        rho(  1,i) = s%Types(s%Basis(i)%Material)%OccupationS
        rho(2:4,i) = s%Types(s%Basis(i)%Material)%OccupationP/3.d0
        rho(5:9,i) = s%Types(s%Basis(i)%Material)%OccupationD/5.d0
        rhod(i)    = s%Types(s%Basis(i)%Material)%OccupationD
      end do
    end if

    call init_Umatrix(mzd,mpd,rhod,s%nAtoms,nOrb)

    return
  end subroutine read_previous_results

  ! This subroutine reads previous band-shifting and magnetization results
  subroutine read_sc_results(err,lsuccess)
    use mod_f90_kind, only: double
    use mod_parameters, only: output, dfttype
    use EnergyIntegration, only: parts
    use mod_magnet, only: rho, mp, mx, my, mz, rhod, mpd, mxd, myd, mzd, hw_count
    use TightBinding, only: nOrb
    use mod_system, only: s => sys
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
      if(((hw_count)==1)) then !.and.(Npl==Npl_i)) then ! Filename in inputcard (1st iteration on loop)
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
    return
  end subroutine read_sc_results

  subroutine calcMagAngle()
    use mod_constants, only: pi
    use mod_magnet, only: mx, my, mz, mabs, &
                          mtheta, mphi, mvec_cartesian, mvec_spherical
    use mod_system, only: s => sys
    use mod_susceptibilities, only: lrot
    implicit none
    integer :: i

    ! Calculating new angles of GS magnetization in units of pi and magnetization vector
    do i = 1,s%nAtoms
      mabs(i)   = sqrt((sum(mx(:,i))**2)+(sum(my(:,i))**2)+(sum(mz(:,i))**2))
      mtheta(i) = acos(sum(mz(:,i))/mabs(i))/pi
      if(abs(mtheta(i))>1.d-8) then
        if(abs(abs(mtheta(i))-1.d0)>1.d-8) then
          mphi(i)   = atan2(sum(my(:,i)),sum(mx(:,i)))/pi
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
    return
  end  subroutine calcMagAngle

  subroutine calcMagneticSelfConsistency()
  !! This subroutine performs the self-consistency
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_parameters, only: output
    use mod_magnet, only: iter, rhod, mxd, myd, mzd
    use mod_mpi_pars, only: calcWorkload, rField
    use adaptiveMesh
    use mod_system, only: s => sys
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
    integer                       :: neq,maxfev,ml,mr,mode,nfev,njev,lwa,ifail=0

    neq = 4*s%nAtoms+1
    allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )

    ! Putting read n and m existing solutions into sc_solu (first guess of the subroutine)
    sc_solu(1:s%nAtoms)              = rhod
    sc_solu(s%nAtoms+1:2*s%nAtoms)   = mxd
    sc_solu(2*s%nAtoms+1:3*s%nAtoms) = myd
    sc_solu(3*s%nAtoms+1:4*s%nAtoms) = mzd
    sc_solu(neq) = s%Ef
    iter  = 1

    if(rField == 0) write(output%unit_loop,"('[self_consistency] Starting self-consistency:')")

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
        factor = 100.d0
        call c05ncf(sc_eqs_old,neq,sc_solu,fvec,mag_tol,maxfev,ml,mr,epsfcn,diag,mode,factor,0,nfev,jac,neq,wa,lwa,qtf,w,ifail)
      else
!         call c05pbf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,neq,mag_tol,wa,lwa,ifail)
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
!         diag(Npl+1:4*Npl) = 100.d0
        factor = 100.d0
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
        factor = 100.d0
        call c05qcf(sc_equations,neq,sc_solu,fvec,mag_tol,maxfev,ml,mr,epsfcn,mode,diag,factor,0,nfev,jac,wa,qtf,iuser,ruser,ifail)
      else
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
!         diag(npl+1:4*npl) = 100.d0
        factor = 100.d0
        call c05rcf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,mag_tol,maxfev,mode,diag,factor,0,nfev,njev,wa,qtf,iuser,ruser,ifail)
      end if
    end if
#endif

    deallocate(sc_solu,diag,qtf,fvec,jac,wa)

    call calcMagAngle()

    return
  end subroutine calcMagneticSelfConsistency

  subroutine calcMagnetization()
    !! Calculates occupation density and magnetization.
    use mod_f90_kind, only: double
    use mod_constants, only: cI, pi, cZero
    use mod_SOC, only: llinearsoc, llineargfsoc
    use EnergyIntegration, only: y, wght
    use mod_system, only: s => sys
    use mod_magnet, only: mx, my, mz, mp, rho, mxd, myd, mzd, mpd, rhod
    use adaptiveMesh, only: bzs, E_k_imag_mesh, activeComm, local_points
    use TightBinding, only: nOrb,nOrb2
    use mod_mpi_pars
    implicit none
    integer  :: i,j, AllocateStatus
    real(double),    dimension(3) :: kp
    real(double),    dimension(:,:), allocatable :: n_orb_u, n_orb_d
    complex(double), dimension(:,:), allocatable :: gdiagud,gdiagdu
    complex(double), dimension(:,:,:,:), allocatable :: gf
    !--------------------- begin MPI vars --------------------
    integer*8 :: ix
    integer :: ncount
    integer :: mu,mup
    real(double) :: weight, ep
    ncount = s%nAtoms * nOrb

    allocate(n_orb_u(s%nAtoms,nOrb), n_orb_d(s%nAtoms,nOrb), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[calcMagnetization] Not enough memory for: n_orb_u, n_orb_d")

    allocate(gdiagud(s%nAtoms,nOrb), gdiagdu(s%nAtoms,nOrb), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[calcMagnetization] Not enough memory for: gdiagdu, gdiagud")

    n_orb_u = 0.d0
    n_orb_d = 0.d0

    gdiagud = cZero
    gdiagdu = cZero

    !$omp parallel default(none) &
    !$omp& private(ix,ep,kp,weight,i,mu,mup,gf,AllocateStatus) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,wght,s,bzs,E_k_imag_mesh,y,n_orb_u,n_orb_d,gdiagud,gdiagdu)
    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat=AllocateStatus)
    if(AllocateStatus /= 0) call AbortProgram("[calcMagnetization] Not enough memory for: gf")
    gf = cZero

    if(llineargfsoc .or. llinearsoc) then
      !$omp do schedule(static) reduction(+:n_orb_u) reduction(+:n_orb_d) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         call greenlineargfsoc(s%Ef,ep,kp,gf)
         do i=1,s%nAtoms
           do mu=1,nOrb
             mup = mu+nOrb
             n_orb_u(i,mu) = n_orb_u(i,mu) + real(gf(mu,mu,i,i)) * weight
             n_orb_d(i,mu) = n_orb_d(i,mu) + real(gf(mup,mup,i,i)) * weight
             gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i) * weight
             gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i) * weight
           end do
         end do
      end do
      !$omp end do
    else
      !$omp do schedule(static) reduction(+:n_orb_u) reduction(+:n_orb_d) reduction(+:gdiagud) reduction(+:gdiagdu)
      do ix = 1, local_points
         ep = y(E_k_imag_mesh(1,ix))
         kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
         weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
         call green(s%Ef,ep,kp,gf)
         do i=1,s%nAtoms
           do mu=1,nOrb
             mup = mu+nOrb
             n_orb_u(i,mu) = n_orb_u(i,mu) + real(gf(mu,mu,i,i)) * weight
             n_orb_d(i,mu) = n_orb_d(i,mu) + real(gf(mup,mup,i,i)) * weight
             gdiagud(i,mu) = gdiagud(i,mu) + gf(mu,mup,i,i) * weight
             gdiagdu(i,mu) = gdiagdu(i,mu) + gf(mup,mu,i,i) * weight
           end do
         end do
      end do
      !$omp end do
    end if

    deallocate(gf)
    !$omp end parallel

    do j=1,s%nAtoms
      mp(:,j)= gdiagdu(j,:) + conjg(gdiagud(j,:))
      mpd(j) = sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9)))
    end do

    call MPI_Allreduce(MPI_IN_PLACE, n_orb_u, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, n_orb_d, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mpd, s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

    n_orb_u = 0.5d0 + n_orb_u/pi
    n_orb_d = 0.5d0 + n_orb_d/pi
    mp      = mp/pi
    mpd     = mpd/pi
    mx      = real(mp)
    my      = aimag(mp)
    mxd     = real(mpd)
    myd     = aimag(mpd)

    do i = 1, s%nAtoms
      rho(:,i) = n_orb_u(i, : ) + n_orb_d(i, : )
      mz (:,i) = n_orb_u(i, : ) - n_orb_d(i, : )
      rhod(i)  = sum(n_orb_u(i,5:9)) + sum(n_orb_d(i,5:9))
      mzd(i)   = sum(n_orb_u(i,5:9)) - sum(n_orb_d(i,5:9))
    end do

    deallocate(n_orb_u, n_orb_d)
    deallocate(gdiagdu, gdiagud)

    return
  end subroutine calcMagnetization

  subroutine calcJacobian(jacobian, N)
    !! Calculated the Jacobian of the spin magnetization
    use mod_f90_kind, only: double
    use mod_constants, only: cI, pi, identorb18, cZero, pauli_dorb, pauli_orb, cOne
    use mod_parameters, only: U, offset, eta
    use mod_SOC, only: llinearsoc, llineargfsoc
    use EnergyIntegration, only: y, wght
    use mod_system, only: s => sys
    use adaptiveMesh, only: local_points, E_k_imag_mesh, bzs, activeComm
    use mod_BrillouinZone, only: realBZ
    use TightBinding, only: nOrb,nOrb2
    use mod_mpi_pars
    implicit none
    integer, intent(in) :: N
    real(double), intent(inout), dimension(N,N) :: jacobian
    complex(double), dimension(:,:), allocatable :: gij,gji,temp,paulitemp
    complex(double), dimension(:,:,:), allocatable :: temp1,temp2
    complex(double), dimension(:,:,:,:), allocatable :: gf,gvg
    integer :: firstPoint, lastPoint

    integer :: i,j
    integer :: AllocateStatus
    integer :: i0,j0,sigma,sigmap

    real(double) :: kp(3), ep
    complex(double) :: weight

    complex(double) :: mhalfU(4,s%nAtoms)
    complex(double), dimension(nOrb2, nOrb2, 4) :: pauli_components1,pauli_components2

    !--------------------- begin MPI vars --------------------
    integer :: ix,mu,ncount2
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount2=N*N

    pauli_components1 = cZero ! Pauli Matrices in Spin-Orbit space
    pauli_components2 = cZero ! drho / dx =
!   Includes d orbitals
    pauli_components1(:,:,1) = identorb18(:,:)
    pauli_components1(:,:,2) = pauli_orb(1,:,:)
    pauli_components1(:,:,3) = pauli_orb(2,:,:)
    pauli_components1(:,:,4) = pauli_orb(3,:,:)
!   Excludes d orbitals
    pauli_components2(5:9, 5:9, 1) = identorb18(5:9, 5:9)           ! sigma_0 (1,1) = 1
    pauli_components2(14:18, 14:18, 1) = identorb18(14:18, 14:18)   ! sigma_0 (2,2) = 1
    pauli_components2(:,:,2) = pauli_dorb(1,:,:)                    ! sigma_1 and so on
    pauli_components2(:,:,3) = pauli_dorb(2,:,:)
    pauli_components2(:,:,4) = pauli_dorb(3,:,:)

    ! Prefactor -U/2 in dH/dm and -U/2 in dH/dn
    do i=1,s%nAtoms
      !mhalfU(1,i) = cOne
      mhalfU(1:4,i) = -0.5d0*U(i+offset)
    end do

    jacobian = 0.d0

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,j,i0,j0,mu,sigma,sigmap,ep,kp,weight,gf,gvg,gij,gji,temp,temp1,temp2,paulitemp) &
    !$omp& shared(llineargfsoc,llinearsoc,local_points,s,realBZ,bzs,E_k_imag_mesh,y,eta,firstPoint,lastPoint,wght,mhalfU,pauli_components1,pauli_components2,jacobian)
    allocate( temp1(nOrb2, nOrb2, 4), &
              temp2(nOrb2, nOrb2, 4), &
              gij(nOrb2,nOrb2), gji(nOrb2,nOrb2), &
              gf(nOrb2, nOrb2, s%nAtoms, s%nAtoms), &
              temp(nOrb2, nOrb2), paulitemp(nOrb2, nOrb2), stat = AllocateStatus)
    if (AllocateStatus/=0) call abortProgram("[sumk_jacobian] Not enough memory for: temp1, temp2, gij, gji, gf, temp")
    temp1     = cZero
    temp2     = cZero
    paulitemp = cZero
    gij       = cZero
    gji       = cZero
    gf        = cZero
    temp      = cZero

    if(llineargfsoc .or. llinearsoc) then
      allocate(gvg(nOrb2, nOrb2, s%nAtoms, s%nAtoms), STAT = AllocateStatus  )
      if (AllocateStatus/=0) call abortProgram("[sumk_jacobian] Not enough memory for: gvg")
      gvg = cZero
    end if

   !$omp do schedule(static) reduction(+:jacobian)
   do ix = 1, local_points
      ep = y(E_k_imag_mesh(1,ix))
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      weight = cmplx(1.d0,0.d0) * wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
      ! Green function on energy Ef + iy, and wave vector kp
      if(llineargfsoc .or. llinearsoc) then
        call greenlinearsoc(s%Ef,ep,kp,gf,gvg)
        gf = gf + gvg
      else
        call green(s%Ef,ep,kp,gf)
      end if

      do j=1,s%nAtoms
        do i=1,s%nAtoms
          gij = gf(:,:,i,j)
          gji = gf(:,:,j,i)

          do sigma = 1,4
            ! temp1 =  pauli*g_ij
            paulitemp = pauli_components1(:,:, sigma)
            call zgemm('n','n',18,18,18,cOne,paulitemp,18,gij,18,cZero,temp,18)
            temp1(:,:, sigma) = temp
          end do

          do sigma = 1,4
            ! temp2 = (-U/2) * sigma* g_ji
            paulitemp = pauli_components2(:,:, sigma)
            call zgemm('n','n',18,18,18,mhalfU(sigma,j),paulitemp,18,gji,18,cZero,temp,18)
            temp2(:,:, sigma) = temp
          end do

          do sigma = 1,4
            do sigmap = 1,4
              ! gdHdxg = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji
              i0 = (sigma-1)*s%nAtoms + i
              j0 = (sigmap-1)*s%nAtoms + j
              gij = temp1(:,:, sigma)
              gji = temp2(:,:, sigmap)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)
              do mu = 1, nOrb2
                if(((mu>=5).and.(mu<=9)).or.((mu>=14).and.(mu<=18))) then
                  jacobian(i0, j0) = jacobian(i0, j0) + real(temp(mu,mu))
                end if

                if(sigma==1) jacobian(4*s%nAtoms+1, j0) = jacobian(4*s%nAtoms+1, j0) + real(temp(mu,mu))
              end do
            end do
          end do

          if(llineargfsoc .or. llinearsoc) then ! non-linear term
            gij = gvg(:,:,i,j)
            gji = gvg(:,:,j,i)

            do sigma = 1,4
              ! temp1 = wkbz* pauli*gvg_ij
              paulitemp = pauli_components1(:,:,sigma)
              call zgemm('n','n',18,18,18,cOne,paulitemp,18,gij,18,cZero,temp,18)
              temp1(:,:,sigma) = temp
            end do

            do sigmap = 1,4
              ! temp2 = (-U/2) * sigma* gvg_ji
              paulitemp = pauli_components2(:,:,sigmap)
              call zgemm('n','n',18,18,18,mhalfU(sigmap,j),paulitemp,18,gji,18,cZero,temp,18)
              temp2(:,:,sigmap) = temp
            end do
            do sigma = 1,4
              do sigmap = 1,4
                ! gdHdxg = temp1*temp2 = wkbz* pauli*gvg_ij*(-U/2)*sigma* gvg_ji
                i0 = (sigma-1)*s%nAtoms + i
                j0 = (sigmap-1)*s%nAtoms + j
                gij = temp1(:,:,sigma)
                gji = temp2(:,:,sigmap)
                call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)
                do mu = 1, nOrb2
                  if(((mu>=5).and.(mu<=9)).or.((mu>=14).and.(mu<=18))) then
                    jacobian(i0, j0) = jacobian(i0, j0) - real(temp(mu,mu)) ! Removing non-linear term
                  end if

                  if(sigma==1) jacobian(4*s%nAtoms+1, j0) = jacobian(4*s%nAtoms+1, j0) - real(temp(mu,mu)) ! Removing non-linear term
                end do
              end do
            end do
          end if ! End linear part
        end do ! End nAtoms i loop
      end do ! End nAtoms j loop
    end do ! End Energy+nkpt loop
    !$omp end do nowait

    !$omp do schedule(static) reduction(+:jacobian)
    do ix = 1, realBZ%workload
      kp = realBZ%kp(1:3,ix)
      weight = cmplx(1.d0,0.d0) * realBZ%w(ix)

      ! Green function at Ef + ieta, and wave vector kp
      if(llineargfsoc .or. llinearsoc) then
        call greenlinearsoc(s%Ef,eta,kp,gf,gvg)
        gf = gf + gvg
      else
        call green(s%Ef,eta,kp,gf)
      end if

      do i=1,s%nAtoms
        gij = gf(:,:,i,i)

        do sigma = 1,4
          ! temp1 =  pauli*g_ii
          paulitemp = pauli_components1(:,:, sigma)
          call zgemm('n','n',18,18,18,weight,paulitemp,18,gij,18,cZero,temp,18)

          i0 = (sigma-1)*s%nAtoms + i
          do mu = 1, nOrb2
            if(((mu>=5).and.(mu<=9)).or.((mu>=14).and.(mu<=18))) then
              jacobian(i0, 4*s%nAtoms+1) = jacobian(i0, 4*s%nAtoms+1) - aimag(temp(mu,mu))
            end if
            if(sigma==1) jacobian(4*s%nAtoms+1, 4*s%nAtoms+1) = jacobian(4*s%nAtoms+1, 4*s%nAtoms+1) - aimag(temp(mu,mu))
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
    do i = 1, 4*s%nAtoms
      jacobian(i,i) = jacobian(i,i) - 1.d0
    end do

    return
  end subroutine calcJacobian

  subroutine calcLGS()
    !! Calculates the ground state charge, magnetization and orbital angular momentum ground state
    use mod_f90_kind, only: double
    use mod_constants, only: cZero,pi
    use mod_System, only: s => sys
    use adaptiveMesh
    use TightBinding, only: nOrb,nOrb2
    use mod_parameters, only: output
    use mod_magnet
    use EnergyIntegration, only: y, wght
    use mod_mpi_pars
    implicit none

    integer :: AllocateStatus
    integer*8 :: ix
    integer :: i,mu,nu,mup,nup
    real(double) :: kp(3)
    complex(double), dimension(:,:,:,:), allocatable :: gf
    complex(double), dimension(:,:,:), allocatable :: gupgd
    real(double) :: weight, ep
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
    if (AllocateStatus/=0) call abortProgram("[L_gs] Not enough memory for: df1,Fint,gf,gfuu,gfud,gfdu,gfdd")

    allocate(gupgd(nOrb, nOrb,s%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[L_gs] Not enough memory for: gupgd")

    if(rField == 0) write(output%unit_loop,"('[L_gs] Calculating Orbital Angular Momentum ground state... ')")

    ! Calculating the number of particles for each spin and orbital using a complex integral

    gupgd  = cZero

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
    !$omp& shared(local_points,s,E_k_imag_mesh,bzs,y,wght,gupgd)

    allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat = AllocateStatus)
    gf = cZero

    !$omp do schedule(static) reduction(+:gupgd)
    do ix = 1, local_points
        kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
        ep = y(E_k_imag_mesh(1,ix))
        weight = bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix)) * wght(E_k_imag_mesh(1,ix))
        !Green function on energy Ef + iy, and wave vector kp
        call green(s%Ef,ep,kp,gf)

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
    lxm = 0.d0
    lym = 0.d0
    lzm = 0.d0

    do nu=5,9
      do mu=5,9
        do i=1,s%nAtoms
          lxpm(i) = lxpm(i) + real(lxp(mu,nu,i)*gupgd(nu,mu,i))
          lypm(i) = lypm(i) + real(lyp(mu,nu,i)*gupgd(nu,mu,i))
          lzpm(i) = lzpm(i) + real(lzp(mu,nu,i)*gupgd(nu,mu,i))
          lxm(i)  = lxm(i)  + real(lx (mu,nu)*gupgd(nu,mu,i))
          lym(i)  = lym(i)  + real(ly (mu,nu)*gupgd(nu,mu,i))
          lzm(i)  = lzm(i)  + real(lz (mu,nu)*gupgd(nu,mu,i))
        end do
      end do
    end do

    ! Calculating angles of GS OAM (in units of pi)
    do i = 1,s%nAtoms
      labs(i)   = sqrt((lxm(i)**2)+(lym(i)**2)+(lzm(i)**2))
      ltheta(i) = acos(lzm(i)/sqrt(lxm(i)**2+lym(i)**2+lzm(i)**2))/pi
      lphi(i)   = atan2(lym(i),lxm(i))/pi
      lpabs(i)  = sqrt((lxpm(i)**2)+(lypm(i)**2)+(lzpm(i)**2))
      lptheta(i)= acos(lzpm(i)/sqrt(lxpm(i)**2+lypm(i)**2+lzpm(i)**2))/pi
      lpphi(i)  = atan2(lypm(i),lxpm(i))/pi
    end do

    deallocate(gupgd)

    return
  end subroutine calcLGS


  subroutine rotate_magnetization_to_field()
  !! Rotate the magnetization to the direction of the field (useful for SOC=F)
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_magnet, only: hw_count, hw_list, hhwx, hhwy, hhwz, &
                          mx, my, mz, mp
    use mod_parameters, only: output
    use mod_System, only: s => sys
    use TightBinding, only: nOrb
    use mod_mpi_pars, only: rField
    implicit none
    integer :: i,j,sign
    real(double) :: mdotb,mabs(nOrb,s%nAtoms)

    if(rField == 0) write(output%unit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    do i = 1, s%nAtoms
      do j=1,nOrb
        mdotb   = hhwx(i)*mx(j,i)+hhwy(i)*my(j,i)+hhwz(i)*mz(j,i)
        sign    = dble(mdotb/abs(mdotb))
        mabs(j,i) = sqrt((mx(j,i)**2)+(my(j,i)**2)+(mz(j,i)**2))
        mx(j,i)   = sign*mabs(j,i)*sin(hw_list(hw_count,2)*pi)*cos(hw_list(hw_count,3)*pi)
        my(j,i)   = sign*mabs(j,i)*sin(hw_list(hw_count,2)*pi)*sin(hw_list(hw_count,3)*pi)
        mz(j,i)   = sign*mabs(j,i)*cos(hw_list(hw_count,2)*pi)
        mp(j,i)   = cmplx(mx(j,i),my(j,i),double)
      end do
    end do

    ! Writing new n and rotated mag to file (without self-consistency)
    if(rField == 0) call write_sc_results()

    ! Writing self-consistency results on screen
    if(rField == 0)  call print_sc_results()

    return
  end subroutine rotate_magnetization_to_field

  ! Writes the self-consistency results on the screen
  subroutine print_sc_results()
    use mod_parameters, only: output
    use mod_system, only: s => sys
    use mod_SOC, only: SOC
    use mod_magnet, only: rho, mvec_cartesian, mp, mvec_spherical, &
                          lxpm, lypm, lzpm, lpphi, lptheta, lxm, lym, lzm, lpabs, labs
    implicit none
    integer :: i

    write(output%unit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    write(output%unit_loop,"(28x,'Ef=',f11.8)") s%Ef
    write(output%unit_loop,"(11x,' *************** Charge density: ***************')")
    do i=1,s%nAtoms
      write(output%unit_loop,"(26x,'N(',i2.0,')=',f11.8)") i, sum(rho(:,i))
    end do
    write(output%unit_loop,"(11x,' *********** Magnetization components: **********')")
    do i=1,s%nAtoms
      write(output%unit_loop,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mvec_cartesian(1,i),i,mvec_cartesian(2,i),i,mvec_cartesian(3,i)
      if(abs(sum(mp(:,i)))/=0) write(output%unit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") mvec_spherical(2,i),mvec_spherical(3,i)
    end do
    if((lGSL).or.(SOC)) then
      write(output%unit_loop,"(11x,' *** Orbital components in local frame:  ***')")
      do i=1,s%nAtoms
        write(output%unit_loop,"(4x,'Lxp(',i2.0,')=',f11.8,4x,'Lyp(',i2.0,')=',f11.8,4x,'Lzp(',i2.0,')=',f11.8)") i,lxpm(i),i,lypm(i),i,lzpm(i)
        if(sqrt(lxpm(i)**2+lypm(i)**2)/=0) write(output%unit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") lptheta(i),lpphi(i)
      end do
      write(output%unit_loop,"(11x,' *** Orbital components in global frame: ***')")
      do i=1,s%nAtoms
        write(output%unit_loop,"(4x,'Lx (',i2.0,')=',f11.8,4x,'Ly (',i2.0,')=',f11.8,4x,'Lz (',i2.0,')=',f11.8)") i,lxm(i),i,lym(i),i,lzm(i)
      end do
      write(output%unit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(output%unit_loop,"(4x,'M (',i2.0,') =',f11.8,4x,'Lp (',i2.0,')=',f11.8,4x,'L (',i2.0,') =',f11.8)") i,mvec_spherical(1,i),i,lpabs(i),i,labs(i)
      end do
    else
      write(output%unit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(output%unit_loop,"(27x,'M (',i2.0,') =',f11.8)") i,mvec_spherical(1,i)
      end do
    end if
    write(output%unit_loop,"('|----------=============================================================----------|')")

    return
  end subroutine print_sc_results


  subroutine write_sc_results()
    !! Writes the self-consistency results into files and broadcasts the scfile for the next iteration.
    use mod_parameters, only: output, dfttype
    use EnergyIntegration, only: parts
    use mod_magnet, only: rho, mx, my, mz
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_mpi_pars
    implicit none
    character(len=30) :: formatvar
    integer :: i,j

    if(rField == 0) then
      ! Writing new results (mx, my, mz and n) to file
      write(output%unit_loop,"('[write_sc_results] Writing new n, mx, my and mz to file...')")
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

    return
  end subroutine write_sc_results


  ! Writes the initial values for the self-consistency
  subroutine print_sc_step(n,mx,my,mz,Ef,fvec)
    use mod_f90_kind, only: double
    use mod_parameters, only: output
    use mod_system, only: s => sys
    use mod_magnet, only: iter
    use mod_mpi_pars
    implicit none
    real(double),dimension(4*s%nAtoms+1),optional :: fvec
    real(double),dimension(s%nAtoms) :: n,mx,my,mz
    real(double)                     :: Ef
    integer :: i

    if(rField==0) then
      if(present(fvec)) then
        do i=1,s%nAtoms
          if(abs(cmplx(mx(i),my(i)))>1.d-10) then
            write(output%unit_loop,"('Plane ',I2,': N(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,n(i),i,mx(i),i,my(i),i,mz(i)
            write(output%unit_loop,"(15x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          else
            write(output%unit_loop,"('Plane ',I2,': N(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,n(i),i,mz(i)
            write(output%unit_loop,"(15x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          end if
        end do
        write(output%unit_loop,"(13x,'Ef=',es16.9)") Ef
        write(output%unit_loop,"(15x,'fvec(',I2,')=',es16.9)") 4*s%nAtoms+1,fvec(4*s%nAtoms+1)
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

    return
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
    use mod_f90_kind, only: double
    use mod_parameters, only: output
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter,rho,mx,my,mz,rhod,mxd,myd,mzd
    use mod_Umatrix, only: update_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag
    integer     ,  intent(inout)        :: iuser(*)
    real(double),  intent(inout)        :: ruser(*)
    real(double),   dimension(N)        :: x,fvec
    real(double),   dimension(N,N)      :: selfconjac
    real(double),   dimension(s%nAtoms) :: mx_in,my_in,mz_in,rho_in
    real(double),   dimension(s%nAtoms) :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms) :: mp_in,mpd_in

    ! Values used in the hamiltonian
    rhod_in = x(           1:  s%nAtoms)
    mxd_in  = x(  s%nAtoms+1:2*s%nAtoms)
    myd_in  = x(2*s%nAtoms+1:3*s%nAtoms)
    mzd_in  = x(3*s%nAtoms+1:4*s%nAtoms)
    mpd_in  = cmplx(mxd_in,myd_in)
    s%Ef    = x(4*s%nAtoms+1)

    call update_Umatrix(mzd_in, mpd_in, rhod_in, s%nAtoms, nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    rho_in = rhod_in + sum(rho(:,:),dim=1) - rhod
    mx_in  = mxd_in  + sum(mx (:,:),dim=1) - mxd
    my_in  = myd_in  + sum(my (:,:),dim=1) - mxd
    mz_in  = mzd_in  + sum(mz (:,:),dim=1) - mxd
    mp_in  = cmplx(mx_in,my_in)

    select case (iflag)
    case(1)
      call calcMagnetization()
      do i = 1, s%nAtoms
        fvec(i           ) = rhod(i) - rhod_in(i)
        fvec(i+1*s%nAtoms) =  mxd(i) -  mxd_in(i)
        fvec(i+2*s%nAtoms) =  myd(i) -  myd_in(i)
        fvec(i+3*s%nAtoms) =  mzd(i) -  mzd_in(i)
      end do
      fvec(4*s%nAtoms+1) = sum(rho) - s%totalOccupation

      call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

      if(lontheflysc) call write_sc_results()
    case(2)
      call calcJacobian(selfconjac, N)
    case default
      write(output%unit,"('[sc_equations_and_jacobian] Problem in self-consistency! iflag = ',I0)") iflag
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select

    iter = iter + 1

    return
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
    use mod_magnet, only: iter,rho,mx,my,mz,rhod,mxd,myd,mzd
    use mod_Umatrix, only: update_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag
    integer     ,   intent(inout)       :: iuser(*)
    real(double),   intent(inout)       :: ruser(*)
    real(double),   dimension(N)        :: x,fvec
    real(double),   dimension(s%nAtoms) :: mx_in,my_in,mz_in,rho_in
    real(double),   dimension(s%nAtoms) :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms) :: mp_in,mpd_in

    iflag=0
    ! Values used in the hamiltonian
    rhod_in = x(           1:  s%nAtoms)
    mxd_in  = x(  s%nAtoms+1:2*s%nAtoms)
    myd_in  = x(2*s%nAtoms+1:3*s%nAtoms)
    mzd_in  = x(3*s%nAtoms+1:4*s%nAtoms)
    mpd_in  = cmplx(mxd_in,myd_in)
    s%Ef    = x(4*s%nAtoms+1)

    call update_Umatrix(mzd_in, mpd_in, rhod_in, s%nAtoms, nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    rho_in = rhod_in + sum(rho(:,:),dim=1) - rhod
    mx_in  = mxd_in  + sum(mx(:,:), dim=1) - mxd
    my_in  = myd_in  + sum(my(:,:), dim=1) - mxd
    mz_in  = mzd_in  + sum(mz(:,:), dim=1) - mxd
    mp_in  = cmplx(mx_in,my_in)

    call calcMagnetization()
    do i = 1, s%nAtoms
        fvec(i           ) = rhod(i) - rhod_in(i)
        fvec(i+1*s%nAtoms) =  mxd(i) -  mxd_in(i)
        fvec(i+2*s%nAtoms) =  myd(i) -  myd_in(i)
        fvec(i+3*s%nAtoms) =  mzd(i) -  mzd_in(i)
    end do
    fvec(4*s%nAtoms+1) = sum(rho) - s%totalOccupation

    call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

    if(lontheflysc) call write_sc_results()

    iter = iter + 1

    return
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
    use mod_f90_kind, only: double
    use mod_parameters, only: output
    use mod_system, only: s => sys
    use TightBinding, only: nOrb
    use mod_magnet, only: iter,rho,mx,my,mz,rhod,mxd,myd,mzd
    use mod_Umatrix, only: update_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag,ldfjac
    real(double),   dimension(N)        :: x,fvec
    real(double),   dimension(ldfjac,N) :: selfconjac
    real(double),   dimension(s%nAtoms) :: mx_in,my_in,mz_in,rho_in
    real(double),   dimension(s%nAtoms) :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms) :: mp_in,mpd_in

    ! Values used in the hamiltonian
    rhod_in = x(           1:  s%nAtoms)
    mxd_in  = x(  s%nAtoms+1:2*s%nAtoms)
    myd_in  = x(2*s%nAtoms+1:3*s%nAtoms)
    mzd_in  = x(3*s%nAtoms+1:4*s%nAtoms)
    mpd_in  = cmplx(mxd_in,myd_in)
    s%Ef    = x(4*s%nAtoms+1)

if(myrank==0) then
  write(*,*) iter
  write(*,*) "Fe", rho(1,1), sum(rho(2:4,1)), sum(rho(5:9,1))
end if

    call update_Umatrix(mzd_in, mpd_in, rhod_in, s%nAtoms, nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    rho_in = rhod_in + sum(rho(:,:),dim=1) - rhod
if(myrank==0) then
  write(*,*) "Fe1", 'd' ,rhod_in(1), 't', sum(rho(:,1),dim=1), 'd',rhod(1),'sp',sum(rho(:,1),dim=1)-rhod(1)
end if
    mx_in  = mxd_in  + sum(mx(:,:), dim=1) - mxd
    my_in  = myd_in  + sum(my(:,:), dim=1) - mxd
    mz_in  = mzd_in  + sum(mz(:,:), dim=1) - mxd
    mp_in  = cmplx(mx_in,my_in)

    flag: select case (iflag)
    case(1)
      call calcMagnetization()
      do i = 1, s%nAtoms
        fvec(i           ) = rhod(i) - rhod_in(i)
        fvec(i+1*s%nAtoms) =  mxd(i) -  mxd_in(i)
        fvec(i+2*s%nAtoms) =  myd(i) -  myd_in(i)
        fvec(i+3*s%nAtoms) =  mzd(i) -  mzd_in(i)
      end do
      fvec(4*s%nAtoms+1) = sum(rho) - s%totalOccupation
if(myrank==0) then
  write(*,*) "Fe2", sum(rho(:,1)), rho_in(1)
  write(*,*) "Fe2", sum(rho), s%totalOccupation
end if

      call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

      if(lontheflysc) call write_sc_results()
    case(2)
      call calcJacobian(selfconjac, N)
    case default
      write(output%unit,"('[sc_eqs_and_jac_old] Problem in self-consistency! iflag = ',I0)") iflag
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select flag

    iter = iter + 1

    return
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
    use mod_magnet, only: iter,rho,mx,my,mz,rhod,mxd,myd,mzd
    use mod_Umatrix, only: update_Umatrix
    use mod_mpi_pars
    implicit none
    integer  :: N,i,iflag
    real(double),   dimension(N)        :: x,fvec
    real(double),   dimension(s%nAtoms) :: mx_in,my_in,mz_in,rho_in
    real(double),   dimension(s%nAtoms) :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms) :: mp_in,mpd_in

    iflag=0
    ! Values used in the hamiltonian
    rhod_in = x(           1:  s%nAtoms)
    mxd_in  = x(  s%nAtoms+1:2*s%nAtoms)
    myd_in  = x(2*s%nAtoms+1:3*s%nAtoms)
    mzd_in  = x(3*s%nAtoms+1:4*s%nAtoms)
    mpd_in  = cmplx(mxd_in,myd_in)
    s%Ef    = x(4*s%nAtoms+1)

    call update_Umatrix(mzd_in, mpd_in, rhod_in, s%nAtoms, nOrb)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,s%Ef)

    rho_in = rhod_in + sum(rho(:,:),dim=1) - rhod
    mx_in  = mxd_in  + sum(mx(:,:), dim=1) - mxd
    my_in  = myd_in  + sum(my(:,:), dim=1) - mxd
    mz_in  = mzd_in  + sum(mz(:,:), dim=1) - mxd
    mp_in  = cmplx(mx_in,my_in)

    call calcMagnetization()
    do i = 1, s%nAtoms
      fvec(i           ) = rhod(i) - rhod_in(i)
      fvec(i+1*s%nAtoms) =  mxd(i) -  mxd_in(i)
      fvec(i+2*s%nAtoms) =  myd(i) -  myd_in(i)
      fvec(i+3*s%nAtoms) =  mzd(i) -  mzd_in(i)
    end do
    fvec(4*s%nAtoms+1) = sum(rho) - s%totalOccupation

    call print_sc_step(rhod,mxd,myd,mzd,s%Ef,fvec)

    if(lontheflysc) call write_sc_results()

    iter = iter + 1

    return
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
    use mod_magnet, only: iter
    use mod_Umatrix, only: update_Umatrix
    implicit none
    integer       :: N,ldfjac,iflag
    real(double)  :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(double),   dimension(s%nAtoms) :: mxd_in,myd_in,mzd_in,rhod_in
    complex(double),dimension(s%nAtoms) :: mpd_in
    !--------------------- begin MPI vars --------------------

    iflag=0
    ! Values used in the hamiltonian
    rhod_in = x(           1:  s%nAtoms)
    mxd_in  = x(  s%nAtoms+1:2*s%nAtoms)
    myd_in  = x(2*s%nAtoms+1:3*s%nAtoms)
    mzd_in  = x(3*s%nAtoms+1:4*s%nAtoms)
    mpd_in  = cmplx(mxd_in,myd_in)
    s%Ef    = x(4*s%nAtoms+1)

    call update_Umatrix(mzd_in, mpd_in, rhod_in, s%nAtoms, nOrb)

    fvec=fvec

    call calcJacobian(selfconjac, N)

    iter = iter + 1

    return
  end subroutine sc_jac_old

end module mod_self_consistency
