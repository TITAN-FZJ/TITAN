module mod_self_consistency
  use mod_kind, only: dp
  implicit none
  integer            :: neq, neq_per_atom
  character(len=300) :: default_file
  real(dp)           :: mag_tol = 1.e-12_dp
  character(len=200) :: scfile = ""
  !! Give a file to start self-consistency
  logical            :: skipsc
  !! Skip self-consistency
  character(len=50)  :: magbasis = ""
  !! Basis to give initial magnetization in 'initialmag' file
  real(dp), allocatable :: initialmag(:,:)
  !! Initial guess for magnetization
  logical :: lselfcon    = .false.
  logical :: lontheflysc = .false.
  logical :: lnojac      = .false.
  logical :: lrotatemag  = .false.
  logical :: lforceoccup = .false.

contains

  subroutine doSelfConsistency()
    use mod_magnet,        only: lp_matrix, mtheta, mphi, lb_matrix, sb_matrix
    use adaptiveMesh,      only: genLocalEKMesh,freeLocalEKMesh
    use mod_BrillouinZone, only: realBZ
    use mod_mpi_pars,      only: rFreq, sFreq, FreqComm, rField, sField, FieldComm
    use mod_parameters,    only: leigenstates,lkpoints,lEf_overwrite,Ef_overwrite
    use mod_System,        only: s => sys
    use mod_expectation,   only: groundstate_L_and_E,expectation_values,expectation_eigenstates_fullhk,calc_GS_L_and_E,calc_GS_L_and_E_fullhk
    use mod_hamiltonian,   only: fullhamiltk,lfullhk
    implicit none
    logical :: lsuccess = .false.

    ! Distribute Energy Integration across all points available
    call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1),lkpoints)
    if(.not.leigenstates) then
      call genLocalEKMesh(s,rField,sField, FieldComm)
    end if

    !---- Checking if full tight-binding hamiltonian can be calculated ----
    if(leigenstates) lfullhk = fullhamiltk(s)
    if(lfullhk) then
      expectation_values => expectation_eigenstates_fullhk
      calc_GS_L_and_E    => calc_GS_L_and_E_fullhk
    end if

    !--------------------------- Self-consistency --------------------------
    ! Trying to read previous densities and Ef from files
    call read_previous_results(lsuccess)
    ! Rotate the magnetization to the direction of the field
    ! (useful for SOC=F)
    ! Only when some previous result was read from file (lsuccess=.true.)
    if(lrotatemag .and. lsuccess) then
      call rotate_magnetization_to_field()
    end if

    ! Overwriting Fermi energy with input value
    if(lEf_overwrite) s%Ef = Ef_overwrite

    ! Performing self-consistency
    if(lselfcon) call calcMagneticSelfConsistency()

    ! Writing new n and mz to file after self-consistency is done
    if(.not. lontheflysc) call write_sc_results()

    ! L matrix in local frame for given quantization direction
    call lp_matrix(mtheta, mphi)

    ! Calculating ground state Orbital Angular Momentum and Band Energy
    call groundstate_L_and_E()

    ! Writing self-consistency results on screen
    if(rField == 0)  call print_sc_results()

    if(.not.leigenstates) call freeLocalEKMesh()

  end subroutine doSelfConsistency

  ! Tries to read n and m if available
  subroutine read_previous_results(lsuccess)
    use mod_kind,              only: dp
    use mod_constants,         only: deg2rad
    use mod_parameters,        only: nOrb, output
    use mod_system,            only: s => sys
    use mod_mpi_pars,          only: rField,abortProgram
    use mod_superconductivity, only: singlet_coupling, lsuperCond
    use mod_Umatrix,           only: init_Umatrix
    use mod_magnet,            only: mx,my,mz,mxd,myd,mzd,mpd,hw_count,hw_list, &
                                     lfield,rho,rhod,rhod0,rho0
    implicit none
    integer             :: i,err, mu
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

      ! Starting magnetization
      allocate(initialmag(s%nAtoms,3))
      if(trim(magbasis)/="") then
        call read_initialmag("initialmag",trim(magbasis),s%nAtoms)
      else
        do i=1,s%nAtoms
          initialmag(i,:) = [0._dp, 0._dp, sign(2.0_dp, hw_list(hw_count,1))]
        end do

        ! Rotating magnetization to field direction
        if(lfield) then
          initialmag(:,1) = sign(2.0_dp, hw_list(hw_count,1))*sin(hw_list(hw_count,2)*deg2rad) * cos(hw_list(hw_count,3)*deg2rad)
          initialmag(:,2) = sign(2.0_dp, hw_list(hw_count,1))*sin(hw_list(hw_count,2)*deg2rad) * sin(hw_list(hw_count,3)*deg2rad)
          initialmag(:,3) = sign(2.0_dp, hw_list(hw_count,1))*cos(hw_list(hw_count,2)*deg2rad)
        end if

      end if

      do i=1,s%nAtoms
        mxd(i) = initialmag(i,1)
        myd(i) = initialmag(i,2)
        mzd(i) = initialmag(i,3)
      end do
      mpd = cmplx(mxd,myd,dp)

      mx = 0._dp
      my = 0._dp
      mz = 0._dp
      ! Initialize hamiltonian using occupations from elemental files and rho=rho0
      rho = rho0
      do i = 1, s%nAtoms
        mx(5:9,i) = 0.2_dp*mxd(i)
        my(5:9,i) = 0.2_dp*myd(i)
        mz(5:9,i) = 0.2_dp*mzd(i)
        rhod(i)   = rhod0(i) !s%Types(s%Basis(i)%Material)%OccupationD
        if(lsuperCond) then
          do mu = 1,nOrb
              singlet_coupling(mu,i) = s%Types(s%Basis(i)%Material)%lambda(mu)
          end do
        end if
      end do
      ! singlet_coupling = 1.0

    end if

    call init_Umatrix(mzd,mpd,rhod,rhod0,rho,rho0,s%nAtoms,nOrb)
  end subroutine read_previous_results

  ! Reads the initial magnetization for nAtoms in 'filename'
  subroutine read_initialmag(filename,basis,nAtoms)
    use mod_kind,       only: dp
    use mod_parameters, only: output
    use mod_constants,  only: deg2rad
    use mod_tools,      only: read_data
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: rField,abortProgram
    implicit none
    character(len=*), intent(in) :: filename,basis
    integer         , intent(in) :: nAtoms
    integer  :: i,err
    real(dp) :: temp(nAtoms,3)

    if(rField == 0) &
    write(output%unit_loop,"('[read_initialmag] Reading initial magnetization from file ',a)") "'" // filename // "'"

    open(unit=321,file=trim(filename),status="old",iostat=err)
    if((rField == 0).and.(err/=0)) call abortProgram("[read_initialmag] File '" // trim(filename) // "'' does not exist!")

    select case(basis)
    case("cartesian","c")
      call read_data(321,nAtoms,3,initialmag)
    case("spherical","s")
      call read_data(321,nAtoms,3,temp)
      do i=1,nAtoms
        initialmag(i,1) = temp(i,1)*sin(temp(i,2)*deg2rad)*cos(temp(i,3)*deg2rad)
        initialmag(i,2) = temp(i,1)*sin(temp(i,2)*deg2rad)*sin(temp(i,3)*deg2rad)
        initialmag(i,3) = temp(i,1)*cos(temp(i,2)*deg2rad)
      end do
    case("bravais","b")
      call read_data(321,nAtoms,3,temp)
      do i=1,nAtoms
        initialmag(i,:) = temp(i,1)*s%a1 + temp(i,2)*s%a2 + temp(i,3)*s%a3
      end do
    case("neighbor","n")
      call read_data(321,nAtoms,1,temp)
      call abortProgram("[read_initialmag] initialmag not implemented for basis = neighbor")
      ! do i=1,nAtoms
      !   initialmag(i,:) = temp(i,1)
      ! end do
    case default
      call abortProgram("[read_initialmag] basis wrongly defined: " // trim(basis))
    end select

    close(unit=321)

  end subroutine read_initialmag

  ! This subroutine reads previous band-shifting and magnetization results
  subroutine read_sc_results(err,lsuccess)
    use mod_kind, only: dp
    use mod_parameters,        only: nOrb, output, dfttype
    use EnergyIntegration,     only: parts
    use mod_system,            only: s => sys
    use mod_superconductivity, only: singlet_coupling
    use mod_magnet,            only: rho, mp, mx, my, mz, rhod, mpd, mxd, myd, mzd, hw_count
    use mod_mpi_pars,          only: rField,FieldComm,ierr,MPI_DOUBLE_PRECISION
    implicit none
    character(len=300)  :: file = ""
    integer,intent(out) :: err
    logical,intent(out) :: lsuccess
    integer             :: i,j
    real(dp)            :: previous_results(5*nOrb,s%nAtoms), previous_Ef

    external :: MPI_Bcast

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
      write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),dfttype,parts,trim(output%BField),trim(output%info),trim(output%SOC),trim(output%suffix)
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
        write(file,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),dfttype,parts,trim(output%BField),trim(output%info),trim(output%SOC),trim(output%suffix)
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
          read(99,fmt=*) (previous_results(j,i), j=1,5*nOrb)
        end do
        read(99,fmt=*) previous_Ef
      end if

      call MPI_Bcast(previous_results,5*nOrb*s%nAtoms,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
      call MPI_Bcast(previous_Ef,1,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

      rho(:,:) = previous_results(       1:  nOrb,:)
      mx (:,:) = previous_results(  nOrb+1:2*nOrb,:)
      my (:,:) = previous_results(2*nOrb+1:3*nOrb,:)
      mz (:,:) = previous_results(3*nOrb+1:4*nOrb,:)
      mp       = cmplx(mx,my,dp)
      singlet_coupling(:,:) = previous_results(4*nOrb+1:5*nOrb,:)

      rhod(:) = sum(rho(5:9,:),dim=1)
      mxd(:)  = sum(mx (5:9,:),dim=1)
      myd(:)  = sum(my (5:9,:),dim=1)
      mzd(:)  = sum(mz (5:9,:),dim=1)
      mpd     = cmplx(mxd,myd,dp)

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
            read(99,fmt=*) (previous_results(j,i), j=1,5*nOrb)
          end do
          read(99,fmt=*) previous_Ef
        end if

        call MPI_Bcast(previous_results,5*nOrb*s%nAtoms,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
        call MPI_Bcast(previous_Ef,1,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

        rho(:,:) = previous_results(       1:  nOrb,:)
        mx (:,:) = previous_results(  nOrb+1:2*nOrb,:)
        my (:,:) = previous_results(2*nOrb+1:3*nOrb,:)
        mz (:,:) = previous_results(3*nOrb+1:4*nOrb,:)
        mp       = cmplx(mx,my,dp)
        singlet_coupling(:,:) = previous_results(4*nOrb+1:5*nOrb,:)

        rhod(:) = sum(rho(5:9,:),dim=1)
        mxd(:)  = sum(mx (5:9,:),dim=1)
        myd(:)  = sum(my (5:9,:),dim=1)
        mzd(:)  = sum(mz (5:9,:),dim=1)
        mpd     = cmplx(mxd,myd,dp)

        call calcMagAngle()

        s%Ef = previous_Ef
        lsuccess = .true. ! Read...
        err = 1           ! ... different parameters
      end if
    end if
    close(99)
  end subroutine read_sc_results

  subroutine calcMagAngle()
    use mod_constants,        only: rad2deg
    use mod_system,           only: s => sys
    use mod_magnet,           only: mx, my, mz, mabs, &
                                    mtheta, mphi, mvec_cartesian, &
                                    mvec_spherical,mtotal_cartesian,mtotal_spherical, lrot
    implicit none
    integer :: i

    ! Calculating new angles of GS magnetization in units of pi and magnetization vector
    do i = 1,s%nAtoms
      mabs(i)   = sqrt((sum(mx(:,i))**2)+(sum(my(:,i))**2)+(sum(mz(:,i))**2))
      if(mabs(i)>1.e-8_dp) then
        mtheta(i) = acos(sum(mz(:,i))/mabs(i))*rad2deg
      else
        mtheta(i) = 0._dp
      end if
      if(abs(mtheta(i))>1.e-8_dp) then
        if(abs(abs(mtheta(i))-180._dp)>1.e-8_dp) then
          mphi(i)   = atan2(sum(my(:,i)),sum(mx(:,i)))*rad2deg
        else
          mphi(i) = 0._dp
        end if
        lrot = .true. ! Susceptibilities need to be rotated
      else
        mphi(i) = 0._dp
      end if
      mvec_cartesian(1,i) = sum(mx(:,i))
      mvec_cartesian(2,i) = sum(my(:,i))
      mvec_cartesian(3,i) = sum(mz(:,i))
      mvec_spherical(1,i) = mabs(i)
      mvec_spherical(2,i) = mtheta(i)
      mvec_spherical(3,i) = mphi(i)
    end do

    mtotal_cartesian(1) = sum(mvec_cartesian(1,:))
    mtotal_cartesian(2) = sum(mvec_cartesian(2,:))
    mtotal_cartesian(3) = sum(mvec_cartesian(3,:))

    mtotal_spherical(1) = sqrt((mtotal_cartesian(1)**2)+(mtotal_cartesian(2)**2)+(mtotal_cartesian(3)**2))
    if(mtotal_spherical(1)>1.e-8_dp) then
      mtotal_spherical(2) = acos(mtotal_cartesian(3)/mtotal_spherical(1))*rad2deg
    else
      mtotal_spherical(2) = 0._dp
    end if
    if((abs(mtotal_spherical(2))>1.e-8_dp).and.(abs(abs(mtotal_spherical(2))-180._dp)>1.e-8_dp)) then
      mtotal_spherical(3)   = atan2(mtotal_cartesian(2),mtotal_cartesian(1))*rad2deg
    else
      mtotal_spherical(3) = 0._dp
    end if

  end  subroutine calcMagAngle

  subroutine calcMagneticSelfConsistency()
  !! This subroutine performs the self-consistency
    use mod_kind,              only: dp
    use mod_parameters,        only: nOrb,output,lfixEf
    use mod_magnet,            only: rho,mxd,myd,mzd
    use mod_mpi_pars,          only: rField
    use mod_system,            only: s => sys
    use mod_dnsqe,             only: dnsqe
    use mod_superconductivity, only: superCond, lsuperCond, singlet_coupling
    implicit none
    real(dp), allocatable :: fvec(:),jac(:,:),wa(:),sc_solu(:)
    real(dp), allocatable :: diag(:),qtf(:)
    integer               :: i,mu,lwa,ifail=0

    neq_per_atom = 8 + (superCond-1)*nOrb
    neq = neq_per_atom*s%nAtoms
    if(.not.lfixEf) neq = neq+1
    allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )

    ! Putting read n and m existing solutions into sc_solu (first guess of the subroutine)
    ! Initialization of the values
    do i = 1, s%nAtoms
      do mu = 5,nOrb
        sc_solu((i-1)*neq_per_atom+(mu-4)) = rho(mu,i)
      end do
      sc_solu((i-1)*neq_per_atom+6) = mxd(i)
      sc_solu((i-1)*neq_per_atom+7) = myd(i)
      sc_solu((i-1)*neq_per_atom+8) = mzd(i)
      if(lsuperCond) then
        do mu = 1, nOrb
          sc_solu((i-1)*neq_per_atom+8+mu) = singlet_coupling(mu,i)
        end do
      end if
    end do
    if(.not.lfixEf) sc_solu(neq_per_atom*s%nAtoms+1) = s%Ef

    if(rField == 0) &
      write(output%unit_loop,"('[self_consistency] Starting self-consistency:')")

    lwa=neq*(3*neq+13)/2
    allocate( wa(lwa) )
    if(lnojac) then
      call dnsqe(sc_eqs,sc_jac,2,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
    else
      call dnsqe(sc_eqs,sc_jac,1,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
    end if
    ifail = ifail-1

    deallocate(sc_solu,diag,qtf,fvec,jac,wa)

    ! Calculating the magnetization in cartesian and spherical coordinates
    call calcMagAngle()
  end subroutine calcMagneticSelfConsistency


  subroutine check_jacobian(neq,x)
    use mod_parameters, only: output,lcheckjac
    use mod_mpi_pars,   only: rField,abortProgram
    use mod_chkder,     only: chkder
    implicit none
    integer :: liw,lw,ifail
    integer , allocatable :: iw(:)
    real(dp), allocatable :: w(:)
    integer , intent(in)  :: neq
    real(dp), intent(in)  :: x(neq)
    real(dp) :: jac(neq,neq)
    real(dp) :: fvec(neq)
    ! integer :: i
    real(dp) :: fvecp(neq),xp(neq),err(neq)

    liw = 1
    lw  = (4+neq)*neq
    allocate(iw(liw),w(lw))

    lcheckjac = .false.

    if(rField == 0) &
    write(output%unit_loop,"('[check_jacobian] Checking Jacobian if Jacobian is correct...')", advance='no')

    call chkder(neq,neq,x,fvec,jac,neq,xp,fvecp,1,err)
    call abortProgram("[check_jacobian] Not implemented for slatec).")

    if(ifail == 0) then
      if(rField == 0) write(output%unit_loop,"(' YES! ')")
    else
      if(rField == 0) write(output%unit_loop,"(' NO! ifail = ',i0)") ifail
    end if

    lcheckjac = .true.

  end subroutine check_jacobian

  subroutine lsqfun(iflag,M,N,x,fvec,selfconjac,ljc,iw,liw,w,lw)
    use mod_kind,              only: dp
    use mod_system,            only: s => sys
    use mod_magnet,            only: rho,rhod,mp,mx,my,mz,mpd,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix,           only: update_Umatrix
    use mod_parameters,        only: nOrb,lfixEf
    use mod_expectation,       only: expectation_values
    use mod_superconductivity, only: lsuperCond, update_singlet_couplings
    ! use mod_mpi_pars
    implicit none
    integer      :: M,N,ljc,i,mu,iflag,lw,liw,iw(liw)
    real(dp)     :: w(lw)
    real(dp),    dimension(N)             :: x,fvec
    real(dp),    dimension(M,ljc)         :: selfconjac
    real(dp),    dimension(nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    real(dp),    dimension(nOrb,s%nAtoms) :: singlet_coupling_in
    real(dp),    dimension(nOrb,s%nAtoms) :: deltas
    complex(dp), dimension(s%nAtoms)      :: mpd_in

    w=w
    iw=iw
    iflag=iflag

    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1,s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*neq_per_atom+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*neq_per_atom+6)
      myd_in(i) = x((i-1)*neq_per_atom+7)
      mzd_in(i) = x((i-1)*neq_per_atom+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      if(lsupercond) then
        do mu = 1,nOrb
          singlet_coupling_in(mu,i) = x((i-1)*(neq_per_atom)+8+mu)
        end do
      end if
    end do
    if(.not.lfixEf) s%Ef    = x(N)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)
    if(lsuperCond) call update_singlet_couplings(s,singlet_coupling_in)

    ! call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,singlet_coupling_in,s%Ef)
    call expectation_values(s,rho,mp,mx,my,mz,deltas)
    do i = 1,s%nAtoms
      rhod(i)   = sum(rho(5:9,i))
      mpd(i)    = sum(mp(5:9,i))
      mxd(i)    = sum(mx(5:9,i))
      myd(i)    = sum(my(5:9,i))
      mzd(i)    = sum(mz(5:9,i))
    end do
    do i = 1,s%nAtoms
      do mu = 5,nOrb
        fvec((i-1)*neq_per_atom+(mu-4)) = rho(mu,i) - rho_in(mu,i)
      end do
      fvec((i-1)*neq_per_atom+6) =  mxd(i) -  mxd_in(i)
      fvec((i-1)*neq_per_atom+7) =  myd(i) -  myd_in(i)
      fvec((i-1)*neq_per_atom+8) =  mzd(i) -  mzd_in(i)
      if(lsupercond) then
        do mu = 1,nOrb
          fvec((i-1)*neq_per_atom+8+mu) = deltas(mu,i) - singlet_coupling_in(mu,i)
        end do
      end if
    end do
    if(.not.lfixEf) fvec(N) = sum(rho) - s%totalOccupation

    call calcJacobian_greenfunction(selfconjac, N)

  end subroutine lsqfun

  subroutine calcJacobian_greenfunction(jacobian,N)
    !! Calculated the Jacobian of the spin magnetization
    use mod_kind,          only: dp,int64
    use mod_constants,     only: pi,ident_norb2,cZero,pauli_dorb,ident_dorb,cOne
    use mod_parameters,    only: nOrb,nOrb2,Un,Um,offset,eta,output,lfixEf
    use mod_SOC,           only: llinearsoc,llineargfsoc
    use EnergyIntegration, only: y,wght
    use mod_System,        only: s => sys
    use adaptiveMesh,      only: local_points,E_k_imag_mesh,bzs,activeComm
    use mod_BrillouinZone, only: realBZ
    use mod_hamiltonian,   only: hamilt_local
    use mod_greenfunction, only: greenlinearsoc,green
    use mod_mpi_pars,      only: rField,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_SUM,ierr,abortProgram
    implicit none
    integer,                       intent(in)    :: N
    real(dp),    dimension(N,N),   intent(inout) :: jacobian
    complex(dp), dimension(:,:),     allocatable :: gij,gji,temp,paulitemp
    complex(dp), dimension(:,:,:),   allocatable :: temp1,temp2
    complex(dp), dimension(:,:,:,:), allocatable :: gf,gvg
    integer(int64) :: ix
    integer     :: AllocateStatus
    integer     :: i,j,mu,nu,sigma,sigmap
    real(dp)    :: kp(3), ep
    complex(dp) :: weight
    complex(dp) :: halfUn(s%nAtoms),halfUm(s%nAtoms)
    complex(dp), dimension(nOrb2, nOrb2, 4) :: pauli_a,pauli_b
    integer :: ncount2

    external :: zgemm,MPI_Allreduce

    if(rField == 0) &
      write(output%unit_loop,"('[calcJacobian_greenfunction] Calculating the Jacobian...')")

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
      halfUn(i) = -0.5_dp*Un(i+offset)
      halfUm(i) = -0.5_dp*Um(i+offset)
    end do

    ! Build local hamiltonian
    if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

    jacobian = 0._dp

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,j,mu,nu,sigma,sigmap,ep,kp,weight,gf,gvg,gij,gji,temp,temp1,temp2,paulitemp) &
    !$omp& shared(llineargfsoc,llinearsoc,lfixEf,local_points,s,nOrb,nOrb2,neq_per_atom,realBZ,bzs,E_k_imag_mesh,y,eta,wght,halfUn,halfUm,pauli_a,pauli_b,jacobian)
    allocate( temp1(nOrb2, nOrb2, 4), &
              temp2(nOrb2, nOrb2, 4), &
              gij(nOrb2,nOrb2), gji(nOrb2,nOrb2), &
              gf(nOrb2, nOrb2, s%nAtoms, s%nAtoms), &
              temp(nOrb2, nOrb2), paulitemp(nOrb2, nOrb2), stat = AllocateStatus)
    if (AllocateStatus/=0) &
      call abortProgram("[calcJacobian_greenfunction] Not enough memory for: temp1, temp2, gij, gji, gf, temp")
    gf        = cZero
    temp      = cZero

    if(llineargfsoc .or. llinearsoc) then
      allocate(gvg(nOrb2, nOrb2, s%nAtoms, s%nAtoms), STAT = AllocateStatus  )
      if (AllocateStatus/=0) &
        call abortProgram("[calcJacobian_greenfunction] Not enough memory for: gvg")
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
            paulitemp(nu  ,nu  ) = paulitemp(nu  ,nu  ) - 2._dp
            paulitemp(nu+9,nu+9) = paulitemp(nu+9,nu+9) - 2._dp
            call zgemm('n','n',18,18,18,halfUn(j),paulitemp,18,gji,18,cZero,temp,18)
            temp2(:,:,1) = temp

            ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

            ! Charge density-charge density part
            gij = temp1(:,:,1)
            gji = temp2(:,:,1)
            call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

            do mu = 1,nOrb
              ! Last line (Total charge neutrality)
              if(.not.lfixEf) jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+(nu-4)) = jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+(nu-4)) + real(temp(mu,mu) + temp(mu+9,mu+9))
              if(mu<5) cycle
              jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+(nu-4)) = jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+(nu-4)) + real(temp(mu,mu) + temp(mu+9,mu+9))
            end do

            ! Magnetic density-charge density part
            do sigma = 2,4
              gij = temp1(:,:,sigma)
              gji = temp2(:,:,1)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

              do mu = 5,nOrb
                jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+(nu-4)) = jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+(nu-4)) + real(temp(mu,mu) + temp(mu+9,mu+9))
              end do
            end do

          end do

          gji = gf(:,:,j,i)

          do sigmap = 2,4
            ! Second product: temp2 = (-U/2) * sigma* g_ji
            paulitemp = pauli_b(:,:,sigmap)
            call zgemm('n','n',18,18,18,halfUm(j),paulitemp,18,gji,18,cZero,temp,18)
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
              if(.not.lfixEf) jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+4+sigmap) = jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+4+sigmap) + real(temp(mu,mu) + temp(mu+9,mu+9))
              if(mu<5) cycle
              jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+4+sigmap) = jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+4+sigmap) + real(temp(mu,mu) + temp(mu+9,mu+9))
            end do
          end do

          ! Magnetic density-magnetic density part
          do sigma = 2,4
            do sigmap = 2,4
              gij = temp1(:,:,sigma)
              gji = temp2(:,:,sigmap)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

              do mu = 5,nOrb
                jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+4+sigmap) = jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+4+sigmap) + real(temp(mu,mu) + temp(mu+9,mu+9))
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
              paulitemp(nu  ,nu  ) = paulitemp(nu  ,nu  ) - 2._dp
              paulitemp(nu+9,nu+9) = paulitemp(nu+9,nu+9) - 2._dp
              call zgemm('n','n',18,18,18,halfUn(j),paulitemp,18,gji,18,cZero,temp,18)
              temp2(:,:, 1) = temp

              ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

              ! Charge density-charge density part
              gij = temp1(:,:, 1)
              gji = temp2(:,:, 1)
              call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

              do mu = 1,nOrb
                ! Last line (Total charge neutrality)
                if(.not.lfixEf) jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+(nu-4)) = jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+(nu-4)) - real(temp(mu,mu) + temp(mu+9,mu+9))
                if(mu<5) cycle
                jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+(nu-4)) = jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+(nu-4)) - real(temp(mu,mu) + temp(mu+9,mu+9))
              end do

              ! Magnetic density-charge density part
              do sigma = 2,4
                gij = temp1(:,:, sigma)
                gji = temp2(:,:, 1)
                call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

                do mu = 5,nOrb
                  jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+(nu-4)) = jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+(nu-4)) - real(temp(mu,mu) + temp(mu+9,mu+9))
                end do
              end do

            end do

            gji = gvg(:,:,j,i)

            do sigmap = 2,4
              ! Second product: temp2 = (-U/2) * sigma* g_ji
              paulitemp = pauli_b(:,:, sigmap)
              call zgemm('n','n',18,18,18,halfUm(j),paulitemp,18,gji,18,cZero,temp,18)
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
                if(.not.lfixEf) jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+4+sigmap) = jacobian(  neq_per_atom*s%nAtoms+1,(j-1)*neq_per_atom+4+sigmap) - real(temp(mu,mu) + temp(mu+9,mu+9))
                if(mu<5) cycle
                jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+4+sigmap) = jacobian((i-1)*neq_per_atom+(mu-4),(j-1)*neq_per_atom+4+sigmap) - real(temp(mu,mu) + temp(mu+9,mu+9))
              end do
            end do

            ! Magnetic density-magnetic density part
            do sigma = 2,4
              do sigmap = 2,4
                gij = temp1(:,:, sigma)
                gji = temp2(:,:, sigmap)
                call zgemm('n','n',18,18,18,weight,gij,18,gji,18,cZero,temp,18)

                do mu = 5,nOrb
                  jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+4+sigmap) = jacobian((i-1)*neq_per_atom+4+sigma,(j-1)*neq_per_atom+4+sigmap) - real(temp(mu,mu) + temp(mu+9,mu+9))
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
    if(.not.lfixEf) then
      !$omp do schedule(static) reduction(+:jacobian)
      do ix = 1, realBZ%workload
        kp = realBZ%kp(1:3,ix)
        weight = cmplx(realBZ%w(ix),0._dp,dp)

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
            jacobian(  neq_per_atom*s%nAtoms+1,neq_per_atom*s%nAtoms+1) = jacobian(  neq_per_atom*s%nAtoms+1,neq_per_atom*s%nAtoms+1) - dimag(temp(mu,mu) + temp(mu+9,mu+9))
            if(mu<5) cycle
            jacobian((i-1)*neq_per_atom+(mu-4),neq_per_atom*s%nAtoms+1) = jacobian((i-1)*neq_per_atom+(mu-4),neq_per_atom*s%nAtoms+1) - dimag(temp(mu,mu) + temp(mu+9,mu+9))
          end do

          do sigma = 2,4
            ! temp1 =  pauli*g_ii
            paulitemp = pauli_a(:,:, sigma)
            call zgemm('n','n',18,18,18,weight,paulitemp,18,gij,18,cZero,temp,18)

            do mu = 5,nOrb
              jacobian((i-1)*neq_per_atom+4+sigma,neq_per_atom*s%nAtoms+1) = jacobian((i-1)*neq_per_atom+4+sigma,neq_per_atom*s%nAtoms+1) - dimag(temp(mu,mu) + temp(mu+9,mu+9))
            end do
          end do

          ! No linear correction is needed since it's a single Green function

        end do ! End nAtoms i loop
      end do ! End nkpt loop
      !$omp end do nowait
    end if

    deallocate(paulitemp)
    deallocate(temp1, temp2, temp)
    deallocate(gf, gij, gji)
    if(allocated(gvg)) deallocate(gvg)
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, jacobian, ncount2, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)

    jacobian = jacobian/pi
    do i = 1, neq_per_atom*s%nAtoms
      jacobian(i,i) = jacobian(i,i) - 1._dp
    end do
  end subroutine calcJacobian_greenfunction

  subroutine rotate_magnetization_to_field()
  !! Rotate the magnetization to the direction of the field (useful for SOC=F)
    use mod_kind,       only: dp
    use mod_constants,  only: deg2rad
    use mod_magnet,     only: lfield,hw_count,hw_list,hhw,mp,mx,my,mz,mpd,mxd,myd,mzd
    use mod_parameters, only: nOrb,output
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: rField
    implicit none
    integer      :: i,j,signal
    real(dp) :: mdotb,mabs(nOrb,s%nAtoms)

    if(rField == 0) &
      write(output%unit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    if(.not.lfield) then
      if(rField == 0) &
        write(output%unit_loop,"('[rotate_magnetization_to_field] Field if OFF! No rotation is done.')")
      return
    end if

    do i = 1, s%nAtoms
      do j=1,nOrb
        mdotb   = hhw(1,i)*mx(j,i)+hhw(2,i)*my(j,i)+hhw(3,i)*mz(j,i)
        signal  = int(mdotb/abs(mdotb))
        mabs(j,i) = sqrt((mx(j,i)**2)+(my(j,i)**2)+(mz(j,i)**2))
        mx  (j,i) = signal*mabs(j,i)*sin(hw_list(hw_count,2)*deg2rad)*cos(hw_list(hw_count,3)*deg2rad)
        my  (j,i) = signal*mabs(j,i)*sin(hw_list(hw_count,2)*deg2rad)*sin(hw_list(hw_count,3)*deg2rad)
        mz  (j,i) = signal*mabs(j,i)*cos(hw_list(hw_count,2)*deg2rad)
        mp  (j,i) = cmplx(mx(j,i),my(j,i),dp)
      end do
    end do

    mxd(:)  = sum(mx (5:9,:),dim=1)
    myd(:)  = sum(my (5:9,:),dim=1)
    mzd(:)  = sum(mz (5:9,:),dim=1)
    mpd     = cmplx(mxd,myd,dp)

    call calcMagAngle()

    ! Writing new n and rotated mag to file (without self-consistency)
    ! if(rField == 0) call write_sc_results()

    ! Writing self-consistency results on screen
    ! if(rField == 0) call print_sc_results()
  end subroutine rotate_magnetization_to_field

  ! Writes the self-consistency results on the screen
  subroutine print_sc_results()
    use mod_parameters,        only: output,leigenstates
    use mod_system,            only: s => sys
    use mod_hamiltonian,       only: energy
    use mod_magnet,            only: rho, mvec_cartesian, mp, mvec_spherical, &
                                     lxm, lym, lzm, ltheta, lphi, labs, iter
    use mod_superconductivity, only: lsuperCond, singlet_coupling
    implicit none
    integer :: i

    write(output%unit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    if(leigenstates) then
      write(output%unit_loop,"(14x,'Ef=',f11.8,4x,'Eband=',f15.7)") s%Ef,energy
    else
      write(output%unit_loop,"(28x,'Ef=',f11.8)") s%Ef
    end if
    write(output%unit_loop,"(11x,' *************** Charge density: ****************')")
    do i=1,s%nAtoms
      write(output%unit_loop,"(a,':',2x,'Ns(',i4.0,')=',f11.8,4x,'Np(',i4.0,')=',f11.8,4x,'Nd(',i4.0,')=',f11.8)") trim(s%Types(s%Basis(i)%Material)%Name),i, rho(1,i),i, sum(rho(2:4,i)),i, sum(rho(5:9,i))
    end do
    if(lsupercond) then
        write(output%unit_loop,"(11x,' ******** Superconducting gap parameter: ********')")
        write(output%unit_loop,"(11x,' *** (Averages of the norms per orbital type) ***')")
        do i=1,s%nAtoms
          write(output%unit_loop,"(a,':',2x,'Ds(',i4.0,')=',f11.8,4x,'Dp(',i4.0,')=',f11.8,4x,'Dd(',i4.0,')=',f11.8)") trim(s%Types(s%Basis(i)%Material)%Name),i, singlet_coupling(1,i),i, sum(singlet_coupling(2:4,i))/3.0,i, sum(singlet_coupling(5:9,i))/5.0
        end do
    end if
    write(output%unit_loop,"(11x,' *********** Magnetization components: **********')")
    if(abs(sum(mp(:,:)))>1.e-8_dp) then
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Mx(',i4.0,')=',f11.8,4x,'My(',i4.0,')=',f11.8,4x,'Mz(',i4.0,')=',f11.8,4x,'theta = ',f11.6,4x,'phi = ',f11.6)") trim(s%Types(s%Basis(i)%Material)%Name),i,mvec_cartesian(1,i),i,mvec_cartesian(2,i),i,mvec_cartesian(3,i),mvec_spherical(2,i),mvec_spherical(3,i)
      end do
    else
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Mx(',i4.0,')=',f11.8,4x,'My(',i4.0,')=',f11.8,4x,'Mz(',i4.0,')=',f11.8)") trim(s%Types(s%Basis(i)%Material)%Name),i,mvec_cartesian(1,i),i,mvec_cartesian(2,i),i,mvec_cartesian(3,i)
      end do
    end if

    write(output%unit_loop,"(11x,' ****** Orbital components in global frame: *****')")
    if(sum(lxm(:)**2+lym(:)**2)>1.e-8_dp) then
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Lx(',i4.0,')=',f11.8,4x,'Ly(',i4.0,')=',f11.8,4x,'Lz(',i4.0,')=',f11.8,4x,'theta = ',f11.6,4x,'phi = ',f11.6)") trim(s%Types(s%Basis(i)%Material)%Name),i,lxm(i),i,lym(i),i,lzm(i),ltheta(i),lphi(i)
      end do
    else
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Lx(',i4.0,')=',f11.8,4x,'Ly(',i4.0,')=',f11.8,4x,'Lz(',i4.0,')=',f11.8)") trim(s%Types(s%Basis(i)%Material)%Name),i,lxm(i),i,lym(i),i,lzm(i)
      end do
    end if
    ! write(output%unit_loop,"(11x,' *** Orbital components in local frame:  ***')")
    ! do i=1,s%nAtoms
    !   write(output%unit_loop,"(4x,'Lxp(',i4.0,')=',f11.8,4x,'Lyp(',i4.0,')=',f11.8,4x,'Lzp(',i4.0,')=',f11.8)") i,lxpm(i),i,lypm(i),i,lzpm(i)
    !   if(sqrt(lxpm(i)**2+lypm(i)**2)/=0) &
    !   write(output%unit_loop,"(12x,'theta = ',f11.6,'  ',4x,'phi = ',f11.6)") lptheta(i),lpphi(i)
    ! end do
    write(output%unit_loop,"(11x,' ******************** Total: ********************')")
    do i=1,s%nAtoms
      write(output%unit_loop,"(a,':',2x,' N(',i4.0,')=',f11.8,4x,' M(',i4.0,')=',f11.8,4x,' L(',i4.0,')=',f11.8)") trim(s%Types(s%Basis(i)%Material)%Name),i,sum(rho(:,i)),i,mvec_spherical(1,i),i,labs(i)
    end do
    write(output%unit_loop,"('|----------===================== (',i4.0,' iterations ) =====================----------|')") iter
  end subroutine print_sc_results


  subroutine write_sc_results()
    !! Writes the self-consistency results into files and broadcasts the scfile for the next iteration.
    use mod_parameters,        only: nOrb, output, dfttype
    use EnergyIntegration,     only: parts
    use mod_magnet,            only: rho, mx, my, mz
    use mod_superconductivity, only: singlet_coupling
    use mod_system,            only: s => sys
    use mod_mpi_pars,          only: rField,MPI_CHARACTER,FieldComm,ierr
    implicit none
    character(len=30) :: formatvar
    integer           :: i,mu

    external :: MPI_Bcast
    if(rField == 0) then
      ! Writing new results (mx, my, mz and n) to file
      write(output%unit_loop,"('[write_sc_results] Writing new n, mx, my, mz and Ef to file...')")
      write(scfile,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,a,'.dat')") output%SOCchar, trim(output%Sites),dfttype,parts,trim(output%BField),trim(output%info),trim(output%SOC),trim(output%suffix)
      open (unit=99,status='replace',file=scfile)

      write(formatvar,fmt="(a,i0,a)") '(',5*nOrb,'(es21.11,2x))'

      do i=1,s%nAtoms
        write(99,fmt=formatvar) (rho(mu,i), mu=1,nOrb),(mx(mu,i), mu=1,nOrb),(my(mu,i), mu=1,nOrb),(mz(mu,i), mu=1,nOrb),(singlet_coupling(mu,i), mu=1,nOrb)
      end do
      write(99,"(es21.11,2x,'! Ef  ')") s%Ef
      write(99,"('! n(1:nOrb), mx(1:nOrb), my(1:nOrb), mz(1:nOrb), delta_sc(1:nOrb) per site ')")
      write(99,"('! Ef ')")

      close(99)
    end if

    call MPI_Bcast(scfile, len(scfile), MPI_CHARACTER, 0, FieldComm,ierr)
  end subroutine write_sc_results


  ! Writes the initial values for the self-consistency
  subroutine print_sc_step(n,mx,my,mz,singlet_coupling,Ef,fvec)
    use mod_kind,              only: dp
    use mod_parameters,        only: nOrb, output
    use mod_system,            only: s => sys
    use mod_magnet,            only: iter
    use mod_superconductivity, only: lsuperCond
    use mod_mpi_pars,          only: rField
    implicit none
    real(dp),dimension(neq),           intent(in), optional :: fvec
    real(dp),dimension(s%nAtoms),      intent(in) :: n,mx,my,mz
    real(dp)                    ,      intent(in) :: Ef
    real(dp),dimension(nOrb,s%nAtoms), intent(in) :: singlet_coupling
    integer :: i,mu
    real(dp) :: fvecsum

    if(rField==0) then
      if(present(fvec)) then
        write(output%unit_loop,"('|------------------------------- Iteration ',i0,' ------------------------------|')") iter
        do i=1,s%nAtoms
          fvecsum = 0.0
          do mu = 5,nOrb
            fvecsum = fvecsum + fvec((i-1)*neq_per_atom+(mu-4))
          end do
          if(abs(cmplx(mx(i),my(i),dp))>1.e-12_dp) then
            write(output%unit_loop,"('Site ',i4,': N(',i4,')=',es16.9,4x,'Mx(',i4,')=',es16.9,4x,'My(',i4,')=',es16.9,4x,'Mz(',i4,')=',es16.9)") i,i,n(i),i,mx(i),i,my(i),i,mz(i)
            write(output%unit_loop,"(15x,'fvec(',i4,')=',es16.9,2x,'fvec(',i4,')=',es16.9,2x,'fvec(',i4,')=',es16.9,2x,'fvec(',i4,')=',es16.9)") i,fvecsum,(i-1)*neq_per_atom+6,fvec((i-1)*neq_per_atom+6),(i-1)*neq_per_atom+7,fvec((i-1)*neq_per_atom+7),(i-1)*neq_per_atom+8,fvec((i-1)*neq_per_atom+8)
          else
            write(output%unit_loop,"('Site ',i4,': N(',i4,')=',es16.9,4x,'Mz(',i4,')=',es16.9)") i,i,n(i),i,mz(i)
            write(output%unit_loop,"(15x,'fvec(',i4,')=',es16.9,2x,'fvec(',i4,')=',es16.9)") i,fvecsum,(i-1)*neq_per_atom+8,fvec((i-1)*neq_per_atom+8)
          end if
        end do
        write(output%unit_loop,"(13x,'Ef=',es16.9)") Ef
        write(output%unit_loop,"(15x,'fvec(',i4,')=',es16.9)") neq,fvec(neq)
      else if(iter == 0) then
        write(output%unit_loop,"('|---------------- Starting charge density, magnetization and Ef ----------------|')")
        do i=1,s%nAtoms
          if(abs(cmplx(mx(i),my(i),dp))>1.e-12_dp) then
            write(output%unit_loop,"('Site ',i4,': N(',i4,')=',es16.9,4x,'Mx(',i4,')=',es16.9,4x,'My(',i4,')=',es16.9,4x,'Mz(',i4,')=',es16.9)") i,i,n(i),i,mx(i),i,my(i),i,mz(i)
          else
            write(output%unit_loop,"('Site ',i4,': N(',i4,')=',es16.9,4x,'Mz(',i4,')=',es16.9)") i,i,n(i),i,mz(i)
          end if
          if(lsuperCond) then
            write(output%unit_loop,"('Site ',i4,': Deltas =',9(es14.7,2x))") i,(singlet_coupling(mu,i),mu=1,nOrb)
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
  subroutine sc_eqs(N,x,fvec,iflag)
    use mod_kind,              only: dp
    use mod_system,            only: s => sys
    use mod_magnet,            only: iter,rho,rhod,mp,mx,my,mz,mpd,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix,           only: update_Umatrix
    use mod_parameters,        only: nOrb,lfixEf
    use mod_expectation,       only: expectation_values
    use mod_superconductivity, only: lsuperCond, update_singlet_couplings
    use mod_mpi_pars
    implicit none
    integer  :: N,i,mu,iflag
    real(dp),    dimension(N)             :: x,fvec
    real(dp),    dimension(nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    real(dp),    dimension(nOrb,s%nAtoms) :: singlet_coupling_in
    real(dp),    dimension(nOrb,s%nAtoms) :: deltas
    complex(dp), dimension(s%nAtoms)      :: mpd_in

    iflag=0

    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1,s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*neq_per_atom+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*neq_per_atom+6)
      myd_in(i) = x((i-1)*neq_per_atom+7)
      mzd_in(i) = x((i-1)*neq_per_atom+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      if(lsupercond) then
        do mu = 1,nOrb
          singlet_coupling_in(mu,i) = x((i-1)*(neq_per_atom)+8+mu)
        end do
      end if
    end do
    if(.not.lfixEf) s%Ef    = x(N)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)
    if(lsuperCond) call update_singlet_couplings(s,singlet_coupling_in)

    call print_sc_step(rhod_in,mxd_in,myd_in,mzd_in,singlet_coupling_in,s%Ef)

    iter = iter + 1

    call expectation_values(s,rho,mp,mx,my,mz,deltas)
    do i = 1,s%nAtoms
      rhod(i)   = sum(rho(5:9,i))
      mpd(i)    = sum(mp(5:9,i))
      mxd(i)    = sum(mx(5:9,i))
      myd(i)    = sum(my(5:9,i))
      mzd(i)    = sum(mz(5:9,i))
    end do
    do i = 1,s%nAtoms
      do mu = 5,nOrb
        fvec((i-1)*neq_per_atom+(mu-4)) = rho(mu,i) - rho_in(mu,i)
      end do
      fvec((i-1)*neq_per_atom+6) =  mxd(i) -  mxd_in(i)
      fvec((i-1)*neq_per_atom+7) =  myd(i) -  myd_in(i)
      fvec((i-1)*neq_per_atom+8) =  mzd(i) -  mzd_in(i)
      if(lsupercond) then
        do mu = 1,nOrb
          fvec((i-1)*neq_per_atom+8+mu) = deltas(mu,i) - singlet_coupling_in(mu,i)
        end do
      end if
    end do
    if(.not.lfixEf) fvec(N) = sum(rho) - s%totalOccupation

    call print_sc_step(rhod,mxd,myd,mzd,deltas,s%Ef,fvec)

    if(lontheflysc) call write_sc_results()

  end subroutine sc_eqs

  ! This subroutine calculates the jacobian of the system of equations
  !     n  - rho_in    = 0
  !     mx - mx_in   = 0
  !     my - my_in   = 0
  !     mz - mz_in   = 0
  !  sum n - n_total = 0
  subroutine sc_jac(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_kind,              only: dp
    use mod_system,            only: s => sys
    use mod_parameters,        only: nOrb,lfixEf
    use mod_magnet,            only: iter,rhod0,rho0,rho
    use mod_Umatrix,           only: update_Umatrix
    use mod_superconductivity, only: lsuperCond, update_singlet_couplings
    implicit none
    integer      :: N,i,mu,ldfjac,iflag
    real(dp)     :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(dp),    dimension(s%nAtoms)      :: mxd_in,myd_in,mzd_in,rhod_in
    real(dp),    dimension(nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(nOrb,s%nAtoms) :: singlet_coupling_in
    complex(dp), dimension(s%nAtoms)      :: mpd_in
    !--------------------- begin MPI vars --------------------

    iflag=0

    ! Values used in the hamiltonian
    rho_in = rho
    do i = 1,s%nAtoms
      do mu = 5,nOrb
        rho_in(mu,i) = x((i-1)*neq_per_atom+(mu-4))
      end do
      rhod_in(i)= sum(rho_in(5:9,i))
      mxd_in(i) = x((i-1)*neq_per_atom+6)
      myd_in(i) = x((i-1)*neq_per_atom+7)
      mzd_in(i) = x((i-1)*neq_per_atom+8)
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      if(lsupercond) then
        do mu = 1,nOrb
          singlet_coupling_in(mu,i) = x((i-1)*(neq_per_atom)+8+mu)
        end do
      end if
    end do
    if(.not.lfixEf) s%Ef    = x(N)

    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)
    if(lsuperCond) call update_singlet_couplings(s,singlet_coupling_in)

    fvec=fvec

    call calcJacobian_greenfunction(selfconjac, N)

    iter = iter + 1

  end subroutine sc_jac

end module mod_self_consistency
