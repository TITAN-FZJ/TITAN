module mod_self_consistency
  use mod_kind, only: dp
  implicit none
  integer               :: neq
  !! Total number of equations on the system
  integer, allocatable  :: neq_per_atom(:)
  !! Number of equations per atom (accumulated)
  character(len=300)    :: default_file
  !! Filename of results to read
  real(dp)              :: mag_tol = 1.e-12_dp
  !! Magnetic tolerance
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
#ifdef _GPU
    use mod_expectation,   only: groundstate_L_and_E,expectation_values,expectation_eigenstates_fullhk_gpu,calc_GS_L_and_E,calc_GS_L_and_E_fullhk_gpu
#else
    use mod_expectation,   only: groundstate_L_and_E,expectation_values,expectation_eigenstates_fullhk,calc_GS_L_and_E,calc_GS_L_and_E_fullhk
#endif
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
#ifdef _GPU
      expectation_values => expectation_eigenstates_fullhk_gpu
      calc_GS_L_and_E    => calc_GS_L_and_E_fullhk_gpu
#else
      expectation_values => expectation_eigenstates_fullhk
      calc_GS_L_and_E    => calc_GS_L_and_E_fullhk
#endif
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
    use mod_parameters,        only: output
    use mod_mpi_pars,          only: rField,abortProgram
    use mod_system,            only: s => sys
#ifdef _GPU
    use mod_superconductivity, only: delta_sc,delta_sc_d,lsuperCond
#else
    use mod_superconductivity, only: delta_sc,lsuperCond
#endif
    use mod_Umatrix,           only: init_Umatrix
    use mod_magnet,            only: mx,my,mz,mxd,myd,mzd,mpd,hw_count,hw_list, &
                                     lfield,rho,rhod,rhod0,rho0
    implicit none
    integer             :: i,err,mu,mud
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

      mx = 0._dp
      my = 0._dp
      mz = 0._dp
      ! Initialize hamiltonian using occupations from elemental files and rho=rho0
      rho = rho0
      atom: do i = 1, s%nAtoms
        ! Initial values for the d-quantities
        mxd(i) = initialmag(i,1)
        myd(i) = initialmag(i,2)
        mzd(i) = initialmag(i,3)
        rhod(i)= rhod0(i)

        do mud=1,s%ndOrb
          mu = s%dOrbs(mud)
          mx(mu,i) = mxd(i)/s%ndOrb
          my(mu,i) = myd(i)/s%ndOrb
          mz(mu,i) = mzd(i)/s%ndOrb
        end do
        if(lsuperCond) then
          do mu = 1,s%nOrb
            delta_sc(mu,i) = s%Types(s%Basis(i)%Material)%lambda(mu)
          end do
        end if
      end do atom
      mpd = cmplx(mxd,myd,dp)

    end if
#ifdef _GPU
    delta_sc_d = delta_sc
#endif

    call init_Umatrix(mzd,mpd,rhod,rhod0,rho,rho0,s)
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
    use mod_kind,              only: dp
    use mod_parameters,        only: output, dfttype
    use EnergyIntegration,     only: parts
    use mod_system,            only: s => sys
    use mod_superconductivity, only: delta_sc
    use mod_magnet,            only: rho, mp, mx, my, mz, rhod, mpd, mxd, myd, mzd, hw_count
    use mod_mpi_pars,          only: rField,FieldComm,ierr,MPI_DOUBLE_PRECISION
    implicit none
    character(len=300)  :: file = ""
    integer,intent(out) :: err
    logical,intent(out) :: lsuccess
    integer             :: i,j,mu,mud
    real(dp)            :: previous_results(5*s%nOrb,s%nAtoms), previous_Ef

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
          read(99,fmt=*) (previous_results(j,i), j=1,5*s%nOrb)
        end do
        read(99,fmt=*) previous_Ef
      end if

      call MPI_Bcast(previous_results,5*s%nOrb*s%nAtoms,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
      call MPI_Bcast(previous_Ef,1,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

      rho(:,:) = previous_results(         1:  s%nOrb,:)
      mx (:,:) = previous_results(  s%nOrb+1:2*s%nOrb,:)
      my (:,:) = previous_results(2*s%nOrb+1:3*s%nOrb,:)
      mz (:,:) = previous_results(3*s%nOrb+1:4*s%nOrb,:)
      mp       = cmplx(mx,my,dp)
      delta_sc(:,:) = previous_results(4*s%nOrb+1:5*s%nOrb,:)

      rhod= 0._dp
      mxd = 0._dp
      myd = 0._dp
      mzd = 0._dp
      do i=1,s%nAtoms
        do mud=1,s%ndOrb
          ! Corresponding orbital of this atom
          mu = s%dOrbs(mud)
          rhod(i) = rhod(i) + rho(mu,i)
          mxd (i) = mxd (i) + mx (mu,i)
          myd (i) = myd (i) + my (mu,i)
          mzd (i) = mzd (i) + mz (mu,i)
        end do
        mpd(i) = cmplx(mxd(i),myd(i),dp)
      end do

      call calcMagAngle()

      s%Ef = previous_Ef
      if(lsuccess) then
        err = 1   ! Read different parameters
      else
        lsuccess   = .true. ! Read same parameters (err=0)
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
    use mod_parameters,        only: output,lfixEf
    use mod_magnet,            only: rho,rhod,mxd,myd,mpd,mzd
    use mod_mpi_pars,          only: rField
    use mod_system,            only: s => sys
    use mod_dnsqe,             only: dnsqe
    use mod_superconductivity, only: lsuperCond,delta_sc
    implicit none
    real(dp), allocatable :: fvec(:),jac(:,:),wa(:),sc_solu(:)
    real(dp), allocatable :: diag(:),qtf(:)
    integer               :: i,lwa,ifail=0

    ! Setting up number of equations (total and accumulated per atom)
    allocate( neq_per_atom(s%nAtoms) )
    neq_per_atom(1) = 0
    do i = 1,s%nAtoms-1
      neq_per_atom(i+1) = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + merge(3,0,abs(s%Basis(i)%Um)>1.e-8_dp) + merge(s%nOrb,0,lsupercond)
    end do
    neq = neq_per_atom(s%nAtoms) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + merge(3,0,abs(s%Basis(i)%Um)>1.e-8_dp) + merge(s%nOrb,0,lsupercond)
    if(.not.lfixEf) neq = neq+1

    if(neq.ne.0) then
      ! Allocating variables dependent on the number of equations
      allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )

      ! Putting read n and m existing solutions into sc_solu (first guess of the subroutine)
      call set_hamiltonian_variables(.false.,s,neq,sc_solu,rho,rhod,mxd,myd,mzd,mpd,delta_sc)

      call print_sc_step(rhod,mpd,mzd,delta_sc,s)

      if(rField == 0) &
        write(output%unit_loop,"('[self_consistency] Starting self-consistency:')")

      ! Performing selfconsistency finding root of non-linear system of equations (SLATEC)
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa) )
      if(lnojac) then
        call dnsqe(sc_eqs,sc_jac,2,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs,sc_jac,1,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      end if
      ifail = ifail-1

      ! Deallocating variables dependent on the number of equations
      deallocate(sc_solu,diag,qtf,fvec,jac,wa)
    else
      if(rField == 0) &
        write(output%unit_loop,"('[self_consistency] Self-consistency not required (neq=0).')")
    end if

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
    use mod_constants,         only: cZero
    use mod_system,            only: s => sys
    use mod_magnet,            only: rho,rhot,rhod,mp,mx,my,mz,mpd,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix,           only: update_Umatrix
    use mod_expectation,       only: expectation_values
    use mod_superconductivity, only: lsuperCond,update_delta_sc
    implicit none
    integer      :: M,N,ljc,i,mu,mud,iflag,lw,liw,iw(liw)
    real(dp)     :: w(lw)
    real(dp),    dimension(N)               :: x,fvec
    real(dp),    dimension(M,ljc)           :: selfconjac
    real(dp),    dimension(s%nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(s%nAtoms)        :: mxd_in,myd_in,mzd_in,rhod_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: delta_sc_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: deltas
    complex(dp), dimension(s%nAtoms)        :: mpd_in

    w=w
    iw=iw
    iflag=iflag

    ! Values used in the hamiltonian
    rho_in = rho ! To store the non-d orbitals into rho_in
    call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mxd_in,myd_in,mzd_in,mpd_in,delta_sc_in)

    ! Update Hubbard term in Hamiltonian
    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s)
    ! Update electron-hole coupling in Hamiltonian
    if(lsuperCond) call update_delta_sc(s,delta_sc_in)

    ! call print_sc_step(rhod_in,mpd_in,mzd_in,delta_sc_in,s)
    call expectation_values(s,rho,mp,mx,my,mz,deltas)

    rhod = 0._dp
    rhot = 0._dp
    mpd  = 0._dp
    mxd  = 0._dp
    myd  = 0._dp
    mzd  = 0._dp
    do i = 1,s%nAtoms
      do mud = 1,s%ndOrb
        mu = s%dOrbs(mud)
        rhod(i) = rhod(i) + rho(mu,i)
        mxd (i) = mxd (i) + mx (mu,i)
        myd (i) = myd (i) + my (mu,i)
        mzd (i) = mzd (i) + mz (mu,i)
      end do
      rhot = rhot + sum(rho(:,i))
      mpd(i) = cmplx(mxd(i),myd(i),dp)
    end do

    ! Setting up linear system of equations:
    call set_system_of_equations(s,N,rhot,rho,mxd,myd,mzd,deltas,&
                                 rho_in,mxd_in,myd_in,mzd_in,delta_sc_in,fvec)


    call calcJacobian_greenfunction(selfconjac, N)

  end subroutine lsqfun

  subroutine calcJacobian_greenfunction(jacobian,N)
    !! Calculated the Jacobian of the spin magnetization
    use mod_kind,          only: dp,int64
    use mod_constants,     only: pi,ident_norb2,cZero,pauli_dorb,ident_dorb,cOne
    use mod_parameters,    only: eta,output,lfixEf
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
    integer     :: i,j,kounti,kountj,mu,nu,mud,nud,sigma,sigmap
    real(dp)    :: kp(3), ep
    complex(dp) :: weight
    complex(dp) :: halfUn(s%nAtoms),halfUm(s%nAtoms)
    complex(dp), dimension(s%nOrb2,s%nOrb2,4) :: pauli_a,pauli_b
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
      halfUn(i) = -0.5_dp*s%Basis(i)%Un
      halfUm(i) = -0.5_dp*s%Basis(i)%Um
    end do

    ! Build local hamiltonian
    if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

    jacobian = 0._dp

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,j,kounti,kountj,mu,mud,nud,nu,sigma,sigmap,ep,kp,weight,gf,gvg,gij,gji,temp,temp1,temp2,paulitemp) &
    !$omp& shared(llineargfsoc,llinearsoc,lfixEf,neq,local_points,s,neq_per_atom,realBZ,bzs,E_k_imag_mesh,y,eta,wght,halfUn,halfUm,pauli_a,pauli_b,jacobian)
    allocate( temp1(s%nOrb2,s%nOrb2,4), &
              temp2(s%nOrb2,s%nOrb2,4), &
              gij(s%nOrb2,s%nOrb2), gji(s%nOrb2,s%nOrb2), &
              gf (s%nOrb2,s%nOrb2, s%nAtoms, s%nAtoms), &
              temp(s%nOrb2,s%nOrb2), paulitemp(s%nOrb2,s%nOrb2), stat = AllocateStatus)
    if (AllocateStatus/=0) &
      call abortProgram("[calcJacobian_greenfunction] Not enough memory for: temp1, temp2, gij, gji, gf, temp")
    gf        = cZero
    temp      = cZero

    if(llineargfsoc .or. llinearsoc) then
      allocate(gvg(s%nOrb2,s%nOrb2,s%nAtoms,s%nAtoms), STAT = AllocateStatus  )
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
        kountj = neq_per_atom(j)
        do i=1,s%nAtoms

          ! First product: temp1 =  pauli*g_ij
          gij = gf(:,:,i,j)
          do sigma = 1,4
            paulitemp = pauli_a(:,:, sigma)
            call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,cOne,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)
            temp1(:,:,sigma) = temp
          end do

          if(abs(s%Basis(j)%Un)>1.e-8_dp) then
            do nu=1,s%ndOrb
              ! Restarting line counter
              kounti = neq_per_atom(i)

              nud = s%dOrbs(nu)
              mud = nud+s%nOrb
              ! Second product: temp2 = (U_j\delta_{mu,nu} -U_j/2) * sigma * g_ji
              gji = gf(:,:,j,i)
              paulitemp = pauli_b(:,:, 1)
              paulitemp(nud,nud) = paulitemp(nud,nud) - 2._dp
              paulitemp(mud,mud) = paulitemp(mud,mud) - 2._dp
              call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,halfUn(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
              temp2(:,:,1) = temp

              ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

              ! Charge density-charge density part
              gij = temp1(:,:,1)
              gji = temp2(:,:,1)
              call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

              if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                do mu=1,s%ndOrb
                  mud = s%dOrbs(mu)
                  jacobian(kounti+mu,kountj+nu) = jacobian(kounti+mu,kountj+nu) + real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                end do
                kounti = neq_per_atom(i) + s%ndOrb
              end if

              ! Last line (Total charge neutrality)
              if(.not.lfixEf) then
                do mu = 1,s%nOrb
                   jacobian(neq,kountj+nu) = jacobian(neq,kountj+nu) + real(temp(mu,mu) + temp(mu+s%nOrb,mu+s%nOrb))
                end do
              end if

              ! Magnetic density-charge density part
              if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                do sigma = 2,4
                  gij = temp1(:,:,sigma)
                  gji = temp2(:,:,1)
                  call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                  do mu=1,s%ndOrb
                    mud = s%dOrbs(mu)
                    jacobian(kounti+sigma-1,kountj+nu) = jacobian(kounti+sigma-1,kountj+nu) + real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                  end do
                end do
                kounti = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
              end if

            end do

            kountj = neq_per_atom(j) + s%ndOrb
          end if ! |Un(j)| > 0

          if(abs(s%Basis(j)%Um)>1.e-8_dp) then

            gji = gf(:,:,j,i)

            do sigmap = 2,4
              ! Second product: temp2 = (-U/2) * sigma* g_ji
              paulitemp = pauli_b(:,:,sigmap)
              call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,halfUm(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
              temp2(:,:,sigmap) = temp
            end do

            ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

            ! Charge density-magnetic density part
            do sigmap = 2,4
              ! Restarting line counter
              kounti = neq_per_atom(i)

              gij = temp1(:,:,1)
              gji = temp2(:,:,sigmap)
              call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

              if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                do mu = 1,s%ndOrb
                  mud = s%dOrbs(mu)
                  jacobian(kounti+mu,kountj+sigmap-1) = jacobian(kounti+mu,kountj+sigmap-1) + real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                end do
                kounti = neq_per_atom(i) + s%ndOrb
              end if

              ! Last line (Total charge neutrality)
              if(.not.lfixEf) then
                do mu = 1,s%nOrb
                  jacobian(neq,kountj+sigmap-1) = jacobian(neq,kountj+sigmap-1) + real(temp(mu,mu) + temp(mu+s%nOrb,mu+s%nOrb))
                end do
              end if

              ! Magnetic density-magnetic density part
              if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                do sigma = 2,4
                  gij = temp1(:,:,sigma)
                  gji = temp2(:,:,sigmap)
                  call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                  do mu=1,s%ndOrb
                    mud = s%dOrbs(mu)
                    jacobian(kounti+sigma-1,kountj+sigmap-1) = jacobian(kounti+sigma-1,kountj+sigmap-1) + real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                  end do
                end do
                kounti = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
              end if
            end do
            kountj = neq_per_atom(j) + merge(s%ndOrb,0,abs(s%Basis(j)%Un)>1.e-8_dp) + 3
          end if ! |Um(j)| > 0

          ! **** Removing non-linear (quadratic) terms: ****
          if(llineargfsoc .or. llinearsoc) then
            ! Restarting column counter
            kountj = neq_per_atom(j)

            gij = gvg(:,:,i,j)

            do sigma = 1,4
              ! First product: temp1 =  pauli*g_ij
              paulitemp = pauli_a(:,:, sigma)
              call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,cOne,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)
              temp1(:,:, sigma) = temp
            end do

            if(abs(s%Basis(j)%Un)>1.e-8_dp) then
              do nu=1,s%ndOrb
                ! Restarting line counter
                kounti = neq_per_atom(i)

                nud = s%dOrbs(nu)
                mud = nud+s%nOrb

                gji = gvg(:,:,j,i)
                ! Second product: temp2 = (-U/2) * sigma* g_ji
                paulitemp = pauli_b(:,:, 1)
                paulitemp(nud,nud) = paulitemp(nud,nud) - 2._dp
                paulitemp(mud,mud) = paulitemp(mud,mud) - 2._dp
                call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,halfUn(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
                temp2(:,:, 1) = temp

                ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

                ! Charge density-charge density part
                gij = temp1(:,:, 1)
                gji = temp2(:,:, 1)
                call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                  do mu=1,s%ndOrb
                    mud = s%dOrbs(mu)
                    jacobian(kounti+mu,kountj+nu) = jacobian(kounti+mu,kountj+nu) - real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                  end do
                  kounti = neq_per_atom(i) + s%ndOrb
                end if

                ! Last line (Total charge neutrality)
                if(.not.lfixEf) then
                  do mu = 1,s%nOrb
                    jacobian(neq,kountj+nu) = jacobian(neq,kountj+nu) - real(temp(mu,mu) + temp(mu+s%nOrb,mu+s%nOrb))
                  end do
                end if

                ! Magnetic density-charge density part
                if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                  do sigma = 2,4
                    gij = temp1(:,:, sigma)
                    gji = temp2(:,:, 1)
                    call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                    do mu=1,s%ndOrb
                      mud = s%dOrbs(mu)
                      jacobian(kounti+sigma-1,kountj+nu) = jacobian(kounti+sigma-1,kountj+nu) - real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                    end do
                  end do
                  kounti = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
                end if

              end do
              kountj = neq_per_atom(j) + s%ndOrb
            end if ! |Un(j)| > 0

            if(abs(s%Basis(j)%Um)>1.e-8_dp) then

              gji = gvg(:,:,j,i)

              do sigmap = 2,4
                ! Second product: temp2 = (-U/2) * sigma* g_ji
                paulitemp = pauli_b(:,:, sigmap)
                call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,halfUm(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
                temp2(:,:, sigmap) = temp
              end do

              ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

              ! Charge density-magnetic density part
              do sigmap = 2,4
                ! Restarting line counter
                kounti = neq_per_atom(i)

                gij = temp1(:,:, 1)
                gji = temp2(:,:, sigmap)
                call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                  do mu = 1,s%ndOrb
                    mud = s%dOrbs(mu)
                    jacobian(kounti+mu,kountj+sigmap-1) = jacobian(kounti+mu,kountj+sigmap-1) - real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                  end do
                  kounti = neq_per_atom(i) + s%ndOrb
                end if

                ! Last line (Total charge neutrality)
                if(.not.lfixEf) then
                  do mu = 1,s%nOrb
                    jacobian(neq,kountj+sigmap-1) = jacobian(neq,kountj+sigmap-1) - real(temp(mu,mu) + temp(mu+s%nOrb,mu+s%nOrb))
                  end do
                end if

                ! Magnetic density-magnetic density part
                if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                  do sigma = 2,4
                    gij = temp1(:,:, sigma)
                    gji = temp2(:,:, sigmap)
                    call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                    do mu=1,s%ndOrb
                      mud = s%dOrbs(mu)
                      jacobian(kounti+sigma-1,kountj+sigmap-1) = jacobian(kounti+sigma-1,kountj+sigmap-1) - real(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
                    end do
                  end do
                  kounti = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
                end if
              end do
              kountj = neq_per_atom(j) + merge(s%ndOrb,0,abs(s%Basis(j)%Un)>1.e-8_dp) + 3
            end if ! |Um(j)| > 0

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
          ! Restarting line counter
          kounti = neq_per_atom(i)

          gij = gf(:,:,i,i)

          ! Charge density per orbital lines
          ! temp1 =  pauli*g_ii
          paulitemp = pauli_a(:,:, 1)
          call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)


          if(abs(s%Basis(i)%Un)>1.e-8_dp) then
            do mu = 1,s%ndOrb
              mud = s%dOrbs(mu)
              jacobian(kounti+mu,neq) = jacobian(kounti+mu,neq) - aimag(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
            end do
            kounti = neq_per_atom(i) + s%ndOrb
          end if

          ! Last line (Total charge neutrality)
          if(.not.lfixEf) then
            do mu = 1,s%nOrb
              jacobian(neq,neq) = jacobian(neq,neq)  - aimag(temp(mu,mu) + temp(mu+s%nOrb,mu+s%nOrb))
            end do
          end if


          ! Magnetic density-last column part
          if(abs(s%Basis(i)%Um)>1.e-8_dp) then
            do sigma = 2,4
              ! temp1 =  pauli*g_ii
              paulitemp = pauli_a(:,:, sigma)
              call zgemm('n','n',s%nOrb2,s%nOrb2,s%nOrb2,weight,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)

              do mu=1,s%ndOrb
                mud = s%dOrbs(mu)
                jacobian(kounti+sigma-1,neq) = jacobian(kounti+sigma-1,neq) - dimag(temp(mud,mud) + temp(mud+s%nOrb,mud+s%nOrb))
              end do
            end do
            kounti = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
          end if

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
    do i = 1, neq-1 
      jacobian(i,i) = jacobian(i,i) - 1._dp
    end do
  end subroutine calcJacobian_greenfunction

  subroutine rotate_magnetization_to_field()
  !! Rotate the magnetization to the direction of the field (useful for SOC=F)
    use mod_kind,       only: dp
    use mod_constants,  only: deg2rad
    use mod_magnet,     only: lfield,hw_count,hw_list,hhw,mp,mx,my,mz,mpd,mxd,myd,mzd
    use mod_parameters, only: output
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: rField
    implicit none
    integer      :: i,mu,mud,signal
    real(dp) :: mdotb,mabs(s%nOrb,s%nAtoms)

    if(rField == 0) &
      write(output%unit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    if(.not.lfield) then
      if(rField == 0) &
        write(output%unit_loop,"('[rotate_magnetization_to_field] Field if OFF! No rotation is done.')")
      return
    end if

    do i = 1, s%nAtoms
      do mu=1,s%nOrb
        mdotb   = hhw(1,i)*mx(mu,i)+hhw(2,i)*my(mu,i)+hhw(3,i)*mz(mu,i)
        signal  = int(mdotb/abs(mdotb))
        mabs(mu,i) = sqrt((mx(mu,i)**2)+(my(mu,i)**2)+(mz(mu,i)**2))
        mx  (mu,i) = signal*mabs(mu,i)*sin(hw_list(hw_count,2)*deg2rad)*cos(hw_list(hw_count,3)*deg2rad)
        my  (mu,i) = signal*mabs(mu,i)*sin(hw_list(hw_count,2)*deg2rad)*sin(hw_list(hw_count,3)*deg2rad)
        mz  (mu,i) = signal*mabs(mu,i)*cos(hw_list(hw_count,2)*deg2rad)
        mp  (mu,i) = cmplx(mx(mu,i),my(mu,i),dp)
      end do
      do mud=1, s%ndOrb
        ! Corresponding orbital of this atom
        mu = s%dOrbs(mud)
        mxd(i)  = mxd(i)  + mx(mu,i)
        myd(i)  = myd(i)  + my(mu,i)
        mzd(i)  = mzd(i)  + mz(mu,i)
      end do
      mpd(i) = cmplx(mxd(i),myd(i),dp)
    end do

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
    use mod_magnet,            only: rho,rhos,rhop,rhod,mvec_cartesian,mp,mvec_spherical, &
                                     lxm,lym,lzm,ltheta,lphi,labs,iter
    use mod_superconductivity, only: lsuperCond,delta_sc
    implicit none
    real(dp), dimension(s%nAtoms) :: deltas,deltap,deltad
    integer  :: i,mu

    write(output%unit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    if(leigenstates) then
      write(output%unit_loop,"(14x,'Ef=',f10.7,4x,'Eband=',f15.7)") s%Ef,energy
    else
      write(output%unit_loop,"(28x,'Ef=',f10.7)") s%Ef
    end if

    write(output%unit_loop,"(11x,' *************** Charge density: ****************')")
    rhos = 0._dp
    rhop = 0._dp
    do i=1,s%nAtoms
      ! Total s-orbital occupation
      do mu=1,s%nsOrb
        rhos(i) = rhos(i) + rho(s%sOrbs(mu),i)
      end do
      ! Total p-orbital occupation
      do mu=1,s%npOrb
        rhop(i) = rhop(i) + rho(s%pOrbs(mu),i)
      end do
      write(output%unit_loop,fmt="(a,':',2x,'Ns=',f10.7,4x,'Np=',f10.7,4x,'Nd=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),rhos(i),rhop(i),rhod(i)
    end do

    if(lsupercond) then
      write(output%unit_loop,"(11x,' ******** Superconducting gap parameter: ********')")
      write(output%unit_loop,"(11x,' *** (Averages of the norms per orbital type) ***')")
      deltas = 0._dp
      deltap = 0._dp
      deltad = 0._dp
      do i=1,s%nAtoms
        ! Average s-orbital superconducting delta
        do mu=1,s%nsOrb
          deltas(i) = deltas(i) + delta_sc(s%sOrbs(mu),i)
        end do
        deltas(i) = deltas(i)/s%nsOrb
        ! Average p-orbital superconducting delta
        do mu=1,s%npOrb
          deltap(i) = deltap(i) + delta_sc(s%pOrbs(mu),i)
        end do
        deltap(i) = deltap(i)/s%npOrb
        ! Average d-orbital superconducting delta
        do mu=1,s%ndOrb
          deltad(i) = deltad(i) + delta_sc(s%dOrbs(mu),i)
        end do
        deltad(i) = deltad(i)/s%ndOrb
        write(output%unit_loop,fmt="(a,':',2x,'Ds=',f10.7,4x,'Dp=',f10.7,4x,'Dd=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),deltas(i),deltap(i),deltad(i)
      end do
    end if

    write(output%unit_loop,"(11x,' *********** Magnetization components: **********')")
    if(abs(sum(mp(:,:)))>1.e-7_dp) then
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Mx=',f10.7,4x,'My=',f10.7,4x,'Mz=',f10.7,4x,'theta = ',f10.5,4x,'phi = ',f10.5)") trim(s%Types(s%Basis(i)%Material)%Name),mvec_cartesian(1,i),mvec_cartesian(2,i),mvec_cartesian(3,i),mvec_spherical(2,i),mvec_spherical(3,i)
      end do
    else
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Mx=',f10.7,4x,'My=',f10.7,4x,'Mz=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),mvec_cartesian(1,i),mvec_cartesian(2,i),mvec_cartesian(3,i)
      end do
    end if

    write(output%unit_loop,"(11x,' ****** Orbital components in global frame: *****')")
    if(sum(lxm(:)**2+lym(:)**2)>1.e-7_dp) then
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Lx=',f10.7,4x,'Ly=',f10.7,4x,'Lz=',f10.7,4x,'theta = ',f10.5,4x,'phi = ',f10.5)") trim(s%Types(s%Basis(i)%Material)%Name),lxm(i),lym(i),lzm(i),ltheta(i),lphi(i)
      end do
    else
      do i=1,s%nAtoms
        write(output%unit_loop,"(a,':',2x,'Lx=',f10.7,4x,'Ly=',f10.7,4x,'Lz=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),lxm(i),lym(i),lzm(i)
      end do
    end if
    ! write(output%unit_loop,"(11x,' *** Orbital components in local frame:  ***')")
    ! do i=1,s%nAtoms
    !   write(output%unit_loop,"(4x,'Lxp(',i4.0,')=',f10.7,4x,'Lyp(',i4.0,')=',f10.7,4x,'Lzp(',i4.0,')=',f10.7)") i,lxpm(i),i,lypm(i),i,lzpm(i)
    !   if(sqrt(lxpm(i)**2+lypm(i)**2)/=0) &
    !   write(output%unit_loop,"(12x,'theta = ',f11.6,'  ',4x,'phi = ',f11.6)") lptheta(i),lpphi(i)
    ! end do
    write(output%unit_loop,"(11x,' ******************** Total: ********************')")
    do i=1,s%nAtoms
      write(output%unit_loop,"(a,':',2x,' N=',f10.7,2x,' |M|=',f10.7,2x,' |L|=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),sum(rho(:,i)),mvec_spherical(1,i),labs(i)
    end do
    write(output%unit_loop,"('|----------===================== (',i4.0,' iterations ) =====================----------|')") iter
  end subroutine print_sc_results


  subroutine write_sc_results()
    !! Writes the self-consistency results into files and broadcasts the scfile for the next iteration.
    use mod_parameters,        only: output,dfttype
    use EnergyIntegration,     only: parts
    use mod_magnet,            only: rho, mx, my, mz
    use mod_superconductivity, only: delta_sc
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

      write(formatvar,fmt="(a,i0,a)") '(',5*s%nOrb,'(es21.11,2x))'

      do i=1,s%nAtoms
        write(99,fmt=formatvar) (rho(mu,i), mu=1,s%nOrb),(mx(mu,i), mu=1,s%nOrb),(my(mu,i), mu=1,s%nOrb),(mz(mu,i), mu=1,s%nOrb),(delta_sc(mu,i), mu=1,s%nOrb)
      end do
      write(99,"(es21.11,2x,'! Ef  ')") s%Ef
      write(99,"('! n(1:nOrb), mx(1:nOrb), my(1:nOrb), mz(1:nOrb), delta_sc(1:nOrb) per site ')")
      write(99,"('! Ef ')")

      close(99)
    end if

    call MPI_Bcast(scfile, len(scfile), MPI_CHARACTER, 0, FieldComm,ierr)
  end subroutine write_sc_results


  ! Writes the initial values for the self-consistency
  subroutine print_sc_step(n,mp,mz,delta_sc,s,fvec)
    use mod_kind,              only: dp
    use mod_parameters,        only: output,lfixEf
    use mod_System,            only: System_type
    use mod_magnet,            only: iter
    use mod_superconductivity, only: lsuperCond
    use mod_mpi_pars,          only: rField
    use mod_tools,             only: RtoS
    implicit none
    type(System_type),                      intent(in) :: s
    real(dp),   dimension(s%nAtoms),        intent(in) :: n,mz
    complex(dp),dimension(s%nAtoms),        intent(in) :: mp
    real(dp),   dimension(s%nOrb,s%nAtoms), intent(in) :: delta_sc
    real(dp),   dimension(neq),             intent(in), optional :: fvec
    character(:), allocatable :: formatvar
    integer  :: i,mu,kount
    real(dp) :: fvecsum

    if(rField/=0) return
    if(present(fvec)) then
      write(output%unit_loop,"('|------------------------------- Iteration ',i0,' ------------------------------|')") iter
      do i=1,s%nAtoms
        formatvar = ""
        if(abs(mp(i))>1.e-8_dp) then
          write(output%unit_loop,"('Site ',i4,': N=',es16.9,9x,'Mx=',es16.9,4x,'My=',es16.9,4x,'Mz=',es16.9)") i,n(i),mp(i)%re,mp(i)%im,mz(i)
        else
          write(output%unit_loop,"('Site ',i4,': N=',es16.9,9x,'Mz=',es16.9)") i,n(i),mz(i)
        end if
        fvecsum = 0._dp
        if(abs(s%Basis(i)%Un)>1.e-8_dp) then
          do mu = 1,s%ndOrb
            fvecsum = fvecsum + fvec(neq_per_atom(i)+mu)
          end do
          formatvar = formatvar // '  fvec(N)= ' // trim(RtoS(fvecsum,"(es16.9)"))
        end if
        kount = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp)
        if(abs(s%Basis(i)%Um)>1.e-8_dp) then
          if(abs(mp(i))>1.e-8_dp) then
            formatvar = formatvar // '  fvec(Mx)= ' // trim(RtoS(fvec(kount+1),"(es16.9)")) // '  fvec(My)= ' // trim(RtoS(fvec(kount+2),"(es16.9)")) // '  fvec(Mz)= ' // trim(RtoS(fvec(kount+3),"(es16.9)"))
          else
            formatvar = formatvar // '  fvec(Mz)= ' // trim(RtoS(fvec(kount+3),"(es16.9)"))
          end if
        end if
        if(len_trim(formatvar)/=0) write(output%unit_loop,"(3x,a)") formatvar
      end do
      write(output%unit_loop,"(13x,'Ef=',es16.9)") s%Ef
      if(.not.lfixEf) write(output%unit_loop,"(8x,'fvec(',i0,')=',es16.9)") neq,fvec(neq)
    else if(iter == 0) then
      write(output%unit_loop,"('|---------------- Starting charge density, magnetization and Ef ----------------|')")
      do i=1,s%nAtoms
        if(abs(mp(i))>1.e-8_dp) then
          write(output%unit_loop,"('Site ',i4,': N=',es16.9,4x,'Mx=',es16.9,4x,'My=',es16.9,4x,'Mz=',es16.9)") i,n(i),mp(i)%re,mp(i)%im,mz(i)
        else
          write(output%unit_loop,"('Site ',i4,': N=',es16.9,4x,'Mz=',es16.9)") i,n(i),mz(i)
        end if
        if(lsuperCond) then
          write(output%unit_loop,"('Site ',i4,': D_sc =',9(es14.7,2x))") i,(delta_sc(mu,i),mu=1,s%nOrb)
        end if
      end do
      write(output%unit_loop,"(13x,'Ef=',es16.9)") s%Ef
    end if

  end subroutine print_sc_step


  subroutine sc_eqs(N,x,fvec,iflag)
  !! This subroutine calculates the self-consistency equations
  !!     n  - rho_in    = 0
  !!     mx - mx_in   = 0
  !!     my - my_in   = 0
  !!     mz - mz_in   = 0
  !!  sum n - n_total = 0
    use mod_kind,              only: dp
    use mod_system,            only: s => sys
    use mod_magnet,            only: iter,rho,rhot,rhod,mp,mx,my,mz,mpd,mxd,myd,mzd,rhod0,rho0
    use mod_Umatrix,           only: update_Umatrix
    use mod_expectation,       only: expectation_values
    use mod_superconductivity, only: lsuperCond, update_delta_sc
    implicit none
    integer  :: N,i,mu,mud,iflag
    real(dp),    dimension(N)               :: x,fvec
    real(dp),    dimension(s%nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(s%nAtoms)        :: mxd_in,myd_in,mzd_in,rhod_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: delta_sc_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: deltas
    complex(dp), dimension(s%nAtoms)        :: mpd_in

    iflag=0

    ! Values used in the hamiltonian
    rho_in = rho ! To store the non-d orbitals into rho_in
    call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mxd_in,myd_in,mzd_in,mpd_in,delta_sc_in)

    ! Update Hubbard term in Hamiltonian
    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s)
    ! Update electron-hole coupling in Hamiltonian
    if(lsuperCond) call update_delta_sc(s,delta_sc_in)


    ! call print_sc_step(rhod_in,mpd_in,mzd_in,delta_sc_in,s)

    iter = iter + 1

    call expectation_values(s,rho,mp,mx,my,mz,deltas)

    rhod = 0._dp
    mpd  = 0._dp
    mxd  = 0._dp
    myd  = 0._dp
    mzd  = 0._dp
    do i = 1,s%nAtoms
      do mud = 1,s%ndOrb
        mu = s%dOrbs(mud)
        rhod(i) = rhod(i) + rho(mu,i)
        mxd (i) = mxd (i) + mx (mu,i)
        myd (i) = myd (i) + my (mu,i)
        mzd (i) = mzd (i) + mz (mu,i)
      end do
      mpd(i) = cmplx(mxd(i),myd(i),dp)
    end do
    rhot = sum(rho(:,:))

    ! Setting up linear system of equations:
    call set_system_of_equations(s,N,rhot,rho,mxd,myd,mzd,deltas,&
                                 rho_in,mxd_in,myd_in,mzd_in,delta_sc_in,fvec)

    call print_sc_step(rhod,mpd,mzd,deltas,s,fvec)

    if(lontheflysc) call write_sc_results()

  end subroutine sc_eqs


  subroutine sc_jac(N,x,fvec,selfconjac,ldfjac,iflag)
  !! This subroutine calculates the jacobian of the system of equations
  !!     n  - rho_in    = 0
  !!     mx - mx_in   = 0
  !!     my - my_in   = 0
  !!     mz - mz_in   = 0
  !!  sum n - n_total = 0
    use mod_kind,              only: dp
    use mod_system,            only: s => sys
    use mod_magnet,            only: iter,rhod0,rho0,rho
    use mod_Umatrix,           only: update_Umatrix
    use mod_superconductivity, only: lsuperCond, update_delta_sc
    implicit none
    integer      :: N,ldfjac,iflag
    real(dp)     :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(dp),    dimension(s%nAtoms)        :: mxd_in,myd_in,mzd_in,rhod_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: delta_sc_in
    complex(dp), dimension(s%nAtoms)        :: mpd_in
    !--------------------- begin MPI vars --------------------

    iflag=0

    ! Values used in the hamiltonian
    rho_in = rho ! To store the non-d orbitals into rho_in
    call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mxd_in,myd_in,mzd_in,mpd_in,delta_sc_in)

    ! Update Hubbard term in Hamiltonian
    call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s)
    ! Update electron-hole coupling in Hamiltonian
    if(lsuperCond) call update_delta_sc(s,delta_sc_in)

    fvec=fvec

    call calcJacobian_greenfunction(selfconjac, N)

    iter = iter + 1

  end subroutine sc_jac

  subroutine set_hamiltonian_variables(set,s,N,x,rho,rhod,mxd,myd,mzd,mpd,delta_sc)
  !! This subroutine transfers the variables from the single array
  !! to the expectation values used in the hamiltonian
  !! Whe set = .false. it does the inverse, to set the initial guess
    use mod_constants,         only: cZero
    use mod_System,            only: System_type
    use mod_superconductivity, only: lsupercond
    use mod_parameters,        only: lfixEf
    implicit none
    logical,                                 intent(in)    :: set
    !! Variable to select which way variables are set
    type(System_type),                       intent(inout) :: s
    !! System quantities
    integer,                                 intent(in)    :: N
    !! Number of equations
    real(dp),    dimension(N),               intent(inout) :: x
    !! Unknowns in single array format
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(inout) :: rho
    !! Variables used in the hamiltonian: orbital-dependent density
    real(dp),    dimension(s%nAtoms),        intent(inout) :: mxd,myd,mzd,rhod
    !! Variables used in the hamiltonian: magnetizatin and density
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(inout) :: delta_sc
    !! Variables used in the hamiltonian: superconducting delta
    complex(dp), dimension(s%nAtoms),        intent(inout) :: mpd
    !! Variables used in the hamiltonian: circular component of magnetic moment

    ! Local variables:
    integer :: i,mu,kount

    if(set) then
      rhod = 0._dp
      mzd  = 0._dp
      mpd  = cZero
      do i = 1,s%nAtoms
        kount = neq_per_atom(i)
        if(abs(s%Basis(i)%Un)>1.e-8_dp) then
          do mu = 1,s%ndOrb
            rho(s%dOrbs(mu),i) = x(kount+mu)
            rhod(i)= rhod(i) + rho(s%dOrbs(mu),i)
          end do
          kount = neq_per_atom(i) + s%ndOrb
        end if
        if(abs(s%Basis(i)%Um)>1.e-8_dp) then
          mxd(i) = x(kount+1)
          myd(i) = x(kount+2)
          mzd(i) = x(kount+3)
          mpd(i) = cmplx(mxd(i),myd(i),dp)
          kount = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
        end if
        if(lsupercond) then
          do mu = 1,s%nOrb
            delta_sc(mu,i) = x(kount+mu)
          end do
        end if
      end do
      if(.not.lfixEf) s%Ef = x(N)
    else
      do i = 1, s%nAtoms
        kount = neq_per_atom(i)
        if(abs(s%Basis(i)%Un)>1.e-8_dp) then
          do mu = 1,s%ndOrb
            x(kount+mu) = rho(s%dOrbs(mu),i)
          end do
          kount = neq_per_atom(i) + s%ndOrb
        end if
        if(abs(s%Basis(i)%Um)>1.e-8_dp) then
          x(kount+1) = mxd(i)
          x(kount+2) = myd(i)
          x(kount+3) = mzd(i)
          kount = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
        end if
        if(lsuperCond) then
          do mu = 1, s%nOrb
            x(kount+mu) = delta_sc(mu,i)
          end do
        end if
      end do
      if(.not.lfixEf) x(neq) = s%Ef
    end if

  end subroutine set_hamiltonian_variables


  subroutine set_system_of_equations(s,N,rhot,rho,mxd,myd,mzd,delta_sc,&
                                     rho_in,mxd_in,myd_in,mzd_in,delta_sc_in,fvec)
  !! This subroutine sets the system of equations to be solved
  !! by the non-linear root finder
    use mod_System,            only: System_type
    use mod_superconductivity, only: lsupercond
    use mod_parameters,        only: lfixEf
    implicit none
    type(System_type),                       intent(in)  :: s
    !! System quantities
    integer,                                 intent(in)  :: N
    !! Number of equations
    real(dp),                                intent(in)  :: rhot
    !! Total density
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(in)  :: rho,rho_in
    !! Orbital-dependent densities (output and input)
    real(dp),    dimension(s%nAtoms),        intent(in)  :: mxd,myd,mzd,mxd_in,myd_in,mzd_in
    !! Variables used in the hamiltonian: magnetizatin and density
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(in)  :: delta_sc,delta_sc_in
    !! Variables used in the hamiltonian: superconducting delta
    real(dp),    dimension(N),               intent(out) :: fvec
    !! Vector holding the equations to find the zeroes

    ! Local variables:
    integer :: i,mu,mud,kount

    do i = 1,s%nAtoms
      kount = neq_per_atom(i)
      if(abs(s%Basis(i)%Un)>1.e-8_dp) then
        do mud = 1,s%ndOrb
          mu = s%dOrbs(mud)
          fvec(kount+mud) = rho(mu,i) - rho_in(mu,i)
        end do
        kount = neq_per_atom(i) + s%ndOrb
      end if
      if(abs(s%Basis(i)%Um)>1.e-8_dp) then
        fvec(kount+1) = mxd(i) - mxd_in(i)
        fvec(kount+2) = myd(i) - myd_in(i)
        fvec(kount+3) = mzd(i) - mzd_in(i)
        kount = neq_per_atom(i) + merge(s%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
      end if
      if(lsuperCond) then
        do mu = 1, s%nOrb
          fvec(kount+mu) = delta_sc(mu,i) - delta_sc_in(mu,i)
        end do
      end if
    end do
    if(.not.lfixEf) fvec(N) = rhot - s%totalOccupation
  end subroutine set_system_of_equations

end module mod_self_consistency
