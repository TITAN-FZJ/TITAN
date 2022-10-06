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
  character(len=500)    :: mixing = ""
  !! Type of mixing used in self-consistency
  character(len=200)    :: scfile = ""
  !! Give a file to start self-consistency
  logical               :: skipsc
  !! Skip self-consistency
  character(len=50)     :: magbasis = ""
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
    use mod_magnet,        only: constr_type,lp_matrix,mtheta,mphi,rho,mx,my,mz,bc
    use adaptiveMesh,      only: genLocalEKMesh,freeLocalEKMesh
    use mod_BrillouinZone, only: realBZ
    use mod_mpi_pars,      only: rFreq, sFreq, FreqComm, rField, sField, FieldComm
    use mod_parameters,    only: leigenstates,lkpoints,lEf_overwrite,Ef_overwrite
    use mod_System,        only: s => sys
#ifdef _GPU
    use mod_expectation,   only: groundstate_L_and_E,expectation_values,expectation_eigenstates_fullhk_gpu,calc_GS_L_and_E,calc_GS_L_and_E_fullhk_gpu,calc_E_dc
#else
    use mod_expectation,   only: groundstate_L_and_E,expectation_values,expectation_eigenstates_fullhk,calc_GS_L_and_E,calc_GS_L_and_E_fullhk,calc_E_dc
#endif
    use mod_hamiltonian,   only: fullhamiltk,lfullhk
    use adaptiveMesh,      only: bzs
    implicit none
    logical :: lsuccess = .false.

    ! Distribute Energy Integration across all points available
    call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1),lkpoints)
    if(.not.leigenstates) then
      call genLocalEKMesh(s,rField,sField, FieldComm,bzs)
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
    if(lselfcon) then
      select case(constr_type)
        case(1)
          !Transversal constraining field in two loops
          call calcSelfConsistency_trans_constr()
        ! case for linear mixing.. not implemented yet
        !case("linear")
        !  call calcSelfConsistency_linear()
          if(.not.lontheflysc) call write_sc_results(rho,mx,my,mz,bc)
        case(2)
          call calcSelfConsistency()
          if(.not.lontheflysc) call write_sc_results(rho,mx,my,mz,bc)
        case default
          call calcSelfConsistency()
          if(.not.lontheflysc) call write_sc_results(rho,mx,my,mz)
      end select
    end if

    ! L matrix in local frame for given quantization direction
    call lp_matrix(mtheta, mphi)

    ! Calculating ground state Orbital Angular Momentum and Band Energy
    call groundstate_L_and_E()

    ! Calculate double counting contributions to the energy
    call calc_E_dc()

    ! Writing self-consistency results on screen
    if(rField == 0)  call print_sc_results()

    if(.not.leigenstates) call freeLocalEKMesh()

  end subroutine doSelfConsistency

  subroutine read_previous_results(lsuccess)
  !> Tries to read n and m if available
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
    use mod_magnet,            only: lconstraining_field,mx,my,mz,mabsd,m_fix,m_fix_abs,mxd,myd,mzd,mzd0,mpd,mpd0,hw_count,hw_list, &
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
    else !  If file does not exist
      if(rField == 0) then
        write(output%unit_loop,"('[read_previous_results] Self-consistency file does not exist:')")
        write(output%unit_loop,"('[read_previous_results] ',a)") trim(default_file)
      end if
      lselfcon = .true.

      ! Starting magnetization
      allocate(initialmag(s%nAtoms,3))
      if(trim(magbasis)/="") then
        call read_initial_values("initialmag",trim(magbasis),s%nAtoms,initialmag)
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
        mabsd(i) = sqrt(mxd(i)**2+myd(i)**2+mzd(i)**2)
        rhod(i)= rhod0(i)

        do mud=1,s%Types(s%Basis(i)%Material)%ndOrb
          mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
          mx(mu,i) = mxd(i)/s%Types(s%Basis(i)%Material)%ndOrb
          my(mu,i) = myd(i)/s%Types(s%Basis(i)%Material)%ndOrb
          mz(mu,i) = mzd(i)/s%Types(s%Basis(i)%Material)%ndOrb
        end do
        if(lsuperCond) then
          do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
            delta_sc(mu,i) = s%Types(s%Basis(i)%Material)%lambda(mu)
          end do
        end if
      end do atom
      mpd = cmplx(mxd,myd,dp)

      if(lconstraining_field) then
        do i=1,s%nAtoms
          m_fix(1,i) = mxd(i)/mabsd(i)
          m_fix(2,i) = myd(i)/mabsd(i)
          m_fix(3,i) = mzd(i)/mabsd(i)
          m_fix_abs(i) = mabsd(i)
        end do
      end if

    end if
#ifdef _GPU
    delta_sc_d = delta_sc
#endif

    call init_Umatrix(mzd,mzd0,mpd,mpd0,rhod,rhod0,rho,rho0,s)
  end subroutine read_previous_results

  subroutine read_initial_values(filename,basis,nAtoms,initial_values)
  !> Reads the initial magnetization for nAtoms in 'filename'
    use mod_kind,       only: dp
    use mod_parameters, only: output
    use mod_constants,  only: deg2rad
    use mod_tools,      only: read_data
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: rField,abortProgram
    use mod_magnet,     only: lconstraining_field
    implicit none
    character(len=*), intent(in)  :: filename,basis
    integer,          intent(in)  :: nAtoms
    real(dp),         intent(out) :: initial_values(nAtoms,3)
    integer  :: i,err
    real(dp) :: temp(nAtoms,3)

    if(rField == 0) &
      write(output%unit_loop,"('[read_initial_values] Reading initial values from file ',a)") "'" // filename // "'"

    open(unit=321,file=trim(filename),status="old",iostat=err)
    ! If file cannot be read or if basis is not given
    if((rField == 0).and.(err/=0)) then
      if(lconstraining_field) write(output%unit_loop,"('[read_initial_values] File ',a,' is required for constraning field!')") "'" // filename // "'"
      call abortProgram("[read_initial_values] File '" // trim(filename) // "'' does not exist!")
    end if

    select case(basis)
    case("cartesian","c")
      call read_data(321,nAtoms,3,initial_values)
    case("spherical","s")
      call read_data(321,nAtoms,3,temp)
      do i=1,nAtoms
        initial_values(i,1) = temp(i,1)*sin(temp(i,2)*deg2rad)*cos(temp(i,3)*deg2rad)
        initial_values(i,2) = temp(i,1)*sin(temp(i,2)*deg2rad)*sin(temp(i,3)*deg2rad)
        initial_values(i,3) = temp(i,1)*cos(temp(i,2)*deg2rad)
      end do
    case("bravais","b")
      call read_data(321,nAtoms,3,temp)
      do i=1,nAtoms
        initial_values(i,:) = temp(i,1)*s%a1 + temp(i,2)*s%a2 + temp(i,3)*s%a3
      end do
    case("neighbor","n")
      call read_data(321,nAtoms,1,temp)
      call abortProgram("[read_initial_values] initial_values not implemented for basis = neighbor")
      ! do i=1,nAtoms
      !   initialmag(i,:) = temp(i,1)
      ! end do
    case default
      call abortProgram("[read_initial_values] basis wrongly defined: " // trim(basis) // ".")
    end select

    close(unit=321)

  end subroutine read_initial_values

  subroutine read_sc_results(err,lsuccess)
  !> This subroutine reads previous band-shifting and magnetization results
    use mod_kind,              only: dp
    use mod_constants,         only: cZero
    use mod_parameters,        only: output,dfttype,dimH
    use EnergyIntegration,     only: parts
    use mod_system,            only: s => sys
    use mod_tools,             only: replaceStr,vec_norm
    use mod_superconductivity, only: delta_sc
    use mod_magnet,            only: rho,mp,mx,my,mz,rhod,mpd,mxd,myd,mzd,hw_count,mabsd,m_fix,m_fix_abs,lconstraining_field,bc
    use mod_mpi_pars,          only: rField,FieldComm,ierr,MPI_DOUBLE_PRECISION
    implicit none
    character(len=300)  :: file="", bcfile=""
    integer,intent(out) :: err
    logical,intent(out) :: lsuccess
    integer             :: i,j,nOrbs,mu,mud,i0,i1
    real(dp)            :: previous_results(5*dimH/2) ! Divide by 2 because of the spins in dimH
    real(dp)            :: previous_Ef

    external :: MPI_Bcast

    if(trim(scfile) /= "") then
      open(unit=99,file=scfile,status="old",iostat=err)
      if(err/=0) then
        if(rField == 0) write(output%unit_loop,"('*** WARNING: Self-consistency file given on input file does not exist! Using default... ***')")
        scfile = " "
      end if
      file = trim(scfile)
      close(99)
    end if

    lsuccess = .false.
    !   Reading previous results (mx, my, mz and n) from files (if available)
    if(trim(scfile)=="") then ! If a filename is not given in inputcard (or does not exist), use the default one
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
        open(unit = 99,file=scfile, status = "old", iostat = err)
        if(err==0 .and. rField==0) then
          write(output%unit_loop,"('[read_sc_results] Using filename given in input file for self-consistency:')")
          write(output%unit_loop,"(a)") trim(scfile)
          file = trim(scfile)
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
            file = trim(scfile)
          end if
          lsuccess   = .true. ! something was read (different parameters)
        end if
      end if
    end if
    if(err==0) then
      if(rField==0) then
        i1 = 0
        do i=1,s%nAtoms
          i0 = i1 + 1
          i1 = i0 + 5*s%Types(s%Basis(i)%Material)%nOrb-1
          read(99,fmt=*) (previous_results(j), j=i0,i1)
        end do
        read(99,fmt=*) previous_Ef
      end if

      call MPI_Bcast(previous_results,size(previous_results),MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
      call MPI_Bcast(previous_Ef,1,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)

      rho = 0._dp
      mx  = 0._dp
      my  = 0._dp
      mz  = 0._dp
      mp  = cZero
      delta_sc = 0._dp
      rhod = 0._dp
      mxd  = 0._dp
      myd  = 0._dp
      mzd  = 0._dp
      
      i1 = 0
      do i=1,s%nAtoms
        i0 = i1 + 1
        i1 = i0 + 5*s%Types(s%Basis(i)%Material)%nOrb-1

        nOrbs = s%Types(s%Basis(i)%Material)%nOrb

        rho(1:nOrbs,i) = previous_results(        i0:  nOrbs+i0-1)
        mx (1:nOrbs,i) = previous_results(  nOrbs+i0:2*nOrbs+i0-1)
        my (1:nOrbs,i) = previous_results(2*nOrbs+i0:3*nOrbs+i0-1)
        mz (1:nOrbs,i) = previous_results(3*nOrbs+i0:4*nOrbs+i0-1)
        mp (1:nOrbs,i) = cmplx(mx(1:nOrbs,i),my(1:nOrbs,i),dp)
        delta_sc(1:nOrbs,i) = previous_results(4*nOrbs+i0:i1)

        do mud=1, s%Types(s%Basis(i)%Material)%ndOrb
          ! Corresponding orbital of this atom
          mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)

          rhod(i) = rhod(i) + rho(mu,i)
          mxd (i) = mxd (i) + mx (mu,i)
          myd (i) = myd (i) + my (mu,i)
          mzd (i) = mzd (i) + mz (mu,i)
        end do
        mabsd(i) = sqrt(mxd(i)**2+myd(i)**2+mzd(i)**2)
        mpd(i)   = cmplx(mxd(i),myd(i),dp)
      end do

      call calcMagAngle()

      s%Ef = previous_Ef

      ! Setting fixed magnetization from initialmag and reading constranining fields
      if(lconstraining_field) then
        ! Starting magnetization
        allocate(initialmag(s%nAtoms,3))
        call read_initial_values("initialmag",trim(magbasis),s%nAtoms,initialmag)

        do i=1,s%nAtoms
          m_fix_abs(i) = vec_norm(initialmag(i,:),3)
          m_fix(:,i) = initialmag(i,:)/m_fix_abs(i)
        end do

        if(rField==0) then
          bcfile=replaceStr( string=file, search="selfconsistency_", substitute="constraniningfield_" )
          open(unit=999,file=bcfile,status="old",iostat=err)
          if(err==0) then
            write(output%unit_loop,"('[read_sc_results] Reading constranining fields from file:')")
            write(output%unit_loop,"(a)") trim(bcfile)
            do i=1,s%nAtoms
              read(999,fmt=*) (bc(mu,i), mu=1,3)
            end do
          else
            write(output%unit_loop,"('[read_sc_results] Cannot read constranining fields from file:')")
            write(output%unit_loop,"(a)") trim(bcfile)
            write(output%unit_loop,"('[read_sc_results] Restarting from bc = 0.0')")
            bc = 0._dp
            lsuccess = .true. ! Marking that parameters are different
          end if
          close(999)
        end if
        call MPI_Bcast(bc,3*s%nAtoms,MPI_DOUBLE_PRECISION,0,FieldComm,ierr)
      end if

      if(lsuccess) then
        err = 1   ! Read different parameters
      else
        lsuccess   = .true. ! Read same parameters (err=0)
      end if
    end if
    close(99)
  end subroutine read_sc_results


  subroutine calcMagAngle()
  !> Calculate magnetization angles and set spherical variables
    use mod_constants,        only: rad2deg
    use mod_system,           only: s => sys
    use mod_magnet,           only: mx,my,mz,mabs,mxd,myd,mzd,mabsd, &
                                    mtheta,mphi,mvec_cartesian,mdvec_cartesian, &
                                    mvec_spherical,mtotal_cartesian,mtotal_spherical,lrot
    implicit none
    integer :: i

    ! Calculating new angles of GS magnetization in units of pi and magnetization vector
    do i = 1,s%nAtoms
      mabs(i)   = sqrt((sum(mx(:,i))**2)+(sum(my(:,i))**2)+(sum(mz(:,i))**2))
      if(mabs(i)>1.e-8_dp) then
        mtheta(i) = merge(acos(sum(mz(:,i))/mabs(i))*rad2deg,0._dp,mabs(i)>1.e-8_dp)
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
      mdvec_cartesian(1,i) = mxd(i)
      mdvec_cartesian(2,i) = myd(i)
      mdvec_cartesian(3,i) = mzd(i)
      mabsd(i) = sqrt((mxd(i)**2)+(myd(i)**2)+(mzd(i)**2))
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


  subroutine calcSelfConsistency()
  !> This subroutine performs the self-consistency
    use mod_kind,              only: dp
    use mod_parameters,        only: output,lfixEf
    use mod_magnet,            only: bc,rho,rhod,rhot,mxd,myd,mpd,mzd,mabsd,mp,mx,my,mz,constr_type,sb_matrix
    use mod_mpi_pars,          only: rField
    use mod_system,            only: s => sys
    use mod_dnsqe,             only: dnsqe
    use mod_superconductivity, only: lsuperCond,delta_sc
    use mod_expectation,       only: expectation_values
#ifdef _GPU
    use nvtx,                  only: nvtxStartRange,nvtxEndRange
#endif
    implicit none
    real(dp), allocatable :: fvec(:),jac(:,:),wa(:),sc_solu(:)
    real(dp), allocatable :: diag(:),qtf(:)
    real(dp), dimension(s%nAtoms) :: in_varx,in_vary,in_varz
    integer               :: i,mu,mud,lwa,ifail=0

    ! Setting up number of equations (total and accumulated per atom)
    allocate( neq_per_atom(s%nAtoms+1) )
    neq_per_atom(1) = 0
    do i = 1,s%nAtoms
      neq_per_atom(i+1) = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + merge(merge(1,3,constr_type==1),0,abs(s%Basis(i)%Um)>1.e-8_dp) + merge(s%Types(s%Basis(i)%Material)%nOrb,0,lsupercond)
    end do
    neq = neq_per_atom(s%nAtoms+1)
    if(.not.lfixEf) neq = neq+1

    if(neq.ne.0) then
      ! Allocating variables dependent on the number of equations
      allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )
      select case(constr_type)
      case(1)
        in_varx(:) = mabsd(:)
      case(2)
        in_varx(:) = bc(1,:)
        in_vary(:) = bc(2,:)
        in_varz(:) = bc(3,:)
        ! Updating field with new constrainig values
        call sb_matrix(s)
      case default
        in_varx(:) = mxd(:)
        in_vary(:) = myd(:)
        in_varz(:) = mzd(:)
      end select

      !!! To be deleted
      ! write(*,*) "Rho = ", rho
      ! write(*,*) "rhod", rhod
      ! write(*,*) "mpd", mpd
      ! write(*,*) "mzd", mzd
      ! write(*,*) "in_varx", in_varx
      ! write(*,*) "in_vary", in_vary
      ! write(*,*) "in_varz", in_varz
      !!!

      ! Putting read n and m existing solutions into sc_solu (first guess of the subroutine)
      call set_hamiltonian_variables(.false.,s,neq,sc_solu,rho,rhod,in_varx,in_vary,in_varz,delta_sc)

      call print_sc_step(rhod,mpd,mzd,delta_sc,s)

      if(rField == 0) &
        write(output%unit_loop,"('[self_consistency] Starting self-consistency:')")

#ifdef _GPU
      ! Starting marker of dnsqe for profiler
      call nvtxStartRange("dnsqe",1)
#endif

      ! Performing selfconsistency finding root of non-linear system of equations (SLATEC)
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa) )
      if(lnojac) then
        call dnsqe(sc_eqs,sc_jac,2,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs,sc_jac,1,neq,sc_solu,fvec,mag_tol,0,ifail,wa,lwa)
      end if
      ifail = ifail-1

#ifdef _GPU
      ! End of dnsqe marker
      call nvtxEndRange
#endif

      ! Deallocating variables dependent on the number of equations
      deallocate(neq_per_atom,sc_solu,diag,qtf,fvec,jac,wa)
    else
      if(rField == 0) &
        write(output%unit_loop,"('[self_consistency] Self-consistency not required (neq=0). Calculating expectations...')")
      call expectation_values(s,rho,mp,mx,my,mz,delta_sc)
      rhod = 0._dp
      rhot = 0._dp
      mpd  = 0._dp
      mxd  = 0._dp
      myd  = 0._dp
      mzd  = 0._dp
      do i = 1,s%nAtoms
        do mud = 1,s%Types(s%Basis(i)%Material)%ndOrb
          mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
          rhod(i) = rhod(i) + rho(mu,i)
          mxd (i) = mxd (i) + mx (mu,i)
          myd (i) = myd (i) + my (mu,i)
          mzd (i) = mzd (i) + mz (mu,i)
        end do
        rhot = rhot + sum(rho(:,i))
        mpd(i) = cmplx(mxd(i),myd(i),dp)
      end do

    end if

    ! Calculating the magnetization in cartesian and spherical coordinates
    call calcMagAngle()

  end subroutine calcSelfConsistency


  subroutine calcSelfConsistency_trans_constr()
  !> This subroutine performs the self-consistency to constrain transversely the magnetic moments
  !> It is done by fixing the transverse contraining fields and converging the length of the magnetization
  !> for that particular field. 
    use mod_kind,       only: dp
    use mod_parameters, only: output
    use mod_magnet,     only: bc,sb_matrix,m_fix,mx,my,mz,iter,cmix
    use mod_system,     only: s => sys
    use mod_tools,      only: transvComponent
    use mod_mpi_pars,   only: rField
    implicit none
    real(dp), dimension(3,s%nAtoms) :: bc_old,m_out_unit,delta_m
    real(dp), dimension(s%nAtoms)   :: m_out_abs
    real(dp) :: rconv 
    integer  :: i,iter_b

    if(rField == 0) &
      write(output%unit_loop,"('[calcSelfConsistency_trans_constr] Starting constraining the transverse components of the magnetic moments...')")

    ! Updating field with new constrainig values
    call sb_matrix(s)

    rconv = 999.e0_dp
    iter_b = 1
    do while(rconv>mag_tol)

      ! Self-consistency for a given transverse constraning field
      iter = 0
      call calcSelfConsistency()

      do i = 1,s%nAtoms
        m_out_abs(i) = sqrt(sum(mx(:,i))**2+sum(my(:,i))**2+sum(mz(:,i))**2)
        m_out_unit(1,i) = sum(mx(:,i))/m_out_abs(i) 
        m_out_unit(2,i) = sum(my(:,i))/m_out_abs(i)
        m_out_unit(3,i) = sum(mz(:,i))/m_out_abs(i)
      end do
      delta_m = m_fix-m_out_unit

      ! Calculating convergency parameter
      rconv = 0.0_dp
      do i = 1,s%nAtoms
        call transvComponent(delta_m(:,i),m_fix(:,i))
        rconv = rconv + dot_product(delta_m(:,i),delta_m(:,i))
      end do

      ! Storing the old value of constraning and setting up new value
      bc_old(:,:) = bc(:,:)
      bc(:,:)     = bc_old(:,:) + cmix*delta_m(:,:)
      if(rField == 0) then
        write(output%unit_loop,"('|################ Iteration ',i4,' of constraining field cycle #################|')") iter_b
        do i = 1,s%nAtoms
          write(output%unit_loop,"(' bc         (',i4,') = ',3es16.8)") i, bc(:,i)
          write(output%unit_loop,"(' delta_m    (',i4,') = ',3es16.8)") i, delta_m(:,i)
          write(output%unit_loop,"(' m_out_unit (',i4,') = ',3es16.8,' |m_out| = ',es16.8)") i, m_out_unit(:,i),m_out_abs(i)
          write(output%unit_loop,"(' m_fix      (',i4,') = ',3es16.8)") i, m_fix(:,i)
        end do
        write(output%unit_loop,"(' Convergence parameter: ',es16.8)") rconv
      end if
      ! Updating the constraninig field in the hamiltonian
      call sb_matrix(s)

      iter_b = iter_b+1
    end do

  end subroutine calcSelfConsistency_trans_constr
  

  subroutine calcSelfConsistency_linear()
  !> Routine not implemented yet! To be done for Filipe :)
  !> This subroutine performs the self-consistency using linear mixing
    use mod_parameters,        only: output
    implicit none

    write(output%unit_loop,"('[calcSelfConsistency_linear] Linear mixing not implemented yet. To be done.')")

  end subroutine calcSelfConsistency_linear


  subroutine check_jacobian(n,x)
  !> This subroutine can be used to check the jacobian implementation
  !> It compares the jacobian subroutine with the one obtained numerically 
  !> by calculating the function is different points
    use mod_parameters, only: output,lcheckjac
    use mod_mpi_pars,   only: rField,abortProgram
    use mod_chkder,     only: chkder
    implicit none
    integer :: liw,lw,ifail
    integer , allocatable :: iw(:)
    real(dp), allocatable :: w(:)
    integer , intent(in)  :: n
    real(dp), intent(in)  :: x(n)
    real(dp) :: jac(n,n)
    real(dp) :: fvec(n)
    ! integer :: i
    real(dp) :: fvecp(n),xp(n),err(n)

    liw = 1
    lw  = (4+n)*n
    allocate(iw(liw),w(lw))

    lcheckjac = .false.

    if(rField == 0) &
    write(output%unit_loop,"('[check_jacobian] Checking Jacobian if Jacobian is correct...')", advance='no')

    call chkder(n,n,x,fvec,jac,n,xp,fvecp,1,err)
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
    use mod_magnet,            only: rho,rhot,rhod,mp,mx,my,mz,mpd,mpd0,mxd,myd,mzd,mzd0,rhod0,rho0
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
    call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mxd_in,myd_in,mzd_in,delta_sc_in)
    do i = 1,s%nAtoms
      mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
    end do

    ! Update Hubbard term in Hamiltonian
    call update_Umatrix(mzd_in,mzd0,mpd_in,mpd0,rhod_in,rhod0,rho_in,rho0,s)
    ! Update electron-hole coupling in Hamiltonian
    if(lsuperCond) call update_delta_sc(s,delta_sc_in)

    call expectation_values(s,rho,mp,mx,my,mz,deltas)

    rhod = 0._dp
    rhot = 0._dp
    mpd  = 0._dp
    mxd  = 0._dp
    myd  = 0._dp
    mzd  = 0._dp
    do i = 1,s%nAtoms
      do mud = 1,s%Types(s%Basis(i)%Material)%ndOrb
        mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
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
    !> Calculated the Jacobian of the spin magnetization
    use mod_kind,          only: dp,int64
    use mod_constants,     only: pi,cZero,cOne
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
    complex(dp), dimension(s%nOrb2,s%nOrb2)      :: gij,gji,temp,paulitemp
    complex(dp), dimension(s%nOrb2,s%nOrb2,4)    :: temp1,temp2
    complex(dp), dimension(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms) :: gf,gvg
    integer(int64) :: ix
    integer     :: i,j,kounti,kountj,mu,nu,mud,nud,sigma,sigmap,nOrb2_i,nOrb2_j
    real(dp)    :: kp(3), ep
    complex(dp) :: weight
    complex(dp) :: halfUn(s%nAtoms),halfUm(s%nAtoms)
    integer :: ncount2

    external :: zgemm,MPI_Allreduce

    if(rField == 0) &
      write(output%unit_loop,"('[calcJacobian_greenfunction] Calculating the Jacobian...')")

    ncount2=N*N

    ! Prefactor -U/2 in dH/dn and dH/dm
    do i=1,s%nAtoms
      halfUn(i) = -0.5_dp*s%Basis(i)%Un
      halfUm(i) = -0.5_dp*s%Basis(i)%Um
    end do

    ! Build local hamiltonian
    if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

    jacobian = 0._dp

    !$omp parallel default(none) &
    !$omp& private(ix,i,j,kounti,kountj,nOrb2_i,nOrb2_j,mu,mud,nud,nu,sigma,sigmap,ep,kp,weight,gf,gvg,gij,gji,temp,temp1,temp2,paulitemp) &
    !$omp& shared(llineargfsoc,llinearsoc,lfixEf,neq,local_points,s,neq_per_atom,realBZ,bzs,E_k_imag_mesh,y,eta,wght,halfUn,halfUm,jacobian)
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
          nOrb2_i = s%Types(s%Basis(i)%Material)%nOrb2
          nOrb2_j = s%Types(s%Basis(j)%Material)%nOrb2

          ! (Re)starting the column counter in each line
          kountj = neq_per_atom(j)

          ! First product: temp1 =  pauli*g_ij
          gij(1:nOrb2_i,1:nOrb2_j) = gf(1:nOrb2_i,1:nOrb2_j,i,j)
          do sigma = 1,4
            paulitemp(1:nOrb2_i,1:nOrb2_i) = merge(s%Types(s%Basis(i)%Material)%pauli_orb(sigma-1,1:nOrb2_i,1:nOrb2_i),s%Types(s%Basis(i)%Material)%pauli_dorb(sigma-1,1:nOrb2_i,1:nOrb2_i),sigma==1) ! sigma_0 includes sp orbitals for the last line - total charge neutrality
            call zgemm('n','n',nOrb2_i,nOrb2_j,nOrb2_i,cOne,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)
            temp1(1:nOrb2_i,1:nOrb2_j,sigma) = temp(1:nOrb2_i,1:nOrb2_j)
          end do

          if(abs(s%Basis(j)%Un)>1.e-8_dp) then
            do nu=1,s%Types(s%Basis(j)%Material)%ndOrb
              ! Restarting line counter
              kounti = neq_per_atom(i)

              nud = s%Types(s%Basis(j)%Material)%dOrbs(nu)
              mud = nud+s%Types(s%Basis(j)%Material)%nOrb
              ! Second product: temp2 = (U_j\delta_{mu,nu} -U_j/2) * pauli * g_ji
              gji(1:nOrb2_j,1:nOrb2_i) = gf(1:nOrb2_j,1:nOrb2_i,j,i)
              paulitemp(1:nOrb2_j,1:nOrb2_j) = s%Types(s%Basis(j)%Material)%pauli_dorb(0,1:nOrb2_j,1:nOrb2_j) ! identity excluding d orbitals from charge part
              paulitemp(nud,nud) = paulitemp(nud,nud) - 2._dp
              paulitemp(mud,mud) = paulitemp(mud,mud) - 2._dp
              call zgemm('n','n',nOrb2_j,nOrb2_i,nOrb2_j,halfUn(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
              temp2(1:nOrb2_j,1:nOrb2_i,1) = temp(1:nOrb2_j,1:nOrb2_i)

              ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*pauli* g_ji

              ! Charge density-charge density part
              gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j,1)
              gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i,1)
              call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

              if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
                  mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                  nud = mud+s%Types(s%Basis(i)%Material)%nOrb

                  jacobian(kounti+mu,kountj+nu) = jacobian(kounti+mu,kountj+nu) + real(temp(mud,mud) + temp(nud,nud))
                end do
                kounti = neq_per_atom(i) + s%ndOrb
              end if

              ! Last line (Total charge neutrality)
              if(.not.lfixEf) then
                do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
                  nud = mu+s%Types(s%Basis(i)%Material)%nOrb
                  jacobian(neq,kountj+nu) = jacobian(neq,kountj+nu) + real(temp(mu,mu) + temp(nud,nud))
                end do
              end if

              ! Magnetic density-charge density part
              if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                do sigma = 2,4
                  gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j,sigma)
                  gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i,1)
                  call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                  do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
                    mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                    nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                    jacobian(kounti+sigma-1,kountj+nu) = jacobian(kounti+sigma-1,kountj+nu) + real(temp(mud,mud) + temp(nud,nud))
                  end do
                end do
                kounti = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
              end if

            end do

            kountj = neq_per_atom(j) + s%Types(s%Basis(j)%Material)%ndOrb
          end if ! |Un(j)| > 0

          if(abs(s%Basis(j)%Um)>1.e-8_dp) then

            gji(1:nOrb2_j,1:nOrb2_i) = gf(1:nOrb2_j,1:nOrb2_i,j,i)

            do sigmap = 2,4
              ! Second product: temp2 = (-U/2) * pauli * g_ji
              paulitemp(1:nOrb2_j,1:nOrb2_j) = s%Types(s%Basis(j)%Material)%pauli_dorb(sigmap-1,1:nOrb2_j,1:nOrb2_j)
              call zgemm('n','n',nOrb2_j,nOrb2_i,nOrb2_j,halfUm(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
              temp2(1:nOrb2_j,1:nOrb2_i,sigmap) = temp(1:nOrb2_j,1:nOrb2_i)
            end do

            ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

            ! Charge density-magnetic density part
            do sigmap = 2,4
              ! Restarting line counter
              kounti = neq_per_atom(i)

              gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j,1)
              gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i,sigmap)
              call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

              if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                do mu = 1,s%Types(s%Basis(i)%Material)%ndOrb
                  mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                  nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                  jacobian(kounti+mu,kountj+sigmap-1) = jacobian(kounti+mu,kountj+sigmap-1) + real(temp(mud,mud) + temp(nud,nud))
                end do
                kounti = neq_per_atom(i) + s%Types(s%Basis(i)%Material)%ndOrb
              end if

              ! Last line (Total charge neutrality)
              if(.not.lfixEf) then
                do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
                  nud = mu+s%Types(s%Basis(i)%Material)%nOrb
                  jacobian(neq,kountj+sigmap-1) = jacobian(neq,kountj+sigmap-1) + real(temp(mu,mu) + temp(nud,nud))
                end do
              end if

              ! Magnetic density-magnetic density part
              if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                do sigma = 2,4
                  gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j,sigma)
                  gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i,sigmap)
                  call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                  do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
                    mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                    nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                    jacobian(kounti+sigma-1,kountj+sigmap-1) = jacobian(kounti+sigma-1,kountj+sigmap-1) + real(temp(mud,mud) + temp(nud,nud))
                  end do
                end do
                kounti = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
              end if
            end do
            kountj = neq_per_atom(j) + merge(s%Types(s%Basis(j)%Material)%ndOrb,0,abs(s%Basis(j)%Un)>1.e-8_dp) + 3
          end if ! |Um(j)| > 0

          ! **** Removing non-linear (quadratic) terms: ****
          if(llineargfsoc .or. llinearsoc) then
            ! Restarting column counter
            kountj = neq_per_atom(j)

            gij(1:nOrb2_i,1:nOrb2_j) = gvg(1:nOrb2_i,1:nOrb2_j,i,j)

            do sigma = 1,4
              ! First product: temp1 =  pauli*g_ij
              paulitemp(1:nOrb2_i,1:nOrb2_i) = merge(s%Types(s%Basis(i)%Material)%pauli_orb(sigma-1,1:nOrb2_i,1:nOrb2_i),s%Types(s%Basis(i)%Material)%pauli_dorb(sigma-1,1:nOrb2_i,1:nOrb2_i),sigma==1) ! sigma_0 includes sp orbitals for the last line - total charge neutrality
              call zgemm('n','n',nOrb2_i,nOrb2_j,nOrb2_i,cOne,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)
              temp1(1:nOrb2_i,1:nOrb2_j,sigma) = temp(1:nOrb2_i,1:nOrb2_j)
            end do

            if(abs(s%Basis(j)%Un)>1.e-8_dp) then
              do nu=1,s%Types(s%Basis(j)%Material)%ndOrb
                ! Restarting line counter
                kounti = neq_per_atom(i)

                nud = s%Types(s%Basis(j)%Material)%dOrbs(nu)
                mud = nud+s%Types(s%Basis(j)%Material)%nOrb

                gji(1:nOrb2_j,1:nOrb2_i) = gvg(1:nOrb2_j,1:nOrb2_i,j,i)
                ! Second product: temp2 = (-U/2) * pauli * g_ji
                paulitemp(1:nOrb2_j,1:nOrb2_j) = s%Types(s%Basis(j)%Material)%pauli_dorb(0,1:nOrb2_j,1:nOrb2_j)
                paulitemp(nud,nud) = paulitemp(nud,nud) - 2._dp
                paulitemp(mud,mud) = paulitemp(mud,mud) - 2._dp
                call zgemm('n','n',nOrb2_j,nOrb2_i,nOrb2_j,halfUn(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
                temp2(1:nOrb2_j,1:nOrb2_i, 1) = temp(1:nOrb2_j,1:nOrb2_i)

                ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji

                ! Charge density-charge density part
                gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j, 1)
                gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i, 1)
                call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                  do mu=1,s%Types(s%Basis(j)%Material)%ndOrb
                    mud = s%Types(s%Basis(j)%Material)%dOrbs(mu)
                    nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                    jacobian(kounti+mu,kountj+nu) = jacobian(kounti+mu,kountj+nu) - real(temp(mud,mud) + temp(nud,nud))
                  end do
                  kounti = neq_per_atom(i) + s%Types(s%Basis(i)%Material)%ndOrb
                end if

                ! Last line (Total charge neutrality)
                if(.not.lfixEf) then
                  do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
                    nud = mu+s%Types(s%Basis(i)%Material)%nOrb
                    jacobian(neq,kountj+nu) = jacobian(neq,kountj+nu) - real(temp(mu,mu) + temp(nud,nud))
                  end do
                end if

                ! Magnetic density-charge density part
                if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                  do sigma = 2,4
                    gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j, sigma)
                    gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i, 1)
                    call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                    do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
                      mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                      nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                      jacobian(kounti+sigma-1,kountj+nu) = jacobian(kounti+sigma-1,kountj+nu) - real(temp(mud,mud) + temp(nud,nud))
                    end do
                  end do
                  kounti = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
                end if

              end do
              kountj = neq_per_atom(j) + s%Types(s%Basis(j)%Material)%ndOrb
            end if ! |Un(j)| > 0

            if(abs(s%Basis(j)%Um)>1.e-8_dp) then

              gji(1:nOrb2_j,1:nOrb2_i) = gvg(1:nOrb2_j,1:nOrb2_i,j,i)

              do sigmap = 2,4
                ! Second product: temp2 = (-U/2) * pauli * g_ji
                paulitemp(1:nOrb2_j,1:nOrb2_j) = s%Types(s%Basis(j)%Material)%pauli_dorb(sigmap-1,1:nOrb2_j,1:nOrb2_j) 
                call zgemm('n','n',nOrb2_j,nOrb2_i,nOrb2_j,halfUm(j),paulitemp,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)
                temp2(1:nOrb2_j,1:nOrb2_i,sigmap) = temp(1:nOrb2_j,1:nOrb2_i)
              end do

              ! Full product:  sigma.g.dHdx.g = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*pauli* g_ji

              ! Charge density-magnetic density part
              do sigmap = 2,4
                ! Restarting line counter
                kounti = neq_per_atom(i)

                gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j,1)
                gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i,sigmap)
                call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                if(abs(s%Basis(i)%Un)>1.e-8_dp) then
                  do mu = 1,s%Types(s%Basis(i)%Material)%ndOrb
                    mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                    nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                    jacobian(kounti+mu,kountj+sigmap-1) = jacobian(kounti+mu,kountj+sigmap-1) - real(temp(mud,mud) + temp(nud,nud))
                  end do
                  kounti = neq_per_atom(i) + s%Types(s%Basis(i)%Material)%ndOrb
                end if

                ! Last line (Total charge neutrality)
                if(.not.lfixEf) then
                  do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
                    nud = mu+s%Types(s%Basis(i)%Material)%nOrb
                    jacobian(neq,kountj+sigmap-1) = jacobian(neq,kountj+sigmap-1) - real(temp(mu,mu) + temp(nud,nud))
                  end do
                end if

                ! Magnetic density-magnetic density part
                if(abs(s%Basis(i)%Um)>1.e-8_dp) then
                  do sigma = 2,4
                    gij(1:nOrb2_i,1:nOrb2_j) = temp1(1:nOrb2_i,1:nOrb2_j,sigma)
                    gji(1:nOrb2_j,1:nOrb2_i) = temp2(1:nOrb2_j,1:nOrb2_i,sigmap)
                    call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,weight,gij,s%nOrb2,gji,s%nOrb2,cZero,temp,s%nOrb2)

                    do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
                      mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                      nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                      jacobian(kounti+sigma-1,kountj+sigmap-1) = jacobian(kounti+sigma-1,kountj+sigmap-1) - real(temp(mud,mud) + temp(nud,nud))
                    end do
                  end do
                  kounti = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
                end if
              end do
              kountj = neq_per_atom(j) + merge(s%Types(s%Basis(j)%Material)%ndOrb,0,abs(s%Basis(j)%Un)>1.e-8_dp) + 3
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
          nOrb2_i = s%Types(s%Basis(i)%Material)%nOrb2
          
          ! Restarting line counter
          kounti = neq_per_atom(i)

          gij(1:nOrb2_i,1:nOrb2_i) = gf(1:nOrb2_i,1:nOrb2_i,i,i)

          ! Charge density per orbital lines
          ! temp1 =  pauli*g_ii
          paulitemp(1:nOrb2_i,1:nOrb2_i) = s%Types(s%Basis(i)%Material)%pauli_orb(0,1:nOrb2_i,1:nOrb2_i) ! sigma_0 includes sp orbitals for the last line - total charge neutrality
          call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_i,weight,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)

          if(abs(s%Basis(i)%Un)>1.e-8_dp) then
            do mu = 1,s%Types(s%Basis(i)%Material)%ndOrb
              mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
              nud = mud+s%Types(s%Basis(i)%Material)%nOrb
              jacobian(kounti+mu,neq) = jacobian(kounti+mu,neq) - aimag(temp(mud,mud) + temp(nud,nud))
            end do
            kounti = neq_per_atom(i) + s%Types(s%Basis(i)%Material)%ndOrb
          end if

          ! Last line (Total charge neutrality)
          if(.not.lfixEf) then
            do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
              nud = mu+s%Types(s%Basis(i)%Material)%nOrb
              jacobian(neq,neq) = jacobian(neq,neq)  - aimag(temp(mu,mu) + temp(nud,nud))
            end do
          end if

          ! Magnetic density-last column part
          if(abs(s%Basis(i)%Um)>1.e-8_dp) then
            do sigma = 2,4
              ! temp1 =  pauli*g_ii
              paulitemp(1:nOrb2_i,1:nOrb2_i) = s%Types(s%Basis(i)%Material)%pauli_dorb(sigma-1,1:nOrb2_i,1:nOrb2_i)
              call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_i,weight,paulitemp,s%nOrb2,gij,s%nOrb2,cZero,temp,s%nOrb2)

              do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
                mud = s%Types(s%Basis(i)%Material)%dOrbs(mu)
                nud = mud+s%Types(s%Basis(i)%Material)%nOrb
                jacobian(kounti+sigma-1,neq) = jacobian(kounti+sigma-1,neq) - dimag(temp(mud,mud) + temp(nud,nud))
              end do
            end do
            kounti = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + 3
          end if

          ! No linear correction is needed since it is a single Green function

        end do ! End nAtoms i loop
      end do ! End nkpt loop
      !$omp end do nowait
    end if
    !$omp end parallel

    call MPI_Allreduce(MPI_IN_PLACE, jacobian, ncount2, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)

    jacobian = jacobian/pi
    do i = 1, neq-1
      jacobian(i,i) = jacobian(i,i) - 1._dp
    end do
  end subroutine calcJacobian_greenfunction


  subroutine rotate_magnetization_to_field()
  !> Rotate the magnetization to the direction of the field (useful for SOC=F)
    use mod_kind,       only: dp
    use mod_constants,  only: deg2rad
    use mod_magnet,     only: lfield,hw_count,hw_list,hhw,mp,mx,my,mz,mpd,mxd,myd,mzd
    use mod_parameters, only: output
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: rField
    implicit none
    integer  :: i,mu,mud,signal
    real(dp) :: mdotb,mabs(s%nOrb,s%nAtoms)

    if(rField == 0) &
      write(output%unit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    if(.not.lfield) then
      if(rField == 0) &
        write(output%unit_loop,"('[rotate_magnetization_to_field] Field if OFF! No rotation is done.')")
      return
    end if

    mxd = 0._dp
    myd = 0._dp
    mzd = 0._dp
    do i = 1, s%nAtoms
      do mu=1,s%Types(s%Basis(i)%Material)%nOrb
        mdotb   = hhw(1,i)*mx(mu,i)+hhw(2,i)*my(mu,i)+hhw(3,i)*mz(mu,i)
        signal  = int(mdotb/abs(mdotb))
        mabs(mu,i) = sqrt((mx(mu,i)**2)+(my(mu,i)**2)+(mz(mu,i)**2))
        mx  (mu,i) = signal*mabs(mu,i)*sin(hw_list(hw_count,2)*deg2rad)*cos(hw_list(hw_count,3)*deg2rad)
        my  (mu,i) = signal*mabs(mu,i)*sin(hw_list(hw_count,2)*deg2rad)*sin(hw_list(hw_count,3)*deg2rad)
        mz  (mu,i) = signal*mabs(mu,i)*cos(hw_list(hw_count,2)*deg2rad)
        mp  (mu,i) = cmplx(mx(mu,i),my(mu,i),dp)
      end do
      do mud=1, s%Types(s%Basis(i)%Material)%ndOrb
        ! Corresponding orbital of this atom
        mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
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

  subroutine print_sc_results()
  !> Writes the self-consistency results on the screen
    use mod_constants,         only: rad2deg
    use mod_parameters,        only: output,leigenstates
    use mod_system,            only: s => sys
    use mod_hamiltonian,       only: energy,energy_dc,energy_dc_n
    use mod_magnet,            only: rho,rhos,rhop,rhod,mvec_cartesian,mp,mvec_spherical, &
                                     lxm,lym,lzm,ltheta,lphi,labs,iter,lconstraining_field,bc
    use mod_superconductivity, only: lsuperCond,delta_sc
    use mod_tools,             only: vec_norm
    implicit none
    real(dp), dimension(s%nAtoms)   :: deltas,deltap,deltad
    real(dp), dimension(3,s%nAtoms) :: bc_spherical
    integer  :: i,mu

    write(output%unit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    if(leigenstates) then
      write(output%unit_loop,"(3x,'Ef=',f13.8,4x,'Eband=',f18.10,4x,'Etotal=',f18.10)") s%Ef,energy,energy+energy_dc+energy_dc_n
      write(output%unit_loop,"(3x,'Double counting:',4x,'Edc_m=',f18.10,5x,'Edc_n=',f18.10)") energy_dc,energy_dc_n
    else
      write(output%unit_loop,"(28x,'Ef=',f10.7)") s%Ef
    end if
    write(output%unit_loop,"(11x,' *************** Charge density: ****************')")
    rhos = 0._dp
    rhop = 0._dp
    do i=1,s%nAtoms
      ! Total s-orbital occupation
      do mu=1,s%Types(s%Basis(i)%Material)%nsOrb
        rhos(i) = rhos(i) + rho(s%Types(s%Basis(i)%Material)%sOrbs(mu),i)
      end do
      ! Total p-orbital occupation
      do mu=1,s%Types(s%Basis(i)%Material)%npOrb
        rhop(i) = rhop(i) + rho(s%Types(s%Basis(i)%Material)%pOrbs(mu),i)
      end do
      write(output%unit_loop,fmt="(a2,':',2x,'Ns=',f10.7,4x,'Np=',f10.7,4x,'Nd=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),rhos(i),rhop(i),rhod(i)
    end do

    if(lsupercond) then
      write(output%unit_loop,"(11x,' ******** Superconducting gap parameter: ********')")
      write(output%unit_loop,"(11x,' *** (Averages of the norms per orbital type) ***')")
      deltas = 0._dp
      deltap = 0._dp
      deltad = 0._dp
      do i=1,s%nAtoms
        ! Average s-orbital superconducting delta
        do mu=1,s%Types(s%Basis(i)%Material)%nsOrb
          deltas(i) = deltas(i) + delta_sc(s%Types(s%Basis(i)%Material)%sOrbs(mu),i)
        end do
        deltas(i) = deltas(i)/s%Types(s%Basis(i)%Material)%nsOrb
        ! Average p-orbital superconducting delta
        do mu=1,s%Types(s%Basis(i)%Material)%npOrb
          deltap(i) = deltap(i) + delta_sc(s%Types(s%Basis(i)%Material)%pOrbs(mu),i)
        end do
        deltap(i) = deltap(i)/s%Types(s%Basis(i)%Material)%npOrb
        ! Average d-orbital superconducting delta
        do mu=1,s%Types(s%Basis(i)%Material)%ndOrb
          deltad(i) = deltad(i) + delta_sc(s%Types(s%Basis(i)%Material)%dOrbs(mu),i)
        end do
        deltad(i) = deltad(i)/s%Types(s%Basis(i)%Material)%ndOrb
        write(output%unit_loop,fmt="(a2,':',2x,'Ds=',f10.7,4x,'Dp=',f10.7,4x,'Dd=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),deltas(i),deltap(i),deltad(i)
      end do
    end if

    write(output%unit_loop,"(11x,' *********** Magnetization components: ***********')")
    if(abs(sum(mp(:,:)))>1.e-7_dp) then
      do i=1,s%nAtoms
        write(output%unit_loop,"(a2,':',2x,'Mx=',f10.7,4x,'My=',f10.7,4x,'Mz=',f10.7,4x,'theta = ',f10.5,4x,'phi = ',f10.5)") trim(s%Types(s%Basis(i)%Material)%Name),mvec_cartesian(1,i),mvec_cartesian(2,i),mvec_cartesian(3,i),mvec_spherical(2,i),mvec_spherical(3,i)
      end do
    else
      do i=1,s%nAtoms
        write(output%unit_loop,"(a2,':',2x,'Mx=',f10.7,4x,'My=',f10.7,4x,'Mz=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),mvec_cartesian(1,i),mvec_cartesian(2,i),mvec_cartesian(3,i)
      end do
    end if
    if(lconstraining_field) then
      write(output%unit_loop,"(11x,' ************** Constraning fields: **************')")
      ! Calculating spherical components of constraning fields
      do i = 1,s%nAtoms
        bc_spherical(1,i) = vec_norm(bc(:,i),3)
        bc_spherical(2,i) = merge(acos(bc(3,i)/bc_spherical(1,i))*rad2deg,0._dp,bc_spherical(1,i)>1.e-8_dp)
        if(abs(bc_spherical(2,i))>1.e-8_dp) then
          if(abs(abs(bc_spherical(2,i))-180._dp)>1.e-8_dp) then
            bc_spherical(3,i) = atan2(bc(2,i),bc(1,i))*rad2deg
          else
            bc_spherical(3,i) = 0._dp
          end if
        else
          bc_spherical(3,i) = 0._dp
        end if
      end do
      if(abs(vec_norm(bc_spherical(2,:),s%nAtoms))>1.e-6_dp) then
        do i=1,s%nAtoms
          write(output%unit_loop,"(a2,':',2x,'Bx=',f10.7,4x,'By=',f10.7,4x,'Bz=',f10.7,4x,'theta = ',f10.5,4x,'phi = ',f10.5)") trim(s%Types(s%Basis(i)%Material)%Name),bc(1,i),bc(2,i),bc(3,i),bc_spherical(2,i),bc_spherical(3,i)
        end do
      else
        do i=1,s%nAtoms
          write(output%unit_loop,"(a2,':',2x,'Bx=',f10.7,4x,'By=',f10.7,4x,'Bz=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),bc(1,i),bc(2,i),bc(3,i)
        end do
      end if
    end if
    write(output%unit_loop,"(11x,' ****** Orbital components in global frame: ******')")
    if(sum(lxm(:)**2+lym(:)**2)>1.e-7_dp) then
      do i=1,s%nAtoms
        write(output%unit_loop,"(a2,':',2x,'Lx=',f10.7,4x,'Ly=',f10.7,4x,'Lz=',f10.7,4x,'theta = ',f10.5,4x,'phi = ',f10.5)") trim(s%Types(s%Basis(i)%Material)%Name),lxm(i),lym(i),lzm(i),ltheta(i),lphi(i)
      end do
    else
      do i=1,s%nAtoms
        write(output%unit_loop,"(a2,':',2x,'Lx=',f10.7,4x,'Ly=',f10.7,4x,'Lz=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),lxm(i),lym(i),lzm(i)
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
      write(output%unit_loop,"(a2,':',2x,' N=',f10.7,2x,' |M|=',f10.7,2x,' |L|=',f10.7)") trim(s%Types(s%Basis(i)%Material)%Name),sum(rho(:,i)),mvec_spherical(1,i),labs(i)
    end do
    write(output%unit_loop,"('|----------===================== (',i4.0,' iterations ) =====================----------|')") iter
  end subroutine print_sc_results


  subroutine write_sc_results(rho,mx,my,mz,bc)
  !> Writes the self-consistency results into files and broadcasts the scfile for the next iteration.
    use mod_parameters,        only: output,dfttype
    use EnergyIntegration,     only: parts
    use mod_superconductivity, only: delta_sc
    use mod_system,            only: s => sys
    use mod_tools,             only: replaceStr
    use mod_mpi_pars,          only: rField,MPI_CHARACTER,FieldComm,ierr
    implicit none
    real(dp),dimension(s%nOrb,s%nAtoms), intent(in) :: rho,mx,my,mz
    !> orbital- and site-dependent charge and magnetic moments
    real(dp),dimension(3,s%nAtoms),      intent(in), optional :: bc
    !> constraining magnetic fields
    character(len=30)  :: formatvar
    character(len=200) :: bcfile
    integer            :: i,mu

    external :: MPI_Bcast
    if(rField == 0) then
      ! Writing new results (mx, my, mz and n) to file
      write(output%unit_loop,"('[write_sc_results] Writing new n, mx, my, mz and Ef to file...')")
      write(scfile,"('./results/',a1,'SOC/selfconsistency/selfconsistency_',a,'_dfttype=',a,'_parts=',i0,a,a,a,a,'.dat')") output%SOCchar, trim(output%Sites),dfttype,parts,trim(output%BField),trim(output%info),trim(output%SOC),trim(output%suffix)
      open (unit=99,status='replace',file=scfile)

      
      do i=1,s%nAtoms
        write(formatvar,fmt="(a,i0,a)") '(',5*s%Types(s%Basis(i)%Material)%nOrb,'(es21.11,2x))'
        write(99,fmt=formatvar) (rho(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb),(mx(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb),(my(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb),(mz(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb),(delta_sc(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb)
      end do
      write(99,"(es21.11,2x,'! Ef  ')") s%Ef
      write(99,"('! n(1:nOrb), mx(1:nOrb), my(1:nOrb), mz(1:nOrb), delta_sc(1:nOrb) per site ')")
      write(99,"('! Ef ')")

      close(99)

      if(present(bc)) then
        ! Writing new results for the constranining fields to file
        write(output%unit_loop,"('[write_sc_results] Writing new bc_x, bc_y, bc_z to file...')")
        bcfile=replaceStr( string=scfile, search="selfconsistency_", substitute="constraniningfield_" )
        open (unit=999,status='replace',file=bcfile)

        do i=1,s%nAtoms
          write(999,fmt="(3(es21.11,2x))") (bc(mu,i), mu=1,3)
        end do

        close(999)
      end if

    end if

    call MPI_Bcast(scfile, len(scfile), MPI_CHARACTER, 0, FieldComm,ierr)
  end subroutine write_sc_results


  subroutine print_sc_step(n,mp,mz,delta_sc,s,fvec)
  !> Writes the intermediate steps or the initial values 
  !> of the self-consistency on screen
    use mod_kind,              only: dp
    use mod_parameters,        only: output,lfixEf
    use mod_System,            only: System_type
    use mod_magnet,            only: iter,lconstraining_field,constr_type,bc
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
          do mu = 1,s%Types(s%Basis(i)%Material)%ndOrb
            fvecsum = fvecsum + fvec(neq_per_atom(i)+mu)
          end do
          formatvar = formatvar // '  fvec(N)= ' // trim(RtoS(fvecsum,"(es16.9)"))
        end if
        kount = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp)
        if(abs(s%Basis(i)%Um)>1.e-8_dp) then
          if(.not.lconstraining_field) then
            if(abs(mp(i))>1.e-8_dp) then
              formatvar = formatvar // '  fvec(Mx)= ' // trim(RtoS(fvec(kount+1),"(es16.9)")) // '  fvec(My)= ' // trim(RtoS(fvec(kount+2),"(es16.9)")) // '  fvec(Mz)= ' // trim(RtoS(fvec(kount+3),"(es16.9)"))
            else
              formatvar = formatvar // '  fvec(Mz)= ' // trim(RtoS(fvec(kount+3),"(es16.9)"))
            end if
          else
            select case(constr_type)
            case(1)
              formatvar = formatvar // '  fvec(|M|)= ' // trim(RtoS(fvec(kount+1),"(es16.9)"))
            case(2)
              formatvar = formatvar // '  fvec(Bx)= ' // trim(RtoS(fvec(kount+1),"(es16.9)")) // '  fvec(By)= ' // trim(RtoS(fvec(kount+2),"(es16.9)")) // '  fvec(Bz)= ' // trim(RtoS(fvec(kount+3),"(es16.9)"))
            end select
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
          if(lconstraining_field) write(output%unit_loop,"(8x,'=>',23x,'Bx=',es16.9,4x,'By=',es16.9,4x,'Bz=',es16.9)") bc(1,i),bc(2,i),bc(3,i)
        else
          write(output%unit_loop,"('Site ',i4,': N=',es16.9,4x,'Mz=',es16.9)") i,n(i),mz(i)
          if(lconstraining_field) write(output%unit_loop,"(8x,'=>',23x,'Bz=',es16.9,4x,'By=',es16.9,4x,'Bz=',es16.9)") bc(1,i),bc(2,i),bc(3,i)
        end if
        if(lsuperCond) then
          write(output%unit_loop,"('Site ',i4,': D_sc =',9(es14.7,2x))") i,(delta_sc(mu,i),mu=1,s%Types(s%Basis(i)%Material)%nOrb)
        end if
      end do
      write(output%unit_loop,"(13x,'Ef=',es16.9)") s%Ef
    end if

  end subroutine print_sc_step


  subroutine sc_eqs(N,x,fvec,iflag)
  !> This subroutine calculates the self-consistency equations
  !>     n  - rho_in  = 0    ! When Un /= 0
  !>     mx - mx_in   = 0
  !>     my - my_in   = 0    ! When Um /= 0
  !>     mz - mz_in   = 0
  !>  sum n - n_total = 0    ! when lfixEf == .false.
  !>
  !> The magnetization equations are substituted by the constraning field
  !> if this option is used.
    use mod_kind,              only: dp
    use mod_parameters,        only: output
    use mod_system,            only: s => sys
    use mod_magnet,            only: iter,maxiter,rho,rhot,rhod,mp,mx,my,mz,mpd,mpd0,mxd,myd,mzd,mzd0,rhod0,rho0,m_fix,m_fix_abs,bc,sb_matrix,mabsd,constr_type
    use mod_Umatrix,           only: update_Umatrix
    use mod_expectation,       only: expectation_values
    use mod_superconductivity, only: lsuperCond,update_delta_sc
    implicit none
    integer  :: N,i,mu,mud,iflag
    real(dp),    dimension(N)               :: x,fvec
    real(dp),    dimension(s%nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(s%nAtoms)        :: mxd_in,myd_in,mzd_in,mabsd_in,rhod_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: delta_sc_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: deltas
    complex(dp), dimension(s%nAtoms)        :: mpd_in

    external :: endTITAN

    iflag=0

    ! Values used in the hamiltonian
    rho_in = rho ! To store the non-d orbitals into rho_in
    select case(constr_type)
    case(1)
      ! myd_in,mzd_in are dummies that are not used
      call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mabsd_in,myd_in,mzd_in,delta_sc_in)
      do i = 1,s%nAtoms
        mxd_in(i) = m_fix(1,i)*mabsd_in(i)
        myd_in(i) = m_fix(2,i)*mabsd_in(i)
        mzd_in(i) = m_fix(3,i)*mabsd_in(i)
        mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      end do
      ! Update Hubbard term in Hamiltonian
      call update_Umatrix(mzd_in,mzd0,mpd_in,mpd0,rhod_in,rhod0,rho_in,rho0,s)
    case(2)
      call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,bc(1,:),bc(2,:),bc(3,:),delta_sc_in)
      do i = 1,s%nAtoms
        mxd_in(i) = m_fix(1,i)*m_fix_abs(i) 
        myd_in(i) = m_fix(2,i)*m_fix_abs(i) 
        mzd_in(i) = m_fix(3,i)*m_fix_abs(i) 
        mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      end do

      ! Updating field with new constrainig values
      call sb_matrix(s)
    case default
      call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mxd_in,myd_in,mzd_in,delta_sc_in)
      do i = 1,s%nAtoms
        mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      end do
      ! Update Hubbard term in Hamiltonian
      call update_Umatrix(mzd_in,mzd0,mpd_in,mpd0,rhod_in,rhod0,rho_in,rho0,s)
    end select

    ! Update electron-hole coupling in Hamiltonian
    if(lsuperCond) call update_delta_sc(s,delta_sc_in)

    iter = iter + 1

    call expectation_values(s,rho,mp,mx,my,mz,deltas)

    rhod = 0._dp
    mpd  = 0._dp
    mxd  = 0._dp
    myd  = 0._dp
    mzd  = 0._dp
    do i = 1,s%nAtoms
      do mud = 1,s%Types(s%Basis(i)%Material)%ndOrb
        mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
        rhod(i) = rhod(i) + rho(mu,i)
        mxd (i) = mxd (i) + mx (mu,i)
        myd (i) = myd (i) + my (mu,i)
        mzd (i) = mzd (i) + mz (mu,i)
        mabsd(i)= sqrt(mxd(i)**2+myd(i)**2+mzd(i)**2)
      end do
      mpd(i) = cmplx(mxd(i),myd(i),dp)
    end do
    rhot = sum(rho(:,:))

    ! Setting up linear system of equations:
    select case(constr_type)
    case(1)
      call set_system_of_equations(s,N,rhot,rho,mabsd,myd,mzd,deltas,&
                                   rho_in,mabsd_in,myd_in,mzd_in,delta_sc_in,fvec)
      if(lontheflysc) call write_sc_results(rho,mx,my,mz,bc)
    case(2)
      call set_system_of_equations(s,N,rhot,rho,sum(mx(:,:),dim=1),sum(my(:,:),dim=1),sum(mz(:,:),dim=1),deltas,&
                                   rho_in,mxd_in,myd_in,mzd_in,delta_sc_in,fvec)
      if(lontheflysc) call write_sc_results(rho,mx,my,mz,bc)
    case default
      call set_system_of_equations(s,N,rhot,rho,mxd,myd,mzd,deltas,&
                                   rho_in,mxd_in,myd_in,mzd_in,delta_sc_in,fvec)
      if(lontheflysc) call write_sc_results(rho,mx,my,mz)
    end select

    call print_sc_step(rhod,mpd,mzd,deltas,s,fvec)

    if(iter>=maxiter) then
      write(output%unit_loop,"('[sc_eqs] Maximum number of iterations reached!')")
      call endTITAN()
    end if
  end subroutine sc_eqs


  subroutine sc_jac(N,x,fvec,selfconjac,ldfjac,iflag)
  !> This subroutine calculates the jacobian of the system of equations
  !>     n  - rho_in    = 0
  !>     mx - mx_in   = 0
  !>     my - my_in   = 0
  !>     mz - mz_in   = 0
  !>  sum n - n_total = 0
    use mod_kind,              only: dp
    use mod_system,            only: s => sys
    use mod_magnet,            only: iter,mzd0,mpd0,rhod0,rho0,rho,m_fix,m_fix_abs,bc,sb_matrix,constr_type
    use mod_Umatrix,           only: update_Umatrix
    use mod_superconductivity, only: lsuperCond, update_delta_sc
    implicit none
    integer      :: i,N,ldfjac,iflag
    real(dp)     :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(dp),    dimension(s%nAtoms)        :: mxd_in,myd_in,mzd_in,mabsd_in,rhod_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: rho_in
    real(dp),    dimension(s%nOrb,s%nAtoms) :: delta_sc_in
    complex(dp), dimension(s%nAtoms)        :: mpd_in
    !--------------------- begin MPI vars --------------------

    iflag=0

    ! Values used in the hamiltonian
    rho_in = rho ! To store the non-d orbitals into rho_in
    select case(constr_type)
    case(1)
      ! myd_in,mzd_in are dummies that are not used
      call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mabsd_in,myd_in,mzd_in,delta_sc_in)
      do i = 1,s%nAtoms
        mpd_in(i) = cmplx(m_fix(1,i)*mabsd_in(i),m_fix(2,i)*mabsd_in(i),dp)
        mzd_in(i) = m_fix(3,i)*mabsd_in(i)
      end do
      ! Update Hubbard term in Hamiltonian
      call update_Umatrix(mzd_in,mzd0,mpd_in,mpd0,rhod_in,rhod0,rho_in,rho0,s)
    case(2)
      call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,bc(1,:),bc(2,:),bc(3,:),delta_sc_in)
      do i = 1,s%nAtoms
        mxd_in(i) = m_fix(1,i)*m_fix_abs(i) 
        myd_in(i) = m_fix(2,i)*m_fix_abs(i) 
        mzd_in(i) = m_fix(3,i)*m_fix_abs(i) 
        mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      end do

      ! Updating field with new constrainig values
      call sb_matrix(s)
    case default
      call set_hamiltonian_variables(.true.,s,N,x,rho_in,rhod_in,mxd_in,myd_in,mzd_in,delta_sc_in)
      do i = 1,s%nAtoms
        mpd_in(i) = cmplx(mxd_in(i),myd_in(i),dp)
      end do
      ! Update Hubbard term in Hamiltonian
      call update_Umatrix(mzd_in,mzd0,mpd_in,mpd0,rhod_in,rhod0,rho_in,rho0,s)
    end select

    ! Update electron-hole coupling in Hamiltonian
    if(lsuperCond) call update_delta_sc(s,delta_sc_in)

    fvec=fvec

    call calcJacobian_greenfunction(selfconjac, N)

    iter = iter + 1

  end subroutine sc_jac

  subroutine set_hamiltonian_variables(set,s,N,x,rho,rhod,mxd,myd,mzd,delta_sc)
  !> This subroutine transfers the variables from the single array
  !> to the expectation values used in the hamiltonian
  !> Whe set = .false. it does the inverse, to set the initial guess
    use mod_System,            only: System_type
    use mod_superconductivity, only: lsupercond
    use mod_parameters,        only: lfixEf
    use mod_magnet,            only: constr_type
    implicit none
    logical,                                 intent(in)    :: set
    !> Variable to select which way variables are set
    type(System_type),                       intent(inout) :: s
    !> System quantities
    integer,                                 intent(in)    :: N
    !> Number of equations
    real(dp),    dimension(N),               intent(inout) :: x
    !> Unknowns in single array format
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(inout) :: rho
    !> Variables used in the hamiltonian: orbital-dependent density
    real(dp),    dimension(s%nAtoms),        intent(inout) :: mxd,myd,mzd,rhod
    !> Variables used in the hamiltonian: magnetizatin and density
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(inout) :: delta_sc
    !> Variables used in the hamiltonian: superconducting delta

    ! Local variables:
    integer :: i,mu,mud,kount

    if(set) then
      rhod = 0._dp
      mxd  = 0._dp
      myd  = 0._dp
      mzd  = 0._dp
      do i = 1,s%nAtoms
        kount = neq_per_atom(i)
        if(abs(s%Basis(i)%Un)>1.e-8_dp) then
          do mud = 1,s%Types(s%Basis(i)%Material)%ndOrb
            mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
            rho(mu,i) = x(kount+mud)
            rhod(i)= rhod(i) + rho(mu,i)
          end do
          kount = neq_per_atom(i) + s%Types(s%Basis(i)%Material)%ndOrb
        end if
        if(abs(s%Basis(i)%Um)>1.e-8_dp) then
          mxd(i) = x(kount+1)
          if(constr_type/=1) then
            myd(i) = x(kount+2)
            mzd(i) = x(kount+3)
          end if
          kount = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + merge(1,3,constr_type==1)
        end if
        if(lsupercond) then
          do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
            delta_sc(mu,i) = x(kount+mu)
          end do
        end if
      end do
      if(.not.lfixEf) s%Ef = x(N)
    else
      do i = 1, s%nAtoms
        kount = neq_per_atom(i)
        if(abs(s%Basis(i)%Un)>1.e-8_dp) then
          do mu = 1,s%Types(s%Basis(i)%Material)%ndOrb
            x(kount+mu) = rho(s%Types(s%Basis(i)%Material)%dOrbs(mu),i)
          end do
          kount = neq_per_atom(i) + s%Types(s%Basis(i)%Material)%ndOrb
        end if
        if(abs(s%Basis(i)%Um)>1.e-8_dp) then
          x(kount+1) = mxd(i)
          if(constr_type/=1) then
            x(kount+2) = myd(i)
            x(kount+3) = mzd(i)
          end if
          kount = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + merge(1,3,constr_type==1)
        end if
        if(lsuperCond) then
          do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
            x(kount+mu) = delta_sc(mu,i)
          end do
        end if
      end do
      if(.not.lfixEf) x(neq) = s%Ef
    end if

  end subroutine set_hamiltonian_variables


  subroutine set_system_of_equations(s,N,rhot,rho,mxd,myd,mzd,delta_sc,&
                                     rho_in,mxd_in,myd_in,mzd_in,delta_sc_in,fvec)
  !> This subroutine sets the system of equations to be solved
  !> by the non-linear root finder
    use mod_System,            only: System_type
    use mod_superconductivity, only: lsupercond
    use mod_parameters,        only: lfixEf
    use mod_magnet,            only: constr_type
    implicit none
    type(System_type),                       intent(in)  :: s
    !> System quantities
    integer,                                 intent(in)  :: N
    !> Number of equations
    real(dp),                                intent(in)  :: rhot
    !> Total density
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(in)  :: rho,rho_in
    !> Orbital-dependent densities (output and input)
    real(dp),    dimension(s%nAtoms),        intent(in)  :: mxd,myd,mzd,mxd_in,myd_in,mzd_in
    !> Variables used in the hamiltonian: magnetizatin and density
    real(dp),    dimension(s%nOrb,s%nAtoms), intent(in)  :: delta_sc,delta_sc_in
    !> Variables used in the hamiltonian: superconducting delta
    real(dp),    dimension(N),               intent(out) :: fvec
    !> Vector holding the equations to find the zeroes

    ! Local variables:
    integer :: i,mu,mud,kount

    do i = 1,s%nAtoms
      kount = neq_per_atom(i)
      if(abs(s%Basis(i)%Un)>1.e-8_dp) then
        do mud = 1,s%Types(s%Basis(i)%Material)%ndOrb
          mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
          fvec(kount+mud) = rho(mu,i) - rho_in(mu,i)
        end do
        kount = neq_per_atom(i) + s%Types(s%Basis(i)%Material)%ndOrb
      end if
      if(abs(s%Basis(i)%Um)>1.e-8_dp) then
        fvec(kount+1) = mxd(i) - mxd_in(i)
        if(constr_type/=1) then
          fvec(kount+2) = myd(i) - myd_in(i)
          fvec(kount+3) = mzd(i) - mzd_in(i)
        end if
        kount = neq_per_atom(i) + merge(s%Types(s%Basis(i)%Material)%ndOrb,0,abs(s%Basis(i)%Un)>1.e-8_dp) + merge(merge(1,3,constr_type==1),0,abs(s%Basis(i)%Um)>1.e-8_dp)
      end if
      if(lsuperCond) then
        do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
          fvec(kount+mu) = delta_sc(mu,i) - delta_sc_in(mu,i)
        end do
      end if
    end do
    if(.not.lfixEf) fvec(N) = rhot - s%totalOccupation
  end subroutine set_system_of_equations

end module mod_self_consistency
