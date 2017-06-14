module mod_self_consistency
  implicit none
  character(len=300)  :: default_file

contains
  ! Tries to read eps1 and m if available - includes hdel, hdelp and hdelm calculations
  subroutine read_previous_results(lsuccess)
    use mod_f90_kind
    use mod_constants
    use mod_magnet
    use mod_parameters
    use mod_system, only: pln_a1, pln_a2, c_nn, pln_cnt
    use mod_mpi_pars, only: myrank_row_hw, ierr, myrank
    use mod_susceptibilities, only: lrot
    implicit none
    integer             :: i,err
    logical,intent(out) :: lsuccess

    lsuccess = .false.
    call read_write_sc_results(0,err,lsuccess)

    if(lsuccess) then
      if(err==0) then ! Same parameters
        if(skipsc) then ! Skip option ON
          if(myrank_row_hw==0) write(outputunit_loop,"('[read_previous_results] Existing results for the same parameters were read. Skipping self-consistency...')")
          lselfcon = .false.
        else ! Skip option OFF
          if(myrank_row_hw==0) write(outputunit_loop,"('[read_previous_results] Existing results for the same parameters were read. Updating values...')")
          lselfcon = .true.
        end if
      else ! Other parameters
        if(myrank_row_hw==0) write(outputunit_loop,"('[read_previous_results] Existing results for other parameters were read. Updating values...')")
        lselfcon = .true.
      end if
      ! Calculating angles of GS magnetization in units of pi and magnetization vector
      do i = 1,Npl
        mabs(i)   = sqrt((mx(i)**2)+(my(i)**2)+(mz(i)**2))
        mtheta(i) = acos(mz(i)/mabs(i))/pi
        if(abs(mtheta(i))>1.d-8) then
          if(abs(abs(mtheta(i))-1.d0)>1.d-8) then
            mphi(i)   = atan2(my(i),mx(i))/pi
          else
            mphi(i) = 0.d0
          end if
          lrot = .true. ! Susceptibilities need to be rotated
        else
          mphi(i) = 0.d0
        end if
        mvec_cartesian(i,1) = mx(i)
        mvec_cartesian(i,2) = my(i)
        mvec_cartesian(i,3) = mz(i)
        mvec_spherical(i,1) = mabs(i)
        mvec_spherical(i,2) = mtheta(i)
        mvec_spherical(i,3) = mphi(i)
      end do
    else !  If file doesn't exist
      if(myrank_row_hw==0) then
        write(outputunit_loop,"('[read_previous_results] Self-consistency file does not exist:')")
        write(outputunit_loop,"('[read_previous_results] ',a)") trim(default_file)
      end if
      lselfcon = .true.
      ! Parameters: center of band, magnetization, exchange split
      eps1 = 0.d0
      if(magaxis == -1) then
        magaxisvec = magaxisvec / sqrt(dot_product(magaxisvec, magaxisvec))
      else if(magaxis == -2) then
        magaxisvec = magaxisvec(1) * pln_a1 + magaxisvec(2) * pln_a2
        magaxisvec = magaxisvec / sqrt(dot_product(magaxisvec, magaxisvec))
      else if(magaxis == -3) then
        magaxisvec = [cos(magaxisvec(1))*sin(magaxisvec(2)), sin(magaxisvec(1))*sin(magaxisvec(2)), cos(magaxisvec(2))]
      else if(magaxis == 0) then
        magaxisvec = [0.d0, 0.d0, sign(1.0d0, hw_list(hw_count,1))]
      else if(magaxis >=1 .and. magaxis <= pln_cnt(1)) then
        magaxisvec(1:3) = c_nn(1:3, magaxis)
      else
        if(myrank.eq.0) write(outputunit,"('[read_previous_results] Unknown magnetization direction!')")
        call MPI_Finalize(ierr)
        stop
      end if

      magaxisvec = magaxisvec * 0.5d0

      mx = magaxisvec(1)
      my = magaxisvec(2)
      mz = magaxisvec(3)


      do i=1,Npl
        if(layertype(i+offset)==2) then
          mx = mx * sign(4.d0,hw_list(hw_count,1))
          my = my * sign(4.d0,hw_list(hw_count,1))
          mz = mz * sign(4.d0,hw_list(hw_count,1))
        endif
      end do

      mp = zero
      if(lfield .and. magaxis == 0) then
        mx = mz*sin(hw_list(hw_count,2)*pi)*cos(hw_list(hw_count,3)*pi)
        my = mz*sin(hw_list(hw_count,2)*pi)*sin(hw_list(hw_count,3)*pi)
        mz = mz*cos(hw_list(hw_count,2)*pi)
      end if
      mp = cmplx(mx,my,double)

      ! Variables used in the hamiltonian
      do i=1,Npl
        hdel(i)   = 0.5d0*U(i+offset)*mz(i)
        hdelp(i)  = 0.5d0*U(i+offset)*mp(i)
      end do
      hdelm = conjg(hdelp)
    end if

    return
  end subroutine read_previous_results

  ! Rotate the magnetization to the direction of the field (useful for SOC=F)
  subroutine rotate_magnetization_to_field()
    use mod_constants
    use mod_magnet
    use mod_parameters
    use mod_mpi_pars, only: myrank_row_hw
    use mod_tight_binding
    implicit none
    integer :: i,err,sign
    logical :: lsuccess
    real(double) :: mdotb

    if(myrank_row_hw==0) write(outputunit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    do i=1,Npl
      mdotb   = hhwx(i)*mx(i)+hhwy(i)*my(i)+hhwz(i)*mz(i)
      sign    = dble(mdotb/abs(mdotb))
      mabs(i) = sqrt(abs(mp(i))**2+(mz(i)**2))
      mx(i)   = sign*mabs(i)*sin(hw_list(hw_count,2)*pi)*cos(hw_list(hw_count,3)*pi)
      my(i)   = sign*mabs(i)*sin(hw_list(hw_count,2)*pi)*sin(hw_list(hw_count,3)*pi)
      mz(i)   = sign*mabs(i)*cos(hw_list(hw_count,2)*pi)
      mp(i)   = cmplx(mx(i),my(i),double)
    end do

    ! Writing new eps1 and rotated mag to file (without self-consistency)
    if(myrank_row_hw==0) call read_write_sc_results(1,err,lsuccess)

    ! Writing self-consistency results on screen
    if(myrank_row_hw==0)  call write_sc_results_on_screen()

    deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1)
    deallocate(mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm)
    deallocate(mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi)
    if(lGSL) deallocate(lxm,lym,lzm,lxpm,lypm,lzpm)
    deallocate(mmlayer,layertype,U,mmlayermag,lambda,npart0)
    deallocate(t0, t0i)

    return
  end subroutine rotate_magnetization_to_field

  ! This subroutine performs the self-consistency
  subroutine do_self_consistency()
    use mod_f90_kind
    use mod_constants, only: pi
    use mod_parameters
    use mod_magnet
    use mod_mpi_pars, only: myrank_row_hw,mpitag
    use mod_dnsqe
    use mod_susceptibilities, only: lrot
    implicit none
    real(double),allocatable      :: fvec(:),jac(:,:),wa(:),sc_solu(:)
    real(double),allocatable      :: diag(:),qtf(:),w(:,:)
    real(double)                  :: epsfcn,factor
#if !defined(_OSX) && !defined(_JUQUEEN)
    real(double)                  :: ruser(1)
    integer                       :: iuser(1)
#endif
    integer                       :: i,neq,maxfev,ml,mr,mode,nfev,njev,lwa,ifail=0

    neq = 4*Npl
    allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )

    ! Putting read eps1 existing solutions into esp1_solu (first guess of the subroutine)
    sc_solu(1:Npl)         = eps1
    sc_solu(Npl+1:2*Npl)   = mx
    sc_solu(2*Npl+1:3*Npl) = my
    sc_solu(3*Npl+1:4*Npl) = mz

    iter  = 1
    mpitag = (Npl-Npl_i)*total_hw_npt1 + hw_count
    if(myrank_row_hw==0) write(outputunit_loop,"('[self_consistency] Starting self-consistency:')")

#if defined(_OSX) || defined(_JUQUEEN)
    if(lslatec) then
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa),w(neq,4) )
      if(lnojac) then
        call dnsqe(sc_eqs_old,sc_jac_old,2,neq,sc_solu,fvec,tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs_old,sc_jac_old,1,neq,sc_solu,fvec,tol,0,ifail,wa,lwa)
      end if
      ifail = ifail-1
    else
      lwa=neq*(neq+1)/2
      allocate( wa(lwa),w(neq,4) )
      if(lnojac) then
!         call c05nbf(sc_equations,neq,sc_solu,fvec,tol,wa,lwa,ifail)
        maxfev = 200*(neq+1)
        ml = neq-1
        mr = neq-1
        epsfcn = 1.d-5
        mode = 1
        factor = 100.d0
        call c05ncf(sc_eqs_old,neq,sc_solu,fvec,tol,maxfev,ml,mr,epsfcn,diag,mode,factor,0,nfev,jac,neq,wa,lwa,qtf,w,ifail)
      else
!         call c05pbf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,neq,tol,wa,lwa,ifail)
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
!         diag(Npl+1:4*Npl) = 100.d0
        factor = 100.d0
        call c05pcf(sc_eqs_and_jac_old,neq,sc_solu,fvec,jac,neq,tol,maxfev,diag,mode,factor,0,nfev,njev,wa,lwa,qtf,w,ifail)
      end if
    end if
    deallocate( w )
#else
    if(lslatec) then
      lwa=neq*(3*neq+13)/2
      allocate( wa(lwa) )
      if(lnojac) then
        call dnsqe(sc_eqs_old,sc_jac_old,2,neq,sc_solu,fvec,tol,0,ifail,wa,lwa)
      else
        call dnsqe(sc_eqs_old,sc_jac_old,1,neq,sc_solu,fvec,tol,0,ifail,wa,lwa)
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
        call c05qcf(sc_equations,neq,sc_solu,fvec,tol,maxfev,ml,mr,epsfcn,mode,diag,factor,0,nfev,jac,wa,qtf,iuser,ruser,ifail)
      else
        maxfev = 100*(neq+1)
        mode = 1
        diag = 1.d0
!         diag(npl+1:4*npl) = 100.d0
        factor = 100.d0
        call c05rcf(sc_equations_and_jacobian,neq,sc_solu,fvec,jac,tol,maxfev,mode,diag,factor,0,nfev,njev,wa,qtf,iuser,ruser,ifail)
      end if
    end if
#endif

    deallocate(sc_solu,diag,qtf,fvec,jac,wa)

    ! Calculating new angles of GS magnetization in units of pi and magnetization vector
    do i = 1,Npl
      mabs(i)   = sqrt((mx(i)**2)+(my(i)**2)+(mz(i)**2))
      mtheta(i) = acos(mz(i)/mabs(i))/pi
      if(abs(mtheta(i))>1.d-8) then
        if(abs(abs(mtheta(i))-1.d0)>1.d-8) then
          mphi(i)   = atan2(my(i),mx(i))/pi
        else
          mphi(i) = 0.d0
        end if
        lrot = .true. ! Susceptibilities need to be rotated
      else
        mphi(i) = 0.d0
      end if
      mvec_cartesian(i,1) = mx(i)
      mvec_cartesian(i,2) = my(i)
      mvec_cartesian(i,3) = mz(i)
      mvec_spherical(i,1) = mabs(i)
      mvec_spherical(i,2) = mtheta(i)
      mvec_spherical(i,3) = mphi(i)
    end do

    return
  end subroutine do_self_consistency

  ! Writes the self-consistency results into files
  subroutine write_sc_results()
    use mod_parameters, only: scfile,hw_count,total_hw_npt1
    use mod_mpi_pars
    implicit none
    integer :: err
    logical :: lsuccess

    if(myrank_row_hw==0) call read_write_sc_results(1,err,lsuccess)
    if(hw_count/=total_hw_npt1) then
      call MPI_Bcast(scfile,len(scfile),MPI_CHARACTER,0,MPI_Comm_Row_hw,ierr)
    else
      scfile = ""
    end if

    return
  end subroutine write_sc_results

  ! Writes the self-consistency results on the screen
  subroutine write_sc_results_on_screen()
    use mod_parameters, only: Npl,outputunit_loop,lGSL
    use mod_mpi_pars
    use mod_magnet
    implicit none
    integer :: i

    write(outputunit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    write(outputunit_loop,"(11x,' *************** Center of d bands: ***************')")
    do i=1,Npl
      write(outputunit_loop,"(26x,'eps1(',i2.0,')=',f11.8)") i,eps1(i)
    end do
    write(outputunit_loop,"(11x,' *********** Magnetization components: **********')")
    do i=1,Npl
      write(outputunit_loop,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mx(i),i,my(i),i,mz(i)
      if(abs(mp(i))/=0) write(outputunit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") mtheta(i),mphi(i)
    end do
    if(lGSL) then
      write(outputunit_loop,"(11x,' *** Orbital components in spin coordinates:  ***')")
      do i=1,Npl
        write(outputunit_loop,"(4x,'Lxp(',i2.0,')=',f11.8,4x,'Lyp(',i2.0,')=',f11.8,4x,'Lzp(',i2.0,')=',f11.8)") i,lxpm(i),i,lypm(i),i,lzpm(i)
        if(sqrt(lxpm(i)**2+lypm(i)**2)/=0) write(outputunit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") lptheta(i),lpphi(i)
      end do
      write(outputunit_loop,"(11x,' *** Orbital components in cubic coordinates: ***')")
      do i=1,Npl
        write(outputunit_loop,"(4x,'Lx (',i2.0,')=',f11.8,4x,'Ly (',i2.0,')=',f11.8,4x,'Lz (',i2.0,')=',f11.8)") i,lxm(i),i,lym(i),i,lzm(i)
      end do
      write(outputunit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,Npl
        write(outputunit_loop,"(4x,'M (',i2.0,') =',f11.8,4x,'Lp (',i2.0,')=',f11.8,4x,'L (',i2.0,') =',f11.8)") i,mabs(i),i,lpabs(i),i,labs(i)
      end do
    else
      write(outputunit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,Npl
        write(outputunit_loop,"(27x,'M (',i2.0,') =',f11.8)") i,mabs(i)
      end do
    end if
    write(outputunit_loop,"('|----------=============================================================----------|')")

    return
  end subroutine write_sc_results_on_screen

  ! This subroutine reads previous band-shifting and magnetization results
  ! and also writes new ones into file
  subroutine read_write_sc_results(iflag,err,lsuccess)
    use mod_constants
    use mod_parameters
    use mod_magnet
    use mod_mpi_pars
    use mod_system, only: nkpt
    implicit none
    character(len=300)  :: file = ""
    character(len=100)  :: fieldpart,socpart,folder,prefix
    character(len=1)    :: SOCc
    integer,intent(in)  :: iflag
    integer,intent(out) :: err
    logical,intent(out) :: lsuccess
    integer             :: i
    real(double)        :: previous_results(Npl,4)
    if((trim(scfile)/="").and.(iflag==0)) then
      open(unit=99,file=scfile,status="old",iostat=err)
      if(err/=0) then
        if(myrank_row_hw==0) write(outputunit_loop,"('*** WARNING: Self-consistency file given on input file does not exist! Using default... ***')")
        scfile = " "
      end if
      close(99)
    end if
    !folder    = "./results/selfconsistency/"
    prefix    = "selfconsistency_"
    fieldpart = ""
    socpart   = ""
    if(lfield) then
      write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
      if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
      if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
      if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
      if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"
    end if
    if(SOC) then
      if(llinearsoc) then
        SOCc = "L"
      else
        SOCc = "T"
      end if
      if(abs(socscale-1.d0)>1.d-6) write(socpart,"('_socscale=',f5.2)") socscale
      if((llineargfsoc).or.(llinearsoc)) socpart = trim(socpart) // "_linearsoc"
    else
      SOCc = "F"
    end if
    lsuccess = .false.
  !   Reading previous results (mx, my, mz and eps1) from files (if available)
    if(iflag==0) then
      if(trim(scfile)=="") then ! If a filename is not given in inputcard (or don't exist), use the default one
        write(file,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts,Utype,trim(fieldpart),nkpt,eta,trim(socpart)
        open(unit=99,file=file,status="old",iostat=err)
        if((err==0).and.(myrank_row_hw==0)) then
          write(outputunit_loop,"('[read_write_sc_results] Self-consistency file already exists. Reading it now...')")
          write(outputunit_loop,"(a)") trim(file)
        else
          default_file = trim(file)
        end if
      else ! If filename in inputcard exists or 2nd+ angular iteration
        if(((hw_count)==1).and.(Npl==Npl_i)) then ! Filename in inputcard (1st iteration on loop)
          open(unit=99,file=scfile,status="old",iostat=err)
          if((err==0).and.(myrank_row_hw==0)) then
            write(outputunit_loop,"('[read_write_sc_results] Using filename given in input file for self-consistency:')")
            write(outputunit_loop,"(a)") trim(scfile)
          end if
        else ! 2nd+ iteration, cheking if default file exists
          write(file,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts,Utype,trim(fieldpart),nkpt,eta,trim(socpart)
          open(unit=99,file=file,status="old",iostat=err)
          if(err==0) then ! Reading file for the same parameters
            if(myrank_row_hw==0) then
              write(outputunit_loop,"('[read_write_sc_results] Self-consistency file already exists. Reading it now...')")
              write(outputunit_loop,"(a)") trim(file)
            end if
          else ! Reading results from previous iteration
            open(unit=99,file=scfile,status="old",iostat=err)
            if((err==0).and.(myrank_row_hw==0)) then
              write(outputunit_loop,"('[read_write_sc_results] Using results from previous iteration as input for self-consistency:')")
              write(outputunit_loop,"(a)") trim(scfile)
            end if
            lsuccess   = .true. ! something was read
          end if
        end if
      end if
      if(err==0) then
        if(myrank_row_hw==0) then
          do i=1,Npl
            read(99,fmt=*) previous_results(i,1)
            read(99,fmt=*) previous_results(i,2)
            read(99,fmt=*) previous_results(i,3)
            read(99,fmt=*) previous_results(i,4)
          end do
        end if
        call MPI_Bcast(previous_results,4*Npl,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
        eps1(:) = previous_results(:,1)
        mx  (:) = previous_results(:,2)
        my  (:) = previous_results(:,3)
        mz  (:) = previous_results(:,4)
        mp  = mx + zi*my
        do i=1,Npl
          hdel(i)   = 0.5d0*U(i+offset)*mz(i)
          hdelp(i)  = 0.5d0*U(i+offset)*mp(i)
        end do
        hdelm = conjg(hdelp)
        if(lsuccess) then
          err = 1   ! Read different parameters
        else
          lsuccess   = .true. ! Read same parameters (err=0)
        end if
      else
        ! If file does not exist, try to read for parts-1
        close(99)
        write(file,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts-1,Utype,trim(fieldpart),nkpt,eta,trim(socpart)
        open(unit=99,file=file,status="old",iostat=err)
        if(err==0) then
          if(myrank_row_hw==0) then
            write(outputunit_loop,"('[read_write_sc_results] Self-consistency file does not exist. Reading results for parts-1 now...')")
            write(outputunit_loop,"('[read_write_sc_results] Updating values obtained for parts-1...')")
            write(outputunit_loop,"(a)") file
          end if
          do i=1,Npl
            read(99,*) eps1(i)
            read(99,*) mx(i)
            read(99,*) my(i)
            read(99,*) mz(i)
          end do
          mp  = mx + zi*my
          do i=1,Npl
            hdel(i)   = 0.5d0*U(i+offset)*mz(i)
            hdelp(i)  = 0.5d0*U(i+offset)*mp(i)
          end do
          hdelm = conjg(hdelp)
          lsuccess = .true. ! Read...
          err = 1           ! ... different parameters
        end if
      end if
      close(99)
!       if((iflag==0).and.(myrank_row_hw==0)) then
!         write(outputunit_loop,"('[read_write_sc_results] File not found:')")
!         write(outputunit_loop,"(a)") trim(scfile)
!       end if
!     Writing new results (mx, my, mz and eps1) and mz to file
    else
      write(outputunit_loop,"('[read_write_sc_results] Writing new eps1, mx, my and mz to file...')")
      write(scfile,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts,Utype,trim(fieldpart),nkpt,eta,trim(socpart)
      open (unit=99,status='replace',file=scfile)
      do i=1,Npl
        write(99,"(es21.11,2x,'! eps1')") eps1(i)
        write(99,"(es21.11,2x,'! mx  ')") mx(i)
        write(99,"(es21.11,2x,'! my  ')") my(i)
        write(99,"(es21.11,2x,'! mz  ')") mz(i)
      end do
      close(99)
    end if
    return
  end subroutine read_write_sc_results

  ! This subroutine calculates the self-consistency equations
  !  n  - n0    = 0
  !  mx - mx_in = 0
  !  my - my_in = 0
  !  mz - mz_in = 0
  ! and the correspondent jacobian
#if !defined(_OSX) && !defined(_JUQUEEN)
  subroutine sc_equations_and_jacobian(N,x,fvec,selfconjac,iuser,ruser,iflag)
    use mod_constants
    use mod_parameters
    use mod_f90_kind
    use mod_generate_epoints
    use mod_magnet
    use mod_tight_binding, only: npart0
    use mod_progress
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag
    integer,           intent(inout) :: iuser(*)
    real(double),      intent(inout) :: ruser(*)
    real(double),dimension(N)        :: x,fvec
    real(double),dimension(N,N)      :: selfconjac
    real(double),dimension(Npl)      :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(Npl)   :: mp_in
    real(double),dimension(N,N)      :: ggr
    real(double),dimension(Npl,9)    :: n_orb_u,n_orb_d,n_orb_t,mag_orb
    real(double),dimension(Npl,9)    :: gdiaguur,gdiagddr
    complex(double),dimension(Npl,9) :: gdiagud,gdiagdu
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount,ncount2
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount=Npl*9
    ncount2=N*N

  ! Values used in the hamiltonian
    eps1  = x(1:Npl)
    mx_in = x(Npl+1:2*Npl)
    my_in = x(2*Npl+1:3*Npl)
    mz_in = x(3*Npl+1:4*Npl)
    mp_in = mx_in+zi*my_in

    do i=1,Npl
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,Npl
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    ix = myrank_row_hw+1
    itask = numprocs ! Number of tasks done initially

    flag: select case (iflag)
    case(1)
      ! Calculating the number of particles for each spin and orbital using a complex integral
      if (myrank_row_hw==0) then ! Process 0 receives all results and send new tasks if necessary
        write(outputunit_loop,"('|------------------- Iteration ',i4,' (densities) ------------------|')") iter
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

        if(lverbose) write(outputunit_loop,"('[sc_equations_and_jacobian1] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
        do i=2,pn1
          if(lverbose) call progress_bar(outputunit_loop,"densities energy points",i,pn1)

          call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),8998+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),7997+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),6996+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          if(lverbose) write(outputunit_loop,"('[sc_equations_and_jacobian1] Point ',i0,' received from ',i0)") i,stat(MPI_SOURCE)

          n_orb_u = n_orb_u + gdiaguur
          n_orb_d = n_orb_d + gdiagddr

          do j=1,Npl
            mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
          end do

          ! If the number of processors is less than the total number of points, sends
          ! the rest of the points to the ones that finish first
          if (itask<pn1) then
            itask = itask + 1
            call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_Comm_Row_hw,ierr)
          else
            call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_Comm_Row_hw,ierr)
          end if
        end do
      else
        ! Other processors calculate each point of the integral and waits for new points
        do
          if(ix>pn1) exit

          ! First and second integrations (in the complex plane)
          call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
          gdiaguur = wght(ix)*gdiaguur
          gdiagddr = wght(ix)*gdiagddr
          gdiagud = wght(ix)*gdiagud
          gdiagdu = wght(ix)*gdiagdu

  !         if(lverbose) write(outputunit_loop,"('[selfconsistencyjacnag1] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
          ! Sending results to process 0
          call MPI_Send(gdiaguur,ncount,MPI_DOUBLE_PRECISION,0,9999+iter+mpitag,MPI_Comm_Row_hw,ierr)
          call MPI_Send(gdiagddr,ncount,MPI_DOUBLE_PRECISION,0,8998+iter+mpitag,MPI_Comm_Row_hw,ierr)
          call MPI_Send(gdiagud,ncount,MPI_DOUBLE_COMPLEX,0,7997+iter+mpitag,MPI_Comm_Row_hw,ierr)
          call MPI_Send(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,0,6996+iter+mpitag,MPI_Comm_Row_hw,ierr)
          ! Receiving new point or signal to exit
          call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_Comm_Row_hw,stat,ierr)
          if(ix==0) exit
        end do
      end if

      ! Send results to all processors
      call MPI_Bcast(n_orb_u,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
      call MPI_Bcast(n_orb_d,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
      call MPI_Bcast(mp,Npl,MPI_DOUBLE_COMPLEX,0,MPI_Comm_Row_hw,ierr)

      n_orb_u = 0.5d0 + n_orb_u/pi
      n_orb_d = 0.5d0 + n_orb_d/pi
      n_orb_t = n_orb_u + n_orb_d
      mag_orb = n_orb_u - n_orb_d
      mp      = mp/pi
      mx      = real(mp)
      my      = aimag(mp)

      do i=1,Npl
        ! Number of particles
        n_t(i) = sum(n_orb_t(i,:))
        fvec(i)   = n_t(i) - npart0(i+offset)
        ! x-component of magnetization
        j = i+Npl
        fvec(j)  = mx(i) - mx_in(i)
        ! y-component of magnetization
        j = j+Npl
        fvec(j)  = my(i) - my_in(i)
        ! z-component of magnetization
        j = j+Npl
        mz(i)    = sum(mag_orb(i,5:9))
        fvec(j)  = mz(i) - mz_in(i)
      end do

      if(myrank_row_hw==0) then
        do i=1,Npl
          if(abs(mp(i))>1.d-10) then
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+Npl,fvec(i+Npl),i+2*Npl,fvec(i+2*Npl),i+3*Npl,fvec(i+3*Npl)
          else
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*Npl,fvec(i+3*Npl)
          end if
        end do
      end if
      if(lontheflysc) call write_sc_results()
    case(2)
      selfconjac = 0.d0
      ! Calculating the jacobian
      if (myrank_row_hw==0) then ! Process 0 receives all results and send new tasks if necessary
        write(outputunit_loop,"('|------------------- Iteration ',i4,' (jacobian) -------------------|')") iter
        call sumk_jacobian(Ef,y(ix),ggr)
        selfconjac = wght(ix)*ggr

        if(lverbose) write(outputunit_loop,"('[sc_equations_and_jacobian2] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
        do i=2,pn1
          if(lverbose) call progress_bar(outputunit_loop,"jacobian energy points",i,pn1)

          call MPI_Recv(ggr,ncount2,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3333+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          if(lverbose) write(outputunit_loop,"('[sc_equations_and_jacobian2] Point ',i0,' received from ',i0)") i,stat(MPI_SOURCE)

          selfconjac = selfconjac + ggr

          ! If the number of processors is less than the total number of points, sends
          ! the rest of the points to the ones that finish first
          if (itask<pn1) then
            itask = itask + 1
            call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_Comm_Row_hw,ierr)
          else
            call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_Comm_Row_hw,ierr)
          end if
        end do
      else
        ! Other processors calculate each point of the integral and waits for new points
        do
          if(ix>pn1) exit

          ! First and second integrations (in the complex plane)
          call sumk_jacobian(Ef,y(ix),ggr)
          ggr = wght(ix)*ggr

  !         if(lverbose) write(outputunit_loop,"('[sc_equations_and_jacobian2] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
          ! Sending results to process 0
          call MPI_Send(ggr,ncount2,MPI_DOUBLE_PRECISION,0,3333+iter+mpitag,MPI_Comm_Row_hw,ierr)
          ! Receiving new point or signal to exit
          call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_Comm_Row_hw,stat,ierr)
          if(ix==0) exit
        end do
      end if

      ! Send results to all processors
      call MPI_Bcast(selfconjac,ncount2,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)

      selfconjac = selfconjac/pi
      do i = Npl+1,4*Npl
        selfconjac(i,i) = selfconjac(i,i) - 1.d0
      end do
    case default
      write(outputunit,"('[sc_equations_and_jacobian] Problem in self-consistency! iflag = ',I0)") iflag
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select flag

    iter = iter + 1

    return
  end subroutine sc_equations_and_jacobian

  ! For a given value of center of band eps1 it calculates the
  ! occupation number and the magnetic moment
  subroutine sc_equations(N,x,fvec,iuser,ruser,iflag)
    use mod_constants
    use mod_parameters
    use mod_f90_kind
    use mod_generate_epoints
    use mod_magnet
    use mod_tight_binding, only: npart0
    use mod_progress
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag
    integer     , intent(inout)         :: iuser(*)
    real(double), intent(inout)         :: ruser(*)
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

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:Npl)
    mx_in = x(Npl+1:2*Npl)
    my_in = x(2*Npl+1:3*Npl)
    mz_in = x(3*Npl+1:4*Npl)
    mp_in = mx_in+zi*my_in
    do i=1,Npl
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,Npl
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    n_orb_u = 0.d0
    n_orb_d = 0.d0

    ix = myrank_row_hw+1
    itask = numprocs ! Number of tasks done initially

    ! Calculating the number of particles for each spin and orbital using a complex integral
    if (myrank_row_hw==0) then ! Process 0 receives all results and send new tasks if necessary
      write(outputunit_loop,"('|------------------- Iteration ',i4,' (densities) ------------------|')") iter
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

      if(lverbose) write(outputunit_loop,"('[sc_equations] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
      do i=2,pn1
        if(lverbose) call progress_bar(outputunit_loop,"densities energy points",i,pn1)

        call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
        call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),9998+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
        call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9997+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
        call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9996+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)

        n_orb_u = n_orb_u + gdiaguur
        n_orb_d = n_orb_d + gdiagddr

        do j=1,Npl
          mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
        end do

        ! If the number of processors is less than the total number of points, sends
        ! the rest of the points to the ones that finish first
        if (itask<pn1) then
          itask = itask + 1
          call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_Comm_Row_hw,ierr)
        else
          call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_Comm_Row_hw,ierr)
        end if
      end do
    else
      ! Other processors calculate each point of the integral and waits for new points
      do
        if(ix>pn1) exit

        ! First and second integrations (in the complex plane)
        call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
        gdiaguur = wght(ix)*gdiaguur
        gdiagddr = wght(ix)*gdiagddr
        gdiagud = wght(ix)*gdiagud
        gdiagdu = wght(ix)*gdiagdu

  !       if(lverbose) write(outputunit_loop,"('[sc_equations] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
        ! Sending results to process 0
        call MPI_Send(gdiaguur,ncount,MPI_DOUBLE_PRECISION,0,9999+iter+mpitag,MPI_Comm_Row_hw,ierr)
        call MPI_Send(gdiagddr,ncount,MPI_DOUBLE_PRECISION,0,9998+iter+mpitag,MPI_Comm_Row_hw,ierr)
        call MPI_Send(gdiagud,ncount,MPI_DOUBLE_COMPLEX,0,9997+iter+mpitag,MPI_Comm_Row_hw,ierr)
        call MPI_Send(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,0,9996+iter+mpitag,MPI_Comm_Row_hw,ierr)
        ! Receiving new point or signal to exit
        call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_Comm_Row_hw,stat,ierr)
        if(ix==0) exit
      end do
    end if

    ! Send results to all processors
    call MPI_Bcast(n_orb_u,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
    call MPI_Bcast(n_orb_d,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
    call MPI_Bcast(mp,Npl,MPI_DOUBLE_COMPLEX,0,MPI_Comm_Row_hw,ierr)

    n_orb_u = 0.5d0 + n_orb_u/pi
    n_orb_d = 0.5d0 + n_orb_d/pi
    n_orb_t = n_orb_u + n_orb_d
    mag_orb = n_orb_u - n_orb_d
    mp      = mp/pi
    mx      = real(mp)
    my      = aimag(mp)

    do i=1,Npl
      ! Number of particles
      n_t(i) = sum(n_orb_t(i,:))
      fvec(i)   = n_t(i) - npart0(i+offset)
      ! x-component of magnetization
      j = i+Npl
      fvec(j)  = mx(i) - mx_in(i)
      ! y-component of magnetization
      j = j+Npl
      fvec(j)  = my(i) - my_in(i)
      ! z-component of magnetization
      j = j+Npl
      mz(i)    = sum(mag_orb(i,5:9))
      fvec(j)  = mz(i) - mz_in(i)
    end do

    if(myrank_row_hw==0) then
      do i=1,Npl
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+Npl,fvec(i+Npl),i+2*Npl,fvec(i+2*Npl),i+3*Npl,fvec(i+3*Npl)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*Npl,fvec(i+3*Npl)
        end if
      end do
    end if

    if(lontheflysc) call write_sc_results()

    iter = iter + 1

    return
  end subroutine sc_equations

#endif

  ! This subroutine calculates the self-consistency equations
  !  n  - n0    = 0
  !  mx - mx_in = 0
  !  my - my_in = 0
  !  mz - mz_in = 0
  ! and the correspondent jacobian
  subroutine sc_eqs_and_jac_old(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_constants
    use mod_parameters
    use mod_f90_kind
    use mod_generate_epoints
    use mod_magnet
    use mod_tight_binding, only: npart0
    use mod_progress
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag
    real(double),dimension(N)        :: x,fvec
    real(double),dimension(ldfjac,N) :: selfconjac
    real(double),dimension(Npl)      :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(Npl)   :: mp_in
    real(double),dimension(N,N)      :: ggr
    real(double),dimension(Npl,9)    :: n_orb_u,n_orb_d,n_orb_t,mag_orb
    real(double),dimension(Npl,9)    :: gdiaguur,gdiagddr
    complex(double),dimension(Npl,9) :: gdiagud,gdiagdu
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount,ncount2
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount=Npl*9
    ncount2=N*N

  ! Values used in the hamiltonian
    eps1  = x(1:Npl)
    mx_in = x(Npl+1:2*Npl)
    my_in = x(2*Npl+1:3*Npl)
    mz_in = x(3*Npl+1:4*Npl)
    mp_in = mx_in+zi*my_in
    do i=1,Npl
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,Npl
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    ix = myrank_row_hw+1
    itask = numprocs ! Number of tasks done initially

    flag: select case (iflag)
    case(1)
      ! Calculating the number of particles for each spin and orbital using a complex integral
      if (myrank_row_hw==0) then ! Process 0 receives all results and send new tasks if necessary
        write(outputunit_loop,"('|------------------- Iteration ',i4,' (densities) ------------------|')") iter
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

        if(lverbose) write(outputunit_loop,"('[sc_eqs_and_jac_old1] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
        do i=2,pn1
          if(lverbose) call progress_bar(outputunit_loop,"densities energy points",i,pn1)

          call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),8998+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),7997+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),6996+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          if(lverbose) write(outputunit_loop,"('[sc_eqs_and_jac_old1] Point ',i0,' received from ',i0)") i,stat(MPI_SOURCE)

          n_orb_u = n_orb_u + gdiaguur
          n_orb_d = n_orb_d + gdiagddr

          do j=1,Npl
            mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
          end do

          ! If the number of processors is less than the total number of points, sends
          ! the rest of the points to the ones that finish first
          if (itask<pn1) then
            itask = itask + 1
            call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_Comm_Row_hw,ierr)
          else
            call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_Comm_Row_hw,ierr)
          end if
        end do
      else
        ! Other processors calculate each point of the integral and waits for new points
        do
          if(ix>pn1) exit

          ! First and second integrations (in the complex plane)
          call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
          gdiaguur = wght(ix)*gdiaguur
          gdiagddr = wght(ix)*gdiagddr
          gdiagud = wght(ix)*gdiagud
          gdiagdu = wght(ix)*gdiagdu

  !         if(lverbose) write(outputunit_loop,"('[selfconsistencyjacnag1] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
          ! Sending results to process 0
          call MPI_Send(gdiaguur,ncount,MPI_DOUBLE_PRECISION,0,9999+iter+mpitag,MPI_Comm_Row_hw,ierr)
          call MPI_Send(gdiagddr,ncount,MPI_DOUBLE_PRECISION,0,8998+iter+mpitag,MPI_Comm_Row_hw,ierr)
          call MPI_Send(gdiagud,ncount,MPI_DOUBLE_COMPLEX,0,7997+iter+mpitag,MPI_Comm_Row_hw,ierr)
          call MPI_Send(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,0,6996+iter+mpitag,MPI_Comm_Row_hw,ierr)
          ! Receiving new point or signal to exit
          call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_Comm_Row_hw,stat,ierr)
          if(ix==0) exit
        end do
      end if

      ! Send results to all processors
      call MPI_Bcast(n_orb_u,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
      call MPI_Bcast(n_orb_d,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
      call MPI_Bcast(mp,Npl,MPI_DOUBLE_COMPLEX,0,MPI_Comm_Row_hw,ierr)

      n_orb_u = 0.5d0 + n_orb_u/pi
      n_orb_d = 0.5d0 + n_orb_d/pi
      n_orb_t = n_orb_u + n_orb_d
      mag_orb = n_orb_u - n_orb_d
      mp      = mp/pi
      mx      = real(mp)
      my      = aimag(mp)

      do i=1,Npl
        ! Number of particles
        n_t(i) = sum(n_orb_t(i,:))
        fvec(i)   = n_t(i) - npart0(i+offset)
        ! x-component of magnetization
        j = i+Npl
        fvec(j)  = mx(i) - mx_in(i)
        ! y-component of magnetization
        j = j+Npl
        fvec(j)  = my(i) - my_in(i)
        ! z-component of magnetization
        j = j+Npl
        mz(i)    = sum(mag_orb(i,5:9))
        fvec(j)  = mz(i) - mz_in(i)
      end do

      if(myrank_row_hw==0) then
        do i=1,Npl
          if(abs(mp(i))>1.d-10) then
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+Npl,fvec(i+Npl),i+2*Npl,fvec(i+2*Npl),i+3*Npl,fvec(i+3*Npl)
          else
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*Npl,fvec(i+3*Npl)
          end if
        end do
      end if
      if(lontheflysc) call write_sc_results()
    case(2)
      selfconjac = 0.d0
      ! Calculating the jacobian
      if (myrank_row_hw==0) then ! Process 0 receives all results and send new tasks if necessary
        write(outputunit_loop,"('|------------------- Iteration ',i4,' (jacobian) -------------------|')") iter
        call sumk_jacobian(Ef,y(ix),ggr)
        selfconjac = wght(ix)*ggr

        if(lverbose) write(outputunit_loop,"('[sc_eqs_and_jac_old2] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
        do i=2,pn1
          if(lverbose) call progress_bar(outputunit_loop,"jacobian energy points",i,pn1)

          call MPI_Recv(ggr,ncount2,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3333+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
          if(lverbose) write(outputunit_loop,"('[sc_eqs_and_jac_old2] Point ',i0,' received from ',i0)") i,stat(MPI_SOURCE)

          selfconjac = selfconjac + ggr

          ! If the number of processors is less than the total number of points, sends
          ! the rest of the points to the ones that finish first
          if (itask<pn1) then
            itask = itask + 1
            call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_Comm_Row_hw,ierr)
          else
            call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_Comm_Row_hw,ierr)
          end if
        end do
      else
        ! Other processors calculate each point of the integral and waits for new points
        do
          if(ix>pn1) exit

          ! First and second integrations (in the complex plane)
          call sumk_jacobian(Ef,y(ix),ggr)
          ggr = wght(ix)*ggr

  !         if(lverbose) write(outputunit_loop,"('[sc_equations_and_jacobian2] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
          ! Sending results to process 0
          call MPI_Send(ggr,ncount2,MPI_DOUBLE_PRECISION,0,3333+iter+mpitag,MPI_Comm_Row_hw,ierr)
          ! Receiving new point or signal to exit
          call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_Comm_Row_hw,stat,ierr)
          if(ix==0) exit
        end do
      end if

      ! Send results to all processors
      call MPI_Bcast(selfconjac,ncount2,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)

      selfconjac = selfconjac/pi
      do i = Npl+1,4*Npl
        selfconjac(i,i) = selfconjac(i,i) - 1.d0
      end do
    case default
      write(outputunit,"('[sc_eqs_and_jac_old] Problem in self-consistency! iflag = ',I0)") iflag
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select flag

    iter = iter + 1

    return
  end subroutine sc_eqs_and_jac_old

  ! For a given value of center of band eps1 it calculates the
  ! occupation number and the magnetic moment
  subroutine sc_eqs_old(N,x,fvec,iflag)
    use mod_constants
    use mod_parameters
    use mod_f90_kind
    use mod_generate_epoints
    use mod_magnet
    use mod_tight_binding, only: npart0
    use mod_progress
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag
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

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:Npl)
    mx_in = x(Npl+1:2*Npl)
    my_in = x(2*Npl+1:3*Npl)
    mz_in = x(3*Npl+1:4*Npl)
    mp_in = mx_in+zi*my_in
    do i=1,Npl
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,Npl
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    n_orb_u = 0.d0
    n_orb_d = 0.d0

    ix = myrank_row_hw+1
    itask = numprocs ! Number of tasks done initially

    ! Calculating the number of particles for each spin and orbital using a complex integral
    if (myrank_row_hw==0) then ! Process 0 receives all results and send new tasks if necessary
      write(outputunit_loop,"('|------------------- Iteration ',i4,' (densities) ------------------|')") iter
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

      if(lverbose) write(outputunit_loop,"('[sc_eqs_old] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
      do i=2,pn1
        if(lverbose) call progress_bar(outputunit_loop,"densities energy points",i,pn1)

        call MPI_Recv(gdiaguur,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,9999+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
        call MPI_Recv(gdiagddr,ncount,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),9998+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
        call MPI_Recv(gdiagud,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9997+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)
        call MPI_Recv(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),9996+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)

        n_orb_u = n_orb_u + gdiaguur
        n_orb_d = n_orb_d + gdiagddr

        do j=1,Npl
          mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
        end do

        ! If the number of processors is less than the total number of points, sends
        ! the rest of the points to the ones that finish first
        if (itask<pn1) then
          itask = itask + 1
          call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_Comm_Row_hw,ierr)
        else
          call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_Comm_Row_hw,ierr)
        end if
      end do
    else
      ! Other processors calculate each point of the integral and waits for new points
      do
        if(ix>pn1) exit

        ! First and second integrations (in the complex plane)
        call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
        gdiaguur = wght(ix)*gdiaguur
        gdiagddr = wght(ix)*gdiagddr
        gdiagud = wght(ix)*gdiagud
        gdiagdu = wght(ix)*gdiagdu

  !       if(lverbose) write(outputunit_loop,"('[sc_equations] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
        ! Sending results to process 0
        call MPI_Send(gdiaguur,ncount,MPI_DOUBLE_PRECISION,0,9999+iter+mpitag,MPI_Comm_Row_hw,ierr)
        call MPI_Send(gdiagddr,ncount,MPI_DOUBLE_PRECISION,0,9998+iter+mpitag,MPI_Comm_Row_hw,ierr)
        call MPI_Send(gdiagud,ncount,MPI_DOUBLE_COMPLEX,0,9997+iter+mpitag,MPI_Comm_Row_hw,ierr)
        call MPI_Send(gdiagdu,ncount,MPI_DOUBLE_COMPLEX,0,9996+iter+mpitag,MPI_Comm_Row_hw,ierr)
        ! Receiving new point or signal to exit
        call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_Comm_Row_hw,stat,ierr)
        if(ix==0) exit
      end do
    end if

    ! Send results to all processors
    call MPI_Bcast(n_orb_u,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
    call MPI_Bcast(n_orb_d,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
    call MPI_Bcast(mp,Npl,MPI_DOUBLE_COMPLEX,0,MPI_Comm_Row_hw,ierr)

    n_orb_u = 0.5d0 + n_orb_u/pi
    n_orb_d = 0.5d0 + n_orb_d/pi
    n_orb_t = n_orb_u + n_orb_d
    mag_orb = n_orb_u - n_orb_d
    mp      = mp/pi
    mx      = real(mp)
    my      = aimag(mp)

    do i=1,Npl
      ! Number of particles
      n_t(i) = sum(n_orb_t(i,:))
      fvec(i)   = n_t(i) - npart0(i+offset)
      ! x-component of magnetization
      j = i+Npl
      fvec(j)  = mx(i) - mx_in(i)
      ! y-component of magnetization
      j = j+Npl
      fvec(j)  = my(i) - my_in(i)
      ! z-component of magnetization
      j = j+Npl
      mz(i)    = sum(mag_orb(i,5:9))
      fvec(j)  = mz(i) - mz_in(i)
    end do

    if(myrank_row_hw==0) then
      do i=1,Npl
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+Npl,fvec(i+Npl),i+2*Npl,fvec(i+2*Npl),i+3*Npl,fvec(i+3*Npl)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*Npl,fvec(i+3*Npl)
        end if
      end do
    end if

    if(lontheflysc) call write_sc_results()

    iter = iter + 1

    return
  end subroutine sc_eqs_old

  subroutine sc_jac_old(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_constants
    use mod_parameters
    use mod_f90_kind
    use mod_generate_epoints
    use mod_magnet
    use mod_progress
    use mod_mpi_pars
    implicit none
    integer       :: N,ldfjac,i,iflag
    real(double)  :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(double),dimension(Npl)         :: mx_in,my_in,mz_in
    complex(double),dimension(Npl)      :: mp_in
    real(double),dimension(4*Npl,4*Npl) :: ggr
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount=16*Npl*Npl

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:Npl)
    mx_in = x(Npl+1:2*Npl)
    my_in = x(2*Npl+1:3*Npl)
    mz_in = x(3*Npl+1:4*Npl)
    mp_in = mx_in+zi*my_in
    do i=1,Npl
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    fvec=fvec

    ix = myrank_row_hw+1
    itask = numprocs ! Number of tasks done initially

    ! Calculating the jacobian
    if (myrank_row_hw==0) then ! Process 0 receives all results and send new tasks if necessary
      write(outputunit_loop,"('|------------------- Iteration ',i4,' (jacobian) -------------------|')") iter
      call sumk_jacobian(Ef,y(ix),ggr)
      selfconjac = wght(ix)*ggr

      if(lverbose) write(outputunit_loop,"('[sc_jac_old] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
      do i=2,pn1
        if(lverbose) call progress_bar(outputunit_loop,"jacobian energy points",i,pn1)
        ! Progress bar

        call MPI_Recv(ggr,ncount,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,3333+iter+mpitag,MPI_Comm_Row_hw,stat,ierr)

        selfconjac = selfconjac + ggr

        ! If the number of processors is less than the total number of points, sends
        ! the rest of the points to the ones that finish first
        if (itask<pn1) then
          itask = itask + 1
          call MPI_Send(itask,1,MPI_INTEGER,stat(MPI_SOURCE),itask,MPI_Comm_Row_hw,ierr)
        else
          call MPI_Send(0,1,MPI_INTEGER,stat(MPI_SOURCE),0,MPI_Comm_Row_hw,ierr)
        end if
      end do
    else
      ! Other processors calculate each point of the integral and waits for new points
      do
        if(ix>pn1) exit

        ! First and second integrations (in the complex plane)
        call sumk_jacobian(Ef,y(ix),ggr)
        ggr = wght(ix)*ggr

  !       if(lverbose) write(outputunit_loop,"('[sc_jac_old] Finished point ',i0,' in rank ',i0,' (',a,')')") ix,myrank_row_hw,trim(host)
        ! Sending results to process 0
        call MPI_Send(ggr,ncount,MPI_DOUBLE_PRECISION,0,3333+iter+mpitag,MPI_Comm_Row_hw,ierr)
        ! Receiving new point or signal to exit
        call MPI_Recv(ix,1,MPI_INTEGER,0,MPI_ANY_TAG,MPI_Comm_Row_hw,stat,ierr)
        if(ix==0) exit
      end do
    end if

    ! Send results to all processors
    call MPI_Bcast(selfconjac,ncount,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)

    selfconjac = selfconjac/pi
    do i = Npl+1,4*Npl
      selfconjac(i,i) = selfconjac(i,i) - 1.d0
    end do

    iter = iter + 1

    return
  end subroutine sc_jac_old

  ! Integration of Green functions over k values to calculate the number of particles
  subroutine sumk_npart(er,ei,gdiaguur,gdiagddr,gdiagud,gdiagdu)
    use mod_f90_kind
    use mod_constants
    use mod_parameters
    use mod_progress
    use mod_mpi_pars
    use mod_system, only: nkpt, kbz, wkbz
  !$  use omp_lib
    implicit none
  !$  integer     :: nthreads,mythread
    integer       :: iz,i,mu,mup
    real(double)  :: kp(3)
    real(double),intent(in)  :: er,ei
    real(double),dimension(Npl,9),intent(out)    :: gdiaguur,gdiagddr
    complex(double),dimension(Npl,9),intent(out) :: gdiagud,gdiagdu
    complex(double),dimension(Npl,Npl,18,18)     :: gf

    gdiaguur= 0.d0
    gdiagddr= 0.d0
    gdiagud = zero
    gdiagdu = zero

  !$omp parallel default(none) &
  !$omp& private(mythread,iz,kp,gf,i,mu,mup) &
  !$omp& shared(llineargfsoc,llinearsoc,lverbose,kbz,wkbz,nkpt,er,ei,gdiaguur,gdiagddr,gdiagud,gdiagdu,Npl,myrank_row_hw,nthreads,outputunit_loop)
  !$  mythread = omp_get_thread_num()
  !$  if((mythread==0).and.(myrank_row_hw==0)) then
  !$    nthreads = omp_get_num_threads()
  !$    write(outputunit_loop,"('[sumk_npart] Number of threads: ',i0)") nthreads
  !$  end if

  !$omp do
    kpoints: do iz=1,nkpt
  !$  if((mythread==0)) then
        if((myrank_row_hw==0).and.(lverbose)) call progress_bar(outputunit_loop,"densities kpoints",iz,nkpt)
  !$   end if

      kp = kbz(:,iz)

      ! Green function on energy Ef + iy, and wave vector kp
      if((llineargfsoc).or.(llinearsoc)) then
        call greenlineargfsoc(er,ei,kp,gf)
      else
        call green(er,ei,kp,gf)
      end if
      !$omp critical
      planes: do i=1,Npl
        orbital_index: do mu=1,9
          mup = mu+9
          gdiaguur(i,mu) = gdiaguur(i,mu) + real(gf(i,i,mu,mu)*wkbz(iz))
          gdiagddr(i,mu) = gdiagddr(i,mu) + real(gf(i,i,mup,mup)*wkbz(iz))
          gdiagud(i,mu) = gdiagud(i,mu) + (gf(i,i,mu,mup)*wkbz(iz))
          gdiagdu(i,mu) = gdiagdu(i,mu) + (gf(i,i,mup,mu)*wkbz(iz))
        end do orbital_index
      end do planes
      !$omp end critical
    end do kpoints
  !$omp end do
  !$omp end parallel
    return
  end subroutine sumk_npart

  ! Integration of Green functions over k vectors to calculate
  ! the jacobian of self-consistency system
  subroutine sumk_jacobian(er,ei,ggr)
    use mod_f90_kind
    use mod_constants
    use mod_parameters
    use mod_progress
    use mod_mpi_pars
    use mod_system, only: nkpt, kbz, wkbz
!$  use omp_lib
    implicit none
!$  integer                  :: nthreads,mythread
    integer                  :: AllocateStatus
    integer                  :: iz,i,j,i0,j0,mu,sigma,sigmap
    real(double)             :: kp(3)
    real(double),intent(in)  :: er,ei
    real(double),dimension(4*Npl,4*Npl),intent(out) :: ggr
    complex(double)                                 :: mhalfU(4,Npl),wkbzc
    complex(double),dimension(4,18,18)              :: pauli_components1,pauli_components2,temp1,temp2
    complex(double),dimension(18,18)                :: gij,gji,temp,paulitemp
    complex(double),dimension(:,:,:,:,:,:),allocatable    :: gdHdxg,gvgdHdxgvg
    complex(double),dimension(Npl,Npl,18,18)        :: gf,gvg

    pauli_components1 = zero
    pauli_components2 = zero
!   Includes d orbitals in the charge component
    pauli_components1(1,:,:)   = identorb18(:,:)
    pauli_components1(2:4,:,:) = pauli_dorb(:,:,:)
!   Excludes d orbitals in the charge component
    pauli_components2(1, 5: 9, 5: 9) = identorb18( 5: 9, 5: 9)
    pauli_components2(1,14:18,14:18) = identorb18(14:18,14:18)
    pauli_components2(2:4,:,:) = pauli_dorb(:,:,:)

  ! Prefactor -U/2 in dH/dm and 1 in dH/deps1
    do j=1,Npl
      mhalfU(1,j) = zum
      mhalfU(2:4,j) = -0.5d0*U(j+offset)
    end do

    ggr    = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,AllocateStatus,errorcode,ierr,iz,kp,wkbzc,gf,gvg,temp,temp1,temp2,gij,gji,paulitemp,gdHdxg,gvgdHdxgvg,i,j,i0,j0,mu,sigma) &
!$omp& shared(llineargfsoc,llinearsoc,lverbose,kbz,wkbz,nkpt,er,ei,ggr,mhalfU,pauli_components1,pauli_components2,Npl,myrank_row_hw,nthreads,outputunit,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if((mythread==0).and.(myrank_row_hw==0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[sumk_jacobian] Number of threads: ',i0)") nthreads
!$  end if
    allocate( gdHdxg(4,4,Npl,Npl,18,18),gvgdHdxgvg(4,4,Npl,Npl,18,18) , STAT = AllocateStatus  )
    if (AllocateStatus/=0) then
      write(outputunit,"('[sumk_jacobian] Not enough memory for: gdHdxg,gvgdHdxgvg')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if

!$omp do reduction(+:ggr)
    kpoints: do iz=1,nkpt
!$  if((mythread==0)) then
        if((myrank_row_hw==0).and.(lverbose)) call progress_bar(outputunit_loop,"jacobian kpoints",iz,nkpt)
!$   end if

      kp = kbz(:,iz)
      wkbzc = cmplx(wkbz(iz),0.d0)

      ! Green function on energy Ef + iy, and wave vector kp
      if((llineargfsoc).or.(llinearsoc)) then
        call greenlinearsoc(er,ei,kp,gf,gvg)
        gf = gf + gvg
      else
        call green(er,ei,kp,gf)
      end if

      do j=1,Npl ; do i=1,Npl
        gij = gf(i,j,:,:)
        gji = gf(j,i,:,:)

        do sigma = 1,4
          ! temp1 =  pauli*g_ij
          paulitemp = pauli_components1(sigma,:,:)
          call zgemm('n','n',18,18,18,zum,paulitemp,18,gij,18,zero,temp,18)
          temp1(sigma,:,:) = temp
        end do

        do sigmap = 1,4
          ! temp2 = (-U/2) * sigma* g_ji
          paulitemp = pauli_components2(sigmap,:,:)
          call zgemm('n','n',18,18,18,mhalfU(sigmap,j),paulitemp,18,gji,18,zero,temp,18)
          temp2(sigmap,:,:) = temp
        end do

        do sigmap = 1,4 ; do sigma = 1,4
          ! gdHdxg = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji
          gij = temp1(sigma,:,:)
          gji = temp2(sigmap,:,:)
          call zgemm('n','n',18,18,18,wkbzc,gij,18,gji,18,zero,temp,18)
          gdHdxg(sigma,sigmap,i,j,:,:) = temp
        end do ; end do


        if((llineargfsoc).or.(llinearsoc)) then ! non-linear term
          gij = gvg(i,j,:,:)
          gji = gvg(j,i,:,:)

          do sigma = 1,4
            ! temp1 = wkbz* pauli*gvg_ij
            paulitemp = pauli_components1(sigma,:,:)
            call zgemm('n','n',18,18,18,zum,paulitemp,18,gij,18,zero,temp,18)
            temp1(sigma,:,:) = temp
          end do

          do sigmap = 1,4
            ! temp2 = (-U/2) * sigma* gvg_ji
            paulitemp = pauli_components2(sigmap,:,:)
            call zgemm('n','n',18,18,18,mhalfU(sigmap,j),paulitemp,18,gji,18,zero,temp,18)
            temp2(sigmap,:,:) = temp
          end do

          do sigmap = 1,4 ; do sigma = 1,4
            ! gdHdxg = temp1*temp2 = wkbz* pauli*gvg_ij*(-U/2)*sigma* gvg_ji
            gij = temp1(sigma,:,:)
            gji = temp2(sigmap,:,:)
            call zgemm('n','n',18,18,18,wkbzc,gij,18,gji,18,zero,temp,18)
            gvgdHdxgvg(sigma,sigmap,i,j,:,:) = temp
          end do ; end do
        end if
      end do ; end do

      ! removing non-linear SOC term
      if((llineargfsoc).or.(llinearsoc)) gdHdxg = gdHdxg - gvgdHdxgvg

      do mu=1,18 ; do j=1,Npl ; do i=1,Npl ; do sigmap=1,4 ; do sigma=1,4
        i0 = (sigma-1)*Npl + i
        j0 = (sigmap-1)*Npl + j
        ! Trace over orbitals and spins of the real part
        ggr(i0,j0) = ggr(i0,j0) + real(gdHdxg(sigma,sigmap,i,j,mu,mu))
      end do ; end do ; end do ; end do ; end do
    end do kpoints
!$omp end do
    deallocate(gdHdxg,gvgdHdxgvg)
!$omp end parallel

    return
  end subroutine sumk_jacobian

end module mod_self_consistency
