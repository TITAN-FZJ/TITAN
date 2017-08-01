module mod_self_consistency
  implicit none
  character(len=300)  :: default_file

contains
  ! Tries to read eps1 and m if available - includes hdel, hdelp and hdelm calculations
  subroutine read_previous_results(lsuccess)
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_magnet, only: mabs, mx, my, mz, mtheta, mphi, &
                          mvec_cartesian, mvec_spherical, &
                          eps1, hw_count, hw_list, lfield, &
                          mp, hdel, hdelm, hdelp
    use mod_parameters, only: skipsc, outputunit_loop, lselfcon, U,&
                              magaxis, magaxisvec, offset, layertype
    use mod_system, only: s => sys
    use mod_mpi_pars, only: myrank_row_hw, myrank, abortProgram
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
      do i = 1,s%nAtoms
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
        continue
      else if(magaxis == -2) then
        magaxisvec = magaxisvec(1) * s%a1 + magaxisvec(2) * s%a2 + magaxisvec(3) * s%a3
      else if(magaxis == -3) then
        magaxisvec = [cos(magaxisvec(2)*pi/180)*sin(magaxisvec(1)*pi/180), sin(magaxisvec(2)*pi/180)*sin(magaxisvec(1)*pi/180), cos(magaxisvec(1)*pi/180)]
      else if(magaxis == 0) then
        magaxisvec = [0.d0, 0.d0, sign(1.0d0, hw_list(hw_count,1))]
      else if(magaxis >=1 .and. magaxis <= s%nAtoms) then
        !magaxisvec(1:3) = c_nn(1:3, magaxis)
        stop "Not Implemented"
      else
        if(myrank.eq.0) call abortProgram("[read_previous_results] Unknown magnetization direction!")
      end if
      magaxisvec = magaxisvec / sqrt(dot_product(magaxisvec, magaxisvec))
      magaxisvec = magaxisvec * 0.5d0

      mx = magaxisvec(1)
      my = magaxisvec(2)
      mz = magaxisvec(3)

      do i=1,s%nAtoms
        if(layertype(i+offset)==2) then
          mx(i) = mx(i) * sign(4.d0,hw_list(hw_count,1))
          my(i) = my(i) * sign(4.d0,hw_list(hw_count,1))
          mz(i) = mz(i) * sign(4.d0,hw_list(hw_count,1))
        endif
      end do

      if(lfield .and. magaxis == 0) then
        mx = mz*sin(hw_list(hw_count,2)*pi)*cos(hw_list(hw_count,3)*pi)
        my = mz*sin(hw_list(hw_count,2)*pi)*sin(hw_list(hw_count,3)*pi)
        mz = mz*cos(hw_list(hw_count,2)*pi)
      end if
      mp = cmplx(mx,my,double)

      ! Variables used in the hamiltonian
      do i=1,s%nAtoms
        hdel(i)   = 0.5d0*U(i+offset)*mz(i)
        hdelp(i)  = 0.5d0*U(i+offset)*mp(i)
      end do
      hdelm = conjg(hdelp)
    end if

    return
  end subroutine read_previous_results

  ! Rotate the magnetization to the direction of the field (useful for SOC=F)
  subroutine rotate_magnetization_to_field()
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_magnet, only: hw_count, hw_list, hhwx, hhwy, hhwz, &
                          mx, my, mz, mabs, mp

    use mod_parameters, only: outputunit_loop
    use mod_System, only: s => sys
    use mod_mpi_pars, only: myrank_row_hw
    implicit none
    integer :: i,err,sign
    logical :: lsuccess
    real(double) :: mdotb

    if(myrank_row_hw==0) write(outputunit_loop,"('[rotate_magnetization_to_field] Rotating previous magnetization to the direction of the field...')")

    do i=1,s%nAtoms
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

    !TODO: Check if it's a good idea to dealloc
    !deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1)
    !deallocate(mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm)
    !deallocate(mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi)
    !if(lGSL) deallocate(lxm,lym,lzm,lxpm,lypm,lzpm)
    !deallocate(mmlayer,layertype,U,mmlayermag,lambda,npart0)
    !deallocate(t0, t0i)

    return
  end subroutine rotate_magnetization_to_field

  ! This subroutine performs the self-consistency
  subroutine do_self_consistency()
    use mod_f90_kind, only: double
    use mod_constants, only: pi
    use mod_parameters, only: outputunit_loop, lslatec, lnojac, tol
    use mod_magnet, only: eps1, mx, my, mz, mabs, mtheta, mphi, mvec_cartesian, mvec_spherical, &
                          hw_count, iter
    use mod_mpi_pars, only: myrank_row_hw,mpitag
    use mod_system, only: s => sys
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

    neq = 4*s%nAtoms
    allocate( sc_solu(neq),diag(neq),qtf(neq),fvec(neq),jac(neq,neq) )

    ! Putting read eps1 existing solutions into esp1_solu (first guess of the subroutine)
    sc_solu(1:s%nAtoms)         = eps1
    sc_solu(s%nAtoms+1:2*s%nAtoms)   = mx
    sc_solu(2*s%nAtoms+1:3*s%nAtoms) = my
    sc_solu(3*s%nAtoms+1:4*s%nAtoms) = mz
    iter  = 1
    !mpitag = (Npl-Npl_i)*total_hw_npt1 + hw_count
    mpitag = hw_count
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
    do i = 1,s%nAtoms
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
    use mod_parameters, only: scfile
    use mod_magnet, only: hw_count,total_hw_npt1
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
    use mod_parameters, only: outputunit_loop,lGSL
    use mod_system, only: s => sys
    !use mod_mpi_pars
    use mod_magnet, only: eps1, mx, my, mz, mp, mphi, mtheta, mabs, &
                          lxpm, lypm, lzpm, lpphi, lptheta, lxm, lym, lzm, lpabs, labs
    implicit none
    integer :: i

    write(outputunit_loop,"('|----------=============== Self-consistent ground state: ===============----------|')")
    write(outputunit_loop,"(11x,' *************** Center of d bands: ***************')")
    do i=1,s%nAtoms
      write(outputunit_loop,"(26x,'eps1(',i2.0,')=',f11.8)") i,eps1(i)
    end do
    write(outputunit_loop,"(11x,' *********** Magnetization components: **********')")
    do i=1,s%nAtoms
      write(outputunit_loop,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mx(i),i,my(i),i,mz(i)
      if(abs(mp(i))/=0) write(outputunit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") mtheta(i),mphi(i)
    end do
    if(lGSL) then
      write(outputunit_loop,"(11x,' *** Orbital components in spin coordinates:  ***')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(4x,'Lxp(',i2.0,')=',f11.8,4x,'Lyp(',i2.0,')=',f11.8,4x,'Lzp(',i2.0,')=',f11.8)") i,lxpm(i),i,lypm(i),i,lzpm(i)
        if(sqrt(lxpm(i)**2+lypm(i)**2)/=0) write(outputunit_loop,"(12x,'theta =',f11.8,' pi',4x,'phi =',f11.8,' pi')") lptheta(i),lpphi(i)
      end do
      write(outputunit_loop,"(11x,' *** Orbital components in cubic coordinates: ***')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(4x,'Lx (',i2.0,')=',f11.8,4x,'Ly (',i2.0,')=',f11.8,4x,'Lz (',i2.0,')=',f11.8)") i,lxm(i),i,lym(i),i,lzm(i)
      end do
      write(outputunit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(4x,'M (',i2.0,') =',f11.8,4x,'Lp (',i2.0,')=',f11.8,4x,'L (',i2.0,') =',f11.8)") i,mabs(i),i,lpabs(i),i,labs(i)
      end do
    else
      write(outputunit_loop,"(11x,' ******************** Total: ********************')")
      do i=1,s%nAtoms
        write(outputunit_loop,"(27x,'M (',i2.0,') =',f11.8)") i,mabs(i)
      end do
    end if
    write(outputunit_loop,"('|----------=============================================================----------|')")

    return
  end subroutine write_sc_results_on_screen

  ! This subroutine reads previous band-shifting and magnetization results
  ! and also writes new ones into file
  subroutine read_write_sc_results(iflag,err,lsuccess)
    use mod_f90_kind, only: double
    use mod_constants, only: zi
    use mod_parameters, only: offset, fieldpart, eta, U,Utype,scfile, outputunit_loop, Npl_folder, dfttype
    use EnergyIntegration, only: parts
    use mod_magnet, only: eps1, hdel, hdelm, hdelp, mp, mz, hw_count, mx, my, mz
    use mod_SOC, only: SOCc, socpart
    use mod_mpi_pars
    use mod_system, only: s => sys
    implicit none
    character(len=300)  :: file = ""
    character(len=100)  :: folder,prefix
    integer,intent(in)  :: iflag
    integer,intent(out) :: err
    logical,intent(out) :: lsuccess
    integer             :: i
    real(double)        :: previous_results(s%nAtoms,4)
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

    lsuccess = .false.
  !   Reading previous results (mx, my, mz and eps1) from files (if available)
    if(iflag==0) then
      if(trim(scfile)=="") then ! If a filename is not given in inputcard (or don't exist), use the default one
        write(file,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts,Utype,trim(fieldpart),s%nkpt,eta,trim(socpart)
        open(unit=99,file=file,status="old",iostat=err)
        if((err==0).and.(myrank_row_hw==0)) then
          write(outputunit_loop,"('[read_write_sc_results] Self-consistency file already exists. Reading it now...')")
          write(outputunit_loop,"(a)") trim(file)
        else
          default_file = trim(file)
        end if
      else ! If filename in inputcard exists or 2nd+ angular iteration
        if(((hw_count)==1)) then !.and.(Npl==Npl_i)) then ! Filename in inputcard (1st iteration on loop)
          open(unit=99,file=scfile,status="old",iostat=err)
          if((err==0).and.(myrank_row_hw==0)) then
            write(outputunit_loop,"('[read_write_sc_results] Using filename given in input file for self-consistency:')")
            write(outputunit_loop,"(a)") trim(scfile)
          end if
        else ! 2nd+ iteration, cheking if default file exists
          write(file,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts,Utype,trim(fieldpart),s%nkpt,eta,trim(socpart)
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
          do i=1,s%nAtoms
            read(99,fmt=*) previous_results(i,1)
            read(99,fmt=*) previous_results(i,2)
            read(99,fmt=*) previous_results(i,3)
            read(99,fmt=*) previous_results(i,4)
          end do
        end if
        call MPI_Bcast(previous_results,4*s%nAtoms,MPI_DOUBLE_PRECISION,0,MPI_Comm_Row_hw,ierr)
        eps1(:) = previous_results(:,1)
        mx  (:) = previous_results(:,2)
        my  (:) = previous_results(:,3)
        mz  (:) = previous_results(:,4)
        mp  = mx + zi*my
        do i=1,s%nAtoms
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
        write(file,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts-1,Utype,trim(fieldpart),s%nkpt,eta,trim(socpart)
        open(unit=99,file=file,status="old",iostat=err)
        if(err==0) then
          if(myrank_row_hw==0) then
            write(outputunit_loop,"('[read_write_sc_results] Self-consistency file does not exist. Reading results for parts-1 now...')")
            write(outputunit_loop,"('[read_write_sc_results] Updating values obtained for parts-1...')")
            write(outputunit_loop,"(a)") file
          end if
          do i=1,s%nAtoms
            read(99,*) eps1(i)
            read(99,*) mx(i)
            read(99,*) my(i)
            read(99,*) mz(i)
          end do
          mp  = mx + zi*my
          do i=1,s%nAtoms
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
      write(scfile,"('./results/',a1,'SOC/selfconsistency/',a,a,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_nkpt=',i0,'_eta=',es8.1,a,'.dat')") SOCc, trim(prefix),trim(Npl_folder),dfttype,parts,Utype,trim(fieldpart),s%nkpt,eta,trim(socpart)
      open (unit=99,status='replace',file=scfile)
      do i=1,s%nAtoms
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
    use mod_f90_kind, only: double
    use mod_constants, only: zi, pi, zero
    use mod_parameters, only: offset, U, outputunit, outputunit_loop, Ef, lverbose, host, lontheflysc
    use EnergyIntegration, only: pn1, y, wght
    use mod_system, only: s => sys
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx, my, mz
    use mod_progress, only: progress_bar
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag
    integer,           intent(inout) :: iuser(*)
    real(double),      intent(inout) :: ruser(*)
    real(double),dimension(N) :: x,fvec
    real(double),dimension(N,N)      :: selfconjac
    real(double),dimension(s%nAtoms)      :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)   :: mp_in
    real(double),dimension(N,N)      :: ggr
    real(double),dimension(s%nAtoms,9)    :: n_orb_u,n_orb_d,n_orb_t,mag_orb
    real(double),dimension(s%nAtoms,9)    :: gdiaguur,gdiagddr
    complex(double),dimension(s%nAtoms,9) :: gdiagud,gdiagdu
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount,ncount2
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount=s%nAtoms*9
    ncount2=N*N

  ! Values used in the hamiltonian
    eps1  = x(           1:  s%nAtoms)
    mx_in = x(  s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in + zi*my_in

    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx_in(i),i,my_in(i),i,mz_in(i)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz_in(i)
        end if
      end do
    end if

    ix = myrank_row_hw+1
    itask = numprocs ! Number of tasks done initially
    select case (iflag)
    case(1)
      n_orb_u = 0.d0
      n_orb_d = 0.d0
      mp = zero

      do while(ix <= pn1)
        call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)

        gdiaguur = wght(ix)*gdiaguur
        gdiagddr = wght(ix)*gdiagddr
        gdiagud = wght(ix)*gdiagud
        gdiagdu = wght(ix)*gdiagdu

        n_orb_u = n_orb_u + gdiaguur
        n_orb_d = n_orb_d + gdiagddr

        do j=1,s%nAtoms
          mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
        end do

        ix = ix + numprocs_row
      end do

      call MPI_Allreduce(MPI_IN_PLACE, n_orb_u, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, n_orb_d, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, mp, s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)

      n_orb_u = 0.5d0 + n_orb_u/pi
      n_orb_d = 0.5d0 + n_orb_d/pi
      n_orb_t = n_orb_u + n_orb_d
      mag_orb = n_orb_u - n_orb_d
      mp      = mp/pi
      mx      = real(mp)
      my      = aimag(mp)

      do i=1,s%nAtoms
        ! Number of particles
        n_t(i) = sum(n_orb_t(i,:))
        fvec(i) = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation !  npart0(i+offset)
        ! x-component of magnetization
        j = i + s%nAtoms
        fvec(j) = mx(i) - mx_in(i)
        ! y-component of magnetization
        j = j + s%nAtoms
        fvec(j) = my(i) - my_in(i)
        ! z-component of magnetization
        j = j+s%nAtoms
        mz(i)    = sum(mag_orb(i,5:9))
        fvec(j)  = mz(i) - mz_in(i)
      end do

      if(myrank_row_hw==0) then
        do i=1,s%nAtoms
          if(abs(mp(i))>1.d-10) then
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          else
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          end if
        end do
      end if
      if(lontheflysc) call write_sc_results()
    case(2)
      ! if((myrank_row_hw==0)) then
      !   write(outputunit_loop,*) "[sumk_jacobian]"
      ! end if

      selfconjac = 0.d0
      do while(ix <= pn1)
        call sumk_jacobian(Ef, y(ix), ggr)
        selfconjac = selfconjac + wght(ix)*ggr
        ix = ix + numprocs_row
      end do

      call MPI_Allreduce(MPI_IN_PLACE, selfconjac, ncount2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
      selfconjac = selfconjac/pi
      do i = s%nAtoms+1, 4*s%nAtoms
        selfconjac(i,i) = selfconjac(i,i) - 1.d0
      end do

    case default
      write(outputunit,"('[sc_equations_and_jacobian] Problem in self-consistency! iflag = ',I0)") iflag
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end select

    iter = iter + 1

    return
  end subroutine sc_equations_and_jacobian

  ! For a given value of center of band eps1 it calculates the
  ! occupation number and the magnetic moment
  subroutine sc_equations(N,x,fvec,iuser,ruser,iflag)
    use mod_constants, only: zi, pi, zero
    use mod_parameters, only: offset, U, outputunit_loop, Ef, lverbose, host, lontheflysc
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use EnergyIntegration, only: pn1, y, wght
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx ,my ,mz
    use mod_progress, only: progress_bar
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag
    integer     , intent(inout)         :: iuser(*)
    real(double), intent(inout)         :: ruser(*)
    real(double),dimension(N)           :: x,fvec
    real(double),dimension(s%nAtoms,9)       :: n_orb_u,n_orb_d,n_orb_t,mag_orb
    real(double),dimension(s%nAtoms)         :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)      :: mp_in
    real(double),dimension(s%nAtoms,9)       :: gdiaguur,gdiagddr
    complex(double),dimension(s%nAtoms,9)    :: gdiagud,gdiagdu
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount
    ncount=s%nAtoms*9
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+zi*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
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
    n_orb_u = 0.d0
    n_orb_d = 0.d0
    mp = zero

    do while(ix <= pn1)
      call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
      gdiaguur = wght(ix)*gdiaguur
      gdiagddr = wght(ix)*gdiagddr
      gdiagud = wght(ix)*gdiagud
      gdiagdu = wght(ix)*gdiagdu

      n_orb_u = n_orb_u + gdiaguur
      n_orb_d = n_orb_d + gdiagddr

      do j=1,s%nAtoms
        mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
      end do

      ix = ix + numprocs_row
    end do

    call MPI_Allreduce(MPI_IN_PLACE, n_orb_u, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, n_orb_d, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp, s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)

    n_orb_u = 0.5d0 + n_orb_u/pi
    n_orb_d = 0.5d0 + n_orb_d/pi
    n_orb_t = n_orb_u + n_orb_d
    mag_orb = n_orb_u - n_orb_d
    mp      = mp/pi
    mx      = real(mp)
    my      = aimag(mp)

    do i=1,s%nAtoms
      ! Number of particles
      n_t(i) = sum(n_orb_t(i,:))
      fvec(i)   = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation !npart0(i+offset)
      ! x-component of magnetization
      j = i+s%nAtoms
      fvec(j)  = mx(i) - mx_in(i)
      ! y-component of magnetization
      j = j+s%nAtoms
      fvec(j)  = my(i) - my_in(i)
      ! z-component of magnetization
      j = j+s%nAtoms
      mz(i)    = sum(mag_orb(i,5:9))
      fvec(j)  = mz(i) - mz_in(i)
    end do

    if(myrank_row_hw==0) then
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
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
    use mod_constants, only: zi, pi, zero
    use mod_parameters, only: offset, U, outputunit_loop, outputunit, Ef, lverbose, host, lontheflysc
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use EnergyIntegration, only: pn1, y, wght
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx, my, mz
    use mod_progress, only: progress_bar
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag,ldfjac
    real(double),dimension(N)        :: x,fvec
    real(double),dimension(ldfjac,N) :: selfconjac
    real(double),dimension(s%nAtoms)      :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)   :: mp_in
    real(double),dimension(N,N)      :: ggr
    real(double),dimension(s%nAtoms,9)    :: n_orb_u,n_orb_d,n_orb_t,mag_orb
    real(double),dimension(s%nAtoms,9)    :: gdiaguur,gdiagddr
    complex(double),dimension(s%nAtoms,9) :: gdiagud,gdiagdu
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount,ncount2
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount=s%nAtoms*9
    ncount2=N*N

  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+zi*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
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
      n_orb_u = 0.d0
      n_orb_d = 0.d0
      mp = zero

      do while(ix <= pn1)
        call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
        gdiaguur = wght(ix)*gdiaguur
        gdiagddr = wght(ix)*gdiagddr
        gdiagud = wght(ix)*gdiagud
        gdiagdu = wght(ix)*gdiagdu

        n_orb_u = n_orb_u + gdiaguur
        n_orb_d = n_orb_d + gdiagddr

        do j=1,s%nAtoms
          mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
        end do

        ix = ix + numprocs_row
      end do

      call MPI_Allreduce(MPI_IN_PLACE, n_orb_u, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, n_orb_d, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, mp, s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)

      n_orb_u = 0.5d0 + n_orb_u/pi
      n_orb_d = 0.5d0 + n_orb_d/pi
      n_orb_t = n_orb_u + n_orb_d
      mag_orb = n_orb_u - n_orb_d
      mp      = mp/pi
      mx      = real(mp)
      my      = aimag(mp)

      do i=1,s%nAtoms
        ! Number of particles
        n_t(i) = sum(n_orb_t(i,:))
        fvec(i)   = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation !npart0(i+offset)
        ! x-component of magnetization
        j = i+s%nAtoms
        fvec(j)  = mx(i) - mx_in(i)
        ! y-component of magnetization
        j = j+s%nAtoms
        fvec(j)  = my(i) - my_in(i)
        ! z-component of magnetization
        j = j+s%nAtoms
        mz(i)    = sum(mag_orb(i,5:9))
        fvec(j)  = mz(i) - mz_in(i)
      end do

      if(myrank_row_hw==0) then
        do i=1,s%nAtoms
          if(abs(mp(i))>1.d-10) then
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          else
            write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
            write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
          end if
        end do
      end if
      if(lontheflysc) call write_sc_results()
    case(2)
      selfconjac = 0.d0
      ! Calculating the jacobian
      do while(ix <= pn1)
        call sumk_jacobian(Ef,y(ix),ggr)
        selfconjac = selfconjac + wght(ix)*ggr

        ix = ix + numprocs_row
      end do
      call MPI_Allreduce(MPI_IN_PLACE, selfconjac, ncount2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)

      selfconjac = selfconjac/pi
      do i = s%nAtoms+1,4*s%nAtoms
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
    use mod_f90_kind, only: double
    use mod_constants, only: zi, pi, zero
    use mod_parameters, only: offset, U, outputunit_loop, Ef, lverbose, host, lontheflysc
    use mod_system, only: s => sys
    use EnergyIntegration, only: pn1, y, wght
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp, mp, mx, my, mz
    use mod_progress, only: progress_bar
    use mod_mpi_pars
    implicit none
    integer  :: N,i,j,iflag
    real(double),dimension(N)           :: x,fvec
    real(double),dimension(s%nAtoms,9)       :: n_orb_u,n_orb_d,n_orb_t,mag_orb
    real(double),dimension(s%nAtoms)         :: n_t,mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)      :: mp_in
    real(double),dimension(s%nAtoms,9)       :: gdiaguur,gdiagddr
    complex(double),dimension(s%nAtoms,9)    :: gdiagud,gdiagdu
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount
    ncount=s%nAtoms*9
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+zi*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    if((myrank_row_hw==0).and.(iter==1)) then
      write(outputunit_loop,"('|---------------- Starting eps1 and magnetization ----------------|')")
      do i=1,s%nAtoms
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
    n_orb_u = 0.d0
    n_orb_d = 0.d0
    mp = zero

    do while(ix <= pn1)
      call sumk_npart(Ef,y(ix),gdiaguur,gdiagddr,gdiagud,gdiagdu)
      gdiaguur = wght(ix)*gdiaguur
      gdiagddr = wght(ix)*gdiagddr
      gdiagud = wght(ix)*gdiagud
      gdiagdu = wght(ix)*gdiagdu

      n_orb_u = n_orb_u + gdiaguur
      n_orb_d = n_orb_d + gdiagddr

      do j=1,s%nAtoms
        mp(j) = mp(j) + (sum(gdiagdu(j,5:9)) + sum(conjg(gdiagud(j,5:9))))
      end do

      ix = ix + numprocs_row
    end do

    call MPI_Allreduce(MPI_IN_PLACE, n_orb_u, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, n_orb_d, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, mp, s%nAtoms, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)

    n_orb_u = 0.5d0 + n_orb_u/pi
    n_orb_d = 0.5d0 + n_orb_d/pi
    n_orb_t = n_orb_u + n_orb_d
    mag_orb = n_orb_u - n_orb_d
    mp      = mp/pi
    mx      = real(mp)
    my      = aimag(mp)

    do i=1,s%nAtoms
      ! Number of particles
      n_t(i) = sum(n_orb_t(i,:))
      fvec(i)   = n_t(i) - s%Types(s%Basis(i)%Material)%Occupation !npart0(i+offset)
      ! x-component of magnetization
      j = i+s%nAtoms
      fvec(j)  = mx(i) - mx_in(i)
      ! y-component of magnetization
      j = j+s%nAtoms
      fvec(j)  = my(i) - my_in(i)
      ! z-component of magnetization
      j = j+s%nAtoms
      mz(i)    = sum(mag_orb(i,5:9))
      fvec(j)  = mz(i) - mz_in(i)
    end do

    if(myrank_row_hw==0) then
      do i=1,s%nAtoms
        if(abs(mp(i))>1.d-10) then
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mx(',I2,')=',es16.9,4x,'My(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mx(i),i,my(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+s%nAtoms,fvec(i+s%nAtoms),i+2*s%nAtoms,fvec(i+2*s%nAtoms),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
        else
          write(outputunit_loop,"('Plane ',I2,': eps1(',I2,')=',es16.9,4x,'Mz(',I2,')=',es16.9)") i,i,eps1(i),i,mz(i)
          write(outputunit_loop,"(10x,'fvec(',I2,')=',es16.9,2x,'fvec(',I2,')=',es16.9)") i,fvec(i),i+3*s%nAtoms,fvec(i+3*s%nAtoms)
        end if
      end do
    end if

    if(lontheflysc) call write_sc_results()

    iter = iter + 1

    return
  end subroutine sc_eqs_old

  subroutine sc_jac_old(N,x,fvec,selfconjac,ldfjac,iflag)
    use mod_f90_kind, only: double
    use mod_constants, only: zi, pi
    use mod_parameters, only: offset, U, outputunit_loop, Ef, lverbose, host
    use mod_system, only: s => sys
    use EnergyIntegration, only: y, wght, pn1
    use mod_magnet, only: iter, eps1, hdel, hdelm, hdelp
    use mod_progress, only: progress_bar
    use mod_mpi_pars
    implicit none
    integer       :: N,ldfjac,i,iflag
    real(double)  :: x(N),fvec(N),selfconjac(ldfjac,N)
    real(double),dimension(s%nAtoms)         :: mx_in,my_in,mz_in
    complex(double),dimension(s%nAtoms)      :: mp_in
    real(double),dimension(4*s%nAtoms,4*s%nAtoms) :: ggr
    !--------------------- begin MPI vars --------------------
    integer :: ix,itask
    integer :: ncount
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^
    ncount=16*s%nAtoms*s%nAtoms

    iflag=0
  ! Values used in the hamiltonian
    eps1  = x(1:s%nAtoms)
    mx_in = x(s%nAtoms+1:2*s%nAtoms)
    my_in = x(2*s%nAtoms+1:3*s%nAtoms)
    mz_in = x(3*s%nAtoms+1:4*s%nAtoms)
    mp_in = mx_in+zi*my_in
    do i=1,s%nAtoms
      hdel(i)   = 0.5d0*U(i+offset)*mz_in(i)
      hdelp(i)  = 0.5d0*U(i+offset)*mp_in(i)
    end do
    hdelm = conjg(hdelp)

    fvec=fvec

    ix = myrank_row_hw+1
    itask = numprocs ! Number of tasks done initially

    selfconjac = 0.d0
    ! Calculating the jacobian
    do while(ix <= pn1)
      call sumk_jacobian(Ef,y(ix),ggr)
      selfconjac = selfconjac + wght(ix)*ggr

      ix = ix + numprocs
    end do
    call MPI_Allreduce(MPI_IN_PLACE, selfconjac, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_Comm_Row_hw, ierr)

    selfconjac = selfconjac/pi
    do i = s%nAtoms+1,4*s%nAtoms
      selfconjac(i,i) = selfconjac(i,i) - 1.d0
    end do

    iter = iter + 1

    return
  end subroutine sc_jac_old

  ! Integration of Green functions over k values to calculate the number of particles
  subroutine sumk_npart(er,ei,gdiaguur,gdiagddr,gdiagud,gdiagdu)
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use mod_constants, only: zero
    use mod_parameters, only: lverbose, outputunit_loop
    use mod_SOC, only: llineargfsoc, llinearsoc
    use mod_progress, only: progress_bar
    use mod_mpi_pars, only: myrank_row_hw, myrank
    use TightBinding, only: nOrb
    !use mod_system, only: nkpt, kbz, wkbz
    !$  use omp_lib
    implicit none
    !$  integer     :: nthreads,mythread
    integer       :: iz,i,mu,mup
    real(double)  :: kp(3)
    real(double),intent(in)  :: er,ei
    real(double),dimension(s%nAtoms,nOrb),intent(out)    :: gdiaguur,gdiagddr
    complex(double),dimension(s%nAtoms,nOrb),intent(out) :: gdiagud,gdiagdu
    complex(double),dimension(s%nAtoms,s%nAtoms,2*nOrb, 2*nOrb)     :: gf
    complex(double),dimension(:,:,:,:),allocatable :: gf_loc

    gdiaguur= 0.d0
    gdiagddr= 0.d0
    gdiagud = zero
    gdiagdu = zero

    !$omp parallel default(none) &
    !$omp& private(mythread,iz,kp,gf,i,mu,mup,gf_loc) &
    !$omp& shared(llineargfsoc,llinearsoc,lverbose,s,er,ei,gdiaguur,gdiagddr,gdiagud,gdiagdu,myrank_row_hw,nthreads,outputunit_loop)
    !$  mythread = omp_get_thread_num()
    !$  if((mythread==0).and.(myrank_row_hw==0)) then
    !$    nthreads = omp_get_num_threads()
    !$    write(outputunit_loop,"('[sumk_npart] Number of threads: ',i0)") nthreads
    !$  end if

    allocate(gf_loc(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb))
    gf_loc = zero

    !$omp do schedule(static)
    do iz=1,s%nkpt
    !!$  if((mythread==0)) then
    !!      if((myrank_row_hw==0).and.(lverbose)) call progress_bar(outputunit_loop,"densities kpoints",iz,s%nkpt)
    !!$   end if

      kp = s%kbz(:,iz)

      ! Green function on energy Ef + iy, and wave vector kp
      if((llineargfsoc).or.(llinearsoc)) then
        call greenlineargfsoc(er,ei,kp,gf)
      else
        call green(er,ei,kp,gf)
      end if
      gf_loc = gf_loc + gf*s%wkbz(iz)
    end do
    !$omp end do

    !$omp critical
    do i=1,s%nAtoms
      do mu=1,nOrb
        mup = mu+nOrb
        gdiaguur(i,mu) = gdiaguur(i,mu) + real(gf_loc(i,i,mu,mu))
        gdiagddr(i,mu) = gdiagddr(i,mu) + real(gf_loc(i,i,mup,mup))
        gdiagud(i,mu) = gdiagud(i,mu) + gf_loc(i,i,mu,mup)
        gdiagdu(i,mu) = gdiagdu(i,mu) + gf_loc(i,i,mup,mu)
      end do
    end do
    !$omp end critical
    deallocate(gf_loc)
    !$omp end parallel
    return
  end subroutine sumk_npart

  ! Integration of Green functions over k vectors to calculate
  ! the jacobian of self-consistency system
  subroutine sumk_jacobian(er,ei,ggr)
    use mod_f90_kind, only: double
    use mod_system, only: s => sys
    use mod_constants, only: identorb18, zero, pauli_dorb, zum
    use mod_parameters, only: offset, U, outputunit_loop, outputunit, lverbose
    use mod_SOC, only: llinearsoc, llineargfsoc
    use TightBinding, only: nOrb
    use mod_progress, only: progress_bar
    use mod_mpi_pars, only: abortProgram, myrank_row_hw
    !$  use omp_lib
    implicit none
    !$  integer                  :: nthreads,mythread
    integer                  :: AllocateStatus
    integer                  :: iz,i,j,i0,j0,mu,sigma,sigmap
    real(double)             :: kp(3)
    real(double),intent(in)  :: er,ei
    real(double),dimension(4*s%nAtoms,4*s%nAtoms),intent(out) :: ggr
    complex(double) :: mhalfU(4,s%nAtoms),wkbzc
    complex(double), dimension(4,2*nOrb, 2*nOrb) :: pauli_components1,pauli_components2,temp1,temp2
    complex(double), dimension(2*nOrb, 2*nOrb) :: gij,gji,temp,paulitemp
    complex(double), dimension(:,:,:,:,:,:), allocatable :: gdHdxg,gvgdHdxgvg
    complex(double), dimension(s%nAtoms,s%nAtoms,2*nOrb, 2*nOrb) :: gf,gvg

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
    do j=1,s%nAtoms
      mhalfU(1,j) = zum
      mhalfU(2:4,j) = -0.5d0*U(j+offset)
    end do

    ggr    = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,AllocateStatus,iz,kp,wkbzc,gf,gvg,temp,temp1,temp2,gij,gji,paulitemp,gdHdxg,gvgdHdxgvg,i,j,i0,j0,mu,sigma) &
!$omp& shared(llineargfsoc,llinearsoc,lverbose,s,er,ei,ggr,mhalfU,pauli_components1,pauli_components2,myrank_row_hw,nthreads,outputunit,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if((mythread==0).and.(myrank_row_hw==0)) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[sumk_jacobian] Number of threads: ',i0)") nthreads
!$  end if
    allocate( gdHdxg(4,4,s%nAtoms,s%nAtoms,18,18),gvgdHdxgvg(4,4,s%nAtoms,s%nAtoms,18,18) , STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumk_jacobian] Not enough memory for: gdHdxg,gvgdHdxgvg")

!$omp do schedule(static), reduction(+:ggr)
    kpoints: do iz=1,s%nkpt
!$  if((mythread==0)) then
        if((myrank_row_hw==0).and.(lverbose)) call progress_bar(outputunit_loop,"jacobian kpoints",iz,s%nkpt)
!$   end if

      kp = s%kbz(:,iz)
      wkbzc = cmplx(s%wkbz(iz), 0.d0)

      ! Green function on energy Ef + iy, and wave vector kp
      if((llineargfsoc).or.(llinearsoc)) then
        call greenlinearsoc(er,ei,kp,gf,gvg)
        gf = gf + gvg
      else
        call green(er,ei,kp,gf)
      end if

      do j=1,s%nAtoms
        do i=1,s%nAtoms
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

          do sigmap = 1,4
            do sigma = 1,4
              ! gdHdxg = temp1*temp2 = wkbz* pauli*g_ij*(-U/2)*sigma* g_ji
              gij = temp1(sigma,:,:)
              gji = temp2(sigmap,:,:)
              call zgemm('n','n',18,18,18,wkbzc,gij,18,gji,18,zero,temp,18)
              gdHdxg(sigma,sigmap,i,j,:,:) = temp
            end do
          end do


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

            do sigmap = 1,4
              do sigma = 1,4
                ! gdHdxg = temp1*temp2 = wkbz* pauli*gvg_ij*(-U/2)*sigma* gvg_ji
                gij = temp1(sigma,:,:)
                gji = temp2(sigmap,:,:)
                call zgemm('n','n',18,18,18,wkbzc,gij,18,gji,18,zero,temp,18)
                gvgdHdxgvg(sigma,sigmap,i,j,:,:) = temp
              end do
            end do
          end if
        end do
      end do

      ! removing non-linear SOC term
      if((llineargfsoc).or.(llinearsoc)) gdHdxg = gdHdxg - gvgdHdxgvg

      do mu=1,18
        do j=1,s%nAtoms
          do i=1,s%nAtoms
            do sigmap=1,4
              do sigma=1,4
                i0 = (sigma-1)*s%nAtoms + i
                j0 = (sigmap-1)*s%nAtoms + j
                ! Trace over orbitals and spins of the real part
                ggr(i0,j0) = ggr(i0,j0) + real(gdHdxg(sigma,sigmap,i,j,mu,mu))
              end do
            end do
          end do
        end do
      end do
    end do kpoints
!$omp end do
    deallocate(gdHdxg,gvgdHdxgvg)
!$omp end parallel

    return
  end subroutine sumk_jacobian

end module mod_self_consistency
