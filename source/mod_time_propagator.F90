module mod_time_propagator
  implicit none
contains

  subroutine time_propagator(s)
    use mod_kind,               only: dp,int32,int64
    use mod_constants,          only: cZero,cI
    use mod_imRK4_parameters,   only: dimH2,step,integration_time,ERR,safe_factor,lelectric, &
                                      hE_0,hw_e,lpulse_e,tau_e,delay_e,lmagnetic,hw1_m,hw_m, &
                                      lpulse_m,tau_m,delay_m,polarization_e,polarization_vec_e, &
                                      polarization_m,polarization_vec_m
    use mod_RK_matrices,        only: A,id,id2,M1,c1,c2,build_identity
    use mod_imRK4,              only: iterate_Zki,calculate_step_error,magnetic_field,vector_potential
    use mod_BrillouinZone,      only: realBZ
    use mod_parameters,         only: dimHsc,output,lprintfieldonly
    use mod_system,             only: System_type
    use mod_expectation,        only: expec_val_n, expec_H_n, expec_L_n
    use mod_Umatrix,            only: update_Umatrix
    use mod_magnet,             only: mzd0,mpd0,rhod0,rho0
    use mod_tools,              only: KronProd,diagonalize,lwork
    use mod_superconductivity,  only: allocate_supercond_variables
    use mod_hamiltonian,        only: calchk,hamilt_local,lfullhk,h0,fullhk
    use mod_checkpoint,         only: save_state,recover_state
    use mod_io,                 only: log_warning
    use mod_time_propagator_io, only: create_time_prop_files,write_header_time_prop,write_field,write_time_prop_files
    use mod_mpi_pars,           only: rFreq,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,FieldComm,ierr
    implicit none
    type(System_type), intent(in)   :: s

    integer(int32)                              :: i,mu,mud,it,n,counter,iter_rej,iter_tot
    real(dp)                                    :: t,tm,t1,t2
    real(dp)                                    :: pinv,h_new,h_old,ERR_old,ERR_kn
    complex(dp), dimension(dimHsc)              :: Yn,Yn_hat,Yn_new,Yn_e,Yn_hat_e,Yn_new_e
    real(dp),    dimension(3)                   :: b_field,A_t,b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2

    logical  :: use_checkpoint
    real(dp) :: t_cp,step_cp

    integer(int64)                          :: ik
    integer(int32)                          :: ncount,ios,ilaenv
    real(dp)                                :: weight, kp(3)
    real(dp),    dimension(s%nOrb,s%nAtoms) :: expec_0, expec_z
    complex(dp), dimension(s%nOrb,s%nAtoms) :: expec_p
    real(dp),    dimension(s%nOrb,s%nAtoms) :: expec_d
    real(dp),    dimension(dimHsc)          :: eval
    complex(dp), dimension(dimHsc,dimHsc)   :: hk,hkev,hamilt_nof
    complex(dp), dimension(dimHsc,dimHsc,realBZ%workload)  :: evec_kn,evec_kn_temp
    real(dp),    dimension(dimHsc,realBZ%workload)         :: eval_kn

    real(dp),    dimension(s%nOrb,s%nAtoms) :: rho_t,mx_t,my_t,mz_t
    complex(dp), dimension(s%nOrb,s%nAtoms) :: mp_t
    real(dp),    dimension(s%nAtoms)        :: rhod_t,mxd_t,myd_t,mzd_t
    complex(dp), dimension(s%nAtoms)        :: mpd_t
    real(dp),    dimension(s%nOrb,s%nAtoms) :: delta_sc_t
    real(dp),    dimension(2,s%nAtoms)      :: lxm,lym,lzm
    real(dp),    dimension(2,s%nAtoms)      :: lxm_t,lym_t,lzm_t
    real(dp)                                :: E_t, E_0
    complex(dp)                             :: exp_eval
   
    external :: MPI_Allreduce,MPI_Barrier,ilaenv,endTITAN

    if(rFreq(1) == 0) &
      write(output%unit_loop,"('CALCULATING TIME-PROPAGATION')")
 
    ! number of elements in the MPI communication
    ncount = s%nOrb*s%nAtoms

    ! Dimensions for RK method
    dimH2  = 2*dimHsc

    allocate( id(dimHsc,dimHsc),id2(dimH2,dimH2),M1(dimH2,dimH2) )

    ! Building identities
    call build_identity(dimHsc,id)
    call build_identity(size(A,1)*dimHsc,id2)

    ! Building matrix M1
    M1 = KronProd(size(A,1),size(A,1),dimHsc,dimHsc,A,id)

    ! Checking for checkpoints
    use_checkpoint = recover_state(rFreq(1),dimHsc,realBZ%workload,t_cp,step_cp,eval_kn,evec_kn)
    if(use_checkpoint) then
      t = t_cp
      step = step_cp
    else
      t = 0._dp
    end if

    ! Creating files and writing headers
    if(.not.(use_checkpoint)) call create_time_prop_files()

    if(lprintfieldonly) then
      if(rFreq(1) == 0) &
        call write_field()
      return
    end if

    ! Getting lwork for diagonalization
    lwork = (ilaenv( 1, 'zhetrd', 'VU', dimHsc, -1, -1, -1 )+1)*dimHsc

    if(rFreq(1) == 0) &
      write(output%unit_loop,"('[time_propagator] Starting propagation from t = ',es9.2)") t

    it = 0       ! Counter of accepted iterations
    iter_tot = 0 ! Counter of total number of iterations (rejected + accepted)

    ! Time propagation over t, kpoints, eigenvectors(Yn) for each k
    t_loop: do while (t <= integration_time)
      t = t + step
      if(rFreq(1) == 0) &
        write(output%unit_loop,"('[time_propagator] Time: ',es10.3,' of ',es10.3)", advance='no') t, integration_time

      ! Build local hamiltonian with m(t) and n(t)
      call hamilt_local(s)

      counter = 0  ! Counter for the calculation of error in the step size for each time t
      iter_rej = 0 ! Counter of rejected steps (for each accepted one)
      ! Propagation loop for a given time t, calculating the optimal step size
      find_step: do

        ! Calculating required time-dependent fields
        b_field  = 0._dp
        b_fieldm = 0._dp
        b_field1 = 0._dp
        b_field2 = 0._dp
        A_t  = 0._dp 
        A_tm = 0._dp 
        A_t1 = 0._dp 
        A_t2 = 0._dp 
        tm = t - step
        t1 = t + step * c1
        t2 = t + step * c2
        if (lmagnetic) then
          call magnetic_field(t,b_field)
          call magnetic_field(tm,b_fieldm)
          call magnetic_field(t1,b_field1)
          call magnetic_field(t2,b_field2)
        end if

        if (lelectric) then
          call vector_potential(t,A_t)
          call vector_potential(tm, A_tm)
          call vector_potential(t1, A_t1)
          call vector_potential(t2, A_t2)
        end if

        rho_t  = 0._dp
        mp_t   = cZero
        mz_t   = 0._dp

        E_t    = 0._dp

        Lxm_t  = 0._dp
        Lym_t  = 0._dp
        Lzm_t  = 0._dp

        ERR    = 0._dp

        delta_sc_t = 0._dp

        !$omp parallel do default(none) schedule(dynamic,1) &
        !$omp& private(Yn_e,Yn_new_e,Yn_hat_e,ik,n,i,kp,weight,hk,hkev,hamilt_nof,exp_eval,ERR_kn,Yn,Yn_new,Yn_hat,eval,expec_0,expec_p,expec_z,expec_d,E_0,lxm,lym,lzm) &
        !$omp& shared(counter,step,s,t,it,id,dimHsc,realBZ,lfullhk,h0,fullhk,evec_kn_temp,evec_kn,eval_kn,use_checkpoint,b_field,A_t,b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2) &
        !$omp& reduction(+:rho_t,mp_t,mz_t,E_t,Lxm_t,Lym_t,Lzm_t,delta_sc_t,ERR)
        kpoints_loop: do ik = 1, realBZ%workload
          kp = realBZ%kp(:,ik)
          weight = realBZ%w(ik)   
          ! Calculating the hamiltonian for a given k-point
          if( lfullhk ) then
            hk = h0 + fullhk(:,:,ik)
          else
            hk = h0 + calchk(s,kp)
          end if

          ! For t=0.0, diagonalize the hamiltonian
          if ((it==0).and.(.not.use_checkpoint)) then
            hkev = hk
            ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
            call diagonalize(dimHsc,hkev,eval)
            eval_kn(:,ik) = eval(:)
            !>>>>> find eigenvalues again
            ! Improve step control by diagonalizing the Hamiltonian after some number of steps:
            ! calling H(t) instead of H(k), at t=0, H(0)= H(k)
          end if

          evs_loop: do n = 1, dimHsc
            exp_eval = exp(-cI*eval_kn(n,ik)*t) 
            if ((it==0).and.(.not.use_checkpoint)) then
              Yn(:)= hkev(:,n)
              !----------------------------------------------------------------------------------------!
              !-------- These steps are to control the error by solving for c^~ instead of c ----------!
              !----------------------------------------------------------------------------------------!
              ! 1- find Yn^~ = Yn * e^(i*En*t/hbar) : for each n and k point (only once).
              ! 2- propagate Yn^~ using H^~ = H(t) - En.
              ! 3- find Yn = Yn^~ * exp( (-i*En*t/hbar) ) at each step.
              ! 4- calculate the error at each step from Yn but propagate with Yn^~.
              !----------------------------------------------------------------------------------------!
              ! Geting the intitial vector Yn^~, setting Yn to Yn^~ for propagation, hbar = 1
              Yn = Yn * conjg(exp_eval)
                
            else
              ! Is the propagated vector Yn^~
              Yn(:) = evec_kn(:,n,ik)
            end if

            ! Calculating the time-dependent Hamiltonian without the time-dependent field
            ! hext_t must be subtracted from hamilt_t
            hamilt_nof = ( eval_kn(n,ik) * id ) - hk

            call iterate_Zki(s,b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2,hamilt_nof,kp,eval_kn(n,ik),step,Yn_new,Yn,Yn_hat)
            ! Note: all iteration outputs correspond to Yn^~ 

            ! Getting Yn's again; Yn = Yn^~ * exp( (-i*En*t/hbar) ), hbar=1, rename Yn to Yn_e
            Yn_e(:)     = Yn     * exp_eval
            Yn_new_e(:) = Yn_new * exp_eval
            Yn_hat_e(:) = Yn_hat * exp_eval

            ! Calculating expectation values of the nth eigenvector
            ! Note: use the Yn_new_e or Yn_nw >>> should give the same result

            ! calculating expectation value of magnetization in eigenvector (n)
            call expec_val_n(s,dimHsc,Yn_new_e,eval_kn(n,ik),expec_0,expec_p,expec_z,expec_d)

            ! calculating expectation value of the T.D Hamiltonian in eigenvector (n)
            call expec_H_n(s,b_field,A_t,hk,kp,Yn_new_e,eval_kn(n,ik),E_0)


            ! calculating expectation value of angular momentum in eigenvector (n)
            call expec_L_n(s, dimHsc, Yn_new_e, eval_kn(n,ik), lxm, lym, lzm)


            rho_t = rho_t + expec_0  * weight 
            mp_t  = mp_t  + expec_p  * cmplx(weight,0._dp,dp)
            mz_t  = mz_t  + expec_z  * weight 

            ! Superconducting order parameter
            delta_sc_t = delta_sc_t + expec_d*weight

            ! Expectation value of the energy
            E_t = E_t + E_0 * weight

            ! Orbital angular momentum
            Lxm_t  = Lxm_t   + lxm  * weight
            Lym_t  = Lym_t   + lym  * weight
            Lzm_t  = Lzm_t   + lzm  * weight

            ! Calculation of the error and the new step size.
            ! Note: use Yn_e, Yn_new_e, Yn_hat_e.
            call calculate_step_error(Yn_e,Yn_new_e,Yn_hat_e,ERR_kn)

            ERR = ERR + ERR_kn * weight

            ! Storing temporary propagated vector before checking if it's Accepted.
            ! Note: save the Yn^~ outputs to be propagated.
            evec_kn_temp(:,n,ik) = Yn_new(:)

          end do evs_loop
        end do kpoints_loop
        !$omp end parallel do

        call MPI_Allreduce(MPI_IN_PLACE, rho_t, ncount    , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, mz_t , ncount    , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, mp_t , ncount    , MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, Lxm_t, 2*s%nAtoms, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, Lym_t, 2*s%nAtoms, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, Lzm_t, 2*s%nAtoms, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, delta_sc_t, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, E_t  , 1         , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, ERR  , 1         , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

        ERR = sqrt(ERR)
        ! Find the new step size h_new
        ! h_new = safe_factor * h_used / (ERR)^(1/p+1) where p = 2*s, safe_factor is some safety factor 
        ! safe_factor = 0.9 * (2*K_max + 1)/ (2*K_max + NEWT) where NEWT is the number of newton iterations
        pinv = 0.5_dp/size(A,1) 
        ! h_new = safe_factor * step * (ERR)**(-1.0/(2.0*s)) !! simpler formula
        if (counter /= 0) then
          h_new = safe_factor * step * (ERR)**(-pinv) * (step/h_old) * (ERR_old/ERR)**(pinv)!! Gustafsson (1994) formula          
        else 
          h_new = safe_factor * step * (ERR)**(-pinv)
        end if 
        ! Tests if the step size is good:
        ! the condition (h_new < safe_factor * step) is equivalent to the condition (ERR > 1)
        ! if ( h_new < 0.9 * step) then
           ! do while ( h_new < step) 
        ! save ERR from iterate_Zki to ERR_old
        ERR_old = ERR

        iter_tot = iter_tot + 1 ! Counter of total number of iterations
        if ( ERR > 1._dp) then
          ! Error is too large! Reject step and
          ! repeat the calculation using h_new
          t = t - step + h_new
          h_old = step
          step = h_new
          ! if(rFreq(1) == 0) write(*,*) "Rejected", t, step, ERR
          iter_rej = iter_rej + 1 ! Counter of rejected steps
          ! this condition seems to mantain a small step size
          ! else if (step <= h_new <= 1.2*step) then
          ! step = step
        else    
          ! Error is low! Accept the step
          ! Update the eignvectors array of dimension ( dimHsc * dimHsc , k )
          evec_kn = evec_kn_temp
          step = h_new
          exit
        end if
          
      end do find_step

      if(rFreq(1) == 0) then
        write(output%unit_loop,"(' (',i0,' rejected iterations)')") iter_rej
        flush(output%unit_loop)
      end if

      ! if(rFreq(1) == 0) write(*,*)  "Accepted", t, step, ERR 

      mx_t = real(mp_t)
      my_t = aimag(mp_t)

      ! obtaining expectation values for d-orbitals
      rhod_t = 0._dp
      mpd_t  = 0._dp
      mxd_t  = 0._dp
      myd_t  = 0._dp
      mzd_t  = 0._dp
      do i = 1, s%nAtoms
        do mud = 1,s%ndOrb
          mu = s%dOrbs(mud)
          rhod_t(i) = rhod_t(i) + rho_t(mu,i)
          mpd_t (i) = mpd_t (i) + mp_t (mu,i)
          mxd_t (i) = mxd_t (i) + mx_t (mu,i)
          myd_t (i) = myd_t (i) + my_t (mu,i)
          mzd_t (i) = mzd_t (i) + mz_t (mu,i)
        end do
      end do
 
      ! Update U-term of the local hamiltonian
      call update_Umatrix(mzd_t,mzd0,mpd_t,mpd0,rhod_t,rhod0,rho_t,rho0,s)

      ! Writing results to file
      if(rFreq(1) == 0) &
        call write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,b_field,A_t,E_t,lxm_t,lym_t,lzm_t) 

      counter = counter + 1
      it = it + 1 ! Counter of accepted iterations

      ! Checking for "save" file to trigger checkpoint
      open(unit=911, file="save", status='old', iostat=ios)
      if(ios==0) then
        close(911)
        call save_state(rFreq(1),dimHsc,realBZ%workload,t,step,eval_kn,evec_kn)
        call MPI_Barrier(FieldComm, ierr)
        if(rFreq(1) == 0) &
          call execute_command_line('rm save')
        call endTITAN()
      end if
    end do t_loop

    ! Creating checkpoint in the last state
    call save_state(rFreq(1),dimHsc,realBZ%workload,t,step,eval_kn,evec_kn)

    deallocate( id,id2,M1 )

    if(rFreq(1) == 0) &
      write(output%unit_loop,"('[time_propagator] Integration time reached. ',i0,' total iterations, with ',i0,' accepted.')") iter_tot,it

    if(lelectric) then
      deallocate(polarization_e,polarization_vec_e,hE_0,hw_e)
      if(lpulse_e) deallocate(tau_e,delay_e) 
    end if

    if(lmagnetic) then
      deallocate(polarization_m,polarization_vec_m,hw1_m,hw_m)
      if(lpulse_m) deallocate(tau_m,delay_m) 
    end if

  end subroutine time_propagator
end module mod_time_propagator

