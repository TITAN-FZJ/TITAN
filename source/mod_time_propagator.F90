module mod_time_propagator
  !! module for time-propagation parameters
  use mod_kind, only: dp
  implicit none
  logical      :: lelectric, lmagnetic, lpulse_e, lpulse_m
  !! Logical variables for choosing which field is applied

  integer      :: npulse_e
  !! Number of electric pulses
  real(dp),     dimension(:) , allocatable :: hE_0
  !! Intensity of electric fields
  real(dp),     dimension(:) , allocatable :: hw_e
  !! Frequency (hwt) of electric fields
  real(dp),     dimension(:) , allocatable :: tau_e
  !! Pulse length of electric fields
  real(dp),     dimension(:) , allocatable :: delay_e
  !! Time delay in electric fields pulses
  character(len=1), dimension(:) , allocatable ::  polarization_e 
  !! Polarization of electric field
  real(dp), dimension(:,:,:) , allocatable ::  polarization_vec_e
  !! Polarization vector, inphase (cos) and out-of-phase (sin) of electric field

  integer      :: npulse_m
  !! Number of magnetic pulses
  real(dp),     dimension(:) , allocatable :: hw1_m
  !! Intensity of magnetic fields
  real(dp),     dimension(:) , allocatable :: hw_m
  !! Frequency (w.t) of magnetic fields
  real(dp),     dimension(:) , allocatable :: tau_m
  !! Pulse length of magnetic fields
  real(dp),     dimension(:) , allocatable :: delay_m
  !! Time delay in magnetic fields pulses
  character(len=1), dimension(:) , allocatable ::  polarization_m
  !! Polarization of magnetic field
  real(dp), dimension(:,:,:) , allocatable ::  polarization_vec_m
  !! Polarization vector, inphase (cos) and out-of-phase (sin) of magnetic field

  real(dp) :: integration_time
  !! Real integration time 
  real(dp) :: step
  !! Step size
  real(dp) :: sc_tol
  !! Time propagation self consistency tolerence 
  real(dp) :: abs_tol, rel_tol, safe_factor
  !! Step size control error(ERR) tolerance
  integer  :: dimH2
  !! Dimension: 2*dimension of the Hamiltonian (2*dimHsc)
  real(dp) :: ERR
  !! Error for the calculation of the step size in time propagation
  real(dp) :: time_conv = 4.84e-5_dp
  !! Conversion of time units to picosecond

  !! Runge-Kutta parameters and Butcher tableu, Pauli matricies,
  !! A_inverse matrix, M1 matrix, identity matricies, coefficients b, d, c and c_avg.
  complex(dp), dimension(:,:), allocatable :: id(:,:),id2(:,:)
  complex(dp), dimension(:,:), allocatable :: M1
  real(dp),    dimension(2,2), parameter   :: A_inverse= reshape([ 3._dp, 0.46410162_dp, -6.46410162_dp, 3._dp ],[ 2,2 ],Order=[ 2,1 ])
  complex(dp), dimension(2,2), parameter   :: A= reshape([ 0.25_dp, -0.03867513_dp, 0.53867513_dp, 0.25_dp ],[ 2,2 ],Order=[ 2,1 ])
  real(dp),    dimension(2)  , parameter   :: b= [ 0.5_dp, 0.5_dp ]
  real(dp)                   , parameter   :: d1= -1.73205081_dp
  real(dp)                   , parameter   :: d2=  1.73205081_dp
  real(dp)                   , parameter   :: d1_hat= 6.4641018_dp
  real(dp)                   , parameter   :: d2_hat= -0.46410171_dp
  real(dp)                   , parameter   :: c1= 0.21132486540518713_dp
  real(dp)                   , parameter   :: c2= 0.5288675134594812_dp
  real(dp)                   , parameter   :: c_avg= (c1+c2)*0.5_dp

contains

  subroutine time_propagator(s)
    use mod_kind,               only: dp,int32,int64
    use mod_constants,          only: cZero,cI
    use mod_BrillouinZone,      only: realBZ
    use mod_parameters,         only: dimHsc,output,lprintfieldonly
    use mod_system,             only: System_type
    use mod_expectation,        only: expec_val_n,expec_H_n,expec_L_n,expec_torque_n
    use mod_Umatrix,            only: update_Umatrix
    use mod_magnet,             only: mzd0,mpd0,rhod0,rho0,mdvec_cartesian,mx,my,mz
    use mod_tools,              only: KronProd,diagonalize,lwork,build_identity
    use mod_superconductivity,  only: allocate_supercond_variables
    use mod_hamiltonian,        only: calchk,hamilt_local,lfullhk,h0,fullhk
    use mod_logging,            only: log_warning
    use mod_mpi_pars,           only: rFreq,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_DOUBLE_COMPLEX,MPI_SUM,FreqComm,FieldComm,ierr
    use mod_torques,            only: SO_torque_operator,xc_torque_operator
    implicit none
    type(System_type), intent(in)   :: s

    integer(int32)                              :: i,mu,mud,it,n,counter,iter_rej,iter_tot
    real(dp)                                    :: t,tm,t1,t2
    real(dp)                                    :: pinv,h_new,h_old,ERR_old,ERR_kn
    complex(dp), dimension(dimHsc)              :: Yn,Yn_hat,Yn_new,Yn_e,Yn_hat_e,Yn_new_e
    complex(dp), dimension(:,:,:), allocatable  :: evec_kn,evec_kn_temp
    real(dp),    dimension(:,:),   allocatable  :: eval_kn
    real(dp),    dimension(3)                   :: b_field,A_t,b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2

    logical  :: use_checkpoint
    real(dp) :: t_cp,step_cp

    integer(int64)                          :: ik
    integer(int32)                          :: ncount,ios,ilaenv
    real(dp)                                :: weight, kp(3)
    real(dp),    dimension(s%nOrb,s%nAtoms) :: expec_0, expec_z
    complex(dp), dimension(s%nOrb,s%nAtoms) :: expec_p
    real(dp),    dimension(s%nOrb,s%nAtoms) :: expec_d
    real(dp),    dimension(:),  allocatable :: eval
    complex(dp), dimension(:,:),allocatable :: hk,hkev,hamilt_nof

    real(dp),    dimension(s%nOrb,s%nAtoms) :: rho_t,mx_t,my_t,mz_t
    real(dp),    dimension(s%nOrb,s%nAtoms) :: mx_t_prev,my_t_prev,mz_t_prev
    complex(dp), dimension(s%nOrb,s%nAtoms) :: mp_t
    real(dp),    dimension(s%nAtoms)        :: rhod_t,mxd_t,myd_t,mzd_t
    complex(dp), dimension(s%nAtoms)        :: mpd_t
    real(dp),    dimension(s%nOrb,s%nAtoms) :: delta_sc_t
    real(dp),    dimension(2,s%nAtoms)      :: lxm,lym,lzm
    real(dp),    dimension(2,s%nAtoms)      :: lxm_t,lym_t,lzm_t

    complex(dp), dimension(s%nOrb2,s%nOrb2,3,s%nAtoms) :: tso_op,txc_op
    real(dp),    dimension(3,2,s%nAtoms)               :: tso,tso_t
    real(dp),    dimension(3,s%nOrb,s%nAtoms)          :: txc,txc_t
    real(dp),    dimension(3,s%nOrb,s%nAtoms)          :: dmdt

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
    allocate( hamilt_nof(dimHsc,dimHsc),hk(dimHsc,dimHsc),hkev(dimHsc,dimHsc),eval(dimHsc),eval_kn(dimHsc,realBZ%workload),evec_kn(dimHsc,dimHsc,realBZ%workload),evec_kn_temp(dimHsc,dimHsc,realBZ%workload) )

    ! Building identities
    call build_identity(dimHsc,id)
    call build_identity(size(A,1)*dimHsc,id2)

    ! Building matrix M1
    M1 = KronProd(size(A,1),size(A,1),dimHsc,dimHsc,A,id)

    ! Checking for checkpoints
    use_checkpoint = recover_state(rFreq(1),s,dimHsc,realBZ%workload,t_cp,step_cp,eval_kn,evec_kn,mx_t,my_t,mz_t)
    if(use_checkpoint) then
      t = t_cp
      step = step_cp
      mx_t_prev = mx_t
      my_t_prev = my_t
      mz_t_prev = mz_t
    else
      t = 0._dp
      ! Storing initial magnetization to use as previous step at t=0
      mx_t_prev = mx
      my_t_prev = my
      mz_t_prev = mz
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

    ! Getting torque operators
    call SO_torque_operator(tso_op)
    call xc_torque_operator(txc_op)


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

        tso_t = 0._dp
        txc_t = 0._dp

        ERR    = 0._dp

        delta_sc_t = 0._dp

        !$omp parallel do default(none) schedule(dynamic,1) &
        !$omp& private(Yn_e,Yn_new_e,Yn_hat_e,ik,n,i,kp,weight,hk,hkev,hamilt_nof,exp_eval,ERR_kn,Yn,Yn_new,Yn_hat,eval,expec_0,expec_p,expec_z,expec_d,E_0,lxm,lym,lzm,tso,txc) &
        !$omp& shared(counter,step,s,t,it,id,tso_op,txc_op,dimHsc,realBZ,lfullhk,h0,fullhk,evec_kn_temp,evec_kn,eval_kn,use_checkpoint,lmagnetic,b_field,lelectric,A_t,b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2) &
        !$omp& reduction(+:rho_t,mp_t,mz_t,E_t,Lxm_t,Lym_t,Lzm_t,delta_sc_t,tso_t,txc_t,ERR)
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

            ! Getting Yn again; Yn = Yn^~ * exp( (-i*En*t/hbar) ), hbar=1, rename Yn to Yn_e
            Yn_e(:)     = Yn     * exp_eval
            Yn_new_e(:) = Yn_new * exp_eval
            Yn_hat_e(:) = Yn_hat * exp_eval

            ! Calculating expectation values of the nth eigenvector
            ! Note: use the Yn_new_e or Yn_nw >>> should give the same result

            ! calculating expectation value of magnetization in eigenvector (n)
            call expec_val_n(s,dimHsc,Yn_new_e,eval_kn(n,ik),expec_0,expec_p,expec_z,expec_d)

            ! calculating expectation value of the T.D Hamiltonian in eigenvector (n)
            call expec_H_n(s,lmagnetic,b_field,lelectric,A_t,hk,kp,Yn_new_e,eval_kn(n,ik),E_0)

            ! calculating expectation value of angular momentum in eigenvector (n)
            call expec_L_n(s,dimHsc,Yn_new_e,eval_kn(n,ik),lxm,lym,lzm)

            ! calculating expectation value of so- and xc-torques in eigenvector (n)
            call expec_torque_n(s,dimHsc,Yn_new_e,eval_kn(n,ik),tso_op,txc_op,tso,txc)

            rho_t = rho_t + expec_0 * weight 
            mp_t  = mp_t  + expec_p * cmplx(weight,0._dp,dp)
            mz_t  = mz_t  + expec_z * weight 

            ! Superconducting order parameter
            delta_sc_t = delta_sc_t + expec_d*weight

            ! Expectation value of the energy
            E_t = E_t + E_0 * weight

            ! Orbital angular momentum
            Lxm_t  = Lxm_t + lxm * weight
            Lym_t  = Lym_t + lym * weight
            Lzm_t  = Lzm_t + lzm * weight

            ! Spin-orbit and xc Torque
            tso_t = tso_t + tso * weight
            txc_t = txc_t + txc * weight

            ! Calculation of the error and the new step size.
            ! Note: use Yn_e, Yn_new_e, Yn_hat_e.
            call calculate_step_error(Yn_e,Yn_new_e,Yn_hat_e,ERR_kn)

            ERR = ERR + ERR_kn * weight

            ! Storing temporary propagated vector before checking if it is accepted.
            ! Note: save the Yn^~ outputs to be propagated.
            evec_kn_temp(:,n,ik) = Yn_new(:)

          end do evs_loop
        end do kpoints_loop
        !$omp end parallel do

        call MPI_Allreduce(MPI_IN_PLACE, rho_t     , ncount           , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, mz_t      , ncount           , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, mp_t      , ncount           , MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, Lxm_t     , 2*s%nAtoms       , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, Lym_t     , 2*s%nAtoms       , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, Lzm_t     , 2*s%nAtoms       , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, tso_t     , 3*2*s%nAtoms     , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, txc_t     , 3*s%nOrb*s%nAtoms, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, delta_sc_t, ncount           , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, E_t       , 1                , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, ERR       , 1                , MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

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

      if(rFreq(1) == 0) &
        write(output%unit_loop,"(' (',i0,' rejected iterations)')") iter_rej

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
        do mud = 1,s%Types(s%Basis(i)%Material)%ndOrb
          mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
          rhod_t(i) = rhod_t(i) + rho_t(mu,i)
          mpd_t (i) = mpd_t (i) + mp_t (mu,i)
          mxd_t (i) = mxd_t (i) + mx_t (mu,i)
          myd_t (i) = myd_t (i) + my_t (mu,i)
          mzd_t (i) = mzd_t (i) + mz_t (mu,i)
        end do
        mdvec_cartesian(1,i) = mxd_t(i)
        mdvec_cartesian(2,i) = myd_t(i)
        mdvec_cartesian(3,i) = mzd_t(i)
      end do

      ! Update U-term of the local hamiltonian
      call update_Umatrix(mzd_t,mzd0,mpd_t,mpd0,rhod_t,rhod0,rho_t,rho0,s)

      ! Updating xc-torque operator
      call xc_torque_operator(txc_op)

      ! Calculating total torque via dM/dt
      dmdt(1,:,:) = (mx_t(:,:)-mx_t_prev(:,:))/step
      dmdt(2,:,:) = (my_t(:,:)-my_t_prev(:,:))/step
      dmdt(3,:,:) = (mz_t(:,:)-mz_t_prev(:,:))/step
      ! Updating previous step with current one
      mx_t_prev = mx_t
      my_t_prev = my_t
      mz_t_prev = mz_t

      ! Writing results to file
      if(rFreq(1) == 0) &
        call write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,b_field,A_t,E_t,lxm_t,lym_t,lzm_t,tso_t,txc_t,dmdt)

      counter = counter + 1
      it = it + 1 ! Counter of accepted iterations

      ! Checking for "save" file to trigger checkpoint
      open(unit=911, file="save", status='old', iostat=ios)
      if(ios==0) then
        close(911)
        call save_state(rFreq(1),s,dimHsc,realBZ%workload,t,step,eval_kn,evec_kn,mx_t,my_t,mz_t)
        call MPI_Barrier(FieldComm, ierr)
        if(rFreq(1) == 0) &
          call execute_command_line('rm save')
        call endTITAN()
      end if
    end do t_loop

    ! Creating checkpoint in the last state
    call save_state(rFreq(1),s,dimHsc,realBZ%workload,t,step,eval_kn,evec_kn,mx_t,my_t,mz_t)

    deallocate( id,id2,M1 )
    deallocate( hamilt_nof,hk,hkev,eval,eval_kn,evec_kn,evec_kn_temp )

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


  subroutine iterate_Zki(s,b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2,hamilt_nof,kp,eval,step,Yn_new,Yn,Yn_hat)
    !! Calculates the vectors Z_ki
    use mod_kind,             only: dp
    use mod_parameters,       only: dimHsc
    use mod_system,           only: System_type
    use mod_tools,            only: vec_norm, LS_solver
    use mod_constants,        only: cZero
    implicit none
    ! define the variables: 
    type(System_type),                     intent(in)    :: s
    real(dp),    dimension(3),             intent(in)    :: b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2
    complex(dp), dimension(dimHsc,dimHsc), intent(in)    :: hamilt_nof
    real(dp),                              intent(in)    :: kp(3)
    real(dp),                              intent(in)    :: eval
    complex(dp), dimension(dimHsc),        intent(inout) :: Yn_new, Yn_hat, Yn
   
    integer                             :: k
    real(dp)                            :: step, error, norm_k_old, norm_k_new, theta_k, eta_k, tol
    complex(dp), dimension(dimH2)       :: deltaZ_k_old, deltaZ_k  
    complex(dp), dimension(dimH2,dimH2) :: M2n
    complex(dp), dimension(dimH2)       :: Fun 
    complex(dp), dimension(dimH2)       :: Z_k 

    ! z = (g - Y)
    ! z = 0.0 is the solution
    ! Initial value of z_k (Taylor expansion point) 
    Z_k   = cZero
    k     = 0
    error = 1._dp
    !
    ! Solve the linear system of equations of size 2 * n ,
    ! where n is the dimension of the hamiltonian:
    !   (I- hA\otimes J)  \Delta z^k = -z^k+ h(A\otimes I)F(z^k)
    !   ----------------               -------------------------
    !        M2n                               Fun
    !
    sc_loop: do while (error >= sc_tol) 
      ! Build M2n matrix
      call M_2n(s,b_fieldm,A_tm,hamilt_nof,kp,eval,Yn,step,M2n)
      ! Build Fun vector
      call Fsystem(b_field1,A_t1,b_field2,A_t2,hamilt_nof,kp,Yn,step,Z_k,Fun)
      ! Solve the Ax=b linear equation (Eq. 24 of the notes)
      call LS_solver(dimH2,M2n,Fun)

      deltaZ_k  = Fun
      
      ! Save old deltaZ(k) to deltaZ(k-1)
      deltaZ_k_old= Z_k
      ! Update Z_k
      Z_k  = Z_k + deltaZ_k

      ! Obtaining new deltaZ(k)
      if (k >= 1) then
        norm_k_old = vec_norm(deltaZ_k_old, dimH2)
        norm_k_new = vec_norm(Z_k, dimH2)
        theta_k    = norm_k_new/norm_k_old
        if (abs(theta_k - 1._dp) < 1.e-15_dp) then 
          error = 0._dp
        else
          eta_k = theta_k/(1._dp-theta_k)
          tol   = abs(norm_k_new - norm_k_old)/min(norm_k_new, norm_k_old)
          error = eta_k*norm_k_new/tol
        end if 
      else
        error = 1._dp
        k     = k+1
      end if 
    end do sc_loop

    ! Find Yn_hat
    Yn_hat = Yn + d1_hat*Z_k(1:dimHsc) + d2_hat*Z_k(dimHsc+1:dimH2)
    ! Find Yn
    Yn_new = Yn + d1*Z_k(1:dimHsc) + d2*Z_k(dimHsc+1:dimH2)

  end subroutine iterate_Zki

  subroutine M_2n(s,b_fieldm,A_tm,hamilt_nof,kp,eval,Yn,step,M2n)
    !! Subroutine to build the matrix M_2n
    !! b_fieldm = b_field(t-step) , A_tm = A_t(t-step)
    use mod_kind,             only: dp
    use mod_parameters,       only: dimHsc
    use mod_system,           only: System_type
    use mod_tools,            only: KronProd
    implicit none
    type(System_type),                       intent(in)  :: s
    real(dp),    dimension(3),               intent(in)  :: b_fieldm,A_tm
    real(dp),                                intent(in)  :: step
    complex(dp), dimension(dimHsc,dimHsc),   intent(in)  :: hamilt_nof
    real(dp),                                intent(in)  :: kp(3)
    real(dp),                                intent(in)  :: eval
    complex(dp), dimension(dimHsc),          intent(in)  :: Yn
    complex(dp), dimension(dimH2,dimH2),     intent(out) :: M2n

    complex(dp), dimension(dimHsc,dimHsc)   :: Jacobian_t
    complex(dp), dimension(dimH2,dimH2) :: Kprod(dimH2,dimH2)

    ! Calculating the Jacobian at previous time
    call build_td_Jacobian(s,b_fieldm,A_tm,hamilt_nof,kp,eval,Yn,Jacobian_t)

    Kprod = KronProd(size(A,1),size(A,1),dimHsc,dimHsc,A,Jacobian_t)

    M2n = id2 - step * Kprod
  end subroutine M_2n

  subroutine Fsystem(b_field1,A_t1,b_field2,A_t2,hamilt_nof,kp,Yn,step,Z_k,Fun) 
    !! subroutine f(Z_k) that builds the right side of the linear system.
    use mod_kind,             only: dp
    use mod_constants,        only: cI,cZero
    use mod_parameters,       only: dimHsc
    use mod_hamiltonian,      only: build_hext
    implicit none
    real(dp),    dimension(3),             intent(in)  :: b_field1,A_t1,b_field2,A_t2
    real(dp),                              intent(in)  :: step
    complex(dp), dimension(dimHsc,dimHsc), intent(in)  :: hamilt_nof
    real(dp),                              intent(in)  :: kp(3)
    complex(dp), dimension(dimHsc),        intent(in)  :: Yn
    complex(dp), dimension(dimH2),         intent(in)  :: Z_k
    complex(dp), dimension(dimH2),         intent(out) :: Fun 
    complex(dp), dimension(dimHsc,dimHsc)  :: hamilt_t,hext_t
    complex(dp), dimension(dimHsc)         :: F1, F2, tempv
    complex(dp), dimension(dimH2)          :: temp
    
    external :: zgemv

    ! Building time dependent hamiltonian
    call build_hext(kp,lmagnetic,b_field1,lelectric,A_t1,hext_t)
    hamilt_t = hamilt_nof - hext_t

    !F1 = -cI * matmul(hamilt_t,Z_k(1:dimHsc) + Yn)
    tempv = Z_k(1:dimHsc)+Yn
    call zgemv('n',dimHsc,dimHsc,-cI,hamilt_t,dimHsc,tempv,1,cZero,F1,1)

    ! Building time dependent hamiltonian
    call build_hext(kp,lmagnetic,b_field2,lelectric,A_t2,hext_t)
    hamilt_t = hamilt_nof - hext_t

    ! F2 = -cI * matmul(hamilt_t,Z_k(dimHsc+1:dimH2) + Yn)
    tempv = Z_k(dimHsc+1:dimH2)+Yn
    call zgemv('n',dimHsc,dimHsc,-cI,hamilt_t,dimHsc,tempv,1,cZero,F2,1)

    Fun = [ F1, F2 ]

    ! Fun = -Z_k + step * matmul(M1,Fun)
    ! First set: temp = step * matmul(M1,Fun)
    call zgemv('n',dimH2,dimH2,cmplx(step,0._dp,dp),M1,dimH2,Fun,1,cZero,temp,1)
    Fun = -Z_k + temp

  end subroutine Fsystem


  subroutine build_td_Jacobian(s,b_field,A_t,hamilt_nof,kp,eval,Yn,Jacobian_t)
    !! build time dependent jacobian for each kp 
    use mod_kind,        only: dp
    use mod_constants,   only: cI
    use mod_parameters,  only: dimHsc
    use mod_system,      only: System_type
    use mod_hamiltonian, only: build_hext
    implicit none
    type(System_type),                     intent(in)  :: s
    real(dp),    dimension(3),             intent(in)  :: b_field,A_t
    complex(dp), dimension(dimHsc,dimHsc), intent(in)  :: hamilt_nof
    real(dp),                              intent(in)  :: kp(3)
    real(dp),                              intent(in)  :: eval
    complex(dp), dimension(dimHsc),        intent(in)  :: Yn
    complex(dp), dimension(dimHsc,dimHsc), intent(out) :: Jacobian_t

    complex(dp), dimension(dimHsc,dimHsc) :: hamilt_t,hext_t,dHdc

    ! Building time dependent hamiltonian
    call build_hext(kp,lmagnetic,b_field,lelectric,A_t,hext_t)
    hamilt_t = hamilt_nof - hext_t

    call build_term_Jacobian(s, eval, Yn, dHdc)

    Jacobian_t = -cI*(hamilt_t + dHdc)

  end subroutine build_td_Jacobian


  subroutine build_term_Jacobian(s,eval,Yn,dHdc)
    !! Calculate the last term of the Jacobian
    !! Given by \sum_j <i| dH/dc^n_k |j> * c_j^n(t)
    use mod_kind,              only: dp
    use mod_parameters,        only: dimHsc,isigmamu2n,eta
    use mod_system,            only: System_type
    use mod_distributions,     only: fd_dist
    use mod_constants,         only: cZero,pauli_mat,pi
    use mod_superconductivity, only: lsupercond
    use mod_mpi_pars,          only: abortProgram
    implicit none
    type(System_type),                     intent(in)  :: s
    real(dp),                              intent(in)  :: eval
    complex(dp), dimension(dimHsc),        intent(in)  :: Yn
    complex(dp), dimension(dimHsc,dimHsc), intent(out) :: dHdc

    integer      :: i,mu,nu,mud,nud,s1,s2,s3,s4,alpha

    if(lsupercond) & 
      call abortProgram("[build_term_Jacobian] Not implemented for superconductivity yet!") ! TODO

    dHdc = cZero
    do i=1,s%nAtoms
      do mud=1,s%Types(s%Basis(i)%Material)%ndOrb
        mu=s%Types(s%Basis(i)%Material)%dOrbs(mud)
        do s1=1,2
          do s2=1,2
            dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,mu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,mu)) + s%Basis(i)%Un*Yn(isigmamu2n(i,s1,mu))*conjg(Yn(isigmamu2n(i,s2,mu)))
            do nud=1,s%Types(s%Basis(i)%Material)%ndOrb
              nu=s%Types(s%Basis(i)%Material)%dOrbs(nud)
              dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) - 0.5_dp*s%Basis(i)%Un*Yn(isigmamu2n(i,s1,mu))*conjg(Yn(isigmamu2n(i,s2,nu)))
              do s3=1,2
                do s4=1,2
                  do alpha = 1,3
                    dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) - 0.5_dp*s%Basis(i)%Um*pauli_mat(s1,s3,alpha)*Yn(isigmamu2n(i,s3,mu))*pauli_mat(s4,s2,alpha)*conjg(Yn(isigmamu2n(i,s4,nu)))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    dHdc = fd_dist(s%Ef, 1._dp/(pi*eta), eval) * dHdc

  end subroutine build_term_Jacobian


  subroutine calculate_step_error(Yn,Yn_new,Yn_hat,ERR_kn)
    !! subroutine to calculate the error in the step size control
    use mod_kind,             only: dp
    use mod_parameters,       only: dimHsc
    implicit none
    complex(dp), dimension(dimHsc), intent(in)  :: Yn, Yn_new, Yn_hat
    real(dp),                     intent(out) :: ERR_kn
    real(dp), dimension(dimHsc)                 :: Eta
    real(dp)                                  :: ERR_i
    integer                                   :: i

    ! Find the difference Eta in the two solutions Yn_new and Yn_hat
    Eta = abs(Yn_new - Yn_hat)
    ! sum over the error components to get the scaled norm (ERR)
    ! ERR= sqrt[ 1\n Sum_i{ abs(Yn - Yn_hat)/( abs_tol + rel_tol * Yn(i) ) }^2 ]
    ERR_kn = 0._dp
    ERR_loop: do i=1, dimHsc
      ERR_i = ( Eta(i)/(abs_tol+ rel_tol* abs(Yn(i)) ) )**2 
      ERR_kn = ERR_i + ERR_kn
    end do ERR_loop
    ERR_kn = ERR_kn/dimHsc
    
  end subroutine calculate_step_error

  subroutine magnetic_field(t,b_field)
    !! Subroutine builds magnetic field B(t)
    !! Using the tame pulse form as the vector potential below
    use mod_kind,             only: dp
    use mod_constants,        only: pi
    implicit none
    real(dp) , intent(in)  :: t
    real(dp) , intent(out) :: b_field(3)
    integer  :: np
    real(dp) :: delay, arg

    b_field(:) = 0._dp
    if(lpulse_m) then
      pulses_m: do np = 1, npulse_m
        if ((t >= delay_m(np)).and.(t <= tau_m(np)+delay_m(np))) then

          delay = 0.5_dp*tau_m(np) + delay_m(np)
          arg   = hw_m(np)*(t-delay)

          b_field(:) = cos(arg)*polarization_vec_m(np,1,:) + sin(arg)*polarization_vec_m(np,2,:)

          ! Cos-squared pulse
          b_field(:) =  b_field(:) * hw1_m(np) * 0.5_dp * ( cos(pi*(t-delay)/tau_m(np)) )**2  
          ! Gaussian pulse
          ! b_pulse = b_pulse * (0.5_dp * hw1_m * exp(-2*log(2*(t-delay/tau_m)**2))
          ! b_pulse = b_pulse * (0.5_dp * hw1_m * exp(-((t-delay)-4._dp*tau_m)**2/tau_m**2))

        end if
      end do pulses_m
    else
      b_field(:) = 0.5_dp*hw1_m(1)*( cos(hw_m(1)*t)*polarization_vec_m(1,1,:)+sin(hw_m(1)*t)*polarization_vec_m(1,2,:) )
    end if
  end subroutine magnetic_field


  subroutine vector_potential(t, A_t)
    !!  Subroutine builds vector potential A(t) = - integral(E(t)dt)
    !! Pulse:
    !!  From paper(DOI: 0.1038/s41567-019-0602-9) the vector potential for a cos^2 pulse is given by: 
    !!  A(t) = (-E_pump/w_pump) * ( cos(pi*t/tau_pump) )^2        * sin(w_pump*t), add delay_e to get:
    !!  A_t  = (-hE_0/hw_e)     * (A cos(pi*(t-delay_e)/tau_e) )^2 * sin(hw_e*t)
    !! center the vector potential at delay_e
    use mod_kind,             only: dp
    use mod_constants,        only: pi
    implicit none 
    real(dp) , intent(in)  :: t
    real(dp) , intent(out) :: A_t(3)
    integer  :: np
    real(dp) :: delay, arg

    A_t(:) = 0._dp
    if(lpulse_e) then
      pulses_e: do np = 1, npulse_e
        if ((t >= delay_e(np)).and.(t <= tau_e(np)+delay_e(np))) then 
          delay = 0.5_dp*tau_e(np) + delay_e(np)
          arg   = hw_e(np)*(t-delay)

          A_t(:) = sin(arg)*polarization_vec_e(np,1,:) - cos(arg)*polarization_vec_e(np,2,:)

          ! Cos-squared pulse:
          A_t(:) = A_t(:) * (-hE_0(np)/hw_e(np)) * ( cos(pi*(t-delay)/tau_e(np)) )**2 
          ! Gaussian pulse:
          ! A_t = A_t * (-hE_0/hw_e) * exp(-2*log(2*(t-delay/tau_e)**2))
        end if
      end do pulses_e
    else
      ! For the electric field cos(wt), vector potential is sin(wt)/w
      A_t(:) = -(hE_0(1)/hw_e(1))*(sin(hw_e(1)*t)*polarization_vec_e(1,1,:)-cos(hw_e(1)*t)*polarization_vec_e(1,2,:)) 
    end if
  end subroutine vector_potential

  subroutine electric_field(t, E_t)
    !!  Subroutine builds Electric field potential E(t) = - dA(t)/dt
    !! Pulse:
    !!  From paper(DOI: 0.1038/s41567-019-0602-9) the vector potential for a cos^2 pulse is given by: 
    !!  A(t) = (-E_pump/w_pump) * ( cos(pi*t/tau_pump) )^2        * sin(w_pump*t), add delay_e to get:
    !!  A_t  = (-hE_0/hw_e)     * (A cos(pi*(t-delay_e)/tau_e) )^2 * sin(hw_e*t)
    !! center the vector potential at delay_e
    use mod_kind,             only: dp
    use mod_constants,        only: pi
    implicit none 
    real(dp) , intent(in)  :: t
    real(dp) , intent(out) :: E_t(3)
    integer  :: np
    real(dp) :: delay, arg

    E_t(:) = 0._dp
    if(lpulse_e) then
      pulses_e: do np = 1, npulse_e
        if ((t >= delay_e(np)).and.(t <= tau_e(np)+delay_e(np))) then 
          delay = 0.5_dp*tau_e(np) + delay_e(np)
          arg   = hw_e(np)*(t-delay)

          E_t(:) = (cos(pi*(t-delay)/tau_e(np)))**2 * (cos(arg)*polarization_vec_e(np,1,:) + sin(arg)*polarization_vec_e(np,2,:)) &
                 - (pi/(tau_e(np)*hw_e(np))) * sin(2*pi*(t-delay)/tau_e(np)) * ( sin(arg)*polarization_vec_e(np,1,:) - cos(arg)*polarization_vec_e(np,2,:) )

          ! Cos-squared pulse:
          E_t(:) = E_t(:) * hE_0(np)
          ! Gaussian pulse:
          ! E_t = E_t * (-hE_0) * exp(-2*log(2*(t-delay/tau_e)**2))
        end if
      end do pulses_e
    else
      ! For the electric field cos(wt), vector potential is sin(wt)/w
      E_t(:) = hE_0(1)*(cos(hw_e(1)*t)*polarization_vec_e(1,1,:)+sin(hw_e(1)*t)*polarization_vec_e(1,2,:)) 
    end if
  end subroutine electric_field



  subroutine save_state(rank,s,dimH,nkpt,rtime,step,eval_kn,evec_kn,mx_t,my_t,mz_t)
  !! This subroutine is to Save current time-propagation state
    use mod_kind,               only: dp, int32, int64
    use mod_parameters,         only: output
    use mod_tools,              only: itos
    use mod_logging,            only: log_warning
    use mod_system,             only: System_type
    implicit none
    integer(int32), intent(in) :: rank
    !! Rank ID of process
    type(System_type), intent(in)   :: s
    !! System derived type containing number of orbitals and number of atoms
    integer(int32), intent(in) :: dimH
    !! Dimension of the hamiltonian and eigenvalues/eigenvectors
    integer(int64), intent(in) :: nkpt
    !! Number of k-points
    real(dp),                               intent(in) :: rtime
    !! Last calculated time instant
    real(dp),                               intent(in) :: step 
    !! Last step-size
    real(dp),    dimension(dimH,nkpt),      intent(in) :: eval_kn
    !! Eigenvalues for all k-points
    complex(dp), dimension(dimH,dimH,nkpt), intent(in) :: evec_kn
    !! Eigenvectors (in columns) for all k-points
    real(dp),    dimension(s%nOrb,s%nAtoms),    intent(in) :: mx_t,my_t,mz_t
    !! Orbital- and site-dependent magnetization components

    ! Local variables:
    integer(int32) :: file_unit = 6101
    !! File unit
    character(len=500) :: output_file
    !! Filename
    character(len=100) :: formatvar
    !! Format variable
    integer(int32) :: i,j,mu
    integer(int64) :: k

    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/checkpoint',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),rank,trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=file_unit, file=output_file, status='replace', form='formatted') 

    call write_header_time_prop(file_unit, "# dim = " // trim(itos(dimH)) // " , nkpt = " // trim(itos(nkpt)) )

    ! Writing iteration, time and step size
    write(unit=file_unit,fmt="(2(es16.9,2x))") rtime, step

    ! Writing eigenvalues and eigenvectors
    write(formatvar,fmt="(a,i0,a)") '(',dimH,'(es16.8e3,2x))'
    do k=1,nkpt
      write(unit=file_unit,fmt=formatvar) (eval_kn(j,k),j=1,dimH)
      do i=1,dimH
        write(unit=file_unit,fmt=formatvar) (evec_kn(i,j,k),j=1,dimH)
      end do
    end do

    ! Writing current orbital-dependent magnetization vector to use for total torque calculation (dM/dt) when recovering
    write(formatvar,fmt="(a,i0,a)") '(',3*s%nOrb*s%nAtoms,'(es16.8e3,2x))'
    write(unit=file_unit,fmt=formatvar) ((mx_t(mu,i),my_t(mu,i),mz_t(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)

    close(file_unit) 

    call log_warning("save_state", "Checkpoint saved successfully.")

  end subroutine save_state


  function recover_state(rank,s,dimH,nkpt,rtime,step,eval_kn,evec_kn,mx_t,my_t,mz_t) result(success)
    !! This subroutine is to recover current time-propagation state
    !! Note that the recover has to use the same MPI setup as used when saving the state, since different files per rank are written
    use mod_kind,               only: dp,int32,int64
    use mod_logging,            only: log_warning
    use mod_parameters,         only: output
    use mod_system,             only: System_type
    implicit none
    integer(int32), intent(in) :: rank
    !! Rank ID of process
    type(System_type), intent(in)   :: s
    !! System derived type containing number of orbitals and number of atoms
    integer(int32), intent(in) :: dimH
    !! Dimension of the hamiltonian and eigenvalues/eigenvectors
    integer(int64), intent(in) :: nkpt
    !! Number of k-points
    real(dp),                               intent(out) :: rtime
    !! Last calculated time instant
    real(dp),                               intent(out) :: step 
    !! Last step-size
    real(dp),    dimension(dimH,nkpt),      intent(out) :: eval_kn
    !! Eigenvalues for all k-points
    complex(dp), dimension(dimH,dimH,nkpt), intent(out) :: evec_kn
    !! Eigenvectors (in columns) for all k-points
    real(dp),    dimension(s%nOrb,s%nAtoms),    intent(out) :: mx_t,my_t,mz_t
    !! Orbital- and site-dependent magnetization components
    logical :: success
    !! Indication of when a state was recovered or not

    ! Local variables:
    integer(int32) :: file_unit = 6102
    !! File unit
    character(len=500) :: output_file
    !! Filename
    character(len=100) :: formatvar
    !! Format variable
    logical :: success_header
    !! Indication of when a state was recovered or not
    integer(int32) :: i,j,mu,stat_check,stat_temp
    integer(int64) :: k

    success = .false.

    ! Try to recover current state (Read file)
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/checkpoint',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),rank,trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=file_unit, file=output_file, status='old', action='read', iostat=stat_check) 
    if (stat_check /= 0) then
      call log_warning("recover_state", "No checkpoint file found.")
      return
    end if

    ! Read header and check if it is the same
    call check_header_time_prop(file_unit,success_header)
    if (.not.success_header) then
      call log_warning("recover_state", "Checkpoint file found, but header differs. Cannot continue from previous point.")
      close(file_unit)
      return
    end if

    ! Reading iteration, time and step size
    read(unit=file_unit,fmt=*, iostat=stat_temp) rtime, step
    if (stat_check /= 0) then
      call log_warning("recover_state", "Checkpoint file with correct header found, but there was a problem reading previous time and step. Cannot continue from previous point.")
      close(file_unit)
      return
    end if

    ! Reading eigenvalues and eigenvectors
    write(formatvar,fmt="(a,i0,a)") '(',dimH,'(es16.8e3,2x))'
    stat_check = 0
    do k=1,nkpt
      read(unit=file_unit,fmt=formatvar, iostat=stat_temp) (eval_kn(j,k),j=1,dimH)
      stat_check = stat_check + stat_temp
      do i=1,dimH
        read(unit=file_unit,fmt=formatvar, iostat=stat_temp) (evec_kn(i,j,k),j=1,dimH)
        stat_check = stat_check + stat_temp
      end do
    end do

    ! Reading current orbital-dependent magnetization vector to use for total torque calculation (dM/dt) when recovering
    write(formatvar,fmt="(a,i0,a)") '(',3*s%nOrb*s%nAtoms,'(es16.8e3,2x))'
    read(unit=file_unit,fmt=formatvar, iostat=stat_temp) ((mx_t(mu,i),my_t(mu,i),mz_t(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)
    stat_check = stat_check + stat_temp

    if (stat_check /= 0) then
      call log_warning("recover_state", "Checkpoint file with correct header found, but there was a problem reading eigenvalues and eigenvectors. Cannot continue from previous point.")
      close(file_unit)
      return
    end if

    ! Read header and check if it is the same
    success = .true.
    call log_warning("recover_state", "Checkpoint recovered.")

    close(file_unit) 


  end function recover_state

  subroutine write_header_time_prop(unit,title_line)
  !! Writing header for previously opened file of unit "unit"
    use mod_kind,            only: int32
    implicit none
    integer(int32),   intent(in) :: unit
    character(len=*), intent(in) :: title_line
    integer(int32) :: i,j

    if(lmagnetic) then
      if(lpulse_m) then
        write(unit=unit, fmt="('# npulse_m = ',i0)") npulse_m
        do i = 1,npulse_m
          write(unit=unit,fmt="('#    pol_m = ( ',es16.9,',',es16.9,',',es16.9,' ) in-phase    ')") (polarization_vec_m(i,1,j), j=1,3)
          write(unit=unit,fmt="('#          = ( ',es16.9,',',es16.9,',',es16.9,' ) out-of-phase')") (polarization_vec_m(i,2,j), j=1,3)
          write(unit=unit,fmt="('#    hw1_m = ',es16.9)") hw1_m(i)
          write(unit=unit,fmt="('#     hw_m = ',es16.9)") hw_m(i)
          write(unit=unit,fmt="('#    tau_m = ',es16.9)") tau_m(i)
          write(unit=unit,fmt="('#  delay_m = ',es16.9)") delay_m(i)
        end do
      else
        write(unit=unit,fmt=  "('#    pol_m = ( ',f6.3,',',f6.3,',',f6.3,' ) in-phase    ')") (polarization_vec_m(1,1,j), j=1,3)
        write(unit=unit,fmt=  "('#          = ( ',f6.3,',',f6.3,',',f6.3,' ) out-of-phase')") (polarization_vec_m(1,2,j), j=1,3)
        write(unit=unit,fmt=  "('#    hw1_m = ',es16.9)") hw1_m(1)
        write(unit=unit,fmt=  "('#     hw_m = ',es16.9)") hw_m(1)
      end if
    end if
    if(lelectric) then
      if(lpulse_e) then
        write(unit=unit,fmt="('# npulse_e = ',i0)") npulse_e
        do i = 1,npulse_e
          write(unit=unit,fmt="('#    pol_e = ( ',es16.9,',',es16.9,',',es16.9,' ) in-phase    ')") (polarization_vec_e(i,1,j), j=1,3)
          write(unit=unit,fmt="('#          = ( ',es16.9,',',es16.9,',',es16.9,' ) out-of-phase')") (polarization_vec_e(i,2,j), j=1,3)
          write(unit=unit,fmt="('#     hE_0 = ',es16.9)") hE_0(i)
          write(unit=unit,fmt="('#     hw_e = ',es16.9)") hw_e(i)
          write(unit=unit,fmt="('#    tau_e = ',es16.9)") tau_e(i)
          write(unit=unit,fmt="('#  delay_e = ',es16.9)") delay_e(i)
        end do
      else
        write(unit=unit,fmt=  "('#    pol_e = ( ',es16.9,',',es16.9,',',es16.9,' ) in-phase    ')") (polarization_vec_e(1,1,j), j=1,3)
        write(unit=unit,fmt=  "('#          = ( ',es16.9,',',es16.9,',',es16.9,' ) out-of-phase')") (polarization_vec_e(1,2,j), j=1,3)
        write(unit=unit,fmt=  "('#     hE_0 = ',es16.9)") hE_0
        write(unit=unit,fmt=  "('#     hw_e = ',es16.9)") hw_e
      end if
    end if

    write(unit=unit, fmt="(a)") title_line

  end subroutine write_header_time_prop


  subroutine check_header_time_prop(unit,success)
  !! Writing header for previously opened file of unit "unit"
    use mod_kind,            only: dp,int32,int64
    use mod_parameters,      only: dimH
    use mod_BrillouinZone,   only: realBZ
    use mod_logging,         only: log_warning
    use mod_tools,           only: itos, rtos
    implicit none
    integer(int32),   intent(in) :: unit
    logical,          intent(out) :: success
    integer(int32) :: i,j

    character(len=100) :: line
    character(len=20), dimension(7) :: title_line
    integer(int32) :: tnpulse_e,tnpulse_m,tdimH
    integer(int64) :: tnkpt
    real(dp)       :: tol=1.e-8_dp
    real(dp)       :: thE_0,thw_e,tpolarization_vec_e(2,3),ttau_e,tdelay_e,thw1_m,thw_m,tpolarization_vec_m(2,3),ttau_m,tdelay_m

    success = .true.

    if(lmagnetic) then
      if(lpulse_m) then
        read(unit=unit, fmt="(a)") line
        read(unit=line(13:), fmt=*) tnpulse_m
        if(tnpulse_m /= npulse_m) then
          call log_warning("check_header_time_prop", "npulse_m from input: " // trim(itos(npulse_m)) // ", npulse_m read: " // trim(itos(tnpulse_m)) ) 
          success = .false.
        end if

        do i = 1,tnpulse_m
          ! In-phase
          read(unit=unit, fmt="(a)") line
          read(unit=line(15:66), fmt=*) (tpolarization_vec_m(1,j), j=1,3)
          ! Out-of-phase
          read(unit=unit, fmt="(a)") line
          read(unit=line(15:66), fmt=*) (tpolarization_vec_m(2,j), j=1,3)
          if(sum(abs(tpolarization_vec_m(:,:)-polarization_vec_m(i,:,:)))  > tol) then
            call log_warning("check_header_time_prop", "polarization_vec_m(" // trim(itos(i)) // ") from input differs from checkpoint.") 
            success = .false.
          end if

          ! Intensity
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) thw1_m
          if(abs(thw1_m-hw1_m(i)) > tol) then
            call log_warning("check_header_time_prop", "hw1_m(" // trim(itos(i)) // ") from input: " // trim(rtos(hw1_m(i),"(es16.9)")) // ", read: " // trim(rtos(thw1_m,"(es16.9)")) // ".") 
            success = .false.
          end if

          ! Frequency
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) thw_m
          if(abs(thw_m-hw_m(i)) > tol) then
            call log_warning("check_header_time_prop", "hw_m(" // trim(itos(i)) // ") from input: " // trim(rtos(hw_m(i),"(es16.9)")) // ", read: " // trim(rtos(thw_m,"(es16.9)")) // ".") 
            success = .false.
          end if

          ! Linewidth
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) ttau_m
          if(abs(ttau_m-tau_m(i)) > tol) then
            call log_warning("check_header_time_prop", "tau_m(" // trim(itos(i)) // ") from input: " // trim(rtos(tau_m(i),"(es16.9)")) // ", read: " // trim(rtos(ttau_m,"(es16.9)")) // ".") 
            success = .false.
          end if

          ! Delay
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) tdelay_m
          if(abs(tdelay_m-delay_m(i)) > tol) then
            call log_warning("check_header_time_prop", "delay_m(" // trim(itos(i)) // ") from input: " // trim(rtos(delay_m(i),"(es16.9)")) // ", read: " // trim(rtos(tdelay_m,"(es16.9)")) // ".") 
            success = .false.
          end if
        end do
      else ! .not. lpulse_m
        ! In-phase
        read(unit=unit, fmt="(a)") line
        read(unit=line(15:66), fmt=*) (tpolarization_vec_m(1,j), j=1,3)
        ! Out-of-phase
        read(unit=unit, fmt="(a)") line
        read(unit=line(15:66), fmt=*) (tpolarization_vec_m(2,j), j=1,3)
        if(sum(abs(tpolarization_vec_m(:,:)-polarization_vec_m(1,:,:)))  > tol) then
          call log_warning("check_header_time_prop", "polarization_vec_m from input differs from checkpoint.") 
          success = .false.
        end if

        ! Intensity
        read(unit=unit, fmt="(a)") line
        read(unit=line(13:), fmt=*) thw1_m
        if(abs(thw1_m-hw1_m(1)) > tol) then
          call log_warning("check_header_time_prop", "hw1_m from input: " // trim(rtos(hw1_m(1),"(es16.9)")) // ", read: " // trim(rtos(thw1_m,"(es16.9)")) // ".") 
          success = .false.
        end if

        ! Frequency
        read(unit=unit, fmt="(a)") line
        read(unit=line(13:), fmt=*) thw_m
        if(abs(thw_m-hw_m(1)) > tol) then
          call log_warning("check_header_time_prop", "hw_m from input: " // trim(rtos(hw_m(1),"(es16.9)")) // ", read: " // trim(rtos(thw_m,"(es16.9)")) // ".") 
          success = .false.
        end if

      end if
    end if
    if(lelectric) then
      if(lpulse_e) then
        read(unit=unit, fmt="(a)") line
        read(unit=line(13:), fmt=*) tnpulse_e
        if(tnpulse_e /= npulse_e) then
          call log_warning("check_header_time_prop", "npulse_e from input: " // trim(itos(npulse_e)) // ", read: " // trim(itos(tnpulse_e)) // ".") 
          success = .false.
        end if

        do i = 1,tnpulse_e
          ! In-phase
          read(unit=unit, fmt="(a)") line
          read(unit=line(15:66), fmt=*) (tpolarization_vec_e(1,j), j=1,3)
          ! Out-of-phase
          read(unit=unit, fmt="(a)") line
          read(unit=line(15:66), fmt=*) (tpolarization_vec_e(2,j), j=1,3)
          if(sum(abs(tpolarization_vec_e(:,:)-polarization_vec_e(i,:,:)))  > tol) then
            call log_warning("check_header_time_prop", "polarization_vec_e(" // trim(itos(i)) // ") from input differs from checkpoint.") 
            success = .false.
          end if

          ! Intensity
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) thE_0
          if(abs(thE_0-hE_0(i)) > tol) then
            call log_warning("check_header_time_prop", "hE_0(" // trim(itos(i)) // ") from input: " // trim(rtos(hE_0(i),"(es16.9)")) // ", read: " // trim(rtos(thE_0,"(es16.9)")) // ".") 
            success = .false.
          end if

          ! Frequency
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) thw_e
          if(abs(thw_e-hw_e(i)) > tol) then
            call log_warning("check_header_time_prop", "hw_e(" // trim(itos(i)) // ") from input: " // trim(rtos(hw_e(i),"(es16.9)")) // ", read: " // trim(rtos(thw_e,"(es16.9)")) // ".") 
            success = .false.
          end if

          ! Linewidth
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) ttau_e
          if(abs(ttau_e-tau_e(i)) > tol) then
            call log_warning("check_header_time_prop", "tau_e(" // trim(itos(i)) // ") from input: " // trim(rtos(tau_e(i),"(es16.9)")) // ", read: " // trim(rtos(ttau_e,"(es16.9)")) // ".") 
            success = .false.
          end if

          ! Delay
          read(unit=unit, fmt="(a)") line
          read(unit=line(13:), fmt=*) tdelay_e
          if(abs(tdelay_e-delay_e(i)) > tol) then
            call log_warning("check_header_time_prop", "delay_e(" // trim(itos(i)) // ") from input: " // trim(rtos(delay_e(i),"(es16.9)")) // ", read: " // trim(rtos(tdelay_e,"(es16.9)")) // ".") 
            success = .false.
          end if
        end do
      else ! .not. lpulse_e
        ! In-phase
        read(unit=unit, fmt="(a)") line
        read(unit=line(15:66), fmt=*) (tpolarization_vec_e(1,j), j=1,3)
        ! Out-of-phase
        read(unit=unit, fmt="(a)") line
        read(unit=line(15:66), fmt=*) (tpolarization_vec_e(2,j), j=1,3)
        if(sum(abs(tpolarization_vec_e(:,:)-polarization_vec_e(1,:,:)))  > tol) then
          call log_warning("check_header_time_prop", "polarization_vec_e from input differs from checkpoint.") 
          success = .false.
        end if

        ! Intensity
        read(unit=unit, fmt="(a)") line
        read(unit=line(13:), fmt=*) thE_0
        if(abs(thE_0-hE_0(1)) > tol) then
          call log_warning("check_header_time_prop", "hE_0 from input: " // trim(rtos(hE_0(1),"(es16.9)")) // ", read: " // trim(rtos(thE_0,"(es16.9)")) // ".") 
          success = .false.
        end if

        ! Frequency
        read(unit=unit, fmt="(a)") line
        read(unit=line(13:), fmt=*) thw_e
        if(abs(thw_e-hw_e(1)) > tol) then
          call log_warning("check_header_time_prop", "hw_e from input: " // trim(rtos(hw_e(1),"(es16.9)")) // ", read: " // trim(rtos(thw_e,"(es16.9)")) // ".") 
          success = .false.
        end if

      end if
    end if

    read(unit=unit, fmt="(a)") line
    read(unit=line, fmt=*) title_line
    ! Dimension of the hamiltonian
    read(unit=title_line(4), fmt=*) tdimH
    if(tdimH /= dimH) then
      call log_warning("check_header_time_prop", "dimH from input: " // trim(itos(dimH)) // ", read: " // trim(itos(tdimH)) // ".") 
      success = .false.
    end if

    ! Number of k-points
    read(unit=title_line(7), fmt=*) tnkpt
    if(tnkpt /= realBZ%workload) then
      call log_warning("check_header_time_prop", "realBZ%workload from input: " // trim(itos(realBZ%workload)) // ", read: " // trim(itos(tnkpt)) // ".") 
      success = .false.
    end if

  end subroutine check_header_time_prop


  subroutine define_time_prop_observables()
  !! This subroutine is used to define the time-propagated observables 
  !! that are written into files
    use mod_parameters, only: output
    implicit none

    if(.not.allocated(output%observable)) then 
      allocate(output%observable(7))
      output%observable(1) = "occupation"
      output%observable(2) = "magnetization"
      output%observable(3) = "AngularMomentum"
      output%observable(4) = "dMdt"
      output%observable(5) = "Tso"
      output%observable(6) = "Txc"
      output%observable(7) = "Energy"
    end if

  end subroutine define_time_prop_observables


  subroutine create_time_prop_files()
  !! subroutine to create files with names and units
  !! parameters: tau_m, tau_e, hw1, hw_m, hw1_m, hw_e, hE_0, integration_time, 
  !! logical parameters: lmagnetic, lpulse_m, lelectric, lpulse_e 
  !! observables: <1>, <sigma>, <L>, <tau>, currents
    use mod_parameters, only: output,lprintfieldonly
    implicit none
    character(len=500) :: output_file
    integer :: i,unit

    call define_time_prop_observables()

    ! Time-dependent fields
    unit = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open(unit=unit,file=trim(output_file), status= 'replace')
    call write_header_time_prop(unit,'#    Time [ps]   ,    field_m x    ,    field_m y    ,    field_m z    ,    field_e x    ,    field_e y    ,    field_e z    ' )
    close(unit)

    if(lprintfieldonly) return

    do i=1,size(output%observable)
      ! for all orbitals
      unit = 5090+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i)),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open(unit=unit,file=trim(output_file), status= 'replace')
      call write_header_time_prop(unit,'#    Time [ps]   , ' // output%observable(i))
      close(unit)

      ! for separate orbitals (only energy is not orbital-dependent))
      if(i<=(size(output%observable)-1)) then
        unit = 5190+i
        write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i))//"_orb",trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open(unit=unit,file=trim(output_file), status= 'replace')
        call write_header_time_prop(unit,'#    Time [ps]   , ' // trim(output%observable(i)) // '_orb')
        close(unit)
      end if
    end do

  end subroutine create_time_prop_files

      
  subroutine open_time_prop_files()
    !! subroutine to open time propagation output files
    use mod_parameters, only: output,missing_files
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=500) :: output_file
    integer :: i,iw,err,errt=0

    call define_time_prop_observables()

    do i=1,size(output%observable)
      ! for all orbitals
      iw = 5090+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i)),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)

      ! for separate orbitals
      if(i<=(size(output%observable)-1)) then
        iw = 5190+i
        write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i))//"_orb",trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)
      end if
    end do

    iw = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
    errt = errt + err
    if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)

    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[open_time_prop_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

  end subroutine open_time_prop_files

  subroutine close_time_prop_files()
    !! subroutine to close time propagation output files
    use mod_parameters, only: output
    implicit none
    integer :: i

    do i=1,size(output%observable)
      ! for all orbitals
      close(5090+i)

      ! for separate orbitals
      if(i<=(size(output%observable)-1)) close(5190+i)
    end do

    close(6090)

    deallocate( output%observable )

  end subroutine close_time_prop_files


  subroutine write_field()
    !! subroutine to write field output files
    use mod_kind,             only: dp
    use mod_parameters,       only: output,missing_files
    use mod_mpi_pars,         only: abortProgram
    implicit none
    character(len=500)     :: output_file
    real(dp), dimension(3) :: field_m, field_e
    real(dp)               :: t, rtime
    integer :: i,iw,err,errt=0

    write(output%unit_loop,"('[write_field] Writing field... ')",advance="no")

    iw = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
    errt = errt + err
    if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)

    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[write_field] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

    t = 0._dp
    t_loop_field: do while (t <= integration_time)
      t = t + step

      ! Calculating time-dependent field
      field_m = 0._dp
      if (lmagnetic) call magnetic_field(t,field_m)
      field_e = 0._dp 
      if (lelectric) call vector_potential(t, field_e)

      rtime = t*time_conv

      write(unit=6090,fmt="(7(es16.9,2x))") rtime, (field_m(i),i=1,3), (field_e(i),i=1,3)
    end do t_loop_field

    close(6090)

    write(output%unit_loop,"('done!')")

  end subroutine write_field


  subroutine write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,field_m,field_e,E_t, Lxm_t,Lym_t,Lzm_t,tso_t,txc_t,dmdt)
    !! subroutine to write in time propagation output files
    use mod_kind,             only: dp
    use mod_system,           only: System_type
    implicit none
    type(System_type),                      intent(in) :: s
    real(dp),                               intent(in) :: t,E_t
    real(dp), dimension(s%nOrb,s%nAtoms),   intent(in) :: rho_t,mx_t,my_t,mz_t
    real(dp), dimension(2,s%nAtoms),        intent(in) :: Lxm_t,Lym_t,Lzm_t
    real(dp), dimension(3,2,s%nAtoms),      intent(in) :: tso_t
    real(dp), dimension(3,s%nOrb,s%nAtoms), intent(in) :: txc_t
    real(dp), dimension(3,s%nOrb,s%nAtoms), intent(in) :: dmdt
    real(dp), dimension(3),                 intent(in) :: field_m,field_e

    integer  :: iw,i,mu,xyz
    real(dp) :: rtime

    call open_time_prop_files()
    
    rtime = t*time_conv

    !---------------========================= Total quantities =========================---------------
    iw = 5090
    ! Site-dependent occupation
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, (sum(rho_t( 1:s%Types(s%Basis(i)%Material)%nOrb,i )),i=1,s%nAtoms)
    ! Site-dependent magnetization
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, (sum(mx_t( 1:s%Types(s%Basis(i)%Material)%nOrb,i )),sum(my_t( 1:s%Types(s%Basis(i)%Material)%nOrb,i )),sum(mz_t( 1:s%Types(s%Basis(i)%Material)%nOrb,i )), sqrt(sum(mx_t( 1:s%Types(s%Basis(i)%Material)%nOrb,i ))**2 + sum(my_t( 1:s%Types(s%Basis(i)%Material)%nOrb,i ))**2 + sum(mz_t( 1:s%Types(s%Basis(i)%Material)%nOrb,i ))**2) ,i=1,s%nAtoms)
    ! Site-dependent OAM
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, (sum(Lxm_t( 1:2,i )), sum(Lym_t(1:2,i )), sum(Lzm_t( 1:2,i )), sqrt(sum(Lxm_t( 1:2,i ))**2 + sum(Lym_t( 1:2,i ))**2 + sum(Lzm_t( 1:2,i ))**2),i=1,s%nAtoms)

    ! dM/dt
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( (sum(dmdt( xyz,1:s%Types(s%Basis(i)%Material)%nOrb,i )), xyz=1,3), i=1,s%nAtoms)
    ! SO-torque
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( (sum(tso_t( xyz,1:2,i )), xyz=1,3), i=1,s%nAtoms)
    ! XC-torque
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( (sum(txc_t( xyz,1:s%Types(s%Basis(i)%Material)%nOrb,i )), xyz=1,3), i=1,s%nAtoms)

    ! Energy
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, E_t

    !----------======================== Orbital-dependent quantities ========================----------
    iw = 5190
    ! Orbital dependent occupation:
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ((rho_t(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)
    ! Orbital dependent magnetization:
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( (mx_t(mu,i),my_t(mu,i),mz_t(mu,i), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)
    ! Orbital dependent OAM, for p (1) and d (2) orbitals
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( (Lxm_t(mu,i), Lym_t(mu,i), Lzm_t(mu,i), mu=1,2), i=1,s%nAtoms)

    ! dM/dt
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( ((dmdt(xyz,mu,i), xyz=1,3), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)
    ! SO-torque
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( ((tso_t(xyz,mu,i), xyz=1,3), mu=1,2), i=1,s%nAtoms)
    ! XC-torque
    iw = iw+1
    write(unit=iw,fmt="(100(es16.9,2x))") rtime, ( ((txc_t(xyz,mu,i), xyz=1,3), mu=1,s%Types(s%Basis(i)%Material)%nOrb), i=1,s%nAtoms)

    !----------------========================= Applied Fields =========================----------------
    ! Time-dependent field
    write(unit=6090,fmt="(7(es16.9,2x))") rtime, (field_m(i),i=1,3), (field_e(i),i=1,3)

    call close_time_prop_files()

  end subroutine write_time_prop_files
end module mod_time_propagator

