module mod_time_propagator
  implicit none
contains

  subroutine time_propagator(s)
    use mod_kind, only: dp, int32, int64
    use mod_constants,         only: cZero, cI
    use mod_imRK4_parameters,  only: dimH2, step, integration_time, ERR, safe_factor, lelectric, &
                                     hE_0, hw_e, lpulse_e, tau_e, delay_e, lmagnetic, hw1_m, hw_m, &
                                     lpulse_m, tau_m, delay_m, polarization_e, polarization_vec_e, &
                                     polarization_m, polarization_vec_m
    use mod_RK_matrices,       only: A, id, id2, M1, build_identity
    use mod_imRK4,             only: iterate_Zki, calculate_step_error, magnetic_field, vector_potential
    use mod_BrillouinZone,     only: realBZ
    use mod_parameters,        only: nOrb,dimH,output,laddresults,lprintfieldonly
    use mod_system,            only: System_type
    use mod_expectation,       only: expec_val_n, expec_H_n, expec_L_n
    use mod_Umatrix,           only: update_Umatrix
    use mod_magnet,            only: rhod0,rho0
    use mod_tools,             only: KronProd
    use mod_superconductivity, only: allocate_supercond_variables
    use mod_hamiltonian,       only: hamiltk,hamilt_local
    use mod_checkpoint,        only: save_state, recover_state
    use mod_io,                only: log_warning
    use mod_time_propagator_io, only: create_time_prop_files, write_header_time_prop, write_field, write_time_prop_files
    use mod_mpi_pars
    implicit none
    type(System_type), intent(in)   :: s

    integer(int32)                              :: i, it, n, counter, iter_rej, iter_tot
    real(dp)                                    :: t, p, h_new, h_old, ERR_old, ERR_kn
    complex(dp), dimension(dimH)                :: Yn, Yn_hat, Yn_new, Yn_e, Yn_hat_e, Yn_new_e
    complex(dp), dimension(:,:,:), allocatable  :: evec_kn,evec_kn_temp
    real(dp),    dimension(:,:),   allocatable  :: eval_kn

    logical  :: use_checkpoint
    real(dp) :: t_cp,step_cp

    integer(int64)                              :: iz
    integer(int32)                              :: info, ncount
    integer(int32)                              :: lwork
    real(dp)                                :: weight, kp(3)
    real(dp),    dimension(nOrb,s%nAtoms)   :: expec_0, expec_z
    complex(dp), dimension(nOrb,s%nAtoms)   :: expec_p
    real(dp),    dimension(nOrb,s%nAtoms)   :: expec_singlet
    real(dp),    dimension(:),  allocatable :: rwork(:), eval(:)
    complex(dp),                allocatable :: work(:), hk(:,:)

    real(dp),    dimension(nOrb,s%nAtoms)   :: rho_t, mx_t, my_t, mz_t
    complex(dp), dimension(nOrb,s%nAtoms)   :: mp_t
    real(dp),    dimension(s%nAtoms)        :: rhod_t, mxd_t, myd_t, mzd_t
    complex(dp), dimension(s%nAtoms)        :: mpd_t
    real(dp),    dimension(nOrb,s%nAtoms)   :: singlet_coupling_t
    real(dp),    dimension(s%nAtoms)        :: lxm, lym, lzm
    real(dp),    dimension(s%nAtoms)        :: lxm_t, lym_t, lzm_t
    real(dp),    dimension(3)               :: field_m, field_e
    real(dp)                                :: E_t, E_0
   
    if(rFreq(1) == 0) then
      write(output%unit_loop,"('CALCULATING TIME-PROPAGATION')")
      ! Creating files and writing headers
      if(.not.laddresults) then
        call create_time_prop_files()
      end if
    end if
 
    if(lprintfieldonly) then
      if(rFreq(1) == 0) &
        call write_field()
      return
    end if

    ! number of elements in the MPI communication
    ncount = nOrb*s%nAtoms

    ! working space for the eigenstate solver
    lwork = 21*dimH

    ! Dimensions for RK method
    dimH2  = 2*dimH

    allocate( id(dimH,dimH),id2(dimH2,dimH2),hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork),eval_kn(dimH,realBZ%workload),evec_kn(dimH,dimH,realBZ%workload),evec_kn_temp(dimH,dimH,realBZ%workload), M1(dimH2,dimH2) )

    ! Build local hamiltonian
    call hamilt_local(s)

    ! Building identities
    call build_identity(dimH,id)
    call build_identity(size(A,1)*dimH,id2)

    ! Building matrix M1
    call KronProd(size(A,1),size(A,1),dimH,dimH,A,id,M1)

    ! Checking for checkpoints
    use_checkpoint = recover_state(rFreq(1),dimH,realBZ%workload,t_cp,step_cp,eval_kn,evec_kn)
    if(use_checkpoint) then
      t = t_cp
      step = step_cp
    else
      t = 0._dp
    end if

    if(rFreq(1) == 0) &
      write(output%unit_loop,"('[time_propagator] Starting propagation from t = ',es9.2)") t

    it = 0       ! Counter of accepted iterations
    iter_tot = 0 ! Counter of total number of iterations (rejected + accepted)

    ! Time propagation over t, kpoints, eigenvectors(Yn) for each k
    t_loop: do while (t <= integration_time)
      t = t + step
      if(rFreq(1) == 0) &
        write(output%unit_loop,"('[time_propagator] Time: ',es10.3,' of ',es10.3)", advance='no') t, integration_time
       
      counter = 0  ! Counter for the calculation of error in the step size for each time t
      iter_rej = 0 ! Counter of rejected steps (for each accepted one)
      ! Propagation loop for a given time t, calculating the optimal step size
      do
        !$omp parallel default(none) &
        !$omp& private(Yn_e,Yn_new_e,Yn_hat_e,iz,n,i,kp,weight,hk,ERR_kn,Yn, Yn_new, Yn_hat, eval,work,rwork,info,expec_0,expec_p,expec_z,expec_singlet,E_0,lxm,lym,lzm) &
        !$omp& shared( ERR,counter,step,s,t,it,dimH,realBZ,evec_kn_temp,evec_kn,eval_kn,use_checkpoint,rho_t,mp_t,mz_t,E_t, Lxm_t,Lym_t,Lzm_t,singlet_coupling_t,lwork)

        rho_t  = 0._dp
        mp_t   = cZero
        mz_t   = 0._dp

        E_t    = 0._dp

        Lxm_t  = 0._dp
        Lym_t  = 0._dp
        Lzm_t  = 0._dp

        ERR    = 0._dp

        singlet_coupling_t = 0._dp
  
        !$omp do reduction(+:rho_t,mp_t,mz_t,E_t,Lxm_t,Lym_t,Lzm_t,singlet_coupling_t,ERR)
        kpoints_loop: do iz = 1, realBZ%workload
          kp = realBZ%kp(1:3,iz)
          weight = realBZ%w(iz)   
          if ((it==0).and.(.not.use_checkpoint)) then
            ! Calculating the hamiltonian for a given k-point
            call hamiltk(s,kp,hk)
            ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
            call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)
            eval_kn(:,iz) = eval(:)
          !>>>>> find eigenvalues again
          ! Improve step control by diagonalizing the Hamiltonian after some number of steps:
          ! calling H(t) instead of hamiltk, at t=0, H(0)= hamiltk
          end if
          evs_loop: do n = 1, dimH
            if ((it==0).and.(.not.use_checkpoint)) then
              Yn(:)= hk(:,n)
    !----------------------------------------------------------------------------------------!
    !-------- These steps are to control the error by solving for c^~ instead of c ----------!
    !----------------------------------------------------------------------------------------!
      ! 1- find Yn^~ = Yn * e^(i*En*t/hbar) : for each n and k point (only once).
      ! 2- propagate Yn^~ using H^~ = H(t) - En.
      ! 3- find Yn = Yn^~ * exp( (-i*En*t/hbar) ) at each step.
      ! 4- calculate the error at each step from Yn but propagate with Yn^~.
    !----------------------------------------------------------------------------------------!
              ! Geting the intitial vector Yn^~, setting Yn to Yn^~ for propagation, hbar = 1
              Yn = Yn * exp(cI*eval_kn(n,iz)*t)
                
            else
              ! Is the propagated vector Yn^~
              Yn(:) = evec_kn(:,n,iz)
               
            end if

            call iterate_Zki(s,t,kp,eval_kn(n,iz),step,Yn_new,Yn,Yn_hat)
            ! Note: all iteration outputs correspond to Yn^~ 

            ! Getting Yn's again; Yn = Yn^~ * exp( (-i*En*t/hbar) ), hbar=1, rename Yn to Yn_e
            Yn_e(:)     = Yn     * exp(-cI*eval_kn(n,iz)*t) 
            Yn_new_e(:) = Yn_new * exp(-cI*eval_kn(n,iz)*t) 
            Yn_hat_e(:) = Yn_hat * exp(-cI*eval_kn(n,iz)*t) 

            ! Calculating expectation values of the nth eigenvector
            ! Note: use the Yn_new_e or Yn_nw >>> should give the same result

            ! calculating expectation value of magnetization in eigenvector (n)
            call expec_val_n(s, dimH, Yn_new_e, eval_kn(n,iz), expec_0, expec_p, expec_z, expec_singlet)

            ! calculating expectation value of the T.D Hamiltonian in eigenvector (n)
            call expec_H_n(s, kp, t, dimH, Yn_new_e, eval_kn(n,iz), E_0)

            ! calculating expectation value of angular momentum in eigenvector (n)
            call expec_L_n(s, dimH, Yn_new_e, eval_kn(n,iz), lxm, lym, lzm)


            rho_t = rho_t + expec_0  * weight 
            mp_t  = mp_t  + expec_p  * weight 
            mz_t  = mz_t  + expec_z  * weight 

            ! Superconducting order parameter
            singlet_coupling_t = singlet_coupling_t + expec_singlet*weight

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
            evec_kn_temp(:,n,iz) = Yn_new(:)

          end do evs_loop
        end do kpoints_loop
        !$omp end do
        !$omp end parallel

        call MPI_Allreduce(MPI_IN_PLACE, rho_t, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, mz_t , ncount, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, mp_t , ncount, MPI_DOUBLE_COMPLEX  , MPI_SUM, FreqComm(1) , ierr)
        call MPI_Allreduce(MPI_IN_PLACE, ERR  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, FreqComm(1) , ierr)

        ERR = sqrt(ERR)
        ! Find the new step size h_new
        ! h_new = safe_factor * h_used / (ERR)^(1/p+1) where p = 2*s, safe_factor is some safety factor 
        ! safe_factor = 0.9 * (2*K_max + 1)/ (2*K_max + NEWT) where NEWT is the number of newton iterations
        p = 2.0*size(A,1) 
        ! h_new = safe_factor * step * (ERR)**(-1.0/(2.0*s)) !! simpler formula
        if (counter /= 0) then
          h_new = safe_factor * step * (ERR)**(-1.0/p) * (step/h_old) * (ERR_old/ERR)**(1.0/p)!! Gustafsson (1994) formula          
        else 
          h_new = safe_factor * step * (ERR)**(-1.0/p)
        end if 
        ! Tests if the step size is good:
        ! the condition (h_new < safe_factor * step) is equivalent to the condition (ERR > 1)
        ! if ( h_new < 0.9 * step) then
           ! do while ( h_new < step) 
        ! save ERR from iterate_Zki to ERR_old
        ERR_old = ERR

        iter_tot = iter_tot + 1 ! Counter of total number of iterations
        if ( ERR > 1._dp) then
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
          ! Update the eignvectors array of dimension ( dimH * dimH , k )
          evec_kn = evec_kn_temp
          step = h_new
          exit
        end if
          
      end do

      if(rFreq(1) == 0) &
        write(output%unit_loop,"(' (',i0,' rejected iterations)')") iter_rej

      ! if(rFreq(1) == 0) write(*,*)  "Accepted", t, step, ERR 

      mx_t = real(mp_t)
      my_t = aimag(mp_t)

      ! obtaining expectation values for d-orbitals
      do i = 1, s%nAtoms
        rhod_t(i) = sum(rho_t(5:9,i))
        mpd_t(i)  = sum(mp_t(5:9,i))
        mxd_t(i)  = sum(mx_t(5:9,i))
        myd_t(i)  = sum(my_t(5:9,i))
        mzd_t(i)  = sum(mz_t(5:9,i))
      end do
 
      call update_Umatrix(mzd_t,mpd_t,rhod_t,rhod0,rho_t,rho0,s%nAtoms,nOrb)

      ! Calculating time-dependent field
      field_m = 0._dp
      if (lmagnetic) call magnetic_field(t,field_m)
      field_e = 0._dp 
      if (lelectric) call vector_potential(t, field_e)

      if(rFreq(1) == 0) &
        call write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,rhod_t, mxd_t, myd_t, mzd_t, field_m, field_e, E_t, lxm_t, lym_t, lzm_t) 

      counter = counter + 1
      it = it + 1 ! Counter of accepted iterations

      ! Checking for "save" file to trigger checkpoint
      open(unit=911, file="save", status='old', iostat=info)
      if(info==0) then
        close(911)
        call save_state(rFreq(1),dimH,realBZ%workload,t,step,eval_kn,evec_kn)
        call MPI_Barrier(FieldComm, ierr)
        if(rFreq(1) == 0) &
          call system('rm save')
      end if
    end do t_loop

    ! Creating checkpoint in the last state
    call save_state(rFreq(1),dimH,realBZ%workload,t,step,eval_kn,evec_kn)

    deallocate( id,id2,hk,rwork,eval,work,eval_kn,evec_kn,evec_kn_temp, M1, output%observable )

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

