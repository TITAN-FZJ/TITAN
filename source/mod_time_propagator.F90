module mod_time_propagator
  implicit none
contains

  subroutine time_propagator(s)
    use mod_f90_kind,         only: double
    use mod_constants,        only: cZero, cI
    use mod_imRK4_parameters, only: dimH2, step, integration_time, ERR, safe_factor, lelectric, hE_0, hw_e, lpulse_e, tau_e, delay_e, lmagnetic, hw1_m, hw_m, lpulse_m, tau_m, delay_m
    use mod_RK_matrices,      only: A, id, id2, M1, build_identity
    use mod_imRK4,            only: iterate_Zki, calculate_step_error, magnetic_pulse_B, evec_potent, hext_t
    use mod_BrillouinZone,    only: realBZ
    use mod_parameters,       only: nOrb,dimH,output,laddresults,lprintfieldonly
    use mod_system,           only: System
    use ElectricField,        only: EshiftBZ,ElectricFieldVector
    use mod_expectation,      only: expec_val_n, expec_H_n, expec_L_n
    use mod_Umatrix,          only: update_Umatrix
    use mod_magnet,           only: rhod0,rho0
    use mod_tools,            only: KronProd
    use mod_mpi_pars
    implicit none
    type(System), intent(in)   :: s

    integer                                         :: i, it, n, counter
    real(double)                                    :: t, p, h_new, h_old, ERR_old, ERR_kn
    complex(double), dimension(dimH)                :: Yn, Yn_hat, Yn_new, Yn_e, Yn_hat_e, Yn_new_e
    complex(double), dimension(:,:,:), allocatable  :: evec_kn,evec_kn_temp
    real(double),    dimension(:,:),   allocatable  :: eval_kn

    integer                                         :: iz, info 
    integer                                         :: lwork
    real(double)                                    :: weight, kp(3)
    real(double),    dimension(nOrb,s%nAtoms)       :: expec_0, expec_z 
    complex(double), dimension(nOrb,s%nAtoms)       :: expec_p
    real(double),    dimension(:),      allocatable :: rwork(:), eval(:)
    complex(double),                    allocatable :: work(:), hk(:,:)

    real(double),    dimension(nOrb,s%nAtoms)   :: rho_t, mx_t, my_t, mz_t
    real(double),    dimension(s%nAtoms)        :: lxm, lxpm, lym, lypm, lzm, lzpm
    real(double),    dimension(s%nAtoms)        :: lxm_t, lxpm_t, lym_t, lypm_t, lzm_t, lzpm_t
    real(double),    dimension(s%nAtoms)        :: lxmd_t, lxpmd_t, lymd_t, lypmd_t, lzmd_t, lzpmd_t
    complex(double), dimension(nOrb,s%nAtoms)   :: mp_t
    real(double),    dimension(s%nAtoms)        :: rhod_t, mxd_t, myd_t, mzd_t
    complex(double), dimension(s%nAtoms)        :: mpd_t
    real(double),    dimension(3)               :: field_m, field_e, A_t
    real(double)                                :: E_t, E_0
   
    if(rFreq(1) == 0) then
      write(output%unit_loop,"('CALCULATING TIME-PROPAGATION')")
      ! Creating files and writing headers
      if(.not.laddresults) then
        call create_time_prop_files()
      end if
    end if
 
    if(lprintfieldonly) then
      call write_field()
      return
    end if

    ! working space for the eigenstate solver
    lwork = 21*dimH

    ! Dimensions for RK method
    dimH2  = 2*dimH

    allocate( id(dimH,dimH),id2(dimH2,dimH2),hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork),eval_kn(dimH,realBZ%workload),evec_kn(dimH,dimH,realBZ%workload),evec_kn_temp(dimH,dimH,realBZ%workload), M1(dimH2,dimH2) )

    ! Building identities
    call build_identity(dimH,id)
    call build_identity(size(A,1)*dimH,id2)

    ! Building matrix M1
    call KronProd(size(A,1),size(A,1),dimH,dimH,A,id,M1)

    ! Time propagation over t, kpoints, eigenvectors(Yn) for each k
    t = 0.d0
    it = 0  
    t_loop: do while (t <= integration_time)
      t = t + step
      if(rField==0) &
        write(output%unit_loop,"('[time_propagator] Time: ',es10.3,' of ',es10.3)") t, integration_time
       
      counter = 0 ! Counter for the calculation of error in the step size for each time t
      ! Propagation loop for a given time t, calculating the optimal step size
      do
        !$omp parallel default(none) &
        !$omp& private(Yn_e,Yn_new_e,Yn_hat_e,iz,n,i,kp,weight,hk,ERR_kn,Yn, Yn_new, Yn_hat, eval,work,rwork,info,expec_0,expec_p,expec_z,E_0,lxm,lxpm,lym,lypm,lzm,lzpm) &
        !$omp& shared( ERR,counter,step,s,t,it,dimH,realBZ,evec_kn_temp,evec_kn,eval_kn,rho_t,mp_t,mz_t,E_t, Lxm_t,Lxpm_t,Lym_t,Lypm_t,Lzm_t,Lzpm_t,lwork,EshiftBZ,ElectricFieldVector)

        rho_t  = 0.d0
        mp_t   = cZero
        mz_t   = 0.d0

        E_t    = 0.d0

        Lxm_t  = 0.d0
        Lxpm_t = 0.d0
        Lym_t  = 0.d0
        Lypm_t = 0.d0
        Lzm_t  = 0.d0
        Lzpm_t = 0.d0

        ERR    = 0.d0
  
        !$omp do reduction(+:rho_t,mp_t,mz_t,ERR)
        kpoints_loop: do iz = 1, realBZ%workload
          kp = realBZ%kp(1:3,iz) + EshiftBZ*ElectricFieldVector
          weight = realBZ%w(iz)   
          if (it==0) then
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
            if (it==0) then
              Yn(:)= hk(:,n)
    !----------------------------------------------------------------------------------------!
    !-------- These steps are to control the error by solving for c^~ instead of c ----------!
    !----------------------------------------------------------------------------------------!
      ! 1- find Yn^~ = Yn * e^(i*En*t/hbar) >>>> for each n and k point (only once).
      ! 2- propagate Yn^~ using H^~ = H(t) - En.
      ! 3- find Yn = Yn^~ * exp( (-i*En*t/hbar) ) at each step.
      ! 4- calculate the error at each step from Yn but propagate with Yn^~.
    !----------------------------------------------------------------------------------------!
    !----------------------------------------------------------------------------------------!
              ! Geting the intitial vector Yn^~, setting Yn to Yn^~ for propagation, hbar = 1
              ! Yn = Yn * exp(cI*eval_kn(n,iz)*t) 
                Yn = Yn * exp(cI*eval_kn(n,iz)*t)
                
            else
              ! Is the propagated vector Yn^~
              Yn(:) = evec_kn(:,n,iz)
               
            end if

            call iterate_Zki(s,t,kp,eval_kn(n,iz),step,Yn_new,Yn,Yn_hat)
            ! Note: all iteration outputs corresponds to Yn^~ 

            ! Getting Yn's again; Yn = Yn^~ * exp( (-i*En*t/hbar) ), hbar=1, rename Yn to Yn_e
            Yn_e(:)     = Yn     * exp(-cI*eval_kn(n,iz)*t) 
            Yn_new_e(:) = Yn_new * exp(-cI*eval_kn(n,iz)*t) 
            Yn_hat_e(:) = Yn_hat * exp(-cI*eval_kn(n,iz)*t) 


            ! Calculating expectation values for the nth eigenvector
            ! Note: use the Yn_new_e or Yn_nw >>> should give the same result

            ! calculating expectation value for magnetization in eigenvector (n)
            call expec_val_n(s, dimH, Yn_new_e, eval_kn(n,iz), expec_0, expec_p, expec_z)

            ! calculating expectation value for the T.D Hamiltonian in eigenvector (n)
            call expec_H_n(s, kp, t, dimH, Yn_new_e, eval_kn(n,iz), E_0)

            ! calculating expectation value of angular momentum in eigenvector (n)
            call expec_L_n(s, dimH, Yn_new_e, eval_kn(n,iz), lxm, lxpm, lym, lypm, lzm, lzpm)


            rho_t = rho_t + expec_0  * weight 
            mp_t  = mp_t  + expec_p  * weight 
            mz_t  = mz_t  + expec_z  * weight 

            E_t = E_t + E_0 * weight

            Lxm_t  = Lxm_t   + lxm  * weight
            Lxpm_t = Lxpm_t  + lxpm * weight
            Lym_t  = Lym_t   + lym  * weight
            Lypm_t = Lypm_t  + lypm * weight
            Lzm_t  = Lzm_t   + lzm  * weight
            Lzpm_t = Lzpm_t  + lzpm * weight
          


            ! Calculation of the error and the new step size.
            ! Note: use Yn_e, Yn_new_e, Yn_hat_e.
            call calculate_step_error(Yn_e,Yn_new_e,Yn_hat_e,ERR_kn)
write(*,*) iz,kp,ERR_kn

            ERR = ERR + ERR_kn * weight

            ! Storing temporary propagated vector before checking if it's Accepted.
            ! Note: save the Yn^~ outputs to be propagated.
            evec_kn_temp(:,n,iz) = Yn_new(:)

          end do evs_loop
        end do kpoints_loop
        !$omp end do
        !$omp end parallel
write(*,*) 'ERR=',ERR
stop
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
        if ( ERR > 1.d0) then
          ! repeat the calculation using h_new
          t = t - step + h_new
          h_old = step
          step = h_new
          write(*,*) "Rejected", t, step, ERR
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

      write(*,*)  "Accepted", t, step, ERR 

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
      field_m = 0.d0
      if (lmagnetic) then 
        if (lpulse_m) then
          if ((t >= delay_m).and.(t <= 8.d0*tau_m+delay_m)) then

            ! ************** update here ******************

            call magnetic_pulse_B(t,field_m)
          end if  
        else
          field_m = [ hw1_m*cos(hw_m*t), hw1_m*sin(hw_m*t), 0.d0 ]
        end if 
      end if 

      field_e = 0.d0 
      if (lelectric) then 
        if (lpulse_e) then
          if ((t >= delay_e).and.(t <= tau_e+delay_e)) then
            call evec_potent(t, A_t)
            field_e = A_t
          end if  
        else
          field_e = [ hE_0*cos(hw_e*t), hE_0*sin(hw_e*t), 0.d0 ]
        end if 
      end if 

      call write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,rhod_t, mxd_t, myd_t, mzd_t, field_m, field_e, E_t, lxm_t, lxpm_t, lym_t, lypm_t, lzm_t, lzpm_t) 

      counter = counter + 1
      it = it + 1

    end do t_loop

  end subroutine time_propagator


  ! Writing header for previously opened file of unit "unit"
  subroutine write_header_time_prop(unit,title_line)
    use mod_f90_kind,         only: double
    use mod_imRK4_parameters, only: lelectric, hE_0, hw_e, lpulse_e, tau_e, delay_e, lmagnetic, hw1_m, hw_m, lpulse_m, tau_m, delay_m
    implicit none
    integer,          intent(in)           :: unit
    character(len=*), intent(in)           :: title_line

    if(lmagnetic) then
      write(unit=unit, fmt="('#   hw1_m = ',es16.9)") hw1_m
      write(unit=unit, fmt="('#    hw_m = ',es16.9)") hw_m
      if(lpulse_m) then
        write(unit=unit, fmt="('#   tau_m = ',es16.9)") tau_m
        write(unit=unit, fmt="('# delay_m = ',es16.9)") delay_m
      end if
    end if
    if(lelectric) then
      write(unit=unit, fmt="('#   hE_0 = ',es16.9)") hE_0
      write(unit=unit, fmt="('#    hw_e = ',es16.9)") hw_e
      if(lpulse_e) then
        write(unit=unit, fmt="('#   tau_e = ',es16.9)") tau_e
        write(unit=unit, fmt="('# delay_e = ',es16.9)") delay_e
      end if
    end if

    write(unit=unit, fmt="(a)") title_line

  end subroutine write_header_time_prop

  ! subroutine to create files with names and units
  ! parameters: tau_m, tau_e, hw1, hw_m, hw1_m, hw_e, hE_0, integration_time, 
  ! logical parameters: lmagnetic, lpulse_m, lelectric, lpulse_e 
  ! observables: <1>, <sigma>, <L>, <tau>, currents
  subroutine create_time_prop_files()
    use mod_parameters,       only: output,lprintfieldonly
    use mod_imRK4_parameters, only: lelectric, lpulse_e, lmagnetic, lpulse_m
    implicit none
    character(len=500) :: output_file
    integer :: i,unit
     
    allocate(output%observable(4))
    output%observable(1) = "occupation"
    output%observable(2) = "magnetization"
    output%observable(3) = "Energy"
    output%observable(4) = "AngularMomentum"

    if(lmagnetic) then
      output%time_field = "_magfield"
      if(lpulse_m) then
        output%time_field = trim(output%time_field) // "pulse"
      end if
    end if
    if(lelectric) then
      output%time_field = trim(output%time_field) // "_efield"
      if(lpulse_e) then
        output%time_field = trim(output%time_field) // "pulse"
      end if
    end if

    ! Time-dependent fields
    unit = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open(unit=unit,file=trim(output_file), status= 'replace')
    call write_header_time_prop(unit,'#    Time [ps]   ,    field_m x   ,    field_m y   ,    field_m z   ,    field_e x   ,    field_e y   ,    field_e z   ' )
    close(unit)

    if(lprintfieldonly) return

    do i=1,size(output%observable)
      ! for all orbitals
      unit = 5090+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i)),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open(unit=unit,file=trim(output_file), status= 'replace')
      call write_header_time_prop(unit,'#    Time [ps]   , ' // output%observable(i))
      close(unit)

      ! for d orbitals
      unit = 5190+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i))//"_d",trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open(unit=unit,file=trim(output_file), status= 'replace')
      call write_header_time_prop(unit,'#    Time [ps]   , ' // output%observable(i) // '_d')
      close(unit)
    end do

  end subroutine create_time_prop_files

      
  ! subroutine to open time propagation output files
  ! is this important or can be inside write_time_prop_files ??? it is also needed inside write
  subroutine open_time_prop_files()
    use mod_parameters, only: output,missing_files
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=500) :: output_file
    integer :: i,iw,err,errt=0

    do i=1,size(output%observable)
      ! for all orbitals
      iw = 5090+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i)),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)

      ! for d orbitals
      iw = 5190+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i))//"_d",trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)
    end do

    iw = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
    errt = errt + err
    if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)

    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[open_time_prop_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

  end subroutine open_time_prop_files

  ! subroutine to close time propagation output files
  subroutine close_time_prop_files()
    use mod_parameters, only: output
    implicit none
    integer :: i

    do i=1,size(output%observable)
      ! for all orbitals
      close(5090+i)

      ! for d orbitals
      close(5190+i)
    end do

    close(6090)

  end subroutine close_time_prop_files


  ! subroutine to close time propagation output files
  subroutine write_field()
    use mod_f90_kind,         only: double
    use mod_parameters,       only: output,missing_files
    use mod_mpi_pars,         only: rFreq,abortProgram
    use mod_imRK4_parameters, only: step, integration_time, lelectric, lpulse_e, delay_e, tau_e, hE_0, hw_e, lmagnetic, lpulse_m, delay_m, tau_m, hw1_m, hw_m
    use mod_imRK4,            only: magnetic_pulse_B, evec_potent
    implicit none
    character(len=500)         :: output_file
    real(double), dimension(3) :: field_m, field_e, A_t
    real(double)               :: t, time
    integer :: i,iw,err,errt=0

    if(rFreq(1) == 0) &
      write(output%unit_loop,"('[write_field] Writing field... ')",advance="no")

    iw = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
    errt = errt + err
    if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)

    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[write_field] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

    t = 0.d0
    t_loop_field: do while (t <= integration_time)
      t = t + step
  
      time = t*4.84d-5

      ! Calculating time-dependent field      
      field_m = 0.d0
      if (lmagnetic) then 
        if (lpulse_m) then
          if ((t >= delay_m).and.(t <= 8.d0*tau_m+delay_m)) then
            call magnetic_pulse_B(t,field_m)
          end if  
        else
          field_m = [ hw1_m*cos(hw_m*t), hw1_m*sin(hw_m*t), 0.d0 ]
        end if 
      end if 

      field_e = 0.d0 
      if (lelectric) then 
        if (lpulse_e) then
          if ((t >= delay_e).and.(t <= tau_e+delay_e)) then
            call evec_potent(t,A_t)
            field_e = A_t
          end if  
        else
          field_e = [ hE_0*cos(hw_e*t), hE_0*sin(hw_e*t), 0.d0 ]
        end if 
      end if 

      write(unit=6090,fmt="(7(es16.9,2x))") time, (field_m(i),i=1,3), (field_e(i),i=1,3)
    end do t_loop_field

    close(6090)

    if(rFreq(1) == 0) &
      write(output%unit_loop,"('done!')")

  end subroutine write_field


  ! subroutine to write in time propagation output files
  subroutine write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,rhod_t, mxd_t, myd_t, mzd_t, field_m, field_e, E_t, Lxm_t, Lxpm_t, Lym_t, Lypm_t, Lzm_t, Lzpm_t) 
    use mod_f90_kind,   only: double
    use mod_system,     only: System
    use mod_parameters, only: nOrb
    implicit none
    type(System),                          intent(in) :: s
    real(double),                          intent(in) :: t, E_t
    real(double), dimension(nOrb,s%nAtoms),intent(in) :: rho_t, mx_t, my_t, mz_t
    real(double), dimension(s%nAtoms)     ,intent(in) :: rhod_t, mxd_t, myd_t, mzd_t, Lxm_t, Lxpm_t, Lym_t, Lypm_t, Lzm_t, Lzpm_t 
    real(double), dimension(3)            ,intent(in) :: field_m, field_e

    integer      :: i
    real(double) :: time

    call open_time_prop_files()
    
    time = t*4.84d-5

    write(unit=5091,fmt="(100(es16.9,2x))") time, (sum(rho_t(:,i)),i=1,s%nAtoms)
    write(unit=5092,fmt="(100(es16.9,2x))") time, (sum(mx_t(:,i)),sum(my_t(:,i)),sum(mz_t(:,i)), sqrt(sum(mx_t(:,i))**2 + sum(my_t(:,i))**2 + sum(mz_t(:,i))**2) ,i=1,s%nAtoms)

    write(unit=5191,fmt="(100(es16.9,2x))") time, (rhod_t(i),i=1,s%nAtoms)
    write(unit=5192,fmt="(100(es16.9,2x))") time, (mxd_t(i),myd_t(i),mzd_t(i),sqrt(mxd_t(i)**2 + myd_t(i)**2 + mzd_t(i)**2),i=1,s%nAtoms)

    write(unit=6090,fmt="(7(es16.9,2x))") time, (field_m(i),i=1,3), (field_e(i),i=1,3)

    write(unit=5093,fmt="(100(es16.9,2x))") time, E_t

    ! write(unit=5094, fmt="(100(es16.9,2x))") time,  sum(Lxm_t(i)), sum(Lym_t(i)), sum(Lzm_t(i)), sum(Lxpm_t(i)),sum(Lypm_t(i)), sum(Lzpm_t(i)), sqrt(sum(Lxm_t(i))**2 + sum(Lym_t(i))**2 + sum(Lzm_t(i))**2 )  
    write(unit=5094, fmt="(100(es16.9,2x))") time, sum(Lxm_t), sum(Lym_t), sum(Lzm_t), sum(Lxpm_t), sum(Lypm_t), sum(Lzpm_t), sqrt(Lxm_t**2 + Lym_t**2 + Lzm_t**2)
   
    call close_time_prop_files()

  end subroutine write_time_prop_files
end module mod_time_propagator

