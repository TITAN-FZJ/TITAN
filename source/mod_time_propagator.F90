module mod_time_propagator
  implicit none
contains

  subroutine time_propagator(s)
    use mod_f90_kind,         only: double
    use mod_constants,        only: cZero
    use mod_imRK4_parameters, only: dimH2, time, step, integration_time, ERR, Delta, lelectric, hw1_e, hw_e, lpulse_e, tau_e, delay_e, lmagnetic, hw1_m, hw_m, lpulse_m, tau_m, delay_m
    use mod_RK_matrices,      only: A, id, id2, M1, build_identity
    use mod_imRK4,            only: iterate_Zki, calculate_step_error, magnetic_pulse_B, electric_pulse_e
    use mod_BrillouinZone,    only: realBZ
    use mod_parameters,       only: dimH
    use mod_system,           only: System
    use TightBinding,         only: nOrb
    use ElectricField,        only: EshiftBZ,ElectricFieldVector
    use mod_expectation,      only: expec_val_n
    use mod_Umatrix,          only: update_Umatrix
    use mod_magnet,           only: rhod0,rho0
    use mod_tools,            only: KronProd
    implicit none
    type(System),         intent(in)  :: s

    integer                           :: i, it, n, counter
    real(double)                      :: t, p, h_new, h_old, ERR_old, ERR_kn
    complex(double), dimension(dimH)  :: Yn, Yn_hat, Yn_new
    complex(double), dimension(:,:,:), allocatable  :: evec_kn,evec_kn_temp
    real(double),    dimension(:,:),   allocatable  :: eval_kn

    integer                                     :: iz, info 
    integer                                     :: lwork
    real(double)                                :: weight, kp(3)
    real(double),    dimension(nOrb,s%nAtoms)   :: expec_0, expec_z
    complex(double), dimension(nOrb,s%nAtoms)   :: expec_p
    real(double),    dimension(:),  allocatable :: rwork(:), eval(:)
    complex(double),                allocatable :: work(:), hk(:,:)

    real(double),    dimension(nOrb,s%nAtoms)   :: rho_t, mx_t, my_t, mz_t
    complex(double), dimension(nOrb,s%nAtoms)   :: mp_t
    real(double),    dimension(s%nAtoms)        :: rhod_t, mxd_t, myd_t, mzd_t
    complex(double), dimension(s%nAtoms)        :: mpd_t
    real(double), dimension(3)                  :: field_m, field_e
    
    ! Creating files and writing headers
    call create_time_prop_files()

    ! Obtaining number of steps for the time loop
    time = int(integration_time/step)

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
      counter = 0 ! Counter for the calculation of error in the step size for each time t

      ! Propagation loop for a given time t, calculating the optimal step size
      do
        !$omp parallel default(none) &
        !$omp& private(iz,n,i,kp,weight,hk,ERR_kn,Yn, Yn_new, Yn_hat, eval,work,rwork,info,expec_0,expec_p,expec_z) &
        !$omp& shared(ERR, counter, step, s,t,it,dimH,realBZ,evec_kn_temp,evec_kn,eval_kn,rho_t,mp_t,mz_t,lwork,EshiftBZ,ElectricFieldVector)

        rho_t = 0.d0
        mp_t  = cZero
        mz_t  = 0.d0
        ERR   = 0.d0

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
          end if
          evs_loop: do n = 1, dimH
            if (it==0) then
              Yn(:)= hk(:,n)
            else
              Yn(:)= evec_kn(:,n,iz)
            end if

            call iterate_Zki(s,t,kp,eval_kn(n,iz),step,Yn_new,Yn,Yn_hat)
             
            ! Calculating expectation values for the nth eigenvector 
            call expec_val_n(s, dimH, Yn_new, eval_kn(n,iz), expec_0, expec_p, expec_z)

            rho_t = rho_t + expec_0 * weight 
            mp_t  = mp_t  + expec_p * weight 
            mz_t  = mz_t  + expec_z * weight 

            ! Calculation of the error and the new step size
            call calculate_step_error(Yn,Yn_new,Yn_hat,ERR_kn)
                                      
            ERR = ERR + ERR_kn * weight

            ! Storing temporary propagated vector before checking if it's Accepted
            evec_kn_temp(:,n,iz) = Yn_new(:)

          end do evs_loop
        end do kpoints_loop
        !$omp end do
        !$omp end parallel

        ERR = sqrt(ERR)
        ! Find the new step size h_new
        ! h_new = delta * h_used / (ERR)^(1/p+1) where p = 2*s, delta is some saftey factor 
        ! delta = 0.9 * (2*K_max + 1)/ (2*K_max + NEWT) where NEWT is the number of newton iterations
        p = 2.0*size(A,1) 
        ! h_new = delta * step * (ERR)**(-1.0/(2.0*s)) !! simpler formula
        if (counter /= 0) then
          h_new = delta * step * (ERR)**(-1.0/p) * (step/h_old) * (ERR_old/ERR)**(1.0/p)!! Gustafsson (1994) formula          
        else 
          h_new = delta * step * (ERR)**(-1.0/p)
        end if 
        ! Tests if the step size is good:
        ! the condition (h_new < delta* step) is equivalent to the condition (ERR > 1)
        ! if ( h_new < 0.9 * step) then
           ! do while ( h_new < step) 
        ! save ERR from iterate_Zki to ERR_old
        ERR_old = ERR
        if ( ERR > 1.d0) then
          ! repeat the calculation using h_new
          t = t - step + h_new
          h_old = step
          step = h_new
          write(*,*) t, step, ERR
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
      do i = 1, s%nAtoms
        rhod_t(i) = sum(rho_t(5:9,i))
        mpd_t(i)  = sum(mp_t(5:9,i))
        mxd_t(i)  = sum(mx_t(5:9,i))
        myd_t(i)  = sum(my_t(5:9,i))
        mzd_t(i)  = sum(mz_t(5:9,i))
      end do

      call update_Umatrix(mzd_t,mpd_t,rhod_t,rhod0,rho_t,rho0,s%nAtoms,nOrb)

      ! Calculating time-dependent field      
      if (lmagnetic) then 
        if (lpulse_m) then
          if ((t >= delay_m).and.(t <= 8.d0*tau_m+delay_m)) then
            call magnetic_pulse_B(t,field_m)
          else
            field_m(1)= 0.d0
            field_m(2)= 0.d0
            field_m(3)= 0.d0 
          end if  
        else
          field_m = [ hw1_m*cos(hw_m*t), hw1_m*sin(hw_m*t), 0.d0 ]
        end if 
      end if 
   
      if (lelectric) then 
        if (lpulse_e) then
          if ((t >= delay_e).and.(t <= 8.d0*tau_e+delay_e)) then
            call electric_pulse_e(t,field_e)
          else
            field_e(1)= 0.d0
            field_e(2)= 0.d0
            field_e(3)= 0.d0 
          end if  
        else
          field_e = [ hw1_e*cos(hw_e*t), hw1_e*sin(hw_e*t), 0.d0 ]
        end if 
      end if 

      call open_time_prop_files()
      call write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,rhod_t, mxd_t, myd_t, mzd_t,field_m,field_e)
      call close_time_prop_files()


      counter = counter + 1
      it = it + 1

    end do t_loop

  end subroutine time_propagator



  ! Writing header for previously opened file of unit "unit"
  subroutine write_header_time_prop(unit,title_line)
    use mod_f90_kind,         only: double
    use mod_imRK4_parameters, only: lelectric, hw1_e, hw_e, lpulse_e, tau_e, delay_e, lmagnetic, hw1_m, hw_m, lpulse_m, tau_m, delay_m

    integer,          intent(in)           :: unit
    character(len=*), intent(in)           :: title_line
    integer :: i

    if(lmagnetic) then
      write(unit=unit, fmt="('#.  hw1_m = ',es16.9)") hw1_m
      write(unit=unit, fmt="('#    hw_m = ',es16.9)") hw_m
      if(lpulse_m) then
        write(unit=unit, fmt="('#   tau_m = ',es16.9)") tau_m
        write(unit=unit, fmt="('# delay_m = ',es16.9)") delay_m
      end if
    end if
    if(lelectric) then
      write(unit=unit, fmt="('#.  hw1_e = ',es16.9)") hw1_e
      write(unit=unit, fmt="('#    hw_e = ',es16.9)") hw_e
      if(lpulse_e) then
        write(unit=unit, fmt="('#   tau_e = ',es16.9)") tau_e
        write(unit=unit, fmt="('# delay_e = ',es16.9)") delay_e
      end if
    end if

    write(unit=unit, fmt="(a)") title_line

  end subroutine write_header_time_prop

  ! subroutine to create files with names and units
  ! parameters: tau_m, tau_e, hw1, hw_m, hw1_m, hw_e, hw1_e, integration_time, 
  ! logical parameters: lmagnetic, lpulse_m, lelectric, lpulse_e 
  ! observables: <1>, <sigma>, <L>, <tau>, currents
  subroutine create_time_prop_files()
    use mod_parameters,       only: output
    use mod_imRK4_parameters, only: lelectric, lpulse_e, lmagnetic, lpulse_m

    character(len=500) :: output_file
    integer :: i,unit
     
    allocate(output%observable(2))
    output%observable(1) = "occupation"
    output%observable(2) = "magnetization"

    if(lmagnetic) then
      output%time_field = "_magfield"
      if(lpulse_m) then
        output%time_field = output%time_field // "pulse"
      end if
    end if
    if(lelectric) then
      output%time_field = output%time_field // "_efield"
      if(lpulse_e) then
        output%time_field = output%time_field // "pulse"
      end if
    end if

    do i=1,size(output%observable)
      ! for all orbitals
      unit = 5090+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i)),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open(unit=unit,file=trim(output_file), status= 'replace')
      call write_header_time_prop(unit,'#      Time      , ' // output%observable(i))
      close(unit)

      ! for d orbitals
      unit = 5190+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i))//"_d",trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open(unit=unit,file=trim(output_file), status= 'replace')
      call write_header_time_prop(unit,'#      Time      , ' // output%observable(i) // '_d')
      close(unit)
    end do

    unit = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open(unit=unit,file=trim(output_file), status= 'replace')
    call write_header_time_prop(unit,'#      Time      ,      field     ' )
    close(unit)

  end subroutine create_time_prop_files

      
  ! subroutine to open time propagation output files
  ! is this important or can be inside write_time_prop_files ??? it is also needed inside write
  subroutine open_time_prop_files()
    use mod_parameters, only: output

    character(len=500) :: output_file
    integer :: i,unit

    do i=1,size(output%observable)
      ! for all orbitals
      unit = 5090+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i)),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open(unit=unit,file=trim(output_file), status= 'old')

      ! for d orbitals
      unit = 5190+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i))//"_d",trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open(unit=unit,file=trim(output_file), status= 'old')
    end do

    unit = 6090
    write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/field',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open(unit=unit,file=trim(output_file), status= 'old')

  end subroutine open_time_prop_files

  ! subroutine to close time propagation output files
  subroutine close_time_prop_files()
    use mod_parameters, only: output
    integer :: i

    do i=1,size(output%observable)
      ! for all orbitals
      close(5090+i)

      ! for d orbitals
      close(5190+i)
    end do

    close(6090)

  end subroutine close_time_prop_files

  ! subroutine to write in time propagation output files
  subroutine write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,rhod_t, mxd_t, myd_t, mzd_t,field_m,field_e)
    use mod_f90_kind, only: double
    use mod_system,   only: System
    use TightBinding, only: nOrb

    type(System),                          intent(in) :: s
    real(double),                          intent(in) :: t
    real(double), dimension(nOrb,s%nAtoms),intent(in) :: rho_t, mx_t, my_t, mz_t
    real(double), dimension(s%nAtoms)     ,intent(in) :: rhod_t, mxd_t, myd_t, mzd_t
    real(double), dimension(3)            ,intent(in) :: field_m, field_e

    integer      :: i
    real(double) :: time

    time = t*6.582d-7
    write(unit=5091,fmt="(100(es16.9,2x))") time, (sum(rho_t(:,i)),i=1,s%nAtoms)
    write(unit=5092,fmt="(100(es16.9,2x))") time, (sum(mx_t(:,i)),sum(my_t(:,i)),sum(mz_t(:,i)), sqrt(sum(mx_t(:,i))**2 + sum(my_t(:,i))**2 + sum(mz_t(:,i))**2) ,i=1,s%nAtoms)

    write(unit=5191,fmt="(100(es16.9,2x))") time, (rhod_t(i),i=1,s%nAtoms)
    write(unit=5192,fmt="(100(es16.9,2x))") time, (mxd_t(i),myd_t(i),mzd_t(i),sqrt(mxd_t(i)**2 + myd_t(i)**2 + mzd_t(i)**2),i=1,s%nAtoms)

    write(unit=5191,fmt="(100(es16.9,2x))") time, (field_m(i),i=1,3), (field_e(i),i=1,3)
  end subroutine write_time_prop_files
end module mod_time_propagator

