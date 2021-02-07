module mod_time_propagator_io
  implicit none
contains

  ! Writing header for previously opened file of unit "unit"
  subroutine write_header_time_prop(unit,title_line)
    use mod_kind,             only: int32
    use mod_imRK4_parameters, only: lelectric, hE_0, hw_e, lpulse_e, npulse_e, polarization_vec_e, tau_e, delay_e, lmagnetic, hw1_m, hw_m, lpulse_m, npulse_m, polarization_vec_m, tau_m, delay_m
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


  ! Writing header for previously opened file of unit "unit"
  subroutine check_header_time_prop(unit,success)
    use mod_kind,             only: dp,int32,int64
    use mod_parameters,       only: dimH
    use mod_BrillouinZone,    only: realBZ
    use mod_io,               only: log_warning
    use mod_tools,            only: itos, rtos
    use mod_imRK4_parameters, only: lelectric, hE_0, hw_e, lpulse_e, npulse_e, polarization_vec_e, tau_e, delay_e, lmagnetic, hw1_m, hw_m, lpulse_m, npulse_m, polarization_vec_m, tau_m, delay_m
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

  ! subroutine to create files with names and units
  ! parameters: tau_m, tau_e, hw1, hw_m, hw1_m, hw_e, hE_0, integration_time, 
  ! logical parameters: lmagnetic, lpulse_m, lelectric, lpulse_e 
  ! observables: <1>, <sigma>, <L>, <tau>, currents
  subroutine create_time_prop_files()
    use mod_parameters, only: output,lprintfieldonly
    implicit none
    character(len=500) :: output_file
    integer :: i,unit
     
    allocate(output%observable(4))
    output%observable(1) = "occupation"
    output%observable(2) = "magnetization"
    output%observable(3) = "Energy"
    output%observable(4) = "AngularMomentum"

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

      ! for separate orbitals
      if(i<=3) then
        unit = 5190+i
        write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i))//"_orb",trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open(unit=unit,file=trim(output_file), status= 'replace')
        call write_header_time_prop(unit,'#    Time [ps]   , ' // trim(output%observable(i)) // '_orb')
        close(unit)
      end if
    end do

  end subroutine create_time_prop_files

      
  ! subroutine to open time propagation output files
  subroutine open_time_prop_files()
    use mod_parameters, only: output,missing_files
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=500) :: output_file
    integer :: i,iw,err,errt=0

    if(.not.allocated(output%observable)) then 
      allocate(output%observable(4))
      output%observable(1) = "occupation"
      output%observable(2) = "magnetization"
      output%observable(3) = "Energy"
      output%observable(4) = "AngularMomentum"
    end if

    do i=1,size(output%observable)
      ! for all orbitals
      iw = 5090+i
      write(output_file,"('./results/',a1,'SOC/',a,'/time_propagation/',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%observable(i)),trim(output%time_field),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=iw, file=trim(output_file), status='old', position='append',form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(output_file)

      ! for separate orbitals
      if(i<=3) then
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

  ! subroutine to close time propagation output files
  subroutine close_time_prop_files()
    use mod_parameters, only: output
    implicit none
    integer :: i

    do i=1,size(output%observable)
      ! for all orbitals
      close(5090+i)

      ! for separate orbitals
      if(i<=3) close(5190+i)
    end do

    close(6090)

    deallocate( output%observable )

  end subroutine close_time_prop_files


  ! subroutine to write field output files
  subroutine write_field()
    use mod_kind,             only: dp
    use mod_parameters,       only: output,missing_files
    use mod_mpi_pars,         only: abortProgram
    use mod_imRK4_parameters, only: step, time_conv, integration_time, lelectric, lmagnetic
    use mod_imRK4,            only: magnetic_field, vector_potential
    implicit none
    character(len=500)     :: output_file
    real(dp), dimension(3) :: field_m, field_e
    real(dp)               :: t, time
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

      time = t*time_conv

      write(unit=6090,fmt="(7(es16.9,2x))") time, (field_m(i),i=1,3), (field_e(i),i=1,3)
    end do t_loop_field

    close(6090)

    write(output%unit_loop,"('done!')")

  end subroutine write_field


  ! subroutine to write in time propagation output files
  subroutine write_time_prop_files(s,t,rho_t,mx_t,my_t,mz_t,field_m,field_e,E_t, Lxm_t,Lym_t,Lzm_t) 
    use mod_kind,             only: dp
    use mod_system,           only: System_type
    use mod_imRK4_parameters, only: time_conv
    implicit none
    type(System_type),                    intent(in) :: s
    real(dp),                             intent(in) :: t,E_t
    real(dp), dimension(s%nOrb,s%nAtoms), intent(in) :: rho_t,mx_t,my_t,mz_t
    real(dp), dimension(2,s%nAtoms),      intent(in) :: Lxm_t,Lym_t,Lzm_t
    real(dp), dimension(3),               intent(in) :: field_m, field_e

    integer  :: i, mu
    real(dp) :: time

    call open_time_prop_files()
    
    time = t*time_conv

    write(unit=5091,fmt="(100(es16.9,2x))") time, (sum(rho_t(:,i)),i=1,s%nAtoms)
    write(unit=5092,fmt="(100(es16.9,2x))") time, (sum(mx_t(:,i)),sum(my_t(:,i)),sum(mz_t(:,i)), sqrt(sum(mx_t(:,i))**2 + sum(my_t(:,i))**2 + sum(mz_t(:,i))**2) ,i=1,s%nAtoms)
    write(unit=5093,fmt="(100(es16.9,2x))") time, E_t
    write(unit=5094,fmt="(100(es16.9,2x))") time, (sum(Lxm_t(:,i)), sum(Lym_t(:,i)), sum(Lzm_t(:,i)), sqrt(sum(Lxm_t(:,i))**2 + sum(Lym_t(:,i))**2 + sum(Lzm_t(:,i))**2),i=1,s%nAtoms)

    ! Orbital dependent occupation:
    write(unit=5191,fmt="(100(es16.9,2x))") time, ((rho_t(mu,i), mu=1,s%nOrb), i=1,s%nAtoms)

    ! Orbital dependent magnetization:        
    write(unit=5192,fmt="(100(es16.9,2x))") time, ( (mx_t(mu,i),my_t(mu,i),mz_t(mu,i), mu=1,s%nOrb), i=1,s%nAtoms)

    ! Orbital dependent OAM, for p (1) and d (2) orbitals
    write(unit=5193,fmt="(100(es16.9,2x))") time, ( (Lxm_t(mu,i), Lym_t(mu,i), Lzm_t(mu,i), mu=1,2), i=1,s%nAtoms)

    write(unit=6090,fmt="(7(es16.9,2x))") time, (field_m(i),i=1,3), (field_e(i),i=1,3)

    call close_time_prop_files()

  end subroutine write_time_prop_files
end module mod_time_propagator_io

