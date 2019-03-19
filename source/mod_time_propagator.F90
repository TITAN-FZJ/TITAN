module mod_time_propagator
  implicit none
contains

  subroutine time_propagator(s)
    use mod_f90_kind,         only: double
    use mod_imRK4_parameters, only: dimH2, time, initialize, step
    use mod_RK_matrices,      only: A, build_identity
    use mod_imRK4,            only: iterate_Zki
    use mod_BrillouinZone,    only: realBZ
    use mod_parameters,       only: output,dimH
    use mod_system,           only: System
    use TightBinding,         only: nOrb,nOrb2
    use ElectricField,        only: EshiftBZ,ElectricFieldVector
    use mod_expectation,      only: expec_val_n
    use mod_Umatrix,          only: update_Umatrix
    implicit none
    character(len=6)                  :: output_file = "output"
    integer                           :: i, j
    real(double)                      :: t
    complex(double), dimension(dimH)  :: Yn
    complex(double), dimension(:,:,:), allocatable  :: evec_kn

    type(System),                                   intent(in)  :: s
    ! real(double),         dimension(nOrb,s%nAtoms), intent(out) :: rho_t, mx_t, my_t, mz_t
    ! complex(double), dimension(nOrb,s%nAtoms),      intent(out) :: mp_t

    integer                                         :: iz, info 
    integer                                         :: lwork
    real(double)                                    :: weight, kp(3)
    real(double),    dimension(nOrb,s%nAtoms)       :: expec_0, expec_z
    complex(double), dimension(nOrb,s%nAtoms)       :: expec_p
    real(double),    dimension(:),  allocatable     :: rwork(:), eval(:)
    complex(double),                    allocatable :: work(:), hk(:,:)

    ! Open file and write headers ******** use un-used unit number ********
    open(unit=11,file=trim(output_file), status= 'replace')
    write(unit=11,fmt=*) '#      Time      ', '        M_x       ', '        M_y       ', '        M_z       ', '      M      '
    

    ! Building identity
    call build_identity(size(A,1)*dimH) 

    ! Obtaining number of steps for the time loop
    time= int(integration_time/step)

    ! Dimensions for RK method
    dimH2  = 2*dimH

    ! working space for the eigenstate solver
    lwork = 21*dimH

    allocate( hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork),evec_kn(dimH,dimH,realBZ%workload) )

    ! Time propagation over t, kpoints, eigenvectors(Yn) for each k
    t_loop: do i = 0, time
      t = i*step
      ! rho_t = 0.d0
      ! mp_t  = zero
      ! mz_t  = 0.d0

      kpoints_loop:    do iz = 1, realBZ%workload
        kp = realBZ%kp(1:3,iz) !+ EshiftBZ*ElectricFieldVector
        weight = realBZ%w(iz)   
        if (i==0) then
          ! Calculating the hamiltonian for a given k-point
          call hamiltk(s,kp,hk)
          ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
          call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)
          ! Save the initial eigenvectors in evec_i matrix 
        end if
        evs_loop: do j=1, dimH
          if (i==0) then
            Yn(:)= hk(:,j)
          else
            Yn(:)= evec_kn(:,j,iz)
          end if
          call iterate_Zki(s, t, kp, Yn)
          ! Calculating expectation values for the nth eigenvector 
          !call expec_val_n(s, dimH, Yn, eval, expec_0, expec_p, expec_z)
          ! TODO: Orbital dependent variables?
          ! rho_t = rho_t + expec_0 * weight 
          ! mp_t  = mp_t + expec_p * weight 
          ! mz_t  = mz_t + expec_z * weight 
          
          ! Update the eignvectors array of dimension ( dimH * dimH , k )
          evec_kn(:,j,iz)= Yn
        end do evs_loop
      end do kpoints_loop
      ! Print outputs
      ! rho_t = rho_t + rho_kn
      ! mp_t  = mp_t  + mp_kn
      ! mz_t  = mz_t  + mz_kn
      ! mx = real(mp)
      ! my = aimag(mp)

      ! ! TODO: who are the orbital dependent variables?
      ! call update_Umatrix(mzd_in,mpd_in,rhod_in,rhod0,rho_in,rho0,s%nAtoms,nOrb)

      ! write(unit=11,fmt="(6(es16.9,2x))") t, mx, my, mz, (mx**2 + my**2 +mz**2)  
    end do t_loop
    close(unit=11)
  end subroutine time_propagator
end module mod_time_propagator


