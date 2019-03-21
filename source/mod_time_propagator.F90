module mod_time_propagator
  implicit none
contains

  subroutine time_propagator(s)
    use mod_f90_kind,         only: double
    use mod_constants,        only: cZero
    use mod_imRK4_parameters, only: dimH2, time, step, integration_time
    use mod_RK_matrices,      only: A, id, id2, M1, build_identity
    use mod_imRK4,            only: iterate_Zki
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

    character(len=9)                  :: output_file = "time_prop"
    integer                           :: i, it, n, m, mu
    real(double)                      :: t
    complex(double), dimension(dimH)  :: Yn
    complex(double), dimension(:,:,:), allocatable  :: evec_kn

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

    ! Open file and write headers
    open(unit=5090,file=trim(output_file), status= 'replace')
    write(unit=5090,fmt=*) '#      Time      ', '        M_x       ', '        M_y       ', '        M_z       ', '      M      '
    
    ! Obtaining number of steps for the time loop
    time = int(integration_time/step)

    ! working space for the eigenstate solver
    lwork = 21*dimH

    ! Dimensions for RK method
    dimH2  = 2*dimH

    allocate( id(dimH,dimH),id2(dimH2,dimH2),hk(dimH,dimH),rwork(3*dimH-2),eval(dimH),work(lwork),evec_kn(dimH,dimH,realBZ%workload), M1(dimH2,dimH2) )

    ! Building identities
    call build_identity(dimH,id)
    call build_identity(size(A,1)*dimH,id2)

    ! Building matrix M1
    call KronProd(size(A,1),size(A,1),dimH,dimH,A,id,M1)

    ! Time propagation over t, kpoints, eigenvectors(Yn) for each k
    t_loop: do it = 0, time
      t = it*step

      !$omp parallel default(none) &
      !$omp& private(iz,n,m,i,mu,kp,weight,hk,Yn,eval,work,rwork,info,expec_0,expec_p,expec_z) &
      !$omp& shared(s,t,it,dimH,realBZ,evec_kn,rho_t,mp_t,mz_t,lwork,EshiftBZ,ElectricFieldVector)

      rho_t = 0.d0
      mp_t  = cZero
      mz_t  = 0.d0

      !$omp do reduction(+:rho_t,mp_t,mz_t)
      kpoints_loop: do iz = 1, realBZ%workload
        kp = realBZ%kp(1:3,iz) + EshiftBZ*ElectricFieldVector
        weight = realBZ%w(iz)   
        if (it==0) then
          ! Calculating the hamiltonian for a given k-point
          call hamiltk(s,kp,hk)
          ! Diagonalizing the hamiltonian to obtain eigenvectors and eigenvalues
          call zheev('V','L',dimH,hk,dimH,eval,work,lwork,rwork,info)
        end if
        evs_loop: do n = 1, dimH
          if (it==0) then
            Yn(:)= hk(:,n)
          else
            Yn(:)= evec_kn(:,n,iz)
          end if
do m=1,dimH
if(abs(Yn(m))>1.d-6) write(100,*) it,t,iz,n,Yn(m)
end do
          call iterate_Zki(s, t, kp, eval(n), Yn)
do m=1,dimH
if(abs(Yn(m))>1.d-6) write(101,*) it,t,iz,n,Yn(m)
end do
          ! Calculating expectation values for the nth eigenvector 
          call expec_val_n(s, dimH, Yn, eval(n), expec_0, expec_p, expec_z)
do i=1,s%nAtoms
do mu=1,9
if(abs(expec_0(mu,i))>1.d-6)  write(102,*) it,t,iz,n,i,mu,expec_0(mu,i)
if(abs(expec_p(mu,i))>1.d-6)  write(103,*) it,t,iz,n,i,mu,expec_p(mu,i)
if(abs(expec_z(mu,i))>1.d-6)  write(104,*) it,t,iz,n,i,mu,expec_z(mu,i)
end do
end do
          rho_t = rho_t + expec_0 * weight 
          mp_t  = mp_t  + expec_p * weight 
          mz_t  = mz_t  + expec_z * weight 
          
          ! Update the eignvectors array of dimension ( dimH * dimH , k )
          evec_kn(:,n,iz)= Yn
        end do evs_loop
      end do kpoints_loop
      !$omp end do
      !$omp end parallel

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

      write(unit=5090,fmt="(100(es16.9,2x))") t, (sum(mx_t(:,i)), sum(my_t(:,i)), sum(mz_t(:,i)), sqrt(sum(mx_t(:,i))**2 + sum(my_t(:,i))**2 + sum(mz_t(:,i))**2), i=1,s%nAtoms)
    end do t_loop
    close(unit=5090)
  end subroutine time_propagator
end module mod_time_propagator


