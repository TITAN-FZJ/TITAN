!> module for imRK4 routines
module mod_imRK4
  implicit none

contains
  ! subroutine to find the vectors Z_ki
  subroutine iterate_Zki(s,t,kp,eval,step,Yn_new,Yn,Yn_hat)
    use mod_f90_kind,         only: double
    use mod_parameters,       only: dimH
    use mod_system,           only: System
    use mod_imRK4_parameters, only: dimH2, sc_tol
    use mod_tools,            only: vec_norm, LS_solver
    use mod_constants,        only: cZero
    use mod_RK_matrices
    implicit none
    ! define the variables: 
    type(System)                     , intent(in)    :: s
    real(double)                     , intent(in)    :: t
    real(double)                     , intent(in)    :: kp(3)
    real(double)                     , intent(in)    :: eval
    complex(double), dimension(dimH) , intent(inout) :: Yn_new, Yn_hat, Yn
   
    integer                                 :: k
    real(double)                            :: step, error, norm_k_old, norm_k_new, theta_k, eta_k, tol
    complex(double), dimension(dimH2)       :: deltaZ_k_old, deltaZ_k  
    complex(double), dimension(dimH2,dimH2) :: M2n
    complex(double), dimension(dimH2)       :: Fun 
    complex(double), dimension(dimH2)       :: Z_k 

    ! z = (g - Y)
    ! z = 0.0 is the solution
    ! Initial value of z_k (Taylor expansion point) 
    Z_k   = cZero
    k     = 0
    error = 1.d0
    !
    ! Solve the linear system of equations of size 2 * n ,
    ! where n is the dimension of the hamiltonian:
    !   (I- hA\otimes J)  \Delta z^k = -z^k+ h(A\otimes I)F(z^k)
    !   ----------------               -------------------------
    !        M2n                               Fun
    !
    sc_loop: do while (error >= sc_tol) 
      ! Build M2n matrix
      call M_2n(s, t, kp, eval, Yn, step, M2n)
      ! Build Fun vector
      call Fsystem(s, t, kp, eval, Yn, step, Z_k, Fun)
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
        if (theta_k == 1.d0) then 
          error = 0.d0
        else
          eta_k = theta_k/(1.d0-theta_k)
          tol   = abs(norm_k_new - norm_k_old)/min(norm_k_new, norm_k_old)
          error = eta_k*norm_k_new/tol
        end if 
      else
        error = 1.d0
        k     = k+1
      end if 
    end do sc_loop

    ! Find Yn_hat
    Yn_hat = Yn + d1_hat*Z_k(1:dimH) + d2_hat*Z_k(dimH+1:dimH2)
    ! Find Yn
    Yn_new = Yn + d1*Z_k(1:dimH) + d2*Z_k(dimH+1:dimH2)
     
  end subroutine iterate_Zki

! Subroutine to build the matrix M_2n
  subroutine M_2n(s, t, kp, eval, Yn, step, M2n)
    use mod_f90_kind,         only: double
    use mod_parameters,       only: dimH
    use mod_system,           only: System
    use mod_imRK4_parameters, only: dimH2
    use mod_RK_matrices,      only: A, id2
    use mod_tools,            only: KronProd
    implicit none
    type(System)                           , intent(in)  :: s
    real(double)                           , intent(in)  :: t, step
    real(double)                           , intent(in)  :: kp(3)
    real(double)                           , intent(in)  :: eval
    complex(double), dimension(dimH)       , intent(in)  :: Yn
    complex(double), dimension(dimH2,dimH2), intent(out) :: M2n

    complex(double), dimension(dimH,dimH)   :: Jacobian_t
    complex(double), dimension(dimH2,dimH2) :: Kprod(dimH2,dimH2)

    ! Calculating the Jacobian at previous time
    call build_td_Jacobian(s, t-step, kp, eval, Yn, Jacobian_t)


    ! TODO: test if size(A) works
    call KronProd(size(A,1),size(A,1),dimH,dimH,A,Jacobian_t,Kprod)
    M2n = id2 - step * Kprod
  end subroutine M_2n

  !> subroutine f(Z_k) that builds the right side of the linear system.
  subroutine Fsystem(s, t, kp, eval, Yn, step, Z_k, Fun) 
    use mod_f90_kind,         only: double
    use mod_constants,        only: cI
    use mod_parameters,       only: dimH
    use mod_imRK4_parameters, only: dimH2
    use mod_system,           only: System
    use mod_RK_matrices,      only: M1, c1, c2
    implicit none
    type(System)                    ,  intent(in)  :: s
    real(double)                    ,  intent(in)  :: t, step, eval
    real(double)                    ,  intent(in)  :: kp(3)
    complex(double), dimension(dimH),  intent(in)  :: Yn
    complex(double), dimension(dimH2), intent(in)  :: Z_k
    complex(double), dimension(dimH2), intent(out) :: Fun 
    real(double)                     :: t1, t2
    complex(double)                  :: hamilt_t(dimH,dimH), hamilt_0(dimH,dimH)
    complex(double), dimension(dimH) :: F1, F2   

    t1 = t + step * c1
    call build_td_hamiltonian(s, t1, kp, eval, hamilt_t, hamilt_0)
    F1 = -cI * matmul( hamilt_t , Z_k(1:dimH) + Yn )

    t2 = t + step * c2
    call build_td_hamiltonian(s, t2, kp, eval, hamilt_t, hamilt_0)
    F2 = -cI * matmul( hamilt_t , Z_k(dimH+1:dimH2) + Yn )

    Fun = [ F1, F2 ]

    Fun = -Z_k + step * matmul(M1,Fun)
  end subroutine Fsystem


  !> build time dependent jacobian for each kp 
  subroutine build_td_Jacobian(s, t, kp, eval, Yn, Jacobian_t)
    use mod_f90_kind,   only: double
    use mod_constants,  only: cI
    use mod_parameters, only: dimH
    use mod_system,     only: System
    ! use mod_Umatrix,   only: hee
    implicit none
    type(System)                          , intent(in)  :: s
    real(double)                          , intent(in)  :: t
    real(double)                          , intent(in)  :: kp(3)
    real(double)                          , intent(in)  :: eval
    complex(double), dimension(dimH)      , intent(in)  :: Yn
    complex(double), dimension(dimH,dimH) , intent(out) :: Jacobian_t

    complex(double), dimension(dimH,dimH)  :: hamilt_t, dHdc, hamilt_0

    call build_td_hamiltonian(s, t, kp, eval, hamilt_t, hamilt_0)
    call build_term_Jacobian(s, eval, Yn, dHdc)

    Jacobian_t = -cI*(hamilt_t + dHdc)

  end subroutine build_td_Jacobian


  !> Calculate the last term of the Jacobian
  !> Given by \sum_j <i| dH/dc^n_k |j> * c_j^n(t)
  subroutine build_term_Jacobian(s, eval, Yn, dHdc)
    use mod_f90_kind,      only: double
    use mod_parameters,    only: nOrb, dimH, Um, Un, isigmamu2n, eta
    use mod_system,        only: System
    use mod_distributions, only: fd_dist
    use mod_constants,     only: cZero, pauli_mat, pi
    ! use mod_Umatrix,   only: hee
    implicit none
    type(System)                          , intent(in)  :: s
    real(double)                          , intent(in)  :: eval
    complex(double), dimension(dimH)      , intent(in)  :: Yn
    complex(double), dimension(dimH,dimH) , intent(out) :: dHdc

    integer      :: i,mu,nu,s1,s2,s3,s4,alpha

    dHdc = cZero
    do i=1,s%nAtoms
      do mu=5,nOrb
        do s1=1,2
          do s2=1,2
            dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,mu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,mu)) + Un(i)*Yn(isigmamu2n(i,s1,mu))*conjg(Yn(isigmamu2n(i,s2,mu)))
            do nu=5,nOrb
              dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) - 0.5d0*Un(i)*Yn(isigmamu2n(i,s1,mu))*conjg(Yn(isigmamu2n(i,s2,nu)))
              do s3=1,2
                do s4=1,2
                  do alpha = 1,3
                    dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) - 0.5d0*Um(i)*pauli_mat(s1,s3,alpha)*Yn(isigmamu2n(i,s3,mu))*pauli_mat(s4,s2,alpha)*conjg(Yn(isigmamu2n(i,s4,nu)))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    dHdc = fd_dist(s%Ef, 1.d0/(pi*eta), eval) * dHdc

  end subroutine build_term_Jacobian



  !> build time dependent Hamiltonian for each kp 
  !> H(t) = hk + hext_t
  subroutine build_td_hamiltonian(s,t,kp,eval,hamilt_t, hamilt_0)
    use mod_f90_kind,    only: double
    use mod_system,      only: System
    use mod_parameters,  only: dimH
    use mod_RK_matrices, only: id
    implicit none
    type(System),    intent(in)  :: s
    real(double),    intent(in)  :: t, eval
    real(double),    intent(in)  :: kp(3)
    complex(double), intent(out) :: hamilt_t(dimH,dimH), hamilt_0(dimH,dimH)

    complex(double), dimension(dimH,dimH)  :: hk,hext_t

    ! Calculating the "ground state" Hamiltonian for a given k-point (with time-dependent expectation values included)
    call hamiltk(s,kp,hk)

    ! Building time dependent hamiltonian
    call build_hext(kp,t,hext_t)

    ! calculating the original hamiltonian H(t) without the eigenvalue term.
    hamilt_0 = hk + hext_t

    ! Calculating the time-dependent Hamiltonian 
    hamilt_t = ( eval * id ) - hamilt_0

    ! Checking if Hamiltonian is hermitian
    ! if( sum(abs(conjg(transpose(hamilt_t))-hamilt_t)) > 1.d-12 ) then
    !  write(*,"('[build_td_hamiltonian] Hamiltonian is not hermitian!')")
    !  stop
    ! end if

  end subroutine build_td_hamiltonian

  !> build time dependent external perturbation Hamiltonian
  !> For a magnetic perturbation: H_ext(t)= S.B(t),  S= Pauli matricies
  !> For an electric perturbation: H_ext(t)= ((P-e*A)^2)/2*m, here only the linear term is implemented.
  subroutine build_hext(kp,t,hext_t)
    use mod_f90_kind,         only: double
    use mod_constants,        only: cI,cZero
    use mod_imRK4_parameters, only: lelectric, lmagnetic
    use mod_System,           only: ia, s => sys
    use mod_parameters,       only: nOrb,nOrb2,dimH
    implicit none
    real(double),    intent(in)  :: kp(3)
    real(double),    intent(in)  :: t
    complex(double), intent(out) :: hext_t(dimH,dimH)

    complex(double)  :: hext(nOrb2,nOrb2), temp(nOrb,nOrb)
    integer          :: i, j, k, mu, nu
    real(double)     :: b_field(3), A_t(3)
    complex(double)  :: kpExp, kpA_t

 
    hext_t = cZero

    hext = cZero
    if(lmagnetic) then
      do mu=1,nOrb
        nu=mu+nOrb
        call magnetic_field(t,b_field)
        hext(mu,mu) = hext(mu,mu) + b_field(3)
        hext(nu,nu) = hext(nu,nu) - b_field(3)
        hext(nu,mu) = hext(nu,mu) + b_field(1)-cI*b_field(2)
        hext(mu,nu) = conjg(hext(nu,mu))
      end do

      do i=1, s%nAtoms
        hext_t(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hext(1:nOrb2,1:nOrb2)
      end do
    end if

    ! The electric field is implemented via a vector potential-renormalization of the hoppings,
    ! as described in Eq. (12) of Ref.:
    ! Electromagnetic fields and dielectric response in empirical tight-binding theory
    ! M. Graf and P. Vogl Phys. Rev. B 51, 4940 (1995)
    if(lelectric) then
      call vector_potential(t, A_t)
      ! Inter-site hopping terms
      !dir$ ivdep:loop
      do k = 1, s%nNeighbors
        j = s%Neighbors(k)%BasisIndex
        ! exp(ik.(R_i-R_j))
        kpExp = exp( cI * dot_product(kp , s%Neighbors(k)%CellVector))

        !dir$ ivdep:loop
        do i = 1, s%nAtoms
          if(.not. s%Neighbors(k)%isHopping(i)) cycle
          kpA_t = exp(-cI * dot_product(A_t, s%Basis(i)%Position(:)-(s%Basis(j)%Position(:)+s%Neighbors(k)%CellVector)))

          temp(1:nOrb,1:nOrb) = s%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp * (kpA_t - 1.d0) ! The -1.d0 term is to discount the usual t(k) term that is already included in H_0
          ! Spin-up
          hext_t(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hext_t(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + temp(1:nOrb,1:nOrb)
          ! Spin-down
          hext_t(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hext_t(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + temp(1:nOrb,1:nOrb)
        end do
      end do

    end if

  end subroutine build_hext

  !> subroutine to calculate the error in the step size control
  subroutine calculate_step_error(Yn,Yn_new,Yn_hat,ERR_kn)
    use mod_f90_kind,         only: double
    use mod_RK_matrices,      only: A
    use mod_parameters,       only: dimH
    use mod_imRK4_parameters, only: abs_tol, rel_tol
    implicit none
    complex(double), dimension(dimH), intent(in) :: Yn, Yn_new, Yn_hat
    real(double), dimension(dimH)                :: Eta
    real(double),                     intent(out):: ERR_kn
    real(double)                                 :: ERR_i
    integer                                      :: i

    ! Find the difference Eta in the two solutions Yn_new and Yn_hat
    Eta = abs(Yn_new - Yn_hat)
    ! sum over the error components to get the scaled norm (ERR)
    ! ERR= sqrt[ 1\n Sum_i{ abs(Yn - Yn_hat)/( abs_tol + rel_tol * Yn(i) ) }^2 ]
    ERR_kn = 0.d0
    ERR_loop: do i=1, dimH
      ERR_i = ( Eta(i)/(abs_tol+ rel_tol* abs(Yn(i)) ) )**2 
      ERR_kn = ERR_i + ERR_kn
    end do ERR_loop
    ERR_kn = ERR_kn/dimH
    
  end subroutine calculate_step_error


 subroutine magnetic_field(t,b_field)
    use mod_f90_kind,         only: double
    use mod_imRK4_parameters, only: hw1_m, hw_m, npulse_m, lpulse_m, tau_m, delay_m, polarization_vec_m
    use mod_constants,        only: pi
    implicit none
    real(double) , intent(in)  :: t
    real(double) , intent(out) :: b_field(3)
    integer      :: np
    real(double) :: delay, arg

    b_field(:) = 0.d0
    if(lpulse_m) then
      pulses_m: do np = 1, npulse_m
        if ((t >= delay_m(np)).and.(t <= tau_m(np)+delay_m(np))) then

          delay = 0.5d0*tau_m(np) - delay_m(np)
          arg   = hw_m(np)*(t-delay)

          b_field(:) = cos(arg)*polarization_vec_m(np,1,:) + sin(arg)*polarization_vec_m(np,2,:)

          ! Cos-squared pulse
          b_field(:) =  b_field(:) * hw1_m(np) * 0.5d0 * ( cos(pi*(t-delay)/tau_m(np)) )**2  
          ! Gaussian pulse
          ! b_pulse = b_pulse * (0.5d0 * hw1_m * exp(-2*log(2*(t-delay/tau_m)**2))
          ! b_pulse = b_pulse * (0.5d0 * hw1_m * exp(-((t-delay)-4.d0*tau_m)**2/tau_m**2))

        end if
      end do pulses_m
    else
      b_field(:) = 0.5d0*hw1_m(1)*( cos(hw_m(1)*t)*polarization_vec_m(1,1,:)+sin(hw_m(1)*t)*polarization_vec_m(1,2,:) )
    end if
  end subroutine magnetic_field


  !> Subroutine builds vector potential A(t) = - integral(E(t)dt)
  ! Pulse:
  !> From paper(DOI: 0.1038/s41567-019-0602-9) the vector potential for a cos^2 pulse is given by: 
  !> A(t) = (-E_pump/w_pump) * ( cos(pi*t/tau_pump) )^2        * sin(w_pump*t), add delay_e to get:
  !> A_t  = (-hE_0/hw_e)     * (A cos(pi*(t-delay_e)/tau_e) )^2 * sin(hw_e*t)
  ! center the vector potential at delay_e
  subroutine vector_potential(t, A_t)
    use mod_f90_kind,         only: double
    use mod_imRK4_parameters, only: hE_0, hw_e, lpulse_e, npulse_e, tau_e, delay_e, polarization_vec_e
    use mod_constants,        only: pi
    implicit none 
    real(double) , intent(in)  :: t
    real(double) , intent(out) :: A_t(3)
    integer      :: np
    real(double) :: delay, arg

    A_t(:) = 0.d0
    if(lpulse_e) then
      pulses_e: do np = 1, npulse_e
        if ((t >= delay_e(np)).and.(t <= tau_e(np)+delay_e(np))) then 
          delay = 0.5d0*tau_e(np) - delay_e(np)
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

  !> Subroutine builds vector potential A(t) = - integral(E(t)dt)
  ! Pulse:
  !> From paper(DOI: 0.1038/s41567-019-0602-9) the vector potential for a cos^2 pulse is given by: 
  !> A(t) = (-E_pump/w_pump) * ( cos(pi*t/tau_pump) )^2        * sin(w_pump*t), add delay_e to get:
  !> A_t  = (-hE_0/hw_e)     * (A cos(pi*(t-delay_e)/tau_e) )^2 * sin(hw_e*t)
  ! center the vector potential at delay_e
  subroutine electric_field(t, E_t)
    use mod_f90_kind,         only: double
    use mod_imRK4_parameters, only: hE_0, hw_e, lpulse_e, npulse_e, tau_e, delay_e, polarization_vec_e
    use mod_constants,        only: pi
    implicit none 
    real(double) , intent(in)  :: t
    real(double) , intent(out) :: E_t(3)
    integer      :: np
    real(double) :: delay, arg

    E_t(:) = 0.d0
    if(lpulse_e) then
      pulses_e: do np = 1, npulse_e
        if ((t >= delay_e(np)).and.(t <= tau_e(np)+delay_e(np))) then 
          delay = 0.5d0*tau_e(np) - delay_e(np)
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

end module mod_imRK4
