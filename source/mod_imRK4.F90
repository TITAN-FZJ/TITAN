module mod_imRK4
  !! module for the implicit 4th-order Runge-Kutta method (imRK4) routines
  implicit none

contains
  subroutine iterate_Zki(s,b_fieldm,A_tm,b_field1,A_t1,b_field2,A_t2,hamilt_nof,kp,eval,step,Yn_new,Yn,Yn_hat)
    !! Calculates the vectors Z_ki
    use mod_kind,             only: dp
    use mod_parameters,       only: dimHsc
    use mod_system,           only: System_type
    use mod_time_propagator,  only: dimH2, sc_tol
    use mod_tools,            only: vec_norm, LS_solver
    use mod_constants,        only: cZero
    use mod_RK_matrices,      only: d1,d2,d1_hat,d2_hat
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
    use mod_time_propagator,  only: dimH2
    use mod_RK_matrices,      only: A, id2
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
    use mod_time_propagator,  only: dimH2
    use mod_RK_matrices,      only: M1
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
    call build_hext(kp,b_field1,A_t1,hext_t)
    hamilt_t = hamilt_nof - hext_t

    !F1 = -cI * matmul(hamilt_t,Z_k(1:dimHsc) + Yn)
    tempv = Z_k(1:dimHsc)+Yn
    call zgemv('n',dimHsc,dimHsc,-cI,hamilt_t,dimHsc,tempv,1,cZero,F1,1)

    ! Building time dependent hamiltonian
    call build_hext(kp,b_field2,A_t2,hext_t)
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
    call build_hext(kp,b_field,A_t,hext_t)
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
    use mod_time_propagator,  only: abs_tol,rel_tol
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
    use mod_time_propagator,  only: hw1_m, hw_m, npulse_m, lpulse_m, tau_m, delay_m, polarization_vec_m
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
    use mod_time_propagator,  only: hE_0, hw_e, lpulse_e, npulse_e, tau_e, delay_e, polarization_vec_e
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
    use mod_time_propagator,  only: hE_0, hw_e, lpulse_e, npulse_e, tau_e, delay_e, polarization_vec_e
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

end module mod_imRK4
