!> module for imRK4 routines
module mod_imRK4
  implicit none

contains
  ! subroutine to find the vectors Z_ki
  subroutine iterate_Zki(s, t, kp, eval, Yn)
    use mod_f90_kind,         only: double
    use mod_parameters,       only: dimH
    use mod_system,           only: System
    use mod_imRK4_parameters, only: step, dimH2, sc_tol
    use mod_tools,            only: vec_norm, LS_solver
    use mod_constants,        only: cZero
    use mod_RK_matrices
    implicit none
    ! define the variables: 
    type(System)                     , intent(in)    :: s
    real(double)                     , intent(in)    :: t
    real(double)                     , intent(in)    :: kp(3)
    real(double)                     , intent(in)    :: eval
    complex(double), dimension(dimH) , intent(inout) :: Yn

    integer                                 :: k 
    real(double)                            :: error, norm_k_old, norm_k_new, theta_k, eta_k, tol
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
      call Fsystem(s, t, kp, Yn, step, Z_k, Fun)
      ! Solve the Ax=b linear equation
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

    Yn = Yn + d1*Z_k(1:dimH) + d2*Z_k(dimH+1:dimH2)
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
  subroutine Fsystem(s, t, kp, Yn, step, Z_k, Fun) 
    use mod_f90_kind,         only: double
    use mod_constants,        only: cI
    use mod_parameters,       only: dimH
    use mod_imRK4_parameters, only: dimH2
    use mod_system,           only: System
    use mod_RK_matrices,      only: M1, c1, c2
    implicit none
    type(System)                    ,  intent(in)  :: s
    real(double)                    ,  intent(in)  :: t, step
    real(double)                    ,  intent(in)  :: kp(3)
    complex(double), dimension(dimH),  intent(in)  :: Yn
    complex(double), dimension(dimH2), intent(in)  :: Z_k
    complex(double), dimension(dimH2), intent(out) :: Fun 
    real(double)                     :: t1, t2
    complex(double)                  :: hamilt_t(dimH,dimH)
    complex(double), dimension(dimH) :: F1, F2   

    t1 = t + step * c1
    call build_td_hamiltonian(s,t1,kp,hamilt_t)
    F1 = -cI * matmul( hamilt_t , Z_k(1:dimH) + Yn )

    t2 = t + step * c2
    call build_td_hamiltonian(s,t2,kp,hamilt_t)
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

    complex(double), dimension(dimH,dimH)  :: hamilt_t,dHdc

    call build_td_hamiltonian(s, t, kp, hamilt_t)
    call build_term_Jacobian(s, eval, Yn, dHdc)

    Jacobian_t = -cI*(hamilt_t + dHdc)

  end subroutine build_td_Jacobian


  !> Calculate the last term of the Jacobian
  !> Given by \sum_j <i| dH/dc^n_k |j> * c_j^n(t)
  subroutine build_term_Jacobian(s, eval, Yn, dHdc)
    use mod_f90_kind,      only: double
    use mod_parameters,    only: dimH, U, isigmamu2n, eta
    use mod_system,        only: System
    use mod_distributions, only: fd_dist
    use TightBinding,      only: nOrb
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
            dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,mu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,mu)) + Yn(isigmamu2n(i,s1,mu))*conjg(Yn(isigmamu2n(i,s2,mu)))
            do nu=5,nOrb
              dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) - 0.5d0*Yn(isigmamu2n(i,s1,mu))*conjg(Yn(isigmamu2n(i,s2,nu)))
              do s3=1,2
                do s4=1,2
                  do alpha = 1,3
                    dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) = dHdc(isigmamu2n(i,s1,mu),isigmamu2n(i,s2,nu)) - 0.5d0*pauli_mat(s1,s3,alpha)*Yn(isigmamu2n(i,s3,mu))*pauli_mat(s4,s2,alpha)*conjg(Yn(isigmamu2n(i,s4,nu)))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      dHdc = U(i) * dHdc
    end do
    dHdc = fd_dist(s%Ef, 1.d0/(pi*eta), eval) * dHdc

  end subroutine build_term_Jacobian



  !> build time dependent Hamiltonian for each kp 
  !> H(t)= S.B(t),  S= Pauli matricies of dimension 18? or dimH
  subroutine build_td_hamiltonian(s,t,kp,hamilt_t)
    use mod_f90_kind,   only: double
    use mod_system,     only: System
    use mod_parameters, only: dimH
    implicit none
    type(System),   intent(in) :: s
    real(double),   intent(in) :: t
    real(double),   intent(in) :: kp(3)
    complex(double), intent(out)  :: hamilt_t(dimH,dimH)
    complex(double), dimension(dimH,dimH)  :: hk

    ! Calculating the "ground state" Hamiltonian for a given k-point (with time-dependent expectation values included)
    call hamiltk(s,kp,hk)

    ! Calculating the time-dependent Hamiltonian
    hamilt_t = hk + hext_t(s%nAtoms,t)

    ! Checking if Hamiltonian is hermitian
    if( sum(abs(conjg(transpose(hamilt_t))-hamilt_t)) > 1.d-12 ) then
      write(*,"('[build_td_hamiltonian] Hamiltonian is not hermitian!')")
      stop
    end if
  end subroutine build_td_hamiltonian

  !> build time dependent external perturbation Hamiltonian
  !> H_ext(t)= S.B(t),  S= Pauli matricies
  function hext_t(nAtoms,t)
    use mod_f90_kind,         only: double
    use mod_constants,        only: cI,cZero
    use mod_imRK4_parameters, only: hw1, hw
    use TightBinding,         only: nOrb,nOrb2
    use mod_System,           only: ia
    use mod_parameters,       only: dimH
    implicit none
    real(double)     :: t
    integer          :: nAtoms
    complex(double)  :: hext_t(dimH,dimH)

    complex(double)  :: hext(nOrb2,nOrb2)
    integer          :: i, mu, nu

    hext = cZero
    do mu=1,nOrb
      nu=mu+nOrb
      hext(nu,mu) = (cos(hw*t) +cI*sin(hw*t))*hw1
      hext(mu,nu) = conjg(hext(nu,mu))
    end do

    hext_t = cZero
    do i=1, nAtoms
      hext_t(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hext(1:nOrb2,1:nOrb2)
    end do

  end function hext_t

end module mod_imRK4

