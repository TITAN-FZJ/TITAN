module mod_hamiltonian
  use mod_kind, only: dp
  implicit none

  logical :: lfullhk = .false.
  !! Logical variable to say when full hk (for every k) is used
  complex(dp), dimension(:,:),   allocatable :: h0
  !! Local hamiltonian
  complex(dp), dimension(:,:,:), allocatable :: fullhk
  !! Full non-local hamiltonian hk (for every k)
#ifdef _GPU
  complex(dp), dimension(:,:),   device, allocatable :: h0_d
  !! Local hamiltonian on GPUs
  complex(dp), dimension(:,:,:), device, allocatable :: fullhk_d
  !! Full non-local hamiltonian hk (for every k) on the GPUs
#endif
  real(dp) :: energy,energy_dc,energy_dc_n,energy_constr
  !! Band energy and double counting part and constraining energy

contains

  ! Deallocate hamiltonian variables
  subroutine deallocate_hamiltonian()
    implicit none

    if(allocated(h0))     deallocate(h0)
    if(allocated(fullhk)) deallocate(fullhk)
#ifdef _GPU
    if(allocated(h0_d))     deallocate(h0_d)
    if(allocated(fullhk_d)) deallocate(fullhk_d)
#endif
  end subroutine deallocate_hamiltonian

  ! Calculate local part of the hamiltonian of the unit cell
  subroutine hamilt_local(s)
    use mod_constants,         only: cZero
    use mod_System,            only: ia,System_type
    use mod_parameters,        only: dimH,dimHsc
    use mod_Umatrix,           only: hee
    use mod_superconductivity, only: lsuperCond,bcs_pairing,delta_sc
    implicit none
    integer :: i
    type(System_type), intent(in) :: s

    if(.not.allocated(h0)) allocate( h0(dimHsc, dimHsc) )

    h0 = cZero

    ! Mouting slab hamiltonian
    ! On-site terms
    do i=1,s%nAtoms
      ! spin-up on-site tight-binding term
      h0(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = s%Types(s%Basis(i)%Material)%onSite(1:s%Types(s%Basis(i)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
      ! spin-down on-site tight-binding term
      h0(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = s%Types(s%Basis(i)%Material)%onSite(1:s%Types(s%Basis(i)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
      ! External magnetic field (orbital + spin) + Electron-electron interaction (Hubbard) + Spin-orbit coupling
      h0(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = h0 (ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                            + s%Basis(i)%lb (1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2) &
                                            + s%Basis(i)%sb (1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2) &
                                            + hee(1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2,i) &
                                            + s%Basis(i)%ls (1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2)
    end do

    ! The form of the superconducting hamiltonian depends on a series of decisions,
    ! such as how to choose the basis after the Bogoliuvob-de Gennes transformation
    ! and how do we define this transformation. In our particular case, we choose to
    ! have a basis (u_up, u_down, v_up, v_down), and the BdG transformation we perform
    ! is c_{i\sigma} = sum_n u_{in\sigma}\gamma_n + v_{in\sigma}\gamma^{n*}\gamma_n^\dagger
    ! The final form of the hamiltonian is something like
    ! | H(k) - E_f        Delta     |
    ! |   Delta^*   -(H(-k) - E_f)* |
    ! roughly. Look at this paper 10.1103/RevModPhys.87.1037 , and to Uriel's thesis to
    ! get a better idea of how to construct this operator

    if(lsuperCond) then
      h0(dimH+1:dimHsc,dimH+1:dimHsc) = -conjg(h0(1:dimH,1:dimH))
      do i = 1,dimH
        h0(     i,     i) = h0(     i,     i) - s%Ef
        h0(dimH+i,dimH+i) = h0(dimH+i,dimH+i) + s%Ef
      end do
      ! Populating the non-diagonal blocks of the hamiltonian. There are several ways to do it.
      call bcs_pairing(s, delta_sc, h0)
    end if

  end subroutine hamilt_local


#ifdef _GPU
  ! Calculate local part of the hamiltonian of the unit cell
  subroutine hamilt_local_gpu(s)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero
    use mod_System,            only: ia_d,System_type
    use mod_parameters,        only: dimH,dimHsc,isigmamu2n_d
    use mod_Umatrix,           only: hee
    use mod_superconductivity, only: lsuperCond,delta_sc_d
    implicit none
    type(System_type), intent(in) :: s
    complex(dp), dimension(s%nOrb ,s%nOrb ,s%nAtoms), device :: onSite_d
    complex(dp), dimension(s%nOrb2,s%nOrb2,s%nAtoms), device :: lb_d,sb_d,hee_d,ls_d
    integer :: i,j,mu,nOrb_i,nOrb2_i
    integer, device :: nOrb_d(s%nAtoms)

    if(.not.allocated(h0_d)) allocate( h0_d(dimHsc, dimHsc) )

    h0_d = cZero

    do i=1,s%nAtoms
      nOrb_d(i) = nOrb_i
      nOrb_i  = s%Types(s%Basis(i)%Material)%nOrb
      nOrb2_i = s%Types(s%Basis(i)%Material)%nOrb2
      onSite_d(:,:,i) = s%Types(s%Basis(i)%Material)%onSite(1:nOrb_i,1:nOrb_i)
      ! s%Basis(i)%sb_d(1:nOrb2_i,1:nOrb2_i) = s%Basis(i)%sb(1:nOrb2_i,1:nOrb2_i)
      ! s%Basis(i)%lb_d(1:nOrb2_i,1:nOrb2_i) = s%Basis(i)%lb(1:nOrb2_i,1:nOrb2_i)
      ! s%Basis(i)%ls_d(1:nOrb2_i,1:nOrb2_i) = s%Basis(i)%ls(1:nOrb2_i,1:nOrb2_i)
      sb_d(1:nOrb2_i,1:nOrb2_i,i) = s%Basis(i)%sb(1:nOrb2_i,1:nOrb2_i)
      lb_d(1:nOrb2_i,1:nOrb2_i,i) = s%Basis(i)%lb(1:nOrb2_i,1:nOrb2_i)
      ls_d(1:nOrb2_i,1:nOrb2_i,i) = s%Basis(i)%ls(1:nOrb2_i,1:nOrb2_i)
    end do
    hee_d = hee


    ! Mouting slab hamiltonian
    ! On-site terms
    !$cuf kernel do <<< *, * >>>
    do i=1,s%nAtoms
      nOrb2_i = s%Types(s%Basis(i)%Material)%nOrb2
      ! spin-up on-site tight-binding term
      h0_d(ia_d(1,i):ia_d(2,i),ia_d(1,i):ia_d(2,i)) = onSite_d(:,:,i)
      ! spin-down on-site tight-binding term
      h0_d(ia_d(3,i):ia_d(4,i),ia_d(3,i):ia_d(4,i)) = onSite_d(:,:,i)
      ! External magnetic field (orbital + spin) + Electron-electron interaction (Hubbard) + Spin-orbit coupling
      h0_d(ia_d(1,i):ia_d(4,i),ia_d(1,i):ia_d(4,i)) = h0_d (ia_d(1,i):ia_d(4,i), ia_d(1,i):ia_d(4,i)) &
                                                    + lb_d (1:nOrb2_i,1:nOrb2_i,i) &
                                                    + sb_d (1:nOrb2_i,1:nOrb2_i,i) &
                                                    + hee_d(1:nOrb2_i,1:nOrb2_i,i) &
                                                    + ls_d (1:nOrb2_i,1:nOrb2_i,i)
    end do

    ! The form of the superconducting hamiltonian depends on a series of decisions,
    ! such as how to choose the basis after the Bogoliuvob-de Gennes transformation
    ! and how do we define this transformation. In our particular case, we choose to
    ! have a basis (u_up, u_down, v_up, v_down), and the BdG transformation we perform
    ! is c_{i\sigma} = sum_n u_{in\sigma}\gamma_n + v_{in\sigma}\gamma^{n*}\gamma_n^\dagger
    ! The final form of the hamiltonian is something like
    ! | H(k) - E_f        Delta     |
    ! |   Delta^*   -(H(-k) - E_f)* |
    ! roughly. Look at this paper 10.1103/RevModPhys.87.1037 , and to Uriel's thesis to
    ! get a better idea of how to construct this operator

    if(lsuperCond) then
      !$cuf kernel do(2) <<< *, * >>>
      do j=1,dimH
        do i=1,dimH
          h0_d(dimH+i,dimH+j) = -conjg(h0_d(i,j))
        end do
      end do

      !$cuf kernel do <<< *, * >>>
      do i = 1,dimH
        h0_d(     i,     i) = h0_d(     i,     i) - s%Ef
        h0_d(dimH+i,dimH+i) = h0_d(dimH+i,dimH+i) + s%Ef
      end do

      ! Populating the non-diagonal blocks of the hamiltonian. There are several ways to do it.
      !$cuf kernel do <<< (*,*), (*,*) >>>
      do i = 1,s%nAtoms
        do mu = 1,nOrb_d(i)
          h0_d(isigmamu2n_d(i,1,mu)     ,isigmamu2n_d(i,2,mu)+dimH) = - cmplx(delta_sc_d(mu,i),0._dp,dp)
          h0_d(isigmamu2n_d(i,2,mu)     ,isigmamu2n_d(i,1,mu)+dimH) =   cmplx(delta_sc_d(mu,i),0._dp,dp)
          h0_d(isigmamu2n_d(i,2,mu)+dimH,isigmamu2n_d(i,1,mu)     ) = - cmplx(delta_sc_d(mu,i),0._dp,dp)
          h0_d(isigmamu2n_d(i,1,mu)+dimH,isigmamu2n_d(i,2,mu)     ) =   cmplx(delta_sc_d(mu,i),0._dp,dp)
        end do
      end do
    end if

  end subroutine hamilt_local_gpu
#endif


  ! Calculate the k-dependent tight-binding hamiltonian of the unit cell
  function calchk(s,kp)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,cI
    use mod_System,            only: ia,ia_sc,System_type
    use mod_parameters,        only: dimHsc
    use mod_superconductivity, only: lsuperCond
    implicit none
    real(dp),          intent(in) :: kp(3)
    type(System_type), intent(in) :: s
    complex(dp), dimension(dimHsc,dimHsc) :: calchk

    integer     :: i,j,k,ia_temp_i,ia_temp_j
    complex(dp) :: tmp(s%nOrb,s%nOrb)
    complex(dp) :: kpExp,mkpExp

    calchk = cZero

    ! Inter-site hopping terms
    do k = 1, s%nNeighbors
      j = s%Neighbors(k)%BasisIndex
      ! exp(ik.(R_i-R_j))
      kpExp = exp(cI * dot_product(kp, s%Neighbors(k)%CellVector))

      if(lsuperCond) mkpExp = conjg(kpExp)

      do i = 1,s%nAtoms
        if(s%Neighbors(k)%isHopping(i)) then
          ! electrons
          tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb) = s%Neighbors(k)%t0i(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb,i) * kpExp
          ! Spin-up
          calchk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = calchk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
          ! Spin-down
          calchk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = calchk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)

          if(lsuperCond) then
            ! holes
            tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb) = s%Neighbors(k)%t0i(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb, i) * mkpExp
            ! Spin-up
            ia_temp_j = ia_sc(3,j) + s%Types(s%Basis(j)%Material)%nOrb - 1
            ia_temp_i = ia_sc(3,i) + s%Types(s%Basis(i)%Material)%nOrb - 1 
            calchk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) = calchk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) - conjg(tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb))
            ! Spin-down
            ia_temp_j = ia_sc(4,j) - s%Types(s%Basis(j)%Material)%nOrb + 1
            ia_temp_i = ia_sc(4,i) - s%Types(s%Basis(i)%Material)%nOrb + 1
            calchk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) = calchk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) - conjg(tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb))
          end if
        end if
      end do
    end do

    ! Test if hamiltonian is Hermitian (to be commented out, uncomment to use it)
    ! do i = ia(1,1), ia(4,s%nAtoms)
    !   do j = i, ia(4,s%nAtoms)
    !     if(abs(calchk(j,i)-conjg(calchk(i,j))) > 1.e-14_dp) then
    !       write(*,"('Hamiltonian not hermitian',i0,2x,i0,2x,es11.4)") i,j,abs(calchk(j,i)-conjg(calchk(i,j)))
    !     end if
    !   end do
    ! end do

  end function calchk


  ! Calculate the k-dependent tight-binding hamiltonian of the unit cell for all k-points
  function fullhamiltk(s) result(success)
    use mod_kind,              only: int64
    use mod_BrillouinZone,     only: realBZ
    use mod_System,            only: System_type
    use mod_parameters,        only: dimHsc,output
    use mod_mpi_pars,          only: rField
    use mod_tools,             only: get_memory
    use mod_progress,          only: write_time
    implicit none
    type(System_type), intent(in) :: s
    integer(int64) :: iz

    logical :: success
    integer :: mem_avail,mem_req

    success = .false.
    ! Getting available memory
    if(.not.get_memory("m",mem_avail)) return

    ! Comparing memory needed with 90% of available memory
    !       nkpt      *  hamilt dim       *cmplx (kb) (mb)
    mem_req = int(dble((dimHsc)**2)*16/1024/1024*dble(realBZ%workload))
    if(mem_req > floor(0.9*mem_avail)) then
      if(rField == 0) &
        write(unit=output%unit, fmt="('[fullhamiltk] Full tight-binding hamiltonian does not fit in 90% of available memory. It must be re-calculated for each k-point.')")
        write(unit=output%unit, fmt="('[fullhamiltk] Required memory: ',i0,'Mb   //  Available memory: ',i0,'Mb')") mem_req, mem_avail
      return
    end if

    if(rField == 0) then
      write(unit=output%unit, fmt="('[fullhamiltk] Full tight-binding hamiltonian fit in memory. Allocating and calculating it.')")
      write(unit=output%unit, fmt="('[fullhamiltk] Required memory: ',i0,'Mb   //  Available memory: ',i0,'Mb')") mem_req, mem_avail
    end if

    allocate(fullhk(dimHsc,dimHsc,realBZ%workload))

    ! Inter-site hopping terms
    !$omp parallel do default(none) schedule(dynamic) &
    !$omp& shared(fullhk,realBZ,dimHsc,s) &
    !$omp& private(iz)
    kloop: do iz = 1,realBZ%workload
      ! Calculating the hamiltonian for a given k-point
      fullhk(:,:,iz) = calchk(s,realBZ%kp(1:3,iz))
    end do kloop
    !$omp end parallel do

#ifdef _GPU
    allocate(fullhk_d(dimHsc,dimHsc,realBZ%workload))
    fullhk_d = fullhk
#endif

    if(rField == 0) call write_time('[fullhamiltk] Finished calculating full Hamiltonian on: ',output%unit)
    success = .true.

  end function fullhamiltk


  ! Calculate hamiltonian of the unit cell
  ! and the spin-orbit coupling contribution separately
  subroutine hamiltklinearsoc(s,kp,hk,vsoc)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,cI
    use mod_system,            only: ia,ia_sc,System_type
    use mod_parameters,        only: dimH,dimHsc
    use mod_Umatrix,           only: hee
    use mod_superconductivity, only: lsuperCond,bcs_pairing,delta_sc
    implicit none
    real(dp),          intent(in) :: kp(3)
    type(System_type), intent(in) :: s

    integer     :: i,j,k,ia_temp_i,ia_temp_j
    complex(dp) :: tmp(s%nOrb,s%nOrb)
    complex(dp) :: kpExp,mkpExp
    complex(dp), dimension(dimHsc,dimHsc), intent(out) :: hk,vsoc

    hk = cZero
    vsoc = cZero

    ! Mouting slab hamiltonian
    ! On-site terms
    do i=1,s%nAtoms
      ! spin-up on-site tight-binding term
      hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = s%Types(s%Basis(i)%Material)%onSite(1:s%Types(s%Basis(i)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
      ! spin-down on-site tight-binding term
      hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = s%Types(s%Basis(i)%Material)%onSite(1:s%Types(s%Basis(i)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
      ! External magnetic field (orbital + spin) + Electron-electron interaction (Hubbard) + Spin-orbit coupling
      hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) =  hk (ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                            + s%Basis(i)%lb (1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2) &
                                            + s%Basis(i)%sb (1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2) &
                                            + hee(1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2,i)
      vsoc(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = s%Basis(i)%ls(1:s%Types(s%Basis(i)%Material)%nOrb2,1:s%Types(s%Basis(i)%Material)%nOrb2)
    end do

    ! The form of the superconducting hamiltonian depends on a series of decisions,
    ! such as how to choose the basis after the Bogoliuvob-de Gennes transformation
    ! and how do we define this transformation. In our particular case, we choose to
    ! have a basis (u_up, u_down, v_up, v_down), and the BdG transformation we perform
    ! is c_{i\sigma} = sum_n u_{in\sigma}\gamma_n + v_{in\sigma}\gamma^{n*}\gamma_n^\dagger
    ! The final form of the hamiltonian is something like
    ! | H(k) - E_f        Delta     |
    ! |   Delta^*   -(H(-k) - E_f)* |
    ! roughly. Look at this paper 10.1103/RevModPhys.87.1037 , and to Uriel's thesis to
    ! get a better idea of how to construct this operator

    if(lsuperCond) then
      hk(dimH+1:dimHsc,dimH+1:dimHsc) = -conjg(hk(1:dimH,1:dimH))
      vsoc(dimH+1:dimHsc,dimH+1:dimHsc) = -conjg(vsoc(1:dimH,1:dimH))
      do i = 1,dimH
        hk(     i,     i) = hk(     i,     i) - s%Ef
        hk(dimH+i,dimH+i) = hk(dimH+i,dimH+i) + s%Ef
      end do
      ! Populating the non-diagonal blocks of the hamiltonian. There are several ways to do it.
      call bcs_pairing(s, delta_sc, hk)
    end if

    ! Inter-site hopping terms
    do k = 1, s%nNeighbors
      j = s%Neighbors(k)%BasisIndex
      ! exp(ik.(R_i-R_j))
      kpExp = exp(cI * dot_product(kp, s%Neighbors(k)%CellVector))

      if(lsuperCond) mkpExp = conjg(kpExp)

      do i = 1,s%nAtoms
        if(s%Neighbors(k)%isHopping(i)) then
          ! electrons
          tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb) = s%Neighbors(k)%t0i(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb,i) * kpExp
          ! Spin-up
          hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
          ! Spin-down
          hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)

          if(lsuperCond) then
            ! holes
            tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb) = s%Neighbors(k)%t0i(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb, i) * mkpExp
            ! Spin-up
            ia_temp_j = ia_sc(3,j) + s%Types(s%Basis(j)%Material)%nOrb - 1
            ia_temp_i = ia_sc(3,i) + s%Types(s%Basis(i)%Material)%nOrb - 1 
            hk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) = hk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) - conjg(tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb))
            ! Spin-down
            ia_temp_j = ia_sc(4,j) - s%Types(s%Basis(j)%Material)%nOrb + 1
            ia_temp_i = ia_sc(4,i) - s%Types(s%Basis(i)%Material)%nOrb + 1
            hk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) = hk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) - conjg(tmp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb))
          end if
        end if
      end do
    end do

  end subroutine hamiltklinearsoc

  !> build time dependent external perturbation Hamiltonian
  !> For a magnetic perturbation: H_ext(t)= S.B(t),  S= Pauli matricies
  !> For an electric perturbation: H_ext(t)= ((P-e*A)^2)/2*m, here only the linear term is implemented.
  subroutine build_hext(kp,b_field,A_t,hext_t)
    use mod_kind,              only: dp
    use mod_constants,         only: cI,cZero
    use mod_imRK4_parameters,  only: lelectric,lmagnetic
    use mod_System,            only: ia,ia_sc,s => sys
    use mod_parameters,        only: dimHsc
    use mod_superconductivity, only: lsuperCond
    implicit none
    real(dp),    intent(in)  :: kp(3)
    real(dp),    intent(in)  :: b_field(3), A_t(3)
    complex(dp), intent(out) :: hext_t(dimHsc,dimHsc)

    complex(dp)  :: hext(s%nOrb2,s%nOrb2), temp(s%nOrb,s%nOrb)
    integer      :: i, j, k, mu, nu, ia_temp_i, ia_temp_j
    complex(dp)  :: kpExp,mkpExp,expA,kpA_t,mkpA_t


    hext_t = cZero

    if(lmagnetic) then
      do i=1, s%nAtoms
        hext = cZero
        do mu=1,s%Types(s%Basis(i)%Material)%nOrb
          nu=mu+s%Types(s%Basis(i)%Material)%nOrb
          hext(mu,mu) = hext(mu,mu) + b_field(3)
          hext(nu,nu) = hext(nu,nu) - b_field(3)
          hext(nu,mu) = hext(nu,mu) + b_field(1)+cI*b_field(2)
          hext(mu,nu) = conjg(hext(nu,mu))
        end do

        hext_t(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hext(1:s%Types(s%Basis(i)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
      end do
    end if

    ! The electric field is implemented via a vector potential-renormalization of the hoppings,
    ! as described in Eq. (12) of Ref.:
    ! Electromagnetic fields and dielectric response in empirical tight-binding theory
    ! M. Graf and P. Vogl Phys. Rev. B 51, 4940 (1995)
    if(lelectric) then
      ! Inter-site hopping terms
      do k = 1, s%nNeighbors
        j = s%Neighbors(k)%BasisIndex
        ! exp(ik.(R_i-R_j))
        kpExp = exp( cI * dot_product(kp,s%Neighbors(k)%CellVector))

        if(lsuperCond) mkpExp = conjg(kpExp)

        do i = 1,s%nAtoms
          if(s%Neighbors(k)%isHopping(i)) then
            expA = ( exp(-cI * dot_product(A_t, s%Basis(i)%Position(:)-(s%Basis(j)%Position(:)+s%Neighbors(k)%CellVector))) - 1._dp) ! The -1._dp term is to discount the usual t(k) term that is already included in H_0
            kpA_t =  kpExp * expA 

            !DIR$ VECTOR ALIGNED
            temp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb) = s%Neighbors(k)%t0i(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb,i) * kpA_t 
            ! Spin-up
            hext_t(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hext_t(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + temp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)
            ! Spin-down
            hext_t(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hext_t(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + temp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb)

            if(lsuperCond) then
              mkpA_t =  mkpExp * expA 
              ! holes
              temp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb) = s%Neighbors(k)%t0i(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb, i) * mkpA_t
              ! Spin-up
              ia_temp_j = ia_sc(3,j) + s%Types(s%Basis(j)%Material)%nOrb - 1
              ia_temp_i = ia_sc(3,i) + s%Types(s%Basis(i)%Material)%nOrb - 1 
              hext_t(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) = hext_t(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) - conjg(temp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb))
              ! Spin-down
              ia_temp_j = ia_sc(4,j) - s%Types(s%Basis(j)%Material)%nOrb + 1
              ia_temp_i = ia_sc(4,i) - s%Types(s%Basis(i)%Material)%nOrb + 1
              hext_t(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) = hext_t(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) - conjg(temp(1:s%Types(s%Basis(j)%Material)%nOrb,1:s%Types(s%Basis(i)%Material)%nOrb))
            end if
          end if

        end do
      end do
    end if

  end subroutine build_hext

end module mod_hamiltonian
