module mod_hamiltonian
  use mod_kind, only: dp
  implicit none
  
  complex(dp), dimension(:,:),   allocatable :: h0
  complex(dp), dimension(:,:,:), allocatable :: fullhk

contains

  ! Deallocate hamiltonian variables
  subroutine deallocate_hamiltonian()
    implicit none

    if(allocated(h0))     deallocate(h0)
    if(allocated(fullhk)) deallocate(fullhk)
  end subroutine deallocate_hamiltonian

  ! Calculate local part of the hamiltonian of the unit cell
  subroutine hamilt_local(sys)
    use mod_constants,  only: cZero
    use mod_System,     only: ia, System_type
    use mod_parameters, only: nOrb,nOrb2,dimH
    use mod_magnet,     only: lb, sb
    use mod_SOC,        only: ls
    use mod_Umatrix,    only: hee
    use mod_superconductivity, only: lsuperCond,supercond,bcs_pairing,singlet_coupling
    implicit none
    integer :: i
    type(System_type), intent(in) :: sys

    if(.not.allocated(h0)) allocate( h0(dimH*supercond, dimH*supercond) )

    h0 = cZero

    ! Mouting slab hamiltonian
    ! On-site terms
    do concurrent (i=1:sys%nAtoms)
      ! spin-up on-site tight-binding term
      h0(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
      ! spin-down on-site tight-binding term
      h0(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
      ! External magnetic field (orbital + spin) + Electron-electron interaction (Hubbard) + Spin-orbit coupling
      h0(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = h0(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                           + lb(1:nOrb2,1:nOrb2,i) + sb(1:nOrb2,1:nOrb2,i) + hee(1:nOrb2,1:nOrb2,i) &
                                           + ls(1:nOrb2,1:nOrb2,i)
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
      h0(dimH+1:2*dimH,dimH+1:2*dimH) = -conjg(h0(1:dimH,1:dimH))
      do concurrent (i = 1:dimH)
        h0(     i,     i) = h0(     i,     i) - sys%Ef
        h0(dimH+i,dimH+i) = h0(dimH+i,dimH+i) + sys%Ef
      end do
      ! Populating the non-diagonal blocks of the hamiltonian. There are several ways to do it.
      call bcs_pairing(sys, singlet_coupling, h0)
    end if


  end subroutine hamilt_local


  ! Calculate the k-dependent tight-binding hamiltonian of the unit cell
  subroutine hamiltk(sys,kp,hk)
    use mod_kind, only: dp
    use mod_constants,  only: cI
    use mod_System,     only: ia,ia_sc,System_type
    use mod_parameters, only: nOrb,dimH
    use mod_superconductivity, only: lsuperCond,supercond
    implicit none
    real(dp), intent(in) :: kp(3)
    type(System_type), intent(in) :: sys
    complex(dp), dimension(dimH*supercond,dimH*supercond), intent(out) :: hk
    integer :: i, j, k , ia_temp_i, ia_temp_j
    complex(dp) :: tmp(nOrb,nOrb)
    complex(dp) :: kpExp,mkpExp

    hk = h0

    ! Inter-site hopping terms
    do k = 1, sys%nNeighbors
      j = sys%Neighbors(k)%BasisIndex
      ! exp(ik.(R_i-R_j))
      kpExp = exp(cI * dot_product(kp, sys%Neighbors(k)%CellVector))

      if(lsuperCond) mkpExp = conjg(kpExp)

      do concurrent (i = 1:sys%nAtoms, sys%Neighbors(k)%isHopping(i))
        ! electrons
        tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
        ! Spin-up
        hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:nOrb,1:nOrb)
        ! Spin-down
        hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:nOrb,1:nOrb)

        if(lsuperCond) then
          ! holes
          tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * mkpExp
          ! Spin-up
          ia_temp_j = ia_sc(3,j) + nOrb - 1
          ia_temp_i = ia_sc(3,i) + nOrb - 1 
          hk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) = hk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i) - conjg(tmp(1:nOrb,1:nOrb))
          ! Spin-down
          ia_temp_j = ia_sc(4,j) - nOrb + 1
          ia_temp_i = ia_sc(4,i) - nOrb + 1
          hk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) = hk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i)) - conjg(tmp(1:nOrb,1:nOrb))
        end if
      end do
    end do

    ! This part is neccessary since hamiltk returns a hamiltonian, hermitian
    ! up to 1e-16 precision. This interfered with some expectation value calculations
    ! causing an artificial anisotropy
    ! With this I'm just making sure that diagonal blocks of the superconducting
    ! Hamiltonian are hermitian, and therefore the total hamiltonian is hermitian.

    hk = (hk + conjg(transpose(hk)))/2._dp

    ! ! Test if hamiltonian is Hermitian (to be commented out, uncomment to use it)
    ! do i = ia(1,1), ia(4,sys%nAtoms)
    !   do j = i, ia(4,sys%nAtoms)
    !     if(abs(hk(j,i)-conjg(hk(i,j))) > 1.e-15_dp) then
    !       write(*,"('Hamiltonian not hermitian',i0,2x,i0,2x,es11.4)") i,j,abs(hk(j,i)-conjg(hk(i,j)))
    !     end if
    !   end do
    ! end do

  end subroutine hamiltk


  ! Calculate the k-dependent tight-binding hamiltonian of the unit cell for all k-points
  function fullhamiltk(sys) result(success)
    use mod_kind, only: dp,int32,int64
    use mod_BrillouinZone, only: realBZ
    use mod_constants, only: cI,cZero
    use mod_System,     only: ia,ia_sc,System_type
    use mod_parameters, only: nOrb,dimH,output
    use mod_mpi_pars, only: rField
    use mod_tools, only: get_memory
    use mod_superconductivity, only: lsuperCond,supercond
    use mod_progress, only: write_time
    implicit none
    type(System_type), intent(in) :: sys
    integer(int32) :: i, j, k, ia_temp_i, ia_temp_j
    integer(int64) :: iz
    complex(dp) :: tmp(nOrb,nOrb)
    complex(dp) :: kpExp,mkpExp

    logical :: success
    integer :: mem_avail,mem_req

    success = .false.
    ! Getting available memory
    if(.not.get_memory("m",mem_avail)) return

    ! Comparing memory needed with 90% of available memory
    !       nkpt      *  hamilt dim       *cmplx (kb) (mb)
    mem_req = int((dimH*supercond)**2* 16/1024/1024*realBZ%workload)
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

    allocate(fullhk(dimH*supercond,dimH*supercond,realBZ%workload))
    fullhk = cZero

    ! Inter-site hopping terms
    !$omp parallel do default(none) shared(fullhk,realBZ,nOrb,ia,sys,lsuperCond,ia_sc) private(i,j,k,iz,ia_temp_i,ia_temp_j,kpExp,mkpExp,tmp) schedule(dynamic)
    do iz = 1,realBZ%workload

      do k = 1, sys%nNeighbors
        j = sys%Neighbors(k)%BasisIndex
        ! exp(ik.(R_i-R_j))
        kpExp = exp(cI * dot_product(realBZ%kp(1:3,iz), sys%Neighbors(k)%CellVector))

        if(lsuperCond) mkpExp = conjg(kpExp)

        do i = 1,sys%nAtoms
          if(.not.sys%Neighbors(k)%isHopping(i)) cycle

          tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
          ! Spin-up
          fullhk(ia(1,j):ia(2,j), ia(1,i):ia(2,i),iz) = fullhk(ia(1,j):ia(2,j), ia(1,i):ia(2,i),iz) + tmp(1:nOrb,1:nOrb)
          ! Spin-down
          fullhk(ia(3,j):ia(4,j), ia(3,i):ia(4,i),iz) = fullhk(ia(3,j):ia(4,j), ia(3,i):ia(4,i),iz) + tmp(1:nOrb,1:nOrb)

          if(lsuperCond) then
            ! holes
            tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * mkpExp
            ! Spin-up
            ia_temp_j = ia_sc(3,j) + nOrb - 1
            ia_temp_i = ia_sc(3,i) + nOrb - 1 
            fullhk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i,iz) = fullhk(ia_sc(3,j):ia_temp_j, ia_sc(3,i):ia_temp_i,iz) - conjg(tmp(1:nOrb,1:nOrb))
            ! Spin-down
            ia_temp_j = ia_sc(4,j) - nOrb + 1
            ia_temp_i = ia_sc(4,i) - nOrb + 1
            fullhk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i),iz) = fullhk(ia_temp_j:ia_sc(4,j), ia_temp_i:ia_sc(4,i),iz) - conjg(tmp(1:nOrb,1:nOrb))
          end if

        end do

      end do
    end do
    !$omp end parallel do

    if(rField == 0) call write_time(output%unit,'[fullhamiltk] Finished calculating full Hamiltonian on: ')
    success = .true.

  end function fullhamiltk


  ! Calculate hamiltonian of the unit cell
  ! and the spin-orbit coupling contribution separately
  subroutine hamiltklinearsoc(sys,kp,hk,vsoc)
    use mod_kind, only: dp
    use mod_constants,  only: cZero, cI
    use mod_system,     only: ia, System_type
    use mod_parameters, only: nOrb,nOrb2
    use mod_magnet,     only: lb, sb
    use mod_SOC,        only: ls
    use mod_Umatrix,    only: hee
    implicit none
    integer :: i, j, k
    real(dp), intent(in) :: kp(3)
    type(System_type), intent(in) :: sys
    complex(dp),dimension(sys%nAtoms*nOrb2,sys%nAtoms*nOrb2),intent(out)  :: hk,vsoc
    complex(dp) :: tmp(nOrb, nOrb)
    complex(dp) :: kpExp

    hk = cZero
    vsoc = cZero

    ! Mouting slab hamiltonian

    ! On-site terms
    do concurrent (i=1:sys%nAtoms)
      ! spin-up on-site tight-binding term
      hk(ia(1,i):ia(2,i), ia(1,i):ia(2,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
      ! spin-down on-site tight-binding term
      hk(ia(3,i):ia(4,i), ia(3,i):ia(4,i)) = sys%Types(sys%Basis(i)%Material)%onSite(1:nOrb,1:nOrb)
      ! External magnetic field (orbital + spin) + Electron-electron interaction (Hubbard)
      hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = hk(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) &
                                           + lb(1:nOrb2,1:nOrb2,i) + sb(1:nOrb2,1:nOrb2,i) + hee(1:nOrb2,1:nOrb2,i)
      ! Spin-orbit coupling term
      vsoc(ia(1,i):ia(4,i), ia(1,i):ia(4,i)) = ls(1:nOrb2,1:nOrb2,i)
    end do

    ! Inter-site hopping terms
    do k = 1, sys%nNeighbors
      j = sys%Neighbors(k)%BasisIndex
      ! exp(ik.(R_i-R_j))
      kpExp = exp(cI * dot_product(kp,sys%Neighbors(k)%CellVector))

      do concurrent (i = 1:sys%nAtoms, sys%Neighbors(k)%isHopping(i))
        tmp(1:nOrb,1:nOrb) = sys%Neighbors(k)%t0i(1:nOrb, 1:nOrb, i) * kpExp
        ! Spin-up
        hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) = hk(ia(1,j):ia(2,j), ia(1,i):ia(2,i)) + tmp(1:nOrb,1:nOrb)
        ! Spin-down
        hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) = hk(ia(3,j):ia(4,j), ia(3,i):ia(4,i)) + tmp(1:nOrb,1:nOrb)
      end do
    end do
  end subroutine hamiltklinearsoc

end module mod_hamiltonian