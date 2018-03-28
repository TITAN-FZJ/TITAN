module expectation

contains

  subroutine calc_initial_Uterms(sys)
    use mod_f90_kind,   only: double
    use mod_constants,  only: cZero
    use mod_System,     only: System,initHamiltkStride
    use TightBinding,   only: nOrb,initTightBinding
    use mod_magnet,     only: l_matrix,lb,sb,allocate_magnet_variables,deallocate_magnet_variables,rho0,rhod0
    use mod_SOC,        only: ls,allocLS
    use adaptiveMesh,   only: generateAdaptiveMeshes,genLocalEKMesh,freeLocalEKMesh
    use mod_parameters, only: kp_in,output
    use mod_polyBasis,  only: polyBasis => read_basis
    use mod_mpi_pars,   only: myrank,FieldComm,rField,sField,abortProgram
    use Lattice,        only: initLattice
    use mod_progress,   only: write_time
    ! use SK_TightBinding
    use mod_BrillouinZone, only: realBZ!,nkpt_x,nkpt_y,nkpt_z
    use EnergyIntegration, only: pn1
    implicit none
    integer :: i,j,mu
    type(System), intent(inout) :: sys
    type(System), allocatable   :: sys0(:)

    if(myrank == 0) &
    call write_time(output%unit,'[calc_initial_Uterms] Start calculating initial density: ')

    allocate(sys0(sys%nTypes))

    allocate(rho0(nOrb,sys%nAtoms),rhod0(sys%nAtoms))

    do i = 1, sys%nTypes
      !------------------ Define the lattice structure -------------------
      call polyBasis(trim(sys%Types(i)%Name), sys0(i))
      if(sys0(i)%nAtoms/=1) call abortProgram("[calc_initial_Uterms] Not implemented for parameter file with more than 1 atom!")
      ! Setting the number of nearest neighbors
      sys0(i)%nStages = sys%nStages
      sys0(i)%relTol  = sys%relTol
      !---------- Generating k points for real axis integration ----------
      if(dot_product(sys0(i)%a3,sys0(i)%a3) == 0.d0) then
        sys0(i)%lbulk = .false.
        realBZ % nkpt_x = kp_in(1)
        realBZ % nkpt_y = kp_in(1)
        realBZ % nkpt_z = 0
      else
        sys0(i)%lbulk = .true.
        realBZ % nkpt_x = kp_in(1)
        realBZ % nkpt_y = kp_in(1)
        realBZ % nkpt_z = kp_in(1)
      end if
      call realBZ % count(sys0(i))

      call initLattice(sys0(i))

! if(myrank==0) call writeLattice(sys0(i))
! stop
      !--- Generating k meshes points for imaginary axis integration -----
      call generateAdaptiveMeshes(sys0(i),pn1)

      !----------- Allocating variables that depend on nAtoms ------------
      call allocate_magnet_variables(sys0(i)%nAtoms, nOrb)
      call allocLS(sys0(i)%nAtoms,nOrb)

      !-------------------- Tight Binding parameters ---------------------
      call initTightBinding(sys0(i))

      !------- Initialize Stride Matrices for hamiltk and dtdksub --------
      call initHamiltkStride(sys0(i)%nAtoms, nOrb)

      !---- L matrix in global frame for given quantization direction ----
      call l_matrix()

      !----- Removing L.B, S.B and L.S matrices from the hamiltonian -----
      lb = cZero
      sb = cZero
      ls = cZero

      ! Distribute Energy Integration across all points available
      call genLocalEKMesh(sys0(i),rField,sField, FieldComm)

      !---------------- Calculating expectation values -------------------
      call calc_expectation_values(sys0(i),sys%Types(i)%rho0,sys%Types(i)%rhod0)

      !------------------------- Test printing ---------------------------
      if(myrank == 0) then
        write(output%unit,"('[calc_initial_Uterms] Calculating initial densities from: ',a)") trim(sys0(i)%Name)
        write(output%unit,"(a,12(2x,es16.9))") trim(sys%Types(i)%Name) , ((sys%Types(i)%rho0(mu,j),mu=1,nOrb),j=1,sys0(i)%nAtoms)
        write(output%unit,"(a,12(2x,es16.9))") trim(sys%Types(i)%Name) , (sys%Types(i)%rhod0(j),j=1,sys0(i)%nAtoms)
      end if

      !------------------------ Freeing memory ---------------------------
      call freeLocalEKMesh()

      !--------------------- Deallocating variables ------------------------
      call deallocate_magnet_variables()
    end do

    ! Transfering from occupations stored on Type to variables used in the hamiltonian
    do i=1,sys%nAtoms
      rho0(:,i) = sys%Types(sys%Basis(i)%Material)%rho0(:,1)
      rhod0(i)  = sys%Types(sys%Basis(i)%Material)%rhod0(1)
    end do

    if(myrank == 0) &
    call write_time(output%unit,'[calc_initial_Uterms] Finished calculating initial density: ')

  end subroutine calc_initial_Uterms


  subroutine calc_expectation_values(sys,rho0,rhod0)
    !! Calculates the expectation values of n_mu^s and n_i/2
    use mod_f90_kind,      only: double
    use mod_constants,     only: cZero,pi
    use mod_System,        only: System
    use TightBinding,      only: nOrb,nOrb2
    use EnergyIntegration, only: y, wght
    use mod_parameters,    only: eta
    use mod_magnet,        only: mzd,mpd,rhod,rho
    use mod_Umatrix
    use adaptiveMesh
    use mod_mpi_pars
    implicit none
    integer      :: AllocateStatus
    integer*8    :: ix
    integer      :: i,mu,mup
    real(double) :: kp(3)
    real(double) :: weight, ep
    type(System)                                     :: sys
    real(double),    dimension(:,:)    , allocatable,intent(out) :: rho0
    real(double),    dimension(:)      , allocatable,intent(out) :: rhod0
    real(double),    dimension(:,:)    , allocatable :: imguu,imgdd
    complex(double), dimension(:,:,:,:), allocatable :: gf
    !--------------------- begin MPI vars --------------------
    integer :: ncount
    ncount=sys%nAtoms*nOrb
    !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

    if(.not.allocated(rho0 )) allocate(rho0(nOrb,sys%nAtoms))
    if(.not.allocated(rhod0)) allocate(rhod0(sys%nAtoms))

    allocate(imguu(nOrb,sys%nAtoms),imgdd(nOrb,sys%nAtoms), stat = AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[calc_expectation_values] Not enough memory for: imguu,imgdd")

    mzd = 0.d0
    mpd = cZero
    do i = 1, sys%nAtoms
      rhod (i) = sys%Types(sys%Basis(i)%Material)%OccupationD
      rhod0(i) = sys%Types(sys%Basis(i)%Material)%OccupationD
    end do
    rho0 = 0.d0
    rho  = rho0
    call init_Umatrix(mzd,mpd,rhod,rhod0,rho,rho0,sys%nAtoms,nOrb)

    imguu = 0.d0
    imgdd = 0.d0

    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,ix,i,mu,mup,kp,ep,weight,gf) &
    !$omp& shared(local_points,sys,E_k_imag_mesh,bzs,eta,y,wght,imguu,imgdd)
    allocate(gf(nOrb2,nOrb2,sys%nAtoms,sys%nAtoms), stat = AllocateStatus)
    gf = cZero

    !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd)
    do ix = 1, local_points
      ep = y(E_k_imag_mesh(1,ix))
      kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
      weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
      !Green function on energy Ef + iy, and wave vector kp
      call green(sys%Ef,ep+eta,sys,kp,gf)
      do i=1,sys%nAtoms
        do mu=1,nOrb
          mup = mu+nOrb
          imguu(mu,i) = imguu(mu,i) + real(gf(mu ,mu ,i,i)) * weight
          imgdd(mu,i) = imgdd(mu,i) + real(gf(mup,mup,i,i)) * weight
        end do
      end do
    end do
    !$omp end do

    deallocate(gf)
    !$omp end parallel
    imguu = imguu / pi
    imgdd = imgdd / pi

    call MPI_Allreduce(MPI_IN_PLACE, imguu, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, imgdd, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)

    rhod0 = 0.d0
    do i=1,sys%nAtoms
      do mu=1,nOrb
        imguu(mu,i) = 0.5d0 + imguu(mu,i)
        imgdd(mu,i) = 0.5d0 + imgdd(mu,i)
        rho0(mu,i)  = imguu(mu,i) + imgdd(mu,i)
        if(mu>=5) rhod0(i) = rhod0(i) + rho0(mu,i)
      end do
    end do

    ! if(rField == 0) write(*,*) rhod0,sum(abs(rho0))
    deallocate(imguu,imgdd)

  end subroutine calc_expectation_values
end module expectation