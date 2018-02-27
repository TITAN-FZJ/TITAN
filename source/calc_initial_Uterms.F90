subroutine calc_initial_Uterms()
  !! Calculates the expectation value \sum_s U^i <c^+_imus c_inus>
  !! and U^i n_i/2
  use mod_f90_kind,      only: double
  use mod_constants,     only: cZero,pi
  use mod_System,        only: s => sys
  use TightBinding,      only: nOrb,nOrb2
  use EnergyIntegration, only: y, wght
  use mod_parameters,    only: eta
  use mod_magnet,        only: mzd,mpd,rhod,rho,rho0,rhod0
  use mod_Umatrix
  use adaptiveMesh
  use mod_mpi_pars
  implicit none
  integer      :: AllocateStatus
  integer*8    :: ix
  integer      :: i,mu,mup
  real(double) :: kp(3)
  real(double) :: weight, ep
  complex(double), dimension(:,:,:,:), allocatable :: gf
  real(double),    dimension(:,:),     allocatable :: imguu,imgdd
  !--------------------- begin MPI vars --------------------
  integer :: ncount
  ncount=s%nAtoms*nOrb
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  allocate(imguu(nOrb,s%nAtoms),imgdd(nOrb,s%nAtoms), stat = AllocateStatus)
  if(AllocateStatus/=0) call abortProgram("[calc_initial_Uterms] Not enough memory for: imguu,imgdd")
  allocate(rho0(nOrb,s%nAtoms),rhod0(s%nAtoms), stat = AllocateStatus)
  if(AllocateStatus/=0) call abortProgram("[calc_initial_Uterms] Not enough memory for: rho0,rhod0")

  mzd = 0.d0
  mpd = cZero
  do i = 1, s%nAtoms
    rhod(i) = s%Types(s%Basis(i)%Material)%OccupationD
    rhod0  (i) = s%Types(s%Basis(i)%Material)%OccupationD
  end do
  rho0 = 0.d0
  rho  = rho0
  call init_Umatrix(mzd,mpd,rhod,rhod0,rho,rho0,s%nAtoms,nOrb)

  imguu = 0.d0
  imgdd = 0.d0

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,i,mu,mup,kp,ep,weight,gf) &
  !$omp& shared(local_points,s,E_k_imag_mesh,bzs,eta,y,wght,imguu,imgdd)
  allocate(gf(nOrb2,nOrb2,s%nAtoms,s%nAtoms), stat = AllocateStatus)
  gf = cZero

  !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd)
  do ix = 1, local_points
    ep = y(E_k_imag_mesh(1,ix))
    kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
    weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
    !Green function on energy Ef + iy, and wave vector kp
    call green(s%Ef,ep+eta,kp,gf)
    do i=1,s%nAtoms
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
  do i=1,s%nAtoms
    do mu=1,nOrb
      imguu(mu,i) = 0.5d0 + imguu(mu,i)
      imgdd(mu,i) = 0.5d0 + imgdd(mu,i)
      if(mu>=5) rhod0(i) = rhod0(i) + real(imguu(mu,i) + imgdd(mu,i))
      rho0(mu,i) = imguu(mu,i) + imgdd(mu,i)
    end do
  end do

  ! if(rField == 0) write(*,*) rhod0,sum(abs(rho0))
  deallocate(imguu,imgdd)

  return
end subroutine calc_initial_Uterms











subroutine calc_initial_Uterms()

  use mod_System, only: System, s => sys
  use SK_TightBinding

  type(System),allocatable :: sys0(:)
  integer :: i

  allocate(sys0(nTypes))
  !------------ Allocating variables that depend on nAtoms -------------
  call allocate_magnet_variables(1, nOrb)
  call allocLS(1,nOrb)
  call allocate_Atom_variables(1) !TODO: Review

  do i=1,nTypes
    sys0(i)%nAtoms = 1
    allocate(sys0(i)%Basis(sys0(i)%nAtoms))
    sys0(i)%nTypes = 1
    allocate(sys0(i)%Types(sys0(i)%nTypes))
    sys0(i)%nStages = s%nStages

    sys0(i)%Basis(1)%Position(:) = [ 0.d0 , 0.d0 , 0.d0 ]
    sys0(i)%Basis(1)%Material = s%Basis(i)%Material

    sys0(i)%a0 = s%Types(s%Basis(i)%Material)%LatticeConstant
    sys0(i)%a1 = s%Types(s%Basis(i)%Material)%a1
    sys0(i)%a2 = s%Types(s%Basis(i)%Material)%a2
    sys0(i)%a3 = s%Types(s%Basis(i)%Material)%a3

    sys0(i)%Ef = s%Types(s%Basis(i)%Material)%FermiLevel

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

    ! allocate ( sys0(i)%Types(1)%onSite(size(s%Types(s%Basis(i)%Material)%onSite(:,:))) )
    ! allocate ( sys0(i)%Types(1)%Hopping(size(s%Types(s%Basis(i)%Material)%Hopping(:,:))) )
    ! allocate ( sys0(i)%Types(1)%Stage(size(s%Types(s%Basis(i)%Material)%Stage(:))) )
    ! sys0(i)%Types(1) = size(s%Types(s%Basis(i)%Material)%onSite(:,1))

    !------------------- Define the lattice structure --------------------
    call initLattice(sys0(i))
    !---- Generating k meshes points for imaginary axis integration ------
    call generateAdaptiveMeshes(sys0(i),pn1)


    !---------------------- Tight Binding parameters -----------------------
    call initTightBinding(sys0(i))

    !---- Initialize Stride Matrices for hamiltk and dtdksub --------------
    call initHamiltkStride(1, nOrb)



    ! Distribute Energy Integration across all points available
    call genLocalEKMesh(sys0(i),rField,sField, FieldComm)

    call calc_expectation_values(sys0(i))

    call freeLocalEKMesh()
  end do

  call deallocate_magnet_variables()


  return
end subroutine calc_initial_Uterms










subroutine calc_expectation_values(sys)
  !! Calculates the expectation value \sum_s U^i <c^+_imus c_inus>
  !! and U^i n_i/2
  use mod_f90_kind,      only: double
  use mod_constants,     only: cZero,tpi
  use mod_System,        only: System
  use TightBinding,      only: nOrb,nOrb2
  use EnergyIntegration, only: y, wght
  use mod_parameters,    only: eta, offset, U
  use mod_magnet,        only: mzd,mpd,rhod,cc,cc0,n0
  use mod_Umatrix
  use adaptiveMesh
  use mod_mpi_pars
  implicit none
  type(System), intent(in) :: sys
  integer      :: AllocateStatus
  integer*8    :: ix
  integer      :: i,mu,nu,mup,nup
  real(double) :: kp(3)
  real(double) :: weight, ep
  complex(double), dimension(:,:,:,:), allocatable :: gf
  complex(double), dimension(:,:,:),   allocatable :: imguu,imgdd
  !--------------------- begin MPI vars --------------------
  integer :: ncount
  ncount=sys%nAtoms*nOrb*nOrb
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  allocate(imguu(nOrb, nOrb,sys%nAtoms),imgdd(nOrb, nOrb,sys%nAtoms), stat = AllocateStatus)
  if(AllocateStatus/=0) call abortProgram("[calc_expectation_values] Not enough memory for: imguu,imgdd")
  allocate(cc(nOrb,nOrb,sys%nAtoms),cc0(nOrb,nOrb,sys%nAtoms),n0(sys%nAtoms), stat = AllocateStatus)
  if(AllocateStatus/=0) call abortProgram("[calc_expectation_values] Not enough memory for: cc,cc0,n0")

  mzd = 0.d0
  mpd = cZero
  do i = 1, sys%nAtoms
    rhod(i) = sys%Types(sys%Basis(i)%Material)%OccupationD
    n0 (i) = 0.5d0*U(i+offset) * sys%Types(sys%Basis(i)%Material)%OccupationD
  end do
  cc0 = cZero
  cc  = cc0
  ! Initialize U matrix with all terms zero
  call init_Umatrix(mzd,mpd,rhod,n0,cc,cc0,sys%nAtoms,nOrb)

  imguu = cZero
  imgdd = cZero

  !$omp parallel default(none) &
  !$omp& private(AllocateStatus,ix,i,mu,nu,mup,nup,kp,ep,weight,gf) &
  !$omp& shared(local_points,s,E_k_imag_mesh,bzs,eta,y,wght,imguu,imgdd)
  allocate(gf(nOrb2,nOrb2,sys%nAtoms,sys%nAtoms), stat = AllocateStatus)
  gf = cZero

  !$omp do schedule(static) reduction(+:imguu) reduction(+:imgdd)
  do ix = 1, local_points
    ep = y(E_k_imag_mesh(1,ix))
    kp = bzs(E_k_imag_mesh(1,ix)) % kp(:,E_k_imag_mesh(2,ix))
    weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix)) % w(E_k_imag_mesh(2,ix))
    !Green function on energy Ef + iy, and wave vector kp
    call green(sys%Ef,ep+eta,s,kp,gf)
    do i=1,sys%nAtoms
    do mu=1,nOrb
      mup = mu+nOrb
      do nu=1,nOrb
      nup = nu+nOrb

      imguu(mu,nu,i) = imguu(mu,nu,i) + ( gf(nu ,mu ,i,i) + conjg(gf(mu ,nu ,i,i)) ) * weight
      imgdd(mu,nu,i) = imgdd(mu,nu,i) + ( gf(nup,mup,i,i) + conjg(gf(mup,nup,i,i)) ) * weight
      end do
    end do
    end do
  end do
  !$omp end do

  deallocate(gf)
  !$omp end parallel
  imguu = imguu / tpi
  imgdd = imgdd / tpi

  call MPI_Allreduce(MPI_IN_PLACE, imguu, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, imgdd, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, activeComm, ierr)

  n0 = 0.d0
  do i=1,sys%nAtoms
    do mu=5,nOrb
      imguu(mu,mu,i) = 0.5d0 + imguu(mu,mu,i)
      imgdd(mu,mu,i) = 0.5d0 + imgdd(mu,mu,i)
      n0(i) = n0(i) + real(imguu(mu,mu,i) + imgdd(mu,mu,i))
    end do
    cc0(:,:,i) = (imguu(:,:,i) + imgdd(:,:,i))
  end do

  deallocate(imguu,imgdd)

  return
end subroutine calc_expectation_values