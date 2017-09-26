! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
subroutine jij_energy(Jij)
  use mod_f90_kind, only: double
  use mod_constants, only: pi, cOne, cZero, pauli_dorb
  use mod_parameters, only: mmlayermag, U, q, mmlayermag, outputunit, nmaglayers, Ef, outputunit
  use EnergyIntegration, only: pn1, y, wght
  use mod_mpi_pars
  use mod_magnet, only: mx,my,mz,mabs
  use mod_system, only: s => sys
  use mod_BrillouinZone, only: BZ
  use TightBinding, only: nOrb,nOrb2


  implicit none
  real(double),dimension(nmaglayers,nmaglayers,3,3) :: Jijint
  real(double),dimension(nmaglayers,nmaglayers,3,3),intent(out) :: Jij

  integer :: ix,iz
  integer :: i,j,mu,nu,alpha
  real(double) :: kp(3), kminusq(3)
  real(double) :: evec(3,nmaglayers)
  real(double) :: Jijk(nmaglayers,nmaglayers,3,3)
  real(double) :: Jijkan(nmaglayers,3,3)
  complex(double), dimension(nmaglayers,3,nOrb2,nOrb2) :: dbxcdm
  complex(double), dimension(nmaglayers,3,3,nOrb2,nOrb2) :: d2bxcdm2
  complex(double), dimension(nmaglayers,nOrb2,nOrb2) :: paulievec
  complex(double), dimension(nOrb2,nOrb2) :: gij, gji, temp1, temp2, paulia, paulib
  complex(double), dimension(nOrb2,nOrb2,s%nAtoms,s%nAtoms) :: gf,gfq
  complex(double) :: weight
  !--------------------- begin MPI vars --------------------
  integer :: start, end, work, remainder
  integer :: ncount
  !^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ! Calculate workload for each MPI process
  if(numprocs <= pn1) then
    remainder = mod(pn1,numprocs)
    if(myrank < remainder) then
      work = ceiling(dble(pn1) / dble(numprocs))
      start = myrank*work + 1
      end = (myrank+1) * work
    else
      work = floor(dble(pn1) / dble(numprocs))
      start = myrank*work + 1 + remainder
      end = (myrank+1) * work + remainder
    end if
  else
    write(outputunit, "('[jij_energy] Too many MPI threads. ',i0,' threads are idle.')") numprocs - pn1
    if(myrank < pn1) then
      start = myrank + 1
      end = myrank + 1
    else
      start = 0
      end = 0
    end if
  end if

  ncount = nmaglayers * nmaglayers * 3 * 3



  do iz = 1, nmaglayers
    ! Unit vector along the direction of the magnetization of each magnetic plane
    evec(:,iz) = [ mx(mmlayermag(iz)), my(mmlayermag(iz)), mz(mmlayermag(iz)) ]/mabs(mmlayermag(iz))

    ! Inner product of pauli matrix in spin and orbital space and unit vector evec
    paulievec(iz,:,:) = pauli_dorb(1,:,:) * evec(1,iz) + pauli_dorb(2,:,:) * evec(2,iz) + pauli_dorb(3,:,:) * evec(3,iz)

    do i = 1, 3
      ! Derivative of Bxc*sigma*evec w.r.t. m_i (Bxc = -U.m/2)
      dbxcdm(iz,i,:,:) = -0.5d0 * U(mmlayermag(iz)) * (pauli_dorb(i,:,:) - (paulievec(iz,:,:)) * evec(i,iz))

      ! Second derivative of Bxc w.r.t. m_i (Bxc = -U.m/2)
      do j=1,3
        d2bxcdm2(iz,i,j,:,:) = evec(i,iz)*pauli_dorb(j,:,:) + pauli_dorb(i,:,:)*evec(j,iz) - 3*paulievec(iz,:,:)*evec(i,iz)*evec(j,iz)
        if(i==j) d2bxcdm2(iz,i,j,:,:) = d2bxcdm2(iz,i,j,:,:) + paulievec(iz,:,:)
        d2bxcdm2(iz,i,j,:,:) = 0.5d0*U(mmlayermag(iz))*d2bxcdm2(iz,i,j,:,:)/(mabs(mmlayermag(iz)))
      end do
    end do
  end do

  ! Calculating the number of particles for each spin and orbital using a complex integral
  Jij = 0.d0

  !$omp parallel default(none) &
  !$omp& private(ix,iz,i,j,mu,nu,alpha,kp,weight,kminusq,gf,gfq,gij,gji,paulia,paulib,temp1,temp2,Jijkan,Jijk,Jijint) &
  !$omp& shared(s,BZ,Ef,y,wght,q,start,end,dbxcdm,d2bxcdm2,mz,nmaglayers,mmlayermag,Jij)

  Jijint = 0.d0

  !$omp do schedule(static) collapse(2)
  do ix = start, end
    do iz = 1, BZ%nkpt
      kp = BZ%kp(1:3,iz)
      weight = BZ%w(iz) * wght(ix)
      ! Green function on energy Ef + iy, and wave vector kp
      call green(Ef,y(ix),kp,gf)

      ! Green function on energy Ef + iy, and wave vector kp-q
      kminusq = kp-q
      call green(Ef,y(ix),kminusq,gfq)

      Jijk   = 0.d0
      Jijkan = 0.d0
      do j = 1,nmaglayers
        do i = 1,nmaglayers
          gij = gf(:,:,mmlayermag(i),mmlayermag(j))
          gji = gfq(:,:,mmlayermag(j),mmlayermag(i))
          do nu = 1,3
            do mu = 1,3
              paulia = dbxcdm(i,mu,:,:)
              paulib = dbxcdm(j,nu,:,:)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,paulia,nOrb2,gij,   nOrb2,cZero,temp1,nOrb2)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp1, nOrb2,paulib,nOrb2,cZero,temp2,nOrb2)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp2, nOrb2,gji,   nOrb2,cZero,temp1,nOrb2)
              ! Trace over orbitals and spins
              do alpha = 1,nOrb2
                Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
              end do
            end do
          end do

          ! Anisotropy (on-site) term
          if(i==j) then
            gij = gf(:,:,mmlayermag(i),mmlayermag(i))
            do nu = 1,3
              do mu = 1,3
                paulia = d2bxcdm2(i,mu,nu,:,:)
                call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,gij,nOrb2,paulia,nOrb2,cZero,temp1,nOrb2)
                ! Trace over orbitals and spins
                do alpha = 1,nOrb2
                  Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
                end do
                Jijk(i,i,mu,nu) = Jijk(i,i,mu,nu) + Jijkan(i,mu,nu)
              end do
            end do
          end if
        end do
      end do

      Jijint = Jijint + Jijk * weight
    end do
  end do
  !$omp end do nowait

  Jijint = -Jijint/pi
  !$omp critical
    Jij = Jij + Jijint
  !$omp end critical
  !$omp end parallel


  call MPI_Allreduce(MPI_IN_PLACE, Jij, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  !Jij = Jij / pi !TODO: Check with filipe if wrong
  return
  end subroutine jij_energy
