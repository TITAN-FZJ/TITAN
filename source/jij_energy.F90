! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
! It uses the convention:
! E = - \sum_ij e_i J_ij e_j /2
! where e_i and e_j are unit vectors along the magnetization direction
! and the sum covers all sites.
subroutine jij_energy(Jij)
  use mod_kind,          only: dp,int32,int64
  use mod_constants,     only: pi, cOne, cZero, pauli_dorb
  use mod_parameters,    only: nOrb2, Um, q, eta
  use EnergyIntegration, only: y, wght
  use mod_magnet,        only: mvec_cartesian,mabs
  use mod_system,        only: s => sys
  use mod_hamiltonian,   only: hamilt_local
  use mod_greenfunction, only: green
  use mod_mpi_pars,      only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_SUM,ierr
  use adaptiveMesh,      only: activeComm,activeRank,local_points,bzs,E_k_imag_mesh
  implicit none
  real(dp),dimension(s%nAtoms,s%nAtoms,3,3)              :: Jijint
  real(dp),dimension(s%nAtoms,s%nAtoms,3,3), intent(out) :: Jij
  integer(int64) :: ix
  integer(int32) :: i,j,mu,nu,alpha
  real(dp) :: kp(3), kminusq(3), ep
  real(dp) :: evec(3,s%nAtoms)
  real(dp) :: Jijk(s%nAtoms,s%nAtoms,3,3)
  real(dp) :: Jijkan(s%nAtoms,3,3)
  real(dp) :: weight
  complex(dp), dimension(s%nAtoms,3,nOrb2,nOrb2)        :: dBxc_dm
  complex(dp), dimension(s%nAtoms,3,3,nOrb2,nOrb2)      :: d2Bxc_dm2
  complex(dp), dimension(s%nAtoms,nOrb2,nOrb2)          :: paulievec
  complex(dp), dimension(nOrb2,nOrb2)                   :: gij, gji, temp1, temp2, paulia, paulib
  complex(dp), dimension(nOrb2,nOrb2,s%nAtoms,s%nAtoms) :: gf,gfq
  integer :: ncount

  external :: zgemm,MPI_Reduce

  ncount = s%nAtoms * s%nAtoms * 3 * 3

  do i = 1, s%nAtoms
    ! Unit vector along the direction of the magnetization of each magnetic plane
    evec(:,i) = [ mvec_cartesian(1,i), mvec_cartesian(2,i), mvec_cartesian(3,i) ]/mabs(i)

    ! Inner product of pauli matrix in spin and orbital space and unit vector evec
    paulievec(i,:,:) = pauli_dorb(1,:,:) * evec(1,i) + pauli_dorb(2,:,:) * evec(2,i) + pauli_dorb(3,:,:) * evec(3,i)

    do mu = 1, 3
      ! Derivative of Bxc*sigma*evec w.r.t. m_i (Bxc = -U.m/2)
      dBxc_dm(i,mu,:,:) = -0.5_dp*Um(i)*(pauli_dorb(mu,:,:) - (paulievec(i,:,:)) * evec(mu,i))
      ! Second derivative of Bxc w.r.t. m_i (Bxc = -U.m/2)
      do nu=1,3
        d2Bxc_dm2(i,mu,nu,:,:) = evec(mu,i)*pauli_dorb(nu,:,:) + pauli_dorb(mu,:,:)*evec(nu,i) - 3*paulievec(i,:,:)*evec(mu,i)*evec(nu,i)
        if(mu==nu) d2Bxc_dm2(i,mu,nu,:,:) = d2Bxc_dm2(i,mu,nu,:,:) + paulievec(i,:,:)
        d2Bxc_dm2(i,mu,nu,:,:) = 0.5_dp*Um(i)*d2Bxc_dm2(i,mu,nu,:,:)/mabs(i)
      end do
    end do
  end do

  ! Build local hamiltonian
  call hamilt_local(s)

  ! Calculating the number of particles for each spin and orbital using a complex integral
  Jij = 0._dp

  !$omp parallel default(none) &
  !$omp& private(ix,i,j,mu,nu,alpha,kp,ep,weight,kminusq,gf,gfq,gij,gji,paulia,paulib,temp1,temp2,Jijkan,Jijk,Jijint) &
  !$omp& shared(s,nOrb2,eta,y,wght,q,local_points,dBxc_dm,d2Bxc_dm2,Jij,bzs,E_k_imag_mesh)

  Jijint = 0._dp

  !$omp do schedule(static)
  do ix = 1, local_points
     ep = y( E_k_imag_mesh(1,ix) )
     kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
     weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix))
     ! Green function on energy Ef + iy, and wave vector kp
      call green(s%Ef,ep+eta,s,kp,gf)

      ! Green function on energy Ef + iy, and wave vector kp-q
      kminusq = kp-q
      call green(s%Ef,ep+eta,s,kminusq,gfq)

      Jijk   = 0._dp
      Jijkan = 0._dp
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          gij = gf(:,:,i,j)
          gji = gfq(:,:,j,i)
          do nu = 1,3
            do mu = 1,3
              paulia = dBxc_dm(i,mu,:,:)
              paulib = dBxc_dm(j,nu,:,:)
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
            gij = gf(:,:,i,i)
            do nu = 1,3
              do mu = 1,3
                paulia = d2Bxc_dm2(i,mu,nu,:,:)
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
  !$omp end do nowait

  Jijint = -Jijint/pi
  !$omp critical
    Jij = Jij + Jijint
  !$omp end critical
  !$omp end parallel

  if(activeRank == 0) then
     call MPI_Reduce(MPI_IN_PLACE, Jij, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, 0, activeComm, ierr)
  else
     call MPI_Reduce(Jij, Jij, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, 0, activeComm, ierr)
  end if

end subroutine jij_energy
