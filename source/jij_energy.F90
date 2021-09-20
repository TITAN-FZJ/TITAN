! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
! It uses the convention:
! E = - \sum_ij e_i J_ij e_j /2
! where e_i and e_j are unit vectors along the magnetization direction
! and the sum covers all sites.
subroutine jij_energy(q,Jij)
  use mod_kind,              only: dp,int32,int64
  use mod_constants,         only: pi,cOne,cZero
  use mod_parameters,        only: eta
  use EnergyIntegration,     only: y,wght
  use mod_magnet,            only: mdvec_cartesian,mabsd
  use mod_system,            only: s => sys
  use mod_hamiltonian,       only: hamilt_local
  use mod_greenfunction,     only: green
  use mod_SOC,               only: llinearsoc,llineargfsoc
  use mod_mpi_pars,          only: MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_SUM,ierr
  use adaptiveMesh,          only: activeComm,activeRank,local_points,bzs,E_k_imag_mesh
  use mod_superconductivity, only: lsupercond
  implicit none
  real(dp), dimension(3),                     intent(in)  :: q
  real(dp), dimension(s%nAtoms,s%nAtoms,3,3), intent(out) :: Jij
  real(dp), dimension(s%nAtoms,s%nAtoms,3,3)              :: Jijint
  integer(int64) :: ix
  integer(int32) :: i,j,mu,nu,alpha,nOrb2_i,nOrb2_j
  real(dp) :: kp(3),kminusq(3),ep
  real(dp) :: evec(3,s%nAtoms)
  real(dp) :: Jijk(s%nAtoms,s%nAtoms,3,3)
  real(dp) :: Jijkan(s%nAtoms,3,3)
  real(dp) :: weight, fermi
  integer(int32) :: po,mo
  complex(dp), dimension(s%nAtoms,3,s%nOrb2,s%nOrb2)        :: dBxc_dm
  complex(dp), dimension(s%nAtoms,3,3,s%nOrb2,s%nOrb2)      :: d2Bxc_dm2
  complex(dp), dimension(s%nAtoms,s%nOrb2,s%nOrb2)          :: paulievec
  complex(dp), dimension(s%nOrb2,s%nOrb2)                   :: gij,gji,temp1,temp2,paulia,paulib,pauli_star
  complex(dp), dimension(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms) :: gf,gfq
  complex(dp), dimension(s%nOrb2,s%nOrb2)                       :: geh, ghe
  integer :: ncount

  external :: zgemm,MPI_Reduce

  fermi = merge(0._dp,s%Ef,lsuperCond)

  ncount = s%nAtoms * s%nAtoms * 3 * 3

  do i = 1, s%nAtoms
    nOrb2_i = s%Types(s%Basis(i)%Material)%nOrb2

    ! Unit vector along the direction of the magnetization of each magnetic plane
    evec(1:3,i) = mdvec_cartesian(1:3,i)/mabsd(i)

    ! Inner product of pauli matrix in spin and orbital space and unit vector evec
    paulievec(i,1:nOrb2_i,1:nOrb2_i) = s%Types(s%Basis(i)%Material)%pauli_dorb(1,1:nOrb2_i,1:nOrb2_i) * evec(1,i) + s%Types(s%Basis(i)%Material)%pauli_dorb(2,1:nOrb2_i,1:nOrb2_i) * evec(2,i) + s%Types(s%Basis(i)%Material)%pauli_dorb(3,1:nOrb2_i,1:nOrb2_i) * evec(3,i)

    do mu = 1, 3
      ! Derivative of Bxc*sigma*evec w.r.t. m_i (Bxc = -U.m/2)
      dBxc_dm(i,mu,1:nOrb2_i,1:nOrb2_i) = -0.5_dp*s%Basis(i)%Um*(s%Types(s%Basis(i)%Material)%pauli_dorb(mu,1:nOrb2_i,1:nOrb2_i) - (paulievec(i,1:nOrb2_i,1:nOrb2_i)) * evec(mu,i))
      ! Second derivative of Bxc w.r.t. m_i (Bxc = -U.m/2)
      do nu=1,3
        d2Bxc_dm2(i,mu,nu,1:nOrb2_i,1:nOrb2_i) = evec(mu,i)*s%Types(s%Basis(i)%Material)%pauli_dorb(nu,1:nOrb2_i,1:nOrb2_i) + s%Types(s%Basis(i)%Material)%pauli_dorb(mu,1:nOrb2_i,1:nOrb2_i)*evec(nu,i) - 3*paulievec(i,1:nOrb2_i,1:nOrb2_i)*evec(mu,i)*evec(nu,i)
        if(mu==nu) d2Bxc_dm2(i,mu,nu,1:nOrb2_i,1:nOrb2_i) = d2Bxc_dm2(i,mu,nu,1:nOrb2_i,1:nOrb2_i) + paulievec(i,1:nOrb2_i,1:nOrb2_i)
        d2Bxc_dm2(i,mu,nu,1:nOrb2_i,1:nOrb2_i) = 0.5_dp*s%Basis(i)%Um*d2Bxc_dm2(i,mu,nu,1:nOrb2_i,1:nOrb2_i)/mabsd(i)
      end do
    end do
  end do

  ! Build local hamiltonian
  ! call hamilt_local(s)
  if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

  ! Calculating Jij using a complex integral
  Jij = 0._dp

  !$omp parallel default(none) &
  !$omp& private(ix,i,j,po,mo,mu,nu,nOrb2_i,nOrb2_j,alpha,kp,ep,weight,kminusq,gf,gfq,gij,gji,paulia,paulib,pauli_star,temp1,temp2,Jijkan,Jijk,Jijint,geh,ghe) &
  !$omp& shared(s,fermi,lsuperCond,eta,y,wght,q,local_points,dBxc_dm,d2Bxc_dm2,Jij,bzs,E_k_imag_mesh,paulievec)

  Jijint = 0._dp

  !$omp do schedule(dynamic)
  do ix = 1, local_points
      ep = y( E_k_imag_mesh(1,ix) )
      kp = bzs( E_k_imag_mesh(1,ix) ) % kp(1:3,E_k_imag_mesh(2,ix))
      weight = wght(E_k_imag_mesh(1,ix)) * bzs(E_k_imag_mesh(1,ix))%w(E_k_imag_mesh(2,ix))
      ! Green function on energy Ef + iy, and wave vector kp
      call green(fermi,ep+eta,s,kp,gf)

      ! Green function on energy Ef + iy, and wave vector kp-q
      kminusq = kp-q
      call green(fermi,ep+eta,s,kminusq,gfq)

      Jijk   = 0._dp
      Jijkan = 0._dp
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          nOrb2_j = s%Types(s%Basis(j)%Material)%nOrb2
          nOrb2_i = s%Types(s%Basis(i)%Material)%nOrb2
          gij = gf (1:nOrb2_i,1:nOrb2_j,i,j)
          gji = gfq(1:nOrb2_j,1:nOrb2_i,j,i)
          if (lsuperCond) then
            geh(1:nOrb2_i,1:nOrb2_j) = gf (1:nOrb2_i,nOrb2_j+1:2*nOrb2_j,i,j)
            ghe(1:nOrb2_j,1:nOrb2_i) = gfq(nOrb2_j+1:2*nOrb2_j,1:nOrb2_i,j,i)
          end if

          do nu = 1,3
            do mu = 1,3
              paulia(1:nOrb2_i,1:nOrb2_i)= dBxc_dm(i,mu,1:nOrb2_i,1:nOrb2_i)
              paulib(1:nOrb2_j,1:nOrb2_j)= dBxc_dm(j,nu,1:nOrb2_j,1:nOrb2_j)
              ! temp1 = paulia * gij = dBxc_dm_i * gij
              call zgemm('n','n',nOrb2_i,nOrb2_j,nOrb2_i,cOne,paulia,s%nOrb2,gij,   s%nOrb2,cZero,temp1,s%nOrb2)
              ! temp2 = temp1 * paulib = paulia * gij * paulib = dBxc_dm_i * gij * dBxc_dm_j
              call zgemm('n','n',nOrb2_i,nOrb2_j,nOrb2_j,cOne,temp1, s%nOrb2,paulib,s%nOrb2,cZero,temp2,s%nOrb2)
              ! temp1 = temp2 * gji = paulia * gij * paulib * gji = dBxc_dm_i * gij * dBxc_dm_j * gji
              call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,cOne,temp2, s%nOrb2,gji,   s%nOrb2,cZero,temp1,s%nOrb2)

              ! Trace over orbitals and spins
              do alpha = 1,nOrb2_i
                Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
              end do
              if ( lsuperCond ) then
                pauli_star(1:nOrb2_j,1:nOrb2_j) = conjg(paulib(1:nOrb2_j,1:nOrb2_j))
                ! temp1 = paulia * geh(ij) = dBxc_dm_i * geh(ij)
                call zgemm('n','n',nOrb2_i,nOrb2_j,nOrb2_i,cOne,paulia,s%nOrb2,geh       ,s%nOrb2,cZero,temp1,s%nOrb2)
                ! temp2 = temp1 * pauli_star = dBxc_dm_i * geh(ij) * (dBxc_dm_j)^*
                call zgemm('n','n',nOrb2_i,nOrb2_j,nOrb2_j,cOne,temp1 ,s%nOrb2,pauli_star,s%nOrb2,cZero,temp2,s%nOrb2)
                ! temp1 = temp2 * ghe(ji) = dBxc_dm_i * geh(ij) * (dBxc_dm_j)^* * ghe(ji)
                call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_j,cOne,temp2 ,s%nOrb2,ghe       ,s%nOrb2,cZero,temp1,s%nOrb2)
                do alpha = 1,nOrb2_i
                  Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
                end do
              end if
            end do
          end do

          ! Anisotropy (on-site) term
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(i==j) then
            gij = gf(1:nOrb2_i,1:nOrb2_i,i,i)
            do nu = 1,3
              do mu = 1,3
                paulia(1:nOrb2_i,1:nOrb2_i) = d2Bxc_dm2(i,mu,nu,1:nOrb2_i,1:nOrb2_i)
                ! temp1 = gii * d2Bxc_dm2_i
                call zgemm('n','n',nOrb2_i,nOrb2_i,nOrb2_i,cOne,gij,s%nOrb2,paulia,s%nOrb2,cZero,temp1,s%nOrb2)
                ! Trace over orbitals and spins
                do alpha = 1,nOrb2_i
                  Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
                end do
                Jijk(i,i,mu,nu) = Jijk(i,i,mu,nu) + Jijkan(i,mu,nu)
              end do
            end do
          end if
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  !call MPI_Allreduce(MPI_IN_PLACE, Jij, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, activeComm, ierr)

end subroutine jij_energy
