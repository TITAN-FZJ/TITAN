! This is the main subroutine to calculate all quantities:
! currents, disturbances, torques, effective fields and susceptibilities
subroutine calculate_all()
  use mod_kind,              only: dp
  use mod_constants,         only: cZero,cOne,cI !, levi_civita
  use mod_parameters,        only: lnodiag,output,laddresults,skip_steps,kount,lhfresponses,&
                                  sigmaimunu2i,emin,deltae,nQvec1,kpoints,dimens,sigmai2i,dimspinAtoms
  use mod_magnet,            only: lfield,hhw,mvec_cartesian,mvec_spherical,lrot
  use mod_SOC,               only: llinearsoc
  use mod_System,            only: s => sys
  use mod_BrillouinZone,     only: realBZ
  use adaptiveMesh,          only: genLocalEKMesh,freeLocalEKMesh
  use mod_prefactors,        only: prefactor,prefactorlsoc,allocate_prefactors,deallocate_prefactors
  use mod_susceptibilities,  only: rottemp,rotmat_i,rotmat_j,&
                                   schitemp,schirot,schi,schihf,&
                                   chiorb,chiorb_hf,chiorb_hflsoc,Umatorb,identt,&
                                   build_identity_and_U_matrix,diagonalize_susceptibilities,&
                                   allocate_susceptibilities,deallocate_susceptibilities,&
                                   create_chi_files,write_susceptibilities
  use mod_torques,           only: torques,total_torques,ntypetorque,&
                                    allocate_torques,deallocate_torques,&
                                    create_torque_files,write_torques
  use mod_disturbances,      only: disturbances,total_disturbances,tchiorbiikl,sdmat,ldmat,&
                                   allocate_disturbances,deallocate_disturbances,&
                                   create_disturbance_files,write_disturbances
  use mod_beff,              only: chiinv,total_Beff,Beff,Beff_cart,&
                                   allocate_beff,deallocate_beff,&
                                   create_beff_files,write_beff
  use mod_progress,          only: write_time
  use mod_rotation_matrices, only: rotation_matrices_chi
  use mod_sumrule,           only: sumrule
  use mod_mpi_pars,          only: rFreq,sFreq,FreqComm,FieldComm,rField,startFreq,endFreq,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,stat,ierr,MPI_SOURCE,MPI_DOUBLE_COMPLEX
  use mod_check_stop,        only: check_stop
  use mod_tools,             only: invers,itos
  use adaptiveMesh,          only: bzs
  !use mod_system, only: n0sc1, n0sc2, n0sc
  !use mod_currents !TODO: Re-Include
  !use mod_sha !TODO: Re-Include
  !use mod_diamagnetic_current
  implicit none
  integer  :: mcount,qcount
  integer  :: i,j,sigma,sigmap,mu,nu,mud,nud!,neighbor
  real(dp) :: e,q(3)!,Icabs

  external :: eintshechi,eintshelinearsoc,eintshechilinearsoc,eintshe,zgemm,MPI_Recv,MPI_Send,MPI_Barrier,sort_all_files

  call allocate_prefactors()
  call allocate_susceptibilities()
  call allocate_disturbances()
  call allocate_beff()
  call allocate_torques()
  ! call allocate_currents() !TODO: Re-Include
  ! call allocate_sha() !TODO: Re-Include
  ! call allocate_idia()

  if(rFreq(1) == 0) then
    write(output%unit_loop,"('CALCULATING PARALLEL CURRENTS, DISTURBANCES, LOCAL SUSCEPTIBILITY, EFFECTIVE FIELDS, SO-TORQUES AND V_DC AS A FUNCTION OF ENERGY')")
    ! Creating files and writing headers
    if(.not.laddresults) then
      call create_chi_files()
      call create_disturbance_files()
      call create_beff_files()
      call create_torque_files()
      ! call create_currents_files() !TODO: Re-Include
      ! call create_sha_files() !TODO: Re-Include
    end if
  end if

  ! Mounting U and identity matrix
  call build_identity_and_U_matrix()

  call genLocalEKMesh(s,rFreq(1),sFreq(1),FreqComm(1),bzs)
  call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))

  ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
  !call OAM_curr_hopping_times_L() !TODO:Re-Include

  ! ! Calculating diamagnetic current contribution in the first 'pn1' processes (final result is on myrank = 0 only)
  ! call calculate_idia()
  ! call MPI_Finalize(ierr)
  ! stop

  if(rField == 0 .and. skip_steps > 0) &
    write(output%unit,"('[calculate_all] Skipping first ',i0,' step(s)...')") skip_steps

  ! Chi wave vector Loop
  do qcount=1,nQvec1
    if(rField==0) &
      write(output%unit_loop,"('[calculate_all] Wave vector Q loop: ',i0,' of ',i0,' points',', Q = [',es10.3,es10.3,es10.3,']')") qcount,nQvec1,(kpoints(i,qcount),i=1,3)
    q = kpoints(:,qcount)
    ! Responses Energy Loop
    do kount = startFreq+skip_steps, endFreq
      e = emin + deltae * (kount-1)
      if(rField==0) &
        write(output%unit_loop,"('[calculate_chi] Starting MPI step ',i0,' of ',i0,':',10x,'E =  ',es10.3)") kount - startFreq - skip_steps + 1, endFreq - startFreq + 1, e

      if(lhfresponses) then
        if(rFreq(1) == 0) &
          write(output%unit_loop,"('[calculate_all] No renormalization will be done. Setting prefactors to identity and calculating HF susceptibilities... ')")
        prefactor = identt
        if(llinearsoc) prefactorlsoc = identt
        call eintshechi(q,e)
        if(rField == 0) &
          call write_time('[calculate_all] Time after susceptibility calculation: ',output%unit_loop)
      else
        if(rField == 0) &
          write(output%unit_loop,"('[calculate_all] Calculating prefactor to use in currents and disturbances calculation. ')")
        if(llinearsoc) then
          call eintshechilinearsoc(q,e) ! Note: chiorb_hflsoc = lambda*dchi_hf/dlambda(lambda=0)
        else
          call eintshechi(q,e)
        end if

        ! Checking sum rule for e=0._dp
        if((abs(e) < 1.e-8_dp).and.(sum(abs(q)) < 1.e-8_dp)) call sumrule(chiorb_hf)

        ! prefactor = (1 + chi_hf*Umat)^-1
        prefactor     = identt
        call zgemm('n','n',dimens,dimens,dimens,cOne,chiorb_hf,dimens,Umatorb,dimens,cOne,prefactor,dimens) ! prefactor = 1+chi_hf*Umat
        call invers(prefactor,dimens)
        if(llinearsoc) then
          prefactorlsoc = identt
          if (.not. allocated(chiorb)) allocate(chiorb(dimens,dimens))
          chiorb = chiorb_hf-chiorb_hflsoc ! the array chiorb will be used as a temporary array
          call zgemm('n','n',dimens,dimens,dimens,cOne,chiorb,dimens,Umatorb,dimens,cOne,prefactorlsoc,dimens) ! prefactorlsoc = 1+chiorb*Umat = 1+(chi_hf + dchi_hf/dlambda)*Umat
          call zgemm('n','n',dimens,dimens,dimens,cOne,prefactor,dimens,prefactorlsoc,dimens,cZero,chiorb,dimens) ! chiorb = prefactor*prefactorlsoc
          call zgemm('n','n',dimens,dimens,dimens,cOne,chiorb,dimens,prefactor,dimens,cZero,prefactorlsoc,dimens) ! prefactorlsoc = chiorb*prefactor = prefactor*prefactorlsoc*prefactor
        end if
        if(rField == 0) &
          call write_time('[calculate_all] Time after prefactor calculation: ',output%unit_loop)
      end if

      ! Start parallelized processes to calculate all quantities for energy e
      if(llinearsoc) then
        call eintshelinearsoc(q,e)
      else
        call eintshe(q,e)
      end if

      if(rField == 0) &
        call write_time('[calculate_all] Time after energy integral: ',output%unit_loop)

      if(rFreq(1) == 0) then

        ! Calculating the full matrix of RPA and HF susceptibilities for energy e
        if(.not.lhfresponses) then

          ! Calculate RPA susceptibility
          if(llinearsoc) then
            ! chiorb = prefactorlsoc*chi_hf + prefactor*chi_hflsoc
            call zgemm('n','n',dimens,dimens,dimens,cOne,prefactorlsoc,dimens,chiorb_hf,dimens,cZero,chiorb,dimens)
            call zgemm('n','n',dimens,dimens,dimens,cOne,prefactor,dimens,chiorb_hflsoc,dimens,cOne,chiorb,dimens)
          else
            ! chiorb = prefactor*chi_hf
            call zgemm('n','n',dimens,dimens,dimens,cOne,prefactor,dimens,chiorb_hf,dimens,cZero,chiorb,dimens) ! (1+chi_hf*Umat)^-1 * chi_hf
          end if

          schi = cZero
          schihf = cZero
          ! Calculating RPA and HF susceptibilities
          do j = 1, s%nAtoms
            do i = 1, s%nAtoms
              do sigmap = 1, 4
                do sigma = 1, 4
                  do nu = 1, s%Types(s%Basis(j)%Material)%nOrb
                    do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
                      schi  (sigmai2i(sigma,i), sigmai2i(sigmap,j)) = schi  (sigmai2i(sigma,i), sigmai2i(sigmap,j)) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                      schihf(sigmai2i(sigma,i), sigmai2i(sigmap,j)) = schihf(sigmai2i(sigma,i), sigmai2i(sigmap,j)) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                    end do
                  end do
                end do
              end do
            end do
          end do

          ! Rotating susceptibilities to the magnetization direction
          if(lrot) then
            do i = 1, s%nAtoms
              call rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,1)
              rotmat_i(:,:,i) = rottemp
              call rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,2)
              rotmat_j(:,:,i) = rottemp
            end do
            do j = 1, s%nAtoms
              do i = 1, s%nAtoms

                rottemp  = rotmat_i(:,:,i)
                do sigma = 1,4
                  do sigmap = 1,4
                    schitemp(sigma,sigmap) = schi(sigmai2i(sigma,i),sigmai2i(sigmap,j))
                  end do
                end do
                call zgemm('n','n',4,4,4,cOne,rottemp,4,schitemp,4,cZero,schirot,4)
                rottemp  = rotmat_j(:,:,j)
                call zgemm('n','n',4,4,4,cOne,schirot,4,rottemp,4,cZero,schitemp,4)
                do sigma = 1, 4
                  do sigmap = 1, 4
                    schi(sigmai2i(sigma,i), sigmai2i(sigmap,j)) = schitemp(sigma,sigmap)
                  end do
                end do
                do sigma = 1, 4
                  do sigmap = 1, 4
                    schitemp(sigma, sigmap) = schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j))
                  end do
                end do
                rottemp  = rotmat_i(:,:,i)
                call zgemm('n','n',4,4,4,cOne,rottemp,4,schitemp,4,cZero,schirot,4)
                rottemp  = rotmat_j(:,:,j)
                call zgemm('n','n',4,4,4,cOne,schirot,4,rottemp,4,cZero,schitemp,4)
                do sigma = 1, 4
                  do sigmap = 1, 4
                    schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schitemp(sigma, sigmap)
                  end do
                end do
              end do
            end do
          end if
        ! Calculating only HF susceptibilities for energy e
        else
          if(llinearsoc) then
            chiorb = chiorb_hf + chiorb_hflsoc
          else
            chiorb = chiorb_hf
          end if

          schihf = cZero
          ! Calculating RPA and HF susceptibilities
          do j = 1, s%nAtoms
            do nu = 1, s%Types(s%Basis(j)%Material)%nOrb
              do i = 1, s%nAtoms
                do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
                  do sigmap = 1, 4
                    do sigma = 1, 4
                      schihf(sigmai2i(sigma,i), sigmai2i(sigmap,j)) = schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                    end do
                  end do
                end do
              end do
            end do
          end do

          ! Rotating susceptibilities to the magnetization direction
          if(lrot) then
            do i = 1, s%nAtoms
              call rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,1)
              rotmat_i(:,:,i) = rottemp
              call rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,2)
              rotmat_j(:,:,i) = rottemp
            end do
            do j = 1, s%nAtoms
              do i = 1, s%nAtoms
                do sigma = 1, 4
                  do sigmap = 1, 4
                    schitemp(sigma, sigmap) = schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j))
                  end do
                end do
                rottemp  = rotmat_i(:,:,i)
                call zgemm('n','n',4,4,4,cOne,rottemp,4,schitemp,4,cZero,schirot,4)
                rottemp  = rotmat_j(:,:,j)
                call zgemm('n','n',4,4,4,cOne,schirot,4,rottemp,4,cZero,schitemp,4)
                do sigma = 1, 4
                  do sigmap = 1, 4
                    schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schitemp(sigma, sigmap)
                  end do
                end do
              end do
            end do
          end if
        end if ! .not.lhfresponses

        ! Calculating inverse susceptibility to use on Beff calculation
        chiinv = cZero
        do i = 1, s%nAtoms
          do j = 1, s%nAtoms
            do sigmap = 1, 4
              do sigma = 1, 4
                do nu = 1, s%Types(s%Basis(j)%Material)%nOrb
                  do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
                    chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))    ! +- , up- , down- , --
                  end do
                end do
              end do
            end do
          end do
        end do
        call invers(chiinv,dimspinAtoms) ! Inverse of the susceptibility chi^(-1)

        disturbances = cZero
        torques      = cZero
        !plane_loop_calculate_all:
        do i=1, s%nAtoms
          ! Spin and charge disturbances
          do mu=1,s%Types(s%Basis(i)%Material)%nOrb
            disturbances(1,i) = disturbances(1,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
            disturbances(2,i) = disturbances(2,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
            disturbances(3,i) = disturbances(3,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))/cI
            disturbances(4,i) = disturbances(4,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
          end do

          ! Spin disturbance matrix to calculate effective field
          sdmat(sigmai2i(1,i)) = disturbances(2,i) + cI*disturbances(3,i) ! +    = x + iy
          sdmat(sigmai2i(2,i)) = disturbances(1,i) + disturbances(4,i)    ! up   = 0 + z
          sdmat(sigmai2i(3,i)) = disturbances(1,i) - disturbances(4,i)    ! down = 0 - z
          sdmat(sigmai2i(4,i)) = disturbances(2,i) - cI*disturbances(3,i) ! -    = x - iy

          ! Orbital angular momentum disturbance in the global frame
          do nu = 1, s%Types(s%Basis(i)%Material)%nOrb
            do mu = 1, s%Types(s%Basis(i)%Material)%nOrb
              ldmat(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
              disturbances(5,i) = disturbances(5,i) + s%Types(s%Basis(i)%Material)%lvec(mu,nu,1)*ldmat(i,mu,nu)
              disturbances(6,i) = disturbances(6,i) + s%Types(s%Basis(i)%Material)%lvec(mu,nu,2)*ldmat(i,mu,nu)
              disturbances(7,i) = disturbances(7,i) + s%Types(s%Basis(i)%Material)%lvec(mu,nu,3)*ldmat(i,mu,nu)
            end do
          end do

          ! Spin-orbit torques (calculated in the spin frame of reference), for the P orbitals
          do nud = 1,s%Types(s%Basis(i)%Material)%npOrb
            nu = s%Types(s%Basis(i)%Material)%pOrbs(nud)
            do mud = 1,s%Types(s%Basis(i)%Material)%npOrb
              mu = s%Types(s%Basis(i)%Material)%pOrbs(mud)
              ! x component: (Ly*Sz - Lz*Sy)/2
              torques(1,1,i) = torques(1,1,i) + s%Types(s%Basis(i)%Material)%LambdaP*(   s%Basis(i)%lpvec(mu,nu,2)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3))) &
                                              + s%Types(s%Basis(i)%Material)%LambdaP*(cI*s%Basis(i)%lpvec(mu,nu,3)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
              ! y component: (Lz*Sx - Lx*Sz)/2
              torques(1,2,i) = torques(1,2,i) + s%Types(s%Basis(i)%Material)%LambdaP*(   s%Basis(i)%lpvec(mu,nu,3)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                              - s%Types(s%Basis(i)%Material)%LambdaP*(   s%Basis(i)%lpvec(mu,nu,1)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)))
              ! z component: (Lx*Sy - Ly*Sx)/2
              torques(1,3,i) = torques(1,3,i) - s%Types(s%Basis(i)%Material)%LambdaP*(cI*s%Basis(i)%lpvec(mu,nu,1)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                              - s%Types(s%Basis(i)%Material)%LambdaP*(   s%Basis(i)%lpvec(mu,nu,2)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
            end do
          end do

          ! Spin-orbit torques (calculated in the spin frame of reference), for the D orbitals
          do nud = 1,s%Types(s%Basis(i)%Material)%ndOrb
            nu = s%Types(s%Basis(i)%Material)%dOrbs(nud)
            do mud = 1,s%Types(s%Basis(i)%Material)%ndOrb
              mu = s%Types(s%Basis(i)%Material)%dOrbs(mud)
              ! x component: (Ly*Sz - Lz*Sy)/2
              torques(1,1,i) = torques(1,1,i) + s%Types(s%Basis(i)%Material)%LambdaD*(   s%Basis(i)%lpvec(mu,nu,2)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3))) &
                                              + s%Types(s%Basis(i)%Material)%LambdaD*(cI*s%Basis(i)%lpvec(mu,nu,3)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
              ! y component: (Lz*Sx - Lx*Sz)/2
              torques(1,2,i) = torques(1,2,i) + s%Types(s%Basis(i)%Material)%LambdaD*(   s%Basis(i)%lpvec(mu,nu,3)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                              - s%Types(s%Basis(i)%Material)%LambdaD*(   s%Basis(i)%lpvec(mu,nu,1)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)))
              ! z component: (Lx*Sy - Ly*Sx)/2
              torques(1,3,i) = torques(1,3,i) - s%Types(s%Basis(i)%Material)%LambdaD*(cI*s%Basis(i)%lpvec(mu,nu,1)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                              - s%Types(s%Basis(i)%Material)%LambdaD*(   s%Basis(i)%lpvec(mu,nu,2)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
            end do
          end do
          torques(1,:,i) = 0.5_dp*torques(1,:,i)

          ! Exchange-correlation torques (calculated in the spin frame of reference)
          torques(2,1,i) = s%Basis(i)%Um*(mvec_cartesian(3,i)*disturbances(3,i)-mvec_cartesian(2,i)*disturbances(4,i))
          torques(2,2,i) = s%Basis(i)%Um*(mvec_cartesian(1,i)*disturbances(4,i)-mvec_cartesian(3,i)*disturbances(2,i))
          torques(2,3,i) = s%Basis(i)%Um*(mvec_cartesian(2,i)*disturbances(2,i)-mvec_cartesian(1,i)*disturbances(3,i))

          ! External torques (calculated in the spin frame of reference)
          if(lfield) then
            torques(3,1,i) = 2._dp*(hhw(2,i)*disturbances(4,i)-hhw(3,i)*disturbances(3,i))
            torques(3,2,i) = 2._dp*(hhw(3,i)*disturbances(2,i)-hhw(1,i)*disturbances(4,i))
            torques(3,3,i) = 2._dp*(hhw(1,i)*disturbances(3,i)-hhw(2,i)*disturbances(2,i))
          end if

          ! ! Calculating spin and charge current for each neighbor TODO : Re-Include
          ! neighbor_loop_calculate_all: do neighbor=n0sc1,n0sc2
          !   ! Charge current
          !   currents(1,neighbor,i) =   ttchiorbiikl(neighbor,sigmai2i(2,i),2)+  ttchiorbiikl(neighbor,sigmai2i(2,i),3)+  ttchiorbiikl(neighbor,sigmai2i(3,i),2)+  ttchiorbiikl(neighbor,sigmai2i(3,i),3)
          !   ! Spin currents
          !   currents(2,neighbor,i) =   ttchiorbiikl(neighbor,sigmai2i(1,i),2)+  ttchiorbiikl(neighbor,sigmai2i(1,i),3)+  ttchiorbiikl(neighbor,sigmai2i(4,i),2)+  ttchiorbiikl(neighbor,sigmai2i(4,i),3)
          !   currents(3,neighbor,i) =   ttchiorbiikl(neighbor,sigmai2i(1,i),2)+  ttchiorbiikl(neighbor,sigmai2i(1,i),3)-  ttchiorbiikl(neighbor,sigmai2i(4,i),2)-  ttchiorbiikl(neighbor,sigmai2i(4,i),3)
          !   currents(4,neighbor,i) =   ttchiorbiikl(neighbor,sigmai2i(2,i),2)+  ttchiorbiikl(neighbor,sigmai2i(2,i),3)-  ttchiorbiikl(neighbor,sigmai2i(3,i),2)-  ttchiorbiikl(neighbor,sigmai2i(3,i),3)
          !   ! Orbital Angular Momentum currents
          !   currents(5,neighbor,i) = Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3)
          !   currents(6,neighbor,i) = Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3)
          !   currents(7,neighbor,i) = Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3)
          !
          ! end do neighbor_loop_calculate_all

        end do

        disturbances(2:4,:) = 0.5_dp*disturbances(2:4,:)
        sdmat    = 0.5_dp*sdmat
        total_disturbances  = sum(disturbances,dim=2)
        total_torques       = sum(torques,dim=3)

        ! currents = currents !TODO:Re-Include
        ! currents(2:4,:,:) =-0.5_dp*currents(2:4,:,:) !TODO:Re-Include
        ! currents(3,:,:)   = currents(3,:,:)/cI !TODO:Re-Include

        ! Total currents for each neighbor direction (Sum of currents over all planes)
        ! total_currents = sum(currents,dim=3)  !TODO:Re-Include

        ! Calculating spin Hall angles
        !call calculate_sha() !TODO: Re-Include

        ! DC spin current (second order)
        ! dc_currents = 0._dp  !TODO:Re-Include
        ! plane_loop_dc_current_all: do i=1,Npl
        !   do j=1,3 ; do mu = 1,3 ; do nu = 1,3
        !     dc_currents(j,i) = dc_currents(j,i) + levi_civita(j,mu,nu)*abs(disturbances(mu+1,i))*abs(disturbances(nu+1,i))*sin(atan2(aimag(disturbances(nu+1,i)),real(disturbances(nu+1,i))) - atan2(aimag(disturbances(mu+1,i)),real(disturbances(mu+1,i))))
        !   end do ; end do ; end do
        ! end do plane_loop_dc_current_all

        ! Effective field in cartesian coordinates
        call zgemm('n','n',dimspinAtoms,1,dimspinAtoms,cOne,chiinv,dimspinAtoms,sdmat,dimspinAtoms,cZero,Beff,dimspinAtoms) ! Beff = chi^(-1)*SD
        ! plane_loop_effective_field_all
        do i = 1, s%nAtoms
          Beff_cart(1,i) =           (Beff(sigmai2i(2,i)) + Beff(sigmai2i(3,i))) ! 0
          Beff_cart(2,i) = 0.5_dp*   (Beff(sigmai2i(1,i)) + Beff(sigmai2i(4,i))) ! x
          Beff_cart(3,i) =-0.5_dp*cI*(Beff(sigmai2i(1,i)) - Beff(sigmai2i(4,i))) ! y
          Beff_cart(4,i) = 0.5_dp*   (Beff(sigmai2i(2,i)) - Beff(sigmai2i(3,i))) ! z
        end do
        total_Beff  = sum(Beff_cart,dim=2)

        ! Sending results to myrank_row_hw = 0 (first process of each field calculation) and writing on files
        if(rFreq(2) == 0) then
          do mcount = 1, sFreq(2)
            if (mcount/=1) then ! Receive all points except the first (that was calculated at myrank_row_hw)
              call MPI_Recv(e                ,1                ,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE ,4000,FreqComm(2),stat,ierr)
              call MPI_Recv(q                ,3                ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),4001,FreqComm(2),stat,ierr)
              call MPI_Recv(mvec_spherical   ,3*s%nAtoms       ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),4100,FreqComm(2),stat,ierr)
              if(.not.lhfresponses) &
                call MPI_Recv(schi             , s%nAtoms*s%nAtoms*16   , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4200,FreqComm(2),stat,ierr)
              call MPI_Recv(schihf             , s%nAtoms*s%nAtoms*16   , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4300,FreqComm(2),stat,ierr)
              call MPI_Recv(Beff_cart          , 4*s%nAtoms             , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4400,FreqComm(2),stat,ierr)
              call MPI_Recv(total_Beff         , 4                      , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4450,FreqComm(2),stat,ierr)
              call MPI_Recv(disturbances       , 7*s%nAtoms             , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4500,FreqComm(2),stat,ierr)
              call MPI_Recv(total_disturbances , 7                      , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4550,FreqComm(2),stat,ierr)
              call MPI_Recv(torques            , ntypetorque*3*s%nAtoms , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4800,FreqComm(2),stat,ierr)
              call MPI_Recv(total_torques      , ntypetorque*3          , MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4850,FreqComm(2),stat,ierr)
              !call MPI_Recv(currents         ,7*n0sc*s%nAtoms  ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4600,FreqComm(2),stat,ierr) !TODO:Re-Include
              !call MPI_Recv(dc_currents      ,3*s%nAtoms       ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),5100,FreqComm(2),stat,ierr) !TODO:Re-Include
              !call MPI_Recv(total_currents   ,7*n0sc           ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),4700,FreqComm(2),stat,ierr) !TODO:Re-Include
              !call MPI_Recv(sha_re           ,4*s%nAtoms       ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),5200,FreqComm(2),stat,ierr) !TODO:Re-Include
              !call MPI_Recv(sha_re_total     ,4                ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),5300,FreqComm(2),stat,ierr) !TODO:Re-Include
              !call MPI_Recv(sha_complex      ,4*s%nAtoms       ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),5400,FreqComm(2),stat,ierr) !TODO:Re-Include
              !call MPI_Recv(sha_complex_total,4                ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),5500,FreqComm(2),stat,ierr) !TODO:Re-Include
            end if

            ! DIAGONALIZING SUSCEPTIBILITY
            if((.not.lhfresponses).and.(.not.lnodiag)) call diagonalize_susceptibilities()

            ! WRITING SUSCEPTIBILITIES
            call write_susceptibilities(qcount,e)

            ! Renormalizing disturbances and currents by the total charge current to neighbor renormnb
            !if(renorm) then !TODO: Re-Include
              ! Obtaining current for renormalization
              !Icabs  = abs(total_currents(1,renormnb)) !TODO: Re-Include

              !rdisturbances   = disturbances/Icabs !TODO: Re-Include
              !rtorques        = torques/Icabs ! TODO: Re-Include
              ! rcurrents       = currents/Icabs  !TODO:Re-Include
              ! rtotal_currents = total_currents/Icabs  !TODO:Re-Include
            !end if !TODO: Re-Include

            ! WRITING DISTURBANCES
            call write_disturbances(e)

            ! WRITING EFFECTIVE FIELDS
            call write_beff(e)

            ! WRITING TORQUES
            call write_torques(e)

            ! WRITING CURRENTS
            !call write_currents(e)  !TODO:Re-Include

            ! WRITING SHA
            !call write_sha(e) !TODO: Re-Include
          end do

          call write_time("[calculate_all] Time after step " // trim(itos(kount)) // ": ",output%unit_loop)

        else
          call MPI_Send(e                  , 1                      , MPI_DOUBLE_PRECISION , 0 , 4000 , FreqComm(2) , ierr)
          call MPI_Send(q                  , 3                      , MPI_DOUBLE_PRECISION , 0 , 4001 , FreqComm(2) , ierr)
          call MPI_Send(mvec_spherical     , 3*s%nAtoms             , MPI_DOUBLE_PRECISION , 0 , 4100 , FreqComm(2) , ierr)
          if(.not.lhfresponses) &
            call MPI_Send(schi             , s%nAtoms*s%nAtoms*16   , MPI_DOUBLE_COMPLEX , 0 , 4200 , FreqComm(2) , ierr)
          call MPI_Send(schihf             , s%nAtoms*s%nAtoms*16   , MPI_DOUBLE_COMPLEX , 0 , 4300 , FreqComm(2) , ierr)
          call MPI_Send(Beff_cart          , 4*s%nAtoms             , MPI_DOUBLE_COMPLEX , 0 , 4400 , FreqComm(2) , ierr)
          call MPI_Send(total_Beff         , 4                      , MPI_DOUBLE_COMPLEX , 0 , 4450 , FreqComm(2) , ierr)
          call MPI_Send(disturbances       , 7*s%nAtoms             , MPI_DOUBLE_COMPLEX , 0 , 4500 , FreqComm(2) , ierr)
          call MPI_Send(total_disturbances , 7                      , MPI_DOUBLE_COMPLEX , 0 , 4550 , FreqComm(2) , ierr)
          call MPI_Send(torques            , ntypetorque*3*s%nAtoms , MPI_DOUBLE_COMPLEX , 0 , 4800 , FreqComm(2) , ierr)
          call MPI_Send(total_torques      , ntypetorque*3          , MPI_DOUBLE_COMPLEX , 0 , 4850 , FreqComm(2) , ierr)
          !call MPI_Send(currents         ,7*n0sc*s%nAtoms  ,MPI_DOUBLE_COMPLEX  ,0,4600,FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(dc_currents      ,3*s%nAtoms       ,MPI_DOUBLE_PRECISION,0,5100,FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(total_currents   ,7*n0sc           ,MPI_DOUBLE_COMPLEX  ,0,4700,FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_re           ,4*s%nAtoms       ,MPI_DOUBLE_PRECISION,0,5200,FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_re_total     ,4                ,MPI_DOUBLE_PRECISION,0,5300,FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_complex      ,4*s%nAtoms       ,MPI_DOUBLE_COMPLEX  ,0,5400,FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_complex_total,4                ,MPI_DOUBLE_COMPLEX  ,0,5500,FreqComm(2),ierr) !TODO:Re-Include
        end if
      end if

      ! Emergency stop
      call MPI_Barrier(FieldComm, ierr)
      call check_stop("stop",0,e)
    end do ! Energy (frequency) loop

  end do ! Wave vector loop

  call freeLocalEKMesh()

  ! Sorting results on files
  if(rField==0) then
    call sort_all_files()
  end if
  call deallocate_prefactors()
  call deallocate_susceptibilities()
  call deallocate_disturbances()
  !call deallocate_currents()  !TODO:Re-Include
  call deallocate_beff()
  call deallocate_torques()
  !call deallocate_sha() !TODO: Re-Include
! call deallocate_idia()

end subroutine calculate_all
