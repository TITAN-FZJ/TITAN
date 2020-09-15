! This is the main subroutine to calculate the fixed-frequency
! (in particular, the dc-limit) response functions
subroutine calculate_dc_limit()
  use mod_kind,              only: dp,int32,int64
  use mod_constants,         only: cZero, cOne, cI
  use mod_parameters,        only: sigmaimunu2i, sigmai2i, dimspinAtoms, Um, offset, lnodiag, output, kount, emin, deltae, nQvec1, kpoints, laddresults, lhfresponses, dimens, skip_steps
  use mod_magnet,            only: lfield, dcfield_dependence, dc_count, dcfield, hw_count, lxp, lyp, lzp, lx, ly, lz, mvec_cartesian, mvec_spherical, hhw, lrot
  use mod_SOC,               only: llinearsoc
  use mod_System,            only: s => sys
  use mod_BrillouinZone,     only: realBZ
  use mod_prefactors,        only: prefactor, prefactorlsoc, allocate_prefactors, deallocate_prefactors
  use mod_progress,          only: write_time
  use adaptiveMesh,          only: genLocalEKMesh, freeLocalEKMesh
  use mod_rotation_matrices, only: rotation_matrices_chi
  use mod_susceptibilities,  only: rottemp, rotmat_i, rotmat_j, &
                                   schitemp, schirot, schi, schihf, &
                                   chiorb, chiorb_hf, chiorb_hflsoc, Umatorb, identt, &
                                   build_identity_and_U_matrix, diagonalize_susceptibilities, &
                                   allocate_susceptibilities, deallocate_susceptibilities, &
                                   create_dc_chi_files, write_dc_susceptibilities
  use mod_torques,           only: torques, total_torques, ntypetorque, &
                                   allocate_torques, deallocate_torques, &
                                   create_dc_torque_files, write_dc_torques
  use mod_disturbances,      only: disturbances, total_disturbances, tchiorbiikl, sdmat, ldmat, &
                                   allocate_disturbances, deallocate_disturbances, &
                                   create_dc_disturbance_files, write_dc_disturbances
  use mod_beff,              only: chiinv, Beff, total_Beff, Beff_cart, &
                                   allocate_beff, deallocate_beff, &
                                   create_dc_beff_files, write_dc_beff
  use mod_check_stop,        only: check_stop
  use mod_tools,             only: invers
  use mod_mpi_pars,          only: rFreq,sFreq,FreqComm,FieldComm,rField,startFreq,endFreq,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,stat,ierr,MPI_SOURCE,MPI_DOUBLE_COMPLEX
  !use mod_system,          only: n0sc1, n0sc2, n0sc
  !use mod_currents !TODO: Re-Include
  !use mod_sha !TODO: Re-Include
  implicit none
  character(len=50) :: time
  integer(int32)    :: mcount,qcount
  integer(int64)    :: count_temp
  integer(int32)    :: i,j,sigma,sigmap,mu,nu,hw_count_temp!,neighbor 
  real(dp)      :: e,q(3),mvec_spherical_temp(3,s%nAtoms)!,Icabs

  external :: eintshechi,eintshelinearsoc,eintshechilinearsoc,eintshe,zgemm,MPI_Recv,MPI_Send,MPI_Barrier,sort_all_files,MPI_Bcast

  dc_count = dc_count + 1

  call genLocalEKMesh(s,rFreq(1),sFreq(1),FreqComm(1))
  call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))

  call allocate_prefactors()
  call allocate_susceptibilities()
  call allocate_disturbances()
  call allocate_beff()
  call allocate_torques()
  !call allocate_currents() !TODO: Re-Include
  !call allocate_sha() !TODO: Re-Include

  if(rFreq(1) == 0) then
    if(emin <= 2.e-6_dp) then
      write(output%unit_loop,"('CALCULATING DC-LIMIT QUANTITIES AS A FUNCTION OF EXTERNAL FIELD (',a,')')") trim(dcfield(dcfield_dependence))
    else
      write(output%unit_loop,"('CALCULATING QUANTITIES AT FIXED FREQUENCY ',es9.2,' AS A FUNCTION OF EXTERNAL FIELD (',a,')')") e,trim(dcfield(dcfield_dependence))
    end if
  end if

  ! Mounting U and identity matrix
  call build_identity_and_U_matrix()

  ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
  !call OAM_curr_hopping_times_L() !TODO: Re-Include

  ! Chi wave vector Loop
  do qcount=1,nQvec1
    if(rField==0) &
      write(output%unit_loop,"('[calculate_all] Wave vector Q loop: ',i0,' of ',i0,' points',', Q = [',es10.3,es10.3,es10.3,']')") qcount,nQvec1,(kpoints(i,qcount),i=1,3)
    q = kpoints(:,qcount)
    ! Chi Energy (frequency) Loop
    do count_temp = startFreq + skip_steps, endFreq + skip_steps
      kount = count_temp
      e = emin + deltae*(kount-1)

      ! Creating files and writing headers
      if((.not.laddresults).and. rField == 0 .and. hw_count == 1) then
        call create_dc_chi_files()
        call create_dc_disturbance_files()
        call create_dc_beff_files()
        call create_dc_torque_files()
        !call create_dc_currents_files() !TODO: Re-Include
        !call create_dc_sha_files() !TODO: Re-Include
      end if

      if(lhfresponses) then
        if(rField == 0) &
          write(output%unit_loop,"('[calculate_dc_limit] No renormalization will be done. Setting prefactors to identity and calculating HF susceptibilities... ')")
        prefactor     = identt
        if(llinearsoc) prefactorlsoc = identt
        call eintshechi(q,e)
        if(rField == 0) &
          call write_time(output%unit_loop,'[calculate_dc_limit] Time after susceptibility calculation: ')
      else
        if(rField == 0) &
          write(output%unit_loop,"('[calculate_dc_limit] Calculating prefactor to use in currents and disturbances calculation. ')")
        if(llinearsoc) then
          call eintshechilinearsoc(q,e) ! Note: chiorb_hflsoc = lambda*dchi_hf/dlambda(lambda=0)
          ! Broadcast chiorb_hflsoc to all processors of the same row
          call MPI_Bcast(chiorb_hflsoc,dimens*dimens,MPI_DOUBLE_COMPLEX,0,FreqComm(1),ierr)
        else
          call eintshechi(q,e)
        end if

        ! Broadcast chiorb_hf to all processors of the same row
        call MPI_Bcast(chiorb_hf,dimens*dimens,MPI_DOUBLE_COMPLEX,0,FreqComm(1),ierr)

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
          call write_time(output%unit_loop,'[calculate_dc_limit] Time after prefactor calculation: ')
      end if

      ! Start parallelized processes to calculate disturbances and currents for energy e
      if(llinearsoc) then
        call eintshelinearsoc(q,e)
      else
        call eintshe(q,e)
      end if

      if(rField == 0) &
        call write_time(output%unit_loop,'[calculate_dc_limit] Time after energy integral: ')

      if(rFreq(1) == 0) then
        if(.not.lhfresponses) then
          ! Calculating the full matrix of RPA and HF susceptibilities for energy e
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
            do nu = 1, 9
              do i = 1, s%nAtoms
                do mu = 1, 9
                  do sigmap = 1, 4
                    do sigma = 1, 4
                      schi  (sigmai2i(sigma,i), sigmai2i(sigmap,j)) = schi(sigmai2i(sigma,i), sigmai2i(sigmap,j))   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
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
        else
          ! Calculating the full matrix of RPA and HF susceptibilities for energy e
          if(llinearsoc) then
            chiorb = chiorb_hf + chiorb_hflsoc
          else
            chiorb = chiorb_hf
          end if

          schihf = cZero
          ! Calculating RPA and HF susceptibilities
          do j = 1, s%nAtoms
            do nu = 1, 9
              do i = 1, s%nAtoms
                do mu = 1, 9
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
        do nu = 1, 9
          do j = 1, s%nAtoms
            do sigmap = 1, 4
              do mu = 1, 9
                do i = 1, s%nAtoms
                  do sigma = 1, 4
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
        do i = 1, s%nAtoms
          ! Spin and charge disturbances
          do mu=1,9
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
          do nu = 1, 9
            do mu = 1, 9
              ldmat(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
              disturbances(5,i) = disturbances(5,i) + lx(mu,nu)*ldmat(i,mu,nu)
              disturbances(6,i) = disturbances(6,i) + ly(mu,nu)*ldmat(i,mu,nu)
              disturbances(7,i) = disturbances(7,i) + lz(mu,nu)*ldmat(i,mu,nu)
            end do
          end do

          ! Spin-orbit torques (calculated in the spin frame of reference), only for the D orbitals
          do nu = 5, 9
            do mu = 5, 9
              ! x component: (Ly*Sz - Lz*Sy)/2
              torques(1,1,i) = torques(1,1,i) + (   lyp(mu,nu,i)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3))) &
                                          + (cI*lzp(mu,nu,i)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
              ! y component: (Lz*Sx - Lx*Sz)/2
              torques(1,2,i) = torques(1,2,i) + (   lzp(mu,nu,i)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                          - (   lxp(mu,nu,i)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)))
              ! z component: (Lx*Sy - Ly*Sx)/2
              torques(1,3,i) = torques(1,3,i) - (cI*lxp(mu,nu,i)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                          - (   lyp(mu,nu,i)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
            end do
          end do
          torques(1,:,i) = 0.5_dp*s%Types(s%Basis(i)%Material)%LambdaD*torques(1,:,i)

          ! Exchange-correlation torques (calculated in the spin frame of reference)
          torques(2,1,i) = Um(i+offset)*(mvec_cartesian(3,i)*disturbances(3,i)-mvec_cartesian(2,i)*disturbances(4,i))
          torques(2,2,i) = Um(i+offset)*(mvec_cartesian(1,i)*disturbances(4,i)-mvec_cartesian(3,i)*disturbances(2,i))
          torques(2,3,i) = Um(i+offset)*(mvec_cartesian(2,i)*disturbances(2,i)-mvec_cartesian(1,i)*disturbances(3,i))

          ! External torques (calculated in the spin frame of reference)
          if(lfield) then
            torques(3,1,i) = 2._dp*(hhw(2,i+offset)*disturbances(4,i)-hhw(3,i+offset)*disturbances(3,i))
            torques(3,2,i) = 2._dp*(hhw(3,i+offset)*disturbances(2,i)-hhw(1,i+offset)*disturbances(4,i))
            torques(3,3,i) = 2._dp*(hhw(1,i+offset)*disturbances(3,i)-hhw(2,i+offset)*disturbances(2,i))
          end if

          ! Calculating spin and charge current for each neighbor
          ! do neighbor = n0sc1, n0sc2 !TODO:Re-Include
          !   ! Charge current
          !   currents(1,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)+ttchiorbiikl  (neighbor,sigmai2i(3,i),2)+ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
          !   ! Spin currents
          !   currents(2,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)+ttchiorbiikl  (neighbor,sigmai2i(4,i),2)+ttchiorbiikl  (neighbor,sigmai2i(4,i),3)
          !   currents(3,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)-ttchiorbiikl  (neighbor,sigmai2i(4,i),2)-ttchiorbiikl  (neighbor,sigmai2i(4,i),3)/cI
          !   currents(4,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)-ttchiorbiikl  (neighbor,sigmai2i(3,i),2)-ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
          !   ! Orbital Angular Momentum currents
          !   currents(5,neighbor,i) = Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3)
          !   currents(6,neighbor,i) = Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3)
          !   currents(7,neighbor,i) = Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3)
          ! end do
        end do

        disturbances(2:4,:) = 0.5_dp*disturbances(2:4,:)
        disturbances = disturbances/e
        torques  = torques/e
        total_disturbances  = sum(disturbances,dim=2)
        total_torques       = sum(torques,dim=3)
        sdmat    = 0.5_dp*sdmat/e
        ! currents = currents/e !TODO:Re-Include
        ! currents(2:4,:,:) =-0.5_dp*currents(2:4,:,:) !TODO:Re-Include
        ! currents(3,:,:)   = currents(3,:,:)/cI !TODO:Re-Include
        ! Total currents for each neighbor direction (Sum of currents over all planes)
        ! total_currents = sum(currents,dim=3) !TODO:Re-Include

        ! Calculating spin Hall angles
        !call calculate_sha() !TODO: Re-Include

        ! DC spin current (second order)
        ! dc_currents = 0._dp  !TODO:Re-Include
        ! do i = 1, s%nAtoms
        !   do j=1,3
        !     do mu = 1, 3
        !       do nu = 1, 3
        !         dc_currents(j,i) = dc_currents(j,i) + levi_civita(j,mu,nu)*abs(disturbances(mu+1,i))*abs(disturbances(nu+1,i))*sin(atan2(aimag(disturbances(nu+1,i)),real(disturbances(nu+1,i))) - atan2(aimag(disturbances(mu+1,i)),real(disturbances(mu+1,i))))/(tpi*e)
        !       end do
        !     end do
        !   end do
        ! end do
        !
        ! Effective field calculation
        call zgemm('n','n',dimspinAtoms,1,dimspinAtoms,cOne,chiinv,dimspinAtoms,sdmat,dimspinAtoms,cZero,Beff,dimspinAtoms) ! Beff = chi^(-1)*SD
        do i = 1, s%nAtoms
          Beff_cart(1,i) =       (Beff(sigmai2i(2,i)) + Beff(sigmai2i(3,i)))    ! 0
          Beff_cart(2,i) = 0.5_dp*(Beff(sigmai2i(1,i)) + Beff(sigmai2i(4,i)))    ! x
          Beff_cart(3,i) = 0.5_dp*(Beff(sigmai2i(1,i)) - Beff(sigmai2i(4,i)))/cI ! y
          Beff_cart(4,i) = 0.5_dp*(Beff(sigmai2i(2,i)) - Beff(sigmai2i(3,i)))    ! z
        end do
        total_Beff  = sum(Beff_cart,dim=2)

        ! Sending results to myrank_row = myrank_col = 0 and writing on file
        if(rFreq(2) == 0) then
          ! Storing rank 0 counter and magnetization direction to recover later
          hw_count_temp       = hw_count
          mvec_spherical_temp = mvec_spherical

          do mcount = 1, sFreq(2)
            if (mcount/=1) then ! Receive all points except the first (that was calculated at myrank_row_hw)
              call MPI_Recv(hw_count          ,1                     ,MPI_INTEGER         ,MPI_ANY_SOURCE ,44000,FreqComm(2),stat,ierr) ! hw_count needed to get correct values of fields on hw_list
              call MPI_Recv(kount             ,1                     ,MPI_INTEGER         ,stat(MPI_SOURCE),44100,FreqComm(2),stat,ierr) ! kount needed to write correct energy on the filename
              call MPI_Recv(q                 ,3                     ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),44101,FreqComm(2),stat,ierr)
              call MPI_Recv(mvec_spherical    ,3*s%nAtoms            ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),44200,FreqComm(2),stat,ierr)
              if(.not.lhfresponses) &
                call MPI_Recv(schi            ,s%nAtoms*s%nAtoms*16  ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44300,FreqComm(2),stat,ierr)
              call MPI_Recv(schihf            ,s%nAtoms*s%nAtoms*16  ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44400,FreqComm(2),stat,ierr)
              call MPI_Recv(Beff_cart         ,4*s%nAtoms            ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44500,FreqComm(2),stat,ierr)
              call MPI_Recv(total_Beff        ,4                     ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44550,FreqComm(2),stat,ierr)
              call MPI_Recv(disturbances      ,7*s%nAtoms            ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44600,FreqComm(2),stat,ierr)
              call MPI_Recv(total_disturbances,7                     ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44650,FreqComm(2),stat,ierr)
              call MPI_Recv(torques           ,ntypetorque*3*s%nAtoms,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44900,FreqComm(2),stat,ierr)
              call MPI_Recv(total_torques     ,ntypetorque*3         ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44950,FreqComm(2),stat,ierr)
              ! call MPI_Recv(currents         ,7*n0sc*s%nAtoms  ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44700,FreqComm(2),stat,ierr) !TODO:Re-Include
              ! call MPI_Recv(dc_currents      ,3*s%nAtoms       ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),45200,FreqComm(2),stat,ierr) !TODO:Re-Include
              ! call MPI_Recv(total_currents   ,7*n0sc           ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44800,FreqComm(2),stat,ierr) !TODO:Re-Include
              ! call MPI_Recv(sha_re           ,4*s%nAtoms       ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),45300,FreqComm(2),stat,ierr) !TODO:Re-Include
              ! call MPI_Recv(sha_re_total     ,4                ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),45400,FreqComm(2),stat,ierr) !TODO:Re-Include
              ! call MPI_Recv(sha_complex      ,4*s%nAtoms       ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),45500,FreqComm(2),stat,ierr) !TODO:Re-Include
              ! call MPI_Recv(sha_complex_total,4                ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),45600,FreqComm(2),stat,ierr) !TODO:Re-Include
            end if

            ! DIAGONALIZING SUSCEPTIBILITY
            if((.not.lhfresponses).and.(.not.lnodiag)) call diagonalize_susceptibilities()

            ! WRITING SUSCEPTIBILITIES
            call write_dc_susceptibilities(qcount)

            ! Renormalizing disturbances and currents by the total charge current to neighbor renormnb
            ! if(renorm) then  !TODO:Re-Include
            !   ! Obtaining current for renormalization
            !   Icabs  = abs(total_currents(1,renormnb))
            !
            !   rdisturbances   = disturbances/Icabs
            !   rtorques        = torques/Icabs
            !   rcurrents       = currents/Icabs
            !   rtotal_currents = total_currents/Icabs
            ! end if

            ! WRITING DISTURBANCES
            call write_dc_disturbances()

            ! WRITING CURRENTS
            !call write_dc_currents() !TODO:Re-Include

            ! Opening B effective files
            call write_dc_beff()

            ! WRITING TORQUES
            call write_dc_torques()

            ! WRITING SHA
            !call write_dc_sha() !TODO:Re-Include
          end do

          write(time,"('[calculate_dc_limit] Time after step ',i0,': ')") dc_count
          call write_time(output%unit_loop,time)

          ! Recovering rank 0 counter and magnetization direction
          hw_count       = hw_count_temp
          mvec_spherical = mvec_spherical_temp
          write(time,"('[calculate_dc_limit] Time after step ',i0,': ')") dc_count
          call write_time(output%unit_loop,time)
        else
          call MPI_Send(hw_count         ,1                ,MPI_INTEGER         ,0,44000, FreqComm(2),ierr)
          call MPI_Send(count_temp       ,1                ,MPI_INTEGER         ,0,44100, FreqComm(2),ierr)
          call MPI_Send(q                ,3                ,MPI_DOUBLE_PRECISION,0,44101, FreqComm(2),ierr)
          call MPI_Send(mvec_spherical   ,3*s%nAtoms       ,MPI_DOUBLE_PRECISION,0,44200, FreqComm(2),ierr)
          if(.not.lhfresponses) &
           call MPI_Send(schi             ,s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX  ,0,44300, FreqComm(2),ierr)
          call MPI_Send(schihf           ,s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX  ,0,44400, FreqComm(2),ierr)
          call MPI_Send(Beff_cart        , 4*s%nAtoms      ,MPI_DOUBLE_COMPLEX  ,0,44500, FreqComm(2),ierr)
          call MPI_Send(total_Beff        ,4                ,MPI_DOUBLE_COMPLEX  ,0,44550, FreqComm(2),ierr)
          call MPI_Send(disturbances     ,7*s%nAtoms       ,MPI_DOUBLE_COMPLEX  ,0,44600, FreqComm(2),ierr)
          call MPI_Send(total_disturbances,7                ,MPI_DOUBLE_COMPLEX  ,0,44650, FreqComm(2),ierr)
          call MPI_Send(torques          ,ntypetorque*3*s%nAtoms,MPI_DOUBLE_COMPLEX  ,0,44900, FreqComm(2),ierr)
          call MPI_Send(total_torques     ,ntypetorque*3    ,MPI_DOUBLE_COMPLEX  ,0,44950, FreqComm(2),ierr)
          !call MPI_Send(currents         ,7*n0sc*s%nAtoms  ,MPI_DOUBLE_COMPLEX  ,0,44700, FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(dc_currents      ,3*s%nAtoms       ,MPI_DOUBLE_PRECISION,0,45200, FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(total_currents   ,7*n0sc           ,MPI_DOUBLE_COMPLEX  ,0,44800, FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_re           ,4*s%nAtoms       ,MPI_DOUBLE_PRECISION,0,45300, FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_re_total     ,4                ,MPI_DOUBLE_PRECISION,0,45400, FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_complex      ,4*s%nAtoms       ,MPI_DOUBLE_COMPLEX  ,0,45500, FreqComm(2),ierr) !TODO:Re-Include
          !call MPI_Send(sha_complex_total,4                ,MPI_DOUBLE_COMPLEX  ,0,45600, FreqComm(2),ierr) !TODO:Re-Include
        end if
      end if

      ! Emergency stop
      call MPI_Barrier(FieldComm, ierr)
      call check_stop("stop",0,e)
    end do ! Energy (frequency) loop

  end do ! Wave vector loop

  ! Sorting results on files
  if(rField==0) then
    call sort_all_files()
  end if
  call deallocate_prefactors()
  call deallocate_susceptibilities()
  call deallocate_disturbances()
  call deallocate_beff()
  call deallocate_torques()
  !call deallocate_currents() !TODO:Re-Include
  !call deallocate_sha()  !TODO:Re-Include
end subroutine calculate_dc_limit
