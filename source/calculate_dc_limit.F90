! This is the main subroutine to calculate the fixed-frequency
! (in particular, the dc-limit) response functions
subroutine calculate_dc_limit()
  use mod_constants,       only: cZero, cOne, cI, levi_civita, tpi
  use mod_parameters, only: sigmaimunu2i, sigmai2i, dimsigmaNpl, U, offset, lnodiag, outputunit, count, emin, deltae, outputunit_loop, npt1, laddresults, lhfresponses, dim, skip_steps
  use mod_magnet,          only: lfield, dcfield_dependence, dc_count, dcfield, hw_count, lxp, lyp, lzp, lx, ly, lz, mvec_cartesian, mvec_spherical, hhwx, hhwy, hhwz
  use mod_SOC, only: llinearsoc
  use mod_System, only: s => sys
  use mod_prefactors, only: prefactor, prefactorlsoc, &
                            allocate_prefactors, deallocate_prefactors
  use mod_susceptibilities, only: lrot, rottemp, rotmat_i, rotmat_j, &
                                  schitemp, schirot, schi, schihf, &
                                  chiorb, chiorb_hf, chiorb_hflsoc, Umatorb, identt, &
                                  build_identity_and_U_matrix, diagonalize_susceptibilities, &
                                  allocate_susceptibilities, deallocate_susceptibilities, &
                                  create_dc_chi_files, write_dc_susceptibilities
  use mod_torques, only: torques, total_torques, ntypetorque, &
                         allocate_torques, deallocate_torques, &
                         create_dc_torque_files, write_dc_torques
  use mod_disturbances, only: disturbances, total_disturbances, tchiorbiikl, sdmat, ldmat, &
                              allocate_disturbances, deallocate_disturbances, &
                              create_dc_disturbance_files, write_dc_disturbances
  use mod_beff, only: chiinv, Beff, total_Beff, Beff_cart, &
                      allocate_beff, deallocate_beff, &
                      create_dc_beff_files, write_dc_beff
  use mod_progress,        only: write_time
  use adaptiveMesh
  use mod_mpi_pars

  !use mod_system,          only: n0sc1, n0sc2, n0sc
  !use mod_currents !TODO: Re-Include
  !use mod_sha !TODO: Re-Include
  implicit none
  character(len=50) :: time
  integer :: mcount
  integer           :: i,j,iw,sigma,sigmap,mu,nu,neighbor,hw_count_temp,count_temp,mpitag2
  real(double)      :: e,Icabs,mvec_spherical_temp(3,s%nAtoms)

  dc_count = dc_count + 1

  call genLocalEKMesh(rFreq(1),sFreq(1),FreqComm(1))

  call allocate_prefactors()
  call allocate_susceptibilities()
  call allocate_disturbances()
  call allocate_beff()
  call allocate_torques()
  !call allocate_currents() !TODO: Re-Include
  !call allocate_sha() !TODO: Re-Include

  if(rFreq(1) == 0) then
    if(emin <= 2.d-6) then
      write(outputunit_loop,"('CALCULATING DC-LIMIT QUANTITIES AS A FUNCTION OF EXTERNAL FIELD (',a,')')") trim(dcfield(dcfield_dependence))
    else
      write(outputunit_loop,"('CALCULATING QUANTITIES AT FIXED FREQUENCY ',es9.2,' AS A FUNCTION OF EXTERNAL FIELD (',a,')')") e,trim(dcfield(dcfield_dependence))
    end if
  end if

  ! Mounting U and identity matrix
  call build_identity_and_U_matrix()

  ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
  !call OAM_curr_hopping_times_L() !TODO: Re-Include

  do count_temp = startFreq + skip_steps, endFreq + skip_steps
    count = count_temp
    e = emin + deltae*(count-1)

    ! Creating files and writing headers
    if((.not.laddresults).and. myrank == 0 .and. hw_count == 1) then
      call create_dc_chi_files()
      call create_dc_disturbance_files()
      call create_dc_beff_files()
      call create_dc_torque_files()
      !call create_dc_currents_files() !TODO: Re-Include
      !call create_dc_sha_files() !TODO: Re-Include
    end if

    if(lhfresponses) then
      if(rField == 0) write(outputunit_loop,"('[calculate_dc_limit] No renormalization will be done. Setting prefactors to identity and calculating HF susceptibilities... ')")
      prefactor     = identt
      if(llinearsoc) prefactorlsoc = identt
      call eintshechi(e)
      if(rField == 0) call write_time(outputunit_loop,'[calculate_dc_limit] Time after susceptibility calculation: ')
    else
      if(rField == 0) write(outputunit_loop,"('[calculate_dc_limit] Calculating prefactor to use in currents and disturbances calculation. ')")
      if(llinearsoc) then
        call eintshechilinearsoc(e) ! Note: chiorb_hflsoc = lambda*dchi_hf/dlambda(lambda=0)
        ! Broadcast chiorb_hflsoc to all processors of the same row
        call MPI_Bcast(chiorb_hflsoc,dim*dim,MPI_DOUBLE_COMPLEX,0,FreqComm(1),ierr)
      else
        call eintshechi(e)
      end if

      ! Broadcast chiorb_hf to all processors of the same row
      call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,FreqComm(1),ierr)

      ! prefactor = (1 + chi_hf*Umat)^-1
      prefactor     = identt
      call zgemm('n','n',dim,dim,dim,cOne,chiorb_hf,dim,Umatorb,dim,cOne,prefactor,dim) ! prefactor = 1+chi_hf*Umat
      call invers(prefactor,dim)
      if(llinearsoc) then
        prefactorlsoc = identt
        if (.not. allocated(chiorb)) allocate(chiorb(dim,dim))
        chiorb = chiorb_hf-chiorb_hflsoc ! the array chiorb will be used as a temporary array
        call zgemm('n','n',dim,dim,dim,cOne,chiorb,dim,Umatorb,dim,cOne,prefactorlsoc,dim) ! prefactorlsoc = 1+chiorb*Umat = 1+(chi_hf + dchi_hf/dlambda)*Umat
        call zgemm('n','n',dim,dim,dim,cOne,prefactor,dim,prefactorlsoc,dim,cZero,chiorb,dim) ! chiorb = prefactor*prefactorlsoc
        call zgemm('n','n',dim,dim,dim,cOne,chiorb,dim,prefactor,dim,cZero,prefactorlsoc,dim) ! prefactorlsoc = chiorb*prefactor = prefactor*prefactorlsoc*prefactor
      end if
      if(rField == 0) call write_time(outputunit_loop,'[calculate_dc_limit] Time after prefactor calculation: ')
    end if

    ! Start parallelized processes to calculate disturbances and currents for energy e
    if(llinearsoc) then
      call eintshelinearsoc(e)
    else
      call eintshe(e)
    end if

    if(rField == 0) call write_time(outputunit_loop,'[calculate_dc_limit] Time after energy integral: ')

    if(rFreq(1) == 0) then

      if(.not.lhfresponses) then
        ! Calculating the full matrix of RPA and HF susceptibilities for energy e
        if(llinearsoc) then
          ! chiorb = prefactorlsoc*chi_hf + prefactor*chi_hflsoc
          call zgemm('n','n',dim,dim,dim,cOne,prefactorlsoc,dim,chiorb_hf,dim,cZero,chiorb,dim)
          call zgemm('n','n',dim,dim,dim,cOne,prefactor,dim,chiorb_hflsoc,dim,cOne,chiorb,dim)
        else
          ! chiorb = prefactor*chi_hf
          call zgemm('n','n',dim,dim,dim,cOne,prefactor,dim,chiorb_hf,dim,cZero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf
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
                    schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                    schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                  end do
                end do
              end do
            end do
          end do
        end do
        ! Rotating susceptibilities to the magnetization direction
        if(lrot) then
          do i = 1, s%nAtoms
            call build_rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,1)
            rotmat_i(:,:,i) = rottemp
            call build_rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,2)
            rotmat_j(:,:,i) = rottemp
          end do
          do j = 1, s%nAtoms
            do i = 1, s%nAtoms
              rottemp  = rotmat_i(:,:,i)
              schitemp = schi(:,:,i,j)
              call zgemm('n','n',4,4,4,cOne,rottemp,4,schitemp,4,cZero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,cOne,schirot,4,rottemp,4,cZero,schitemp,4)
              schi(:,:,i,j) = schitemp

              rottemp  = rotmat_i(:,:,i)
              schitemp = schihf(:,:,i,j)
              call zgemm('n','n',4,4,4,cOne,rottemp,4,schitemp,4,cZero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,cOne,schirot,4,rottemp,4,cZero,schitemp,4)
              schihf(:,:,i,j) = schitemp
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
                    schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                  end do
                end do
              end do
            end do
          end do
        end do
        ! Rotating susceptibilities to the magnetization direction
        if(lrot) then
          do i = 1, s%nAtoms
            call build_rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,1)
            rotmat_i(:,:,i) = rottemp
            call build_rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,2)
            rotmat_j(:,:,i) = rottemp
          end do
          do j = 1, s%nAtoms
            do i = 1, s%nAtoms
              rottemp  = rotmat_i(:,:,i)
              schitemp = schihf(:,:,i,j)
              call zgemm('n','n',4,4,4,cOne,rottemp,4,schitemp,4,cZero,schirot,4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n','n',4,4,4,cOne,schirot,4,rottemp,4,cZero,schitemp,4)
              schihf(:,:,i,j) = schitemp
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
      call invers(chiinv,dimsigmaNpl) ! Inverse of the susceptibility chi^(-1)

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

        ! Spin-orbit torques (calculated in the spin frame of reference)
        do nu = 1, 9
          do mu = 1, 9
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
        torques(1,:,i) = 0.5d0*s%Types(s%Basis(i)%Material)%Lambda*torques(1,:,i)

        ! Exchange-correlation torques (calculated in the spin frame of reference)
        torques(2,1,i) = U(i+offset)*(mvec_cartesian(3,i)*disturbances(3,i)-mvec_cartesian(2,i)*disturbances(4,i))
        torques(2,2,i) = U(i+offset)*(mvec_cartesian(1,i)*disturbances(4,i)-mvec_cartesian(3,i)*disturbances(2,i))
        torques(2,3,i) = U(i+offset)*(mvec_cartesian(2,i)*disturbances(2,i)-mvec_cartesian(1,i)*disturbances(3,i))

        ! External torques (calculated in the spin frame of reference)
        if(lfield) then
          torques(3,1,i) = 2.d0*(hhwy(i+offset)*disturbances(4,i)-hhwz(i+offset)*disturbances(3,i))
          torques(3,2,i) = 2.d0*(hhwz(i+offset)*disturbances(2,i)-hhwx(i+offset)*disturbances(4,i))
          torques(3,3,i) = 2.d0*(hhwx(i+offset)*disturbances(3,i)-hhwy(i+offset)*disturbances(2,i))
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

      disturbances(2:4,:) = 0.5d0*disturbances(2:4,:)
      disturbances = disturbances/e
      torques  = torques/e
      total_disturbances  = sum(disturbances,dim=2)
      total_torques       = sum(torques,dim=3)
      sdmat    = 0.5d0*sdmat/e
      ! currents = currents/e !TODO:Re-Include
      ! currents(2:4,:,:) =-0.5d0*currents(2:4,:,:) !TODO:Re-Include
      ! currents(3,:,:)   = currents(3,:,:)/cI !TODO:Re-Include
      ! Total currents for each neighbor direction (Sum of currents over all planes)
      ! total_currents = sum(currents,dim=3) !TODO:Re-Include

      ! Calculating spin Hall angles
      !call calculate_sha() !TODO: Re-Include

      ! DC spin current (second order)
      ! dc_currents = 0.d0  !TODO:Re-Include
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
      call zgemm('n','n',dimsigmaNpl,1,dimsigmaNpl,cOne,chiinv,dimsigmaNpl,sdmat,dimsigmaNpl,cZero,Beff,dimsigmaNpl) ! Beff = chi^(-1)*SD
      do i = 1, s%nAtoms
        Beff_cart(1,i) =       (Beff(sigmai2i(2,i)) + Beff(sigmai2i(3,i)))    ! 0
        Beff_cart(2,i) = 0.5d0*(Beff(sigmai2i(1,i)) + Beff(sigmai2i(4,i)))    ! x
        Beff_cart(3,i) = 0.5d0*(Beff(sigmai2i(1,i)) - Beff(sigmai2i(4,i)))/cI ! y
        Beff_cart(4,i) = 0.5d0*(Beff(sigmai2i(2,i)) - Beff(sigmai2i(3,i)))    ! z
      end do
      total_Beff  = sum(Beff_cart,dim=2)

      ! Sending results to myrank_row = myrank_col = 0 and writing on file
      if(rFreq(2) == 0) then
        ! Storing rank 0 counter and magnetization direction to recover later
        hw_count_temp       = hw_count
        mvec_spherical_temp = mvec_spherical

        do mcount = 1, sFreq(2)
          if (mcount/=1) then ! Receive all points except the first (that was calculated at myrank_row_hw)
            call MPI_Recv(hw_count         ,1                , MPI_INTEGER         ,MPI_ANY_SOURCE  ,44000,FreqComm(2),stat,ierr) ! hw_count needed to get correct values of fields on hw_list
            call MPI_Recv(count            ,1                , MPI_INTEGER         ,stat(MPI_SOURCE),44100,FreqComm(2),stat,ierr) ! count needed to write correct energy on the filename
            call MPI_Recv(mvec_spherical   ,3*s%nAtoms       , MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),44200,FreqComm(2),stat,ierr)
            if(.not.lhfresponses) &
               call MPI_Recv(schi             ,s%nAtoms*s%nAtoms*16, MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44300,FreqComm(2),stat,ierr)
            call MPI_Recv(schihf           ,s%nAtoms*s%nAtoms*16, MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44400,FreqComm(2),stat,ierr)
            call MPI_Recv(Beff_cart        ,4*s%nAtoms      ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44500,FreqComm(2),stat,ierr)
            call MPI_Recv(total_Beff        ,4                ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44550,FreqComm(2),stat,ierr)
            call MPI_Recv(disturbances     ,7*s%nAtoms       ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44600,FreqComm(2),stat,ierr)
            call MPI_Recv(total_disturbances,7                ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44650,FreqComm(2),stat,ierr)
            call MPI_Recv(torques          ,ntypetorque*3*s%nAtoms,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44900,FreqComm(2),stat,ierr)
            call MPI_Recv(total_torques     ,ntypetorque*3    ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44950,FreqComm(2),stat,ierr)
            ! call MPI_Recv(currents         ,7*n0sc*s%nAtoms  ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44700,FreqComm(2),stat,ierr)  !TODO:Re-Include
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
          call write_dc_susceptibilities()

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
        call write_time(outputunit_loop,time)

        ! Emergency stop
        open(unit=911, file="stop", status='old', iostat=iw)
        if(iw==0) then
          close(911)
            write(outputunit,"('[calculate_dc_limit] Emergency ""stop"" file found! Stopping after step ',i0,'...')") dc_count
            call system ('rm stop')
            write(outputunit,"('[calculate_dc_limit] (""stop"" file deleted!)')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        ! Recovering rank 0 counter and magnetization direction
        hw_count       = hw_count_temp
        mvec_spherical = mvec_spherical_temp
        write(time,"('[calculate_dc_limit] Time after step ',i0,': ')") dc_count
        call write_time(outputunit_loop,time)
     else
        call MPI_Send(hw_count         ,1                ,MPI_INTEGER         ,0,44000, FreqComm(2),ierr)
        call MPI_Send(count_temp       ,1                ,MPI_INTEGER         ,0,44100, FreqComm(2),ierr)
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
  end do

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
  return
end subroutine calculate_dc_limit
