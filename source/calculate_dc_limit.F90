! This is the main subroutine to calculate the fixed-frequency
! (in particular, the dc-limit) response functions
subroutine calculate_dc_limit()
  use mod_parameters
  use mod_constants
  use mod_tight_binding, only: lambda
  use mod_mpi_pars
  use mod_magnet
  use mod_progress
  use mod_prefactors
  use mod_susceptibilities
  use mod_disturbances
  use mod_currents
  use mod_beff
  use mod_torques
  use mod_sha
  implicit none
  character(len=50) :: time
  integer           :: i,j,iw,sigma,sigmap,mu,nu,neighbor,hw_count_temp,count_temp,mpitag2
  real(double)      :: e,Icabs,mvec_spherical_temp(Npl,3)

  dc_count = dc_count + 1

  call allocate_susceptibilities()
  call allocate_disturbances()
  call allocate_currents()
  call allocate_beff()
  call allocate_torques()
  call allocate_sha()
  if(myrank_row.eq.0) then
    if(emin.le.2.d-6) then
      write(outputunit_loop,"('CALCULATING DC-LIMIT QUANTITIES AS A FUNCTION OF EXTERNAL FIELD (',a,')')") trim(dcfield(dcfield_dependence))
    else
      write(outputunit_loop,"('CALCULATING QUANTITIES AT FIXED FREQUENCY ',es9.2,' AS A FUNCTION OF EXTERNAL FIELD (',a,')')") e,trim(dcfield(dcfield_dependence))
    end if
  end if

  ! Mounting U and identity matrix
  call build_identity_and_U_matrix()

  ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
  call allocate_prefactors()
  call OAM_curr_hopping_times_L()

  dclimit_energy_loop: do count_temp=1,npt1
    count = count_temp
    e = emin + deltae*(count-1)
    mpitag2 = (Npl-Npl_i)*npt1 + count
    mpitag  = (hw_count-1)*(Npl_f-Npl_i)*npt1 + mpitag2

    ! Creating files and writing headers
    if((.not.laddresults).and.(myrank.eq.0).and.(hw_count.eq.1)) then
      call openclose_dc_chi_files(0)
      call openclose_dc_disturbance_files(0)
      call openclose_dc_currents_files(0)
      call openclose_dc_beff_files(0)
      call openclose_dc_torque_files(0)
      call openclose_dc_sha_files(0)
    end if

    if(lhfresponses) then
      if(myrank_row_hw.eq.0) write(outputunit_loop,"('[calculate_dc_limit] No renormalization will be done. Setting prefactors to identity and calculating HF susceptibilities... ')")
      prefactor     = identt
      if(llinearsoc) prefactorlsoc = identt
      call eintshechi(e)
      if(myrank_row_hw.eq.0) call write_time(outputunit_loop,'[calculate_dc_limit] Time after susceptibility calculation: ')
    else
      if(myrank_row_hw.eq.0) write(outputunit_loop,"('[calculate_dc_limit] Calculating prefactor to use in currents and disturbances calculation. ')")
      if(llinearsoc) then
        call eintshechilinearsoc(e) ! Note: chiorb_hflsoc = lambda*dchi_hf/dlambda(lambda=0)
        ! Broadcast chiorb_hflsoc to all processors of the same row
        call MPI_Bcast(chiorb_hflsoc,dim*dim,MPI_DOUBLE_COMPLEX,0,MPI_Comm_Row,ierr)
      else
        call eintshechi(e)
      end if

      ! Broadcast chiorb_hf to all processors of the same row
      call MPI_Bcast(chiorb_hf,dim*dim,MPI_DOUBLE_COMPLEX,0,MPI_Comm_Row,ierr)

      ! prefactor = (1 + chi_hf*Umat)^-1
      prefactor     = identt
      call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zum,prefactor,dim) ! prefactor = 1+chi_hf*Umat
      call invers(prefactor,dim)
      if(llinearsoc) then
        prefactorlsoc = identt
        if (.not. allocated(chiorb)) allocate(chiorb(dim,dim))
        chiorb = chiorb_hf-chiorb_hflsoc ! the array chiorb will be used as a temporary array
        call zgemm('n','n',dim,dim,dim,zum,chiorb,dim,Umatorb,dim,zum,prefactorlsoc,dim) ! prefactorlsoc = 1+chiorb*Umat = 1+(chi_hf + dchi_hf/dlambda)*Umat
        call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,prefactorlsoc,dim,zero,chiorb,dim) ! chiorb = prefactor*prefactorlsoc
        call zgemm('n','n',dim,dim,dim,zum,chiorb,dim,prefactor,dim,zero,prefactorlsoc,dim) ! prefactorlsoc = chiorb*prefactor = prefactor*prefactorlsoc*prefactor
      end if
      if(myrank_row_hw.eq.0) call write_time(outputunit_loop,'[calculate_dc_limit] Time after prefactor calculation: ')
    end if

    ! Start parallelized processes to calculate disturbances and currents for energy e
    call eintshe(e)

    if(myrank_row_hw.eq.0) call write_time(outputunit_loop,'[calculate_dc_limit] Time after energy integral: ')

    if(myrank_row.eq.0) then

      if(.not.lhfresponses) then
        ! Calculating the full matrix of RPA and HF susceptibilities for energy e
        if(llinearsoc) then
          ! chiorb = prefactorlsoc*chi_hf + prefactor*chi_hflsoc
          call zgemm('n','n',dim,dim,dim,zum,prefactorlsoc,dim,chiorb_hf,dim,zero,chiorb,dim)
          call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hflsoc,dim,zum,chiorb,dim)
        else
          ! chiorb = prefactor*chi_hf
          call zgemm('n','n',dim,dim,dim,zum,prefactor,dim,chiorb_hf,dim,zero,chiorb,dim) ! (1+chi_hf*Umat)^-1 * chi_hf
        end if

        schi = zero
        schihf = zero
        ! Calculating RPA and HF susceptibilities
        calculate_susceptibility_dclimit: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
          schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
          schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
        end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_dclimit
        ! Rotating susceptibilities to the magnetization direction
        if(lrot) then
          do i=1,Npl
            call build_rotation_matrices_chi(mtheta(i),mphi(i),rottemp,1)
            rotmat_i(:,:,i) = rottemp
            call build_rotation_matrices_chi(mtheta(i),mphi(i),rottemp,2)
            rotmat_j(:,:,i) = rottemp
          end do
          rotate_susceptibility_dclimit: do j=1,Npl ; do i=1,Npl
            rottemp  = rotmat_i(:,:,i)
            schitemp = schi(:,:,i,j)
            call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
            rottemp  = rotmat_j(:,:,j)
            call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
            schi(:,:,i,j) = schitemp

            rottemp  = rotmat_i(:,:,i)
            schitemp = schihf(:,:,i,j)
            call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
            rottemp  = rotmat_j(:,:,j)
            call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
            schihf(:,:,i,j) = schitemp
          end do ; end do rotate_susceptibility_dclimit
        end if
      else
        ! Calculating the full matrix of RPA and HF susceptibilities for energy e
        if(llinearsoc) then
          chiorb = chiorb_hf + chiorb_hflsoc
        else
          chiorb = chiorb_hf
        end if

        schihf = zero
        ! Calculating RPA and HF susceptibilities
        calculate_hfsusceptibility_dclimit: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
          schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
        end do ; end do ; end do ; end do ; end do ; end do calculate_hfsusceptibility_dclimit
        ! Rotating susceptibilities to the magnetization direction
        if(lrot) then
          do i=1,Npl
            call build_rotation_matrices_chi(mtheta(i),mphi(i),rottemp,1)
            rotmat_i(:,:,i) = rottemp
            call build_rotation_matrices_chi(mtheta(i),mphi(i),rottemp,2)
            rotmat_j(:,:,i) = rottemp
          end do
          rotate_hfsusceptibility_dclimit: do j=1,Npl ; do i=1,Npl
            rottemp  = rotmat_i(:,:,i)
            schitemp = schihf(:,:,i,j)
            call zgemm('n','n',4,4,4,zum,rottemp,4,schitemp,4,zero,schirot,4)
            rottemp  = rotmat_j(:,:,j)
            call zgemm('n','n',4,4,4,zum,schirot,4,rottemp,4,zero,schitemp,4)
            schihf(:,:,i,j) = schitemp
          end do ; end do rotate_hfsusceptibility_dclimit
        end if

      end if ! .not.lhfresponses

      ! Calculating inverse susceptibility to use on Beff calculation
      chiinv = zero
      calculate_susceptibility_Beff_dclimit: do nu=1,9 ; do j=1,Npl ; do sigmap=1,4 ; do mu=1,9 ; do i=1,Npl ; do sigma=1,4
        chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = chiinv(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))    ! +- , up- , down- , --
      end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_Beff_dclimit
      call invers(chiinv,dimsigmaNpl) ! Inverse of the susceptibility chi^(-1)

      disturbances = zero
      torques      = zero
      plane_loop_calculate_dclimit: do i=1,Npl
        ! Spin and charge disturbances
        do mu=1,9
          disturbances(1,i) = disturbances(1,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
          disturbances(2,i) = disturbances(2,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))
          disturbances(3,i) = disturbances(3,i) + (tchiorbiikl(sigmaimunu2i(1,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,mu),3))/zi
          disturbances(4,i) = disturbances(4,i) + (tchiorbiikl(sigmaimunu2i(2,i,mu,mu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,mu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,mu),3))
        end do

        ! Spin disturbance matrix to calculate effective field
        sdmat(sigmai2i(1,i)) = disturbances(2,i) + zi*disturbances(3,i) ! +    = x + iy
        sdmat(sigmai2i(2,i)) = disturbances(1,i) + disturbances(4,i)    ! up   = 0 + z
        sdmat(sigmai2i(3,i)) = disturbances(1,i) - disturbances(4,i)    ! down = 0 - z
        sdmat(sigmai2i(4,i)) = disturbances(2,i) - zi*disturbances(3,i) ! -    = x - iy

        ! Orbital angular momentum disturbance
        do nu=1,9; do mu=1,9
          ldmat(i,mu,nu) = tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)
          disturbances(5,i) = disturbances(5,i) + lxp(mu,nu)*ldmat(i,mu,nu)
          disturbances(6,i) = disturbances(6,i) + lyp(mu,nu)*ldmat(i,mu,nu)
          disturbances(7,i) = disturbances(7,i) + lzp(mu,nu)*ldmat(i,mu,nu)
        end do; end do

        ! Spin-orbit torques (calculated in the spin frame of reference)
        do nu=1,9; do mu=1,9
          ! x component: (Ly*Sz - Lz*Sy)/2
          torques(1,1,i) = torques(1,1,i) + (   lyp(mu,nu)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3))) &
                                      + (zi*lzp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
          ! y component: (Lz*Sx - Lx*Sz)/2
          torques(1,2,i) = torques(1,2,i) + (   lzp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                      - (   lxp(mu,nu)*(tchiorbiikl(sigmaimunu2i(2,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(2,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(3,i,mu,nu),3)))
          ! z component: (Lx*Sy - Ly*Sx)/2
          torques(1,3,i) = torques(1,3,i) - (zi*lxp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)-tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3))) &
                                      - (   lyp(mu,nu)*(tchiorbiikl(sigmaimunu2i(1,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(1,i,mu,nu),3)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),2)+tchiorbiikl(sigmaimunu2i(4,i,mu,nu),3)))
        end do; end do
        torques(1,:,i) = 0.5d0*lambda(i+1)*torques(1,:,i)

        ! Exchange-correlation torques (calculated in the spin frame of reference)
        torques(2,1,i) = U(i+1)*(mz(i)*disturbances(3,i)-my(i)*disturbances(4,i))
        torques(2,2,i) = U(i+1)*(mx(i)*disturbances(4,i)-mz(i)*disturbances(2,i))
        torques(2,3,i) = U(i+1)*(my(i)*disturbances(2,i)-mx(i)*disturbances(3,i))

        ! External torques (calculated in the spin frame of reference)
        if(lfield) then
          torques(3,1,i) = 2.d0*(hhwy(i+1)*disturbances(4,i)-hhwz(i+1)*disturbances(3,i))
          torques(3,2,i) = 2.d0*(hhwz(i+1)*disturbances(2,i)-hhwx(i+1)*disturbances(4,i))
          torques(3,3,i) = 2.d0*(hhwx(i+1)*disturbances(3,i)-hhwy(i+1)*disturbances(2,i))
        end if

        ! Calculating spin and charge current for each neighbor
        neighbor_loop_calculate_dclimit: do neighbor=n0sc1,n0sc2
          ! Charge current
          currents(1,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)+ttchiorbiikl  (neighbor,sigmai2i(3,i),2)+ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
          ! Spin currents
          currents(2,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)+ttchiorbiikl  (neighbor,sigmai2i(4,i),2)+ttchiorbiikl  (neighbor,sigmai2i(4,i),3)
          currents(3,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(1,i),2)+ttchiorbiikl  (neighbor,sigmai2i(1,i),3)-ttchiorbiikl  (neighbor,sigmai2i(4,i),2)-ttchiorbiikl  (neighbor,sigmai2i(4,i),3)/zi
          currents(4,neighbor,i) = ttchiorbiikl  (neighbor,sigmai2i(2,i),2)+ttchiorbiikl  (neighbor,sigmai2i(2,i),3)-ttchiorbiikl  (neighbor,sigmai2i(3,i),2)-ttchiorbiikl  (neighbor,sigmai2i(3,i),3)
          ! Orbital Angular Momentum currents
          currents(5,neighbor,i) = Lxttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lxttchiorbiikl(neighbor,sigmai2i(3,i),3)
          currents(6,neighbor,i) = Lyttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lyttchiorbiikl(neighbor,sigmai2i(3,i),3)
          currents(7,neighbor,i) = Lzttchiorbiikl(neighbor,sigmai2i(2,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(2,i),3)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),2)+Lzttchiorbiikl(neighbor,sigmai2i(3,i),3)
        end do neighbor_loop_calculate_dclimit
      end do plane_loop_calculate_dclimit
      disturbances(2:4,:) = 0.5d0*disturbances(2:4,:)
      disturbances = disturbances/e
      torques  = torques/e
      sdmat    = 0.5d0*sdmat/e
      currents = currents/e
      currents(2:4,:,:) =-0.5d0*currents(2:4,:,:)
      currents(3,:,:)   = currents(3,:,:)/zi
      ! Total currents for each neighbor direction (Sum of currents over all planes)
      total_currents = sum(currents,dim=3)

      ! Calculating spin Hall angles
      call calculate_sha()

      ! DC spin current (second order)
      dc_currents = 0.d0
      plane_loop_dc_current_all: do i=1,Npl
        do j=1,3 ; do mu = 1,3 ; do nu = 1,3
          dc_currents(j,i) = dc_currents(j,i) + levi_civita(j,mu,nu)*abs(disturbances(mu+1,i))*abs(disturbances(nu+1,i))*sin(atan2(aimag(disturbances(nu+1,i)),real(disturbances(nu+1,i))) - atan2(aimag(disturbances(mu+1,i)),real(disturbances(mu+1,i))))/(tpi*e)
        end do ; end do ; end do
      end do plane_loop_dc_current_all

      ! Effective field calculation
      call zgemm('n','n',dimsigmaNpl,1,dimsigmaNpl,zum,chiinv,dimsigmaNpl,sdmat,dimsigmaNpl,zero,Beff,dimsigmaNpl) ! Beff = chi^(-1)*SD
      plane_loop_effective_field_dclimit: do i=1,Npl
        Beff_cart(sigmai2i(1,i)) =       (Beff(sigmai2i(2,i)) + Beff(sigmai2i(3,i)))    ! 0
        Beff_cart(sigmai2i(2,i)) = 0.5d0*(Beff(sigmai2i(1,i)) + Beff(sigmai2i(4,i)))    ! x
        Beff_cart(sigmai2i(3,i)) = 0.5d0*(Beff(sigmai2i(1,i)) - Beff(sigmai2i(4,i)))/zi ! y
        Beff_cart(sigmai2i(4,i)) = 0.5d0*(Beff(sigmai2i(2,i)) - Beff(sigmai2i(3,i)))    ! z
      end do plane_loop_effective_field_dclimit

      ! Sending results to myrank_row = myrank_col = 0 and writing on file
      if(myrank_col.eq.0) then
        ! Storing rank 0 counter and magnetization direction to recover later
        hw_count_temp       = hw_count
        mvec_spherical_temp = mvec_spherical

        MPI_points_dclimit: do mcount=1,MPIpts_hw
          if (mcount.ne.1) then ! Receive all points except the first (that was calculated at myrank_row_hw)
            call MPI_Recv(hw_count         ,1                ,MPI_INTEGER         ,MPI_ANY_SOURCE  ,44000+mpitag2,MPI_Comm_Col,stat,ierr) ! hw_count needed to get correct values of fields on hw_list
            call MPI_Recv(count            ,1                ,MPI_INTEGER         ,stat(MPI_SOURCE),44100+mpitag2,MPI_Comm_Col,stat,ierr) ! count needed to write correct energy on the filename
            call MPI_Recv(mvec_spherical   ,3*Npl            ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),44200+mpitag2,MPI_Comm_Col,stat,ierr)
            if(.not.lhfresponses) &
            call MPI_Recv(schi             ,Npl*Npl*16       ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44300+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(schihf           ,Npl*Npl*16       ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44400+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(Beff_cart        ,dimsigmaNpl      ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44500+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(disturbances     ,7*Npl            ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44600+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(currents         ,7*n0sc*Npl       ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44700+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(total_currents   ,7*n0sc           ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44800+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(torques          ,ntypetorque*3*Npl,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),44900+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(dc_currents      ,3*Npl            ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),45200+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(sha_re           ,4*Npl            ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),45300+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(sha_re_total     ,4                ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),45400+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(sha_complex      ,4*Npl            ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),45500+mpitag2,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(sha_complex_total,4                ,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),45600+mpitag2,MPI_Comm_Col,stat,ierr)
          end if

          ! DIAGONALIZING SUSCEPTIBILITY
          if((.not.lhfresponses).and.(.not.lnodiag)) call diagonalize_susceptibilities()

          ! WRITING SUSCEPTIBILITIES
          ! Opening chi and diag files
          call openclose_dc_chi_files(1)
          ! Writing susceptibilities
          call write_dc_susceptibilities()
          ! Closing chi and diag files
          call openclose_dc_chi_files(2)

          ! Renormalizing disturbances and currents by the total charge current to neighbor renormnb
          if(renorm) then
            ! Obtaining current for renormalization
            Icabs  = abs(total_currents(1,renormnb))

            rdisturbances   = disturbances/Icabs
            rtorques        = torques/Icabs
            rcurrents       = currents/Icabs
            rtotal_currents = total_currents/Icabs
          end if

          ! WRITING DISTURBANCES
          ! Opening disturbance files
          call openclose_dc_disturbance_files(1)
          ! Writing disturbances
          call write_dc_disturbances()
          ! Closing disturbance files
          call openclose_dc_disturbance_files(2)

          ! WRITING CURRENTS
          ! Opening current files
          call openclose_dc_currents_files(1)
          ! Writing currents
          call write_dc_currents()
          ! Closing current files
          call openclose_dc_currents_files(2)

          ! Opening B effective files
          call openclose_dc_beff_files(1)
          ! Writing effective fields
          call write_dc_beff()
          ! Closing B effective files
          call openclose_dc_beff_files(2)

          ! WRITING TORQUES
          ! Opening torque files
          call openclose_dc_torque_files(1)
          ! Writing torques
          call write_dc_torques()
          ! Closing torque files
          call openclose_dc_torque_files(2)

          ! WRITING SHA
          ! Opening SHA files
          call openclose_dc_sha_files(1)
          ! Writing SHA
          call write_dc_sha()
          ! Closing SHA files
          call openclose_dc_sha_files(2)
        end do MPI_points_dclimit

        write(time,"('[calculate_dc_limit] Time after step ',i0,': ')") dc_count
        call write_time(outputunit_loop,time)

        ! Emergency stop
        open(unit=911, file="stop", status='old', iostat=iw)
        if(iw.eq.0) then
          close(911)
            write(outputunit,"('[calculate_dc_limit] Emergency ""stop"" file found! Stopping after step ',i0,'...')") dc_count
            call system ('rm stop')
            write(outputunit,"('[calculate_dc_limit] (""stop"" file deleted!)')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
        ! Recovering rank 0 counter and magnetization direction
        hw_count       = hw_count_temp
        mvec_spherical = mvec_spherical_temp
      else
        if(myrank_row_hw.eq.0) write(outputunit_loop,"('[calculate_dc_limit] Sending results to process zero... ')")
        call MPI_Send(hw_count         ,1                ,MPI_INTEGER         ,0,44000+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(count_temp       ,1                ,MPI_INTEGER         ,0,44100+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(mvec_spherical   ,3*Npl            ,MPI_DOUBLE_PRECISION,0,44200+mpitag2,MPI_Comm_Col,ierr)
        if(.not.lhfresponses) &
        call MPI_Send(schi             ,Npl*Npl*16       ,MPI_DOUBLE_COMPLEX  ,0,44300+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(schihf           ,Npl*Npl*16       ,MPI_DOUBLE_COMPLEX  ,0,44400+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(Beff_cart        ,dimsigmaNpl      ,MPI_DOUBLE_COMPLEX  ,0,44500+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(disturbances     ,7*Npl            ,MPI_DOUBLE_COMPLEX  ,0,44600+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(currents         ,7*n0sc*Npl       ,MPI_DOUBLE_COMPLEX  ,0,44700+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(total_currents   ,7*n0sc           ,MPI_DOUBLE_COMPLEX  ,0,44800+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(torques          ,ntypetorque*3*Npl,MPI_DOUBLE_COMPLEX  ,0,44900+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(dc_currents      ,3*Npl            ,MPI_DOUBLE_PRECISION,0,45200+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(sha_re           ,4*Npl            ,MPI_DOUBLE_PRECISION,0,45300+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(sha_re_total     ,4                ,MPI_DOUBLE_PRECISION,0,45400+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(sha_complex      ,4*Npl            ,MPI_DOUBLE_COMPLEX  ,0,45500+mpitag2,MPI_Comm_Col,ierr)
        call MPI_Send(sha_complex_total,4                ,MPI_DOUBLE_COMPLEX  ,0,45600+mpitag2,MPI_Comm_Col,ierr)
        write(time,"('[calculate_dc_limit] Time after step ',i0,': ')") dc_count
        call write_time(outputunit_loop,time)
      end if
    end if
  end do dclimit_energy_loop

  ! Sorting results on files
  if(myrank_row_hw.eq.0) then
    call sort_all_files()
  end if
  call deallocate_prefactors()
  call deallocate_susceptibilities()
  call deallocate_disturbances()
  call deallocate_currents()
  call deallocate_beff()
  call deallocate_torques()
  call deallocate_sha()
  return
end subroutine calculate_dc_limit