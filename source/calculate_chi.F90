! This is the main subroutine to calculate the susceptibilities
subroutine calculate_chi()
  use mod_kind,              only: dp
  use mod_constants,         only: cZero,cOne,StoC,CtoS
  use mod_parameters,        only: kount,emin,deltae,nQvec1,kpoints,dimens,sigmaimunu2i,output,lhfresponses,lnodiag,laddresults,skip_steps,sigmai2i
  use mod_magnet,            only: mvec_spherical,lrot!,lfield
  use mod_alpha,             only: create_alpha_files,write_alpha
  use mod_system,            only: s => sys
  use mod_BrillouinZone,     only: realBZ
  use mod_tools,             only: itos,invers
  use adaptiveMesh,          only: genLocalEKMesh,freeLocalEKMesh
  use mod_progress,          only: write_time
  use TorqueTorqueResponse,  only: calcTTResponse,create_TTR_files,allocTTResponse
  use TorqueSpinResponse,    only: calcTSResponse,create_TSR_files,allocTSResponse
  use mod_rotation_matrices, only: rotation_matrices_chi
  use mod_SOC,               only: SOC
  use mod_Coupling,          only: get_J_K_from_chi
  use mod_susceptibilities,  only: identt,Umatorb,schi,schihf,schiLS,schiSL,schiLL,schirot,rotmat_i,&
                                   rotmat_j,rottemp,schitemp,chiorb_hf,chiorb,&
                                   build_identity_and_U_matrix,diagonalize_susceptibilities,&
                                   create_chi_files,write_susceptibilities,&
                                   allocate_susceptibilities,deallocate_susceptibilities
  use mod_sumrule,           only: sumrule
  use mod_mpi_pars,          only: rFreq,sFreq,FreqComm,FieldComm,rField,startFreq,endFreq,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,stat,ierr,MPI_SOURCE,MPI_DOUBLE_COMPLEX
  use mod_check_stop,        only: check_stop
  use adaptiveMesh,          only: bzs
  implicit none
  integer  :: mcount,qcount
  integer  :: i,j,sigma,sigmap,mu,nu,gama,xi,p
  real(dp) :: e,q(3)
  complex(dp), dimension(:,:), allocatable :: temp

  external :: eintshechi,zgemm,MPI_Recv,MPI_Send,MPI_Barrier,sort_all_files

  call allocate_susceptibilities()
  call allocTTResponse(s%ndAtoms)
  call allocTSResponse(s%ndAtoms)
  call genLocalEKMesh(s,rFreq(1), sFreq(1), FreqComm(1),bzs)
  call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))

  if(rFreq(1) == 0) allocate(temp(dimens,dimens))
  if(rField == 0) then
    write(output%unit_loop,"('CALCULATING LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
    ! write(outputunit_loop,"('Qx = ',es10.3,', Qz = ',es10.3)") q(1),q(3)
    ! Creating files and writing headers
    if(.not.laddresults) then
      call create_chi_files()
      call create_alpha_files()
      call create_TTR_files()
      call create_TSR_files()
    end if
  end if

  ! Mounting U and identity matrix
  call build_identity_and_U_matrix()

  if(rField == 0 .and. skip_steps > 0) &
    write(output%unit_loop,"('[calculate_chi] Skipping first ',i0,' step(s)...')") skip_steps

  ! Chi wave vector Loop
  do qcount=1,nQvec1
    if(rField==0) &
      write(output%unit_loop,"('[calculate_chi] Wave vector Q loop: ',i0,' of ',i0,' points',', Q = [',es10.3,es10.3,es10.3,']')") qcount,nQvec1,(kpoints(i,qcount),i=1,3)
    q = kpoints(:,qcount)
    ! Chi Energy (frequency) Loop
    do kount = startFreq+skip_steps, endFreq
      e = emin + deltae * (kount-1)

      if(rField==0) &
        write(output%unit_loop,"('[calculate_chi] Starting MPI step ',i0,' of ',i0,':',10x,'E =  ',es10.3)") kount - startFreq - skip_steps + 1, endFreq - startFreq + 1, e

      ! Start parallelized processes to calculate chiorb_hf and chiorbi0_hf for energy e
      call eintshechi(q,e)

      ! Time now
      if(rField == 0) &
        call write_time('[calculate_chi] Time after calculating chi HF: ',output%unit_loop)

      ! Checking sum rule for e=0._dp
      if((abs(e) < 1.e-8_dp).and.(sum(abs(q)) < 1.e-8_dp)) call sumrule(chiorb_hf)

      ! From here on all other processes except for rFreq(1) == 0 are idle :/
      if(rFreq(1)==0) then
        ! (1 + chi_hf*Umat)^-1
        temp = identt
        call zgemm('n','n',dimens,dimens,dimens,cOne,chiorb_hf,dimens,Umatorb,dimens,cOne,temp,dimens)
        ! temp = identt + temp
        call invers(temp,dimens)
        call zgemm('n','n',dimens,dimens,dimens,cOne,temp,dimens,chiorb_hf,dimens,cZero,chiorb,dimens)

        schi   = cZero
        schihf = cZero
        ! Calculating RPA and HF susceptibilities in global frame of reference
        do j=1, s%nAtoms
          do i=1, s%nAtoms
            do sigmap=1, 4
              do sigma=1, 4
                do nu=1, s%Types(s%Basis(j)%Material)%nOrb
                  do mu=1, s%Types(s%Basis(i)%Material)%nOrb
                    schi  (sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schi  (sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb   (sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                    schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                  end do
                end do
              end do
            end do
          end do
        end do

        if(SOC) then
          schiLS = cZero
          schiSL = cZero
          schiLL = cZero
          ! Calculating 3x3 <<L,S>>, <<S,L>> and <<L,L>> responses in cartesian components (x,y,z)
          do j=1, s%nAtoms
            do i=1, s%nAtoms
              do sigmap=1, 3
                do sigma=1, 3
                  do nu=1, s%Types(s%Basis(j)%Material)%nOrb
                    do mu=1, s%Types(s%Basis(i)%Material)%nOrb
                      do gama=1, s%Types(s%Basis(i)%Material)%nOrb
                        do p=1, 4
                          schiLS(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schiLS(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + s%Types(s%Basis(i)%Material)%lvec(mu,gama,sigma)*( chiorb(sigmaimunu2i(2,i,mu,gama),sigmaimunu2i(p,j,nu,nu)   ) + chiorb(sigmaimunu2i(3,i,mu,gama),sigmaimunu2i(p,j,nu,nu)   ) )*CtoS(p,sigmap+1)
                        end do ! p
                        do xi=1, s%Types(s%Basis(j)%Material)%nOrb
                          schiLL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schiLL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + s%Types(s%Basis(i)%Material)%lvec(mu,gama,sigma)*( chiorb(sigmaimunu2i(2,i,mu,gama),sigmaimunu2i(2,j,nu,xi)) + chiorb(sigmaimunu2i(2,i,mu,gama),sigmaimunu2i(3,j,nu,xi)) + chiorb(sigmaimunu2i(3,i,mu,gama),sigmaimunu2i(2,j,nu,xi)) + chiorb(sigmaimunu2i(3,i,mu,gama),sigmaimunu2i(3,j,nu,xi)) )*s%Types(s%Basis(i)%Material)%lvec(nu,xi,sigmap)
                        end do ! xi
                      end do ! gama

                      do xi=1, s%Types(s%Basis(j)%Material)%nOrb
                        do p=1, 4
                          schiSL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schiSL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) +      StoC(sigma+1,p)*( chiorb(sigmaimunu2i(p,i,mu,mu)   ,sigmaimunu2i(2,j,nu,xi)) + chiorb(sigmaimunu2i(p,i,mu,mu)   ,sigmaimunu2i(3,j,nu,xi)) )*s%Types(s%Basis(i)%Material)%lvec(nu,xi,sigmap)
                        end do ! p
                      end do ! gama
                    end do !mu
                  end do !nu
                end do !sigma
              end do !sigmap
            end do !i
          end do !j
        end if !SOC

        ! Rotating susceptibilities to the magnetization direction (local frame of reference)
        if(lrot) then
          do i=1, s%nAtoms
            ! Left rotation matrix
            call rotation_matrices_chi(180._dp-mvec_spherical(2,i),180._dp+mvec_spherical(3,i),rottemp,1)
            rotmat_i(:,:,i) = rottemp
            ! Right rotation matrix
            call rotation_matrices_chi(180._dp-mvec_spherical(2,i),180._dp+mvec_spherical(3,i),rottemp,2)
            rotmat_j(:,:,i) = rottemp
          end do
          do j=1, s%nAtoms
            do i=1, s%nAtoms
              rottemp  = rotmat_i(:,:,i)
              do sigma = 1, 4
                do sigmap = 1, 4
                  schitemp(sigma,sigmap) = schi(sigmai2i(sigma,i),sigmai2i(sigmap,j))
                end do
              end do
              call zgemm('n', 'n', 4, 4, 4, cOne, rottemp, 4, schitemp, 4, cZero, schirot, 4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n', 'n', 4, 4, 4, cOne, schirot, 4, rottemp, 4, cZero, schitemp, 4)
              do sigma = 1, 4
                do sigmap = 1, 4
                  schi(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schitemp(sigma,sigmap)
                end do
              end do

              rottemp  = rotmat_i(:,:,i)
              do sigma = 1, 4
                do sigmap = 1, 4
                  schitemp(sigma,sigmap) = schihf(sigmai2i(sigma,i), sigmai2i(sigmap,j))
                end do
              end do
              call zgemm('n', 'n', 4, 4, 4, cOne, rottemp, 4, schitemp, 4, cZero, schirot, 4)
              rottemp  = rotmat_j(:,:,j)
              call zgemm('n', 'n', 4, 4, 4, cOne, schirot, 4, rottemp, 4, cZero, schitemp, 4)
              do sigma = 1, 4
                do sigmap = 1, 4
                  schihf(sigmai2i(sigma,i), sigmai2i(sigmap,j)) = schitemp(sigma,sigmap)
                end do
              end do
            end do
          end do
        end if

        ! Sending results to myrank_row = myrank_col = 0 and writing on file
        if(rFreq(2) == 0) then
          do mcount=1, sFreq(2)
            if (mcount/=1) then
              call MPI_Recv(e,     1,                   MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE  ,1000,FreqComm(2),stat,ierr)
              call MPI_Recv(q,     3,                   MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1200,FreqComm(2),stat,ierr)
              call MPI_Recv(schi,  s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),1300,FreqComm(2),stat,ierr)
              call MPI_Recv(schihf,s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX  ,stat(MPI_SOURCE),1400,FreqComm(2),stat,ierr)
            end if

            ! DIAGONALIZING SUSCEPTIBILITY
            if((.not.lhfresponses).and.(.not.lnodiag)) call diagonalize_susceptibilities()

            ! TODO: Fix 
            ! if(.not.lfield) then
            !   call calcTTResponse(e)
            !   call calcTSResponse(e)
            ! end if

            ! WRITING GILBERT DAMPING
            if(abs(e)>1.e-8_dp) then
              call write_alpha(e)
            else if(sum(abs(q)) < 1.e-8_dp) then
              call get_J_K_from_chi()
            end if

            ! WRITING RPA AND HF SUSCEPTIBILITIES
            call write_susceptibilities(qcount,e)
          end do

          call write_time("[calculate_chi] Time after step " // trim(itos(kount)) // ": ",output%unit_loop)

        else
          call MPI_Send(e,     1,                   MPI_DOUBLE_PRECISION,0,1000, FreqComm(2),ierr)
          call MPI_Send(q,     3,                   MPI_DOUBLE_PRECISION,0,1100, FreqComm(2),ierr)
          call MPI_Send(schi,  s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX,  0,1200, FreqComm(2),ierr)
          call MPI_Send(schihf,s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX,  0,1300, FreqComm(2),ierr)
        end if
      end if

      ! Emergency stop
      call MPI_Barrier(FieldComm, ierr)
      call check_stop("stop",0,e)
    end do ! Energy (frequency) loop

  end do ! Wave vector loop

!   ! Sorting results on files
  if(rField == 0) call sort_all_files()

  call freeLocalEKMesh()

  call deallocate_susceptibilities()
  if(rFreq(1) == 0) deallocate(temp)

end subroutine calculate_chi
