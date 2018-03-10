! This is the main subroutine to calculate the susceptibilities
subroutine calculate_chi()
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cOne, StoC, CtoS
  use mod_parameters, only:  count, emin, deltae, dim, sigmaimunu2i, output, lnodiag,laddresults, skip_steps, sigmai2i
  use mod_mpi_pars
  use mod_magnet, only: lfield,mvec_spherical,lvec
  use mod_susceptibilities, only: identt, Umatorb, schi, schihf, schiLS, schiSL, schiLL, schirot, rotmat_i, &
       rotmat_j, rottemp, schitemp, lrot, chiorb_hf, chiorb, &
       build_identity_and_U_matrix, diagonalize_susceptibilities, &
       create_chi_files, write_susceptibilities, &
       allocate_susceptibilities, deallocate_susceptibilities
  use mod_alpha, only: create_alpha_files, write_alpha, allocate_alpha, deallocate_alpha
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use mod_BrillouinZone, only: realBZ
  use mod_tools, only: itos
  use adaptiveMesh, only: genLocalEKMesh, freeLocalEKMesh
  use mod_progress, only: write_time
  use TorqueTorqueResponse, only: calcTTResponse, create_TTR_files, allocTTResponse
  use TorqueSpinResponse, only: calcTSResponse, create_TSR_files, allocTSResponse
  use mod_rotation_matrices, only: rotation_matrices_chi
  use mod_sumrule
  implicit none
  character(len=50) :: time
  integer           :: mcount
  integer           :: i,j,sigma,sigmap,mu,nu,gamma,xi,p
  real(double)      :: e
  complex(double), dimension(:,:),   allocatable :: temp
  call allocate_susceptibilities()
  call allocate_alpha()
  call allocTTResponse(s%nAtoms)
  call allocTSResponse(s%nAtoms)
  call genLocalEKMesh(s,rFreq(1), sFreq(1), FreqComm(1))
  call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))


  if(rFreq(1) == 0) allocate(temp(dim,dim))
  if(rField == 0) then
     write(output%unit_loop,"('CALCULATING LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
     ! write(outputunit_loop,"('Qx = ',es10.3,', Qz = ',es10.3)") q(1),q(3)
     ! Creating files and writing headers
     if(.not.laddresults) then
        call create_chi_files()
        call create_alpha_files()
        call create_TTR_files
        call create_TSR_files()
     end if
  end if

  ! Mounting U and identity matrix
  call build_identity_and_U_matrix()

  if(rField == 0 .and. skip_steps > 0) write(output%unit_loop,"('[calculate_chi] Skipping first ',i0,' step(s)...')") skip_steps

  ! Chi Energy Loop
  do count = startFreq+skip_steps, endFreq+skip_steps
     e = emin + deltae * (count-1)
     if(rField==0) write(output%unit_loop,"('[calculate_chi] Starting MPI step ',i0,' of ',i0)") count - startFreq - skip_steps + 1, endFreq - startFreq + 1

     ! Start parallelized processes to calculate chiorb_hf and chiorbi0_hf for energy e
     call eintshechi(e)

     ! Time now
     if(rField == 0)  call write_time(output%unit_loop,'[calculate_chi] Time after calculating chi HF: ')

     ! Checking sum rule for e=0.d0
     if(e == 0.d0) call sumrule(chiorb_hf)

     ! From here on all other processes except for rFreq(1) == 0 are idle :/
     if(rFreq(1)==0) then
        ! (1 + chi_hf*Umat)^-1
        temp = identt
        call zgemm('n','n',dim,dim,dim,cOne,chiorb_hf,dim,Umatorb,dim,cOne,temp,dim)
        ! temp = identt + temp
        call invers(temp,dim)
        call zgemm('n','n',dim,dim,dim,cOne,temp,dim,chiorb_hf,dim,cZero,chiorb,dim)

        schi   = cZero
        schihf = cZero
        ! Calculating RPA and HF susceptibilities in global frame of reference
        do j=1, s%nAtoms
           do i=1, s%nAtoms
              do sigmap=1, 4
                 do sigma=1, 4
                    do nu=1, nOrb
                       do mu=1, nOrb
                          schi  (sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schi  (sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb   (sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                          schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schihf(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
                       end do
                    end do
                 end do
              end do
           end do
        end do

        schiLS = cZero
        schiSL = cZero
        schiLL = cZero
        ! Calculating 3x3 <<L,S>>, <<S,L>> and <<L,L>> responses in cartesian components (x,y,z)
        do j=1, s%nAtoms
           do i=1, s%nAtoms
              do sigmap=1, 3
                 do sigma=1, 3
                    do nu=1, nOrb
                       do mu=1, nOrb
                          do gamma=1, nOrb
                            do p=1, 4
                              schiLS(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schiLS(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + lvec(mu,gamma,sigma)*( chiorb(sigmaimunu2i(2,i,mu,gamma),sigmaimunu2i(p,j,nu,nu)) + chiorb(sigmaimunu2i(3,i,mu,gamma),sigmaimunu2i(p,j,nu,nu)) )*CtoS(p,sigmap+1)
                              schiSL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schiSL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + StoC(sigma+1,p)*( chiorb(sigmaimunu2i(p,i,mu,mu),sigmaimunu2i(2,j,nu,gamma)) + chiorb(sigmaimunu2i(p,i,mu,mu),sigmaimunu2i(3,j,nu,gamma)) )*lvec(nu,gamma,sigmap)
                            end do
                            do xi=1, nOrb
                              schiLL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) = schiLL(sigmai2i(sigma,i),sigmai2i(sigmap,j)) + lvec(mu,nu,sigma)*( chiorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) + chiorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) + chiorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) + chiorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) )*lvec(gamma,xi,sigmap)
                            end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do


        ! Rotating susceptibilities to the magnetization direction (local frame of reference)
        if(lrot) then
           do i=1, s%nAtoms
              call rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,1)
              rotmat_i(:,:,i) = rottemp
              call rotation_matrices_chi(mvec_spherical(2,i),mvec_spherical(3,i),rottemp,2)
              rotmat_j(:,:,i) = rottemp
           end do
           do j=1, s%nAtoms
              do i=1, s%nAtoms
                 rottemp  = rotmat_i(:,:,i)
                 do sigma = 1, 4
                   do sigmap = 1, 4
                     schitemp(sigma,sigmap) = schi(sigmai2i(sigmap,j),sigmai2i(sigma,i))
                   end do
                 end do
                 call zgemm('n', 'n', 4, 4, 4, cOne, rottemp, 4, schitemp, 4, cZero, schirot, 4)
                 rottemp  = rotmat_j(:,:,j)
                 call zgemm('n', 'n', 4, 4, 4, cOne, schirot, 4, rottemp, 4, cZero, schitemp, 4)
                 do sigma = 1, 4
                   do sigmap = 1, 4
                     schi(sigmai2i(sigmap,j),sigmai2i(sigma,i)) = schitemp(sigma,sigmap)
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
                 call MPI_Recv(e,     1,                   MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,1000,FreqComm(2),stat,ierr)
                 call MPI_Recv(schi,  s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1100,FreqComm(2),stat,ierr)
                 call MPI_Recv(schihf,s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1200,FreqComm(2),stat,ierr)
              end if

              ! DIAGONALIZING SUSCEPTIBILITY
              if((.not.lhfresponses).and.(.not.lnodiag)) call diagonalize_susceptibilities()

              if(.not.lfield) then
                ! Gonna be stupidly slow :/
                call calcTTResponse(e)
                call calcTSResponse(e)
              end if
              ! WRITING GILBERT DAMPING
              if(e/=0) call write_alpha(e)

              ! WRITING RPA AND HF SUSCEPTIBILITIES
              call write_susceptibilities(e)
           end do

           write(time,"('[calculate_chi] Time after step ',i0,': ')") count
           call write_time(output%unit_loop,time)

        else
           call MPI_Send(e,     1,                   MPI_DOUBLE_PRECISION,0,1000, FreqComm(2),ierr)
           call MPI_Send(schi,  s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX,  0,1100, FreqComm(2),ierr)
           call MPI_Send(schihf,s%nAtoms*s%nAtoms*16,MPI_DOUBLE_COMPLEX,  0,1200, FreqComm(2),ierr)
        end if
     end if
     call MPI_Barrier(FieldComm, ierr)
  end do

  ! Sorting results on files
  if(rField == 0) call sort_all_files()

  call freeLocalEKMesh()

  call deallocate_susceptibilities()
  call deallocate_alpha()

  if(rFreq(1) == 0) deallocate(temp)

  return
end subroutine calculate_chi
