! This is the main subroutine to calculate the susceptibilities
subroutine calculate_chi()
  use mod_parameters
  use mod_constants
  use mod_mpi_pars
  use mod_magnet, only: mtheta,mphi
  use mod_alpha, only: open_alpha_files, close_alpha_files, write_alpha
  use mod_progress
  use mod_susceptibilities
  implicit none
  character(len=50) :: time
  integer           :: i,j,iw,sigma,sigmap,mu,nu
  real(double)      :: e
  complex(double), dimension(:,:),   allocatable :: temp
  call allocate_susceptibilities()
  if(myrank_row==0) allocate(temp(dim,dim))
  if(myrank==0) then
    write(outputunit_loop,"('CALCULATING LOCAL SUSCEPTIBILITY AS A FUNCTION OF ENERGY')")
    ! write(outputunit_loop,"('Qx = ',es10.3,', Qz = ',es10.3)") q(1),q(3)
    ! Creating files and writing headers
    if(.not.laddresults) then
      call openclose_chi_files(0)
    end if
    call open_alpha_files()
  end if

  ! Mounting U and identity matrix
  call build_identity_and_U_matrix()

  if((myrank==0).and.(skip_steps>0)) write(outputunit_loop,"('[calculate_chi] Skipping first ',i0,' step(s)...')") skip_steps

  chi_energy_loop: do count=1+skip_steps,MPIsteps
    mpitag = (Npl-Npl_i)*total_hw_npt1*MPIsteps + (hw_count-1)*MPIsteps + count
    e = emin + deltae*myrank_col + MPIdelta*(count-1)
    if(myrank==0) then
      ! write(outputunit_loop,"(i0,' of ',i0,' points',', e = ',es10.3,' in myrank_col ',i0)") ((count-1)*MPIpts+myrank_col+1),npt1,e,myrank_col
      write(outputunit_loop,"('[calculate_chi] Starting MPI step ',i0,' of ',i0)") count,MPIsteps
    end if

    ! Start parallelized processes to calculate chiorb_hf and chiorbi0_hf for energy e
    call eintshechi(e)

    if(myrank_row==0) then
      ! (1 + chi_hf*Umat)^-1
      temp = identt
      call zgemm('n','n',dim,dim,dim,zum,chiorb_hf,dim,Umatorb,dim,zum,temp,dim)
      ! temp = identt + temp
      call invers(temp,dim)
      call zgemm('n','n',dim,dim,dim,zum,temp,dim,chiorb_hf,dim,zero,chiorb,dim)

      schi = zero
      schihf = zero
      ! Calculating RPA and HF susceptibilities
      calculate_susceptibility_chi: do j=1,Npl ; do nu=1,9 ; do i=1,Npl ; do mu=1,9 ; do sigmap=1,4 ; do sigma=1,4
        schi  (sigma,sigmap,i,j) = schi(sigma,sigmap,i,j)   + chiorb(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
        schihf(sigma,sigmap,i,j) = schihf(sigma,sigmap,i,j) + chiorb_hf(sigmaimunu2i(sigma,i,mu,mu),sigmaimunu2i(sigmap,j,nu,nu))
      end do ; end do ; end do ; end do ; end do ; end do calculate_susceptibility_chi
      ! Rotating susceptibilities to the magnetization direction
      if(lrot) then
        do i=1,Npl
          call build_rotation_matrices_chi(mtheta(i),mphi(i),rottemp,1)
          rotmat_i(:,:,i) = rottemp
          call build_rotation_matrices_chi(mtheta(i),mphi(i),rottemp,2)
          rotmat_j(:,:,i) = rottemp
        end do
        rotate_susceptibility_chi: do j=1,Npl ; do i=1,Npl
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
        end do ; end do rotate_susceptibility_chi
      end if

      ! Sending results to myrank_row = myrank_col = 0 and writing on file
      if(myrank_col==0) then
        MPI_points_chi: do mcount=1,MPIpts
          if (mcount/=1) then
            call MPI_Recv(e,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,1000,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(schi,Npl*Npl*16,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1100,MPI_Comm_Col,stat,ierr)
            call MPI_Recv(schihf,Npl*Npl*16,MPI_DOUBLE_COMPLEX,stat(MPI_SOURCE),1200,MPI_Comm_Col,stat,ierr)
          end if

          ! DIAGONALIZING SUSCEPTIBILITY
          if(.not.lnodiag) call diagonalize_susceptibilities()

          ! WRITING GILBERT DAMPING
          call write_alpha(e)
          ! WRITING RPA AND HF SUSCEPTIBILITIES
          ! Opening chi and diag files
          call openclose_chi_files(1)
          ! Writing susceptibilities
          call write_susceptibilities(e)
          ! Closing chi and diag files
          call openclose_chi_files(2)
        end do MPI_points_chi

        write(time,"('[calculate_chi] Time after step ',i0,': ')") count
        call write_time(outputunit_loop,time)

        ! Emergency stop
        open(unit=911, file="stop", status='old', iostat=iw)
        if(iw==0) then
          close(911)
          write(outputunit,"('[calculate_chi] Emergency ""stop"" file found! Stopping after step ',i0,'...')") count
          call system ('rm stop')
          write(outputunit,"('[calculate_chi] (""stop"" file deleted!)')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      else
        call MPI_Send(e,1,MPI_DOUBLE_PRECISION,0,1000,MPI_Comm_Col,ierr)
        call MPI_Send(schi,Npl*Npl*16,MPI_DOUBLE_COMPLEX,0,1100,MPI_Comm_Col,ierr)
        call MPI_Send(schihf,Npl*Npl*16,MPI_DOUBLE_COMPLEX,0,1200,MPI_Comm_Col,ierr)
      end if
    end if
  end do chi_energy_loop

  ! Sorting results on files
  if(myrank==0) then
    call close_alpha_files()
    call sort_all_files()
  end if

  call deallocate_susceptibilities()
  if(myrank_row==0) then
    deallocate(temp)
  end if

  return
end subroutine calculate_chi
