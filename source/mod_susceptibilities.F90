module mod_susceptibilities
  use mod_f90_kind, only: double
  implicit none
  ! Spin and orbital susceptibilities
  complex(double), dimension(:,:), allocatable :: schi, schihf
  complex(double), dimension(:,:), allocatable :: schiLS, schiSL, schiLL
  ! Full response functions
  complex(double), dimension(:,:), allocatable :: chiorb_hf,chiorb_hflsoc,chiorb
  complex(double), dimension(:,:), allocatable :: identt,Umatorb
  ! Rotation of spin susceptibilities
  logical :: lrot = .false.
  complex(double), dimension(:,:,:), allocatable :: rotmat_i, rotmat_j
  complex(double), dimension(:,:),   allocatable :: rottemp, schitemp, schirot

  ! Susceptibility diagonalization
  integer :: lwork
  real(double),    dimension(:),   allocatable :: rwork
  complex(double), dimension(:),   allocatable :: eval, work
  complex(double), dimension(:,:), allocatable :: chimag,evecl,evecr
#ifdef _JUQUEEN
  integer      :: ilo,ihi
  real(double) :: abnrm
  real(double), dimension(:), allocatable :: dscale, rconde, rcondv
#endif


contains

  subroutine allocate_susceptibilities()
    !! This subroutine allocates variables related to the susceptibility calculation
    use mod_f90_kind,   only: double
    use mod_parameters, only: dim
    use mod_SOC,        only: llinearsoc
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: rFreq,abortProgram
    implicit none
    integer :: AllocateStatus

    if(rFreq(1) == 0) then
       allocate( schi(4*s%nAtoms, 4*s%nAtoms), &
            schihf(4*s%nAtoms,4*s%nAtoms), &
            chiorb(dim,dim), STAT = AllocateStatus )
       if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: schi,schihf,chiorb")

       allocate( schiLS(3*s%nAtoms, 3*s%nAtoms), &
            schiSL(3*s%nAtoms,3*s%nAtoms), &
            schiLL(3*s%nAtoms,3*s%nAtoms), STAT = AllocateStatus )
       if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: schiLS, schiSL, schiLL")

       if(lrot) then
          allocate( rotmat_i(4,4,s%nAtoms), &
               rotmat_j(4,4,s%nAtoms), &
               rottemp(4,4), &
               schitemp(4,4), &
               schirot(4,4), stat = AllocateStatus )
          if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: rotmat_i,rotmat_j,rottemp,schitemp,schirot")
       end if

       if(rFreq(2) == 0) then
          if(s%nAtoms>1) then
             lwork = 33*s%nAtoms
             allocate( chimag(s%nAtoms,s%nAtoms), &
                  rwork(2*s%nAtoms), &
                  eval(s%nAtoms), &
                  evecl(1,s%nAtoms), &
                  evecr(s%nAtoms,s%nAtoms), &
                  work(lwork), STAT = AllocateStatus )
             if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: chimag, rwork, eval, evecl, evecr, work")
#ifdef _JUQUEEN
             allocate( dscale(s%nAtoms), &
                  rconde(s%nAtoms), &
                  rcondv(s%nAtoms), STAT = AllocateStatus )
             if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: dscale, rconde, rcondv")
#endif
          end if
       end if
    end if

    allocate( chiorb_hf(dim,dim), &
         Umatorb(dim,dim), &
         identt(dim,dim), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: chiorb_hf,Umatorb,identt")

    if(llinearsoc) then
       allocate( chiorb_hflsoc(dim,dim), STAT = AllocateStatus  )
       if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: chiorb_hflsoc")
    end if

  end subroutine allocate_susceptibilities

  subroutine deallocate_susceptibilities()
    !! This subroutine deallocates variables related to the susceptibility calculation
    implicit none

    if(allocated(schi)) deallocate(schi)
    if(allocated(schihf)) deallocate(schihf)
    if(allocated(schiLS)) deallocate(schiLS)
    if(allocated(schiSL)) deallocate(schiSL)
    if(allocated(schiLL)) deallocate(schiLL)
    if(allocated(rotmat_i)) deallocate(rotmat_i)
    if(allocated(rotmat_j)) deallocate(rotmat_j)
    if(allocated(rottemp)) deallocate(rottemp)
    if(allocated(schitemp)) deallocate(schitemp)
    if(allocated(schirot)) deallocate(schirot)
    if(allocated(chimag)) deallocate(chimag)
    if(allocated(rwork)) deallocate(rwork)
    if(allocated(eval)) deallocate(eval)
    if(allocated(evecl)) deallocate(evecl)
    if(allocated(evecr)) deallocate(evecr)
    if(allocated(work)) deallocate(work)

#ifdef _JUQUEEN
    if(allocated(dscale)) deallocate(dscale)
    if(allocated(rconde)) deallocate(rconde)
    if(allocated(rcondv)) deallocate(rcondv)
#endif

    if(allocated(chiorb)) deallocate(chiorb)
    if(allocated(chiorb_hf)) deallocate(chiorb_hf)
    if(allocated(chiorb_hflsoc)) deallocate(chiorb_hflsoc)
    if(allocated(Umatorb)) deallocate(Umatorb)
    if(allocated(identt)) deallocate(identt)

  end subroutine deallocate_susceptibilities

  subroutine build_identity_and_U_matrix()
    !! Mounts U and identity matrix
    use mod_parameters, only: offset, U, sigmaimunu2i, dim
    use mod_constants, only: cZero, cOne
    use mod_system, only: s => sys
    implicit none
    integer :: xi,gamma,nu,mu,i

    Umatorb = cZero
    do xi=5,9
       do gamma=5,9
          do nu=5,9
             do mu=5,9
                do i=1,s%nAtoms
                   if((mu/=nu).or.(gamma/=xi)) cycle
                   Umatorb(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,i,gamma,xi)) = cmplx(U(i+offset),0.d0)
                   Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,i,gamma,xi)) = cmplx(U(i+offset),0.d0)
                   Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,i,gamma,xi)) = cmplx(U(i+offset),0.d0)
                   Umatorb(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,i,gamma,xi)) = cmplx(U(i+offset),0.d0)
                end do
             end do
          end do
       end do
    end do
    do xi=5,9
       do gamma=5,9
          do nu=5,9
             do mu=5,9
                do i=1,s%nAtoms
                   if((mu/=xi).or.(nu/=gamma)) cycle
                   Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,i,gamma,xi)) = Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,i,gamma,xi))-cmplx(U(i+offset),0.d0)
                   Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,i,gamma,xi)) = Umatorb(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,i,gamma,xi))-cmplx(U(i+offset),0.d0)
                   Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,i,gamma,xi)) = Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,i,gamma,xi))-cmplx(U(i+offset),0.d0)
                   Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,i,gamma,xi)) = Umatorb(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,i,gamma,xi))-cmplx(U(i+offset),0.d0)
                end do
             end do
          end do
       end do
    end do

    identt     = cZero
    do i=1,dim
       identt(i,i) = cOne
    end do

  end subroutine build_identity_and_U_matrix


  ! This subroutine diagonalize the transverse susceptibility
  subroutine diagonalize_susceptibilities()
    use mod_f90_kind,   only: double
    use mod_system,     only: s => sys
    use mod_parameters, only: sigmai2i
    use mod_tools,      only: itos
    use mod_mpi_pars,   only: abortProgram
    implicit none
    integer  :: i,j,ifail

    if(s%nAtoms>1) then
       do i=1, s%nAtoms
          do j=1, s%nAtoms
             chimag(i,j) = schi(sigmai2i(1,i), sigmai2i(1,j))
          end do
       end do

#ifdef _JUQUEEN
       call zgeevx('N','N','V','N',s%nAtoms,chimag,s%nAtoms,eval,evecl,1,evecr,s%nAtoms,ilo,ihi,dscale,abnrm,rconde,rcondv,work,lwork,rwork,ifail)
#else
       call zgeev('N','V',s%nAtoms,chimag,s%nAtoms,eval,evecl,1,evecr,s%nAtoms,work,lwork,rwork,ifail)
#endif

       if(ifail/=0) then
          call abortProgram("[diagonalize_susceptibilities] Problem with diagonalization. ifail = " // trim(itos(ifail)))
          ! else
          !   write(outputunit,"('[diagonalize_susceptibilities] optimal lwork = ',i0,' lwork = ',i0)") work(1),lwork
       end if
    end if ! s%nAtoms

  end subroutine diagonalize_susceptibilities

  subroutine create_chi_files()
    !! This subroutine creates all the files needed for the susceptibilities
    use mod_parameters, only: output, lnodiag, lhfresponses
    use mod_system,     only: s => sys
    use mod_mpi_pars, only: abortProgram
    implicit none

    character(len=500)  :: varm
    integer :: i,j,iw

     do j=1,s%nAtoms
        do i=1,s%nAtoms
           iw = 1050 + (j-1) * s%nAtoms + i
           ! RPA SUSCEPTIBILITIES
           if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/chi_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),i,j,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
              open (unit=iw, file=varm, status='replace', form='formatted')
              write(unit=iw, fmt="('#     energy    ,  ((real[chi(j,i)], imag[chi(j,i)], j=1,4),i=1,4)')")
              close(unit=iw)
           end if
           iw = iw+1000
           ! HF SUSCEPTIBILITIES
           write(varm,"('./results/',a1,'SOC/',a,'/HF/chihf_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),i,j,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
           open (unit=iw, file=varm, status='replace', form='formatted')
           write(unit=iw, fmt="('#     energy    ,  ((real[chihf(j,i)], imag[chihf(j,i)], j=1,4),i=1,4)  ')")
           close(unit=iw)
        end do
     end do
    ! RPA DIAGONALIZATION
    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       write(varm,"('./results/',a1,'SOC/',a,'/RPA/chi_eval',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
       open (unit=19900, file=varm,status='replace', form='formatted')
       write(unit=19900,fmt="('#     energy    ,  real part of 1st eigenvalue  ,  imaginary part of 1st eigenvalue  ,  real part of 2nd eigenvalue  ,  imaginary part of 2nd eigenvalue  , ... ')")
       close (unit=19900)
       do i=1,s%nAtoms
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/chi_evec',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
          open (unit=19900+i, file=varm,status='replace', form='formatted')
          write(unit=19900+i,fmt="('#     energy    ,  real part of 1st component  ,  imaginary part of 1st component  ,  real part of 2nd component  ,  imaginary part of 2nd component  , ...   ')")
          close (unit=19900+i)
       end do
    end if

  end subroutine create_chi_files

  subroutine open_chi_files()
    !! This subroutine opens all the files needed for the susceptibilities
    use mod_parameters, only: output, lnodiag, lhfresponses, missing_files
    use mod_System,     only: s => sys
    use mod_mpi_pars, only: abortProgram
    implicit none

    character(len=500)  :: varm
    integer :: i,j,iw,err,errt=0

     do j=1,s%nAtoms
        do i=1,s%nAtoms
           iw = 1050 + (j-1) * s%nAtoms + i

           ! RPA SUSCEPTIBILITIES
           if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/chi_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),i,j,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
              open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
              errt = errt + err
              if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
           end if
           iw = iw+1000
           ! HF SUSCEPTIBILITIES
           write(varm,"('./results/',a1,'SOC/',a,'/HF/chihf_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),i,j,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
           open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
           errt = errt + err
           if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
        end do
     end do
    ! RPA DIAGONALIZATION
    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       write(varm,"('./results/',a1,'SOC/',a,'/RPA/chi_eval',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
       open (unit=19900, file=varm, status='old', position='append', form='formatted', iostat=err)
       errt = errt + err
       if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
       do i=1,s%nAtoms
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/chi_evec',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),i,trim(output%Energy),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
          open (unit=19900+i, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
       end do
    end if
    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openclose_chi_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

  end subroutine open_chi_files

  subroutine close_chi_files()
    !! This subroutine closes all the files needed for the susceptibilities
    use mod_parameters, only: lnodiag, lhfresponses
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none
    integer :: i,j,iw

     do j=1,s%nAtoms
        do i=1,s%nAtoms
           iw = 1050 + (j-1)*s%nAtoms + i
           ! RPA SUSCEPTIBILITIES
           if(.not.lhfresponses) close(unit=iw)
           iw = iw+1000
           ! HF SUSCEPTIBILITIES
           close(unit=iw)
        end do
     end do
    ! RPA DIAGONALIZATION
    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       close (unit=19900)
       do i=1,s%nAtoms
          close (unit=19900+i)
       end do
    end if

  end subroutine close_chi_files

  ! This subroutine write all the susceptibilities into files
  ! (already opened with openclose_chi_files(1))
  ! Some information may be written on the screen
  subroutine write_susceptibilities(e)
    use mod_f90_kind, only: double
    use mod_parameters, only: lwriteonscreen, lhfresponses, lnodiag, output, sigmai2i
    use mod_system, only: s => sys
    implicit none
    character(len=100)      :: varm
    integer                 :: i,j,iw,p,q
    real(double),intent(in) :: e

    call open_chi_files()

    if(lwriteonscreen) write(output%unit_loop,"(' #################  Susceptibilities:  #################')")
    do j=1, s%nAtoms
       do i=1, s%nAtoms
           iw = 1050 + (j-1)*s%nAtoms + i
           if(.not.lhfresponses) then
              write(unit=iw,fmt="(33(es16.9,2x))") e, ((real(schi(sigmai2i(q,i),sigmai2i(p,j))), aimag(schi(sigmai2i(q,i),sigmai2i(p,j))), q = 1, 4), p = 1, 4)
              if(i == j .and. lwriteonscreen) write(output%unit_loop,"('E = ',es11.4,', Plane: ',i0,' Chi+- = (',es16.9,') + i(',es16.9,')')") e,i,real(schi(sigmai2i(1,i),sigmai2i(1,j))),aimag(schi(sigmai2i(1,i),sigmai2i(1,j)))
           end if

           iw = iw+1000
           write(unit=iw,fmt="(33(es16.9,2x))") e, ((real(schihf(sigmai2i(q,i),sigmai2i(p,j))), aimag(schihf(sigmai2i(q,i),sigmai2i(p,j))), q = 1, 4), p = 1, 4)
       end do
    end do

    if(s%nAtoms > 1 .and. (.not. lhfresponses) .and. (.not. lnodiag)) then
       write(varm,fmt="('(',i0,'(es16.9,2x))')") 2*s%nAtoms+1
       write(unit=19900,fmt=varm) e,(real(eval(i)),aimag(eval(i)),i=1,s%nAtoms)
       do i=1,s%nAtoms
          write(unit=19900+i,fmt=varm) e,(real(evecr(j,i)),aimag(evecr(j,i)),j=1,s%nAtoms)
       end do
    end if

    call close_chi_files()

  end subroutine write_susceptibilities

  subroutine create_dc_chi_files()
    !! This subroutine creates all the files needed for the dc-limit susceptibilities
    use mod_parameters, only: count, output, lhfresponses, lnodiag
    use mod_magnet,     only: dcprefix, dcfield_dependence, dcfield, dc_header
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=500)  :: varm
    integer :: i,j,iw

     do j=1,s%nAtoms
        do i=1,s%nAtoms
           iw = 10500 + (j-1)*s%nAtoms + i
           ! RPA SUSCEPTIBILITIES
           if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a,'chi_',a,'_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
              open (unit=iw, file=varm, status='replace', form='formatted')
              write(unit=iw, fmt="('#',a,'  real part of chi ',a,'  ,  imaginary part of chi ',a,'  ,  amplitude of chi ',a,' ')") trim(dc_header)
              close(unit=iw)
           end if
           iw = iw+1000
           ! HF SUSCEPTIBILITIES
           write(varm,"('./results/',a1,'SOC/',a,'/HF/',a,'chihf_',a,'_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
           open (unit=iw, file=varm, status='replace', form='formatted')
           write(unit=iw, fmt="('#',a,'  real part of chi ',a,' HF ,  imaginary part of chi ',a,' HF  ,  amplitude of chi ',a,' HF ')") trim(dc_header)
           close(unit=iw)
        end do
     end do
    ! RPA DIAGONALIZATION
    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a,'chi_eval_',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
       open (unit=199000, file=varm,status='replace', form='formatted')
       write(unit=199000,fmt="('#',a,'  real part of 1st eigenvalue  ,  imaginary part of 1st eigenvalue  ,  real part of 2nd eigenvalue  ,  imaginary part of 2nd eigenvalue  , ...')") trim(dc_header)
       close (unit=199000)
       do i=1,s%nAtoms
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a,'chi_evec',i0,'_',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
          open (unit=199000+i, file=varm,status='replace', form='formatted')
          write(unit=199000+i,fmt="('#',a,'  real part of 1st component  ,  imaginary part of 1st component  ,  real part of 2nd component  ,  imaginary part of 2nd component  , ...')") trim(dc_header)
          close (unit=199000+i)
       end do
    end if

  end subroutine create_dc_chi_files

  subroutine open_dc_chi_files()
    !! This subroutine opens all the files needed for the dc-limit susceptibilities
    use mod_parameters, only: count, output, lhfresponses, lnodiag, missing_files
    use mod_magnet,     only: dcprefix, dcfield_dependence, dcfield
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=500)  :: varm
    integer :: i,j,iw,err,errt=0

     do j=1,s%nAtoms
        do i=1,s%nAtoms
           iw = 10500 + (j-1)*s%nAtoms + i
           ! RPA SUSCEPTIBILITIES
           if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a,'chi_',a,'_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
              open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
              errt = errt + err
              if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
           end if
           iw = iw+1000
           ! HF SUSCEPTIBILITIES
           write(varm,"('./results/',a1,'SOC/',a,'/HF/',a,'chihf_',a,'_',i0,'_',i0,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
           open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
           errt = errt + err
           if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
        end do
     end do
    ! RPA DIAGONALIZATION
    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a,'chi_eval_',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
       open (unit=199000, file=varm, status='old', position='append', form='formatted', iostat=err)
       errt = errt + err
       if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
       do i=1,s%nAtoms
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a,'chi_evec',i0,'_',a,a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(output%dcBField),trim(output%SOC),trim(output%suffix)
          open (unit=199000+i, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
       end do
    end if
    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openclose_dc_chi_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))


  end subroutine open_dc_chi_files

  subroutine close_dc_chi_files()
    !! This subroutine closes all the files needed for the dc-limit susceptibilities
    use mod_parameters, only: lhfresponses, lnodiag
    use mod_system,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none
    integer :: i,j,iw

     do j=1,s%nAtoms
        do i=1,s%nAtoms
           ! RPA SUSCEPTIBILITIES
           iw = 10500 + (j-1) * s%nAtoms+i
           if(.not.lhfresponses) close(unit=iw)
           ! HF SUSCEPTIBILITIES
           iw = iw+1000
           close(unit=iw)
        end do
     end do
    ! RPA DIAGONALIZATION
    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       close (unit=199000)
       do i=1,s%nAtoms
          close (unit=199000+i)
       end do
    end if

  end subroutine close_dc_chi_files


  ! This subroutine write all the susceptibilities into files
  ! (already opened with openclose_chi_files(1))
  ! Some information is also written on the screen
  subroutine write_dc_susceptibilities()
    use mod_f90_kind, only: double
    use mod_parameters, only: lwriteonscreen, output, lhfresponses, lnodiag, sigmai2i
    use mod_magnet, only: hw_count, dc_fields, dcfield_dependence, dcfield
    use mod_system, only: s => sys
    implicit none
    character(len=100)      :: varm
    integer                 :: i,j,iw,p,q

    call open_dc_chi_files()

    if(lwriteonscreen) write(output%unit_loop,"(' #################  Susceptibilities:  #################')")
    do j=1,s%nAtoms
      do i=1,s%nAtoms
       iw = 10500 + (j-1)*s%nAtoms + i
       if(.not.lhfresponses) then
          write(unit=iw,fmt="(a,2x,32(es16.9,2x))") trim(dc_fields(hw_count)) , ((real(schi(sigmai2i(q,i),sigmai2i(p,j))), aimag(schi(sigmai2i(q,i),sigmai2i(p,j))), q = 1, 4), p = 1, 4)
          if( i == j .and. lwriteonscreen ) write(output%unit_loop,"(a,' = ',a,', Plane: ',i0,' Chi+- = (',es16.9,') + i(',es16.9,')')") trim(dcfield(dcfield_dependence)),trim(dc_fields(hw_count)),i,real(schi(sigmai2i(1,i),sigmai2i(1,j))),aimag(schi(sigmai2i(1,i),sigmai2i(1,j)))
       end if

       iw = iw+1000
       write(unit=iw,fmt="(a,2x,32(es16.9,2x))") trim(dc_fields(hw_count)) , ((real(schi(sigmai2i(q,i),sigmai2i(p,j))), aimag(schi(sigmai2i(q,i),sigmai2i(p,j))), q = 1, 4), p = 1, 4)
      end do
    end do

    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       write(varm,fmt="(a,i0,a)") '(a,2x,',2*s%nAtoms,'(es16.9,2x))'
       write(unit=199000,fmt=varm) trim(dc_fields(hw_count)) , (real(eval(i)) , aimag(eval(i)),i=1,s%nAtoms)
       do i=1,s%nAtoms
          write(unit=199000+i,fmt=varm) trim(dc_fields(hw_count)) , (real(evecr(j,i)) , aimag(evecr(j,i)),j=1,s%nAtoms)
       end do
    end if ! s%nAtoms

    call close_dc_chi_files
  end subroutine write_dc_susceptibilities

  ! This subroutine sorts susceptibilities files
  subroutine sort_susceptibilities()
    use mod_f90_kind, only: double
    use mod_parameters, only: itype, lhfresponses, lnodiag
    use mod_tools, only: sort_file
    use mod_system, only: s => sys
    implicit none
    integer :: i,j,iw,idc=1

    ! Opening chi and diag files
    if(itype==9) then
       idc=10
       call open_dc_chi_files()
    else
       call open_chi_files()
    end if

    do j=1, s%nAtoms
       do i=1, s%nAtoms
         iw = 1050*idc + (j-1)*s%nAtoms + i
         if(.not.lhfresponses) then
            call sort_file(iw,.true.)
         end if

         iw = iw+1000
         call sort_file(iw,.true.)
       end do
    end do

    if((s%nAtoms>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
       call sort_file(19900*idc,.true.)
       do i=1,s%nAtoms
          call sort_file(19900*idc+i,.true.)
       end do
    end if

    ! Closing chi and diag files
    if(itype==9) then
       call close_dc_chi_files()
    else
       call close_chi_files()
    end if

  end subroutine sort_susceptibilities

end module mod_susceptibilities
