module mod_susceptibilities
  use mod_f90_kind, only: double
  implicit none
  ! Spin susceptibilities
  complex(double), dimension(:,:,:,:), allocatable   :: schi, schihf, shirhoz
  ! Rotation of spin susceptibilities
  logical                       :: lrot = .false.
  complex(double), dimension(:,:,:), allocatable   :: rotmat_i, rotmat_j
  complex(double), dimension(:,:), allocatable     :: rottemp, schitemp, schirot
  ! Susceptibility diagonalization
  integer                       :: lwork
  real(double),allocatable      :: rwork(:)
  complex(double),allocatable   :: eval(:),work(:)
  complex(double), dimension(:,:),allocatable       :: chimag,evecl,evecr
#ifdef _JUQUEEN
  integer :: ilo,ihi
  real(double) :: abnrm
  real(double),allocatable      :: dscale(:),rconde(:),rcondv(:)
#endif
  ! Full response functions
  complex(double), dimension(:,:),   allocatable :: chiorb_hf,chiorb_hflsoc,chiorb
  complex(double), dimension(:,:),   allocatable :: identt,Umatorb
contains

  ! This subroutine allocates variables related to the susceptibility calculation
  subroutine allocate_susceptibilities()
    use mod_f90_kind, only: double
    use mod_parameters, only: dim, nmaglayers
    use mod_SOC, only: llinearsoc
    use mod_system, only: s => sys
    use mod_mpi_pars, only: abortProgram, myrank_row, myrank_col
    implicit none
    integer           :: AllocateStatus

    if(myrank_row==0) then
      allocate( schi(4,4,s%nAtoms, s%nAtoms),schihf(4,4,s%nAtoms, s%nAtoms),chiorb(dim,dim), STAT = AllocateStatus )
      if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: schi,schihf,chiorb")
      if(lrot) allocate( rotmat_i(4,4,s%nAtoms),rotmat_j(4,4,s%nAtoms),rottemp(4,4),schitemp(4,4),schirot(4,4) )

      if(myrank_col==0) then
        if(nmaglayers>1) then
          lwork = 33*nmaglayers
          allocate( chimag(nmaglayers,nmaglayers),rwork(2*nmaglayers),eval(nmaglayers),evecl(1,nmaglayers),evecr(nmaglayers,nmaglayers),work(lwork) )
#ifdef _JUQUEEN
          allocate( dscale(nmaglayers),rconde(nmaglayers),rcondv(nmaglayers) )
#endif
        end if
      end if
    end if
    allocate( chiorb_hf(dim,dim),Umatorb(dim,dim),identt(dim,dim), STAT = AllocateStatus )
    if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: chiorb_hf,Umatorb,identt")
    if(llinearsoc) then
      allocate( chiorb_hflsoc(dim,dim), STAT = AllocateStatus  )
      if (AllocateStatus/=0) call abortProgram("[allocate_susceptibilities] Not enough memory for: chiorb_hflsoc")
    end if

    return
  end subroutine allocate_susceptibilities

  ! Mounts U and identity matrix
  subroutine build_identity_and_U_matrix()
    use mod_parameters, only: offset, U, Utype, layertype, sigmaimunu2i, dim
    use mod_constants, only: zero, zum
    use mod_system, only: s => sys
    implicit none
    integer :: xi,gamma,nu,mu,i

    Umatorb = zero
    if(Utype/=0) then
      do xi=5,9
        do gamma=5,9
          do nu=5,9
            do mu=5,9
              do i=1,s%nAtoms
                if((Utype==1).and.(layertype(i+offset)/=2)) cycle
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
                if((Utype==1).and.(layertype(i+offset)/=2)) cycle
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
    end if

    identt     = zero
    do i=1,dim
     identt(i,i) = zum
    end do

    return
  end subroutine build_identity_and_U_matrix

  ! This subroutine deallocates variables related to the susceptibility calculation
  subroutine deallocate_susceptibilities()
    use mod_f90_kind, only: double
    use mod_parameters, only: nmaglayers
    use mod_SOC, only: llinearsoc
    use mod_mpi_pars, only: myrank_row, myrank_col
    implicit none

    if(myrank_row==0) then
      deallocate(schi,schihf)
      if(lrot) deallocate(rotmat_i,rotmat_j,rottemp,schitemp,schirot)
      if(myrank_col==0) then
        if(nmaglayers>1) then
          deallocate(chimag,rwork,eval,evecl,evecr,work)
#ifdef _JUQUEEN
          deallocate(dscale,rconde,rcondv)
#endif
        end if
      end if
    end if
    if(allocated(chiorb)) deallocate(chiorb)
    deallocate(chiorb_hf,Umatorb,identt)
    if(llinearsoc) deallocate(chiorb_hflsoc)

    return
  end subroutine deallocate_susceptibilities

  ! This subroutine diagonalize the transverse susceptibility
  subroutine diagonalize_susceptibilities()
    use mod_f90_kind, only: double
    use mod_parameters, only: nmaglayers,mmlayermag,outputunit
    use mod_mpi_pars, only: MPI_Abort, MPI_COMM_WORLD, errorcode, ierr
    implicit none
    integer  :: i,j,ifail

    if(nmaglayers>1) then
      do i=1, nmaglayers
        do j=1, nmaglayers
          chimag(i,j) = schi(1,1,mmlayermag(i)-1,mmlayermag(j)-1)
        end do
      end do

#ifdef _JUQUEEN
      call zgeevx('N','N','V','N',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,ilo,ihi,dscale,abnrm,rconde,rcondv,work,lwork,rwork,ifail)
#else
      call zgeev('N','V',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,work,lwork,rwork,ifail)
#endif

      if(ifail/=0) then
        write(outputunit,"('[diagonalize_susceptibilities] Problem with diagonalization. ifail = ',i0)") ifail
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      ! else
      !   write(outputunit,"('[diagonalize_susceptibilities] optimal lwork = ',i0,' lwork = ',i0)") work(1),lwork
      end if
    end if ! nmaglayers

    return
  end subroutine diagonalize_susceptibilities

  ! This subroutine opens and closes all the files needed for the susceptibilities
  subroutine openclose_chi_files(iflag)
    use mod_parameters, only: fieldpart, nmaglayers, Npl_folder, suffix, eta, lnodiag, ltestcharge, lhfresponses, Utype
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use mod_mpi_pars, only: abortProgram
    use EnergyIntegration, only: strEnergyParts
    implicit none

    character(len=500)  :: varm
    character(len=2)    :: spin(4)
    integer :: i,j,sigma,iw,iflag,err,errt=0

    spin(1) = "pm"
    spin(2) = "um"
    spin(3) = "dm"
    spin(4) = "mm"

    if(iflag==0) then
      do sigma=1,4
        do j=1,s%nAtoms
          do i=1,s%nAtoms
            iw = 1000 + (sigma-1) * s%nAtoms * s%nAtoms + (j-1) * s%nAtoms + i
            ! RPA SUSCEPTIBILITIES
            if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a2,'/chi_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
              open (unit=iw, file=varm, status='replace', form='formatted')
              write(unit=iw, fmt="('#     energy    ,  real part of chi ',a,'  ,  imaginary part of chi ',a,'  ,  amplitude of chi ',a,'  ')") spin(sigma),spin(sigma),spin(sigma)
              close(unit=iw)
            end if
            iw = iw+1000
            ! HF SUSCEPTIBILITIES
            write(varm,"('./results/',a1,'SOC/',a,'/HF/',a2,'/chihf_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
            open (unit=iw, file=varm, status='replace', form='formatted')
            write(unit=iw, fmt="('#     energy    ,  real part of chi ',a,' HF ,  imaginary part of chi ',a,' HF  ,  amplitude of chi ',a,' HF  ')") spin(sigma),spin(sigma),spin(sigma)
            close(unit=iw)
          end do
        end do
      end do
      ! RPA DIAGONALIZATION
      if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
        write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/chi_eval',a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
        open (unit=1990, file=varm,status='replace', form='formatted')
        write(unit=1990,fmt="('#     energy    ,  real part of 1st eigenvalue  ,  imaginary part of 1st eigenvalue  ,  real part of 2nd eigenvalue  ,  imaginary part of 2nd eigenvalue  , ... ')")
        close (unit=1990)
        do i=1,nmaglayers
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/chi_evec',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=1990+i, file=varm,status='replace', form='formatted')
          write(unit=1990+i,fmt="('#     energy    ,  real part of 1st component  ,  imaginary part of 1st component  ,  real part of 2nd component  ,  imaginary part of 2nd component  , ...   ')")
          close (unit=1990+i)
        end do
      end if
      if(ltestcharge) then
        write(varm, "('./results/',a1,'SOC/',a,'/RPA/testcharge',a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
        open (unit=17964, file=varm, status='replace', form='formatted', iostat=err)
        close(unit=17964)
      end if

    else if(iflag==1) then
      do sigma=1,4
        do j=1,s%nAtoms
          do i=1,s%nAtoms
            iw = 1000 + (sigma-1) * s%nAtoms * s%nAtoms + (j-1) * s%nAtoms + i

            ! RPA SUSCEPTIBILITIES
            if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a2,'/chi_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
              open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
              errt = errt + err
            end if
            iw = iw+1000
            ! HF SUSCEPTIBILITIES
            write(varm,"('./results/',a1,'SOC/',a,'/HF/',a2,'/chihf_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
            open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
            errt = errt + err
          end do
        end do
      end do
      ! RPA DIAGONALIZATION
      if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
        write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/chi_eval',a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
        open (unit=1990, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        do i=1,nmaglayers
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/chi_evec',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),i,trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=1990+i, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end do
      end if
      if(ltestcharge) then
        write(varm, "('./results/',a1,'SOC/',a,'/RPA/testcharge',a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(strEnergyParts),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
        open (unit=17964, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end if
      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_chi_files] Some file(s) do(es) not exist! Stopping before starting calculations...")

    else
      do sigma=1,4
        do j=1,s%nAtoms
          do i=1,s%nAtoms
            iw = 1000+(sigma-1)*s%nAtoms*s%nAtoms+(j-1)*s%nAtoms+i
            ! RPA SUSCEPTIBILITIES
            if(.not.lhfresponses) close(unit=iw)
            iw = iw+1000
            ! HF SUSCEPTIBILITIES
            close(unit=iw)
          end do
        end do
      end do
      ! RPA DIAGONALIZATION
      if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
        close (unit=1990)
        do i=1,nmaglayers
          close (unit=1990+i)
        end do
      end if
      if(ltestcharge) then
        close(unit=17964)
      end if
    end if

    return
  end subroutine openclose_chi_files

  ! This subroutine write all the susceptibilities into files
  ! (already opened with openclose_chi_files(1))
  ! Some information may be written on the screen
  subroutine write_susceptibilities(e)
    use mod_f90_kind, only: double
    use mod_parameters, only: lwriteonscreen, ltestcharge, lhfresponses, lnodiag, outputunit_loop, nmaglayers
    use mod_system, only: s => sys
    implicit none
    character(len=100)      :: varm
    integer                 :: i,j,iw,sigma
    real(double),intent(in) :: e

    if(lwriteonscreen) write(outputunit_loop,"(' #################  Susceptibilities:  #################')")
    do j=1, s%nAtoms
      do i=1, s%nAtoms
        do sigma=1, 4
          iw = 1000+(sigma-1)*s%nAtoms*s%nAtoms+(j-1)*s%nAtoms+i
          if(.not.lhfresponses) then
            write(unit=iw,fmt="(4(es16.9,2x))") e,real(schi(sigma,1,i,j)),aimag(schi(sigma,1,i,j)),abs(schi(sigma,1,i,j))
            if(((sigma==1).and.(i==j)).and.(lwriteonscreen)) write(outputunit_loop,"('E = ',es11.4,', Plane: ',i0,' Chi+- = (',es16.9,') + i(',es16.9,')')") e,i,real(schi(sigma,1,i,j)),aimag(schi(sigma,1,i,j))
          end if

          iw = iw+1000
          write(unit=iw,fmt="(4(es16.9,2x))") e,real(schihf(sigma,1,i,j)),aimag(schihf(sigma,1,i,j)),abs(schihf(sigma,1,i,j))
        end do
      end do
    end do

    if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
      write(varm,fmt="(a,i0,a)") '(',2*nmaglayers+1,'(es16.9,2x))'
      write(unit=1990,fmt=varm) e,(real(eval(i)),aimag(eval(i)),i=1,nmaglayers)
      do i=1,nmaglayers
        write(unit=1990+i,fmt=varm) e,(real(evecr(j,i)),aimag(evecr(j,i)),j=1,nmaglayers)
      end do
    end if ! nmaglayers

    if(ltestcharge) then
      do j=1, s%nAtoms
        do i=1, s%nAtoms
          write(unit=17964, fmt="(es16.9,' ',i0,' ',i0,' ',es16.9,' ',es16.9)") e, i, j, real(schi(2,2,i,j)+schi(3,2,i,j)-schi(2,3,i,j)-schi(3,3,i,j)), aimag(schi(2,2,i,j)+schi(3,2,i,j)-schi(2,3,i,j)-schi(3,3,i,j))
        end do
      end do
    end if

    return
  end subroutine write_susceptibilities

  ! This subroutine opens and closes all the files needed for the dc-limit susceptibilities
  subroutine openclose_dc_chi_files(iflag, count)
    use mod_parameters, only: dcfieldpart, lhfresponses, Npl_folder, eta, Utype, suffix, nmaglayers, lnodiag
    use mod_magnet, only: dcprefix, dcfield_dependence, dcfield, dc_header
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use mod_mpi_pars, only: abortProgram
    use EnergyIntegration, only: strEnergyParts
    implicit none
    integer, intent(in) :: iflag
    integer, intent(in) :: count
    character(len=500)  :: varm
    character(len=2)    :: spin(4)
    integer :: i,j,sigma,iw,err,errt=0

    spin(1) = "pm"
    spin(2) = "um"
    spin(3) = "dm"
    spin(4) = "mm"

    if(iflag==0) then
      do sigma=1,4
        do j=1,s%nAtoms
          do i=1,s%nAtoms
            iw = 10000+(sigma-1)*s%nAtoms*s%nAtoms+(j-1)*s%nAtoms+i
            ! RPA SUSCEPTIBILITIES
            if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a2,'/',a,'chi_',a,'_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
              open (unit=iw, file=varm, status='replace', form='formatted')
              write(unit=iw, fmt="('#',a,'  real part of chi ',a,'  ,  imaginary part of chi ',a,'  ,  amplitude of chi ',a,' ')") trim(dc_header),spin(sigma),spin(sigma),spin(sigma)
              close(unit=iw)
            end if
            iw = iw+1000
            ! HF SUSCEPTIBILITIES
            write(varm,"('./results/',a1,'SOC/',a,'/HF/',a2,'/',a,'chihf_',a,'_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
            open (unit=iw, file=varm, status='replace', form='formatted')
            write(unit=iw, fmt="('#',a,'  real part of chi ',a,' HF ,  imaginary part of chi ',a,' HF  ,  amplitude of chi ',a,' HF ')") trim(dc_header),spin(sigma),spin(sigma),spin(sigma)
            close(unit=iw)
          end do
        end do
      end do
      ! RPA DIAGONALIZATION
      if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
        write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/',a,'chi_eval_',a,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
        open (unit=19900, file=varm,status='replace', form='formatted')
        write(unit=19900,fmt="('#',a,'  real part of 1st eigenvalue  ,  imaginary part of 1st eigenvalue  ,  real part of 2nd eigenvalue  ,  imaginary part of 2nd eigenvalue  , ...')") trim(dc_header)
        close (unit=19900)
        do i=1,nmaglayers
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/',a,'chi_evec',i0,'_',a,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(dcprefix(count)),i,trim(dcfield(dcfield_dependence)),trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
          open (unit=19900+i, file=varm,status='replace', form='formatted')
          write(unit=19900+i,fmt="('#',a,'  real part of 1st component  ,  imaginary part of 1st component  ,  real part of 2nd component  ,  imaginary part of 2nd component  , ...')") trim(dc_header)
          close (unit=19900+i)
        end do
      end if

    else if(iflag==1) then
      do sigma=1,4
        do j=1,s%nAtoms
          do i=1,s%nAtoms
            iw = 10000+(sigma-1)*s%nAtoms*s%nAtoms+(j-1)*s%nAtoms+i
            ! RPA SUSCEPTIBILITIES
            if(.not.lhfresponses) then
              write(varm,"('./results/',a1,'SOC/',a,'/RPA/',a2,'/',a,'chi_',a,'_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
              open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
              errt = errt + err
            end if
            iw = iw+1000
            ! HF SUSCEPTIBILITIES
            write(varm,"('./results/',a1,'SOC/',a,'/HF/',a2,'/',a,'chihf_',a,'_',i0,'_',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),spin(sigma),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),i,j,trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
            open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
            errt = errt + err
          end do
        end do
      end do
      ! RPA DIAGONALIZATION
      if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
        write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/',a,'chi_eval_',a,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(dcprefix(count)),trim(dcfield(dcfield_dependence)),trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
        open (unit=19900, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        do i=1,nmaglayers
          write(varm,"('./results/',a1,'SOC/',a,'/RPA/pm/',a,'chi_evec',i0,'_',a,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(dcprefix(count)),i,trim(dcfield(dcfield_dependence)),trim(strEnergyParts),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(suffix)
          open (unit=19900+i, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end do
      end if
      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_dc_chi_files] Some file(s) do(es) not exist! Stopping before starting calculations...")

    else
      do sigma=1,4
        do j=1,s%nAtoms
          do i=1,s%nAtoms
            ! RPA SUSCEPTIBILITIES
            iw = 10000+(sigma-1) * s%nAtoms * s%nAtoms + (j-1) * s%nAtoms+i
            if(.not.lhfresponses) close(unit=iw)
            ! HF SUSCEPTIBILITIES
            iw = iw+1000
            close(unit=iw)
          end do
        end do
      end do
      ! RPA DIAGONALIZATION
      if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
        close (unit=19900)
        do i=1,nmaglayers
          close (unit=19900+i)
        end do
      end if
    end if

    return
  end subroutine openclose_dc_chi_files

  ! This subroutine write all the susceptibilities into files
  ! (already opened with openclose_chi_files(1))
  ! Some information is also written on the screen
  subroutine write_dc_susceptibilities()
    use mod_f90_kind, only: double
    use mod_parameters, only: lwriteonscreen, outputunit_loop, lhfresponses, nmaglayers, lnodiag
    use mod_magnet, only: hw_count, dc_fields, dcfield_dependence, dcfield
    use mod_system, only: s => sys
    implicit none
    character(len=100)      :: varm
    integer                 :: i,j,iw,sigma

    if(lwriteonscreen) write(outputunit_loop,"(' #################  Susceptibilities:  #################')")
    do j=1,s%nAtoms ;  do i=1,s%nAtoms ; do sigma=1,4
      iw = 10000+(sigma-1)*s%nAtoms*s%nAtoms+(j-1)*s%nAtoms+i
      if(.not.lhfresponses) then
        write(unit=iw,fmt="(a,2x,3(es16.9,2x))") trim(dc_fields(hw_count)) , real(schi(sigma,1,i,j)) , aimag(schi(sigma,1,i,j)) , abs(schi(sigma,1,i,j))
        if((sigma==1).and.(i==j).and.(lwriteonscreen)) write(outputunit_loop,"(a,' = ',a,', Plane: ',i0,' Chi+- = (',es16.9,') + i(',es16.9,')')") trim(dcfield(dcfield_dependence)),trim(dc_fields(hw_count)),i,real(schi(sigma,1,i,j)),aimag(schi(sigma,1,i,j))
      end if

      iw = iw+1000
      write(unit=iw,fmt="(a,2x,3(es16.9,2x))") trim(dc_fields(hw_count)) , real(schihf(sigma,1,i,j)) , aimag(schihf(sigma,1,i,j)) , abs(schihf(sigma,1,i,j))
    end do ; end do ; end do

    if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
      write(varm,fmt="(a,i0,a)") '(a,2x,',2*nmaglayers,'(es16.9,2x))'
      write(unit=19900,fmt=varm) trim(dc_fields(hw_count)) , (real(eval(i)) , aimag(eval(i)),i=1,nmaglayers)
      do i=1,nmaglayers
        write(unit=19900+i,fmt=varm) trim(dc_fields(hw_count)) , (real(evecr(j,i)) , aimag(evecr(j,i)),j=1,nmaglayers)
      end do
    end if ! nmaglayers

    return
  end subroutine write_dc_susceptibilities

  ! This subroutine sorts susceptibilities files
  subroutine sort_susceptibilities(count)
    use mod_f90_kind, only: double
    use mod_parameters, only: itype, lhfresponses, nmaglayers, lnodiag
    use mod_tools, only: sort_file
    use mod_system, only: s => sys
    implicit none
    integer, optional, intent(in) :: count
    integer :: i,j,iw,sigma,idc=1

    ! Opening chi and diag files
    if(itype==9) then
      idc=10
      call openclose_dc_chi_files(1, count)
    else
      call openclose_chi_files(1)
    end if

    do j=1, s%nAtoms
      do i=1, s%nAtoms
        do sigma=1, 4
          iw = 1000*idc+(sigma-1)*s%nAtoms*s%nAtoms+(j-1)*s%nAtoms+i
          if(.not.lhfresponses) then
            call sort_file(iw,.true.)
          end if

          iw = iw+1000
          call sort_file(iw,.true.)
        end do
      end do
    end do

    if((nmaglayers>1).and.(.not.lhfresponses).and.(.not.lnodiag)) then
      call sort_file(1990*idc,.true.)
      do i=1,nmaglayers
        call sort_file(1990*idc+i,.true.)
      end do
    end if ! nmaglayers

    ! Closing chi and diag files
    if(itype==9) then
      call openclose_dc_chi_files(2, count)
    else
      call openclose_chi_files(2)
    end if

    return
  end subroutine sort_susceptibilities

end module mod_susceptibilities
