module mod_susceptibilities
  use mod_f90_kind
  implicit none
  ! Spin susceptibilities
  complex(double),allocatable   :: schi(:,:,:,:),schihf(:,:,:,:)
  ! Rotation of spin susceptibilities
  logical                       :: lrot = .false.
  complex(double),allocatable   :: rotmat_i(:,:,:),rotmat_j(:,:,:),rottemp(:,:),schitemp(:,:),schirot(:,:)
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

contains

  ! This subroutine allocates variables related to the susceptibility calculation
  subroutine allocate_susceptibilities()
    use mod_f90_kind
    use mod_parameters, only: Npl,dim,nmaglayers,llinearsoc
    use mod_mpi_pars
    implicit none
    integer           :: AllocateStatus

    if(myrank_row.eq.0) then
      allocate( schi(4,4,Npl,Npl),schihf(4,4,Npl,Npl),chiorb(dim,dim), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        if(myrank.eq.0) write(*,"('[allocate_susceptibilities] Not enough memory for: schi,schihf,chiorb')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      if(lrot) allocate( rotmat_i(4,4,Npl),rotmat_j(4,4,Npl),rottemp(4,4),schitemp(4,4),schirot(4,4) )

      if(myrank_col.eq.0) then
        if(nmaglayers.gt.1) then
          lwork = 33*nmaglayers
          allocate( chimag(nmaglayers,nmaglayers),rwork(2*nmaglayers),eval(nmaglayers),evecl(1,nmaglayers),evecr(nmaglayers,nmaglayers),work(lwork) )
#ifdef _JUQUEEN
          allocate( dscale(nmaglayers),rconde(nmaglayers),rcondv(nmaglayers) )
#endif
        end if
      end if
    end if
    allocate( chiorb_hf(dim,dim), STAT = AllocateStatus )
    if (AllocateStatus.ne.0) then
      if(myrank.eq.0) write(*,"('[allocate_susceptibilities] Not enough memory for: chiorb_hf')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    if(llinearsoc) then
      allocate( chiorb_hflsoc(dim,dim), STAT = AllocateStatus  )
      if (AllocateStatus.ne.0) then
        write(*,"('[allocate_susceptibilities] Not enough memory for: chiorb_hflsoc')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if

    return
  end subroutine allocate_susceptibilities

  ! This subroutine allocates variables related to the susceptibility calculation
  subroutine deallocate_susceptibilities()
    use mod_f90_kind
    use mod_parameters, only: nmaglayers,llinearsoc
    use mod_mpi_pars
    implicit none

    if(myrank_row.eq.0) then
      deallocate(schi,schihf)
      if(lrot) deallocate(rotmat_i,rotmat_j,rottemp,schitemp,schirot)
      if(myrank_col.eq.0) then
        if(nmaglayers.gt.1) then
          deallocate(rwork,eval,evecl,evecr,work)
#ifdef _JUQUEEN
          deallocate(dscale,rconde,rcondv)
#endif
        end if
      end if
    end if
    if(allocated(chiorb)) deallocate(chiorb)
    deallocate(chiorb_hf)
    if(llinearsoc) deallocate(chiorb_hflsoc)

    return
  end subroutine deallocate_susceptibilities

  ! This subroutine diagonalize the transverse susceptibility
  subroutine diagonalize_susceptibilities()
    use mod_f90_kind
    use mod_parameters, only: nmaglayers,mmlayermag
    use mod_mpi_pars
    implicit none
    integer  :: i,j,ifail

    if(nmaglayers.gt.1) then

      do i=1,nmaglayers ; do j=1,nmaglayers
        chimag(i,j) = schi(1,1,mmlayermag(i)-1,mmlayermag(j)-1)
      end do ; end do

#ifdef _JUQUEEN
      call zgeevx('N','N','V','N',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,ilo,ihi,dscale,abnrm,rconde,rcondv,work,lwork,rwork,ifail)
#else
      call zgeev('N','V',nmaglayers,chimag,nmaglayers,eval,evecl,1,evecr,nmaglayers,work,lwork,rwork,ifail)
#endif

      if(ifail.ne.0) then
        write(*,*) '[diagonalize_susceptibilities] Problem with diagonalization. ifail = ',ifail
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
!         else
!           write(*,*) ' optimal lwork = ',work(1),' lwork = ',lwork
      end if
    end if ! nmaglayers

    return
  end subroutine diagonalize_susceptibilities

  ! This subroutine opens and closes all the files needed for the susceptibilities
  subroutine openclose_chi_files(iflag)
    use mod_parameters
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=2)    :: spin(4)
    character(len=1)    :: SOCc
    integer :: i,j,sigma,iw,iflag,err,errt=0

    spin(1) = "pm"
    spin(2) = "um"
    spin(3) = "dm"
    spin(4) = "mm"

    fieldpart = ""
    socpart   = ""
    if(SOC) then
      if(llinearsoc) then
        SOCc = "L"
      else
        SOCc = "T"
      end if
      write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
    else
      SOCc = "F"
    end if
    if(FIELD) then
      write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hwa,hwt,hwp
      if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
    end if

    if(iflag.eq.0) then
      do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
        ! RPA SUSCEPTIBILITIES
        iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/RPA/',a2,'/chi_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   energy  ,  real part of chi ',a,'  ,  imaginary part of chi ',a,'  ,  amplitude of chi ',a,'  ')") spin(sigma),spin(sigma),spin(sigma)
        close(unit=iw)
        ! HF SUSCEPTIBILITIES
        iw = iw+1000
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/HF/',a2,'/chihf_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   energy  ,  real part of chi ',a,' HF ,  imaginary part of chi ',a,' HF  ,  amplitude of chi ',a,' HF  ')") spin(sigma),spin(sigma),spin(sigma)
        close(unit=iw)
      end do ; end do ; end do
      ! RPA DIAGONALIZATION
      if(nmaglayers.gt.1) then
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/RPA/pm/chi_eval_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=1990, file=varm,status='unknown', form='formatted')
        write(unit=1990,fmt="('#   energy  ,  real part of 1st eigenvalue  ,  imaginary part of 1st eigenvalue  ,  real part of 2nd eigenvalue  ,  imaginary part of 2nd eigenvalue  , ... ')")
        close (unit=1990)
        do i=1,nmaglayers
          write(varm,"('./results/',a1,'SOC/',i0,'Npl/RPA/pm/chi_evec',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
          open (unit=1990+i, file=varm,status='unknown', form='formatted')
          write(unit=1990+i,fmt="('#   energy  ,  real part of 1st component  ,  imaginary part of 1st component  ,  real part of 2nd component  ,  imaginary part of 2nd component  , ...   ')")
          close (unit=1990+i)
        end do
      end if

    else if(iflag.eq.1) then
      do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
        ! RPA SUSCEPTIBILITIES
        iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/RPA/',a2,'/chi_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        ! HF SUSCEPTIBILITIES
        iw = iw+1000
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/HF/',a2,'/chihf_',i0,'_',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,spin(sigma),i,j,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do ; end do
      ! RPA DIAGONALIZATION
      if(nmaglayers.gt.1) then
        write(varm,"('./results/',a1,'SOC/',i0,'Npl/RPA/pm/chi_eval_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=1990, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        do i=1,nmaglayers
          write(varm,"('./results/',a1,'SOC/',i0,'Npl/RPA/pm/chi_evec',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,Npl,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart)
          open (unit=1990+i, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end do
      end if
      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(*,"(a,i0,a)") "[openclose_chi_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        stop
      end if

    else
      do sigma=1,4 ; do j=1,Npl ; do i=1,Npl
        ! RPA SUSCEPTIBILITIES
        iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i
        close(unit=iw)
        ! HF SUSCEPTIBILITIES
        iw = iw+1000
        close(unit=iw)
      end do ; end do ; end do
      ! RPA DIAGONALIZATION
      if(nmaglayers.gt.1) then
        close (unit=1990)
        do i=1,nmaglayers
          close (unit=1990+i)
        end do
      end if
    end if

    return
  end subroutine openclose_chi_files

  ! This subroutine write all the susceptibilities into files
  ! (already opened with openclose_chi_files(1))
  ! Some information is also written on the screen
  subroutine write_susceptibilities(e)
    use mod_f90_kind
    use mod_parameters, only: Npl,nmaglayers
    implicit none
    character(len=100)      :: varm
    integer                 :: i,j,iw,sigma
    real(double),intent(in) :: e

    write(*,"(' #################  Susceptibilities:  #################')")
    do j=1,Npl ;  do i=1,Npl ; do sigma=1,4
      iw = 1000+(sigma-1)*Npl*Npl+(j-1)*Npl+i

      write(unit=iw,fmt="(4(es16.9,2x))") e,real(schi(sigma,1,i,j)),aimag(schi(sigma,1,i,j)),abs(schi(sigma,1,i,j))
      if((sigma.eq.1).and.(i.eq.j)) write(*,"('E = ',es11.4,', Plane: ',i0,' Chi+- = (',es16.9,') + i(',es16.9,')')") e,i,real(schi(sigma,1,i,j)),aimag(schi(sigma,1,i,j))

      iw = iw+1000
      write(unit=iw,fmt="(4(es16.9,2x))") e,real(schihf(sigma,1,i,j)),aimag(schihf(sigma,1,i,j)),abs(schihf(sigma,1,i,j))
    end do ; end do ; end do

    if(nmaglayers.gt.1) then
      write(varm,fmt="(a,i0,a)") '(',2*nmaglayers+1,'(es16.9,2x))'
      write(unit=1990,fmt=varm) e,(real(eval(i)),aimag(eval(i)),i=1,nmaglayers)
      do i=1,nmaglayers
        write(unit=1990+i,fmt=varm) e,(real(evecr(j,i)),aimag(evecr(j,i)),j=1,nmaglayers)
      end do
    end if ! nmaglayers

    return
  end subroutine write_susceptibilities

end module mod_susceptibilities