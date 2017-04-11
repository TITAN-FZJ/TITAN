module mod_disturbances
  use mod_f90_kind
  implicit none
  ! Disturbances and renormalized disturbances: chd,sdx,sdy,sdz,ldx,ldy,ldz
  complex(double),allocatable   :: disturbances(:,:),rdisturbances(:,:)
  ! Full response functions
  complex(double), dimension(:,:),   allocatable :: tchiorbiikl
  complex(double), dimension(:,:,:), allocatable :: ldmat
  complex(double),dimension(:),allocatable       :: sdmat
contains

  ! This subroutine allocates variables related to the disturbance calculation
  subroutine allocate_disturbances()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_prefactors,only: prefactor,prefactorlsoc
    use mod_parameters, only: Npl,renorm,llinearsoc,dim,dimsigmaNpl,outputunit
    implicit none
    integer           :: AllocateStatus

    if(myrank_row==0) then
      allocate( disturbances(7,Npl),sdmat(dimsigmaNpl),ldmat(Npl,9,9), STAT = AllocateStatus )
      if (AllocateStatus/=0) then
        write(outputunit,"('[allocate_disturbances] Not enough memory for: disturbances,sdmat,ldmat')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      if(renorm) then
        allocate( rdisturbances(7,Npl), STAT = AllocateStatus )
        if (AllocateStatus/=0) then
          write(outputunit,"('[allocate_disturbances] Not enough memory for: rdisturbances')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if
    allocate( tchiorbiikl(dim,4), STAT = AllocateStatus  )
    if (AllocateStatus/=0) then
      write(outputunit,"('[allocate_disturbances] Not enough memory for: tchiorbiikl')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    if (.not. allocated(prefactor)) then
      allocate(prefactor(dim,dim), STAT = AllocateStatus  )
      if (AllocateStatus/=0) then
        write(outputunit,"('[allocate_disturbances] Not enough memory for: prefactor')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if
    if (.not. allocated(prefactorlsoc)) then
      if(llinearsoc) then
        allocate(prefactorlsoc(dim,dim), STAT = AllocateStatus  )
        if (AllocateStatus/=0) then
          write(outputunit,"('[allocate_disturbances] Not enough memory for: prefactorlsoc')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if

    return
  end subroutine allocate_disturbances

  ! This subroutine deallocates variables related to the disturbance calculation
  subroutine deallocate_disturbances()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_prefactors,only: prefactor,prefactorlsoc
    use mod_parameters, only: renorm
    implicit none

    if(myrank_row==0) then
      deallocate(disturbances,sdmat,ldmat)
      if(renorm) deallocate(rdisturbances)
    end if
    deallocate(tchiorbiikl)
    if (allocated(prefactor)) deallocate(prefactor)
    if (allocated(prefactorlsoc)) deallocate(prefactorlsoc)

    return
  end subroutine deallocate_disturbances

  ! This subroutine opens and closes all the files needed for the disturbances
  subroutine openclose_disturbance_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=5)    :: folder(7)
    character(len=2)    :: filename(7)
    character(len=1)    :: SOCc
    integer :: i,j,iw,iflag,err,errt=0

    fieldpart = ""
    socpart   = ""
    if(SOC) then
      if(llinearsoc) then
        SOCc = "L"
      else
        SOCc = "T"
      end if
      if(abs(socscale-1.d0)>1.d-6) write(socpart,"('_socscale=',f5.2)") socscale
    else
      SOCc = "F"
    end if
    if(lfield) then
      write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
      if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
      if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
      if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
      if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"
    end if

    folder(1) = "CD"
    folder(2) = "SD"
    folder(3) = "SD"
    folder(4) = "SD"
    folder(5) = "LD"
    folder(6) = "LD"
    folder(7) = "LD"
    if(lhfresponses) then
      do i=1,7
        folder(i) = trim(folder(i)) // "_HF"
      end do
    end if

    filename(1) = "Cd"
    filename(2) = "Sx"
    filename(3) = "Sy"
    filename(4) = "Sz"
    filename(5) = "Lx"
    filename(6) = "Ly"
    filename(7) = "Lz"

    if(iflag==0) then
      do i=1,Npl ; do j=1,7
        iw = 3000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/r',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,'_renormnb=',i0,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,renormnb,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do
    else if (iflag==1) then
      do i=1,Npl ; do j=1,7
        iw = 3000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/r',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,'_renormnb=',i0,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),filename(j),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,renormnb,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do
      ! Stop if some file does not exist
      if(errt/=0) then
        write(outputunit,"(a,i0,a)") "[openclose_disturbance_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      do i=1,Npl ; do j=1,7
        iw = 3000+(i-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do
    end if

    return
  end subroutine openclose_disturbance_files

  ! This subroutine write all the disturbances into files
  ! (already opened with openclose_disturbance_files(1))
  ! Some information may also be written on the screen
  subroutine write_disturbances(e)
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm,outputunit_loop,lwriteonscreen
    use mod_magnet, only: mvec_spherical
    implicit none
    integer  :: i,iw
    real(double),intent(in) :: e

    if(lwriteonscreen) write(outputunit_loop,"(' ################# Disturbances: #################')")
    ! Writing Spin, Charge and Orbital disturbances
    do i=1,Npl
      if(lwriteonscreen) then
        write(outputunit_loop,"('|--------------- Energy = ',es11.4,' , Plane: ',i0,' ---------------|')") e,i

        write(outputunit_loop,"('     Cd  = (',es16.9,') + i(',es16.9,')')") real(disturbances(1,i)),aimag(disturbances(1,i))
        write(outputunit_loop,"(' abs(Cd) = ',es16.9)") abs(disturbances(1,i))
        write(outputunit_loop,"('atan(Cd) = ',es16.9)") atan2(aimag(disturbances(1,i)),real(disturbances(1,i)))

        write(outputunit_loop,"('     Sdx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(2,i)),aimag(disturbances(2,i))
        write(outputunit_loop,"(' abs(Sdx) = ',es16.9)") abs(disturbances(2,i))
        write(outputunit_loop,"('atan(Sdx) = ',es16.9)") atan2(aimag(disturbances(2,i)),real(disturbances(2,i)))

        write(outputunit_loop,"('     Sdy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(3,i)),aimag(disturbances(3,i))
        write(outputunit_loop,"(' abs(Sdy) = ',es16.9)") abs(disturbances(3,i))
        write(outputunit_loop,"('atan(Sdy) = ',es16.9)") atan2(aimag(disturbances(3,i)),real(disturbances(3,i)))

        write(outputunit_loop,"('     Sdz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(4,i)),aimag(disturbances(4,i))
        write(outputunit_loop,"(' abs(Sdz) = ',es16.9)") abs(disturbances(4,i))
        write(outputunit_loop,"('atan(Sdz) = ',es16.9)") atan2(aimag(disturbances(4,i)),real(disturbances(4,i)))

        write(outputunit_loop,"('     Ldx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(5,i)),aimag(disturbances(5,i))
        write(outputunit_loop,"(' abs(Ldx) = ',es16.9)") abs(disturbances(5,i))
        write(outputunit_loop,"('atan(Ldx) = ',es16.9)") atan2(aimag(disturbances(5,i)),real(disturbances(5,i)))

        write(outputunit_loop,"('     Ldy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(6,i)),aimag(disturbances(6,i))
        write(outputunit_loop,"(' abs(Ldy) = ',es16.9)") abs(disturbances(6,i))
        write(outputunit_loop,"('atan(Ldy) = ',es16.9)") atan2(aimag(disturbances(6,i)),real(disturbances(6,i)))

        write(outputunit_loop,"('     Ldz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(7,i)),aimag(disturbances(7,i))
        write(outputunit_loop,"(' abs(Ldz) = ',es16.9)") abs(disturbances(7,i))
        write(outputunit_loop,"('atan(Ldz) = ',es16.9)") atan2(aimag(disturbances(7,i)),real(disturbances(7,i)))
      end if

      ! Writing charge disturbance
      iw = 3000+(i-1)*7
      write(unit=iw+1,fmt="(7(es16.9,2x))") e,abs(disturbances(1,i)),real(disturbances(1,i)),aimag(disturbances(1,i)),atan2(aimag(disturbances(1,i)),real(disturbances(1,i))),real(disturbances(1,i))/abs(disturbances(1,i)),aimag(disturbances(1,i))/abs(disturbances(1,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing x-component spin disturbance
      write(unit=iw+2,fmt="(7(es16.9,2x))") e,abs(disturbances(2,i)),real(disturbances(2,i)),aimag(disturbances(2,i)),atan2(aimag(disturbances(2,i)),real(disturbances(2,i))),real(disturbances(2,i))/abs(disturbances(2,i)),aimag(disturbances(2,i))/abs(disturbances(2,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing y-component spin disturbance
      write(unit=iw+3,fmt="(7(es16.9,2x))") e,abs(disturbances(3,i)),real(disturbances(3,i)),aimag(disturbances(3,i)),atan2(aimag(disturbances(3,i)),real(disturbances(3,i))),real(disturbances(3,i))/abs(disturbances(3,i)),aimag(disturbances(3,i))/abs(disturbances(3,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing z-component spin disturbance
      write(unit=iw+4,fmt="(7(es16.9,2x))") e,abs(disturbances(4,i)),real(disturbances(4,i)),aimag(disturbances(4,i)),atan2(aimag(disturbances(4,i)),real(disturbances(4,i))),real(disturbances(4,i))/abs(disturbances(4,i)),aimag(disturbances(4,i))/abs(disturbances(4,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)

      ! Writing x-component orbital disturbance
      write(unit=iw+5,fmt="(7(es16.9,2x))") e,abs(disturbances(5,i)),real(disturbances(5,i)),aimag(disturbances(5,i)),atan2(aimag(disturbances(5,i)),real(disturbances(5,i))),real(disturbances(5,i))/abs(disturbances(5,i)),aimag(disturbances(5,i))/abs(disturbances(5,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing y-component orbital disturbance
      write(unit=iw+6,fmt="(7(es16.9,2x))") e,abs(disturbances(6,i)),real(disturbances(6,i)),aimag(disturbances(6,i)),atan2(aimag(disturbances(6,i)),real(disturbances(6,i))),real(disturbances(6,i))/abs(disturbances(6,i)),aimag(disturbances(6,i))/abs(disturbances(6,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing z-component orbital disturbance
      write(unit=iw+7,fmt="(7(es16.9,2x))") e,abs(disturbances(7,i)),real(disturbances(7,i)),aimag(disturbances(7,i)),atan2(aimag(disturbances(7,i)),real(disturbances(7,i))),real(disturbances(7,i))/abs(disturbances(7,i)),aimag(disturbances(7,i))/abs(disturbances(7,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)

      ! Writing renormalized disturbances
      if(renorm) then
        ! Writing renormalized charge disturbance
        write(unit=iw+1001,fmt="(7(es16.9,2x))") e,abs(rdisturbances(1,i)),real(rdisturbances(1,i)),aimag(rdisturbances(1,i)),atan2(aimag(rdisturbances(1,i)),real(rdisturbances(1,i))),real(rdisturbances(1,i))/abs(rdisturbances(1,i)),aimag(rdisturbances(1,i))/abs(rdisturbances(1,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized x-component spin disturbance
        write(unit=iw+1002,fmt="(7(es16.9,2x))") e,abs(rdisturbances(2,i)),real(rdisturbances(2,i)),aimag(rdisturbances(2,i)),atan2(aimag(rdisturbances(2,i)),real(rdisturbances(2,i))),real(rdisturbances(2,i))/abs(rdisturbances(2,i)),aimag(rdisturbances(2,i))/abs(rdisturbances(2,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized y-component spin disturbance
        write(unit=iw+1003,fmt="(7(es16.9,2x))") e,abs(rdisturbances(3,i)),real(rdisturbances(3,i)),aimag(rdisturbances(3,i)),atan2(aimag(rdisturbances(3,i)),real(rdisturbances(3,i))),real(rdisturbances(3,i))/abs(rdisturbances(3,i)),aimag(rdisturbances(3,i))/abs(rdisturbances(3,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized z-component spin disturbance
        write(unit=iw+1004,fmt="(7(es16.9,2x))") e,abs(rdisturbances(4,i)),real(rdisturbances(4,i)),aimag(rdisturbances(4,i)),atan2(aimag(rdisturbances(4,i)),real(rdisturbances(4,i))),real(rdisturbances(4,i))/abs(rdisturbances(4,i)),aimag(rdisturbances(4,i))/abs(rdisturbances(4,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)

        ! Writing renormalized x-component orbital disturbance
        write(unit=iw+1005,fmt="(7(es16.9,2x))") e,abs(rdisturbances(5,i)),real(rdisturbances(5,i)),aimag(rdisturbances(5,i)),atan2(aimag(rdisturbances(5,i)),real(rdisturbances(5,i))),real(rdisturbances(5,i))/abs(rdisturbances(5,i)),aimag(rdisturbances(5,i))/abs(rdisturbances(5,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized y-component orbital disturbance
        write(unit=iw+1006,fmt="(7(es16.9,2x))") e,abs(rdisturbances(6,i)),real(rdisturbances(6,i)),aimag(rdisturbances(6,i)),atan2(aimag(rdisturbances(6,i)),real(rdisturbances(6,i))),real(rdisturbances(6,i))/abs(rdisturbances(6,i)),aimag(rdisturbances(6,i))/abs(rdisturbances(6,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized z-component orbital disturbance
        write(unit=iw+1007,fmt="(7(es16.9,2x))") e,abs(rdisturbances(7,i)),real(rdisturbances(7,i)),aimag(rdisturbances(7,i)),atan2(aimag(rdisturbances(7,i)),real(rdisturbances(7,i))),real(rdisturbances(7,i))/abs(rdisturbances(7,i)),aimag(rdisturbances(7,i))/abs(rdisturbances(7,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      end if
    end do

    return
  end subroutine write_disturbances

  ! This subroutine opens and closes all the files needed for the dc-limit disturbances
  subroutine openclose_dc_disturbance_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=5)    :: folder(7)
    character(len=2)    :: filename(7)
    character(len=1)    :: SOCc
    integer :: i,j,iw,iflag,err,errt=0

    fieldpart = ""
    socpart   = ""
    if(SOC) then
      if(llinearsoc) then
        SOCc = "L"
      else
        SOCc = "T"
      end if
      if(abs(socscale-1.d0)>1.d-6) write(socpart,"('_socscale=',f5.2)") socscale
    else
      SOCc = "F"
    end if

    if(dcfield_dependence/=7) then
      if((dcfield_dependence/=1).and.(dcfield_dependence/=4).and.(dcfield_dependence/=5)) write(fieldpart,"(a,'_hwa=',es9.2)") trim(fieldpart),hwa
      if((dcfield_dependence/=2).and.(dcfield_dependence/=4).and.(dcfield_dependence/=6)) write(fieldpart,"(a,'_hwt=',f5.2)") trim(fieldpart),hwt
      if((dcfield_dependence/=3).and.(dcfield_dependence/=5).and.(dcfield_dependence/=6)) write(fieldpart,"(a,'_hwp=',f5.2)") trim(fieldpart),hwp
    end if
    if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
    if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
    if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
    if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"

    folder(1) = "CD"
    folder(2) = "SD"
    folder(3) = "SD"
    folder(4) = "SD"
    folder(5) = "LD"
    folder(6) = "LD"
    folder(7) = "LD"
    if(lhfresponses) then
      do i=1,7
        folder(i) = trim(folder(i)) // "_HF"
      end do
    end if

    filename(1) = "Cd"
    filename(2) = "Sx"
    filename(3) = "Sy"
    filename(4) = "Sz"
    filename(5) = "Lx"
    filename(6) = "Ly"
    filename(7) = "Lz"

    if(iflag==0) then
      do i=1,Npl ; do j=1,7
        iw = 30000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,' imag part of ',a,' , real part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'r',a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,'_renormnb=',i0,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,renormnb,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#',a,' imag part of ',a,' , real part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do
    else if (iflag==1) then
      do i=1,Npl ; do j=1,7
        iw = 30000+(i-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'r',a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_EFp=',es8.1,'_EFt=',es8.1,'_renormnb=',i0,a,'.dat')") SOCc,trim(Npl_folder),trim(folder(j)),trim(dcprefix(count)),filename(j),trim(dcfield(dcfield_dependence)),i,parts,parts3,nkpt,eta,Utype,trim(fieldpart),trim(socpart),EFp,Eft,renormnb,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do
      ! Stop if some file does not exist
      if(errt/=0) then
        write(outputunit,"(a,i0,a)") "[openclose_dc_disturbance_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      do i=1,Npl ; do j=1,7
        iw = 30000+(i-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do
    end if

    return
  end subroutine openclose_dc_disturbance_files

  ! This subroutine write all the dc-limit disturbances into files
  ! (already opened with openclose_dc_disturbance_files(1))
  ! Some information may also be written on the screen
  subroutine write_dc_disturbances()
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm,dcfield,dcfield_dependence,dc_fields,hw_count,outputunit_loop,lwriteonscreen
    use mod_magnet, only: mvec_spherical
    implicit none
    integer  :: i,iw

    if(lwriteonscreen) write(outputunit_loop,"(' ################# Disturbances: #################')")
    ! Writing Spin, Charge and Orbital disturbances
    do i=1,Npl
      if(lwriteonscreen) then
        write(outputunit_loop,"('|--------------- ',a,' = ',a,' , Plane: ',i0,' ---------------|')") trim(dcfield(dcfield_dependence)),trim(dc_fields(hw_count)),i

        write(outputunit_loop,"('     Cd  = (',es16.9,') + i(',es16.9,')')") real(disturbances(1,i)),aimag(disturbances(1,i))
        write(outputunit_loop,"(' abs(Cd) = ',es16.9)") abs(disturbances(1,i))
        write(outputunit_loop,"('atan(Cd) = ',es16.9)") atan2(aimag(disturbances(1,i)),real(disturbances(1,i)))

        write(outputunit_loop,"('     Sdx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(2,i)),aimag(disturbances(2,i))
        write(outputunit_loop,"(' abs(Sdx) = ',es16.9)") abs(disturbances(2,i))
        write(outputunit_loop,"('atan(Sdx) = ',es16.9)") atan2(aimag(disturbances(2,i)),real(disturbances(2,i)))

        write(outputunit_loop,"('     Sdy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(3,i)),aimag(disturbances(3,i))
        write(outputunit_loop,"(' abs(Sdy) = ',es16.9)") abs(disturbances(3,i))
        write(outputunit_loop,"('atan(Sdy) = ',es16.9)") atan2(aimag(disturbances(3,i)),real(disturbances(3,i)))

        write(outputunit_loop,"('     Sdz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(4,i)),aimag(disturbances(4,i))
        write(outputunit_loop,"(' abs(Sdz) = ',es16.9)") abs(disturbances(4,i))
        write(outputunit_loop,"('atan(Sdz) = ',es16.9)") atan2(aimag(disturbances(4,i)),real(disturbances(4,i)))

        write(outputunit_loop,"('     Ldx  = (',es16.9,') + i(',es16.9,')')") real(disturbances(5,i)),aimag(disturbances(5,i))
        write(outputunit_loop,"(' abs(Ldx) = ',es16.9)") abs(disturbances(5,i))
        write(outputunit_loop,"('atan(Ldx) = ',es16.9)") atan2(aimag(disturbances(5,i)),real(disturbances(5,i)))

        write(outputunit_loop,"('     Ldy  = (',es16.9,') + i(',es16.9,')')") real(disturbances(6,i)),aimag(disturbances(6,i))
        write(outputunit_loop,"(' abs(Ldy) = ',es16.9)") abs(disturbances(6,i))
        write(outputunit_loop,"('atan(Ldy) = ',es16.9)") atan2(aimag(disturbances(6,i)),real(disturbances(6,i)))

        write(outputunit_loop,"('     Ldz  = (',es16.9,') + i(',es16.9,')')") real(disturbances(7,i)),aimag(disturbances(7,i))
        write(outputunit_loop,"(' abs(Ldz) = ',es16.9)") abs(disturbances(7,i))
        write(outputunit_loop,"('atan(Ldz) = ',es16.9)") atan2(aimag(disturbances(7,i)),real(disturbances(7,i)))
      end if

      ! Writing charge disturbance
      iw = 30000+(i-1)*7
      write(unit=iw+1,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(disturbances(1,i)) , real(disturbances(1,i)) , atan2(aimag(disturbances(1,i)),real(disturbances(1,i))) , real(disturbances(1,i))/abs(disturbances(1,i)) , aimag(disturbances(1,i))/abs(disturbances(1,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing x-component spin disturbance
      write(unit=iw+2,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(disturbances(2,i)) , real(disturbances(2,i)) , atan2(aimag(disturbances(2,i)),real(disturbances(2,i))) , real(disturbances(2,i))/abs(disturbances(2,i)) , aimag(disturbances(2,i))/abs(disturbances(2,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing y-component spin disturbance
      write(unit=iw+3,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(disturbances(3,i)) , real(disturbances(3,i)) , atan2(aimag(disturbances(3,i)),real(disturbances(3,i))) , real(disturbances(3,i))/abs(disturbances(3,i)) , aimag(disturbances(3,i))/abs(disturbances(3,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing z-component spin disturbance
      write(unit=iw+4,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(disturbances(4,i)) , real(disturbances(4,i)) , atan2(aimag(disturbances(4,i)),real(disturbances(4,i))) , real(disturbances(4,i))/abs(disturbances(4,i)) , aimag(disturbances(4,i))/abs(disturbances(4,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)

      ! Writing x-component orbital disturbance
      write(unit=iw+5,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(disturbances(5,i)) , real(disturbances(5,i)) , atan2(aimag(disturbances(5,i)),real(disturbances(5,i))) , real(disturbances(5,i))/abs(disturbances(5,i)) , aimag(disturbances(5,i))/abs(disturbances(5,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing y-component orbital disturbance
      write(unit=iw+6,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(disturbances(6,i)) , real(disturbances(6,i)) , atan2(aimag(disturbances(6,i)),real(disturbances(6,i))) , real(disturbances(6,i))/abs(disturbances(6,i)) , aimag(disturbances(6,i))/abs(disturbances(6,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      ! Writing z-component orbital disturbance
      write(unit=iw+7,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(disturbances(7,i)) , real(disturbances(7,i)) , atan2(aimag(disturbances(7,i)),real(disturbances(7,i))) , real(disturbances(7,i))/abs(disturbances(7,i)) , aimag(disturbances(7,i))/abs(disturbances(7,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)

      ! Writing renormalized disturbances
      if(renorm) then
        ! Writing renormalized charge disturbance
        write(unit=iw+1001,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rdisturbances(1,i)) , real(rdisturbances(1,i)) , atan2(aimag(rdisturbances(1,i)),real(rdisturbances(1,i))) , real(rdisturbances(1,i))/abs(rdisturbances(1,i)) , aimag(rdisturbances(1,i))/abs(rdisturbances(1,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized x-component spin disturbance
        write(unit=iw+1002,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rdisturbances(2,i)) , real(rdisturbances(2,i)) , atan2(aimag(rdisturbances(2,i)),real(rdisturbances(2,i))) , real(rdisturbances(2,i))/abs(rdisturbances(2,i)) , aimag(rdisturbances(2,i))/abs(rdisturbances(2,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized y-component spin disturbance
        write(unit=iw+1003,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rdisturbances(3,i)) , real(rdisturbances(3,i)) , atan2(aimag(rdisturbances(3,i)),real(rdisturbances(3,i))) , real(rdisturbances(3,i))/abs(rdisturbances(3,i)) , aimag(rdisturbances(3,i))/abs(rdisturbances(3,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized z-component spin disturbance
        write(unit=iw+1004,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rdisturbances(4,i)) , real(rdisturbances(4,i)) , atan2(aimag(rdisturbances(4,i)),real(rdisturbances(4,i))) , real(rdisturbances(4,i))/abs(rdisturbances(4,i)) , aimag(rdisturbances(4,i))/abs(rdisturbances(4,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)

        ! Writing renormalized x-component orbital disturbance
        write(unit=iw+1005,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rdisturbances(5,i)) , real(rdisturbances(5,i)) , atan2(aimag(rdisturbances(5,i)),real(rdisturbances(5,i))) , real(rdisturbances(5,i))/abs(rdisturbances(5,i)) , aimag(rdisturbances(5,i))/abs(rdisturbances(5,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized y-component orbital disturbance
        write(unit=iw+1006,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rdisturbances(6,i)) , real(rdisturbances(6,i)) , atan2(aimag(rdisturbances(6,i)),real(rdisturbances(6,i))) , real(rdisturbances(6,i))/abs(rdisturbances(6,i)) , aimag(rdisturbances(6,i))/abs(rdisturbances(6,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
        ! Writing renormalized z-component orbital disturbance
        write(unit=iw+1007,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(rdisturbances(7,i)) , real(rdisturbances(7,i)) , atan2(aimag(rdisturbances(7,i)),real(rdisturbances(7,i))) , real(rdisturbances(7,i))/abs(rdisturbances(7,i)) , aimag(rdisturbances(7,i))/abs(rdisturbances(7,i)) , mvec_spherical(i,2) , mvec_spherical(i,3)
      end if
    end do

    return
  end subroutine write_dc_disturbances

  ! This subroutine sorts disturbance files
  subroutine sort_disturbances()
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm,itype
    use mod_tools, only: sort_file
    implicit none
    integer :: i,j,iw,idc=1

    ! Opening disturbance files
    if(itype==9) then
      idc=10
      call openclose_dc_disturbance_files(1)
    else
      call openclose_disturbance_files(1)
    end if

    do i=1,Npl
      ! Sorting disturbance files
      iw = 3000*idc+(i-1)*7
      do j=1,7
        call sort_file(iw+j,.true.)
      end do

      ! Sorting renormalized disturbances
      if(renorm) then
        do j=1001,1007
          call sort_file(iw+j,.true.)
        end do
      end if
    end do

    ! Closing disturbance files
    if(itype==9) then
      call openclose_dc_disturbance_files(2)
    else
      call openclose_disturbance_files(2)
    end if

    return
  end subroutine sort_disturbances

end module mod_disturbances