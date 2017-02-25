module mod_beff
  use mod_f90_kind
  implicit none
  ! Effective field
  complex(double),dimension(:,:),allocatable :: chiinv
  complex(double),dimension(:)  ,allocatable :: Beff,Beff_cart
contains

  ! This subroutine allocates variables related to the effective field calculation
  subroutine allocate_beff()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters, only: dimsigmaNpl,outputunit
    implicit none
    integer :: AllocateStatus

    if(myrank_row.eq.0) then
      allocate( Beff(dimsigmaNpl),Beff_cart(dimsigmaNpl),chiinv(dimsigmaNpl,dimsigmaNpl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        write(outputunit,"('[allocate_beff] Not enough memory for: Beff,Beff_cart,chiinv')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if

    return
  end subroutine allocate_beff

  ! This subroutine deallocates variables related to the effective field calculation
  subroutine deallocate_beff()
    use mod_f90_kind
    use mod_mpi_pars
    implicit none

    if(myrank_row.eq.0) then
      deallocate(Beff,Beff_cart,chiinv)
    end if

    return
  end subroutine deallocate_beff

  ! This subroutine opens and closes all the files needed for the effective field
  subroutine openclose_beff_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=7)    :: folder
    character(len=4)    :: filename
    character(len=1)    :: direction(4),SOCc
    integer :: i,sigma,iw,iflag,err,errt=0

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
    if(lfield) then
      write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
      if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
      if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
      if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
      if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"
    end if

    folder = "Beff"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename = "Beff"

    direction(1) = "0"
    direction(2) = "x"
    direction(3) = "y"
    direction(4) = "z"

    if(iflag.eq.0) then
      do sigma=1,4 ; do i=1,Npl
        iw = 8000+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma)
        close(unit=iw)
      end do ; end do
    else if (iflag.eq.1) then
      do sigma=1,4 ; do i=1,Npl
        iw = 8000+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do
      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(outputunit,"(a,i0,a)") "[openclose_beff_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      do sigma=1,4 ; do i=1,Npl
        iw = 8000+(sigma-1)*Npl+i
        close(unit=iw)
      end do ; end do
    end if

    return
  end subroutine openclose_beff_files

  ! This subroutine write all the effective fields into files
  ! (already opened with openclose_beff_files(1))
  ! Some information may also be written on the screen
  subroutine write_beff(e)
    use mod_f90_kind
    use mod_parameters, only: Npl,sigmai2i
    use mod_magnet, only: mvec_spherical
    implicit none
    integer      :: i,iw,sigma
    real(double) :: phase,sine,cosine
    real(double),intent(in) :: e

    do sigma=1,4 ; do i=1,Npl
      iw = 8000+(sigma-1)*Npl+i

      if(abs(Beff_cart(sigmai2i(sigma,i))).ge.1.d-10) then
        phase  = atan2(aimag(Beff_cart(sigmai2i(sigma,i))),real(Beff_cart(sigmai2i(sigma,i))))
        sine   = real(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i)))
        cosine = aimag(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i)))
      else
        phase  = 0.d0
        sine   = 0.d0
        cosine = 0.d0
      end if

      write(unit=iw,fmt="(9(es16.9,2x))") e , abs(Beff_cart(sigmai2i(sigma,i))) , real(Beff_cart(sigmai2i(sigma,i))) , aimag(Beff_cart(sigmai2i(sigma,i))) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)
    end do ; end do

    return
  end subroutine write_beff

  ! This subroutine opens and closes all the files needed for the effective field
  subroutine openclose_dc_beff_files(iflag)
    use mod_parameters
    use mod_mpi_pars
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=7)    :: folder
    character(len=4)    :: filename
    character(len=1)    :: direction(4),SOCc
    integer :: i,sigma,iw,iflag,err,errt=0

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

    if(dcfield_dependence.ne.7) then
      if((dcfield_dependence.ne.1).and.(dcfield_dependence.ne.4).and.(dcfield_dependence.ne.5)) write(fieldpart,"(a,'_hwa=',es9.2)") trim(fieldpart),hwa
      if((dcfield_dependence.ne.2).and.(dcfield_dependence.ne.4).and.(dcfield_dependence.ne.6)) write(fieldpart,"(a,'_hwt=',f5.2)") trim(fieldpart),hwt
      if((dcfield_dependence.ne.3).and.(dcfield_dependence.ne.5).and.(dcfield_dependence.ne.6)) write(fieldpart,"(a,'_hwp=',f5.2)") trim(fieldpart),hwp
    end if
    if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
    if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
    if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
    if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"

    folder = "Beff"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename = "Beff"

    direction(1) = "0"
    direction(2) = "x"
    direction(3) = "y"
    direction(4) = "z"

    if(iflag.eq.0) then
      do sigma=1,4 ; do i=1,Npl
        iw = 80000+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#',a,' imaginary part of ',a,a,' , real part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma)
        close(unit=iw)
      end do ; end do
    else if (iflag.eq.1) then
      do sigma=1,4 ; do i=1,Npl
        iw = 80000+(sigma-1)*Npl+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename),direction(sigma),trim(dcfield(dcfield_dependence)),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do
      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(outputunit,"(a,i0,a)") "[openclose_dc_beff_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    else
      do sigma=1,4 ; do i=1,Npl
        iw = 80000+(sigma-1)*Npl+i
        close(unit=iw)
      end do ; end do
    end if

    return
  end subroutine openclose_dc_beff_files

  ! This subroutine write all the effective fields into files
  ! (already opened with openclose_dc_beff_files(1))
  ! Some information may also be written on the screen
  subroutine write_dc_beff()
    use mod_f90_kind
    use mod_parameters, only: Npl,sigmai2i,dc_fields,hw_count
    use mod_magnet, only: mvec_spherical
    implicit none
    integer      :: i,iw,sigma
    real(double) :: phase,sine,cosine

    do sigma=1,4 ; do i=1,Npl
      iw = 80000+(sigma-1)*Npl+i

      if(abs(Beff_cart(sigmai2i(sigma,i))).ge.1.d-10) then
        phase  = atan2(aimag(Beff_cart(sigmai2i(sigma,i))),real(Beff_cart(sigmai2i(sigma,i))))
        sine   = real(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i)))
        cosine = aimag(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i)))
      else
        phase  = 0.d0
        sine   = 0.d0
        cosine = 0.d0
      end if

      write(unit=iw,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , aimag(Beff_cart(sigmai2i(sigma,i))) , real(Beff_cart(sigmai2i(sigma,i))) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)
    end do ; end do

    return
  end subroutine write_dc_beff

  ! This subroutine sorts effective field files
  subroutine sort_beff()
    use mod_f90_kind
    use mod_parameters, only: Npl,itype
    use mod_tools, only: sort_file
    implicit none
    integer :: i,sigma,iw,idc=1

    ! Opening effective field files
    if(itype.eq.9) then
      idc=10
      call openclose_dc_beff_files(1)
    else
      call openclose_beff_files(1)
    end if

    do sigma=1,4 ; do i=1,Npl
      iw = 8000*idc+(sigma-1)*Npl+i
      ! Sorting effective field files
      call sort_file(iw,.true.)
    end do ; end do

    ! Closing effective field files
    if(itype.eq.9) then
      call openclose_dc_beff_files(2)
    else
      call openclose_beff_files(2)
    end if

    return
  end subroutine sort_beff

end module mod_beff