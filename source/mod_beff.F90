module mod_beff
  use mod_f90_kind
  implicit none
  ! Effective field
  complex(double),dimension(:,:),allocatable :: chiinv
  complex(double),dimension(:),allocatable   :: Beff,Beff_cart
contains

  ! This subroutine allocates variables related to the effective field calculation
  subroutine allocate_beff()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_parameters, only: dimsigmaNpl
    implicit none
    integer           :: AllocateStatus

    if(myrank_row.eq.0) then
      allocate( Beff(dimsigmaNpl),Beff_cart(dimsigmaNpl),chiinv(dimsigmaNpl,dimsigmaNpl), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        if(myrank.eq.0) write(*,"('[allocate_beff] Not enough memory for: Beff,Beff_cart,chiinv')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if

    return
  end subroutine allocate_beff

  ! This subroutine allocates variables related to the effective field calculation
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
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=7)    :: folder
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
    if(FIELD) then
      write(fieldpart,"('_hwx=',e8.1,'_hwy=',e8.1,'_hwz=',e8.1)") hwx,hwy,hwz
      if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
    end if

    folder = "Beff"
    if(lnorenorm) folder = trim(folder) // "_HF"

    direction(1) = "0"
    direction(2) = "x"
    direction(3) = "y"
    direction(4) = "z"

    if(iflag.eq.0) then
      do sigma=1,4 ; do i=1,Npl
        iw = 7000+(sigma-1)*Npl+i
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/Beff',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   energy  , amplitude of Beff_',a,' , real part of Beff_',a,' , imaginary part of Beff_',a,' , phase of Beff_',a,' , cosine of Beff_',a,'  ,  sine of ',a,'  ')") direction(sigma),direction(sigma),direction(sigma),direction(sigma),direction(sigma),direction(sigma)
        close(unit=iw)
      end do ; end do
    else if (iflag.eq.1) then
      do sigma=1,4 ; do i=1,Npl
        iw = 7000+(sigma-1)*Npl+i
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/Beff',a,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',e8.1,'_Utype=',i0,a,a,'_dirEfield=',A,'.dat')") SOCc,Npl,trim(folder),direction(sigma),i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
      end do ; end do
      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(*,"(a,i0,a)") "[openclose_beff_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        stop
      end if
    else
      do sigma=1,4 ; do i=1,Npl
        iw = 7000+(sigma-1)*Npl+i
        close(unit=iw)
      end do ; end do
    end if

    return
  end subroutine openclose_beff_files

  ! This subroutine write all the effective fields into files
  ! (already opened with openclose_beff_files(1))
  ! Some information is also written on the screen
  subroutine write_beff(e)
    use mod_f90_kind
    use mod_parameters, only: Npl,sigmai2i
    implicit none
    integer  :: i,iw,sigma
    real(double),intent(in) :: e

    do sigma=1,4 ; do i=1,Npl
      iw = 7000+(sigma-1)*Npl+i

      write(unit=iw,fmt="(7(e16.9,2x))") e,abs(Beff_cart(sigmai2i(sigma,i))),real(Beff_cart(sigmai2i(sigma,i))),aimag(Beff_cart(sigmai2i(sigma,i))),atan2(aimag(Beff_cart(sigmai2i(sigma,i))),real(Beff_cart(sigmai2i(sigma,i)))),real(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i))),aimag(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i)))
    end do ; end do

    return
  end subroutine write_beff

end module mod_beff