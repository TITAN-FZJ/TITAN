module mod_beff
  use mod_f90_kind, only: double
  implicit none
  ! Effective field
  complex(double),dimension(:,:), allocatable :: chiinv
  complex(double),dimension(:), allocatable :: Beff,Beff_cart

contains

  ! This subroutine allocates variables related to the effective field calculation
  subroutine allocate_beff()
    use mod_f90_kind, only: double
    use mod_mpi_pars, only: abortProgram, myrank_row
    use mod_parameters, only: dimsigmaNpl,outputunit
    implicit none
    integer :: AllocateStatus

    if(myrank_row==0) then
      allocate( Beff(dimsigmaNpl),Beff_cart(dimsigmaNpl),chiinv(dimsigmaNpl,dimsigmaNpl), STAT = AllocateStatus )
      if (AllocateStatus /= 0) call abortProgram("[allocate_beff] Not enough memory for: Beff,Beff_cart,chiinv")
    end if

    return
  end subroutine allocate_beff

  ! This subroutine deallocates variables related to the effective field calculation
  subroutine deallocate_beff()
    use mod_mpi_pars, only: myrank_row
    implicit none

    if(myrank_row==0) then
      deallocate(Beff,Beff_cart,chiinv)
    end if

    return
  end subroutine deallocate_beff

  ! This subroutine opens and closes all the files needed for the effective field
  subroutine openclose_beff_files(iflag)
    use mod_parameters, only: suffix, fieldpart, energypart, lhfresponses, Npl_folder, eta, Utype
    use mod_mpi_pars, only: abortProgram
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use electricfield, only: strElectricField
    implicit none

    character(len=500)  :: varm
    character(len=7)    :: folder
    character(len=4)    :: filename
    character(len=1)    :: direction(4)
    integer :: i,sigma,iw,iflag,err,errt=0

    folder = "Beff"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename = "Beff"

    direction(1) = "0"
    direction(2) = "x"
    direction(3) = "y"
    direction(4) = "z"

    if(iflag==0) then
      do sigma=1,4 ; do i=1, s%nAtoms
        iw = 8000+(sigma-1)*s%nAtoms+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename),direction(sigma),i,trim(energypart),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,a,' , real part of ',a,a,' , imaginary part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  ')") trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma)
        close(unit=iw)
      end do ; end do
    else if (iflag==1) then
      do sigma=1,4 ; do i=1,s%nAtoms
        iw = 8000+(sigma-1)*s%nAtoms+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(filename),direction(sigma),i,trim(energypart),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
      end do ; end do
      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_beff_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
    else
      do sigma=1,4 ; do i=1,s%nAtoms
        iw = 8000+(sigma-1)*s%nAtoms+i
        close(unit=iw)
      end do ; end do
    end if

    return
  end subroutine openclose_beff_files

  ! This subroutine write all the effective fields into files
  ! (already opened with openclose_beff_files(1))
  ! Some information may also be written on the screen
  subroutine write_beff(e)
    use mod_f90_kind, only: double
    use mod_parameters, only: sigmai2i
    use mod_magnet, only: mvec_spherical
    use mod_system, only: s => sys
    implicit none
    integer      :: i,iw,sigma
    real(double) :: phase,sine,cosine
    real(double),intent(in) :: e

    do sigma=1,4
      do i=1,s%nAtoms
        iw = 8000+(sigma-1)*s%nAtoms+i

        if(abs(Beff_cart(sigmai2i(sigma,i)))>=1.d-10) then
          phase  = atan2(aimag(Beff_cart(sigmai2i(sigma,i))),real(Beff_cart(sigmai2i(sigma,i))))
          sine   = real(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i)))
          cosine = aimag(Beff_cart(sigmai2i(sigma,i)))/abs(Beff_cart(sigmai2i(sigma,i)))
        else
          phase  = 0.d0
          sine   = 0.d0
          cosine = 0.d0
        end if

        write(unit=iw,fmt="(9(es16.9,2x))") e , abs(Beff_cart(sigmai2i(sigma,i))) , real(Beff_cart(sigmai2i(sigma,i))) , aimag(Beff_cart(sigmai2i(sigma,i))) , phase , sine , cosine , mvec_spherical(i,2) , mvec_spherical(i,3)
      end do
    end do

    return
  end subroutine write_beff

  ! This subroutine opens and closes all the files needed for the effective field
  subroutine openclose_dc_beff_files(iflag)
    use mod_parameters, only: lhfresponses,dcfieldpart, energypart, suffix, Npl_folder, Utype, eta
    use mod_magnet, only: dcprefix, dc_header, dcfield, dcfield_dependence
    use mod_mpi_pars, only: abortProgram
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use electricfield, only: strElectricField
    implicit none

    character(len=500)  :: varm
    character(len=7)    :: folder
    character(len=4)    :: filename
    character(len=1)    :: direction(4)
    integer :: i,sigma,iw,iflag,err,errt=0, count

    folder = "Beff"
    if(lhfresponses) folder = trim(folder) // "_HF"
    filename = "Beff"

    direction(1) = "0"
    direction(2) = "x"
    direction(3) = "y"
    direction(4) = "z"

    if(iflag==0) then
      do sigma=1,4 ; do i=1,s%nAtoms
        iw = 80000+(sigma-1)*s%nAtoms+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(energypart),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,' imaginary part of ',a,a,' , real part of ',a,a,' , phase of ',a,a,' , cosine of ',a,a,'  ,  sine of ',a,a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma),trim(filename),direction(sigma)
        close(unit=iw)
      end do ; end do
    else if (iflag==1) then
      do sigma=1,4 ; do i=1,s%nAtoms
        iw = 80000+(sigma-1)*s%nAtoms+i
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,a,'_',a,'_pos=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,a,'.dat')") SOCc,trim(Npl_folder),trim(folder),trim(dcprefix(count)),trim(filename),direction(sigma),trim(dcfield(dcfield_dependence)),i,trim(energypart),s%nkpt,eta,Utype,trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err.ne.0) missing_files = trim(missing_files) // " " // trim(varm)
      end do ; end do
      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_dc_beff_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
    else
      do sigma=1,4 ; do i=1,s%nAtoms
        iw = 80000+(sigma-1)*s%nAtoms+i
        close(unit=iw)
      end do ; end do
    end if

    return
  end subroutine openclose_dc_beff_files

  ! This subroutine write all the effective fields into files
  ! (already opened with openclose_dc_beff_files(1))
  ! Some information may also be written on the screen
  subroutine write_dc_beff()
    !use mod_f90_kind
    use mod_parameters, only: sigmai2i
    use mod_magnet, only: mvec_spherical, dc_fields, hw_count
    use mod_system, only: s => sys
    implicit none
    integer      :: i,iw,sigma
    real(double) :: phase,sine,cosine

    do sigma=1,4 ; do i=1,s%nAtoms
      iw = 80000+(sigma-1)*s%nAtoms+i

      if(abs(Beff_cart(sigmai2i(sigma,i)))>=1.d-10) then
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
    !use mod_f90_kind
    use mod_parameters, only: Npl,itype
    use mod_tools, only: sort_file
    implicit none
    integer :: i,sigma,iw,idc=1

    ! Opening effective field files
    if(itype==9) then
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
    if(itype==9) then
      call openclose_dc_beff_files(2)
    else
      call openclose_beff_files(2)
    end if

    return
  end subroutine sort_beff

end module mod_beff
