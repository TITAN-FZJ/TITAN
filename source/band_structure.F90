!   Calculates magnetic LDOS
subroutine band_structure()
  use mod_f90_kind
  use mod_constants, only: pi,sq2,tpi
  use mod_parameters
  use mod_mpi_pars, only: mpitag
  implicit none
  character(len=400) :: varm
  character(len=50)  :: fieldpart,socpart
  character(len=1)   :: SOCc
  integer            :: j,ifail
  integer                       :: lwork,dimbs
  real(double)                  :: kmin(3),kmax(3),deltak(3)
  real(double),allocatable      :: rwork(:),kpoints(:,:)
  complex(double),allocatable   :: eval(:),evecl(:,:),evecr(:,:),work(:)
  complex(double),allocatable   :: hk(:,:)

  dimbs = (Npl+2)*18
  lwork = 33*dimbs
  allocate( hk(dimbs,dimbs),rwork(2*dimbs),eval(dimbs),evecl(1,dimbs),evecr(1,dimbs),work(lwork) )

  write(outputunit_loop,"('CALCULATING THE BAND STRUCTURE')")

  select case (lattice)
  case("bcc110")
    band_struct_direction_bcc110: select case (kdirection)
    case ("GH")
      kmin = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
      kmax = [ 0.d0 , 0.d0 , 1.5d0*pi/a0 ]      ! H
    case ("HP")
      kmin = [ 0.d0 , 0.d0 , 1.5d0*pi/a0 ]      ! H
      kmax = [ pi*sq2/a0 , 0.d0 , 0.5d0*pi/a0 ] ! P
    case ("PN")
      kmin = [ pi*sq2/a0 , 0.d0 , 0.5d0*pi/a0 ] ! P
      kmax = [ pi*sq2/a0 , 0.d0 , 0.d0 ]        ! N
    case ("NG")
      kmin = [ pi*sq2/a0 , 0.d0 , 0.d0 ]        ! N
      kmax = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
    case ("HG")
      kmin = [ 0.d0 , 0.d0 , 1.5d0*pi/a0 ]      ! H
      kmax = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
    case ("PH")
      kmin = [ pi*sq2/a0 , 0.d0 , 0.5d0*pi/a0 ] ! P
      kmax = [ 0.d0 , 0.d0 , 1.5d0*pi/a0 ]      ! H
    case ("NP")
      kmin = [ pi*sq2/a0 , 0.d0 , 0.d0 ]        ! N
      kmax = [ pi*sq2/a0 , 0.d0 , 0.5d0*pi/a0 ] ! P
    case ("GN")
      kmin = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
      kmax = [ pi*sq2/a0 , 0.d0 , 0.d0 ]        ! N
    case default
      write(outputunit_loop,"('[band_structure] Choose one of the following directions in k space:')")
      write(outputunit_loop,"('[band_structure] GH, HP, PN, NG (or the opposite HG, PH, NP, GN)')")
      stop
    end select band_struct_direction_bcc110
    ! Transformation to cartesian axis
    kmin(1) = kmin(1)*0.5d0*sq2
    kmin(2) = kmin(1)
    kmin(3) = kmin(3)

    kmax(1) = kmax(1)*0.5d0*sq2
    kmax(2) = kmax(1)
    kmax(3) = kmax(3)
  case("fcc100")
    band_struct_direction_fcc100: select case (kdirection)
    case ("GM")
      kmin = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
      kmax = [ 1.d0 , 0.d0 , 0.d0 ]*tpi/a0      ! M
    case ("MX")
      kmin = [ 1.d0 , 0.d0 , 0.d0 ]*tpi/a0      ! M
      kmax = [ 1.d0 , 1.d0 , 0.d0 ]*pi/a0       ! X
    case ("XG")
      kmin = [ 1.d0 , 1.d0 , 0.d0 ]*pi/a0       ! X
      kmax = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
    case ("GX")
      kmin = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
      kmax = [ 1.d0 , 1.d0 , 0.d0 ]*pi/a0       ! X
    case ("XM")
      kmin = [ 1.d0 , 1.d0 , 0.d0 ]*pi/a0       ! X
      kmax = [ 1.d0 , 0.d0 , 0.d0 ]*tpi/a0      ! M
    case ("MG")
      kmin = [ 1.d0 , 0.d0 , 0.d0 ]*tpi/a0      ! M
      kmax = [ 0.d0 , 0.d0 , 0.d0 ]             ! Gamma
    case default
      write(outputunit_loop,"('[band_structure] Choose one of the following directions in k space:')")
      write(outputunit_loop,"('[band_structure] GM, MX, XG (or the opposite GX, XM, MG)')")
      stop
    end select band_struct_direction_fcc100

    ! Transformation to cartesian axis
    kmin(3) = kmin(1)+kmin(2) ! using kmin(3) as a temporary variable
    kmin(2) =-kmin(1)+kmin(2)
    kmin(1) = kmin(3)
    kmin(3) = 0.d0

    kmax(3) = kmax(1)+kmax(2) ! using kmax(3) as a temporary variable
    kmax(2) =-kmax(1)+kmax(2)
    kmax(1) = kmax(3)
    kmax(3) = 0.d0

  end select

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

  write(varm,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a2,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),kdirection,nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=666+mpitag, file=varm,status='unknown')

  deltak = (kmax - kmin)/npts
  allocate( kpoints(npt1,3) )
  do count=1,npt1
    kpoints(count,:) = kmin + (count-1)*deltak
  end do

  band_structure_loop: do count=1,npt1
    write(outputunit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") count,npt1,dble((count-1.d0)/npts)
    call hamiltk(kpoints(count,:),hk)

    call zgeev('N','N',dimbs,hk,dimbs,eval,evecl,1,evecr,1,work,lwork,rwork,ifail)
    if(ifail.ne.0) then
      write(outputunit_loop,"('[band_structure] Problem with diagonalization. ifail = ',i0)") ifail
      stop
!       else
!         write(outputunit_loop,"('[band_structure] optimal lwork = ',i0,' lwork = ',i0)") work(1),lwork
    end if
    ! Transform energy to eV if runoption is on
    eval = eval
    write(666+mpitag,'(1000(es16.8))') dble((count-1.d0)/npts),(real(eval(j)),j=1,dimbs)
  end do band_structure_loop
  close(666+mpitag)
  deallocate(hk,rwork,eval,evecl,evecr,work)

  return
end subroutine band_structure
