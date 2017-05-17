!   Calculates magnetic LDOS
subroutine band_structure()
  use mod_f90_kind
  use mod_constants, only: pi,sq2,tpi
  use mod_parameters
  use mod_mpi_pars, only: mpitag
  use mod_system, only: a0, nkpt
  implicit none
  character(len=400) :: varm
  character(len=50)  :: fieldpart,socpart
  character(len=1)   :: SOCc
  integer            :: i, j,ifail
  integer                       :: lwork,dimbs
  real(double)                  :: kmin(3),kmax(3), dir(3), deltak
  real(double),allocatable      :: rwork(:),kpoints(:,:)
  complex(double),allocatable   :: eval(:),evecl(:,:),evecr(:,:),work(:)
  complex(double),allocatable   :: hk(:,:)
  real(double), dimension(:,:), allocatable :: band_points
  real(double) :: total_length

  dimbs = (Npl_total)*18
  lwork = 33*dimbs
  allocate( hk(dimbs,dimbs),rwork(2*dimbs),eval(dimbs),evecl(1,dimbs),evecr(1,dimbs),work(lwork) )

  write(outputunit_loop,"('CALCULATING THE BAND STRUCTURE')")

  call read_band_points(band_points)
  total_length = 0.d0
  do i = 1, band_cnt - 1
    total_length = total_length + sqrt(dot_product(band_points(:,i)-band_points(:,i+1), band_points(:,i)-band_points(:,i+1)))
  end do

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

  write(varm,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a2,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),kdirection,nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=666+mpitag, file=varm,status='replace')

  deltak = total_length / npts
  allocate(kpoints(3,npt1))

  total_length = 0.d0
  count = 0
  i = 0
  do count = 1, npt1
    if(deltak*count >= total_length) then
      i = i + 1
      total_length = total_length + sqrt(dot_product(band_points(:,i)-band_points(:,i+1), band_points(:,i)-band_points(:,i+1)))
      j = 0
      dir = (band_points(:,i+1) - band_points(:,i))
      dir = dir / sqrt(dot_product(dir,dir))
    end if
    kpoints(:,count) = band_points(:,i) + dir * j * deltak
    j = j + 1
    print *, kpoints(:,count)
  end do

  ! deltak = (kmax - kmin)/npts
  ! allocate( kpoints(npt1,3) )
  ! do count=1,npt1
  !   kpoints(count,:) = kmin + (count-1)*deltak
  ! end do

  band_structure_loop: do count=1,npt1
    write(outputunit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") count,npt1,dble((count-1.d0)/npts)
    call hamiltk(kpoints(:,count),hk)

    call zgeev('N','N',dimbs,hk,dimbs,eval,evecl,1,evecr,1,work,lwork,rwork,ifail)
    if(ifail/=0) then
      write(outputunit_loop,"('[band_structure] Problem with diagonalization. ifail = ',i0)") ifail
      stop
!       else
!         write(outputunit_loop,"('[band_structure] optimal lwork = ',i0,' lwork = ',i0)") work(1),lwork
    end if
    ! Transform energy to eV if runoption is on
    eval = eval
    write(666+mpitag,'(1000(es16.8))') dble((count-1.d0)*deltak), (real(eval(j)),j=1,dimbs)
  end do band_structure_loop
  close(666+mpitag)
  deallocate(hk,rwork,eval,evecl,evecr,work)

  return
end subroutine band_structure

subroutine read_band_points(kbands)
  use mod_parameters, only: bands, band_cnt
  use mod_system, only: b1, b2, b3
  use mod_f90_kind, only: double
  implicit none
  real(double), dimension(:,:), allocatable, intent(out) :: kbands
  character(len=40) :: band_file = "kbands"
  character(len=200) :: line
  integer :: ios, line_count, i, j
  logical :: found
  type :: band_point
    character(len=5) :: name
    real(double), dimension(3) :: kp
  end type
  type(band_point), dimension(:), allocatable :: kband

  line_count = 0
  open(unit = 666999, file = trim(band_file), status='old', iostat=ios)
  if(ios /= 0) then
    !TODO: Call upon hell
  endif

  ! Count non commented lines
  read(unit=666999, fmt='(A)', iostat=ios) line
  do while (ios == 0)
    if(len_trim(line) > 0 .and. len_trim(line) < 200 .and. index(adjustl(line), "#") /= 1) line_count = line_count + 1
    read(unit=666999, fmt='(A)', iostat=ios) line
  end do
  rewind(666999)
  print *, line_count

  allocate(kband(line_count))
  i = 0
  read(unit=666999, fmt='(A)', iostat=ios) line
  do while (ios == 0)
    if(len_trim(line) <= 0 .or. len_trim(line) >= 200 .or. index(adjustl(line),"#") == 1) cycle
    if(index(line, "#") > 0) line = line(1:index(line, "#")-1)
    i = i + 1
    read(unit = line, fmt = *, iostat=ios) kband(i)%name, (kband(i)%kp(j), j = 1,3)
    read(unit=666999, fmt='(A)', iostat=ios) line
  end do

  allocate(kbands(3,band_cnt))
  do i = 1, band_cnt
    found = .false.
    do j = 1, line_count
      if(trim(kband(j)%name) == bands(i)) then
        found = .true.
        kbands(1:3,i) = kband(j)%kp(1) * b1 + kband(j)%kp(2) * b2 + kband(j)%kp(3) * b3
        exit
      endif
    end do
    if(.not. found) then
      print *, "read_bands went wrong"
      stop
    endif
  end do

  close(666999)
end subroutine read_band_points
