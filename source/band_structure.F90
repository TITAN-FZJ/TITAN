subroutine read_band_points(kbands, b1, b2, b3)
  use mod_f90_kind, only: double
  use mod_parameters, only: bands, band_cnt
  implicit none

  real(double), dimension(:,:), allocatable, intent(out) :: kbands
  real(double), dimension(3), intent(in) :: b1, b2, b3
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
      if(trim(kband(j)%name) == trim(bands(i))) then
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


!   Calculates magnetic LDOS
subroutine band_structure(s)
  use mod_f90_kind, only: double
  use mod_parameters, only: fieldpart, outputunit_loop, Ef, eta, Npl_folder, npts, npt1, Utype, bands, band_cnt
  use mod_SOC, only: SOCc, socpart
  use mod_mpi_pars, only: mpitag
  use mod_system, only: System
  implicit none
  type(System), intent(in) :: s
  character(len=400) :: varm
  integer :: i, j, ifail, count
  integer :: lwork,dimbs
  real(double) :: dir(3), deltak
  real(double), allocatable :: rwork(:),kpoints(:,:)
  complex(double), allocatable :: eval(:),evecl(:,:),evecr(:,:),work(:)
  complex(double), allocatable :: hk(:,:)
  character(len=20) :: kdirection
  real(double), dimension(:,:), allocatable :: band_points
  real(double) :: total_length

    interface
      subroutine read_band_points(kbands, b1, b2, b3)
        use mod_f90_kind, only: double
        real(double), dimension(:,:), allocatable, intent(out) :: kbands
        real(double), dimension(3), intent(in) :: b1, b2, b3
      end subroutine
    end interface

  dimbs = (s%nAtoms)*18
  lwork = 33*dimbs
  allocate( hk(dimbs,dimbs),rwork(2*dimbs),eval(dimbs),evecl(1,dimbs),evecr(1,dimbs),work(lwork) )

  write(outputunit_loop,"('CALCULATING THE BAND STRUCTURE')")

  call read_band_points(band_points, s%b1, s%b2, s%b3)
  total_length = 0.d0
  do i = 1, band_cnt - 1
    total_length = total_length + sqrt(dot_product(band_points(:,i)-band_points(:,i+1), band_points(:,i)-band_points(:,i+1)))
  end do

  write(kdirection,*) (trim(adjustl(bands(i))), i = 1,band_cnt)
  write(varm,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),trim(adjustl(kdirection)),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=666+mpitag, file=varm,status='replace')
  write(unit=666+mpitag, fmt=*) band_cnt, npts
  write(unit=666+mpitag, fmt=*) Ef
  deltak = total_length / npts
  allocate(kpoints(3,npt1))

  total_length = 0.d0
  count = 0
  i = 0
  do count = 1, npt1
    if(deltak*count >= total_length .and. i < band_cnt - 1) then
      i = i + 1
      write(unit=666+mpitag, fmt=*) bands(i), total_length
      total_length = total_length + sqrt(dot_product(band_points(1:3,i)-band_points(1:3,i+1), band_points(1:3,i)-band_points(1:3,i+1)))
      j = 0
      dir = (band_points(1:3,i+1) - band_points(1:3,i))
      dir = dir / sqrt(dot_product(dir,dir))
    end if
    kpoints(:,count) = band_points(:,i) + dir * j * deltak
    j = j + 1
  end do
  i = i+1
  write(unit=666+mpitag, fmt=*) bands(i), total_length

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
    ! else
    !   write(outputunit_loop,"('[band_structure] optimal lwork = ',i0,' lwork = ',i0)") work(1),lwork
    end if
    ! Transform energy to eV if runoption is on
    eval = eval
    write(666+mpitag,'(1000(es16.8))') dble((count-1.d0)*deltak), (real(eval(j)),j=1,dimbs)
  end do band_structure_loop
  close(666+mpitag)
  deallocate(hk,rwork,eval,evecl,evecr,work)

  return
end subroutine band_structure
