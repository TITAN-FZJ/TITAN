!   Calculates band structure
subroutine band_structure(s)
  use mod_kind,              only: dp
  use mod_parameters,        only: output, kdirection, nQvec, nQvec1, kpoints, bsfile, wsfile, deltak
  use mod_system,            only: System_type
  use mod_tools,             only: cross,diagonalize
  use mod_mpi_pars,          only: rField
  use mod_superconductivity, only: superCond
  use mod_io,                only: write_header
  use mod_hamiltonian,       only: hamilt_local,hamiltk
  implicit none
  type(System_type), intent(in) :: s
  integer :: i, kount, f_unit=666, w_unit=667, n
  integer :: dimbs
  real(dp), dimension(:), allocatable :: eval
  complex(dp), allocatable :: hk(:,:)
  character(len=30) :: formatvar1,formatvar2

  external :: zheev
  
  if(rField == 0) write(output%unit_loop,"('CALCULATING THE BAND STRUCTURE')")

  dimbs = (s%nAtoms)*18*superCond
  allocate( hk(dimbs,dimbs),eval(dimbs) )

  call hamilt_local(s)

  ! Band structure file
  write(bsfile,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
  open (unit=f_unit, file=bsfile, status='replace')
  call write_header(f_unit,"# dble((kount-1._dp)*deltak), (eval(i),i=1,dimbs)",s%Ef)

  write(wsfile,"('./results/',a1,'SOC/',a,'/BS/weights_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
  open (unit=w_unit, file=wsfile, status='replace')
  call write_header(w_unit,"# ((abs(hk(i,n))**2,i=1,dimbs),n=1,dimbs)",s%Ef)

  write(unit=f_unit, fmt="(a,2x,i3)") "# dimbs ",dimbs
  write(formatvar1,fmt="(a,i0,a)") '(',1+dimbs,'(es16.8e3,2x))'
  write(formatvar2,fmt="(a,i0,a)") '(',dimbs*dimbs,'(es10.3e3,2x))'

  do kount=1,nQvec1
    write(output%unit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") kount,nQvec1,dble((kount-1._dp)/nQvec)

    call hamiltk(s,kpoints(:,kount),hk)
  
    call diagonalize(dimbs,hk,eval)

    write(unit=f_unit,fmt=formatvar1) dble((kount-1._dp)*deltak), (eval(i),i=1,dimbs)
    write(unit=w_unit,fmt=formatvar2) ((abs(hk(i,n))**2,i=1,dimbs),n=1,dimbs)
  end do

  close(f_unit)
  close(w_unit)
  deallocate(hk,eval)

end subroutine band_structure
