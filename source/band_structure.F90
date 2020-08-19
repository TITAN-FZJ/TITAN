!   Calculates band structure
subroutine band_structure(s)

  use mod_kind, only: dp
  use mod_parameters,        only: output, kdirection, nQvec, nQvec1, kpoints, bsfile, wsfile, deltak
  use mod_system,            only: System
  use mod_tools,             only: cross
  use mod_mpi_pars,          only: rField
  use mod_superconductivity, only: superCond, hamiltk_sc, lsuperCond
  use mod_io,                only: write_header
  use mod_hamiltonian,       only: hamiltk,hamilt_local,h0
  implicit none
  type(System), intent(in) :: s
  integer :: i, info, count, f_unit=666, w_unit=667, n
  integer :: lwork,dimbs
  real(dp), dimension(:), allocatable :: rwork,eval
  complex(dp), allocatable :: work(:),hk(:,:)
  character(len=30) :: formatvar1,formatvar2

  if(rField == 0) write(output%unit_loop,"('CALCULATING THE BAND STRUCTURE')")

  dimbs = (s%nAtoms)*18*superCond
  lwork = 2*dimbs-1
  allocate( hk(dimbs,dimbs),rwork(3*dimbs-2),eval(dimbs),work(lwork) )

  call hamilt_local(s)

  ! Band structure file
  write(bsfile,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
  open (unit=f_unit, file=bsfile, status='replace')
  call write_header(f_unit,"# dble((count-1._dp)*deltak), (eval(i),i=1,dimbs)",s%Ef)

  write(wsfile,"('./results/',a1,'SOC/',a,'/BS/weights_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
  open (unit=w_unit, file=wsfile, status='replace')
  call write_header(w_unit,"# ((abs(hk(i,n))**2,i=1,dimbs),n=1,dimbs)",s%Ef)

  write(unit=f_unit, fmt="(a,2x,i3)") "# dimbs ",dimbs
  write(formatvar1,fmt="(a,i0,a)") '(',1+dimbs,'(es16.8e3,2x))'
  write(formatvar2,fmt="(a,i0,a)") '(',dimbs*dimbs,'(es16.8e3,2x))'

  do count=1,nQvec1
    write(output%unit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") count,nQvec1,dble((count-1._dp)/nQvec)

    if(lsuperCond) then
        call hamiltk_sc(s,kpoints(:,count),hk)
    else
        call hamiltk(s,kpoints(:,count),hk)
    end if

    call zheev('V','U',dimbs,hk,dimbs,eval,work,lwork,rwork,info)

    if(info/=0) then
      write(output%unit_loop,"('[band_structure] Problem with diagonalization. info = ',i0)") info
      stop
    end if

    write(unit=f_unit,fmt=formatvar1) dble((count-1._dp)*deltak), (eval(i),i=1,dimbs)
    write(unit=w_unit,fmt=formatvar2) ((abs(hk(i,n))**2,i=1,dimbs),n=1,dimbs)
  end do

  close(f_unit)
  close(w_unit)
  deallocate(h0,hk,rwork,eval,work)

end subroutine band_structure
