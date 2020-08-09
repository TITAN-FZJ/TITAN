!   Calculates band structure
subroutine band_structure(s)

  use mod_f90_kind,          only: double
  use mod_parameters,        only: output, kdirection, nQvec, nQvec1, kpoints, bsfile, deltak
  use mod_system,            only: System
  use mod_tools,             only: cross
  use mod_mpi_pars,          only: rField
  use mod_superconductivity, only: superCond, hamiltk_sc, lsuperCond
  use mod_io,                only: write_header
  use mod_hamiltonian,       only: hamiltk
  implicit none
  type(System), intent(in) :: s
  integer :: i, info, count, f_unit=666, n
  integer :: lwork,dimbs
  real(double), dimension(:), allocatable :: rwork,eval
  complex(double), allocatable :: work(:),hk(:,:)
  character(len=30) :: formatvar

  dimbs = (s%nAtoms)*18*superCond
  lwork = 2*dimbs-1
  allocate( hk(dimbs,dimbs),rwork(3*dimbs-2),eval(dimbs),work(lwork) )

  if(rField == 0) write(output%unit_loop,"('CALCULATING THE BAND STRUCTURE')")

  write(bsfile,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
  open (unit=f_unit, file=bsfile, status='replace')
  call write_header(f_unit,"# dble((count-1.d0)*deltak), (eval(i),i=1,dimbs)",s%Ef)

  write(unit=f_unit, fmt="(a,2x,i3)") "# dimbs ",dimbs
  write(formatvar,fmt="(a,i0,a)") '(',1+dimbs*(dimbs+1),'(es16.8e3,2x))'


  do count=1,nQvec1
    write(output%unit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") count,nQvec1,dble((count-1.d0)/nQvec)

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

    write(unit=f_unit,fmt=formatvar) dble((count-1.d0)*deltak), (eval(i),i=1,dimbs), ((abs(hk(i,n))**2,i=1,dimbs),n=1,dimbs)
  end do

  close(f_unit)
  deallocate(hk,rwork,eval,work)

end subroutine band_structure
