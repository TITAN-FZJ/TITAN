!   Calculates band structure
subroutine band_structure(s)
  use mod_f90_kind,          only: double
  use mod_parameters,        only: output, nQvec, nQvec1, kpoints, bsfile, deltak
  use mod_system,            only: System
  use mod_tools,             only: cross
  use mod_mpi_pars,          only: rField
  use mod_superconductivity, only: superCond, hamiltk_sc, lsuperCond
  implicit none
  type(System), intent(in) :: s
  integer :: j, info, count
  integer :: lwork,dimbs
  real(double), dimension(:), allocatable :: rwork,eval
  complex(double), allocatable :: work(:),hk(:,:)

  dimbs = (s%nAtoms)*18*superCond
  lwork = 2*dimbs-1
  allocate( hk(dimbs,dimbs),rwork(3*dimbs-2),eval(dimbs),work(lwork) )

  if(rField == 0) write(output%unit_loop,"('CALCULATING THE BAND STRUCTURE')")

  open (unit=666, file=bsfile, status='old')
  do count=1,nQvec1
    write(output%unit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") count,nQvec1,dble((count-1.d0)/nQvec)

    if(lsuperCond) then
        call hamiltk_sc(s,kpoints(:,count),hk)
    else
        call hamiltk(s,kpoints(:,count),hk)
    end if

    call zheev('N','U',dimbs,hk,dimbs,eval,work,lwork,rwork,info)

    if(info/=0) then
      write(output%unit_loop,"('[band_structure] Problem with diagonalization. info = ',i0)") info
      stop
    end if
    ! Transform energy to eV if runoption is on
    ! eval = eval
    write(unit=666,fmt='(1000(es16.8))') dble((count-1.d0)*deltak), (eval(j),j=1,dimbs)
  end do
  close(666)
  deallocate(hk,rwork,eval,work)

end subroutine band_structure
