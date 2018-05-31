!   Calculates band structure
subroutine band_structure(s)
  use mod_f90_kind,      only: double
  use mod_parameters,    only: output, nQvec, nQvec1, kpoints, bsfile, deltak
  use mod_system,        only: System
  use TightBinding,      only: nOrb
  use mod_tools,         only: cross
  use mod_mpi_pars,      only: rField
  implicit none
  type(System), intent(in) :: s
  integer :: j, ifail, count
  integer :: lwork,dimbs
  real(double), allocatable :: rwork(:)
  complex(double), allocatable :: eval(:),evecl(:,:),evecr(:,:),work(:)
  complex(double), allocatable :: hk(:,:)

  dimbs = (s%nAtoms)*18
  lwork = 33*dimbs
  allocate( hk(dimbs,dimbs),rwork(2*dimbs),eval(dimbs),evecl(1,dimbs),evecr(1,dimbs),work(lwork) )

  if(rField == 0) write(output%unit_loop,"('CALCULATING THE BAND STRUCTURE')")

  open (unit=666, file=bsfile, status='old')
  do count=1,nQvec1
    write(output%unit_loop,"('[band_structure] ',i0,' of ',i0,' points',', i = ',es10.3)") count,nQvec1,dble((count-1.d0)/nQvec)
    call hamiltk(s,kpoints(:,count),hk)

    call zgeev('N','N',dimbs,hk,dimbs,eval,evecl,1,evecr,1,work,lwork,rwork,ifail)
    if(ifail/=0) then
      write(output%unit_loop,"('[band_structure] Problem with diagonalization. ifail = ',i0)") ifail
      stop
    end if
    ! Transform energy to eV if runoption is on
    ! eval = eval
    write(unit=666,fmt='(1000(es16.8))') dble((count-1.d0)*deltak), (real(eval(j)),j=1,dimbs)
  end do
  close(666)
  deallocate(hk,rwork,eval,evecl,evecr,work)

end subroutine band_structure
