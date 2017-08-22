! This subroutine calculates LDOS
subroutine ldos()
  use mod_f90_kind, only: double
  use mod_parameters, only: outputunit_loop, fieldpart, Npl_folder, eta, Utype, npt1, emin, deltae
  use mod_SOC, only: SOCc, socpart
  use mod_system, only: s => sys
  use mod_mpi_pars, only: mpitag
  use TightBinding, only: nOrb
  implicit none
  character(len=400) :: varm
  integer            :: i,iw, count
  real(double)       :: e
  real(double),dimension(:,:),allocatable :: ldosu,ldosd

  allocate(ldosu(s%nAtoms, nOrb),ldosd(s%nAtoms,nOrb))

  write(outputunit_loop,"('CALCULATING LDOS')")

  ! Opening files
  ! LDOS
  do i=1,s%nAtoms
    iw = 1000+(mpitag-1)*s%nAtoms*2 + (i-1)*2 + 1
    write(varm,"('./results/',a1,'SOC/',a,'/LDOS/ldosu_layer',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),i,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=iw, file=varm,status='unknown')
    iw = iw+1
    write(varm,"('./results/',a1,'SOC/',a,'/LDOS/ldosd_layer',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),i,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=iw, file=varm,status='unknown')
  end do

  ldos_energy_loop: do count=1,npt1
    e = emin + (count-1)*deltae
    write(outputunit_loop,"('[ldos] ',i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e

    call ldos_energy(e,ldosu,ldosd)

    ! Writing into files
    ! LDOS
    ldos_writing_plane_loop: do i=1,s%nAtoms
        iw = 1000+(mpitag-1)*s%nAtoms*2 + (i-1)*2 + 1
        write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
        iw = iw+1
        write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
    end do ldos_writing_plane_loop

  end do ldos_energy_loop

  ! Closing files
  do i=1,s%nAtoms
    iw = 1000+(mpitag-1)*s%nAtoms*2 + (i-1)*2 + 1
    close (iw)
    iw = iw+1
    close (iw)
  end do

  deallocate(ldosu,ldosd)

  return
end subroutine ldos
