! This subroutine calculates LDOS
subroutine ldos()
  use mod_f90_kind,      only: double
  use mod_parameters,    only: nOrb, output, nEner1, emin, deltae,laddresults
  use mod_system,        only: s => sys
  use mod_BrillouinZone, only: realBZ
  use mod_LDOS
  use mod_mpi_pars
  implicit none
  integer :: i, j
  real(double) :: e

  call allocateLDOS()
  call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))

  ! Opening files
  if(rField == 0) then
    write(output%unit_loop,"('CALCULATING LDOS')")
    if(.not.laddresults) call createLDOSFiles()
    call openLDOSFiles()
  end if

  do i = startFreq, endFreq
     e = emin + (i-1)*deltae
     if(rFreq(1) == 0) write(output%unit_loop,"('[ldos] ',i0,' of ',i0,' points',', e = ',es10.3)") i,nEner1,e
     call ldos_energy(e,ldosu,ldosd)

     if(rFreq(1) == 0) then
        if(rFreq(2) == 0) then
           do j = 1, sFreq(2)
             if (j /= 1) then
               call MPI_Recv(e,     1            ,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE  ,1000,FreqComm(2),stat,ierr)
               call MPI_Recv(ldosd, s%nAtoms*nOrb,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1100,FreqComm(2),stat,ierr)
               call MPI_Recv(ldosu, s%nAtoms*nOrb,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1200,FreqComm(2),stat,ierr)
               end if

               ! Writing into files
               call writeLDOS(e)
            end do
         else
            call MPI_Send(e,     1            ,MPI_DOUBLE_PRECISION,0,1000,FreqComm(2),stat,ierr)
            call MPI_Send(ldosd, s%nAtoms*nOrb,MPI_DOUBLE_PRECISION,0,1100,FreqComm(2),stat,ierr)
            call MPI_Send(ldosu, s%nAtoms*nOrb,MPI_DOUBLE_PRECISION,0,1200,FreqComm(2),stat,ierr)
         end if
      end if
      call MPI_Barrier(FieldComm, ierr)
   end do

  call deallocateLDOS()

  ! Closing files
  if(rField == 0) call closeLDOSFiles()

end subroutine ldos
