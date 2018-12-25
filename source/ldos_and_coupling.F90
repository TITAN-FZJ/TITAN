! This subroutine calculates LDOS and coupling as a function of energy
subroutine ldos_and_coupling()
  use mod_f90_kind, only: double
  use mod_parameters, only: output, emin, deltae, nEner1, skip_steps
  use mod_system, only: s => sys
  use mod_LDOS
  use mod_Coupling
  use TightBinding, only: nOrb
  use mod_BrillouinZone, only: realBZ
  use mod_mpi_pars
  implicit none
  integer :: i, j, mu, count, ncount,ncount2,ncount3
  real(double) :: e

  ncount  = s%nAtoms*nOrb
  ncount2 = s%nAtoms*s%nAtoms
  ncount3 = ncount*9

  if(rField == 0) write(output%unit_loop,"('CALCULATING LDOS AND EXCHANGE INTERACTIONS AS A FUNCTION OF ENERGY')")

  call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))

  ! Opening files
  if(rField == 0) then
     ! LDOS
     call openLDOSFiles()
     ! Exchange interactions
     call openCouplingFiles()
  end if

  call allocateLDOS()
  call allocateCoupling()

  do count = startFreq + skip_steps, endFreq + skip_steps
    e = emin + (count-1) * deltae
    if(rFreq(1) == 0) write(output%unit_loop,"('[ldos_and_coupling] ',i0,' of ',i0,' points',', e = ',es10.3)") count,nEner1,e

    call ldos_jij_energy(e,ldosu,ldosd,Jij)

    if(rFreq(1) == 0) then
       do i = 1, s%nAtoms
         do j = 1, s%nAtoms
            trJij(i,j)    = 0.5d0*(Jij(i,j,1,1) + Jij(i,j,2,2))
            Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
            Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
            do mu = 1, 3
              Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
            end do
         end do
       end do


       if(rFreq(2) == 0) then
         do i = 1, sFreq(2)
            if (i /= 1) then
               call MPI_Recv(e,     1       ,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE  ,1000,FreqComm(2),stat,ierr)
               call MPI_Recv(ldosd, ncount  ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1100,FreqComm(2),stat,ierr)
               call MPI_Recv(ldosu, ncount  ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1200,FreqComm(2),stat,ierr)
               call MPI_Recv(trJij, ncount2 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1300,FreqComm(2),stat,ierr)
               call MPI_Recv(Jij,   ncount3 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1400,FreqComm(2),stat,ierr)
               call MPI_Recv(Jijs,  ncount3 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1500,FreqComm(2),stat,ierr)
               call MPI_Recv(Jija,  ncount3 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1600,FreqComm(2),stat,ierr)
            end if

            ! Writing into files
            call writeLDOS(e)
            ! Exchange interactions
            call writeCoupling(e)
         end do
       else
          call MPI_Recv(e,     1       ,MPI_DOUBLE_PRECISION,0,1000,FreqComm(2),stat,ierr)
          call MPI_Recv(ldosd, ncount  ,MPI_DOUBLE_PRECISION,0,1100,FreqComm(2),stat,ierr)
          call MPI_Recv(ldosu, ncount  ,MPI_DOUBLE_PRECISION,0,1200,FreqComm(2),stat,ierr)
          call MPI_Recv(trJij, ncount2 ,MPI_DOUBLE_PRECISION,0,1300,FreqComm(2),stat,ierr)
          call MPI_Recv(Jij,   ncount3 ,MPI_DOUBLE_PRECISION,0,1400,FreqComm(2),stat,ierr)
          call MPI_Recv(Jijs,  ncount3 ,MPI_DOUBLE_PRECISION,0,1500,FreqComm(2),stat,ierr)
          call MPI_Recv(Jija,  ncount3 ,MPI_DOUBLE_PRECISION,0,1600,FreqComm(2),stat,ierr)
       end if
    end if
    call MPI_Barrier(FieldComm, ierr)

  end do

  call deallocateLDOS()
  call deallocateCoupling()

  ! Closing files
  if(rField == 0) then
     call closeLDOSFiles()
     call closeCouplingFiles()
  end if

end subroutine ldos_and_coupling
