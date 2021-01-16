! This subroutine may be used to do tests on the program
subroutine debugging()
  use, intrinsic :: iso_fortran_env
  use mod_constants
  use mod_parameters
  use mod_magnet
  use mod_prefactors
  use mod_generate_epoints
  use mod_mpi_pars
  use mod_progress
  use mod_diamagnetic_current
  implicit none

  if(myrank==0) write(outputunit,"('[debugging] Starting to debug...')")

!   integer   :: i,i0,i1,j,j0,j1
!   real(dp)  :: kp(3)
!   complex(dp),dimension((Npl+2)*18,(Npl+2)*18)  :: hk,test,test2
!   complex(dp),dimension(Npl+2,Npl+2,18,18)  :: gf

!   write(*,*) 'Hamiltonian:'
  mz    = 2._dp
  mp  = cZero
  ! Variables used in the hamiltonian
  eps1 = 0._dp
!   lambda=0._dp
  ! call calculate_idia()



!   Finalizing program
  if(myrank==0) call write_time('[main] Finished on: ',outputunit)
  call MPI_Finalize(ierr)
  if ((ierr/=0).and.(myrank==0)) write(outputunit,"('[main] Something went wrong in the parallelization! ierr = ',i0)") ierr
  stop

end subroutine debugging
