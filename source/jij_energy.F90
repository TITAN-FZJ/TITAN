! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
subroutine jij_energy(Jij)
  use mod_f90_kind, only: double
  use mod_constants, only: pi
  use mod_parameters, only: lverbose, host, outputunit, nmaglayers, Ef, pn1
  use mod_generate_epoints, only: y, wght
  use mod_mpi_pars
  use mod_progress
  implicit none
  integer         :: i
  real(double),dimension(nmaglayers,nmaglayers,3,3) :: Jijint
  real(double),dimension(nmaglayers,nmaglayers,3,3),intent(out) :: Jij
!--------------------- begin MPI vars --------------------
  integer :: ix,itask,ncount
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  ncount=nmaglayers*nmaglayers*3*3


  ix = myrank+1
  itask = numprocs ! Number of tasks done initially

  ! Calculating the number of particles for each spin and orbital using a complex integral
  Jij = 0.d0
  do while( ix <= pn1)
    call sumk_jij(Ef,y(ix),Jijint)
    Jij = Jij + Jijint * wght(ix)
    ix = ix + numprocs
  end do

  call MPI_Allreduce(MPI_IN_PLACE, Jij, ncount, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  Jij = Jij / pi
  return
  end subroutine jij_energy
