subroutine lattice_definitions()
  use mod_f90_kind, only: double
  use mod_constants, only: pi
  use mod_mpi_pars
  use mod_system, only: l_nn, r_nn, c_nn, pln_a1, pln_a2, pln_normal, pln_cnt
  use mod_parameters
  use mod_tools, only: cross_unit, is_parallel
  implicit none
  integer :: i
  real(double),dimension(3) :: versor_Eperp

  if(dirEfield == -1) then
    dirEfieldvec = dirEfieldvec / sqrt(dot_product(dirEfieldvec, dirEfieldvec))
    EFp = atan(dirEfieldvec(2) / dirEfieldvec(1))
    EFt = acos(dirEfieldvec(3))
  else if(dirEfield == -2) then
    dirEfieldvec = dirEfieldvec(1) * pln_a1 + dirEfieldvec(2) * pln_a2
    dirEfieldvec = dirEfieldvec / sqrt(dot_product(dirEfieldvec, dirEfieldvec))
    EFp = atan(dirEfieldvec(2) / dirEfieldvec(1))
    EFt = acos(dirEfieldvec(3))
  else if(dirEfield == -3) then
    dirEfieldvec = [cos(EFp*pi/180)*sin(EFt*pi/180), sin(EFp*pi/180)*sin(EFt*pi/180), cos(EFt*pi/180)]
  else if(dirEfield >=1 .and. dirEfield <= pln_cnt(1)) then
    dirEfieldvec = c_nn(:, dirEfield)
  else
    if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Unknown E-Field direction!')")
    call MPI_Finalize(ierr)
    stop
  end if

  ! Calculating out-of-plane unit vector

  ! Calculating vector in-plane perpendicular to electric field direction
  if(0.d0 == dot_product(pln_normal, pln_normal)) then
    versor_Eperp = 0.d0
    do i = l_nn(1,1), l_nn(1,2) - 1
      if(.not. is_parallel(r_nn(:,i), dirEfieldvec)) then
        versor_Eperp = cross_unit(r_nn(:,i), dirEfieldvec)
        exit
      end if
    end do
    if(0.d0 == dot_product(versor_Eperp, versor_Eperp)) then
      !TODO: proper error message
      if(0 == myrank) write(outputunit,*) "[lattice_defintions] Error dirEfieldvec"
      call MPI_Finalize(ierr)
    end if
  else
    versor_Eperp = cross_unit(pln_normal, dirEfieldvec)
  end if


  ! Storing neighbors and angles for longitudinal and transverse currents
  sha_longitudinal = 0
  sha_transverse   = 0
  longitudinal_neighbors = 0
  transverse_neighbors = 0
  do i=l_nn(1,1), l_nn(1,2)-1
    long_cos(i) = dot_product(c_nn(:,i), dirEfieldvec)
    if(long_cos(i)>1.d-8) then
      longitudinal_neighbors = longitudinal_neighbors+1
      sha_longitudinal(longitudinal_neighbors) = i
    end if
    transv_cos(i) = dot_product(c_nn(:,i),versor_Eperp)
    if(transv_cos(i)>1.d-8) then
      transverse_neighbors = transverse_neighbors+1
      sha_transverse(transverse_neighbors) = i
    end if
  end do
  if(0 == myrank) then
    write(outputunit,"('[lattice_definitions] Longitudinal neighbors: ',10(i0,2x))") (sha_longitudinal(i),i=1,longitudinal_neighbors)
    write(outputunit,"('[lattice_definitions]   Transverse neighbors: ',10(i0,2x))") (sha_transverse(i),i=1,transverse_neighbors)
  end if

  return
end subroutine lattice_definitions
