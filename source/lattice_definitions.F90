subroutine lattice_definitions()
  use mod_f90_kind
  use mod_constants
  use mod_mpi_pars
  use mod_lattice
  use mod_parameters
  use mod_tools, only: cross_unit
  integer :: i,j,err
  real(double),dimension(3) :: vec1,vec2



  if(dirEfield == -1) then
    dirEfieldvec = dirEfieldvec / sqrt(dot_product(dirEfieldvec, dirEfieldvec))
    EFp = atan(dirEfieldvec(2) / dirEfieldvec(1))
    EFt = acos(dirEfieldvec(3))
  else if(dirEfield == -2) then
    dirEfieldvec = dirEfieldvec(1) * a1_pln + dirEfieldvec(2) * a2_pln
    dirEfieldvec = dirEfieldvec / sqrt(dot_product(dirEfieldvec, dirEfieldvec))
    EFp = atan(dirEfieldvec(2) / dirEfieldvec(1))
    EFt = acos(dirEfieldvec(3))
  else if(dirEfield == -3) then 
    dirEfieldvec = [cos(EFp)*sin(EFt), sin(EFp)*sin(EFt), cos(EFt)]
  else if(dirEfield >=1 .and. dirEfield <= n0) then
    dirEfieldvec = c0(dirEfield,:)
  else
    if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Unknown E-Field direction!')")
    call MPI_Finalize(ierr)
    stop
  end if

  ! Calculating out-of-plane unit vector
  vec1 = r0(1,:)
  vec2 = r0(2,:)
  versor_oop   = cross_unit(vec1,vec2)
  ! Calculating vector in-plane perpendicular to electric field direction
  versor_Eperp = cross_unit(versor_oop,dirEfieldvec)

  ! Storing neighbors and angles for longitudinal and transverse currents
  sha_longitudinal = 0
  sha_transverse   = 0
  longitudinal_neighbors = 0
  transverse_neighbors = 0
  do i=1,n0
    long_cos(i) = dot_product(c0(i,:),dirEfieldvec)
    if(long_cos(i)>1.d-8) then
      longitudinal_neighbors = longitudinal_neighbors+1
      sha_longitudinal(longitudinal_neighbors) = i
    end if
    transv_cos(i) = dot_product(c0(i,:),versor_Eperp)
    if(transv_cos(i)>1.d-8) then
      transverse_neighbors = transverse_neighbors+1
      sha_transverse(transverse_neighbors) = i
    end if
  end do
  if(myrank==0) then
    write(outputunit,"('[lattice_definitions] Longitudinal neighbors: ',10(i0,2x))") (sha_longitudinal(i),i=1,longitudinal_neighbors)
    write(outputunit,"('[lattice_definitions]   Transverse neighbors: ',10(i0,2x))") (sha_transverse(i),i=1,transverse_neighbors)
  end if

  return
end subroutine lattice_definitions
