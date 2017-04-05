subroutine lattice_definitions()
  use mod_f90_kind
  use mod_constants
  use mod_mpi_pars
  use mod_lattice
  use mod_parameters
  use mod_tools, only: cross_unit
  integer :: i,j,err
  real(double),dimension(3) :: vec1,vec2

  lattice_dependent_definitions: select case (lattice)
  case("bcc110")
!----------- Direction of applied in-plane electric field --------------
    read(dirEfield,fmt=*,iostat=err) i
    if(err==0) then
      direction_E_field_bcc110_neighbors: select case (i)
      case (1:6)   !    In plane neighbors:
        dirEfieldvec = c0(i,:)
      case default !    Other direction:
        dirEfieldvec = [1.d0/(sq2*sq3),1.d0/(sq2*sq3),-sq2/sq3] ! In-plane, Perpendicular to the 1st n.n.
      end select direction_E_field_bcc110_neighbors
    else
      direction_E_field_bcc110_axis: select case (dirEfield)
      case ("L")
    !   In-plane long axis:
        dirEfieldvec = [hsq2,hsq2,0.d0]
      case ("S")
    !   In-plane short axis
        dirEfieldvec = [0.d0,0.d0,1.d0]
      case ("O")
    !   Other direction:
    !     dirEfieldvec = [1.d0/sq3,1.d0/sq3,1.d0/sq3] ! In the direction of the 1st n.n.
        dirEfieldvec = [1.d0/(sq2*sq3),1.d0/(sq2*sq3),-sq2/sq3] ! In-plane, Perpendicular to the 1st n.n.
      end select direction_E_field_bcc110_axis
    end if
  case("fcc100")
!----------- Direction of applied in-plane electric field --------------
    read(unit=dirEfield,fmt=*) i
    direction_E_field_fcc100: select case (i)
    case (1:8)   !    In plane neighbors:
      dirEfieldvec = c0(i,:)
    case default !    Other direction:
      dirEfieldvec = [1.d0,0.d0,0.d0] ! In-plane
    end select direction_E_field_fcc100
  case("fcc111")
!----------- Direction of applied in-plane electric field --------------
    read(unit=dirEfield,fmt=*) j
    direction_E_field_fcc111: select case (j)
!       case (1:8)   !    In plane neighbors:
!         dirEfieldvec = c0(j,:)
!       case default !    Other direction:
!         dirEfieldvec = [1.d0,0.d0,0.d0] ! In-plane
    end select direction_E_field_fcc111
  end select lattice_dependent_definitions

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