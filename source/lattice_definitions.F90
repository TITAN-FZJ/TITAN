subroutine lattice_definitions()
  use mod_f90_kind
  use mod_constants
  use mod_mpi_pars
  use mod_parameters
  use mod_lattice
  integer :: i,j,err

  lattice_dependent_definitions: select case (lattice)
  case("bcc110")
!----------- Direction of applied in-plane electric field --------------
    read(dirEfield,fmt=*,iostat=err) i
    if(err.eq.0) then
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
!------------------- Spin quantization direction -----------------------
    magnetization_axis_bcc110: select case (magaxis)
    case ("L")
  !   In-plane long axis: (x - short axis; y - out-of-plane; z - long axis)
      phi   = pi/4.d0 ! around z (counter-clockwise from the top X->Y)
      theta = pi/2.d0 ! around y' (counter-clockwise from the top Z->X')
    case ("P")
  !   Out-of-plane: (x - short axis; y - long axis; z - out-of-plane)
      phi   = -pi/4.d0 ! around z (counter-clockwise from the top X->Y)
      theta = pi/2.d0 ! around y' (counter-clockwise from the top Z->X')
    case ("S")
  !   In-plane short axis: (x - long axis; y - out-of-plane; z - short axis)
      phi   = pi/4.d0 ! around z (counter-clockwise from the top X->Y)
      theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
    case default
      if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Choose a correct magnetization axis!')")
      call MPI_Finalize(ierr)
      stop
  !     phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
  !     theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
    end select magnetization_axis_bcc110
!-------------- Longitudinal and transverse DC voltage -----------------
    read(dirEfield,fmt=*,iostat=err) i
    if(err.eq.0) then
      vdc_definitions_bcc110_neighbors: select case (i)
      case (1:6) ! If the electric field is applied along a first n.n.
        lvdc = .false.
        vdcneighbor(1) = i ! Longitudinal neighbor
        vdcneighbor(2) = 0 ! No transverse neighbor
        ! There is a longitudinal neighbor, but it's not possible to map
        ! the transverse spin disturbance without rotation
        ! Magnetization components used in Vdc calculations
        ! x -> direction of the electric field
        ! y -> in-plane transverse to direction of the field
        ! z -> out-of-plane
        ! mapped in the magnetization reference frame
        select case (magaxis)
        case ("P")
          mvdcvector(1) = 2  ! x -> y (field)
          mvdcvector(2) = 1  ! y -> x
          mvdcvector(3) = 3  ! z -> z (out-of-plane)
        case ("L")
          mvdcvector(1) = 3  ! x -> z (field)
          mvdcvector(2) = 1  ! y -> x
          mvdcvector(3) = 2  ! z -> y (out-of-plane)
        case ("S")
          mvdcvector(1) = 1  ! x -> x (field)
          mvdcvector(2) = 3  ! y -> z
          mvdcvector(3) = 2  ! z -> y (out-of-plane)
        end select
      case default
    !   Other direction:
        if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Vdc not implemented for this electric field!')")
        lvdc = .false.
        vdcneighbor(1) = 0
        vdcneighbor(2) = 0
      end select vdc_definitions_bcc110_neighbors
    else
      vdc_definitions_bcc110_axis: select case (dirEfield)
      case ("L") ! If the electric field is applied along the long axis
        lvdc = .true.
        vdcneighbor(1) = 0 ! No longitudinal neighbor
        vdcneighbor(2) = 5 ! Transverse neighbor
        ! Magnetization components used in Vdc calculations
        ! x -> direction of the electric field
        ! y -> in-plane transverse to direction of the field
        ! z -> out-of-plane
        ! mapped in the magnetization reference frame
        select case (magaxis)
        case ("P")
          mvdcvector(1) = 2  ! x -> y (field)
          mvdcvector(2) = 1  ! y -> x
          mvdcvector(3) = 3  ! z -> z (out-of-plane)
        case ("L")
          mvdcvector(1) = 3  ! x -> z (field)
          mvdcvector(2) = 1  ! y -> x
          mvdcvector(3) = 2  ! z -> y (out-of-plane)
        case ("S")
          mvdcvector(1) = 1  ! x -> x (field)
          mvdcvector(2) = 3  ! y -> z
          mvdcvector(3) = 2  ! z -> y (out-of-plane)
        end select
      case ("S")
        lvdc = .true.
        vdcneighbor(1) = 5 ! Longitudinal neighbor
        vdcneighbor(2) = 0 ! No transverse neighbor
        select case (magaxis)
        case ("P")
          mvdcvector(1) = 1  ! x -> x (field)
          mvdcvector(2) = 2  ! y -> y
          mvdcvector(3) = 3  ! z -> z (out-of-plane)
        case ("L")
          mvdcvector(1) = 1  ! x -> x (field)
          mvdcvector(2) = 3  ! y -> z
          mvdcvector(3) = 2  ! z -> y (out-of-plane)
        case ("S")
          mvdcvector(1) = 3  ! x -> z (field)
          mvdcvector(2) = 1  ! y -> x
          mvdcvector(3) = 2  ! z -> y (out-of-plane)
        end select
      case ("O")
    !   Other direction:
        if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Vdc not implemented for this electric field!')")
        lvdc = .false.
        vdcneighbor(1) = 0
        vdcneighbor(2) = 0
      end select vdc_definitions_bcc110_axis
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
!------------------- Spin quantization direction -----------------------
    read(unit=magaxis,fmt=*) j
    magnetization_axis_fcc100: select case (j)
    case (1:4)
  !   In-plane 1st n.n.: (x - out-of-plane)
      phi   = dble(2*j-1)*pi/4.d0 ! around z (counter-clockwise from the top X->Y)
      theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
    case (5:8)
  !   In-plane 2nd n.n.: (x - out-of-plane)
      phi   = dble(j-5)*pi/2.d0   ! around z (counter-clockwise from the top X->Y)
      theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
    case (9)
  !   Out-of-plane: (x in the direction of dirEfield)
      select case (i)             ! around z (counter-clockwise from the top X->Y)
      case (1:4)
        phi = dble(2*i-1)*pi/4.d0
      case (5:8)
        phi = dble(i-5)*pi/2.d0
      case default
        phi = 0.d0
      end select
      theta = 0.d0                ! around y' (counter-clockwise from the top Z->X')
    case default
      if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Choose a correct magnetization axis!')")
      call MPI_Finalize(ierr)
      stop
  !     phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
  !     theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
    end select magnetization_axis_fcc100
!-------------- Longitudinal and transverse DC voltage -----------------
    vdc_definitions_fcc100: select case (i)
    case (1:4) ! If the electric field is applied along the long axis
      lvdc = .true.
      vdcneighbor(1) = i          ! Longitudinal neighbor
      vdcneighbor(2) = mod(i,4)+1 ! Transverse neighbor
      ! Magnetization components used in Vdc calculations
      ! x -> direction of the electric field
      ! y -> in-plane transverse to direction of the field
      ! z -> out-of-plane
      ! mapped in the magnetization reference frame
      if(j.eq.9) then
        mvdcvector(1) = 1  ! x -> x (field)
        mvdcvector(2) = 2  ! y -> y
        mvdcvector(3) = 3  ! z -> z (out-of-plane)
      else if(j.eq.i) then
        mvdcvector(1) = 3  ! x -> z (field)
        mvdcvector(2) = 2  ! y -> y
        mvdcvector(3) = 1  ! z -> x (out-of-plane)
      else
        if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Vdc not implemented for this direction of magnetization!')")
        lvdc = .false.
        vdcneighbor(1) = 0 ! Longitudinal neighbor?
        vdcneighbor(2) = 0 ! Transverse neighbor?
      end if
    case (5:8)
      lvdc = .true.
      vdcneighbor(1) = i            ! Longitudinal neighbor
      vdcneighbor(2) = mod(i-4,4)+5 ! Transverse neighbor
      if(j.eq.9) then
        mvdcvector(1) = 1  ! x -> x (field)
        mvdcvector(2) = 2  ! y -> y
        mvdcvector(3) = 3  ! z -> z (out-of-plane)
      else if(j.eq.i) then
        mvdcvector(1) = 3  ! x -> z (field)
        mvdcvector(2) = 2  ! y -> y
        mvdcvector(3) = 1  ! z -> x (out-of-plane)
      else
        if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Vdc not implemented for this direction of magnetization!')")
        lvdc = .false.
        vdcneighbor(1) = 0 ! Longitudinal neighbor?
        vdcneighbor(2) = 0 ! Transverse neighbor?
      end if
    case default
  !   Other direction:
      if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Vdc not implemented for this electric field!')")
      vdcneighbor(1) = 0
      vdcneighbor(2) = 0
    end select vdc_definitions_fcc100
  case("fcc111")
!------------------- Spin quantization direction -----------------------
    read(unit=magaxis,fmt=*) j
    magnetization_axis_fcc111: select case (j)
!         case (1:4)
!       !   In-plane 1st n.n.:
!           phi   = dble(2*j-1)*pi/4.d0 ! around z (counter-clockwise from the top X->Y)
!           theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
!         case (5:8)
!       !   In-plane 2nd n.n.:
!           phi   = dble(j-5)*pi/2.d0   ! around z (counter-clockwise from the top X->Y)
!           theta = pi/2.d0             ! around y' (counter-clockwise from the top Z->X')
!         case (9)
!       !   Out-of-plane:
!           phi   = pi/4.d0             ! around z (counter-clockwise from the top X->Y)
!           theta = 0.d0                ! around y' (counter-clockwise from the top Z->X')
    case default
      if(myrank.eq.0) write(outputunit,"('[lattice_definitions] Choose a correct magnetization axis!')")
      call MPI_Finalize(ierr)
      stop
  !     phi   = 0.d0 ! around z (counter-clockwise from the top X->Y)
  !     theta = 0.d0 ! around y' (counter-clockwise from the top Z->X')
    end select magnetization_axis_fcc111
!----------- Direction of applied in-plane electric field --------------
    read(unit=dirEfield,fmt=*) j
    direction_E_field_fcc111: select case (j)
!       case (1:8)   !    In plane neighbors:
!         dirEfieldvec = c0(j,:)
!       case default !    Other direction:
!         dirEfieldvec = [1.d0,0.d0,0.d0] ! In-plane
    end select direction_E_field_fcc111
  end select lattice_dependent_definitions

!-------------- Check if current is calculated for Vdc -----------------
  if((vdcneighbor(1).ne.0).and.((vdcneighbor(1).lt.n0sc1).or.(vdcneighbor(1).gt.n0sc2))) then
    if(myrank.eq.0) then
      write(outputunit,"('[lattice_definitions] Invalid longitudinal neighbor for Vdc calculation: ',i0,'!')") vdcneighbor(1)
      write(outputunit,"('[lattice_definitions] Current calculated between ',i0,' and ',i0,'.')") n0sc1,n0sc2
    end if
    call MPI_Finalize(ierr)
    stop
  end if
  if((vdcneighbor(2).ne.0).and.((vdcneighbor(2).lt.n0sc1).or.(vdcneighbor(2).gt.n0sc2))) then
    if(myrank.eq.0) then
      write(outputunit,"('[lattice_definitions] Invalid transverse neighbor for Vdc calculation: ',i0,'!')") vdcneighbor(2)
      write(outputunit,"('[lattice_definitions] Current calculated between ',i0,' and ',i0,'.')") n0sc1,n0sc2
    end if
    call MPI_Finalize(ierr)
    stop
  end if

  return
end subroutine lattice_definitions