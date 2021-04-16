!-------------------------------------------------------------------------------
! TITAN - Time-dependent description of Itinerant electrons: Transport and Angular momentum properties of Nanostructures
!-------------------------------------------------------------------------------
!
! MODULE: mod_dft
!
!> @author
!> Filipe Guimaraes, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> This module deals with the input parameters from DFT calculations
!> initially from PAOFLOW
!
! REVISION HISTORY:
! 25 February 2020 - Initial Version
!-------------------------------------------------------------------------------
module mod_dft
  implicit none
contains

  subroutine readHamiltonian(s,filename)
  !! Reading element file, including all the parameters
    use mod_kind,      only: dp
    use mod_system,    only: System_type
    use AtomTypes,     only: NeighborAtom
    use mod_io,        only: log_error
    use mod_tools,     only: ItoS,StoI,StoR,StoArray,next_line
    implicit none
    type(System_type),   intent(inout) :: s
    character(len=*),    intent(in)    :: filename

    character(len=100) :: title

    integer  :: f_unit = 9954
    integer  :: nOrb
    !! Number of orbitals read from file
    integer  :: nNeighbors
    !! Number of neighbors
    integer  :: nCells,nCells_per_dim
    !! Total number of cells and per dimension
    integer  :: dim_ini(3)=1, dim_end(3)=1, pos(3), orb(2)
    !! Dimensions and temporary index variables
    real(dp) :: hop(2)
    !! Temporary variable to store complex hoppings
    integer, allocatable :: weight(:)
    !! Weight of each cell
    complex(dp), dimension(:,:,:,:),  allocatable :: hoppings
    !! Hopping parameters read from the file
    type(NeighborAtom), dimension(:), allocatable :: list
    !! Temporary list to store neighbor information
    integer  :: x,y,z,i,j,i0,i1,mu,nu,ios,wgt_per_line = 15

    ! Opening file
    open(f_unit, file=trim(filename), status='old', iostat=ios)
    if(ios /= 0) call log_error("readHamiltonian", "Error occured when trying to read file " // trim(filename))

    ! Read title
    read(f_unit, fmt='(A)', iostat=ios) title
    title = trim(adjustl(title))

    ! Read number of orbitals * atoms in the unit cell
    nOrb = StoI(next_line("readHamiltonian",f_unit,"number of orbitals"))
    nOrb = int(nOrb/s%nAtoms)
    if(nOrb /= s%nOrb) call log_error("readHamiltonian", "Number of orbitals selected in input (" // trim(itos(s%nOrb)) // ") is different than in hamiltonian file (" // trim(itos(nOrb)) // ")")

    ! Read number of Cells
    nCells = StoI(next_line("readHamiltonian",f_unit,"number of cells"))

    ! Number of Cells per dimension
    nCells_per_dim = ceiling(nCells**(1.d0/s%isysdim))
    do i = 1,s%isysdim
      dim_ini(i) = 1
      dim_end(i) = nCells_per_dim
    end do

    ! Allocate array for all atoms in all unit cells
    allocate(list(nCells*s%nAtoms))
    do i = 1,nCells*s%nAtoms
      allocate(list(i)%t0i(nOrb,nOrb,s%nAtoms))
      allocate(list(i)%isHopping(s%nAtoms))
      list(i)%t0i = 0._dp
      list(i)%isHopping = .false.
    end do
    do i = 1, s%nTypes
      allocate(s%Types(i)%onSite(nOrb,nOrb))
    end do

    ! Read weight of each cell
    allocate(weight(nCells))
    j = floor(dble(nCells/wgt_per_line))
    do i = 1,j
      i0 = (i-1)*wgt_per_line+1
      i1 = i0+wgt_per_line-1
      weight(i0:i1) = StoI(next_line("readHamiltonian",f_unit,"weights"),wgt_per_line)
    end do
    i1 = wgt_per_line*j
    weight(i1+1:i1+mod(nCells,wgt_per_line)) = StoI(next_line("readHamiltonian",f_unit,"weights"),mod(nCells,wgt_per_line)-1)    

    ! Allocating temporary hopping matrix
    allocate(hoppings(s%nAtoms,s%nAtoms,nOrb,nOrb))


    ! "nNeighbors" is added up until all atoms in all unit cells, nCells*s%nAtoms
    nNeighbors = 0
    do x = dim_ini(1),dim_end(1)
      do y = dim_ini(2),dim_end(2)
        do z = dim_ini(3),dim_end(3)
          do i = 1,s%nAtoms

            ! "nNeighbors" is the current atom
            nNeighbors = nNeighbors + 1

            do mu = 1,nOrb
              ! Reading hamiltonian for a given unit cell
              do j = 1,s%nAtoms
                do nu = 1,nOrb
                  ! Reading information from file
                  read(f_unit, fmt=*, iostat=ios) pos(1), pos(2), pos(3), orb(1), orb(2), hop(1), hop(2)
                  ! hoppings(i,j,mu,nu) = cmplx(hop(1),hop(2),dp)
                  if((pos(1)==0).and.(pos(2)==0).and.(pos(3)==0).and.(i==j)) then
                    s%Types(s%Basis(i)%Material)%onSite(mu,nu) = cmplx(hop(1),hop(2),dp)
                    cycle
                  end if

                  list(nNeighbors)%isHopping(j) = .true.
                  list(nNeighbors)%t0i(mu,nu,j) = cmplx(hop(1),0.e0_dp,dp)
                end do !nu
              end do !j
            end do ! mu

            ! list(nNeighbors)%isHopping(i) = .true.

            ! Determine the indices of the unit-cell containing atom "size"
            list(nNeighbors)%Cell(1) = pos(1)
            list(nNeighbors)%Cell(2) = pos(2)
            list(nNeighbors)%Cell(3) = pos(3)
            ! Cell position is R_i = i*a1 + j*a2 + k*a3
            list(nNeighbors)%CellVector = list(nNeighbors)%Cell(1) * s%a1 &
                                        + list(nNeighbors)%Cell(2) * s%a2 &
                                        + list(nNeighbors)%Cell(3) * s%a3
            ! Atom position is r = R_i + r_j
            list(nNeighbors)%Position = s%Basis(i)%Position + list(nNeighbors)%CellVector
            list(nNeighbors)%BasisIndex = i
            list(nNeighbors)%Material = s%Basis(i)%Material
          end do ! i
        end do ! z
      end do ! y
    end do ! x

    close(f_unit)
  
    ! Defining number of neighbors (total number of atoms)
    ! TODO: If cell is to be cut, a temporary array must be created here, 
    ! and the exact number of neighbors can be counted in the loop below
    s%nNeighbors = nCells*s%nAtoms
    if(allocated(s%Neighbors)) deallocate(s%Neighbors)
    allocate(s%Neighbors(s%nNeighbors))

    s%Neighbors(1:s%nNeighbors) = list(1:s%nNeighbors)

  end subroutine readHamiltonian

end module mod_dft