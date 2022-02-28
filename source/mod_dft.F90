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
    use mod_kind,      only: dp,sp,int32,int64
    use mod_system,    only: System_type
    use AtomTypes,     only: NeighborAtom
    use mod_constants, only: cZero
    use mod_io,        only: log_error
    use mod_tools,     only: ItoS,RtoS,StoI,StoR,StoArray,next_line,get_string_size
    implicit none
    type(System_type),   intent(inout) :: s
    character(len=*),    intent(in)    :: filename

    character(len=100) :: title,temp
    !! Title of the hamiltonian file, and temporary variable
    logical :: tb
    !! Flag to indicate which kind of hamiltonian file to read 
    !! (.true. -> wannier tb or .false. -> wannier/paoflow)
    integer  :: f_unit = 9954
    !! File unit to use
    integer  :: nOrb,cnt
    !! Number of orbitals read from file
    integer  :: nNeighbors
    !! Number of neighbors
    integer  :: nCells
    !! Total number of cells
    integer  :: pos(3), orb(2)
    !! Temporary index variables
    real(dp) :: hop(2)
    !! Temporary variable to store complex hoppings
    integer, allocatable :: weight(:)
    !! Weight of each cell
    type(NeighborAtom), dimension(:), allocatable :: list
    !! Temporary list to store neighbor information
    integer  :: cell,i,j,i0,i1,mu,nu,ios,wgt_per_line = 15

    ! Opening file
    open(f_unit, file=trim(filename), status='old', iostat=ios)
    if(ios /= 0) call log_error("readHamiltonian", "Error occured when trying to read file " // trim(filename))

    ! Read title
    read(f_unit, fmt='(A)', iostat=ios) title

    ! Read total number of orbitals / number of wannier functions per unit cell
    temp = next_line("readHamiltonian",f_unit,"number of orbitals or first Bravais vector")
    cnt = get_string_size(temp)
    if(cnt==1) then
      tb = .false.
    else if(cnt==3) then
      tb = .true.
      ! Bravais vectors are not used from here
      temp = next_line("readHamiltonian",f_unit,"second Bravais vector")
      temp = next_line("readHamiltonian",f_unit,"third Bravais vector")
      ! Reading number of orbitals
      temp = next_line("readHamiltonian",f_unit,"number of orbitals")
    else
      call log_error("readHamiltonian", "Invalid hamiltonian file. First line contains" // itos(cnt) // "elements (should be 1 or 3).")
    end if
    nOrb = StoI(temp)
    if(nOrb /= s%total_nOrb) call log_error("readHamiltonian", "Total number of orbitals from input (" // trim(itos(s%total_nOrb)) // ") is different than in hamiltonian file (" // trim(itos(nOrb)) // ")")

    ! Read number of Cells / number of Wigner-Seitz points
    nCells = StoI(next_line("readHamiltonian",f_unit,"number of cells"))

    ! Allocate array for all atoms in all unit cells to store degeneracies of Wigner-Seitz points
    allocate(list(nCells*s%nAtoms))
    do i = 1,nCells*s%nAtoms
      allocate(list(i)%isHopping(s%nAtoms))
      list(i)%isHopping = .false.
    end do
    do i = 1, s%nTypes
      allocate(s%Types(i)%onSite(s%Types(i)%nOrb,s%Types(i)%nOrb))
    end do

    ! Read weight of each cell / degeneracies of Wigner-Seitz points
    allocate(weight(nCells))
    j = floor(dble(nCells/wgt_per_line))
    do i = 1,j
      i0 = (i-1)*wgt_per_line+1
      i1 = i0+wgt_per_line-1
      weight(i0:i1) = StoI(next_line("readHamiltonian",f_unit,"weights"),wgt_per_line)
    end do
    if (mod(nCells,wgt_per_line)/=0) then
      i1 = wgt_per_line*j
      weight(i1+1:i1+mod(nCells,wgt_per_line)) = StoI(next_line("readHamiltonian",f_unit,"weights"),mod(nCells,wgt_per_line))    
    end if

    ! "nNeighbors" is added up until all atoms in all unit cells, nCells*s%nAtoms
    nNeighbors = 0
    do cell = 1,nCells
      if(tb) &
        read(f_unit, fmt=*, iostat=ios) pos(1), pos(2), pos(3)
      do j = 1,s%nAtoms
        
        ! "nNeighbors" is the current atom (composed by {cell,i})
        nNeighbors = nNeighbors + 1
        
        do nu = 1,s%Types(s%Basis(j)%Material)%nOrb
          ! Reading hamiltonian for a given unit cell
          do i = 1,s%nAtoms
            if(.not.allocated(list(nNeighbors)%t0i)) then
              ! Allocating with the maximum number of orbitals, to be able to store
              ! the largest matrix
              allocate(list(nNeighbors)%t0i(s%nOrb,s%nOrb,s%nAtoms))
              list(nNeighbors)%t0i = cZero
            end if

            do mu = 1,s%Types(s%Basis(i)%Material)%nOrb
              ! Reading information from file
              if(tb) then
                read(f_unit, fmt=*, iostat=ios) orb(1), orb(2), hop(1), hop(2)
              else
                read(f_unit, fmt=*, iostat=ios) pos(1), pos(2), pos(3), orb(1), orb(2), hop(1), hop(2)
              end if
              if((pos(1)==0).and.(pos(2)==0).and.(pos(3)==0).and.(i==j)) then
                if((mu==nu).and.(hop(2)>1.0e-12_dp)) &
                  call log_error("readHamiltonian", "On-site, on-orbital term for i = j = " // trim(itos(i)) // ", mu = nu = " // trim(itos(mu)) // " is not real: Im(H) = " // trim(rtos(hop(2),"(f7.2)")))

                s%Types(s%Basis(i)%Material)%onSite(mu,nu) = cmplx(hop(1),hop(2),dp)/weight(cell)
                cycle
              end if

              list(nNeighbors)%isHopping(i) = .true.
              list(nNeighbors)%t0i(nu,mu,i) = cmplx(hop(1),hop(2),dp)/weight(cell)
            end do !mu
          end do !i
        end do ! nu

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
        list(nNeighbors)%Position = s%Basis(j)%Position + list(nNeighbors)%CellVector

        list(nNeighbors)%BasisIndex = j
        list(nNeighbors)%Material = s%Basis(j)%Material
      end do ! j
    end do ! cell

    close(f_unit)
  
    ! Defining number of neighbors (total number of atoms)
    ! TODO: If number of neighbors is to be cut, a temporary array must be created here, 
    ! and the exact number of neighbors can be counted in a loop below
    s%nNeighbors = nCells*s%nAtoms
    if(allocated(s%Neighbors)) deallocate(s%Neighbors)
    allocate(s%Neighbors(s%nNeighbors))

    s%Neighbors(1:s%nNeighbors) = list(1:s%nNeighbors)

  end subroutine readHamiltonian

end module mod_dft
