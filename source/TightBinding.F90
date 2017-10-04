module TightBinding
  use mod_f90_kind, only: double
  implicit none

  integer :: tbmode ! (1) Slater-Koster, (2) DFT

  integer, parameter :: nOrb = 9
  integer, parameter :: nOrb2 = 18
  logical, dimension(9) :: Orbitals

  integer :: fermi_layer   ! Maybe remove it (read from outside)
  character(len=20), dimension(:), allocatable :: layers
contains

  subroutine initTightBinding(s)
    use mod_system, only: System
    use SK_TightBinding, only: SKParam => get_parameter
    implicit none
    type(System), intent(inout) :: s
    if(tbmode == 1) then
      call SKParam(s, fermi_layer, nOrb, Orbitals)
    else if(tbmode == 2) then
      stop "Not Implemented"
      !call get_DFT_hopping()
    end if
  end subroutine initTightBinding

end module TightBinding