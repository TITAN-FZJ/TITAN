! This subroutine builds the U vector from the elemental files
subroutine build_U(sys)
  use mod_parameters, only: U
  use mod_System,     only: System
  implicit none
  type(System), intent(in) :: sys
  integer :: i

  do i=1,sys%nAtoms
    U(i) = sys%Types(sys%Basis(i)%Material)%U
  end do
end subroutine build_U
