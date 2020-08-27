! This subroutine builds the U vector from the elemental files
subroutine build_U(sys)
  use mod_parameters, only: Un, Um
  use mod_System,     only: System_type
  implicit none
  type(System_type), intent(in) :: sys
  integer :: i
  
  do i=1,sys%nAtoms
    Un(i) = sys%Types(sys%Basis(i)%Material)%Un
  end do
  
  do i=1,sys%nAtoms
    Um(i) = sys%Types(sys%Basis(i)%Material)%Um
  end do
end subroutine build_U
