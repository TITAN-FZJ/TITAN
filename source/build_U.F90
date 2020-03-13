! This subroutine builds the U vector from the elemental files
subroutine build_U(sys)
  use mod_parameters, only: U, Un, Um
  use mod_System,     only: System
  implicit none
  type(System), intent(in) :: sys
  integer :: i

  do i=1,sys%nAtoms
    U(i) = sys%Types(sys%Basis(i)%Material)%U
  end do
  
  do i=1,sys%nAtoms
    Un(i) = sys%Types(sys%Basis(i)%Material)%Un
  end do
  
  do i=1,sys%nAtoms
    Um(i) = sys%Types(sys%Basis(i)%Material)%Um
  end do
end subroutine build_U
