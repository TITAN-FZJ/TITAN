! Calculate the derivative of H0 of a slab containing
! Npl layers (ES not included)
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!        1     2     3       Npl-2 Npl-1  Npl
!         <-S-> <S-1>           <S-1> <-S->
subroutine dtdksub(kp,dtdk)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_tight_binding
  implicit none
  integer     :: i
  real(double), intent(in)  :: kp(3)
  complex(double),dimension(9,9)    :: h00,h01,h10
  complex(double),dimension(Npl,Npl,9,9),intent(out)  :: dtdk

  dtdk = zero

! Mouting derivative of slab's hamiltonian
  do i=1,Npl
    call helpdtdk(kp,h00,h01,h10,i+1)

    dtdk(i,i,:,:) = h00

    if (i.ne.Npl) then
      dtdk(i,i+1,:,:) = h01(:,:)
      dtdk(i+1,i,:,:) = h10(:,:)
    end if
  end do

  return
end subroutine dtdksub