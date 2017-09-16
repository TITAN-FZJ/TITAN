!------------------------------------------------------------------------------------!
! TITAN - Time-dependent Transport and Angular momentum properties of Nanostructures !
!------------------------------------------------------------------------------------!
!
! MODULE: mod_monoBasis
!
!> @author
!> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> Generate polyatomic Basis from Layer information
!
! REVISION HISTORY:
! 06 July 2017 - Initial Version
!------------------------------------------------------------------------------------!

module mod_monoBasis
  use mod_f90_kind, only: double
  implicit none

contains

  subroutine create_basis()
    ! TODO
    ! Take cZero position
    ! Generate bubble around it
    ! Take closest out of plane atom
    ! Write current shift current cZero by closest atom and write to array
    ! Find closest out of plane atoms, if more than one, find the one closest to all prior basis atoms
  end subroutine create_basis

end module mod_monoBasis
