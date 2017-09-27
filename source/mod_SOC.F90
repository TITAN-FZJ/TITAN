module mod_SOC
use mod_f90_kind, only: double
implicit none

!========================================================================================!
logical :: SOC                                          !< Turn on/off SOC
logical :: llineargfsoc = .false.
logical :: llinearsoc = .false.                         !< Linear SOC
real(double) :: socscale = 1.d0                         !< Rescale of SOC parameter
character(len=50)  :: socpart = ""
character(len=1)   :: SOCc = ""

complex(double), dimension(:,:), allocatable :: ls


contains

  subroutine initLS(theta, phi, nOrb)
    use mod_f90_kind, only: double
    use mod_constants, only: cI, sq3, cZero
    implicit none
    real(double), intent(in) :: theta, phi
    integer, intent(in) :: nOrb

    if(nOrb /= 9) stop "LS Matrix only implemented for nOrb = 9"

    if(allocated(ls)) deallocate(ls)
    allocate(ls(2*nOrb, 2*nOrb))

    ! the spin-orbit matrix
    ls = cZero

    if(SOC) then
      ! diagonal in spin
      ! p-block
      !   ls( 2, 3) = -0.5d0*cI*cos(theta)
      !   ls( 2, 4) =  0.5d0*cI*sin(theta)*sin(phi)
      !   ls( 3, 2) =  conjg(ls(2,3))
      !   ls( 3, 4) = -0.5d0*cI*sin(theta)*cos(phi)
      !   ls( 4, 2) =  conjg(ls(2,4))
      !   ls( 4, 3) =  conjg(ls(3,4))

      !   ls(11:13,11:13) = -ls(2:4,2:4)

      ! d-block
      ls( 5, 6) =  0.5d0*cI*sin(theta)*sin(phi)
      ls( 5, 7) = -0.5d0*cI*sin(theta)*cos(phi)
      ls( 5, 8) =  cI*cos(theta)
      ls( 6, 5) =  conjg(ls(5,6))
      ls( 6, 7) =  0.5d0*cI*cos(theta)
      ls( 6, 8) = -0.5d0*cI*sin(theta)*cos(phi)
      ls( 6, 9) = -0.5d0*cI*sq3*sin(theta)*cos(phi)
      ls( 7, 5) =  conjg(ls(5,7))
      ls( 7, 6) =  conjg(ls(6,7))
      ls( 7, 8) = -0.5d0*cI*sin(theta)*sin(phi)
      ls( 7, 9) =  0.5d0*cI*sq3*sin(theta)*sin(phi)
      ls( 8, 5) =  conjg(ls(5,8))
      ls( 8, 6) =  conjg(ls(6,8))
      ls( 8, 7) =  conjg(ls(7,8))
      ls( 9, 6) =  conjg(ls(6,9))
      ls( 9, 7) =  conjg(ls(7,9))

      ls(14:18,14:18) = -ls(5:9,5:9)

      ! off-diagonal in spin
      ! p-block
      !   ls( 2,12) =  0.5d0*cI*sin(theta)
      !   ls( 2,13) =  0.5d0*(cos(phi)+cI*cos(theta)*sin(phi))
      !   ls( 3,11) = -0.5d0*cI*sin(theta)
      !   ls( 3,13) =  0.5d0*(sin(phi)-cI*cos(theta)*cos(phi))
      !   ls( 4,11) = -0.5d0*(cos(phi)+cI*cos(theta)*sin(phi))
      !   ls( 4,12) = -0.5d0*(sin(phi)-cI*cos(theta)*cos(phi))

      !   ls(11:13,2:4) = transpose(conjg(ls(2:4,11:13)))

      ! d-block
      ls( 5,15) =  0.5d0*(cos(phi)+cI*cos(theta)*sin(phi))
      ls( 5,16) =  0.5d0*(sin(phi)-cI*cos(theta)*cos(phi))
      ls( 5,17) = -cI*sin(theta)
      ls( 6,14) = -0.5d0*(cos(phi)+cI*cos(theta)*sin(phi))
      ls( 6,16) = -0.5d0*cI*sin(theta)
      ls( 6,17) =  0.5d0*(sin(phi)-cI*cos(theta)*cos(phi))
      ls( 6,18) =  0.5d0*sq3*(sin(phi)-cI*cos(theta)*cos(phi))
      ls( 7,14) =  -0.5d0*(sin(phi)-cI*cos(theta)*cos(phi))
      ls( 7,15) =  0.5d0*cI*sin(theta)
      ls( 7,17) = -0.5d0*(cos(phi)+cI*cos(theta)*sin(phi))
      ls( 7,18) =  0.5d0*sq3*(cos(phi)+cI*cos(theta)*sin(phi))
      ls( 8,14) =  cI*sin(theta)
      ls( 8,15) = -0.5d0*(sin(phi)-cI*cos(theta)*cos(phi))
      ls( 8,16) =  0.5d0*(cos(phi)+cI*cos(theta)*sin(phi))
      ls( 9,15) = -0.5d0*sq3*(sin(phi)-cI*cos(theta)*cos(phi))
      ls( 9,16) = -0.5d0*sq3*(cos(phi)+cI*cos(theta)*sin(phi))

      ls(14:18,5:9) = transpose(conjg(ls(5:9,14:18)))

      !   l_x = ls(14:18, 5: 9)+ls( 5: 9,14:18)
      !   l_y = -cI*(ls(14:18, 5: 9)-ls( 5: 9,14:18))
      !   l_z = 2.d0*ls( 5: 9, 5: 9)
    end if
    return
  end subroutine initLS

end module mod_SOC
