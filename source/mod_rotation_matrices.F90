! Rotation matrices
module mod_rotation_matrices
  
contains

  subroutine rotation_matrices_chi(theta_in,phi_in,rotmat,iflag)
    ! Build rotation matrices for the spin-spin susceptibility in the
    ! +,up,down,- basis: chi' = A.chi.B
    ! Due to the form of the susceptibility, the left
    ! rotation matrix (A) is different than the right one (B).
    ! iflag = 1 calculates A
    ! iflag = 2 calculates B
    ! Angles should be given in degrees
    use mod_kind, only: dp
    use mod_constants, only: cI,deg2rad
    implicit none
    integer,                        intent(in)  :: iflag
    real(dp),                   intent(in)  :: theta_in,phi_in
    real(dp)                                :: theta,phi,onepluscos,oneminuscos,sintheta
    complex(dp)                             :: expphi,expmphi
    complex(dp),dimension(4,4), intent(out) :: rotmat

    theta = theta_in*deg2rad
    phi   = phi_in*deg2rad

    expphi      = exp(cI*phi)
    expmphi     = conjg(expphi)
    onepluscos  = 1._dp + cos(theta)
    oneminuscos = 1._dp - cos(theta)
    sintheta    = sin(theta)

    if(iflag==1) then ! left matrix
      rotmat(1,1) = expmphi*onepluscos
      rotmat(1,2) =-sintheta
      rotmat(1,3) = sintheta
      rotmat(1,4) =-expphi*oneminuscos

      rotmat(2,1) = expmphi*sintheta
      rotmat(2,2) = onepluscos
      rotmat(2,3) = oneminuscos
      rotmat(2,4) = expphi*sintheta

      rotmat(3,1) =-expmphi*sintheta
      rotmat(3,2) = oneminuscos
      rotmat(3,3) = onepluscos
      rotmat(3,4) =-expphi*sintheta

      rotmat(4,1) =-expmphi*oneminuscos
      rotmat(4,2) =-sintheta
      rotmat(4,3) = sintheta
      rotmat(4,4) = expphi*onepluscos
    else ! right matrix
      rotmat(1,1) = expphi*onepluscos
      rotmat(2,1) =-sintheta
      rotmat(3,1) = sintheta
      rotmat(4,1) =-expmphi*oneminuscos

      rotmat(1,2) = expphi*sintheta
      rotmat(2,2) = onepluscos
      rotmat(3,2) = oneminuscos
      rotmat(4,2) = expmphi*sintheta

      rotmat(1,3) =-expphi*sintheta
      rotmat(2,3) = oneminuscos
      rotmat(3,3) = onepluscos
      rotmat(4,3) =-expmphi*sintheta

      rotmat(1,4) =-expphi*oneminuscos
      rotmat(2,4) =-sintheta
      rotmat(3,4) = sintheta
      rotmat(4,4) = expmphi*onepluscos
    end if

    rotmat = 0.5_dp*rotmat

  end subroutine rotation_matrices_chi

  ! Rotation matrix of an angle theta (in degrees) around y axis
  subroutine rotation_matrix_ry(theta,ry)
    use mod_kind, only: dp
    use mod_constants, only: deg2rad
    implicit none
    real(dp),                intent(in)  :: theta
    real(dp),dimension(3,3), intent(out) :: ry

    ry = 0._dp
    ry(1,1) = cos(theta*deg2rad)
    ry(1,3) = sin(theta*deg2rad)
    ry(2,2) = 1._dp
    ry(3,1) =-sin(theta*deg2rad)
    ry(3,3) = cos(theta*deg2rad)

  end subroutine rotation_matrix_ry

  ! Rotation matrix of an angle phi (in degrees)  around z axis
  subroutine rotation_matrix_rz(phi,rz)
    use mod_kind, only: dp
    use mod_constants, only: deg2rad
    implicit none
    real(dp),                intent(in)  :: phi
    real(dp),dimension(3,3), intent(out) :: rz

    rz = 0._dp
    rz(1,1) = cos(phi*deg2rad)
    rz(1,2) =-sin(phi*deg2rad)
    rz(2,1) = sin(phi*deg2rad)
    rz(2,2) = cos(phi*deg2rad)
    rz(3,3) = 1._dp

  end subroutine rotation_matrix_rz
end module mod_rotation_matrices