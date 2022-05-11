module mod_rotation_matrices
!! Module containing rotation matrices

contains

  subroutine rotation_matrices_chi(theta_in,phi_in,rotmat,iflag)
    !! Build rotation matrices for the spin-spin susceptibility in the
    !! +,up,down,- basis: chi_r = A.chi.B
    !! Due to the form of the susceptibility, the left
    !! rotation matrix (A) is different than the right one (B).
    !! iflag = 1 calculates A
    !! iflag = 2 calculates B
    !! Angles should be given in degrees
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

  subroutine rotation_matrix_ry(theta,ry)
    !! Rotation matrix of an angle theta (in degrees) around y axis
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

  subroutine rotation_matrix_rz(phi,rz)
    !! Rotation matrix of an angle phi (in degrees)  around z axis
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

  subroutine rot_rodrigues(u0,u1,rotmat)
    !! Rotation matrix rotating vector u0 to vector u1
    !! Rodrigues rotation formula
    use mod_kind, only: dp
    implicit none
    real(dp), intent(in)  :: u0(3)
    real(dp), intent(in)  :: u1(3)
    real(dp), intent(out) :: rotmat(3,3)
  ! ---------------------------------------
    real(dp), parameter :: tol = 1.e-6_dp, rtol = 1.e-8_dp
    real(dp)  :: axis(3), u0len, u1len, u01len, cosa, sina, cross(3,3), ui(3), fac
    integer   :: i
  
  ! Lengths of the vectors (if not unit vectors)
    u0len = sqrt(dot_product(u0,u0))
    u1len = sqrt(dot_product(u1,u1))
    !if (u0len < tol .or. u1len < tol) stop 'rotvec: check u0 and u1'
    cosa = dot_product(u0,u1)/(u0len*u1len)
  ! Rotation axis: axis = u0 x u1 / | u0 x u1 |
    axis(1) = u0(2)*u1(3) - u0(3)*u1(2)
    axis(2) = u0(3)*u1(1) - u0(1)*u1(3)
    axis(3) = u0(1)*u1(2) - u0(2)*u1(1)
  ! chop small numbers
    where (abs(axis) < tol) axis = 0.e0_dp
    u01len = sqrt(dot_product(axis,axis))
  ! ----------------------------------------------------------------------
    if (u01len > rtol) then
  !   noncollinear: no problem
      sina = u01len/(u0len*u1len)
      axis = axis/u01len
      fac = 1.e0_dp - cosa
    else
  !   collinear: handle zeros
      if (cosa > 0.e0_dp) then
  !     cosa = +1 => no rotation
        cosa = +1.e0_dp; sina = 0.e0_dp
        axis(:) = 0.e0_dp
        fac = 0.e0_dp
      else
  !     cosa = -1 => inversion
        cosa = -1.e0_dp; sina = 0.e0_dp
        fac = 2.e0_dp
  !     choose a rotation axis
        i = minloc(abs(u0),dim=1)
        ui(:) = 0.e0_dp
        ui(i) = 1.e0_dp
        axis(1) = u0(2)*ui(3) - u0(3)*ui(2)
        axis(2) = u0(3)*ui(1) - u0(1)*ui(3)
        axis(3) = u0(1)*ui(2) - u0(2)*ui(1)
        u01len = sqrt(dot_product(axis,axis))
        axis(:) = axis(:)/u01len
      end if
    end if
  ! ----------------------------------------------------------------------
  ! cross product matrix
    cross(:,:) = 0.e0_dp
    cross(2,1) =  axis(3); cross(1,2) = -axis(3)
    cross(3,1) = -axis(2); cross(1,3) =  axis(2)
    cross(3,2) =  axis(1); cross(2,3) = -axis(1)
  ! rotation matrix
    rotmat(:,:) = sina*cross(:,:) + fac*matmul(cross,cross)
    do i=1,3
      rotmat(i,i) = rotmat(i,i) + 1.e0_dp
    end do
  ! All done!
  end subroutine rot_rodrigues

end module mod_rotation_matrices