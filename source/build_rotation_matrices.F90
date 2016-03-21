! Build rotation matrices for the spin-spin susceptibility in the
! +,up,down,- basis: chi' = A.chi.B
! Due to the form of the susceptibility, the left
! rotation matrix (A) is different than the right one (B).
! iflag = 1 calculates A
! iflag = 2 calculates B
subroutine build_rotation_matrices(theta,phi,rotmat,iflag)
  use mod_f90_kind
  use mod_constants, only: zi
  implicit none
  integer,                        intent(in)  :: iflag
  real(double),                   intent(in)  :: theta,phi
  real(double)                                :: onepluscos,oneminuscos,sintheta
  complex(double)                             :: expphi,expmphi
  complex(double),dimension(4,4), intent(out) :: rotmat

  expphi      = exp(zi*phi)
  expmphi     = conjg(expphi)
  onepluscos  = 1.d0 + cos(theta)
  oneminuscos = 1.d0 - cos(theta)
  sintheta    = sin(theta)

  if(iflag.eq.1) then
    rotmat(1,1) = expmphi*onepluscos
    rotmat(1,2) =-expmphi*sintheta
    rotmat(1,3) =-rotmat(1,2)
    rotmat(1,4) =-expmphi*oneminuscos

    rotmat(2,1) = sintheta
    rotmat(2,2) = onepluscos
    rotmat(2,3) = oneminuscos
    rotmat(2,4) = sintheta

    rotmat(3,1) =-sintheta
    rotmat(3,2) = oneminuscos
    rotmat(3,3) = onepluscos
    rotmat(3,4) =-sintheta

    rotmat(4,1) = conjg(rotmat(1,4))
    rotmat(4,2) = conjg(rotmat(1,2))
    rotmat(4,3) = conjg(rotmat(1,3))
    rotmat(4,4) = conjg(rotmat(1,1))
  else
    rotmat(1,4) =-expmphi*oneminuscos
    rotmat(2,4) =-expmphi*sintheta
    rotmat(3,4) =-rotmat(2,4)
    rotmat(4,4) = expmphi*onepluscos

    rotmat(1,3) =-sintheta
    rotmat(2,3) = oneminuscos
    rotmat(3,3) = onepluscos
    rotmat(4,3) =-sintheta

    rotmat(1,2) = sintheta
    rotmat(2,2) = onepluscos
    rotmat(3,2) = oneminuscos
    rotmat(4,2) = sintheta

    rotmat(1,1) = conjg(rotmat(4,4))
    rotmat(2,1) = conjg(rotmat(2,4))
    rotmat(3,1) = conjg(rotmat(3,4))
    rotmat(4,1) = conjg(rotmat(1,4))
  end if

  rotmat = 0.5d0*rotmat

  return
end subroutine build_rotation_matrices
