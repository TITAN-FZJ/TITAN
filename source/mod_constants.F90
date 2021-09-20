module mod_constants
  !! Contains commonly used constants
  use mod_kind, only: dp
  implicit none
  complex(dp), parameter :: cZero=cmplx(0._dp,0._dp,dp)
  !! Complex scalar cZero
  complex(dp), parameter :: cOne=cmplx(1._dp,0._dp,dp)
  !! Complex scalar one
  complex(dp), parameter :: cI=cmplx(0._dp,1._dp,dp)
  !! Complex scalar i
  real(dp), parameter :: pi = 4._dp * atan(1._dp)
  !! Pi
  real(dp), parameter :: tpi = 2._dp * pi
  !! 2*Pi
  real(dp), parameter :: sq2 = sqrt(2._dp)
  !! Square root of 2
  real(dp), parameter :: hsq2 = 0.5_dp * sq2
  !! (Square root of 2) / 2
  real(dp), parameter :: sq3 = sqrt(3._dp)
  !! Square root of 3
  real(dp), parameter :: deg2rad = pi/180._dp
  !! degrees to radians
  real(dp), parameter :: rad2deg = 180._dp/pi
  !! radians to degrees
  complex(dp) :: pauli_mat(2,2,0:5)
  !! Identity and pauli matrices  (0,x,y,z,+,-)
#ifdef _GPU
  complex(dp), device :: pauli_mat_d(2,2,0:5)
  !! Identity and pauli matrices  (0,x,y,z,+,-) on the GPUs
#endif
  real(dp)    :: levi_civita(3,3,3)
  !! Levi Civita Tensor

  complex(dp), dimension(4,4) :: StoC = reshape([ cmplx(0.0_dp,0._dp,dp), cmplx(0.5_dp,0._dp,dp),  cmplx(0.0_dp,-0.5_dp,dp),  cmplx( 0.0_dp,0._dp,dp), &
                                                  cmplx(1.0_dp,0._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx(0.0_dp, 0.0_dp,dp),  cmplx( 0.5_dp,0._dp,dp), &
                                                  cmplx(1.0_dp,0._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx(0.0_dp, 0.0_dp,dp),  cmplx(-0.5_dp,0._dp,dp), &
                                                  cmplx(0.0_dp,0._dp,dp), cmplx(0.5_dp,0._dp,dp),  cmplx(0.0_dp, 0.5_dp,dp),  cmplx( 0.0_dp,0._dp,dp)], [4,4])
  !! Transformation matrix spin (+,up,down,-) to cartesian (0,x,y,z)

  complex(dp), dimension(4,4) :: CtoS = reshape([ cmplx(0.0_dp,0._dp,dp), cmplx(0.5_dp,0._dp,dp),  cmplx( 0.5_dp, 0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp), &
                                                  cmplx(1.0_dp,0._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp),  cmplx( 1.0_dp, 0._dp,dp), &
                                                  cmplx(0.0_dp,1._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp),  cmplx( 0.0_dp,-1._dp,dp), &
                                                  cmplx(0.0_dp,0._dp,dp), cmplx(1.0_dp,0._dp,dp),  cmplx(-1.0_dp, 0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp)], [4,4])
  !! Transformation matrix cartesian (0,x,y,z) to spin (+,up,down,-)
contains

  subroutine define_constants(s)
    use mod_System, only: System_type
    implicit none
    type(System_type), intent(inout) :: s
    integer :: i,mu,nu,mud,nOrb,nOrb2
    
    ! All four Pauli matrices in spin space
    ! identity
    pauli_mat(1,:,0) = [cOne, cZero]
    pauli_mat(2,:,0) = [cZero, cOne]
    ! sigma_x
    pauli_mat(1,:,1) = [cZero,cOne]
    pauli_mat(2,:,1) = [cOne, cZero]
    ! sigma_y
    pauli_mat(1,:,2) = [cZero,-cI]
    pauli_mat(2,:,2) = [cI,cZero]
    ! sigma_z
    pauli_mat(1,:,3) = [cOne,cZero]
    pauli_mat(2,:,3) = [cZero,-cOne]
    ! sigma_+
    pauli_mat(:,:,4) = pauli_mat(:,:,1) + cI*pauli_mat(:,:,2)
    ! sigma_-
    pauli_mat(:,:,5) = pauli_mat(:,:,1) - cI*pauli_mat(:,:,2)

#ifdef _GPU
    pauli_mat_d = pauli_mat
#endif

    levi_civita(1,1,:) = [ 0._dp , 0._dp , 0._dp ]
    levi_civita(1,2,:) = [ 0._dp , 0._dp , 1._dp ]
    levi_civita(1,3,:) = [ 0._dp ,-1._dp , 0._dp ]
    levi_civita(2,1,:) = [ 0._dp , 0._dp ,-1._dp ]
    levi_civita(2,2,:) = [ 0._dp , 0._dp , 0._dp ]
    levi_civita(2,3,:) = [ 1._dp , 0._dp , 0._dp ]
    levi_civita(3,1,:) = [ 0._dp , 1._dp , 0._dp ]
    levi_civita(3,2,:) = [-1._dp , 0._dp , 0._dp ]
    levi_civita(3,3,:) = [ 0._dp , 0._dp , 0._dp ]

    ! Element-dependent matrices
    do i = 1, s%nTypes
      nOrb = s%Types(i)%nOrb
      nOrb2 = 2*nOrb
      allocate( s%Types(i)%ident_norb (nOrb,nOrb),&
                s%Types(i)%ident_norb2(nOrb2,nOrb2),&
                s%Types(i)%ident_dorb (nOrb2,nOrb2),&
                s%Types(i)%pauli_orb  (0:5,nOrb2,nOrb2),&
                s%Types(i)%pauli_dorb (0:5,nOrb2,nOrb2)&
                )

      s%Types(i)%ident_norb  = cZero
      s%Types(i)%ident_norb2 = cZero
      s%Types(i)%ident_dorb  = cZero
      s%Types(i)%pauli_orb   = cZero
      s%Types(i)%pauli_dorb  = cZero

      ! identity and Pauli matrices for given element (nOrb x nOrb matrix)
      do mu = 1,nOrb
        s%Types(i)%ident_norb (mu,mu) = cOne
        s%Types(i)%ident_norb2(mu,mu) = cOne
        nu = mu + s%Types(i)%nOrb
        s%Types(i)%ident_norb2(nu,nu) = cOne

        ! pauli matrix x
        s%Types(i)%pauli_orb(1,mu,nu) = cOne
        s%Types(i)%pauli_orb(1,nu,mu) = cOne
        ! pauli matrix y
        s%Types(i)%pauli_orb(2,mu,nu) = -cI
        s%Types(i)%pauli_orb(2,nu,mu) = cI
        ! pauli matrix z
        s%Types(i)%pauli_orb(3,mu,mu) = cOne
        s%Types(i)%pauli_orb(3,nu,nu) = -cOne
      end do
      ! identity and Pauli matrices for d-orbitals only:
      do mud=1,s%Types(i)%ndOrb
        mu = s%Types(i)%dOrbs(mud)

        ! Identity
        s%Types(i)%ident_dorb(mu,mu)   = cOne
        s%Types(i)%ident_dorb(nu,nu)   = cOne
        ! pauli matrix x
        s%Types(i)%pauli_dorb(1,mu,nu) = cOne
        s%Types(i)%pauli_dorb(1,nu,mu) = cOne
        ! pauli matrix y
        s%Types(i)%pauli_dorb(2,mu,nu) = -cI
        s%Types(i)%pauli_dorb(2,nu,mu) = cI
        ! pauli matrix z
        s%Types(i)%pauli_dorb(3,mu,mu) = cOne
        s%Types(i)%pauli_dorb(3,nu,nu) = -cOne
      end do
      ! Idendity as Pauli matrix 0
      s%Types(i)%pauli_orb (0,:,:) = s%Types(i)%ident_norb2(:,:)
      s%Types(i)%pauli_dorb(0,:,:) = s%Types(i)%ident_dorb (:,:)
      ! Pauli + and - matrices
      s%Types(i)%pauli_orb (4,:,:) = s%Types(i)%pauli_orb (1,:,:) + cI * s%Types(i)%pauli_orb (2,:,:)
      s%Types(i)%pauli_orb (5,:,:) = s%Types(i)%pauli_orb (1,:,:) - cI * s%Types(i)%pauli_orb (2,:,:)
      s%Types(i)%pauli_dorb(4,:,:) = s%Types(i)%pauli_dorb(1,:,:) + cI * s%Types(i)%pauli_dorb(2,:,:)
      s%Types(i)%pauli_dorb(5,:,:) = s%Types(i)%pauli_dorb(1,:,:) - cI * s%Types(i)%pauli_dorb(2,:,:)
    end do

  end subroutine define_constants

  integer function delta(i,j)
    implicit none
    integer :: i,j
    delta = 0
    if(i==j) delta = 1
  end function delta
end module mod_constants
