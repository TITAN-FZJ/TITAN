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
  real(dp)    :: levi_civita(3,3,3)
  !! Levi Civita Tensor
  complex(dp), dimension(:,:), allocatable :: ident_norb
  !! Identity in orbital space
  complex(dp), dimension(:,:), allocatable :: ident_norb2
  !! Identity in spin and orbital space
  complex(dp), dimension(:,:), allocatable :: ident_dorb
  !! Identity in spin and orbital space (only non-zero on d-orbitals)
  complex(dp) :: pauli_mat(2,2,0:5)
  !! Identity and pauli matrices  (0,x,y,z,+,-)
#ifdef _GPU
  complex(dp), device :: pauli_mat_d(2,2,0:5)
  !! Identity and pauli matrices  (0,x,y,z,+,-) on the GPUs
#endif
  complex(dp), dimension(:,:,:), allocatable :: pauli_orb
  !! Pauli matrices in spin and orbital space (x,y,z,+,-)
  complex(dp), dimension(:,:,:), allocatable :: pauli_dorb
  !! Pauli matrices in spin and orbital space (x,y,z,+,-) (only non-zero on d-orbitals)

  complex(dp), dimension(4,4) :: StoC = reshape([cmplx(0.0_dp,0._dp,dp), cmplx(0.5_dp,0._dp,dp),  cmplx(0.0_dp,-0.5_dp,dp),  cmplx( 0.0_dp,0._dp,dp), &
                                                 cmplx(1.0_dp,0._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx(0.0_dp, 0.0_dp,dp),  cmplx( 0.5_dp,0._dp,dp), &
                                                 cmplx(1.0_dp,0._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx(0.0_dp, 0.0_dp,dp),  cmplx(-0.5_dp,0._dp,dp), &
                                                 cmplx(0.0_dp,0._dp,dp), cmplx(0.5_dp,0._dp,dp),  cmplx(0.0_dp, 0.5_dp,dp),  cmplx( 0.0_dp,0._dp,dp)], [4,4])
  !! Transformation matrix spin (+,up,down,-) to cartesian (0,x,y,z)

  complex(dp), dimension(4,4) :: CtoS = reshape([cmplx(0.0_dp,0._dp,dp), cmplx(0.5_dp,0._dp,dp),  cmplx( 0.5_dp, 0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp), &
                                                 cmplx(1.0_dp,0._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp),  cmplx( 1.0_dp, 0._dp,dp), &
                                                 cmplx(0.0_dp,1._dp,dp), cmplx(0.0_dp,0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp),  cmplx( 0.0_dp,-1._dp,dp), &
                                                 cmplx(0.0_dp,0._dp,dp), cmplx(1.0_dp,0._dp,dp),  cmplx(-1.0_dp, 0._dp,dp),  cmplx( 0.0_dp, 0._dp,dp)], [4,4])
  !! Transformation matrix cartesian (0,x,y,z) to spin (+,up,down,-)
contains
  subroutine allocate_constants(nOrb)
    implicit none
    integer, intent(in) :: nOrb

    allocate(ident_norb(nOrb,nOrb),ident_norb2(2*nOrb,2*nOrb),ident_dorb(2*nOrb,2*nOrb),pauli_orb(5,2*nOrb,2*nOrb),pauli_dorb(5,2*nOrb,2*nOrb))

  end subroutine allocate_constants


  subroutine deallocate_constants()
    implicit none

    deallocate(ident_norb,ident_norb2,ident_dorb,pauli_orb,pauli_dorb)

  end subroutine deallocate_constants

  subroutine define_constants(s)
    use mod_System, only: System_type
    implicit none
    type(System_type), intent(in) :: s
    integer :: i,mu,nu,mud,nud

    ident_norb  = cZero
    ident_norb2 = cZero
    do i=1,s%nOrb
      ident_norb(i,i) = cOne
    end do
    do i=1,s%nOrb2
      ident_norb2(i,i) = cOne
    end do

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

! Pauli matrices in spin and orbital space
    pauli_dorb = cZero
    pauli_orb = cZero
    ident_dorb = cZero
    do mu = 1,s%nOrb
      nu = mu+s%nOrb

      ! pauli matrix x
      pauli_orb(1,mu,nu) = cOne
      pauli_orb(1,nu,mu) = cOne
      ! pauli matrix y
      pauli_orb(2,mu,nu) = -cI
      pauli_orb(2,nu,mu) = cI
      ! pauli matrix z
      pauli_orb(3,mu,mu) = cOne
      pauli_orb(3,nu,nu) = -cOne
    end do

    do mu = 1,s%ndOrb
      mud = s%dOrbs(mu)
      nud = mud+s%nOrb

      ! Identity
      ident_dorb(mud,mud)   = cOne
      ident_dorb(nud,nud)   = cOne
      ! pauli matrix x
      pauli_dorb(1,mud,nud) = cOne
      pauli_dorb(1,nud,mud) = cOne
      ! pauli matrix y
      pauli_dorb(2,mud,nud) = -cI
      pauli_dorb(2,nud,mud) = cI
      ! pauli matrix z
      pauli_dorb(3,mud,mud) = cOne
      pauli_dorb(3,nud,nud) = -cOne
    end do
    pauli_orb(4,:,:) = pauli_orb(1,:,:) + cI*pauli_orb(2,:,:)
    pauli_orb(5,:,:) = pauli_orb(1,:,:) - cI*pauli_orb(2,:,:)
    pauli_dorb(4,:,:) = pauli_dorb(1,:,:) + cI*pauli_dorb(2,:,:)
    pauli_dorb(5,:,:) = pauli_dorb(1,:,:) - cI*pauli_dorb(2,:,:)

    levi_civita(1,1,:) = [ 0._dp , 0._dp , 0._dp ]
    levi_civita(1,2,:) = [ 0._dp , 0._dp , 1._dp ]
    levi_civita(1,3,:) = [ 0._dp ,-1._dp , 0._dp ]
    levi_civita(2,1,:) = [ 0._dp , 0._dp ,-1._dp ]
    levi_civita(2,2,:) = [ 0._dp , 0._dp , 0._dp ]
    levi_civita(2,3,:) = [ 1._dp , 0._dp , 0._dp ]
    levi_civita(3,1,:) = [ 0._dp , 1._dp , 0._dp ]
    levi_civita(3,2,:) = [-1._dp , 0._dp , 0._dp ]
    levi_civita(3,3,:) = [ 0._dp , 0._dp , 0._dp ]

  end subroutine define_constants

  integer function delta(i,j)
    implicit none
    integer :: i,j
    delta = 0
    if(i==j) delta = 1
  end function delta
end module mod_constants
