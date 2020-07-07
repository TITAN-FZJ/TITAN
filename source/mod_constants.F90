module mod_constants
  !! Contains commonly used constants
  use mod_f90_kind
  implicit none
  complex(double), parameter :: cZero=(0.d0,0.d0)
  !! Complex scalar cZero
  complex(double), parameter :: cOne=(1.d0,0.d0)
  !! Complex scalar one
  complex(double), parameter :: cI=(0.d0,1.d0)
  !! Complex scalar i
  real(double), parameter :: pi = 4.d0 * atan(1.d0)
  !! Pi
  real(double), parameter :: tpi = 2.d0 * pi
  !! 2*Pi
  real(double), parameter :: sq2 = sqrt(2.d0)
  !! Square root of 2
  real(double), parameter :: hsq2 = 0.5d0 * sq2
  !! (Square root of 2) / 2
  real(double), parameter :: sq3 = sqrt(3.d0)
  !! Square root of 3
  real(double), parameter :: deg2rad = pi/180
  !! degrees to radians
  real(double), parameter :: rad2deg = 180/pi
  !! radians to degrees
  real(double)    :: levi_civita(3,3,3)
  !! Levi Civita Tensor
  complex(double) :: ident_norb(9,9)
  !! Identity in orbital space
  complex(double) :: ident_norb2(18,18)
  !! Identity in spin and orbital space
  complex(double) :: ident_dorb(18,18)
  !! Identity in spin and orbital space (only non-zero on d-orbitals)
  complex(double) :: pauli_mat(2,2,0:5)
  !! Identity and pauli matrices  (0,x,y,z,+,-)
  complex(double) :: pauli_orb(5,18,18)
  !! Pauli matrices in spin and orbital space (x,y,z,+,-)
  complex(double) :: pauli_dorb(5,18,18)
  !! Pauli matrices in spin and orbital space (x,y,z,+,-) (only non-zero on d-orbitals)

  complex(double), dimension(4,4) :: StoC = reshape([cmplx(0.0d0,0.d0,double), cmplx(0.5d0,0.d0,double),  cmplx(0.0d0,-0.5d0,double),  cmplx( 0.0d0,0.d0,double), &
                                                     cmplx(1.0d0,0.d0,double), cmplx(0.0d0,0.d0,double),  cmplx(0.0d0, 0.0d0,double),  cmplx( 0.5d0,0.d0,double), &
                                                     cmplx(1.0d0,0.d0,double), cmplx(0.0d0,0.d0,double),  cmplx(0.0d0, 0.0d0,double),  cmplx(-0.5d0,0.d0,double), &
                                                     cmplx(0.0d0,0.d0,double), cmplx(0.5d0,0.d0,double),  cmplx(0.0d0, 0.5d0,double),  cmplx( 0.0d0,0.d0,double)], [4,4])
  !! Transformation matrix spin (+,up,down,-) to cartesian (0,x,y,z)

  complex(double), dimension(4,4) :: CtoS = reshape([cmplx(0.0d0,0.d0,double), cmplx(0.5d0,0.d0,double),  cmplx( 0.5d0, 0.d0,double),  cmplx( 0.0d0, 0.d0,double), &
                                                     cmplx(1.0d0,0.d0,double), cmplx(0.0d0,0.d0,double),  cmplx( 0.0d0, 0.d0,double),  cmplx( 1.0d0, 0.d0,double), &
                                                     cmplx(0.0d0,1.d0,double), cmplx(0.0d0,0.d0,double),  cmplx( 0.0d0, 0.d0,double),  cmplx( 0.0d0,-1.d0,double), &
                                                     cmplx(0.0d0,0.d0,double), cmplx(1.0d0,0.d0,double),  cmplx(-1.0d0, 0.d0,double),  cmplx( 0.0d0, 0.d0,double)], [4,4])
  !! Transformation matrix cartesian (0,x,y,z) to spin (+,up,down,-)
contains
  subroutine define_constants()
    implicit none
    integer :: i,mu,nu

    ident_norb  = cZero
    ident_norb2 = cZero
    do i=1,9
      ident_norb(i,i) = cOne
    end do
    do i=1,18
      ident_norb2(i,i) = cOne
    end do

    ! All four Pauli matrices in spin space
    pauli_mat(:,1,0) = [cOne, cZero]
    pauli_mat(:,2,0) = [cZero, cOne]
    pauli_mat(:,1,1) = [cZero,cOne]
    pauli_mat(:,2,1) = [cOne, cZero]
    pauli_mat(:,1,2) = [cZero,cI]
    pauli_mat(:,2,2) = [-cI,cZero]
    pauli_mat(:,1,3) = [cOne,cZero]
    pauli_mat(:,2,3) = [cZero,-cOne]

    pauli_mat(:,:,4) = pauli_mat(:,:,1) + cI*pauli_mat(:,:,2)
    pauli_mat(:,:,5) = pauli_mat(:,:,1) - cI*pauli_mat(:,:,2)


! Pauli matrices in spin and orbital space
    pauli_dorb = cZero
    ident_dorb = cZero
    do mu = 1,9
      nu = mu+9

      ! pauli matrix x
      pauli_orb(1,mu,nu) = cOne
      pauli_orb(1,nu,mu) = cOne
      ! pauli matrix y
      pauli_orb(2,mu,nu) = -cI
      pauli_orb(2,nu,mu) = cI
      ! pauli matrix z
      pauli_orb(3,mu,mu) = cOne
      pauli_orb(3,nu,nu) = -cOne

      if (mu<5) cycle     ! Identity and Pauli matrices for d orbitals only
      ! Identity
      ident_dorb(mu,mu)   = cOne
      ident_dorb(nu,nu)   = cOne
      ! pauli matrix x
      pauli_dorb(1,mu,nu) = cOne
      pauli_dorb(1,nu,mu) = cOne
      ! pauli matrix y
      pauli_dorb(2,mu,nu) = -cI
      pauli_dorb(2,nu,mu) = cI
      ! pauli matrix z
      pauli_dorb(3,mu,mu) = cOne
      pauli_dorb(3,nu,nu) = -cOne
    end do
    pauli_orb(4,:,:) = pauli_orb(1,:,:) + cI*pauli_orb(2,:,:)
    pauli_orb(5,:,:) = pauli_orb(1,:,:) - cI*pauli_orb(2,:,:)
    pauli_dorb(4,:,:) = pauli_dorb(1,:,:) + cI*pauli_dorb(2,:,:)
    pauli_dorb(5,:,:) = pauli_dorb(1,:,:) - cI*pauli_dorb(2,:,:)

    levi_civita(1,1,:) = [ 0.d0 , 0.d0 , 0.d0 ]
    levi_civita(1,2,:) = [ 0.d0 , 0.d0 , 1.d0 ]
    levi_civita(1,3,:) = [ 0.d0 ,-1.d0 , 0.d0 ]
    levi_civita(2,1,:) = [ 0.d0 , 0.d0 ,-1.d0 ]
    levi_civita(2,2,:) = [ 0.d0 , 0.d0 , 0.d0 ]
    levi_civita(2,3,:) = [ 1.d0 , 0.d0 , 0.d0 ]
    levi_civita(3,1,:) = [ 0.d0 , 1.d0 , 0.d0 ]
    levi_civita(3,2,:) = [-1.d0 , 0.d0 , 0.d0 ]
    levi_civita(3,3,:) = [ 0.d0 , 0.d0 , 0.d0 ]

!   ! Pauli matrices in spin and orbital space
!     pauli  = cZero
!     paulid = cZero
!     do mu = 1,9
!       nu = mu+9
!       ! identity
!       pauli(1,mu,mu) = cOne
!       pauli(1,nu,nu) = cOne
!       if (mu<5) cycle     ! Pauli matrices for d orbitals only
!       ! paulid matrix x
!       pauli(2,mu,nu) = cOne
!       pauli(2,nu,mu) = cOne
!       ! paulid matrix y
!       pauli(3,mu,nu) = -cI
!       pauli(3,nu,mu) = cI
!       ! paulid matrix z
!       pauli(4,mu,mu) = cOne
!       pauli(4,nu,nu) = -cOne

!       ! identity
!       paulid(1,mu,mu) = cOne
!       paulid(1,nu,nu) = cOne
!     end do
!     paulid(2:4,:,:) = pauli(2:4,:,:)

  end subroutine define_constants

  integer function delta(i,j)
    implicit none
    integer :: i,j
    delta = 0
    if(i==j) delta = 1
  end function delta
end module mod_constants
