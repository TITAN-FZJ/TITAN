module mod_constants
  !! Contains commonly used constants
  use mod_f90_kind
  implicit none
  complex(double), parameter :: zero=(0.d0,0.d0)
  !! Complex scalar zero
  complex(double), parameter :: zum=(1.d0,0.d0)
  !! Complex scalar one
  complex(double), parameter :: zi=(0.d0,1.d0)
  !! Complex scalar i
  real(double), parameter :: pi = 4.d0 * atan(1.d0)
  !! Pi
  real(double), parameter :: tpi = 2.d0 * pi
  !! 2*Pi
  real(double), parameter :: sq2 = sqrt(2.d0)
  !! Square root of pi
  real(double), parameter :: hsq2 = 0.5d0 * sq2
  !! (Square root of pi) / 2
  real(double), parameter :: sq3 = sqrt(3.d0)
  !! Square root of 3
  real(double)    :: levi_civita(3,3,3)
  !! Levi Civita Tensor
  complex(double) :: identorb18(18,18)
  !! Identity
  complex(double) :: identorb9(9,9)
  complex(double) :: pauli_orb(3,18,18)
  complex(double) :: pauli_dorb(3,18,18)

contains
  subroutine define_constants()
    implicit none
    integer :: i,mu,nu

    identorb9  = zero
    identorb18 = zero
    do i=1,9
      identorb9(i,i) = zum
    end do
    do i=1,18
      identorb18(i,i) = zum
    end do

! Pauli matrices in spin and orbital space
    pauli_dorb = zero
    do mu = 1,9
      nu = mu+9

      ! pauli matrix x
      pauli_orb(1,mu,nu) = zum
      pauli_orb(1,nu,mu) = zum
      ! pauli matrix y
      pauli_orb(2,mu,nu) = -zi
      pauli_orb(2,nu,mu) = zi
      ! pauli matrix z
      pauli_orb(3,mu,mu) = zum
      pauli_orb(3,nu,nu) = -zum

      if (mu<5) cycle     ! Pauli matrices for d orbitals only
      ! pauli matrix x
      pauli_dorb(1,mu,nu) = zum
      pauli_dorb(1,nu,mu) = zum
      ! pauli matrix y
      pauli_dorb(2,mu,nu) = -zi
      pauli_dorb(2,nu,mu) = zi
      ! pauli matrix z
      pauli_dorb(3,mu,mu) = zum
      pauli_dorb(3,nu,nu) = -zum
    end do

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
!     pauli  = zero
!     paulid = zero
!     do mu = 1,9
!       nu = mu+9
!       ! identity
!       pauli(1,mu,mu) = zum
!       pauli(1,nu,nu) = zum
!       if (mu<5) cycle     ! Pauli matrices for d orbitals only
!       ! paulid matrix x
!       pauli(2,mu,nu) = zum
!       pauli(2,nu,mu) = zum
!       ! paulid matrix y
!       pauli(3,mu,nu) = -zi
!       pauli(3,nu,mu) = zi
!       ! paulid matrix z
!       pauli(4,mu,mu) = zum
!       pauli(4,nu,nu) = -zum

!       ! identity
!       paulid(1,mu,mu) = zum
!       paulid(1,nu,nu) = zum
!     end do
!     paulid(2:4,:,:) = pauli(2:4,:,:)

    return
  end subroutine define_constants
end module mod_constants
