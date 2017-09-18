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

    identorb9  = cZero
    identorb18 = cZero
    do i=1,9
      identorb9(i,i) = cOne
    end do
    do i=1,18
      identorb18(i,i) = cOne
    end do

! Pauli matrices in spin and orbital space
    pauli_dorb = cZero
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

      if (mu<5) cycle     ! Pauli matrices for d orbitals only
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

    return
  end subroutine define_constants
end module mod_constants
