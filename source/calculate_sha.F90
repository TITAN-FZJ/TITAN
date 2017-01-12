! This subroutine calculates the spin hall angle
! given by
! SHA = (2.d0*e/\hbar) J_S/J_C

subroutine calculate_sha()
  use mod_parameters
  use mod_constants, only: zero
  use mod_magnet
  use mod_progress
  use mod_currents
  use mod_lattice
  use mod_sha
  implicit none
  integer      :: i,j,neighbor

  sha_re = 0.d0
  sha_complex = zero
  sha_re_total = 0.d0
  sha_complex_total = zero
  do neighbor=n0sc1,n0sc2
    ! Summing charge currents flowing on longitudinal direction
    if(any(neighbor.eq.sha_longitudinal(1:longitudinal_neighbors))) then
      do i=1,Npl
        sha_re(1,i) = sha_re(1,i) + real(currents(1,neighbor,i))
        sha_complex(1,i) = sha_complex(1,i) + currents(1,neighbor,i)
      end do
      sha_re_total(1) = sha_re_total(1) + real(total_currents(1,neighbor))
      sha_complex_total(1) = sha_complex_total(1) + total_currents(1,neighbor)
    end if

    ! Summing spin currents flowing on transverse direction
    if(any(neighbor.eq.sha_transverse(1:transverse_neighbors))) then
      do j=2,4
        do i=1,Npl
          sha_re(j,i) = sha_re(j,i) + real(currents(j,neighbor,i))
          sha_complex(j,i) = sha_complex(j,i) + currents(j,neighbor,i)
        end do
        sha_re_total(j) = sha_re_total(j) + real(total_currents(j,neighbor))
        sha_complex_total(j) = sha_complex_total(j) + total_currents(j,neighbor)
      end do
    end if
  end do
  ! Calculating real and imaginary part of spin Hall angles
  do i=1,Npl
    sha_re(2:4,i) = 2.d0*sha_re(2:4,i)/sha_re(1,i)
    sha_complex(2:4,i) = 2.d0*sha_complex(2:4,i)/sha_complex(1,i)
  end do
  sha_re_total(2:4) = 2.d0*sha_re_total(2:4)/sha_re_total(1)
  sha_complex_total(2:4) = 2.d0*sha_complex_total(2:4)/sha_complex_total(1)

  return
end subroutine calculate_sha