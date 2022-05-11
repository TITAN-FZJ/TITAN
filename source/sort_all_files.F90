subroutine sort_all_files()
  !! This subroutine is a wrapper for the sorting calls depending on the itype
  use mod_parameters,       only: itype,kount,nEner1,output
  use mod_magnet,           only: hw_count
  use mod_susceptibilities, only: sort_susceptibilities
  use mod_disturbances,     only: sort_disturbances
  use mod_beff,             only: sort_beff
  use mod_torques,          only: sort_torques
  use mod_LDOS,             only: sortLDOS
  !use mod_currents!TODO: Re-Include
  !use mod_sha !TODO: Re-Include
  implicit none

  write(output%unit_loop,"('[sort_all_files] Sorting files... ')", advance='no')
  select case (itype)
  case (2)
    ! SORTING LDOS
    call sortLDOS()
  case (7)
    ! SORTING SUSCEPTIBILITIES
    call sort_susceptibilities()
  case (8)
    ! SORTING SUSCEPTIBILITIES
    call sort_susceptibilities()
    ! SORTING DISTURBANCES
    call sort_disturbances()
    ! SORTING EFFECTIVE FIELDS
    call sort_beff()
    ! SORTING TORQUES
    call sort_torques()
    ! SORTING CURRENTS
    !call sort_currents() !TODO: Re-Include
    ! SORTING SHA
    !call sort_sha() !TODO: Re-Include
  case (9)
    if(hw_count==1) then
      do kount = 1, nEner1
        ! SORTING SUSCEPTIBILITIES
        call sort_susceptibilities()
        ! SORTING DISTURBANCES
        call sort_disturbances()
        ! SORTING EFFECTIVE FIELDS
        call sort_beff()
        ! SORTING TORQUES
        call sort_torques()
        ! SORTING CURRENTS
        !call sort_currents() !TODO: Re-Include
        ! SORTING SHA
        !call sort_sha() !TODO: Re-Include
      end do
    end if
  end select
  write(output%unit_loop,"('Done! All files sorted!')")

end subroutine sort_all_files
