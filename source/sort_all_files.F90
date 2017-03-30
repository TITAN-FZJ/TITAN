! Sort all files
subroutine sort_all_files()
  use mod_parameters, only: itype,hw_count,count,npt1,outputunit_loop
  use mod_susceptibilities
  use mod_disturbances
  use mod_currents
  use mod_beff
  use mod_torques
  use mod_sha
  implicit none

  write(outputunit_loop,"('[sort_all_files] Sorting files... ',$)")
  select case (itype)
  case (7)
    ! SORTING SUSCEPTIBILITIES
    call sort_susceptibilities()
  case (8)
    ! SORTING SUSCEPTIBILITIES
    call sort_susceptibilities()
    ! SORTING DISTURBANCES
    call sort_disturbances()
    ! SORTING CURRENTS
    call sort_currents()
    ! SORTING EFFECTIVE FIELDS
    call sort_beff()
    ! SORTING TORQUES
    call sort_torques()
    ! SORTING SHA
    call sort_sha()
  case (9)
    if(hw_count.eq.1) then
      energy_loop: do count=1,npt1
        ! SORTING SUSCEPTIBILITIES
        call sort_susceptibilities()
        ! SORTING DISTURBANCES
        call sort_disturbances()
        ! SORTING CURRENTS
        call sort_currents()
        ! SORTING EFFECTIVE FIELDS
        call sort_beff()
        ! SORTING TORQUES
        call sort_torques()
        ! SORTING SHA
        call sort_sha()
      end do energy_loop
    end if
  end select
  write(outputunit_loop,"('Done! All files sorted!')")

  return
end subroutine sort_all_files