! Sort all files
subroutine sort_all_files()
  use mod_parameters, only: itype,count,npt1,outputunit_loop
  use mod_susceptibilities
  !use mod_disturbances!TODO: Re-Include
  !use mod_currents!TODO: Re-Include
  !use mod_beff!TODO: Re-Include
  !use mod_torques!TODO: Re-Include
  !use mod_sha !TODO: Re-Include
  use mod_magnet, only: hw_count
  implicit none

  write(outputunit_loop,"('[sort_all_files] Sorting files... ')", advance='no')
  select case (itype)
  case (7)
    ! SORTING SUSCEPTIBILITIES
    call sort_susceptibilities()
  case (8)
    ! SORTING SUSCEPTIBILITIES
    call sort_susceptibilities()
    ! SORTING DISTURBANCES
    !call sort_disturbances() !TODO: Re-Include
    ! SORTING CURRENTS
    !call sort_currents() !TODO: Re-Include
    ! SORTING EFFECTIVE FIELDS
    !call sort_beff() !TODO: Re-Include
    ! SORTING TORQUES
    !call sort_torques() !TODO: Re-Include
    ! SORTING SHA
    !call sort_sha() !TODO: Re-Include
  case (9)
    if(hw_count==1) then
      energy_loop: do count=1,npt1
        ! SORTING SUSCEPTIBILITIES
        call sort_susceptibilities(count)
        ! SORTING DISTURBANCES
        !call sort_disturbances() !TODO: Re-Include
        ! SORTING CURRENTS
        !call sort_currents() !TODO: Re-Include
        ! SORTING EFFECTIVE FIELDS
        !call sort_beff() !TODO: Re-Include
        ! SORTING TORQUES
        !call sort_torques() !TODO: Re-Include
        ! SORTING SHA
        !call sort_sha() !TODO: Re-Include
      end do energy_loop
    end if
  end select
  write(outputunit_loop,"('Done! All files sorted!')")

  return
end subroutine sort_all_files
