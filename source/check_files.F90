! Check if files exist
subroutine check_files()
  use mod_parameters, only: itype,npt1,outputunit
  use mod_susceptibilities, only: openclose_chi_files, openclose_dc_chi_files
  !use mod_disturbances, only: openclose_disturbance_files, openclose_dc_disturbance_files !TODO: Re-Include
  !use mod_currents, only: openclose_currents_files, openclose_dc_currents_files !TODO: Re-Include
  !use mod_beff, only: openclose_beff_files, openclose_dc_beff_files !TODO: Re-Include
  !use mod_torques, only: openclose_torque_files, openclose_dc_torque_files !TODO: Re-Include
  !use mod_sha, only: openclose_sha_files, openclose_dc_sha_files !TODO: Re-Include
  use mod_magnet, only: hw_count
  implicit none
  integer :: count
  write(outputunit,"('[check_files] Checking if files exist...')")
  select case (itype)
  case (7)
    call openclose_chi_files(1)
    call openclose_chi_files(2)
  case (8)
    call openclose_chi_files(1)
    call openclose_chi_files(2)
    !call openclose_disturbance_files(1) !TODO: Re-Include
    !call openclose_disturbance_files(2) !TODO: Re-Include
    !call openclose_currents_files(1) !TODO: Re-Include
    !call openclose_currents_files(2) !TODO: Re-Include
    !call openclose_beff_files(1) !TODO: Re-Include
    !call openclose_beff_files(2) !TODO: Re-Include
    !call openclose_torque_files(1) !TODO: Re-Include
    !call openclose_torque_files(2) !TODO: Re-Include
    !call openclose_sha_files(1) !TODO: Re-Include
    !call openclose_sha_files(2) !TODO: Re-Include
  case (9)
    if(hw_count==1) then
      do count=1, npt1
        call openclose_dc_chi_files(1, count)
        call openclose_dc_chi_files(2, count)
        !call openclose_dc_disturbance_files(1) !TODO: Re-Include
        !call openclose_dc_disturbance_files(2) !TODO: Re-Include
        !call openclose_dc_currents_files(1) !TODO: Re-Include
        !call openclose_dc_currents_files(2) !TODO: Re-Include
        !call openclose_dc_beff_files(1) !TODO: Re-Include
        !call openclose_dc_beff_files(2) !TODO: Re-Include
        !call openclose_dc_torque_files(1) !TODO: Re-Include
        !call openclose_dc_torque_files(2) !TODO: Re-Include
        !call openclose_dc_sha_files(1) !TODO: Re-Include
        !call openclose_dc_sha_files(2) !TODO: Re-Include
      end do
    end if
  end select
  write(outputunit,"('[check_files] All files exist!')")

  return
end subroutine check_files
