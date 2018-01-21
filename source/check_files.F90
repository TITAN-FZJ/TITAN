! Check if files exist
subroutine check_files()
  use mod_parameters, only: itype,npt1,output
  use mod_susceptibilities, only: open_chi_files, close_chi_files, &
                                  open_dc_chi_files, close_dc_chi_files
  use mod_disturbances, only: open_disturbance_files, close_disturbance_files, &
                              open_dc_disturbance_files, close_dc_disturbance_files
  use mod_beff, only: open_beff_files, close_beff_files, &
                      open_dc_beff_files, close_dc_beff_files
  use mod_torques, only: open_torque_files, close_torque_files, &
                         open_dc_torque_files, close_dc_torque_files
  use mod_magnet, only: hw_count
  use mod_LDOS, only: openLDOSFiles,closeLDOSFiles
  !use mod_currents, only: openclose_currents_files, openclose_dc_currents_files !TODO: Re-Include
  !use mod_sha, only: openclose_sha_files, openclose_dc_sha_files !TODO: Re-Include
  implicit none
  integer :: count
  write(output%unit,"('[check_files] Checking if files exist...')")
  select case (itype)
  case (2)
    call openLDOSFiles()
    call closeLDOSFiles()
  case (7)
    call open_chi_files()
    call close_chi_files()
  case (8)
    call open_chi_files()
    call close_chi_files()
    call open_disturbance_files()
    call close_disturbance_files()
    call open_beff_files()
    call close_beff_files()
    call open_torque_files()
    call close_torque_files()
    !call open_currents_files() !TODO: Re-Include
    !call close_currents_files() !TODO: Re-Include
    !call open_sha_files(1) !TODO: Re-Include
    !call close_sha_files(2) !TODO: Re-Include
  case (9)
    if(hw_count==1) then
      do count=1, npt1
        call open_dc_chi_files()
        call close_dc_chi_files()
        call open_dc_disturbance_files()
        call close_dc_disturbance_files()
        call open_dc_beff_files()
        call close_dc_beff_files()
        call open_dc_torque_files()
        call close_dc_torque_files()
        !call open_dc_currents_files(1) !TODO: Re-Include
        !call close_dc_currents_files(2) !TODO: Re-Include
        !call open_dc_sha_files(1) !TODO: Re-Include
        !call close_dc_sha_files(2) !TODO: Re-Include
      end do
    end if
  end select
  write(output%unit,"('[check_files] All files exist!')")

  return
end subroutine check_files
