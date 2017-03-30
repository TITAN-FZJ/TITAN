! Check if files exist
subroutine check_files()
  use mod_parameters, only: itype,hw_count,count,npt1,outputunit
  use mod_susceptibilities
  use mod_disturbances
  use mod_currents
  use mod_beff
  use mod_torques
  use mod_sha
  implicit none

  write(outputunit,"('[check_files] Checking if files exist...')")
  select case (itype)
  case (7)
    call openclose_chi_files(1)
    call openclose_chi_files(2)
  case (8)
    call openclose_chi_files(1)
    call openclose_chi_files(2)
    call openclose_disturbance_files(1)
    call openclose_disturbance_files(2)
    call openclose_currents_files(1)
    call openclose_currents_files(2)
    call openclose_beff_files(1)
    call openclose_beff_files(2)
    call openclose_torque_files(1)
    call openclose_torque_files(2)
    call openclose_sha_files(1)
    call openclose_sha_files(2)
  case (9)
    if(hw_count.eq.1) then
      energy_loop: do count=1,npt1
        call openclose_dc_chi_files(1)
        call openclose_dc_chi_files(2)
        call openclose_dc_disturbance_files(1)
        call openclose_dc_disturbance_files(2)
        call openclose_dc_currents_files(1)
        call openclose_dc_currents_files(2)
        call openclose_dc_beff_files(1)
        call openclose_dc_beff_files(2)
        call openclose_dc_torque_files(1)
        call openclose_dc_torque_files(2)
        call openclose_dc_sha_files(1)
        call openclose_dc_sha_files(2)
      end do energy_loop
    end if
  end select
  write(outputunit,"('[check_files] All files exist!')")

  return
end subroutine check_files