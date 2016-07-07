! Check if files exist
subroutine check_files()
  use mod_parameters, only: itype,lhfresponses,hw_count
  use mod_susceptibilities
  use mod_disturbances
  use mod_currents
  use mod_beff
  use mod_torques
  implicit none

  select case (itype)
  case (7)
    call openclose_chi_files(1)
    call openclose_chi_files(2)
  case (8)
    if(.not.lhfresponses) then
      call openclose_chi_files(1)
      call openclose_chi_files(2)
    end if
    call openclose_disturbance_files(1)
    call openclose_disturbance_files(2)
    call openclose_currents_files(1)
    call openclose_currents_files(2)
    call openclose_beff_files(1)
    call openclose_beff_files(2)
    call openclose_torque_files(1)
    call openclose_torque_files(2)
  case (9)
    if(hw_count.eq.1) then
      if(.not.lhfresponses) then
        call openclose_dc_chi_files(1)
        call openclose_dc_chi_files(2)
      end if
      call openclose_dc_disturbance_files(1)
      call openclose_dc_disturbance_files(2)
      call openclose_dc_currents_files(1)
      call openclose_dc_currents_files(2)
      call openclose_dc_beff_files(1)
      call openclose_dc_beff_files(2)
      call openclose_dc_torque_files(1)
      call openclose_dc_torque_files(2)
    end if
  end select

  return
end subroutine check_files