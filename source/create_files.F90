! Create files with headers
subroutine create_files()
  use mod_parameters, only: itype,outputunit,count,npt1
  use mod_susceptibilities
  use mod_disturbances
  use mod_currents
  use mod_beff
  use mod_torques
  use mod_sha
  implicit none

  select case (itype)
  case (7)
    call openclose_chi_files(0)
    write(outputunit,"('[create_files] Susceptibilities files created/overwritten!')")
  case (8)
    call openclose_chi_files(0)
    call openclose_disturbance_files(0)
    call openclose_currents_files(0)
    call openclose_beff_files(0)
    call openclose_torque_files(0)
    call openclose_sha_files(0)
    write(outputunit,"('[create_files] Susceptibilities, disturbances, current, SHA, effective field, torque and DC voltage files created/overwritten!')")
  case (9)
    energy_loop: do count=1,npt1
      call openclose_dc_chi_files(0)
      call openclose_dc_disturbance_files(0)
      call openclose_dc_currents_files(0)
      call openclose_dc_beff_files(0)
      call openclose_dc_torque_files(0)
      call openclose_dc_sha_files(0)
    end do energy_loop
    write(outputunit,"('[create_files] Susceptibilities, disturbances, current, SHA, effective field, torque and DC voltage files created/overwritten!')")
  case default
    write(outputunit,"('[create_files] No files to create for selected option! (itype = ',i0,')')") itype
  end select

  return
end subroutine create_files