! Create files with headers
subroutine create_files()
  use mod_parameters, only: itype,output,count,npt1
  use mod_susceptibilities
  use mod_disturbances !TODO: Re-Include
  use mod_beff !TODO: Re-Include
  use mod_torques !TODO: Re-Include
  !use mod_currents !TODO: Re-Include
  !use mod_sha !TODO: Re-Include
  implicit none

  select case (itype)
  case (7)
    call create_chi_files()
    write(output%unit,"('[create_files] Susceptibilities files created/overwritten!')")
  case (8)
    call create_chi_files()
    call create_disturbance_files() !TODO: Re-Include
    call create_beff_files() !TODO: Re-Include
    call create_torque_files() !TODO: Re-Include
    !call create_currents_files() !TODO: Re-Include
    !call create_sha_files() !TODO: Re-Include
    write(output%unit,"('[create_files] Susceptibilities, disturbances, current, SHA, effective field, torque and DC voltage files created/overwritten!')")
  case (9)
    do count = 1, npt1
      call create_dc_chi_files()
      call create_dc_disturbance_files() !TODO: Re-Include
      call create_dc_beff_files() !TODO: Re-Include
      call create_dc_torque_files() !TODO: Re-Include
      !call create_dc_currents_files() !TODO: Re-Include
      !call create_dc_sha_files() !TODO: Re-Include
    end do
    write(output%unit,"('[create_files] Susceptibilities, disturbances, current, SHA, effective field, torque and DC voltage files created/overwritten!')")
  case default
    write(output%unit,"('[create_files] No files to create for selected option! (itype = ',i0,')')") itype
  end select

end subroutine create_files
