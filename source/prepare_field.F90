! This subroutine sets up external magnetic fields and related loop
subroutine prepare_field()
  use mod_parameters
  implicit none
  integer       :: i

  if(lfield) then
    if (hwa_npts.eq.0) then
      hwa_f = hwa_i
      hwa_npts = 1
    end if
    if((abs(hwa_i)+abs(hwa_f)).gt.1.d-8) then
!         if(myrank.eq.0) then
!           write(outputunit,*) "hwa_i = ", hwa_i
!           write(outputunit,*) "hwa_f = ", hwa_f
!           write(outputunit,*) "hwa_npts = ", hwa_npts
!           write(outputunit,*) "hwa_npt1 = ", hwa_npt1
!           write(outputunit,*) "hwt_i = ", hwt_i
!           write(outputunit,*) "hwt_f = ", hwt_f
!           write(outputunit,*) "hwt_npts = ", hwt_npts
!           write(outputunit,*) "hwt_npt1 = ", hwt_npt1
!           write(outputunit,*) "hwp_i = ", hwp_i
!           write(outputunit,*) "hwp_f = ", hwp_f
!           write(outputunit,*) "hwp_npts = ", hwp_npts
!           write(outputunit,*) "hwp_npt1 = ", hwp_npt1
!         end if
      if (hwt_npts.eq.0) then
        hwt_f = hwt_i
        hwt_npts = 1
      end if
      if (hwp_npts.eq.0) then
        hwp_f = hwp_i
        hwp_npts = 1
      end if
      ! External field angular loops steps
      hwa_s = (hwa_f - hwa_i)/hwa_npts
      if(abs(hwa_s).le.1.d-10) hwa_npt1 = 1
      hwt_s = (hwt_f - hwt_i)/hwt_npts
      if(abs(hwt_s).le.1.d-10) hwt_npt1 = 1
      hwp_s = (hwp_f - hwp_i)/hwp_npts
      if(abs(hwp_s).le.1.d-10) hwp_npt1 = 1
    else ! hwa_i and hwa_f = 0
      ! Cartesian coordinates on spin system of reference
      hwa_i   = sqrt(hwx**2+hwy**2+hwz**2)
      hwa_f   = hwa_i
      if(abs(hwa_i).lt.1.d-8) then
        lfield = .false.
      else
        hwt_i    = acos(hwz/hwa_i)
        hwt_s    = 0.d0
        hwt_npt1 = 1
        hwp_i    = atan2(hwy,hwx)
        hwp_s    = 0.d0
        hwp_npt1 = 1
      end if
    end if
  else ! lfield
    hwa_i    = 0.d0
    hwa_f    = 0.d0
    hwa_s    = 0.d0
    hwa_npt1 = 1
    hwt_i    = 0.d0
    hwt_f    = 0.d0
    hwt_s    = 0.d0
    hwt_npt1 = 1
    hwp_i    = 0.d0
    hwp_f    = 0.d0
    hwp_s    = 0.d0
    hwp_npt1 = 1
  end if

  ! Total number of points in the loops
  total_hw_npt1 = hwa_npt1*hwt_npt1*hwp_npt1
  if(total_hw_npt1.eq.1) skip_steps_hw = 0
  allocate(hw_list(total_hw_npt1,3))

  ! Creating list of magnetic fields (in spherical coordinates)
  i = 1
  hw_intensity: do hwa_count=1,hwa_npt1
    hwa = hwa_i + (hwa_count-1)*hwa_s
    hw_theta: do hwt_count=1,hwt_npt1
      hwt = hwt_i + (hwt_count-1)*hwt_s
      hw_phi: do hwp_count=1,hwp_npt1
        hwp = hwp_i + (hwp_count-1)*hwp_s

        hw_list(i,:) = [ hwa , hwt , hwp ]
        i = i+1

      end do hw_phi
    end do hw_theta
  end do hw_intensity

  return
end subroutine prepare_field