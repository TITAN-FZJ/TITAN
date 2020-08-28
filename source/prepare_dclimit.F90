! This subroutine sets up external magnetic fields and related loop
subroutine prepare_dclimit()
  use, intrinsic :: iso_fortran_env
  use mod_parameters
  use mod_mpi_pars, only: myrank,ierr
  implicit none
  real(dp) :: e

  ! Inicial checks
  if(.not.lfield) then
    if(myrank==0) write(outputunit,"('[prepare_dclimit] External Field is off! Calculation of dc-limit needs external field dependence!')")
    call MPI_Finalize(ierr)
    stop
  end if

  if(total_hw_npt1==1) then
    if(myrank==0) write(outputunit,"('[prepare_dclimit] dc-limit calculation needs variation of one field variable (abs, theta or phi)!')")
    call MPI_Finalize(ierr)
    stop
  end if

  ! Allocating variable to write value of the fields on file
  allocate(dc_fields(total_hw_npt1),dcprefix(nEner1))

  ! Prefix for filenames containing the frequency/energy
  do kount=1,nEner1
    e = emin + deltae*(kount-1)
    if(e<2.e-6_dp) then
      write(dcprefix(kount),fmt="('dc')")
    else
      write(dcprefix(kount),fmt="('hw=',es8.1,'_')") e
    end if
  end do

  if((hwa_npt1>1).and.(hwt_npt1==1).and.(hwp_npt1==1)) then
    dcfield_dependence = 1
    dc_header = "      hwa       ,"
    do hw_count=1,total_hw_npt1
      write(dc_fields(hw_count),fmt="(es16.9,2x)") hw_list(hw_count,1)
    end do
  else if((hwa_npt1==1).and.(hwt_npt1>1).and.(hwp_npt1==1)) then
    dcfield_dependence = 2
    dc_header = "      hwt       ,"
    do hw_count=1,total_hw_npt1
      write(dc_fields(hw_count),fmt="(es16.9,2x)") hw_list(hw_count,2)
    end do
  else if((hwa_npt1==1).and.(hwt_npt1==1).and.(hwp_npt1>1)) then
    dcfield_dependence = 3
    dc_header = "      hwp       ,"
    do hw_count=1,total_hw_npt1
      write(dc_fields(hw_count),fmt="(es16.9,2x)") hw_list(hw_count,3)
    end do
  else if((hwa_npt1>1).and.(hwt_npt1>1).and.(hwp_npt1==1)) then
    dcfield_dependence = 4
    dc_header = "      hwa       ,      hwt       ,"
    do hw_count=1,total_hw_npt1
      write(dc_fields(hw_count),fmt="(2(es16.9,2x))") hw_list(hw_count,1),hw_list(hw_count,2)
    end do
  else if((hwa_npt1>1).and.(hwt_npt1==1).and.(hwp_npt1>1)) then
    dcfield_dependence = 5
    dc_header = "      hwa       ,      hwp       ,"
    do hw_count=1,total_hw_npt1
      write(dc_fields(hw_count),fmt="(2(es16.9,2x))") hw_list(hw_count,1),hw_list(hw_count,3)
    end do
  else if((hwa_npt1==1).and.(hwt_npt1>1).and.(hwp_npt1>1)) then
    dcfield_dependence = 6
    dc_header = "      hwt       ,      hwp       ,"
    do hw_count=1,total_hw_npt1
      write(dc_fields(hw_count),fmt="(2(es16.9,2x))") hw_list(hw_count,2),hw_list(hw_count,3)
    end do
  else if((hwa_npt1>1).and.(hwt_npt1>1).and.(hwp_npt1>1)) then
    dcfield_dependence = 7
    dc_header = "      hwa       ,      hwt       ,      hwp       ,"
    do hw_count=1,total_hw_npt1
      write(dc_fields(hw_count),fmt="(3(es16.9,2x))") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
    end do
  end if

end subroutine prepare_dclimit
