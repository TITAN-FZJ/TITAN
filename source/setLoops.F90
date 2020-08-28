! Set energy (frequency) and wave vector loops
subroutine setLoops(s)
  use mod_kind, only: dp
  use mod_parameters, only: itype, nEner, nEner1, emin, emax, deltae, nQvec, nQvec1, deltak, kdirection, band_points, band_cnt, bands, kpoints, partial_length
  use mod_system,     only: System_type
  use mod_tools,      only: vecDist
  implicit none
  type(System_type), intent(in)  :: s
  integer                   :: kount,i,j
  real(dp)              :: dir(3),total_length
    interface
      subroutine read_band_points(kbands, a0, b1, b2, b3)
        use mod_kind, only: dp
        real(dp), dimension(:,:), allocatable, intent(out) :: kbands
        real(dp), dimension(3), intent(in) :: b1, b2, b3
        real(dp),               intent(in) :: a0
      end subroutine
    end interface

  ! Energy loop step
  if(nEner < 1) then
    deltae = 0
  else
    deltae = (emax - emin)/nEner
  end if
  if(deltae<=1.e-15_dp) nEner1 = 1

  ! Wave vector loop
  if((itype == 4).or.((itype >= 7).and.(itype <= 9))) then
    if(nQvec1==1) then
      deltak = 0._dp
      allocate(kpoints(3,1))

      if(band_cnt == 1) then
        call read_band_points(band_points,s%a0,s%b1,s%b2,s%b3)
        kpoints(:,1) = band_points(:,1)
      else
        kpoints(:,1) = [0._dp,0._dp,0._dp]
      end if

    else ! If a path is given
      ! Reading different points
      call read_band_points(band_points,s%a0,s%b1,s%b2,s%b3)

      ! Calculating total length
      total_length = 0._dp
      do i = 1, band_cnt - 1
        total_length = total_length + vecDist(band_points(:,i),band_points(:,i+1))
      end do

      ! Writing path in symmetry points
      write(kdirection,*) (trim(adjustl(bands(i))), i = 1,band_cnt)

      ! Calculating step size between each k-point
      deltak = total_length / nQvec
      allocate( kpoints(3,nQvec1), partial_length(band_cnt) )

      ! Obtaining the k-points in the chosen path and the partial lenghts
      i = 0
      partial_length(1) = 0._dp
      do kount = 1, nQvec1
        ! if( (deltak*kount >= total_length1).and.(i < band_cnt - 1) ) then
        if( (deltak*kount >= sum(partial_length(1:i+1))).and.(i < band_cnt - 1) ) then
          i = i + 1
          partial_length(i+1) = vecDist(band_points(1:3,i),band_points(1:3,i+1))
          j = 0
          dir = (band_points(1:3,i+1) - band_points(1:3,i))
          dir = dir / partial_length(i+1)
        end if
        kpoints(:,kount) = band_points(:,i) + dir * j * deltak
        j = j + 1
      end do

      ! write(*,*) sum(partial_length),total_length
    end if ! nQvec
  end if ! itype

end subroutine setLoops



! Deallocate loop variables
subroutine deallocateLoops()
  use mod_parameters, only: itype, kpoints, partial_length
  use mod_magnet,     only: hw_list
  implicit none

  if(allocated(hw_list)) deallocate(hw_list)

  ! Wave vector loop
  if((itype == 4).or.((itype >= 7).and.(itype <= 9))) then
    if(allocated(kpoints)) deallocate(kpoints)
    if(allocated(partial_length)) deallocate(partial_length)
  end if ! itype

end subroutine deallocateLoops



