! Set energy (frequency) and wave vector loops
subroutine setLoops(s)
  use mod_f90_kind,   only: double
  use mod_parameters, only: itype, nEner, nEner1, emin, emax, deltae, nQvec, nQvec1, deltak, kdirection, band_points, band_cnt, bands, kpoints, partial_length
  use mod_system,     only: System
  use mod_tools,      only: vecDist
  implicit none
  type(System), intent(in)  :: s
  integer                   :: count,i,j
  real(double)              :: dir(3),total_length,total_length1
    interface
      subroutine read_band_points(kbands, a0, b1, b2, b3)
        use mod_f90_kind, only: double
        real(double), dimension(:,:), allocatable, intent(out) :: kbands
        real(double), dimension(3), intent(in) :: b1, b2, b3
        real(double),               intent(in) :: a0
      end subroutine
    end interface

  ! Energy loop step
  if(nEner < 1) then
    deltae = 0
  else
    deltae = (emax - emin)/nEner
  end if
  if(deltae<=1.d-14) nEner1 = 1

  ! Wave vector loop
  if((itype == 4).or.((itype >= 7).and.(itype <= 9))) then
    if(nQvec1==1) then
      deltak = 0.d0
      allocate(kpoints(3,1))

      if(band_cnt == 1) then
        call read_band_points(band_points,s%a0,s%b1,s%b2,s%b3)
        kpoints(:,1) = band_points(:,1)
      else
        kpoints(:,1) = [0.d0,0.d0,0.d0]
      end if

    else ! If a path is given
      ! Reading different points
      call read_band_points(band_points,s%a0,s%b1,s%b2,s%b3)

      ! Calculating total length
      total_length = 0.d0
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
      partial_length(1) = 0.d0
      do count = 1, nQvec1
        ! if( (deltak*count >= total_length1).and.(i < band_cnt - 1) ) then
        if( (deltak*count >= sum(partial_length(1:i+1))).and.(i < band_cnt - 1) ) then
          i = i + 1
          partial_length(i+1) = vecDist(band_points(1:3,i),band_points(1:3,i+1))
          j = 0
          dir = (band_points(1:3,i+1) - band_points(1:3,i))
          dir = dir / partial_length(i+1)
        end if
        kpoints(:,count) = band_points(:,i) + dir * j * deltak
        j = j + 1
      end do

      ! write(*,*) sum(partial_length),total_length
    end if ! nQvec
  end if ! itype

end subroutine setLoops