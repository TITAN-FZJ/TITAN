! Set energy (frequency) and wave vector loops
subroutine setLoops(s)
  use mod_f90_kind,   only: double
  use mod_parameters, only: itype, output, nEner, nEner1, emin, emax, deltae, nQvec, nQvec1, deltak, band_points, band_cnt, bands, kpoints, bsfile
  use mod_system,     only: System
  implicit none
  type(System), intent(in) :: s
  character(len=20)        :: kdirection
  integer                  :: count,i,j
  real(double)             :: dir(3),total_length

  ! Energy loop step
  if(nEner < 1) then
    deltae = 0
  else
    deltae = (emax - emin)/nEner
  end if
  if(deltae<=1.d-14) nEner1 = 1

  ! Wave vector loop
  if((itype == 4).or.(itype == 7)) then

    if(nQvec1==1) then
      deltak = 0.d0
      allocate(kpoints(3,1))

      if(band_cnt == 1) then
        call read_band_points(band_points, s%b1,s%b2,s%b3)
        kpoints(:,1) = band_points(:,1)
      else
        kpoints(:,1) = [0.d0,0.d0,0.d0]
      end if

    else ! If a path is given
      call read_band_points(band_points, s%b1,s%b2,s%b3)
      total_length = 0.d0
      do i = 1, band_cnt - 1
        total_length = total_length + sqrt(dot_product(band_points(:,i)-band_points(:,i+1), band_points(:,i)-band_points(:,i+1)))
      end do

      if(itype == 4) then
        write(kdirection,*) (trim(adjustl(bands(i))), i = 1,band_cnt)
        write(bsfile,"('./results/',a1,'SOC/',a,'/BS/bandstructure_kdir=',a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(adjustl(kdirection)),trim(output%info),trim(output%BField),trim(output%SOC)
        open (unit=666, file=bsfile,status='replace')
        write(unit=666, fmt=*) band_cnt, nQvec
        write(unit=666, fmt=*) s%Ef
      end if

      deltak = total_length / nQvec
      allocate(kpoints(3,nQvec1))

      total_length = 0.d0
      count = 0
      i = 0
      do count = 1, nQvec1
        if(deltak*count >= total_length .and. i < band_cnt - 1) then
          i = i + 1
          if(itype == 4) write(unit=666, fmt=*) bands(i), total_length
          total_length = total_length + sqrt(dot_product(band_points(1:3,i)-band_points(1:3,i+1), band_points(1:3,i)-band_points(1:3,i+1)))
          j = 0
          dir = (band_points(1:3,i+1) - band_points(1:3,i))
          dir = dir / sqrt(dot_product(dir,dir))
        end if
        kpoints(:,count) = band_points(:,i) + dir * j * deltak
        j = j + 1
      end do
      i = i+1
      if(itype == 4) then
        write(unit=666, fmt=*) bands(i), total_length
        close(unit=666)
      end if
    end if ! nQvec

  end if ! itype

end subroutine setLoops