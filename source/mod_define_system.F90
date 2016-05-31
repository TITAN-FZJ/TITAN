! This module chooses the parameters for a slab using set1 and set2
! read in the input file
! The first principles parameters are listed in the module
! mod_tight_binding.F90
!
! The Fermi level is determined by the parameter of the first layer
!
module mod_define_system
  use mod_f90_kind
  implicit none

contains

  subroutine define_system()
    use mod_parameters, only: Npl,mmlayer,set1,set2
    use mod_mpi_pars, only: myrank
  implicit none
  integer :: i,j,ios,set1g,set2g

  if((set1.eq.9).or.(set2.eq.9)) then
    open(unit=112358, file="system.dat", status='old', iostat=ios)
    if(myrank.eq.0) then
      if (ios.ne.0) stop "[define_system] File 'system.dat' does not exist. "
    end if
    read(unit=112358,fmt=*,iostat=ios) (mmlayer(i),i=1,Npl+2)
    if(myrank.eq.0) then
      if (ios.ne.0) stop "[define_system] Problem defining system from file 'system.dat'. "
    end if
!     write(*,"(20(2x,i0))") (mmlayer(i),i=1,Npl+2)
!     stop
    return
  end if

  set1g = (set1-1)*6
  if(Npl.eq.1) set2 = set1
  set2g = (set2-1)*6

  i = ceiling((Npl+2)/2.d0)
  do j=1,i
    if(j.le.6) then
      mmlayer(j) = j + set1g
    else !bulk
      mmlayer(j) = 6 + set1g
    end if
  end do
  do j=i+1,Npl+2
    mmlayer(j) = Npl + 3 - j + set2g
    if(mmlayer(j).gt.(6+set2g)) then ! bulk
      mmlayer(j) = 6 + set2g
    end if
  end do

!   write(*,"(20(2x,i0))") (mmlayer(i),i=1,Npl+2)
!   stop

  return
  end subroutine define_system
end module mod_define_system