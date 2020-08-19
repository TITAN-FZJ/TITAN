! This module chooses the parameters for a slab using set1 and set2
! read in the input file
! The first principles parameters are listed in the module
! mod_tight_binding.F90
!
! The Fermi level is determined by the parameter of the first layer
!
module mod_define_system
  implicit none

contains

  subroutine define_system()
  use mod_parameters
  use mod_mpi_pars
  use mod_System,   only: s => sys
  implicit none
  integer :: i,j,ios,set1g,set2g

  if((set1==9).or.(set2==9)) then
    open(unit=112358, file="system.dat", status='old', iostat=ios)
    if (ios/=0) then
      if(myrank==0) write(output%unit_loop,"('[define_system] File ""system.dat"" does not exist. ')")
      call MPI_Finalize(ierr)
      stop
    end if
    read(unit=112358,fmt=*,iostat=ios) (mmlayer(i),i=1,Npl+2)
    if (ios/=0) then
      if(myrank==0) write(output%unit_loop,"('[define_system] Problem defining system from file ""system.dat"".')")
      call MPI_Finalize(ierr)
      stop
    end if
!     write(output%unit_loop,"(20(2x,i0))") (mmlayer(i),i=1,Npl+2)
!     stop
      end if

  set1g = (set1-1)*6
  if(Npl_input==1) set2 = set1
  set2g = (set2-1)*6

  ! i = ceiling((Npl_input+2)/2._dp)
  ! do j=1,i
  !   if(j<=6) then
  !     mmlayer(j) = j + set1g
  !   else !bulk
  !     mmlayer(j) = 6 + set1g
  !   end if
  ! end do
  ! do j=i+1,Npl_input+2
  !   mmlayer(j) = Npl_input + 3 - j + set2g
  !   if(mmlayer(j)>(6+set2g)) then ! bulk
  !     mmlayer(j) = 6 + set2g
  !   end if
  ! end do

!  write(strSites,fmt="(i0,'Npl')") Npl_input

  ! do i=1,naddlayers
  !   mmlayer(Npl_input+i+1) = addlayers(i)
  !   write(strSites,fmt="(a,'_',i0)") trim(strSites),addlayers(i)
  ! end do

!   write(output%unit_loop,"(20(2x,i0))") (mmlayer(i),i=1,Npl+2)
!   stop

  end subroutine define_system
end module mod_define_system
