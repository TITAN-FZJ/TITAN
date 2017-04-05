! This subroutine may be used to do tests on the program
subroutine debugging()
  use mod_f90_kind
  use mod_constants
  use mod_lattice
  use mod_parameters
  use mod_magnet
  use mod_tight_binding
  use mod_prefactors
  use mod_generate_kpoints
  use mod_generate_epoints
  use mod_mpi_pars
  use mod_progress
  use mod_diamagnetic_current
  implicit none

  if(myrank==0) write(outputunit,"('[debugging] Starting to debug...')")

!   integer   :: i,i0,i1,j,j0,j1
!   real(double)  :: kp(3)
!   complex(double),dimension((Npl+2)*18,(Npl+2)*18)  :: hk,test,test2
!   complex(double),dimension(Npl+2,Npl+2,18,18)  :: gf

!   write(*,*) 'Hamiltonian:'
  mz    = 2.d0
  mp  = zero
  mm  = conjg(mp)
  ! Variables used in the hamiltonian
  eps1 = 0.d0
  hdel  = 0.d0
  hdelp = zero
  hdelm = zero
  hdel(1:Npl)   = 0.5d0*U*mz
  hdelp(1:Npl)  = 0.5d0*U*mp
  hdelm(1:Npl)  = 0.5d0*U*mm
!   lambda=0.d0
  ! call calculate_idia()

!   write(*,*) 'Real space hoppings:'
!   write(*,*) '*********************'
!   write(*,*) sum(abs(t00(1,:,:,:))),sum(t00(1,:,:,:))
!   write(*,*) sum(abs(t00(7,:,:,:))),sum(t00(7,:,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(t00(2,:,:,:))),sum(t00(2,:,:,:))
!   write(*,*) sum(abs(t00(6,:,:,:))),sum(t00(6,:,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(t00(3,:,:,:))),sum(t00(3,:,:,:))
!   write(*,*) sum(abs(t00(5,:,:,:))),sum(t00(5,:,:,:))
!   write(*,*) '*********************'
!   write(*,*) sum(abs(t01(1,:,:,:))),sum(t01(1,:,:,:))
!   write(*,*) sum(abs(t01(6,:,:,:))),sum(t01(6,:,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(t01(2,:,:,:))),sum(t01(2,:,:,:))
!   write(*,*) sum(abs(t01(5,:,:,:))),sum(t01(5,:,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(t01(3,:,:,:))),sum(t01(3,:,:,:))
!   write(*,*) sum(abs(t01(4,:,:,:))),sum(t01(4,:,:,:))
!   write(*,*) '*********************'

!   write(*,*) 'Occupations:'
!   write(*,*) '*********************'
!   write(*,*) npart0(1,1),npart0(1,2),npart0(1,3)
!   write(*,*) npart0(7,1),npart0(7,2),npart0(7,3)
!   write(*,*) '---------------------'
!   write(*,*) npart0(2,1),npart0(2,2),npart0(2,3)
!   write(*,*) npart0(6,1),npart0(6,2),npart0(6,3)
!   write(*,*) '---------------------'
!   write(*,*) npart0(3,1),npart0(3,2),npart0(3,3)
!   write(*,*) npart0(5,1),npart0(5,2),npart0(5,3)
!   write(*,*) '*********************'

!   write(*,*) 'SOC parameters:'
!   write(*,*) '*********************'
!   write(*,*) lambda(1)
!   write(*,*) lambda(7)
!   write(*,*) '---------------------'
!   write(*,*) lambda(2)
!   write(*,*) lambda(6)
!   write(*,*) '---------------------'
!   write(*,*) lambda(3)
!   write(*,*) lambda(5)
!   write(*,*) '*********************'


!   write(*,*) 'Hamiltonian:'
!   mz    = 0.5d0
!   mp  = zero
!   mm  = conjg(mp)
!   ! Variables used in the hamiltonian
!   eps1 = 0.d0
!   hdel  = 0.d0
!   hdelp = zero
!   hdelm = zero
!   hdel(2:Npl+1)   = 0.5d0*U*mz
!   hdelp(2:Npl+1)  = 0.5d0*U*mp
!   hdelm(2:Npl+1)  = 0.5d0*U*mm
!   lambda=0.d0

!   kp = 0.d0!kbz(1,:)
!   call hamiltk(kp,hk)

!   test = transpose(conjg(hk))
!   test2 = hk-test
!   if(sum(abs((test2)))>1.d-10) then
!     write(*,*) 'Hamiltonian is not hermitian!'
!   else
!     write(*,*) 'Hamiltonian is hermitian'
!   end if

!   test = transpose(hk)
!   test2 = hk-test
!   if(sum(abs((test2)))>1.d-10) then
!     write(*,*) 'Hamiltonian is not symmetric!'
!     do i=1,(Npl+2)*18 ; do j=1,(Npl+2)*18
!       if(abs(test2(i,j))>1.d-10) write(*,*) i,j,test2(i,j)
!     end do ; end do
!   else
!     write(*,*) 'Hamiltonian is symmetric'
!   end if


!   write(*,*) '*********************'
!   i = 1
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 7
!   j0 = (j-1)*18+1
!   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,i0:i1))),sum(hk(i0:i1,i0:i1))
!   write(*,*) '11'
!   do i=1,4
!     write(*,'(i4,100f8.5)') i, (hk(i,j),j=1,4)
!   end do

!   write(*,*) '12'
!   do i=1,4
!     write(*,'(i4,100f8.5)') i, (hk(i,j),j=19,22)
!   end do

!   write(*,*) '21'
!   do i=19,22
!     write(*,'(i4,100f8.5)') i, (hk(i,j),j=1,4)
!   end do

!   write(*,*) '22'
!   do i=19,22
!     write(*,'(i4,100f8.5)') i, (hk(i,j),j=19,22)
!   end do

!   write(*,*) '66'
!   do i=91,94 ; do j=91,94
!     write(*,*) i,j, hk(i,j)
!   end do ; end do

!   write(*,*) '67'
!   do i=91,94 ; do j=109,112
!     write(*,*) i,j, hk(i,j)
!   end do ; end do

!   write(*,*) '76'
!   do i=109,112 ; do j=91,94
!     write(*,*) i,j, hk(i,j)
!   end do ; end do

!   write(*,*) '77'
!   do i=109,112 ; do j=109,112
!     write(*,*) i,j, hk(i,j)
!   end do ; end do

! !   stop
!   do i=19,22 ; do j=19,22
! !     i0 = (i-1)*9+1
! !     i1 = i0+8
! !     j0 = (j-1)*9+1
! !     j1 = j0+8
!     write(*,*) i,j, hk(i,j)
!   end do ; end do
! !   stop
!   write(*,*) sum(abs(hk(j0:j1,j0:j1))),sum(hk(j0:j1,j0:j1))
!   write(*,*) '---------------------'
!   i = 2
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 6
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,i0:i1))),sum(hk(i0:i1,i0:i1))
!   write(*,*) sum(abs(hk(j0:j1,j0:j1))),sum(hk(j0:j1,j0:j1))
!   write(*,*) '---------------------'
!   i = 3
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 5
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,i0:i1))),sum(hk(i0:i1,i0:i1))
!   write(*,*) sum(abs(hk(j0:j1,j0:j1))),sum(hk(j0:j1,j0:j1))
!   write(*,*) '*********************'
!   i = 1
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 2
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
!   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
!   i = 6
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 7
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
!   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
!   write(*,*) '---------------------'
!   i = 2
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 3
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
!   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
!   i = 5
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 6
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
!   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
!   write(*,*) '---------------------'
!   i = 3
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 4
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
!   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
!   i = 4
!   i0 = (i-1)*18+1
!   i1 = i0+17
!   j = 5
!   j0 = (j-1)*18+1
!   j1 = j0+17
!   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
!   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
!   write(*,*) '*********************'

!   do i=1,Npl+2 ; do j=1,Npl+2
!     if (abs(i-j)>plnn) then
!       i0 = (i-1)*18+1
!       i1 = i0+17
!       j0 = (j-1)*18+1
!       j1 = j0+17
!       if(sum(abs(hk(i0:i1,j0:j1)))>1.d-10) then
!         write(*,*) 'Non zero values out-of tri-diagonal part!'
!       end if
!     end if
!   end do ; end do

! !   do i=1,(Npl+2)*18 ; do j=1,(Npl+2)*18
! !     if((mod(i,9)==1).and.(mod(j,9)==1)) hk(i,j) = zero !ss

! !     if((mod(i,9)==1).and.(mod(j,9)==2)) hk(i,j) = zero !sp
! !     if((mod(i,9)==1).and.(mod(j,9)==3)) hk(i,j) = zero !sp
! !     if((mod(i,9)==1).and.(mod(j,9)==4)) hk(i,j) = zero !sp

! !     if((mod(i,9)==2).and.(mod(j,9)==1)) hk(i,j) = zero !ps
! !     if((mod(i,9)==3).and.(mod(j,9)==1)) hk(i,j) = zero !ps
! !     if((mod(i,9)==4).and.(mod(j,9)==1)) hk(i,j) = zero !ps

! !     if((mod(i,9)==2).and.(mod(j,9)==2)) hk(i,j) = zero !pp
! !     if((mod(i,9)==2).and.(mod(j,9)==3)) hk(i,j) = zero !pp
! !     if((mod(i,9)==2).and.(mod(j,9)==4)) hk(i,j) = zero !pp

! !     if((mod(i,9)==3).and.(mod(j,9)==2)) hk(i,j) = zero !pp
! !     if((mod(i,9)==3).and.(mod(j,9)==3)) hk(i,j) = zero !pp
! !     if((mod(i,9)==3).and.(mod(j,9)==4)) hk(i,j) = zero !pp

! !     if((mod(i,9)==4).and.(mod(j,9)==2)) hk(i,j) = zero !pp
! !     if((mod(i,9)==4).and.(mod(j,9)==3)) hk(i,j) = zero !pp
! !     if((mod(i,9)==4).and.(mod(j,9)==4)) hk(i,j) = zero !pp

! !     if((mod(i,9)==1).and.(mod(j,9)==5)) hk(i,j) = zero !sd
! !     if((mod(i,9)==1).and.(mod(j,9)==6)) hk(i,j) = zero !sd
! !     if((mod(i,9)==1).and.(mod(j,9)==7)) hk(i,j) = zero !sd
! !     if((mod(i,9)==1).and.(mod(j,9)==8)) hk(i,j) = zero !sd
! !     if((mod(i,9)==1).and.(mod(j,9)==0)) hk(i,j) = zero !sd

! !     if((mod(i,9)==5).and.(mod(j,9)==1)) hk(i,j) = zero !ds
! !     if((mod(i,9)==6).and.(mod(j,9)==1)) hk(i,j) = zero !ds
! !     if((mod(i,9)==7).and.(mod(j,9)==1)) hk(i,j) = zero !ds
! !     if((mod(i,9)==8).and.(mod(j,9)==1)) hk(i,j) = zero !ds
! !     if((mod(i,9)==0).and.(mod(j,9)==1)) hk(i,j) = zero !ds

! !     if((mod(i,9)==2).and.(mod(j,9)==5)) hk(i,j) = zero !pd
! !     if((mod(i,9)==2).and.(mod(j,9)==6)) hk(i,j) = zero !pd
! !     if((mod(i,9)==2).and.(mod(j,9)==7)) hk(i,j) = zero !pd
! !     if((mod(i,9)==2).and.(mod(j,9)==8)) hk(i,j) = zero !pd
! !     if((mod(i,9)==2).and.(mod(j,9)==0)) hk(i,j) = zero !pd

! !     if((mod(i,9)==3).and.(mod(j,9)==5)) hk(i,j) = zero !pd
! !     if((mod(i,9)==3).and.(mod(j,9)==6)) hk(i,j) = zero !pd
! !     if((mod(i,9)==3).and.(mod(j,9)==7)) hk(i,j) = zero !pd
! !     if((mod(i,9)==3).and.(mod(j,9)==8)) hk(i,j) = zero !pd
! !     if((mod(i,9)==3).and.(mod(j,9)==0)) hk(i,j) = zero !pd

! !     if((mod(i,9)==4).and.(mod(j,9)==5)) hk(i,j) = zero !pd
! !     if((mod(i,9)==4).and.(mod(j,9)==6)) hk(i,j) = zero !pd
! !     if((mod(i,9)==4).and.(mod(j,9)==7)) hk(i,j) = zero !pd
! !     if((mod(i,9)==4).and.(mod(j,9)==8)) hk(i,j) = zero !pd
! !     if((mod(i,9)==4).and.(mod(j,9)==0)) hk(i,j) = zero !pd

! !     if((mod(i,9)==5).and.(mod(j,9)==2)) hk(i,j) = zero !dp
! !     if((mod(i,9)==6).and.(mod(j,9)==2)) hk(i,j) = zero !dp
! !     if((mod(i,9)==7).and.(mod(j,9)==2)) hk(i,j) = zero !dp
! !     if((mod(i,9)==8).and.(mod(j,9)==2)) hk(i,j) = zero !dp
! !     if((mod(i,9)==0).and.(mod(j,9)==2)) hk(i,j) = zero !dp

! !     if((mod(i,9)==5).and.(mod(j,9)==3)) hk(i,j) = zero !dp
! !     if((mod(i,9)==6).and.(mod(j,9)==3)) hk(i,j) = zero !dp
! !     if((mod(i,9)==7).and.(mod(j,9)==3)) hk(i,j) = zero !dp
! !     if((mod(i,9)==8).and.(mod(j,9)==3)) hk(i,j) = zero !dp
! !     if((mod(i,9)==0).and.(mod(j,9)==3)) hk(i,j) = zero !dp

! !     if((mod(i,9)==5).and.(mod(j,9)==4)) hk(i,j) = zero !dp
! !     if((mod(i,9)==6).and.(mod(j,9)==4)) hk(i,j) = zero !dp
! !     if((mod(i,9)==7).and.(mod(j,9)==4)) hk(i,j) = zero !dp
! !     if((mod(i,9)==8).and.(mod(j,9)==4)) hk(i,j) = zero !dp
! !     if((mod(i,9)==0).and.(mod(j,9)==4)) hk(i,j) = zero !dp

! !     if((mod(i,9)==5).and.(mod(j,9)==5)) hk(i,j) = zero !dd
! !     if((mod(i,9)==5).and.(mod(j,9)==6)) hk(i,j) = zero !dd
! !     if((mod(i,9)==5).and.(mod(j,9)==7)) hk(i,j) = zero !dd
! !     if((mod(i,9)==5).and.(mod(j,9)==8)) hk(i,j) = zero !dd
! !     if((mod(i,9)==5).and.(mod(j,9)==0)) hk(i,j) = zero !dd

! !     if((mod(i,9)==6).and.(mod(j,9)==5)) hk(i,j) = zero !dd
! !     if((mod(i,9)==6).and.(mod(j,9)==6)) hk(i,j) = zero !dd
! !     if((mod(i,9)==6).and.(mod(j,9)==7)) hk(i,j) = zero !dd
! !     if((mod(i,9)==6).and.(mod(j,9)==8)) hk(i,j) = zero !dd
! !     if((mod(i,9)==6).and.(mod(j,9)==0)) hk(i,j) = zero !dd

! !     if((mod(i,9)==7).and.(mod(j,9)==5)) hk(i,j) = zero !dd
! !     if((mod(i,9)==7).and.(mod(j,9)==6)) hk(i,j) = zero !dd
! !     if((mod(i,9)==7).and.(mod(j,9)==7)) hk(i,j) = zero !dd
! !     if((mod(i,9)==7).and.(mod(j,9)==8)) hk(i,j) = zero !dd
! !     if((mod(i,9)==7).and.(mod(j,9)==0)) hk(i,j) = zero !dd

! !     if((mod(i,9)==8).and.(mod(j,9)==5)) hk(i,j) = zero !dd
! !     if((mod(i,9)==8).and.(mod(j,9)==6)) hk(i,j) = zero !dd
! !     if((mod(i,9)==8).and.(mod(j,9)==7)) hk(i,j) = zero !dd
! !     if((mod(i,9)==8).and.(mod(j,9)==8)) hk(i,j) = zero !dd
! !     if((mod(i,9)==8).and.(mod(j,9)==0)) hk(i,j) = zero !dd

! !     if((mod(i,9)==0).and.(mod(j,9)==5)) hk(i,j) = zero !dd
! !     if((mod(i,9)==0).and.(mod(j,9)==6)) hk(i,j) = zero !dd
! !     if((mod(i,9)==0).and.(mod(j,9)==7)) hk(i,j) = zero !dd
! !     if((mod(i,9)==0).and.(mod(j,9)==8)) hk(i,j) = zero !dd
! !     if((mod(i,9)==0).and.(mod(j,9)==0)) hk(i,j) = zero !dd
! !   end do ; end do



!   write(*,*) 'Green function:'
!   lambda=0.d0
!   kp = 0.d0 !kbz(1,:)
!   call green_es(-1.d0,0.d0,kp,gf)

!   write(*,*) '*********************'
!   write(*,*) sum(abs(gf(1,1,:,:))),sum(gf(1,1,:,:))
!   write(*,*) sum(abs(gf(7,7,:,:))),sum(gf(7,7,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(gf(2,2,:,:))),sum(gf(2,2,:,:))
!   write(*,*) sum(abs(gf(6,6,:,:))),sum(gf(6,6,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(gf(3,3,:,:))),sum(gf(3,3,:,:))
!   write(*,*) sum(abs(gf(5,5,:,:))),sum(gf(5,5,:,:))
!   write(*,*) '*********************'
!   write(*,*) sum(abs(gf(1,2,:,:))),sum(gf(1,2,:,:))
!   write(*,*) sum(abs(gf(6,7,:,:))),sum(gf(6,7,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(gf(2,3,:,:))),sum(gf(2,3,:,:))
!   write(*,*) sum(abs(gf(5,6,:,:))),sum(gf(5,6,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(gf(3,4,:,:))),sum(gf(3,4,:,:))
!   write(*,*) sum(abs(gf(4,5,:,:))),sum(gf(4,5,:,:))
!   write(*,*) '*********************'
!   write(*,*) sum(abs(gf(1,3,:,:))),sum(gf(1,3,:,:))
!   write(*,*) sum(abs(gf(5,7,:,:))),sum(gf(5,7,:,:))
!   write(*,*) '---------------------'
!   write(*,*) sum(abs(gf(2,4,:,:))),sum(gf(2,4,:,:))
!   write(*,*) sum(abs(gf(4,6,:,:))),sum(gf(4,6,:,:))
!   write(*,*) '*********************'

! !   test = transpose(conjg(hk))
! !   test2 = hk-test
! !   if(sum(abs((test2)))>1.d-10) then
! !     write(*,*) 'Green function is not hermitian!'
! !   else
! !     write(*,*) 'Green function is hermitian'
! !   end if

! !   test = transpose(hk)
! !   test2 = hk-test
! !   if(sum(abs((test2)))>1.d-10) then
! !     write(*,*) 'Green function is not symmetric!'
! !     do i=1,(Npl+2)*18 ; do j=1,(Npl+2)*18
! !       if(abs(test2(i,j))>1.d-10) write(*,*) i,j,test2(i,j)
! !     end do ; end do
! !   else
! !     write(*,*) 'Green function is symmetric'
! !   end if


! !   write(*,*) '*********************'
! !   i = 1
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 7
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) 'Diagonais'
! !   do i=1,18
! !     write(*,*) '11, i = ',i,hk(i,i)
! !     write(*,*) '77, i = ',i,hk((j-1)*18+i,(j-1)*18+i)
! !   end do
! !   write(*,*) sum(abs(hk(i0:i1,i0:i1))),sum(hk(i0:i1,i0:i1))
! !   write(*,*) sum(abs(hk(j0:j1,j0:j1))),sum(hk(j0:j1,j0:j1))
! !   test(1:18,1:18) = hk(i0:i1,i0:i1) - hk(j0:j1,j0:j1)
! !   if(sum(abs(test(1:18,1:18)))>1.d-10) then
! !     write(*,*) '11 e 77 diferentes'
! !     do i=1,18 ; do j=1,18
! !       if(abs(test(i,j))>1.d-10) write(*,*) i,j,test(i,j)
! !     end do ; end do
! !   end if

! !   write(*,*) '---------------------'
! !   i = 2
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 6
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,i0:i1))),sum(hk(i0:i1,i0:i1))
! !   write(*,*) sum(abs(hk(j0:j1,j0:j1))),sum(hk(j0:j1,j0:j1))
! !   write(*,*) '---------------------'
! !   i = 3
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 5
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,i0:i1))),sum(hk(i0:i1,i0:i1))
! !   write(*,*) sum(abs(hk(j0:j1,j0:j1))),sum(hk(j0:j1,j0:j1))
! !   write(*,*) '*********************'
! !   i = 1
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 2
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
! !   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
! !   i = 6
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 7
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
! !   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
! !   write(*,*) '---------------------'
! !   i = 2
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 3
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
! !   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
! !   i = 5
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 6
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
! !   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
! !   write(*,*) '---------------------'
! !   i = 3
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 4
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
! !   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
! !   i = 4
! !   i0 = (i-1)*18+1
! !   i1 = i0+17
! !   j = 5
! !   j0 = (j-1)*18+1
! !   j1 = j0+17
! !   write(*,*) sum(abs(hk(i0:i1,j0:j1))),sum(hk(i0:i1,j0:j1))
! !   write(*,*) sum(abs(hk(j0:j1,i0:i1))),sum(hk(j0:j1,i0:i1))
! !   write(*,*) '*********************'

!   Finalizing program
  if(myrank==0) call write_time(outputunit,'[main] Finished on: ')
  call MPI_Finalize(ierr)
  if ((ierr/=0).and.(myrank==0)) write(outputunit,"('[main] Something went wrong in the parallelization! ierr = ',i0)") ierr
  stop

  return
end subroutine debugging