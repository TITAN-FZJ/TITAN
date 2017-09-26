!-------------------------------------------------------------------------------
! TITAN - Time-dependent description of Itinerant electrons: Transport and Angular momentum properties of Nanostructures
!-------------------------------------------------------------------------------
!
! MODULE: mod_system
!
!> @author
!> Jens Renè Suckert, PGI-1/IAS-1, Peter-Grünberg-Institut,FZ Jülich
!
! DESCRIPTION:
!> Describes the system in real and reciprocal space.
!> Contains function to generate all neighbors and the Brillouin Zone for a layered system or bulk.
!
! REVISION HISTORY:
! 28 April 2017 - Initial Version
!-------------------------------------------------------------------------------
module mod_system
  use mod_f90_kind, only: double
  use AtomTypes, only: BasisAtom, NeighborAtom, AtomType
  implicit none

  type :: System
    character(len=200) :: Name = ""
    real(double), dimension(3) :: a1, a2, a3
    real(double), dimension(3) :: b1, b2, b3
    real(double) :: a0
    logical :: lbulk = .false.

    integer :: nAtoms = 0
    type(BasisAtom), dimension(:), allocatable :: Basis
    integer :: nNeighbors
    type(NeighborAtom), dimension(:), allocatable :: Neighbors
    integer :: nStages = 0
    real(double), dimension(:,:), allocatable :: Distances ! size (nStages, nAtoms)
    integer :: nTypes = 0
    type(AtomType), dimension(:), allocatable :: Types
    integer :: nkpt = 0
    real(double), dimension(:,:), allocatable :: kbz ! size (3,nkpt)
    real(double), dimension(:), allocatable :: wkbz
  end type System

  type(System) :: sys

  integer :: n0sc1 !< first neighbor to calculate the in-plane spin and charge current
  integer :: n0sc2 !< last neighbor to calculate the in-plane spin and charge current
  integer :: n0sc  !< Number of neighbors to calculate currents
  real(double), dimension(3) :: pln_normal

  integer, dimension(:,:), allocatable :: ia



 contains

   subroutine initHamiltkStride(nAtoms, nOrb)
     implicit none
     integer, intent(in) :: nAtoms, nOrb
     integer :: i
     allocate(ia(4,nAtoms))
     do i = 1, nAtoms
        ia(1,i) = (i-1) * 2 * nOrb + 1
        ia(2,i) = ia(1,i) + nOrb - 1
        ia(3,i) = ia(1,i) + nOrb
        ia(4,i) = ia(3,i) + nOrb - 1
     end do
     return
   end subroutine initHamiltkStride
!
!   subroutine generate_neighbors()
!     use mod_f90_kind, only: double
!     use mod_tools, only: is_parallel
!     implicit none
!
!     real(double), dimension(:,:), allocatable :: cube_pos
!     real(double), dimension(:,:), allocatable :: cube_dist
!     real(double), dimension(3) :: tmp_pos
!     real(double), dimension(2) :: tmp_dist
!     integer :: i,j,k,l,m
!     logical :: flag
!     allocate( cube_pos(3,(4*nstages+1)**3-1))
!     allocate(cube_dist(2,(4*nstages+1)**3-1))
!     allocate(stages_list(4*nstages))
!     allocate(planes_list(4*nstages))
!
!     stages_list = 0.d0
!     planes_list = 0.d0
!
!     ! Generate all neighbors in cuboid
!     l = 0
!     do i = -2*nstages, 2*nstages
!        do j = -2*nstages, 2*nstages
!           do k = -2*nstages, 2*nstages
!              if(0 == i .and. 0 == j .and. 0 == k) cycle
!
!              tmp_pos  = i*a1 + j*a2 + k*a3
!              tmp_dist(1) = sqrt(dot_product(tmp_pos, tmp_pos))
!              tmp_dist(2) = dot_product(tmp_pos, pln_normal)
!
!              if(tmp_dist(2) > 1d-9) cycle
!
!              m = l
!              l = l + 1
!
!              do while(1 <= m)
!                 if(abs(tmp_dist(2)) - abs(cube_dist(2,m)) > 1.0d-9 ) exit
!                 if(abs(abs(tmp_dist(2)) - abs(cube_dist(2,m))) < 1.0d-9 .and. cube_dist(1,m) <= tmp_dist(1)) exit
!                 cube_pos(:, m+1) = cube_pos(:, m)
!                 cube_dist(:,m+1) = cube_dist(:,m)
!                 m = m - 1
!              end do
!              cube_pos(1:3, m+1) = tmp_pos(1:3)
!              cube_dist(1:2,m+1) = tmp_dist(1:2)
!           end do
!        end do
!     end do
!     ! Find first nstages next neighbor distances
!     j = 0
!     do i = 1, l
!        flag = .false.
!        do k = 1, j
!           if(abs(cube_dist(1,i) - stages_list(k)) <= 1.0d-9) flag = .true.
!        end do
!        if(.not. flag) then
!           if(j > nstages .and. cube_dist(1,i) > stages_list(nstages)) cycle
!           m = j
!           j = j + 1
!           do while(1 <= m)
!              if(cube_dist(1,i) - stages_list(m) > 1d-9) exit
!              stages_list(m+1) = stages_list(m)
!              m = m - 1
!           end do
!           stages_list(m+1) = cube_dist(1,i)
!        end if
!     end do
!     ! Find amount of available planes
!     npln = 0
!     j = 0
!     do i = 1, l
!        if(cube_dist(1,i) - stages_list(nstages) > 1d-9) cycle
!        flag = .false.
!        do k = 1, npln
!           if(abs(cube_dist(2,i) - planes_list(npln)) <= 1d-9) flag = .true.
!        end do
!        if(.not. flag) then
!           npln = npln + 1
!           planes_list(npln) = cube_dist(2,i)
!        end if
!     end do
!     ! Count final number of elements
!     ncount = 0
!     do i = 1, l
!        if(planes_list(npln) - cube_dist(2,i) > 1d-9) cycle
!        if(cube_dist(1,i) - stages_list(nstages) > 1d-9) cycle
!        ncount = ncount + 1
!     end do
!
!     allocate(r_nn(3,ncount))
!     allocate(c_nn(3,ncount))
!     allocate(l_nn(nstages+1, npln+1))
!
!     l_nn = 1
!     m = 0
!     do i = 1, l
!        if(planes_list(npln) - cube_dist(2,i) > 1d-9) cycle
!        if(cube_dist(1,i) - stages_list(nstages) > 1d-9) cycle
!
!        do j = 1, npln
!           if(abs(cube_dist(2,i) - planes_list(j)) <= 1d-9) exit
!        end do
!        do k = 1, nstages
!           if(abs(cube_dist(1,i) - stages_list(k)) <= 1d-9) exit
!        end do
!
!        l_nn(k+1:, j) = l_nn(k+1:, j) + 1
!        if(j < npln+1) l_nn(:, j+1:) = l_nn(:, j+1:) + 1
!
!        m = m + 1
!        r_nn(1:3, m) = cube_pos(1:3, i)
!        c_nn(1:3, m) = cube_pos(1:3, i) / cube_dist(1,i)
!     end do
!
!     allocate(pln_cnt(npln))
!     do i = 1, npln
!       pln_cnt(i) = l_nn(1,i+1)-l_nn(1,i)
!     end do
!
!     if(dot_product(pln_normal, pln_normal) > 1d-9) then
!        if(l_nn(1,1) - l_nn(1,2) == 0) stop !TODO: Error message
!        pln_a1 = r_nn(:,l_nn(1,1))
!        do i = l_nn(1,1)+1, l_nn(1,2) - 1
!           if(.not. is_parallel(pln_a1, r_nn(1:3,i))) then
!              pln_a2 = r_nn(1:3,i)
!              exit
!           end if
!        end do
!        if(dot_product(pln_a2, pln_a2) < 1d-9) stop !TODO: Error message
!     end if
!
!     if(n0sc1 <= 0 .or. n0sc1 > pln_cnt(1)) n0sc1 = 1
!     if(n0sc2 <= 0 .or. n0sc1 < n0sc1) n0sc2 = n0sc1
!     if(n0sc2 > pln_cnt(1)) n0sc2 = pln_cnt(1)
!     n0sc = n0sc2 - n0sc2 + 1
!
!
!     deallocate(cube_pos)
!     deallocate(cube_dist)
!   end subroutine generate_neighbors
!
!   subroutine generate_kpoints()
!     implicit none
!
!     if(1d-9 > dot_product(pln_normal, pln_normal)) then
!        call generate_kpoints_bulk()
!     else
!        call generate_kpoints_plane()
!     end if
!     return
!   end subroutine generate_kpoints
!
!   subroutine generate_kpoints_plane()
!     use mod_f90_kind,  only: double
!     use mod_tools,     only: cross
!     use mod_constants, only: tpi
!     implicit none
!     real(double), dimension(3,4) :: BZ
!     real(double), dimension(3)   :: diff
!     real(double), dimension(:,:), allocatable :: extrakbz
!     real(double), dimension(:), allocatable :: extrawkbz
!     real(double), dimension(:,:), allocatable :: inikbz
!     real(double), dimension(:), allocatable :: iniwkbz
!     real(double) :: smallest_dist, distance, ini_smallest_dist, vol
!     integer :: j,l,m, smallest_index, numextrakbz, nkpt_perdim
!
!     allocate(extrakbz(3,nkpt*10))
!     allocate(extrawkbz(nkpt*10))
!     nkpt_perdim = ceiling(sqrt(dble(nkpt)))
!     nkpt_perdim = nkpt_perdim
!     nkpt = nkpt_perdim**2
!
!     allocate(iniwkbz(nkpt), inikbz(3,nkpt))
!
!     vol = tpi / dot_product(pln_normal, cross(pln_a1, pln_a2))
!     b1 = vol * cross(pln_a1, pln_normal)
!     b2 = vol * cross(pln_normal, pln_a2)
!
!     BZ(:,1) = 0.d0
!     BZ(:,2) = b1
!     BZ(:,3) = b2
!     BZ(:,4) = b1 + b2
!
!     !Generate k-points in the paralelogram determined by b1 and b2
!     iniwkbz = 1.d0
!     m = 0
!     do l = 1, nkpt_perdim;
!        do j = 1, nkpt_perdim !Run over the entire BZ
!           m = m + 1
!           inikbz(:, m) = ( dble(l-1)*b1 + dble(j-1)*b2 ) / dble(nkpt_perdim)
!        end do
!     end do
!
!     !Translate the k-points to the 1st BZ.
!     !10*|b1+b2|, bigger than the distance of any genarated kpoint
!     ini_smallest_dist = 10.d0 * sqrt(dot_product(b1+b2, b1+b2))
!     numextrakbz = 0
!     !Run over all the kpoints generated initially.
!     do l = 1, nkpt
!        smallest_dist = ini_smallest_dist
!        m = 0
!        !Checks to each of the 4 BZ's the kpoint belongs by checking
!        ! to each BZ it's closer.
!        do j = 1, 4
!           diff = inikbz(:,l) - BZ(:,j)
!           distance = sqrt(dot_product(diff, diff))
!           if(distance < smallest_dist) then
!              smallest_dist = distance
!              smallest_index = j
!           end if
!        end do
!        !Checks if the kpoint is in the border between two or more
!        ! BZ's. If yes, create a clone of it to translate later into
!        ! the 1st BZ.
!        do j = 1, 4
!           diff = inikbz(:, l) - BZ(:, j)
!           distance = sqrt(dot_product(diff, diff))
!           if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. j /= smallest_index ) then
!              m = m + 1
!              numextrakbz = numextrakbz + 1
!              extrakbz(:, numextrakbz) = inikbz(:, l) - BZ(:,j)
!           end if
!        end do
!        if(m /= 0) then
!           !The weight of the kpoint in the border is shared with
!           ! its clones.
!           iniwkbz(l) = 1.d0 / dble(m+1)
!           extrawkbz(numextrakbz - m + 1 : numextrakbz) = 1.d0 / dble(m+1)
!        end if
!        !Translate the kpoint to the 1st BZ
!        inikbz(:, l) = inikbz(:, l) - BZ(:, smallest_index)
!     end do
!
!     allocate( wkbz(nkpt + numextrakbz), kbz(3, nkpt + numextrakbz) )
!     !The final array of kpoints will be the initial one by the clones
!     ! of the ones in the border between BZ's.
!     kbz (:,1:nkpt) = inikbz (:,:)
!     wkbz(  1:nkpt) = iniwkbz(:)
!     if(numextrakbz/=0) then
!        kbz(:, nkpt + 1 : nkpt + numextrakbz) = extrakbz(:, 1 : numextrakbz)
!        wkbz(nkpt + 1 : nkpt + numextrakbz) = extrawkbz(1 : numextrakbz)
!     end if
!     wkbz=wkbz/dble(nkpt)
!     nkpt=nkpt + numextrakbz
!     return
!   end subroutine generate_kpoints_plane
!
!
!   subroutine generate_kpoints_bulk()
!     use mod_f90_kind,   only: double
!     use mod_tools,      only: cross
!     use mod_constants,  only: tpi
!
!     implicit none
!     real(double) :: vol
!     real(double), dimension(3,8) :: BZ
!     real(double) :: smallest_dist, distance, ini_smallest_dist
!     real(double), dimension(3) ::  diff
!     real(double), allocatable :: inikbz(:,:), iniwkbz(:)
!     real(double), allocatable, dimension(:,:) :: extrakbz
!     real(double), allocatable, dimension(:) :: extrawkbz
!     integer :: l,j,m, k, smallest_index, numextrakbz
!     integer :: nkpt_perdim !n. of k point per dimension
!
!     allocate(extrakbz(3,nkpt*10), extrawkbz(nkpt*10))
!
!     nkpt_perdim=ceiling((dble(nkpt))**(1.d0/3.d0))
!     nkpt_perdim=nkpt_perdim
!     nkpt = nkpt_perdim**3
!     allocate( iniwkbz(nkpt), inikbz(3, nkpt) )
!
!     vol = tpi / dot_product(a1, cross(a2,a3))
!     b1 = vol * cross(a2,a3)
!     b2 = vol * cross(a3,a1)
!     b3 = vol * cross(a1,a2)
!
!     BZ(:,1)=0.d0
!     BZ(:,2)=b1
!     BZ(:,3)=b2
!     BZ(:,4)=b1+b2
!     BZ(:,5)=b3
!     BZ(:,6)=b1+b3
!     BZ(:,7)=b2+b3
!     BZ(:,8)=b1+b2+b3
!
!     !Generate k-points in the paralelogram determined by b1 and b2
!     iniwkbz = 1.d0
!     m = 0
!     do l = 1, nkpt_perdim
!        do j = 1, nkpt_perdim
!           do k = 1, nkpt_perdim !Run over the entire BZ
!              m = m + 1
!              inikbz(:, m)= ( dble(l-1)*b1 + dble(j-1)*b2 + dble(k-1)*b3)/ dble(nkpt_perdim)
!           end do
!        end do
!     end do
!
!     !Translate the k-points to the 1st BZ.
!     !10*|b1+b2|, bigger than the distance of any genarated kpoint
!     ini_smallest_dist=10.d0*sqrt(dot_product(b1+b2,b1+b2))
!     numextrakbz=0
!     !Run over all the kpoints generated initially.
!     do l=1, nkpt
!        smallest_dist=ini_smallest_dist
!        m=0
!        !Checks to each of the 4 BZ's the kpoint belongs by checking
!        ! to each BZ it's closer.
!        do j=1, 8
!           diff=inikbz(:,l)-BZ(:,j)
!           distance=sqrt(dot_product(diff,diff))
!           if(distance<smallest_dist) then
!              smallest_dist=distance
!              smallest_index=j
!           end if
!        end do
!        !Checks if the kpoint is in the border between two or more
!        ! BZ's. If yes, create a clone of it to translate later into
!        ! the 1st BZ.
!        do j=1, 8
!           diff=inikbz(:,l)-BZ(:,j)
!           distance=sqrt(dot_product(diff,diff))
!           if( ( abs(distance-smallest_dist) < 1.d-12 ) .and. &
!                j/=smallest_index ) then
!              m=m+1
!              numextrakbz=numextrakbz+1
!              extrakbz(:,numextrakbz)=inikbz(:,l)-BZ(:,j)
!           end if
!        end do
!        if(m/=0) then
!           !The weight of the kpoint in the border is shared with
!           ! its clones.
!           iniwkbz(l)=1.d0/dble(m+1)
!           extrawkbz(numextrakbz-m+1:numextrakbz)=1.d0/dble(m+1)
!        end if
!        !Translate the kpoint to the 1st BZ
!        inikbz(:,l)=inikbz(:,l)-BZ(:,smallest_index)
!     end do
!
!     allocate( wkbz(nkpt+numextrakbz),kbz(3,nkpt+numextrakbz) )
!     !The final array of kpoints will be the initial one by the clones
!     ! of the ones in the border between BZ's.
!     kbz (:,1:nkpt)=inikbz (:,:)
!     wkbz(1:nkpt)  =iniwkbz(:)
!     if(numextrakbz/=0) then
!        kbz(:,nkpt+1:nkpt+numextrakbz)=extrakbz(:,1:numextrakbz)
!        wkbz(nkpt+1:nkpt+numextrakbz) =extrawkbz(1:numextrakbz)
!     end if
!     wkbz=wkbz/dble(nkpt)
!     nkpt=nkpt + numextrakbz
!
!     deallocate(extrakbz, extrawkbz)
!
!     return
!   end subroutine generate_kpoints_bulk
!
!   subroutine write_neighbors_to_file()
!     implicit none
!     integer :: i, j, k
!
!     open(unit = 5675, file='positions', status='replace')
!     write(5675, "('#In plane basis: ',3(f12.9,2x),' ',3(f12.9,2x))") pln_a1, pln_a2
!
!     write(5675, *) "# Number   NN-Stage  NN-Distance  Position   Dir. Cosine"
!     do i=1,npln
!        write(5675, "('# Plane ',i0,': ',f12.9)") i, planes_list(i)
!        do j = 1, nstages
!           do k = l_nn(j,i), l_nn(j+1,i)-1
!              write(5675, "('',I0,' ',I0,' ',f12.9,' ',3(f12.9,2x),' ',3(f12.9,2x),'')") k, j, stages_list(j), r_nn(:,k), c_nn(:,k)
!           end do
!        end do
!     end do
!
!     close(unit = 5675)
!
!     return
!   end subroutine write_neighbors_to_file
!
!   subroutine write_kpoints_to_file()
!     implicit none
!     integer :: i
!
!     !open (unit=2222, file='kpoints2d',status='replace')
!     !write(unit=2222,fmt="(a)") ' #      kx            ky            wk'
!     open (unit=3333, file='kpoints',status='replace')
!     write(unit=3333,fmt="(a)") ' #      kx            ky            kz            wk'
!     do i=1,nkpt
!        !write(unit=2222,fmt="(3(f12.9,2x))") kbz2d(i,1),kbz2d(i,2),wkbz(i)
!        write(unit=3333,fmt="(4(f12.9,2x))") kbz(:,i),wkbz(i)
!     end do
!     !close(2222)
!     close(3333)
!
!     return
!   end subroutine write_kpoints_to_file


end module mod_system
