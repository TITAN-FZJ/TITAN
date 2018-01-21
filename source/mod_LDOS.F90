module mod_LDOS
   use mod_f90_kind, only: double
   implicit none
   character(len=4), private :: folder = "LDOS"
   character(len=6), dimension(4), private :: filename = ["ldosu","ldosd","lmdosu", "lmdosd"]

   real(double),dimension(:,:),allocatable :: ldosu,ldosd

contains
   subroutine allocateLDOS()
      use mod_System, only: s => sys
      use TightBinding, only: nOrb
      implicit none

      if(allocated(ldosu)) deallocate(ldosu)
      if(allocated(ldosd)) deallocate(ldosd)
      allocate(ldosu(s%nAtoms,nOrb))
      allocate(ldosd(s%nAtoms,nOrb))

      return
   end subroutine allocateLDOS

   subroutine deallocateLDOS()
      implicit none
      if(allocated(ldosu)) deallocate(ldosu)
      if(allocated(ldosd)) deallocate(ldosd)
      return
   end subroutine deallocateLDOS

   subroutine openLDOSFiles()
      use mod_parameters, only: output
      use mod_System, only: s => sys
      implicit none
      character(len=400) :: varm
      integer            :: i, iw, j

      do i = 1, s%nAtoms
         do j = 1, size(filename)
            iw = 1000 + (i-1) * size(filename) + j
            write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_site=',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(j)),i,trim(output%info),trim(output%BField),trim(output%SOC)
            open (unit=iw, file=varm,status='replace')
            write(unit=iw, fmt="('#   energy      ,  LDOS SUM        ,  LDOS S          ,  LDOS P          ,  LDOS EG         ,  LDOS T2G        ,')")
         end do
      end do
      return
   end subroutine openLDOSFiles

   subroutine closeLDOSFiles()
      use mod_System, only: s => sys
      implicit none
      integer i,j, iw
      do i = 1, s%nAtoms
         do j = 1, size(filename)
            iw = 1000 + (i-1) * size(filename) + j
            close(iw)
         end do
      end do
      return
   end subroutine closeLDOSFiles

   subroutine writeLDOS(e)
      use mod_f90_kind, only: double
      use mod_System, only: s => sys
      implicit none
      real(double), intent(in) :: e
      integer :: i, iw, j

      do i = 1, s%nAtoms
          iw = 1000 + (i-1) * 2 + 1
          write(unit=iw,fmt="( 6(es16.9,2x))") e, sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:7)),sum(ldosu(i,8:9))
          iw = iw + 1
          write(unit=iw,fmt="( 6(es16.9,2x))") e, sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:7)),sum(ldosd(i,8:9))
          iw = iw + 1
          write(unit=iw,fmt="(10(es16.9,2x))") e, (ldosu(i,j), j = 1,9) !sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
          iw = iw + 1
          write(unit=iw,fmt="(10(es16.9,2x))") e, (ldosd(i,j), j = 1,9) !sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
      end do

   end subroutine writeLDOS

end module mod_LDOS
