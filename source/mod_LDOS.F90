module mod_LDOS
  use mod_f90_kind, only: double
  implicit none
  character(len=4), private :: folder = "LDOS"
  character(len=6), dimension(4), private :: filename = ["ldosu ","ldosd ","lmdosu", "lmdosd"]

  real(double),dimension(:,:),allocatable :: ldosu,ldosd

contains
  subroutine allocateLDOS()
    use mod_System,     only: s => sys
    use mod_parameters, only: nOrb
    use mod_superconductivity, only: superCond
    implicit none

    if(allocated(ldosu)) deallocate(ldosu)
    if(allocated(ldosd)) deallocate(ldosd)
    allocate(ldosu(s%nAtoms,nOrb*superCond))
    allocate(ldosd(s%nAtoms,nOrb*superCond))

  end subroutine allocateLDOS

  subroutine deallocateLDOS()
    implicit none
    if(allocated(ldosu)) deallocate(ldosu)
    if(allocated(ldosd)) deallocate(ldosd)
  end subroutine deallocateLDOS

  subroutine createLDOSFiles()
    use mod_parameters,        only: output
    use mod_System,            only: s => sys
    use mod_io,                only: write_header
    use mod_superconductivity, only: lsuperCond
    implicit none
    character(len=400) :: varm,title(size(filename))
    integer            :: i, iw, j

    title(1:2) = "#    energy      ,  LDOS SUM       ,  LDOS S         ,  LDOS P         ,  LDOS D         "
    title(3:4) = "#    energy      ,  LDOS S         ,  LDOS PX        ,  LDOS PY        ,  LDOS PZ        ,  LDOS DXY       ,  LDOS DYZ       ,  LDOS DZX       ,  LDOS DX2       ,  LDOS DZ2       "

    do i = 1, s%nAtoms
      do j = 1, size(filename)
        iw = 1000 + (i-1) * size(filename) + j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_site=',i0,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(j)),i,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open (unit=iw, file=varm,status='replace', form='formatted')
        !                             if lsuperCond is true then it writes 0.0, otherwise it writes s%Ef
        call write_header(iw,title(j),merge(0.0,s%Ef,lsuperCond))
        ! write(unit=iw, fmt="('#   energy      ,  LDOS SUM        ,  LDOS S          ,  LDOS P          ,  LDOS T2G        ,  LDOS EG         ')")
        close(unit=iw)
      end do
    end do
  end subroutine createLDOSFiles

  subroutine openLDOSFiles()
    use mod_parameters, only: output,missing_files
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=400) :: varm
    integer            :: i, iw, j, err, errt=0

    do i = 1, s%nAtoms
      do j = 1, size(filename)
        iw = 1000 + (i-1) * size(filename) + j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_site=',i0,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(j)),i,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
      end do
    end do
    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openLDOSFiles] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

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
  end subroutine closeLDOSFiles

  subroutine writeLDOS(e)
    use mod_f90_kind, only: double
    use mod_System,   only: s => sys
    use mod_superconductivity, only: lsuperCond
    implicit none
    real(double), intent(in) :: e
    integer :: i, iw, j

    if(lsuperCond) then
        do i = 1, s%nAtoms
           iw = 1000 + (i-1) * size(filename) + 1
           write(unit=iw,fmt="( 8(es16.9,2x))") e, sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9)), ldosu(i,10),sum(ldosu(i,11:13)),sum(ldosu(i,14:18))
           iw = iw + 1
           write(unit=iw,fmt="( 8(es16.9,2x))") e, sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9)), ldosd(i,10),sum(ldosd(i,11:13)),sum(ldosd(i,14:18))
           iw = iw + 1
           write(unit=iw,fmt="(19(es16.9,2x))") e, (ldosu(i,j), j = 1,18) !sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="(19(es16.9,2x))") e, (ldosd(i,j), j = 1,18) !sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
        end do
    else
        do i = 1, s%nAtoms
           iw = 1000 + (i-1) * size(filename) + 1
           write(unit=iw,fmt="( 5(es16.9,2x))") e, sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="( 5(es16.9,2x))") e, sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="(10(es16.9,2x))") e, (ldosu(i,j), j = 1,9) !sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="(10(es16.9,2x))") e, (ldosd(i,j), j = 1,9) !sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
        end do
    end if

  end subroutine writeLDOS

  ! This subroutine sorts LDOS files
  subroutine sortLDOS()
    use mod_f90_kind, only: double
    use mod_tools,    only: sort_file
    use mod_system,   only: s => sys
    implicit none
    integer :: i,j,iw

    ! Opening LDOS files
    call openLDOSFiles()

    do i = 1, s%nAtoms
      do j = 1, size(filename)
        iw = 1000 + (i-1) * size(filename) + j
        call sort_file(iw)
      end do
    end do

    ! Closing LDOS files
    call closeLDOSFiles()

  end subroutine sortLDOS

end module mod_LDOS
