module mod_currents
  use, intrinsic :: iso_fortran_env
  implicit none
  ! Currents, total currents and renormalized currents: Ich, Isx, Isy, Isz, Ilx, Ily,Ilz
  complex(dp),allocatable :: currents(:,:,:)
  !! Currents
  complex(dp),allocatable :: total_currents(:,:)
  !! Total Currents
  complex(dp),allocatable :: rcurrents(:,:,:)
  !! Renormalized Currents
  complex(dp),allocatable :: rtotal_currents(:,:)
  !! Total renormalized currents

  real(dp),allocatable      :: dc_currents(:,:)
  ! Full response functions
  complex(dp), dimension(:,:,:), allocatable :: ttchiorbiikl
  !! Currents response function
  complex(dp), dimension(:,:,:), allocatable :: Lxttchiorbiikl
  !! Orbital Current Lx response function
  complex(dp), dimension(:,:,:), allocatable :: Lyttchiorbiikl
  !! Orbital Current Ly response function
  complex(dp), dimension(:,:,:), allocatable :: Lzttchiorbiikl
  !! Orbital Current Lz response function
contains

  ! This subroutine allocates variables related to the current calculation
  subroutine allocate_currents()
    use, intrinsic :: iso_fortran_env
    use mod_mpi_pars
    use mod_prefactors,only: prefactor,prefactorlsoc
    use mod_parameters, only: Npl,dimspinAtoms,renorm,llinearsoc,dim,outputunit
    use mod_system, only: n0sc1, n0sc2
    implicit none
    integer           :: AllocateStatus

    if(myrank_row==0) then
      allocate( currents       (7, n0sc1:n0sc2, Npl), &
                total_currents (7, n0sc1:n0sc2), &
                dc_currents    (3, Npl), STAT = AllocateStatus )
      if (AllocateStatus/=0) then
         write(outputunit,"('[allocate_currents] Not enough memory for: currents,total_currents,dc_currents')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      if(renorm) then
        allocate( rcurrents       (7, n0sc1:n0sc2, Npl), &
                  rtotal_currents (7, n0sc1:n0sc2), STAT = AllocateStatus )
        if (AllocateStatus/=0) then
          write(outputunit,"('[allocate_currents] Not enough memory for: rcurrents,rtotal_currents')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if
    allocate( ttchiorbiikl   (n0sc1:n0sc2, dimspinAtoms,4), &
              Lxttchiorbiikl (n0sc1:n0sc2, dimspinAtoms,4), &
              Lyttchiorbiikl (n0sc1:n0sc2, dimspinAtoms,4), &
              Lzttchiorbiikl (n0sc1:n0sc2, dimspinAtoms,4), STAT = AllocateStatus  )
    if (AllocateStatus/=0) then
      write(outputunit,"('[allocate_currents] Not enough memory for: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    if (.not. allocated(prefactor)) then
      allocate(prefactor(dim,dim), STAT = AllocateStatus  )
      if (AllocateStatus/=0) then
        write(outputunit,"('[allocate_currents] Not enough memory for: prefactor')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if
    if (.not. allocated(prefactorlsoc)) then
      if(llinearsoc) then
        allocate(prefactorlsoc(dim,dim), STAT = AllocateStatus  )
        if (AllocateStatus/=0) then
          write(outputunit,"('[allocate_currents] Not enough memory for: prefactorlsoc')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if

  end subroutine allocate_currents

  ! This subroutine allocates variables related to the current calculation
  subroutine deallocate_currents()
    use, intrinsic :: iso_fortran_env
    use mod_mpi_pars
    use mod_prefactors,only: prefactor,prefactorlsoc
    use mod_parameters, only: renorm
    implicit none

    if(myrank_row==0) then
      deallocate(currents,total_currents,dc_currents)
      if(renorm) deallocate(rcurrents,rtotal_currents)
    end if
    deallocate(ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)
    if (allocated(prefactor)) deallocate(prefactor)
    if (allocated(prefactorlsoc)) deallocate(prefactorlsoc)

  end subroutine deallocate_currents

  ! This subroutine opens and closes all the files needed for the currents
  subroutine openclose_currents_files(iflag)
    use mod_parameters, only: fieldpart, output
    use mod_mpi_pars
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: nkpt, n0sc1, n0sc2
    use electricfield, only: strElectricField
    implicit none

    character(len=500)  :: varm
    character(len=5)    :: folder(7)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,neighbor,iflag,err,errt=0


    folder(1) = "CC"
    folder(2) = "SC"
    folder(3) = "SC"
    folder(4) = "SC"
    folder(5) = "LC"
    folder(6) = "LC"
    folder(7) = "LC"
    if(lhfresponses) then
      do i=1,7
        folder(i) = trim(folder(i)) // "_HF"
      end do
    end if

    filename(1) = "Ich"
    filename(2) = "Isx"
    filename(3) = "Isy"
    filename(4) = "Isz"
    filename(5) = "Ilx"
    filename(6) = "Ily"
    filename(7) = "Ilz"

    if(iflag==0) then
      ! Header for currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/rprll',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do ; end do
      ! Header for total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 7000+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,'_neighbor=',i0,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/rprll',a,'_neighbor=',i0,'_total',a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#     energy    , amplitude of ',a,' , real part of ',a,' , imag part of ',a,' ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do
      ! Header for DC spin current
      do i=1,Npl ; do j=2,4
        iw = 8100+(i-1)*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'pumpdc_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(3)),filename(j),i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#      energy     ,    ',a,'pumpdc   ')") filename(j)
        close(unit=iw)
      end do ; end do

    else if(iflag==1) then
      ! Currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/rprll',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end if
      end do ; end do ; end do

      ! Total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 7000+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/prll',a,'_neighbor=',i0,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/rprll',a,'_neighbor=',i0,'_total',a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),filename(j),neighbor,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end if
      end do ; end do

      ! DC spin current
      do i=1,Npl ; do j=2,4
        iw = 8100+(i-1)*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'pumpdc_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(3)),filename(j),i,trim(output%Energy),trim(output%info),trim(fieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
      end do ; end do

      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_currents_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
    else
      ! Closing currents per plane per neighbor files
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do ; end do

      ! Closing total currents for each neighbor direction files
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 7000+(neighbor-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do

      ! Closing DC spin currents
      do i=1,Npl ; do j=2,4
        iw = 8100+(i-1)*3+j-1
        close(unit=iw)
      end do ; end do
    end if

  end subroutine openclose_currents_files

  ! This subroutine write all the currents into files
  ! (already opened with openclose_currents_files(1))
  ! Some information is also written on the screen
  subroutine write_currents(e)
    use, intrinsic :: iso_fortran_env
    use mod_parameters, only: renorm,outputunit_loop,lwriteonscreen
    use mod_magnet, only: mvec_spherical,mtotal_spherical
    use mod_system, only: n0sc1, n0sc2
    implicit none
    integer  :: neighbor,i,j,iw
    real(dp),intent(in) :: e

    if(lwriteonscreen) write(outputunit_loop,"(' ################# Currents: #################')")
    do neighbor=n0sc1,n0sc2
      if(lwriteonscreen) then
        write(outputunit_loop,"('|--------- Neighbor: ',i0,' , Energy = ',es11.4,' ---------|')") neighbor,e

        write(outputunit_loop,"('     Ich = (',es16.9,') + i(',es16.9,')')") real(total_currents(1,neighbor)),dimag(total_currents(1,neighbor))
        write(outputunit_loop,"(' abs(Ich) = ',es16.9)") abs(total_currents(1,neighbor))
        write(outputunit_loop,"('atan(Ich) = ',es16.9)") atan2(dimag(total_currents(1,neighbor)),real(total_currents(1,neighbor)))

        write(outputunit_loop,"('     Isx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(2,neighbor)),dimag(total_currents(2,neighbor))
        write(outputunit_loop,"(' abs(Isx) = ',es16.9)") abs(total_currents(2,neighbor))
        write(outputunit_loop,"('atan(Isx) = ',es16.9)") atan2(dimag(total_currents(2,neighbor)),real(total_currents(2,neighbor)))

        write(outputunit_loop,"('     Isy  = (',es16.9,') + i(',es16.9,')')") real(total_currents(3,neighbor)),dimag(total_currents(3,neighbor))
        write(outputunit_loop,"(' abs(Isy) = ',es16.9)") abs(total_currents(3,neighbor))
        write(outputunit_loop,"('atan(Isy) = ',es16.9)") atan2(dimag(total_currents(3,neighbor)),real(total_currents(3,neighbor)))

        write(outputunit_loop,"('     Isz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(4,neighbor)),dimag(total_currents(4,neighbor))
        write(outputunit_loop,"(' abs(Isz) = ',es16.9)") abs(total_currents(4,neighbor))
        write(outputunit_loop,"('atan(Isz) = ',es16.9)") atan2(dimag(total_currents(4,neighbor)),real(total_currents(4,neighbor)))

        write(outputunit_loop,"('     Ilx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(5,neighbor)),dimag(total_currents(5,neighbor))
        write(outputunit_loop,"(' abs(Ilx) = ',es16.9)") abs(total_currents(5,neighbor))
        write(outputunit_loop,"('atan(Ilx) = ',es16.9)") atan2(dimag(total_currents(5,neighbor)),real(total_currents(5,neighbor)))

        write(outputunit_loop,"('     Ily  = (',es16.9,') + i(',es16.9,')')") real(total_currents(6,neighbor)),dimag(total_currents(6,neighbor))
        write(outputunit_loop,"(' abs(Ily) = ',es16.9)") abs(total_currents(6,neighbor))
        write(outputunit_loop,"('atan(Ily) = ',es16.9)") atan2(dimag(total_currents(6,neighbor)),real(total_currents(6,neighbor)))

        write(outputunit_loop,"('     Ilz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(7,neighbor)),dimag(total_currents(7,neighbor))
        write(outputunit_loop,"(' abs(Ilz) = ',es16.9)") abs(total_currents(7,neighbor))
        write(outputunit_loop,"('atan(Ilz) = ',es16.9)") atan2(dimag(total_currents(7,neighbor)),real(total_currents(7,neighbor)))
      end if

      ! Total currents for each neighbor direction
      ! Writing charge current
      iw = 7000+(neighbor-1)*7
      write(unit=iw+1,fmt="(9(es16.9,2x))") e , abs(total_currents(1,neighbor)) , real(total_currents(1,neighbor)) , dimag(total_currents(1,neighbor)) , atan2(dimag(total_currents(1,neighbor)),real(total_currents(1,neighbor))) , real(total_currents(1,neighbor))/abs(total_currents(1,neighbor)) , dimag(total_currents(1,neighbor))/abs(total_currents(1,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

      ! Writing x-component spin current
      write(unit=iw+2,fmt="(9(es16.9,2x))") e , abs(total_currents(2,neighbor)) , real(total_currents(2,neighbor)) , dimag(total_currents(2,neighbor)) , atan2(dimag(total_currents(2,neighbor)),real(total_currents(2,neighbor))) , real(total_currents(2,neighbor))/abs(total_currents(2,neighbor)) , dimag(total_currents(2,neighbor))/abs(total_currents(2,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing y-component spin current
      write(unit=iw+3,fmt="(9(es16.9,2x))") e , abs(total_currents(3,neighbor)) , real(total_currents(3,neighbor)) , dimag(total_currents(3,neighbor)) , atan2(dimag(total_currents(3,neighbor)),real(total_currents(3,neighbor))) , real(total_currents(3,neighbor))/abs(total_currents(3,neighbor)) , dimag(total_currents(3,neighbor))/abs(total_currents(3,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing z-component spin current
      write(unit=iw+4,fmt="(9(es16.9,2x))") e , abs(total_currents(4,neighbor)) , real(total_currents(4,neighbor)) , dimag(total_currents(4,neighbor)) , atan2(dimag(total_currents(4,neighbor)),real(total_currents(4,neighbor))) , real(total_currents(4,neighbor))/abs(total_currents(4,neighbor)) , dimag(total_currents(4,neighbor))/abs(total_currents(4,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

      ! Writing x-component orbital angular momentum current
      write(unit=iw+5,fmt="(9(es16.9,2x))") e , abs(total_currents(5,neighbor)) , real(total_currents(5,neighbor)) , dimag(total_currents(5,neighbor)) , atan2(dimag(total_currents(5,neighbor)),real(total_currents(5,neighbor))) , real(total_currents(5,neighbor))/abs(total_currents(5,neighbor)) , dimag(total_currents(5,neighbor))/abs(total_currents(5,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing y-component orbital angular momentum current
      write(unit=iw+6,fmt="(9(es16.9,2x))") e , abs(total_currents(6,neighbor)) , real(total_currents(6,neighbor)) , dimag(total_currents(6,neighbor)) , atan2(dimag(total_currents(6,neighbor)),real(total_currents(6,neighbor))) , real(total_currents(6,neighbor))/abs(total_currents(6,neighbor)) , dimag(total_currents(6,neighbor))/abs(total_currents(6,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing z-component orbital angular momentum current
      write(unit=iw+7,fmt="(9(es16.9,2x))") e , abs(total_currents(7,neighbor)) , real(total_currents(7,neighbor)) , dimag(total_currents(7,neighbor)) , atan2(dimag(total_currents(7,neighbor)),real(total_currents(7,neighbor))) , real(total_currents(7,neighbor))/abs(total_currents(7,neighbor)) , dimag(total_currents(7,neighbor))/abs(total_currents(7,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

      ! Writing renormalized currents
      if(renorm) then
        ! Writing renormalized charge current
        write(unit=iw+1001,fmt="(9(es16.9,2x))") e , abs(rtotal_currents(1,neighbor)) , real(rtotal_currents(1,neighbor)) , dimag(rtotal_currents(1,neighbor)) , atan2(dimag(rtotal_currents(1,neighbor)),real(rtotal_currents(1,neighbor))) , real(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor)) , dimag(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

        ! Writing renormalized x-component spin current
        write(unit=iw+1002,fmt="(9(es16.9,2x))") e , abs(rtotal_currents(2,neighbor)) , real(rtotal_currents(2,neighbor)) , dimag(rtotal_currents(2,neighbor)) , atan2(dimag(rtotal_currents(2,neighbor)),real(rtotal_currents(2,neighbor))) , real(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor)) , dimag(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing renormalized y-component spin current
        write(unit=iw+1003,fmt="(9(es16.9,2x))") e , abs(rtotal_currents(3,neighbor)) , real(rtotal_currents(3,neighbor)) , dimag(rtotal_currents(3,neighbor)) , atan2(dimag(rtotal_currents(3,neighbor)),real(rtotal_currents(3,neighbor))) , real(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor)) , dimag(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing renormalized z-component spin current
        write(unit=iw+1004,fmt="(9(es16.9,2x))") e , abs(rtotal_currents(4,neighbor)) , real(rtotal_currents(4,neighbor)) , dimag(rtotal_currents(4,neighbor)) , atan2(dimag(rtotal_currents(4,neighbor)),real(rtotal_currents(4,neighbor))) , real(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor)) , dimag(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

        ! Writing x-component orbital angular momentum current
        write(unit=iw+1005,fmt="(9(es16.9,2x))") e , abs(rtotal_currents(5,neighbor)) , real(rtotal_currents(5,neighbor)) , dimag(rtotal_currents(5,neighbor)) , atan2(dimag(rtotal_currents(5,neighbor)),real(rtotal_currents(5,neighbor))) , real(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor)) , dimag(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing y-component orbital angular momentum current
        write(unit=iw+1006,fmt="(9(es16.9,2x))") e , abs(rtotal_currents(6,neighbor)) , real(rtotal_currents(6,neighbor)) , dimag(rtotal_currents(6,neighbor)) , atan2(dimag(rtotal_currents(6,neighbor)),real(rtotal_currents(6,neighbor))) , real(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor)) , dimag(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing z-component orbital angular momentum current
        write(unit=iw+1007,fmt="(9(es16.9,2x))") e , abs(rtotal_currents(7,neighbor)) , real(rtotal_currents(7,neighbor)) , dimag(rtotal_currents(7,neighbor)) , atan2(dimag(rtotal_currents(7,neighbor)),real(rtotal_currents(7,neighbor))) , real(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor)) , dimag(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      end if

      ! Writing currents per plane
      do i=1,Npl
        ! Writing charge current
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7
        write(unit=iw+1,fmt="(9(es16.9,2x))") e , abs(currents(1,neighbor,i)) , real(currents(1,neighbor,i)) , dimag(currents(1,neighbor,i)) , atan2(dimag(currents(1,neighbor,i)),real(currents(1,neighbor,i))) , real(currents(1,neighbor,i))/abs(currents(1,neighbor,i)) , dimag(currents(1,neighbor,i))/abs(currents(1,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing x-component spin current
        write(unit=iw+2,fmt="(9(es16.9,2x))") e , abs(currents(2,neighbor,i)) , real(currents(2,neighbor,i)) , dimag(currents(2,neighbor,i)) , atan2(dimag(currents(2,neighbor,i)),real(currents(2,neighbor,i))) , real(currents(2,neighbor,i))/abs(currents(2,neighbor,i)) , dimag(currents(2,neighbor,i))/abs(currents(2,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing y-component spin current
        write(unit=iw+3,fmt="(9(es16.9,2x))") e , abs(currents(3,neighbor,i)) , real(currents(3,neighbor,i)) , dimag(currents(3,neighbor,i)) , atan2(dimag(currents(3,neighbor,i)),real(currents(3,neighbor,i))) , real(currents(3,neighbor,i))/abs(currents(3,neighbor,i)) , dimag(currents(3,neighbor,i))/abs(currents(3,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing z-component spin current
        write(unit=iw+4,fmt="(9(es16.9,2x))") e , abs(currents(4,neighbor,i)) , real(currents(4,neighbor,i)) , dimag(currents(4,neighbor,i)) , atan2(dimag(currents(4,neighbor,i)),real(currents(4,neighbor,i))) , real(currents(4,neighbor,i))/abs(currents(4,neighbor,i)) , dimag(currents(4,neighbor,i))/abs(currents(4,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing x-component orbital angular momentum current
        write(unit=iw+5,fmt="(9(es16.9,2x))") e , abs(currents(5,neighbor,i)) , real(currents(5,neighbor,i)) , dimag(currents(5,neighbor,i)) , atan2(dimag(currents(5,neighbor,i)),real(currents(5,neighbor,i))) , real(currents(5,neighbor,i))/abs(currents(5,neighbor,i)) , dimag(currents(5,neighbor,i))/abs(currents(5,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing y-component orbital angular momentum current
        write(unit=iw+6,fmt="(9(es16.9,2x))") e , abs(currents(6,neighbor,i)) , real(currents(6,neighbor,i)) , dimag(currents(6,neighbor,i)) , atan2(dimag(currents(6,neighbor,i)),real(currents(6,neighbor,i))) , real(currents(6,neighbor,i))/abs(currents(6,neighbor,i)) , dimag(currents(6,neighbor,i))/abs(currents(6,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing z-component orbital angular momentum current
        write(unit=iw+7,fmt="(9(es16.9,2x))") e , abs(currents(7,neighbor,i)) , real(currents(7,neighbor,i)) , dimag(currents(7,neighbor,i)) , atan2(dimag(currents(7,neighbor,i)),real(currents(7,neighbor,i))) , real(currents(7,neighbor,i))/abs(currents(7,neighbor,i)) , dimag(currents(7,neighbor,i))/abs(currents(7,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing renormalized currents
        if(renorm) then
          ! Writing renormalized charge current
          write(unit=iw+1001,fmt="(9(es16.9,2x))") e , abs(rcurrents(1,neighbor,i)) , real(rcurrents(1,neighbor,i)) , dimag(rcurrents(1,neighbor,i)) , atan2(dimag(rcurrents(1,neighbor,i)),real(rcurrents(1,neighbor,i))) , real(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i)) , dimag(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

          ! Writing renormalized x-component spin current
          write(unit=iw+1002,fmt="(9(es16.9,2x))") e , abs(rcurrents(2,neighbor,i)) , real(rcurrents(2,neighbor,i)) , dimag(rcurrents(2,neighbor,i)) , atan2(dimag(rcurrents(2,neighbor,i)),real(rcurrents(2,neighbor,i))) , real(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i)) , dimag(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing renormalized y-component spin current
          write(unit=iw+1003,fmt="(9(es16.9,2x))") e , abs(rcurrents(3,neighbor,i)) , real(rcurrents(3,neighbor,i)) , dimag(rcurrents(3,neighbor,i)) , atan2(dimag(rcurrents(3,neighbor,i)),real(rcurrents(3,neighbor,i))) , real(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i)) , dimag(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing renormalized z-component spin current
          write(unit=iw+1004,fmt="(9(es16.9,2x))") e , abs(rcurrents(4,neighbor,i)) , real(rcurrents(4,neighbor,i)) , dimag(rcurrents(4,neighbor,i)) , atan2(dimag(rcurrents(4,neighbor,i)),real(rcurrents(4,neighbor,i))) , real(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i)) , dimag(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

          ! Writing x-component orbital angular momentum current
          write(unit=iw+1005,fmt="(9(es16.9,2x))") e , abs(rcurrents(5,neighbor,i)) , real(rcurrents(5,neighbor,i)) , dimag(rcurrents(5,neighbor,i)) , atan2(dimag(rcurrents(5,neighbor,i)),real(rcurrents(5,neighbor,i))) , real(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i)) , dimag(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing y-component orbital angular momentum current
          write(unit=iw+1006,fmt="(9(es16.9,2x))") e , abs(rcurrents(6,neighbor,i)) , real(rcurrents(6,neighbor,i)) , dimag(rcurrents(6,neighbor,i)) , atan2(dimag(rcurrents(6,neighbor,i)),real(rcurrents(6,neighbor,i))) , real(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i)) , dimag(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing z-component orbital angular momentum current
          write(unit=iw+1007,fmt="(9(es16.9,2x))") e , abs(rcurrents(7,neighbor,i)) , real(rcurrents(7,neighbor,i)) , dimag(rcurrents(7,neighbor,i)) , atan2(dimag(rcurrents(7,neighbor,i)),real(rcurrents(7,neighbor,i))) , real(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i)) , dimag(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        end if
      end do
    end do

    ! Writing DC spin currents
    do i=1,Npl ; do j=1,3
      iw = 8100+(i-1)*3+j
      write(unit=iw,fmt="(4(es16.9,2x))") e , dc_currents(j,i) , mvec_spherical(2,i) , mvec_spherical(3,i)
    end do ; end do

  end subroutine write_currents

  ! This subroutine opens and closes all the files needed for the currents
  subroutine openclose_dc_currents_files(iflag)
    use mod_parameters, only: dcfieldpart, output
    use mod_mpi_pars
    use mod_SOC, only: SOCc, socpart
    use mod_system, only:nkpt, n0sc1, n0sc2
    use electricfield, only: strElectricField
    implicit none

    character(len=500)  :: varm
    character(len=5)    :: folder(7)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,neighbor,iflag,err,errt=0

    folder(1) = "CC"
    folder(2) = "SC"
    folder(3) = "SC"
    folder(4) = "SC"
    folder(5) = "LC"
    folder(6) = "LC"
    folder(7) = "LC"
    if(lhfresponses) then
      do i=1,7
        folder(i) = trim(folder(i)) // "_HF"
      end do
    end if

    filename(1) = "Ich"
    filename(2) = "Isx"
    filename(3) = "Isy"
    filename(4) = "Isz"
    filename(5) = "Ilx"
    filename(6) = "Ily"
    filename(7) = "Ilz"

    if(iflag==0) then
      ! Header for currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,'_',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,' real part of ',a,' ,  imag part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'rprll',a,'_',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#',a,' real part of ',a,' ,  imag part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do ; end do
      ! Header for total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 70000+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,'_',a,'_neighbor=',i0,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,' real part of ',a,' ,  imag part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'rprll',a,'_',a,'_neighbor=',i0,'_total',a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='replace', form='formatted')
          write(unit=iw, fmt="('#',a,' real part of ',a,' ,  imag part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do

      ! Header for DC spin current
      do i=1,Npl ; do j=2,4
        iw = 81000+(i-1)*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'pumpdc_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(3)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='replace', form='formatted')
        write(unit=iw, fmt="('#',a,'    ',a,'pumpdc    , mag angle theta , mag angle phi  ')") trim(dc_header),filename(j)
        close(unit=iw)
      end do ; end do

    else if(iflag==1) then
      ! Currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,'_',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'rprll',a,'_',a,'_neighbor=',i0,'_pos=',i0,a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end if
      end do ; end do ; end do

      ! Total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 70000+(neighbor-1)*7+j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'prll',a,'_',a,'_neighbor=',i0,'_total',a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'rprll',a,'_',a,'_neighbor=',i0,'_total',a,a,a,a,a,'_renormnb=',i0,a,'.dat')") SOCc,trim(strSites),trim(folder(j)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),neighbor,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),renormnb,trim(suffix)
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
        end if
      end do ; end do

      ! DC spin current
      do i=1,Npl ; do j=2,4
        iw = 81000+(i-1)*3+j-1
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,a,'pumpdc_',a,'_pos=',i0,a,a,a,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(3)),trim(dcprefix(kount)),filename(j),trim(dcfield(dcfield_dependence)),i,trim(output%Energy),trim(output%info),trim(dcfieldpart),trim(socpart),trim(strElectricField),trim(suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
      end do ; end do

      ! Stop if some file does not exist
      if(errt/=0) call abortProgram("[openclose_dc_currents_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_LINE('A') // trim(missing_files))
    else
      ! Closing currents per plane per neighbor files
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do ; end do

      ! Closing total currents for each neighbor direction files
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 70000+(neighbor-1)*7+j
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          close(unit=iw)
        end if
      end do ; end do

      ! Closing DC spin currents
      do i=1,Npl ; do j=2,4
        iw = 81000+(i-1)*3+j-1
        close(unit=iw)
      end do ; end do
    end if

  end subroutine openclose_dc_currents_files

  ! This subroutine write all the currents in the dc limit into files
  ! (already opened with openclose_dc_currents_files(1))
  ! Some information is also written on the screen

  ! NOTE: Usually, the dclimit is obtained from the imaginary part of the
  ! response function. But we have used the i of the current opecator to
  ! change it to the real part in such a way that the phase we calculate
  ! is directly related to the applied electric field.
  subroutine write_dc_currents()
    use, intrinsic :: iso_fortran_env
    use mod_parameters
    use mod_magnet, only: mvec_spherical,mtotal_spherical
    use mod_system, only: n0sc1, n0sc2
    implicit none
    integer  :: neighbor,i,j,iw

    if(lwriteonscreen) write(outputunit_loop,"(' ################# Currents: #################')")
    do neighbor=n0sc1,n0sc2
      if(lwriteonscreen) then
        write(outputunit_loop,"('|--------- Neighbor: ',i0,' , ',a,' = ',a,' ---------|')") neighbor,trim(dcfield(dcfield_dependence)),trim(dc_fields(hw_count))

        write(outputunit_loop,"('     Ich = (',es16.9,') + i(',es16.9,')')") real(total_currents(1,neighbor)),dimag(total_currents(1,neighbor))
        write(outputunit_loop,"(' abs(Ich) = ',es16.9)") abs(total_currents(1,neighbor))
        write(outputunit_loop,"('atan(Ich) = ',es16.9)") atan2(dimag(total_currents(1,neighbor)),real(total_currents(1,neighbor)))

        write(outputunit_loop,"('     Isx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(2,neighbor)),dimag(total_currents(2,neighbor))
        write(outputunit_loop,"(' abs(Isx) = ',es16.9)") abs(total_currents(2,neighbor))
        write(outputunit_loop,"('atan(Isx) = ',es16.9)") atan2(dimag(total_currents(2,neighbor)),real(total_currents(2,neighbor)))

        write(outputunit_loop,"('     Isy  = (',es16.9,') + i(',es16.9,')')") real(total_currents(3,neighbor)),dimag(total_currents(3,neighbor))
        write(outputunit_loop,"(' abs(Isy) = ',es16.9)") abs(total_currents(3,neighbor))
        write(outputunit_loop,"('atan(Isy) = ',es16.9)") atan2(dimag(total_currents(3,neighbor)),real(total_currents(3,neighbor)))

        write(outputunit_loop,"('     Isz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(4,neighbor)),dimag(total_currents(4,neighbor))
        write(outputunit_loop,"(' abs(Isz) = ',es16.9)") abs(total_currents(4,neighbor))
        write(outputunit_loop,"('atan(Isz) = ',es16.9)") atan2(dimag(total_currents(4,neighbor)),real(total_currents(4,neighbor)))

        write(outputunit_loop,"('     Ilx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(5,neighbor)),dimag(total_currents(5,neighbor))
        write(outputunit_loop,"(' abs(Ilx) = ',es16.9)") abs(total_currents(5,neighbor))
        write(outputunit_loop,"('atan(Ilx) = ',es16.9)") atan2(dimag(total_currents(5,neighbor)),real(total_currents(5,neighbor)))

        write(outputunit_loop,"('     Ily  = (',es16.9,') + i(',es16.9,')')") real(total_currents(6,neighbor)),dimag(total_currents(6,neighbor))
        write(outputunit_loop,"(' abs(Ily) = ',es16.9)") abs(total_currents(6,neighbor))
        write(outputunit_loop,"('atan(Ily) = ',es16.9)") atan2(dimag(total_currents(6,neighbor)),real(total_currents(6,neighbor)))

        write(outputunit_loop,"('     Ilz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(7,neighbor)),dimag(total_currents(7,neighbor))
        write(outputunit_loop,"(' abs(Ilz) = ',es16.9)") abs(total_currents(7,neighbor))
        write(outputunit_loop,"('atan(Ilz) = ',es16.9)") atan2(dimag(total_currents(7,neighbor)),real(total_currents(7,neighbor)))
      end if

      ! Total currents for each neighbor direction
      ! Writing charge current
      iw = 70000+(neighbor-1)*7
      write(unit=iw+1,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_currents(1,neighbor)) , dimag(total_currents(1,neighbor)) , atan2(dimag(total_currents(1,neighbor)),real(total_currents(1,neighbor))) , real(total_currents(1,neighbor))/abs(total_currents(1,neighbor)) , dimag(total_currents(1,neighbor))/abs(total_currents(1,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

      ! Writing x-component spin current
      write(unit=iw+2,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_currents(2,neighbor)) , dimag(total_currents(2,neighbor)) , atan2(dimag(total_currents(2,neighbor)),real(total_currents(2,neighbor))) , real(total_currents(2,neighbor))/abs(total_currents(2,neighbor)) , dimag(total_currents(2,neighbor))/abs(total_currents(2,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing y-component spin current
      write(unit=iw+3,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_currents(3,neighbor)) , dimag(total_currents(3,neighbor)) , atan2(dimag(total_currents(3,neighbor)),real(total_currents(3,neighbor))) , real(total_currents(3,neighbor))/abs(total_currents(3,neighbor)) , dimag(total_currents(3,neighbor))/abs(total_currents(3,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing z-component spin current
      write(unit=iw+4,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_currents(4,neighbor)) , dimag(total_currents(4,neighbor)) , atan2(dimag(total_currents(4,neighbor)),real(total_currents(4,neighbor))) , real(total_currents(4,neighbor))/abs(total_currents(4,neighbor)) , dimag(total_currents(4,neighbor))/abs(total_currents(4,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

      ! Writing x-component orbital angular momentum current
      write(unit=iw+5,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_currents(5,neighbor)) , dimag(total_currents(5,neighbor)) , atan2(dimag(total_currents(5,neighbor)),real(total_currents(5,neighbor))) , real(total_currents(5,neighbor))/abs(total_currents(5,neighbor)) , dimag(total_currents(5,neighbor))/abs(total_currents(5,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing y-component orbital angular momentum current
      write(unit=iw+6,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_currents(6,neighbor)) , dimag(total_currents(6,neighbor)) , atan2(dimag(total_currents(6,neighbor)),real(total_currents(6,neighbor))) , real(total_currents(6,neighbor))/abs(total_currents(6,neighbor)) , dimag(total_currents(6,neighbor))/abs(total_currents(6,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      ! Writing z-component orbital angular momentum current
      write(unit=iw+7,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(total_currents(7,neighbor)) , dimag(total_currents(7,neighbor)) , atan2(dimag(total_currents(7,neighbor)),real(total_currents(7,neighbor))) , real(total_currents(7,neighbor))/abs(total_currents(7,neighbor)) , dimag(total_currents(7,neighbor))/abs(total_currents(7,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

      ! Writing renormalized currents
      if(renorm) then
        ! Writing renormalized charge current
        write(unit=iw+1001,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rtotal_currents(1,neighbor)) , dimag(rtotal_currents(1,neighbor)) , atan2(dimag(rtotal_currents(1,neighbor)),real(rtotal_currents(1,neighbor))) , real(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor)) , dimag(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

        ! Writing renormalized x-component spin current
        write(unit=iw+1002,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rtotal_currents(2,neighbor)) , dimag(rtotal_currents(2,neighbor)) , atan2(dimag(rtotal_currents(2,neighbor)),real(rtotal_currents(2,neighbor))) , real(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor)) , dimag(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing renormalized y-component spin current
        write(unit=iw+1003,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rtotal_currents(3,neighbor)) , dimag(rtotal_currents(3,neighbor)) , atan2(dimag(rtotal_currents(3,neighbor)),real(rtotal_currents(3,neighbor))) , real(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor)) , dimag(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing renormalized z-component spin current
        write(unit=iw+1004,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rtotal_currents(4,neighbor)) , dimag(rtotal_currents(4,neighbor)) , atan2(dimag(rtotal_currents(4,neighbor)),real(rtotal_currents(4,neighbor))) , real(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor)) , dimag(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)

        ! Writing x-component orbital angular momentum current
        write(unit=iw+1005,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rtotal_currents(5,neighbor)) , dimag(rtotal_currents(5,neighbor)) , atan2(dimag(rtotal_currents(5,neighbor)),real(rtotal_currents(5,neighbor))) , real(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor)) , dimag(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing y-component orbital angular momentum current
        write(unit=iw+1006,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rtotal_currents(6,neighbor)) , dimag(rtotal_currents(6,neighbor)) , atan2(dimag(rtotal_currents(6,neighbor)),real(rtotal_currents(6,neighbor))) , real(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor)) , dimag(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
        ! Writing z-component orbital angular momentum current
        write(unit=iw+1007,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rtotal_currents(7,neighbor)) , dimag(rtotal_currents(7,neighbor)) , atan2(dimag(rtotal_currents(7,neighbor)),real(rtotal_currents(7,neighbor))) , real(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor)) , dimag(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor)) , mtotal_spherical(2) , mtotal_spherical(3)
      end if

      ! Writing currents per plane
      do i=1,Npl
        ! Writing charge current
        iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7
        write(unit=iw+1,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(currents(1,neighbor,i)) , dimag(currents(1,neighbor,i)) , atan2(dimag(currents(1,neighbor,i)),real(currents(1,neighbor,i))) , real(currents(1,neighbor,i))/abs(currents(1,neighbor,i)) , dimag(currents(1,neighbor,i))/abs(currents(1,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing x-component spin current
        write(unit=iw+2,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(currents(2,neighbor,i)) , dimag(currents(2,neighbor,i)) , atan2(dimag(currents(2,neighbor,i)),real(currents(2,neighbor,i))) , real(currents(2,neighbor,i))/abs(currents(2,neighbor,i)) , dimag(currents(2,neighbor,i))/abs(currents(2,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing y-component spin current
        write(unit=iw+3,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(currents(3,neighbor,i)) , dimag(currents(3,neighbor,i)) , atan2(dimag(currents(3,neighbor,i)),real(currents(3,neighbor,i))) , real(currents(3,neighbor,i))/abs(currents(3,neighbor,i)) , dimag(currents(3,neighbor,i))/abs(currents(3,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing z-component spin current
        write(unit=iw+4,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(currents(4,neighbor,i)) , dimag(currents(4,neighbor,i)) , atan2(dimag(currents(4,neighbor,i)),real(currents(4,neighbor,i))) , real(currents(4,neighbor,i))/abs(currents(4,neighbor,i)) , dimag(currents(4,neighbor,i))/abs(currents(4,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing x-component orbital angular momentum current
        write(unit=iw+5,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(currents(5,neighbor,i)) , dimag(currents(5,neighbor,i)) , atan2(dimag(currents(5,neighbor,i)),real(currents(5,neighbor,i))) , real(currents(5,neighbor,i))/abs(currents(5,neighbor,i)) , dimag(currents(5,neighbor,i))/abs(currents(5,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing y-component orbital angular momentum current
        write(unit=iw+6,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(currents(6,neighbor,i)) , dimag(currents(6,neighbor,i)) , atan2(dimag(currents(6,neighbor,i)),real(currents(6,neighbor,i))) , real(currents(6,neighbor,i))/abs(currents(6,neighbor,i)) , dimag(currents(6,neighbor,i))/abs(currents(6,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        ! Writing z-component orbital angular momentum current
        write(unit=iw+7,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(currents(7,neighbor,i)) , dimag(currents(7,neighbor,i)) , atan2(dimag(currents(7,neighbor,i)),real(currents(7,neighbor,i))) , real(currents(7,neighbor,i))/abs(currents(7,neighbor,i)) , dimag(currents(7,neighbor,i))/abs(currents(7,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

        ! Writing renormalized currents
        if(renorm) then
          ! Writing renormalized charge current
          write(unit=iw+1001,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rcurrents(1,neighbor,i)) , dimag(rcurrents(1,neighbor,i)) , atan2(dimag(rcurrents(1,neighbor,i)),real(rcurrents(1,neighbor,i))) , real(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i)) , dimag(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

          ! Writing renormalized x-component spin current
          write(unit=iw+1002,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rcurrents(2,neighbor,i)) , dimag(rcurrents(2,neighbor,i)) , atan2(dimag(rcurrents(2,neighbor,i)),real(rcurrents(2,neighbor,i))) , real(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i)) , dimag(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing renormalized y-component spin current
          write(unit=iw+1003,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rcurrents(3,neighbor,i)) , dimag(rcurrents(3,neighbor,i)) , atan2(dimag(rcurrents(3,neighbor,i)),real(rcurrents(3,neighbor,i))) , real(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i)) , dimag(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing renormalized z-component spin current
          write(unit=iw+1004,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rcurrents(4,neighbor,i)) , dimag(rcurrents(4,neighbor,i)) , atan2(dimag(rcurrents(4,neighbor,i)),real(rcurrents(4,neighbor,i))) , real(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i)) , dimag(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)

          ! Writing x-component orbital angular momentum current
          write(unit=iw+1005,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rcurrents(5,neighbor,i)) , dimag(rcurrents(5,neighbor,i)) , atan2(dimag(rcurrents(5,neighbor,i)),real(rcurrents(5,neighbor,i))) , real(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i)) , dimag(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing y-component orbital angular momentum current
          write(unit=iw+1006,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rcurrents(6,neighbor,i)) , dimag(rcurrents(6,neighbor,i)) , atan2(dimag(rcurrents(6,neighbor,i)),real(rcurrents(6,neighbor,i))) , real(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i)) , dimag(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
          ! Writing z-component orbital angular momentum current
          write(unit=iw+1007,fmt="(a,2x,7(es16.9,2x))") trim(dc_fields(hw_count)) , real(rcurrents(7,neighbor,i)) , dimag(rcurrents(7,neighbor,i)) , atan2(dimag(rcurrents(7,neighbor,i)),real(rcurrents(7,neighbor,i))) , real(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i)) , dimag(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i)) , mvec_spherical(2,i) , mvec_spherical(3,i)
        end if
      end do
    end do

    ! Writing DC spin currents
    do i=1,Npl ; do j=1,3
      iw = 81000+(i-1)*3+j
      write(unit=iw,fmt="(a,2x,3(es16.9,2x))") trim(dc_fields(hw_count)) , dc_currents(j,i) , mvec_spherical(2,i) , mvec_spherical(3,i)
    end do ; end do

  end subroutine write_dc_currents

  ! This subroutine sorts current files
  subroutine sort_currents()
    use, intrinsic :: iso_fortran_env
    use mod_parameters, only: Npl,renorm,itype
    use mod_tools, only: sort_file
    use mod_system, only: n0sc1, n0sc2
    implicit none
    integer  :: neighbor,i,j,iw,idc=1

    ! Opening current files
    if(itype==9) then
      idc=10
      call openclose_dc_currents_files(1)
    else
      call openclose_currents_files(1)
    end if

    do neighbor=n0sc1,n0sc2
      ! Sorting total current files
      iw = 7000*idc+(neighbor-1)*7
      do j=1,7
        call sort_file(iw+j)
      end do

      ! Sorting renormalized total currents
      if(renorm) then
        do j=1001,1007
          call sort_file(iw+j)
        end do
      end if

      ! Sorting currents per plane
      do i=1,Npl
        ! Sorting current files
        iw = 5000*idc+(i-1)*n0sc2*7+(neighbor-1)*7
        do j=1,7
          call sort_file(iw+j)
        end do

        ! Sorting renormalized currents
        if(renorm) then
          do j=1001,1007
            call sort_file(iw+j)
          end do
        end if
      end do
    end do

    ! Sorting DC spin currents
    do i=1,Npl ; do j=1,3
      iw = 8100*idc+(i-1)*3+j
      call sort_file(iw)
    end do ; end do

    ! Closing current files
    if(itype==9) then
      call openclose_dc_currents_files(2)
    else
      call openclose_currents_files(2)
    end if

  end subroutine sort_currents

end module mod_currents
