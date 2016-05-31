module mod_currents
  use mod_f90_kind
  implicit none
  ! Currents, total currents and renormalized currents: Ich, Isx, Isy, Isz, Ilx, Ily,Ilz
  complex(double),allocatable   :: currents(:,:,:),total_currents(:,:),rcurrents(:,:,:),rtotal_currents(:,:)
  ! Full response functions
  complex(double), dimension(:,:,:), allocatable :: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl
contains

  ! This subroutine allocates variables related to the current calculation
  subroutine allocate_currents()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_prefactors,only: prefactor,prefactorlsoc
    use mod_parameters, only: n0sc1,n0sc2,Npl,dimsigmaNpl,renorm,llinearsoc,dim
    implicit none
    integer           :: AllocateStatus

    if(myrank_row.eq.0) then
      allocate( currents(7,n0sc1:n0sc2,Npl),total_currents(7,n0sc1:n0sc2), STAT = AllocateStatus )
      if (AllocateStatus.ne.0) then
        if(myrank.eq.0) write(*,"('[allocate_currents] Not enough memory for: currents,total_currents')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
      if(renorm) then
        allocate( rcurrents(7,n0sc1:n0sc2,Npl),rtotal_currents(7,n0sc1:n0sc2), STAT = AllocateStatus )
        if (AllocateStatus.ne.0) then
          if(myrank.eq.0) write(*,"('[allocate_currents] Not enough memory for: rcurrents,rtotal_currents')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if
    allocate( ttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lxttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lyttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4),Lzttchiorbiikl(n0sc1:n0sc2,dimsigmaNpl,4), STAT = AllocateStatus  )
    if (AllocateStatus.ne.0) then
      write(*,"('[allocate_currents] Not enough memory for: ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl')")
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
    end if
    if (.not. allocated(prefactor)) then
      allocate(prefactor(dim,dim), STAT = AllocateStatus  )
      if (AllocateStatus.ne.0) then
        write(*,"('[allocate_currents] Not enough memory for: prefactor')")
        call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
      end if
    end if
    if (.not. allocated(prefactorlsoc)) then
      if(llinearsoc) then
        allocate(prefactorlsoc(dim,dim), STAT = AllocateStatus  )
        if (AllocateStatus.ne.0) then
          write(*,"('[allocate_currents] Not enough memory for: prefactorlsoc')")
          call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
        end if
      end if
    end if

    return
  end subroutine allocate_currents

  ! This subroutine allocates variables related to the current calculation
  subroutine deallocate_currents()
    use mod_f90_kind
    use mod_mpi_pars
    use mod_prefactors,only: prefactor,prefactorlsoc
    use mod_parameters, only: renorm
    implicit none

    if(myrank_row.eq.0) then
      deallocate(currents,total_currents)
      if(renorm) deallocate(rcurrents,rtotal_currents)
    end if
    deallocate(ttchiorbiikl,Lxttchiorbiikl,Lyttchiorbiikl,Lzttchiorbiikl)
    if (allocated(prefactor)) deallocate(prefactor)
    if (allocated(prefactorlsoc)) deallocate(prefactorlsoc)

    return
  end subroutine deallocate_currents

  ! This subroutine opens and closes all the files needed for the currents
  subroutine openclose_currents_files(iflag)
    use mod_parameters
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=1)    :: SOCc
    character(len=5)    :: folder(7)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,neighbor,iflag,err,errt=0

    fieldpart = ""
    socpart   = ""
    if(SOC) then
      if(llinearsoc) then
        SOCc = "L"
      else
        SOCc = "T"
      end if
      write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
    else
      SOCc = "F"
    end if
    if(FIELD) then
      write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hwa,hwt,hwp
      if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
    end if

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

    if(iflag.eq.0) then
      ! Header for currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/prll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/rprll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do ; end do
      ! Header for total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 7000+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/prll',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/rprll',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#   energy  ,  amplitude of ',a,'  ,  real part of ',a,'  ,  imaginary part of ',a,'  ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  ')") filename(j),filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do
    else if(iflag.eq.1) then
      ! Currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/prll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/rprll',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do ; end do

      ! Total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 7000+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/prll',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/rprll',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do

      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(*,"(a,i0,a)") "[openclose_currents_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        stop
      end if
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
    end if

    return
  end subroutine openclose_currents_files

  ! This subroutine write all the currents into files
  ! (already opened with openclose_currents_files(1))
  ! Some information is also written on the screen
  subroutine write_currents(e)
    use mod_f90_kind
    use mod_parameters, only: Npl,renorm,n0sc1,n0sc2
    implicit none
    integer  :: neighbor,i,iw
    real(double),intent(in) :: e

    write(*,"(' ################# Currents: #################')")
    do neighbor=n0sc1,n0sc2
      write(*,"('|--------- Neighbor: ',i0,' , Energy = ',es11.4,' ---------|')") neighbor,e

      write(*,"('     Ich = (',es16.9,') + i(',es16.9,')')") real(total_currents(1,neighbor)),aimag(total_currents(1,neighbor))
      write(*,"(' abs(Ich) = ',es16.9)") abs(total_currents(1,neighbor))
      write(*,"('atan(Ich) = ',es16.9)") atan2(aimag(total_currents(1,neighbor)),real(total_currents(1,neighbor)))

      write(*,"('     Isx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(2,neighbor)),aimag(total_currents(2,neighbor))
      write(*,"(' abs(Isx) = ',es16.9)") abs(total_currents(2,neighbor))
      write(*,"('atan(Isx) = ',es16.9)") atan2(aimag(total_currents(2,neighbor)),real(total_currents(2,neighbor)))

      write(*,"('     Isy  = (',es16.9,') + i(',es16.9,')')") real(total_currents(3,neighbor)),aimag(total_currents(3,neighbor))
      write(*,"(' abs(Isy) = ',es16.9)") abs(total_currents(3,neighbor))
      write(*,"('atan(Isy) = ',es16.9)") atan2(aimag(total_currents(3,neighbor)),real(total_currents(3,neighbor)))

      write(*,"('     Isz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(4,neighbor)),aimag(total_currents(4,neighbor))
      write(*,"(' abs(Isz) = ',es16.9)") abs(total_currents(4,neighbor))
      write(*,"('atan(Isz) = ',es16.9)") atan2(aimag(total_currents(4,neighbor)),real(total_currents(4,neighbor)))

      write(*,"('     Ilx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(5,neighbor)),aimag(total_currents(5,neighbor))
      write(*,"(' abs(Ilx) = ',es16.9)") abs(total_currents(5,neighbor))
      write(*,"('atan(Ilx) = ',es16.9)") atan2(aimag(total_currents(5,neighbor)),real(total_currents(5,neighbor)))

      write(*,"('     Ily  = (',es16.9,') + i(',es16.9,')')") real(total_currents(6,neighbor)),aimag(total_currents(6,neighbor))
      write(*,"(' abs(Ily) = ',es16.9)") abs(total_currents(6,neighbor))
      write(*,"('atan(Ily) = ',es16.9)") atan2(aimag(total_currents(6,neighbor)),real(total_currents(6,neighbor)))

      write(*,"('     Ilz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(7,neighbor)),aimag(total_currents(7,neighbor))
      write(*,"(' abs(Ilz) = ',es16.9)") abs(total_currents(7,neighbor))
      write(*,"('atan(Ilz) = ',es16.9)") atan2(aimag(total_currents(7,neighbor)),real(total_currents(7,neighbor)))

      ! Total currents for each neighbor direction
      ! Writing charge current
      iw = 7000+(neighbor-1)*7
      write(unit=iw+1,fmt="(7(es16.9,2x))") e,abs(total_currents(1,neighbor)),real(total_currents(1,neighbor)),aimag(total_currents(1,neighbor)),atan2(aimag(total_currents(1,neighbor)),real(total_currents(1,neighbor))),real(total_currents(1,neighbor))/abs(total_currents(1,neighbor)),aimag(total_currents(1,neighbor))/abs(total_currents(1,neighbor))

      ! Writing x-component spin current
      write(unit=iw+2,fmt="(7(es16.9,2x))") e,abs(total_currents(2,neighbor)),real(total_currents(2,neighbor)),aimag(total_currents(2,neighbor)),atan2(aimag(total_currents(2,neighbor)),real(total_currents(2,neighbor))),real(total_currents(2,neighbor))/abs(total_currents(2,neighbor)),aimag(total_currents(2,neighbor))/abs(total_currents(2,neighbor))
      ! Writing y-component spin current
      write(unit=iw+3,fmt="(7(es16.9,2x))") e,abs(total_currents(3,neighbor)),real(total_currents(3,neighbor)),aimag(total_currents(3,neighbor)),atan2(aimag(total_currents(3,neighbor)),real(total_currents(3,neighbor))),real(total_currents(3,neighbor))/abs(total_currents(3,neighbor)),aimag(total_currents(3,neighbor))/abs(total_currents(3,neighbor))
      ! Writing z-component spin current
      write(unit=iw+4,fmt="(7(es16.9,2x))") e,abs(total_currents(4,neighbor)),real(total_currents(4,neighbor)),aimag(total_currents(4,neighbor)),atan2(aimag(total_currents(4,neighbor)),real(total_currents(4,neighbor))),real(total_currents(4,neighbor))/abs(total_currents(4,neighbor)),aimag(total_currents(4,neighbor))/abs(total_currents(4,neighbor))

      ! Writing x-component orbital angular momentum current
      write(unit=iw+5,fmt="(7(es16.9,2x))") e,abs(total_currents(5,neighbor)),real(total_currents(5,neighbor)),aimag(total_currents(5,neighbor)),atan2(aimag(total_currents(5,neighbor)),real(total_currents(5,neighbor))),real(total_currents(5,neighbor))/abs(total_currents(5,neighbor)),aimag(total_currents(5,neighbor))/abs(total_currents(5,neighbor))
      ! Writing y-component orbital angular momentum current
      write(unit=iw+6,fmt="(7(es16.9,2x))") e,abs(total_currents(6,neighbor)),real(total_currents(6,neighbor)),aimag(total_currents(6,neighbor)),atan2(aimag(total_currents(6,neighbor)),real(total_currents(6,neighbor))),real(total_currents(6,neighbor))/abs(total_currents(6,neighbor)),aimag(total_currents(6,neighbor))/abs(total_currents(6,neighbor))
      ! Writing z-component orbital angular momentum current
      write(unit=iw+7,fmt="(7(es16.9,2x))") e,abs(total_currents(7,neighbor)),real(total_currents(7,neighbor)),aimag(total_currents(7,neighbor)),atan2(aimag(total_currents(7,neighbor)),real(total_currents(7,neighbor))),real(total_currents(7,neighbor))/abs(total_currents(7,neighbor)),aimag(total_currents(7,neighbor))/abs(total_currents(7,neighbor))

      ! Writing renormalized currents
      if(renorm) then
        ! Writing renormalized charge current
        write(unit=iw+1001,fmt="(7(es16.9,2x))") e,abs(rtotal_currents(1,neighbor)),real(rtotal_currents(1,neighbor)),aimag(rtotal_currents(1,neighbor)),atan2(aimag(rtotal_currents(1,neighbor)),real(rtotal_currents(1,neighbor))),real(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor)),aimag(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor))

        ! Writing renormalized x-component spin current
        write(unit=iw+1002,fmt="(7(es16.9,2x))") e,abs(rtotal_currents(2,neighbor)),real(rtotal_currents(2,neighbor)),aimag(rtotal_currents(2,neighbor)),atan2(aimag(rtotal_currents(2,neighbor)),real(rtotal_currents(2,neighbor))),real(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor)),aimag(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor))
        ! Writing renormalized y-component spin current
        write(unit=iw+1003,fmt="(7(es16.9,2x))") e,abs(rtotal_currents(3,neighbor)),real(rtotal_currents(3,neighbor)),aimag(rtotal_currents(3,neighbor)),atan2(aimag(rtotal_currents(3,neighbor)),real(rtotal_currents(3,neighbor))),real(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor)),aimag(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor))
        ! Writing renormalized z-component spin current
        write(unit=iw+1004,fmt="(7(es16.9,2x))") e,abs(rtotal_currents(4,neighbor)),real(rtotal_currents(4,neighbor)),aimag(rtotal_currents(4,neighbor)),atan2(aimag(rtotal_currents(4,neighbor)),real(rtotal_currents(4,neighbor))),real(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor)),aimag(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor))

        ! Writing x-component orbital angular momentum current
        write(unit=iw+1005,fmt="(7(es16.9,2x))") e,abs(rtotal_currents(5,neighbor)),real(rtotal_currents(5,neighbor)),aimag(rtotal_currents(5,neighbor)),atan2(aimag(rtotal_currents(5,neighbor)),real(rtotal_currents(5,neighbor))),real(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor)),aimag(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor))
        ! Writing y-component orbital angular momentum current
        write(unit=iw+1006,fmt="(7(es16.9,2x))") e,abs(rtotal_currents(6,neighbor)),real(rtotal_currents(6,neighbor)),aimag(rtotal_currents(6,neighbor)),atan2(aimag(rtotal_currents(6,neighbor)),real(rtotal_currents(6,neighbor))),real(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor)),aimag(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor))
        ! Writing z-component orbital angular momentum current
        write(unit=iw+1007,fmt="(7(es16.9,2x))") e,abs(rtotal_currents(7,neighbor)),real(rtotal_currents(7,neighbor)),aimag(rtotal_currents(7,neighbor)),atan2(aimag(rtotal_currents(7,neighbor)),real(rtotal_currents(7,neighbor))),real(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor)),aimag(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor))
      end if

      ! Writing currents per plane
      do i=1,Npl
        ! Writing charge current
        iw = 5000+(i-1)*n0sc2*7+(neighbor-1)*7
        write(unit=iw+1,fmt="(7(es16.9,2x))") e,abs(currents(1,neighbor,i)),real(currents(1,neighbor,i)),aimag(currents(1,neighbor,i)),atan2(aimag(currents(1,neighbor,i)),real(currents(1,neighbor,i))),real(currents(1,neighbor,i))/abs(currents(1,neighbor,i)),aimag(currents(1,neighbor,i))/abs(currents(1,neighbor,i))

        ! Writing x-component spin current
        write(unit=iw+2,fmt="(7(es16.9,2x))") e,abs(currents(2,neighbor,i)),real(currents(2,neighbor,i)),aimag(currents(2,neighbor,i)),atan2(aimag(currents(2,neighbor,i)),real(currents(2,neighbor,i))),real(currents(2,neighbor,i))/abs(currents(2,neighbor,i)),aimag(currents(2,neighbor,i))/abs(currents(2,neighbor,i))
        ! Writing y-component spin current
        write(unit=iw+3,fmt="(7(es16.9,2x))") e,abs(currents(3,neighbor,i)),real(currents(3,neighbor,i)),aimag(currents(3,neighbor,i)),atan2(aimag(currents(3,neighbor,i)),real(currents(3,neighbor,i))),real(currents(3,neighbor,i))/abs(currents(3,neighbor,i)),aimag(currents(3,neighbor,i))/abs(currents(3,neighbor,i))
        ! Writing z-component spin current
        write(unit=iw+4,fmt="(7(es16.9,2x))") e,abs(currents(4,neighbor,i)),real(currents(4,neighbor,i)),aimag(currents(4,neighbor,i)),atan2(aimag(currents(4,neighbor,i)),real(currents(4,neighbor,i))),real(currents(4,neighbor,i))/abs(currents(4,neighbor,i)),aimag(currents(4,neighbor,i))/abs(currents(4,neighbor,i))

        ! Writing x-component orbital angular momentum current
        write(unit=iw+5,fmt="(7(es16.9,2x))") e,abs(currents(5,neighbor,i)),real(currents(5,neighbor,i)),aimag(currents(5,neighbor,i)),atan2(aimag(currents(5,neighbor,i)),real(currents(5,neighbor,i))),real(currents(5,neighbor,i))/abs(currents(5,neighbor,i)),aimag(currents(5,neighbor,i))/abs(currents(5,neighbor,i))
        ! Writing y-component orbital angular momentum current
        write(unit=iw+6,fmt="(7(es16.9,2x))") e,abs(currents(6,neighbor,i)),real(currents(6,neighbor,i)),aimag(currents(6,neighbor,i)),atan2(aimag(currents(6,neighbor,i)),real(currents(6,neighbor,i))),real(currents(6,neighbor,i))/abs(currents(6,neighbor,i)),aimag(currents(6,neighbor,i))/abs(currents(6,neighbor,i))
        ! Writing z-component orbital angular momentum current
        write(unit=iw+7,fmt="(7(es16.9,2x))") e,abs(currents(7,neighbor,i)),real(currents(7,neighbor,i)),aimag(currents(7,neighbor,i)),atan2(aimag(currents(7,neighbor,i)),real(currents(7,neighbor,i))),real(currents(7,neighbor,i))/abs(currents(7,neighbor,i)),aimag(currents(7,neighbor,i))/abs(currents(7,neighbor,i))

        ! Writing renormalized currents
        if(renorm) then
          ! Writing renormalized charge current
          write(unit=iw+1001,fmt="(7(es16.9,2x))") e,abs(rcurrents(1,neighbor,i)),real(rcurrents(1,neighbor,i)),aimag(rcurrents(1,neighbor,i)),atan2(aimag(rcurrents(1,neighbor,i)),real(rcurrents(1,neighbor,i))),real(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i)),aimag(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i))

          ! Writing renormalized x-component spin current
          write(unit=iw+1002,fmt="(7(es16.9,2x))") e,abs(rcurrents(2,neighbor,i)),real(rcurrents(2,neighbor,i)),aimag(rcurrents(2,neighbor,i)),atan2(aimag(rcurrents(2,neighbor,i)),real(rcurrents(2,neighbor,i))),real(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i)),aimag(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i))
          ! Writing renormalized y-component spin current
          write(unit=iw+1003,fmt="(7(es16.9,2x))") e,abs(rcurrents(3,neighbor,i)),real(rcurrents(3,neighbor,i)),aimag(rcurrents(3,neighbor,i)),atan2(aimag(rcurrents(3,neighbor,i)),real(rcurrents(3,neighbor,i))),real(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i)),aimag(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i))
          ! Writing renormalized z-component spin current
          write(unit=iw+1004,fmt="(7(es16.9,2x))") e,abs(rcurrents(4,neighbor,i)),real(rcurrents(4,neighbor,i)),aimag(rcurrents(4,neighbor,i)),atan2(aimag(rcurrents(4,neighbor,i)),real(rcurrents(4,neighbor,i))),real(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i)),aimag(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i))

          ! Writing x-component orbital angular momentum current
          write(unit=iw+1005,fmt="(7(es16.9,2x))") e,abs(rcurrents(5,neighbor,i)),real(rcurrents(5,neighbor,i)),aimag(rcurrents(5,neighbor,i)),atan2(aimag(rcurrents(5,neighbor,i)),real(rcurrents(5,neighbor,i))),real(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i)),aimag(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i))
          ! Writing y-component orbital angular momentum current
          write(unit=iw+1006,fmt="(7(es16.9,2x))") e,abs(rcurrents(6,neighbor,i)),real(rcurrents(6,neighbor,i)),aimag(rcurrents(6,neighbor,i)),atan2(aimag(rcurrents(6,neighbor,i)),real(rcurrents(6,neighbor,i))),real(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i)),aimag(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i))
          ! Writing z-component orbital angular momentum current
          write(unit=iw+1007,fmt="(7(es16.9,2x))") e,abs(rcurrents(7,neighbor,i)),real(rcurrents(7,neighbor,i)),aimag(rcurrents(7,neighbor,i)),atan2(aimag(rcurrents(7,neighbor,i)),real(rcurrents(7,neighbor,i))),real(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i)),aimag(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i))
        end if
      end do
    end do

    return
  end subroutine write_currents


  ! This subroutine opens and closes all the files needed for the currents
  subroutine openclose_dclimit_currents_files(iflag)
    use mod_parameters
    implicit none

    character(len=500)  :: varm
    character(len=50)   :: fieldpart,socpart
    character(len=1)    :: SOCc
    character(len=5)    :: folder(7)
    character(len=3)    :: filename(7)
    integer :: i,j,iw,neighbor,iflag,err,errt=0

    fieldpart = ""
    socpart   = ""
    if(SOC) then
      if(llinearsoc) then
        SOCc = "L"
      else
        SOCc = "T"
      end if
      write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
    else
      SOCc = "F"
    end if
    if(FIELD) then
      if(dcfield_dependence.ne.1) write(fieldpart,"(a,'_hwa=',es9.2)") trim(fieldpart),hwa
      if(dcfield_dependence.ne.2) write(fieldpart,"(a,'_hwt=',f5.2)") trim(fieldpart),hwt
      if(dcfield_dependence.ne.3) write(fieldpart,"(a,'_hwp=',f5.2)") trim(fieldpart),hwp
      if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
    end if

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

    if(iflag.eq.0) then
      ! Header for currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcprll',a,'_',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   field  ,  imaginary part of ',a,'  ,   real part of ',a,'   ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcrprll',a,'_',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#   field  ,  imaginary part of ',a,'  ,   real part of ',a,'   ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do ; end do
      ! Header for total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 70000+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcprll',a,'_',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='unknown', form='formatted')
        write(unit=iw, fmt="('#   field  ,  imaginary part of ',a,'  ,   real part of ',a,'   ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") filename(j),filename(j),filename(j),filename(j),filename(j)
        close(unit=iw)

        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcrprll',a,'_',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='unknown', form='formatted')
          write(unit=iw, fmt="('#   field  ,  imaginary part of ',a,'  ,   real part of ',a,'   ,  phase of ',a,'  ,  cosine of ',a,'  ,  sine of ',a,'  , mag angle theta , mag angle phi  ')") filename(j),filename(j),filename(j),filename(j),filename(j)
          close(unit=iw)
        end if
      end do ; end do
    else if(iflag.eq.1) then
      ! Currents per plane per neighbor
      do i=1,Npl ; do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcprll',a,'_',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcrprll',a,'_',a,'_neighbor=',i0,'_pos=',i0,'_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,i,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do ; end do

      ! Total currents for each neighbor direction
      do neighbor=n0sc1,n0sc2 ; do j=1,7
        iw = 70000+(neighbor-1)*7+j
        write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcprll',a,'_',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(renorm) then
          iw = iw+1000
          write(varm,"('./results/SOC=',a1,'/Npl=',i0,'/',a,'/dcrprll',a,'_',a,'_neighbor=',i0,'_total_parts=',i0,'_parts3=',i0,'_ncp=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'_dirEfield=',a,'_renormnb=',i0,'.dat')") SOCc,Npl,trim(folder(j)),filename(j),dcfield(dcfield_dependence),neighbor,parts,parts3,ncp,eta,Utype,trim(fieldpart),trim(socpart),dirEfield,renormnb
          open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
        end if
      end do ; end do

      ! Stop if some file does not exist
      if(errt.ne.0) then
        write(*,"(a,i0,a)") "[openclose_dclimit_currents_files] Some file(s) do(es) not exist! Stopping before starting calculations..."
        stop
      end if
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
    end if

    return
  end subroutine openclose_dclimit_currents_files

  ! This subroutine write all the currents into files
  ! (already opened with openclose_currents_files(1))
  ! Some information is also written on the screen
  subroutine write_dclimit_currents(dcfield_variable)
    use mod_f90_kind
    use mod_constants, only: pi
    use mod_parameters, only: Npl,renorm,n0sc1,n0sc2,mmlayermag
    use mod_magnet, only: mtheta,mphi
    implicit none
    integer  :: neighbor,i,iw
    real(double),intent(in) :: dcfield_variable

    write(*,"(' ################# Currents: #################')")
    do neighbor=n0sc1,n0sc2
      write(*,"('|--------- Neighbor: ',i0,' , Field = ',es11.4,' ---------|')") neighbor,dcfield_variable

      write(*,"('     Ich = (',es16.9,') + i(',es16.9,')')") real(total_currents(1,neighbor)),aimag(total_currents(1,neighbor))
      write(*,"(' abs(Ich) = ',es16.9)") abs(total_currents(1,neighbor))
      write(*,"('atan(Ich) = ',es16.9)") atan2(aimag(total_currents(1,neighbor)),real(total_currents(1,neighbor)))

      write(*,"('     Isx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(2,neighbor)),aimag(total_currents(2,neighbor))
      write(*,"(' abs(Isx) = ',es16.9)") abs(total_currents(2,neighbor))
      write(*,"('atan(Isx) = ',es16.9)") atan2(aimag(total_currents(2,neighbor)),real(total_currents(2,neighbor)))

      write(*,"('     Isy  = (',es16.9,') + i(',es16.9,')')") real(total_currents(3,neighbor)),aimag(total_currents(3,neighbor))
      write(*,"(' abs(Isy) = ',es16.9)") abs(total_currents(3,neighbor))
      write(*,"('atan(Isy) = ',es16.9)") atan2(aimag(total_currents(3,neighbor)),real(total_currents(3,neighbor)))

      write(*,"('     Isz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(4,neighbor)),aimag(total_currents(4,neighbor))
      write(*,"(' abs(Isz) = ',es16.9)") abs(total_currents(4,neighbor))
      write(*,"('atan(Isz) = ',es16.9)") atan2(aimag(total_currents(4,neighbor)),real(total_currents(4,neighbor)))

      write(*,"('     Ilx  = (',es16.9,') + i(',es16.9,')')") real(total_currents(5,neighbor)),aimag(total_currents(5,neighbor))
      write(*,"(' abs(Ilx) = ',es16.9)") abs(total_currents(5,neighbor))
      write(*,"('atan(Ilx) = ',es16.9)") atan2(aimag(total_currents(5,neighbor)),real(total_currents(5,neighbor)))

      write(*,"('     Ily  = (',es16.9,') + i(',es16.9,')')") real(total_currents(6,neighbor)),aimag(total_currents(6,neighbor))
      write(*,"(' abs(Ily) = ',es16.9)") abs(total_currents(6,neighbor))
      write(*,"('atan(Ily) = ',es16.9)") atan2(aimag(total_currents(6,neighbor)),real(total_currents(6,neighbor)))

      write(*,"('     Ilz  = (',es16.9,') + i(',es16.9,')')") real(total_currents(7,neighbor)),aimag(total_currents(7,neighbor))
      write(*,"(' abs(Ilz) = ',es16.9)") abs(total_currents(7,neighbor))
      write(*,"('atan(Ilz) = ',es16.9)") atan2(aimag(total_currents(7,neighbor)),real(total_currents(7,neighbor)))

      ! Total currents for each neighbor direction
      ! Writing charge current
      iw = 70000+(neighbor-1)*7
      write(unit=iw+1,fmt="(8(es16.9,2x))") dcfield_variable,aimag(total_currents(1,neighbor)),real(total_currents(1,neighbor)),atan2(aimag(total_currents(1,neighbor)),real(total_currents(1,neighbor))),real(total_currents(1,neighbor))/abs(total_currents(1,neighbor)),aimag(total_currents(1,neighbor))/abs(total_currents(1,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi

      ! Writing x-component spin current
      write(unit=iw+2,fmt="(8(es16.9,2x))") dcfield_variable,aimag(total_currents(2,neighbor)),real(total_currents(2,neighbor)),atan2(aimag(total_currents(2,neighbor)),real(total_currents(2,neighbor))),real(total_currents(2,neighbor))/abs(total_currents(2,neighbor)),aimag(total_currents(2,neighbor))/abs(total_currents(2,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
      ! Writing y-component spin current
      write(unit=iw+3,fmt="(8(es16.9,2x))") dcfield_variable,aimag(total_currents(3,neighbor)),real(total_currents(3,neighbor)),atan2(aimag(total_currents(3,neighbor)),real(total_currents(3,neighbor))),real(total_currents(3,neighbor))/abs(total_currents(3,neighbor)),aimag(total_currents(3,neighbor))/abs(total_currents(3,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
      ! Writing z-component spin current
      write(unit=iw+4,fmt="(8(es16.9,2x))") dcfield_variable,aimag(total_currents(4,neighbor)),real(total_currents(4,neighbor)),atan2(aimag(total_currents(4,neighbor)),real(total_currents(4,neighbor))),real(total_currents(4,neighbor))/abs(total_currents(4,neighbor)),aimag(total_currents(4,neighbor))/abs(total_currents(4,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi

      ! Writing x-component orbital angular momentum current
      write(unit=iw+5,fmt="(8(es16.9,2x))") dcfield_variable,aimag(total_currents(5,neighbor)),real(total_currents(5,neighbor)),atan2(aimag(total_currents(5,neighbor)),real(total_currents(5,neighbor))),real(total_currents(5,neighbor))/abs(total_currents(5,neighbor)),aimag(total_currents(5,neighbor))/abs(total_currents(5,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
      ! Writing y-component orbital angular momentum current
      write(unit=iw+6,fmt="(8(es16.9,2x))") dcfield_variable,aimag(total_currents(6,neighbor)),real(total_currents(6,neighbor)),atan2(aimag(total_currents(6,neighbor)),real(total_currents(6,neighbor))),real(total_currents(6,neighbor))/abs(total_currents(6,neighbor)),aimag(total_currents(6,neighbor))/abs(total_currents(6,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
      ! Writing z-component orbital angular momentum current
      write(unit=iw+7,fmt="(8(es16.9,2x))") dcfield_variable,aimag(total_currents(7,neighbor)),real(total_currents(7,neighbor)),atan2(aimag(total_currents(7,neighbor)),real(total_currents(7,neighbor))),real(total_currents(7,neighbor))/abs(total_currents(7,neighbor)),aimag(total_currents(7,neighbor))/abs(total_currents(7,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi

      ! Writing renormalized currents
      if(renorm) then
        ! Writing renormalized charge current
        write(unit=iw+1001,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rtotal_currents(1,neighbor)),real(rtotal_currents(1,neighbor)),atan2(aimag(rtotal_currents(1,neighbor)),real(rtotal_currents(1,neighbor))),real(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor)),aimag(rtotal_currents(1,neighbor))/abs(rtotal_currents(1,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi

        ! Writing renormalized x-component spin current
        write(unit=iw+1002,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rtotal_currents(2,neighbor)),real(rtotal_currents(2,neighbor)),atan2(aimag(rtotal_currents(2,neighbor)),real(rtotal_currents(2,neighbor))),real(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor)),aimag(rtotal_currents(2,neighbor))/abs(rtotal_currents(2,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
        ! Writing renormalized y-component spin current
        write(unit=iw+1003,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rtotal_currents(3,neighbor)),real(rtotal_currents(3,neighbor)),atan2(aimag(rtotal_currents(3,neighbor)),real(rtotal_currents(3,neighbor))),real(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor)),aimag(rtotal_currents(3,neighbor))/abs(rtotal_currents(3,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
        ! Writing renormalized z-component spin current
        write(unit=iw+1004,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rtotal_currents(4,neighbor)),real(rtotal_currents(4,neighbor)),atan2(aimag(rtotal_currents(4,neighbor)),real(rtotal_currents(4,neighbor))),real(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor)),aimag(rtotal_currents(4,neighbor))/abs(rtotal_currents(4,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi

        ! Writing x-component orbital angular momentum current
        write(unit=iw+1005,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rtotal_currents(5,neighbor)),real(rtotal_currents(5,neighbor)),atan2(aimag(rtotal_currents(5,neighbor)),real(rtotal_currents(5,neighbor))),real(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor)),aimag(rtotal_currents(5,neighbor))/abs(rtotal_currents(5,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
        ! Writing y-component orbital angular momentum current
        write(unit=iw+1006,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rtotal_currents(6,neighbor)),real(rtotal_currents(6,neighbor)),atan2(aimag(rtotal_currents(6,neighbor)),real(rtotal_currents(6,neighbor))),real(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor)),aimag(rtotal_currents(6,neighbor))/abs(rtotal_currents(6,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
        ! Writing z-component orbital angular momentum current
        write(unit=iw+1007,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rtotal_currents(7,neighbor)),real(rtotal_currents(7,neighbor)),atan2(aimag(rtotal_currents(7,neighbor)),real(rtotal_currents(7,neighbor))),real(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor)),aimag(rtotal_currents(7,neighbor))/abs(rtotal_currents(7,neighbor)),mtheta(mmlayermag(1)-1)*pi,mphi(mmlayermag(1)-1)*pi
      end if

      ! Writing currents per plane
      do i=1,Npl
        ! Writing charge current
        iw = 50000+(i-1)*n0sc2*7+(neighbor-1)*7
        write(unit=iw+1,fmt="(8(es16.9,2x))") dcfield_variable,aimag(currents(1,neighbor,i)),real(currents(1,neighbor,i)),atan2(aimag(currents(1,neighbor,i)),real(currents(1,neighbor,i))),real(currents(1,neighbor,i))/abs(currents(1,neighbor,i)),aimag(currents(1,neighbor,i))/abs(currents(1,neighbor,i)),mtheta(i)*pi,mphi(i)*pi

        ! Writing x-component spin current
        write(unit=iw+2,fmt="(8(es16.9,2x))") dcfield_variable,aimag(currents(2,neighbor,i)),real(currents(2,neighbor,i)),atan2(aimag(currents(2,neighbor,i)),real(currents(2,neighbor,i))),real(currents(2,neighbor,i))/abs(currents(2,neighbor,i)),aimag(currents(2,neighbor,i))/abs(currents(2,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
        ! Writing y-component spin current
        write(unit=iw+3,fmt="(8(es16.9,2x))") dcfield_variable,aimag(currents(3,neighbor,i)),real(currents(3,neighbor,i)),atan2(aimag(currents(3,neighbor,i)),real(currents(3,neighbor,i))),real(currents(3,neighbor,i))/abs(currents(3,neighbor,i)),aimag(currents(3,neighbor,i))/abs(currents(3,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
        ! Writing z-component spin current
        write(unit=iw+4,fmt="(8(es16.9,2x))") dcfield_variable,aimag(currents(4,neighbor,i)),real(currents(4,neighbor,i)),atan2(aimag(currents(4,neighbor,i)),real(currents(4,neighbor,i))),real(currents(4,neighbor,i))/abs(currents(4,neighbor,i)),aimag(currents(4,neighbor,i))/abs(currents(4,neighbor,i)),mtheta(i)*pi,mphi(i)*pi

        ! Writing x-component orbital angular momentum current
        write(unit=iw+5,fmt="(8(es16.9,2x))") dcfield_variable,aimag(currents(5,neighbor,i)),real(currents(5,neighbor,i)),atan2(aimag(currents(5,neighbor,i)),real(currents(5,neighbor,i))),real(currents(5,neighbor,i))/abs(currents(5,neighbor,i)),aimag(currents(5,neighbor,i))/abs(currents(5,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
        ! Writing y-component orbital angular momentum current
        write(unit=iw+6,fmt="(8(es16.9,2x))") dcfield_variable,aimag(currents(6,neighbor,i)),real(currents(6,neighbor,i)),atan2(aimag(currents(6,neighbor,i)),real(currents(6,neighbor,i))),real(currents(6,neighbor,i))/abs(currents(6,neighbor,i)),aimag(currents(6,neighbor,i))/abs(currents(6,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
        ! Writing z-component orbital angular momentum current
        write(unit=iw+7,fmt="(8(es16.9,2x))") dcfield_variable,aimag(currents(7,neighbor,i)),real(currents(7,neighbor,i)),atan2(aimag(currents(7,neighbor,i)),real(currents(7,neighbor,i))),real(currents(7,neighbor,i))/abs(currents(7,neighbor,i)),aimag(currents(7,neighbor,i))/abs(currents(7,neighbor,i)),mtheta(i)*pi,mphi(i)*pi

        ! Writing renormalized currents
        if(renorm) then
          ! Writing renormalized charge current
          write(unit=iw+1001,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rcurrents(1,neighbor,i)),real(rcurrents(1,neighbor,i)),atan2(aimag(rcurrents(1,neighbor,i)),real(rcurrents(1,neighbor,i))),real(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i)),aimag(rcurrents(1,neighbor,i))/abs(rcurrents(1,neighbor,i)),mtheta(i)*pi,mphi(i)*pi

          ! Writing renormalized x-component spin current
          write(unit=iw+1002,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rcurrents(2,neighbor,i)),real(rcurrents(2,neighbor,i)),atan2(aimag(rcurrents(2,neighbor,i)),real(rcurrents(2,neighbor,i))),real(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i)),aimag(rcurrents(2,neighbor,i))/abs(rcurrents(2,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
          ! Writing renormalized y-component spin current
          write(unit=iw+1003,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rcurrents(3,neighbor,i)),real(rcurrents(3,neighbor,i)),atan2(aimag(rcurrents(3,neighbor,i)),real(rcurrents(3,neighbor,i))),real(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i)),aimag(rcurrents(3,neighbor,i))/abs(rcurrents(3,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
          ! Writing renormalized z-component spin current
          write(unit=iw+1004,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rcurrents(4,neighbor,i)),real(rcurrents(4,neighbor,i)),atan2(aimag(rcurrents(4,neighbor,i)),real(rcurrents(4,neighbor,i))),real(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i)),aimag(rcurrents(4,neighbor,i))/abs(rcurrents(4,neighbor,i)),mtheta(i)*pi,mphi(i)*pi

          ! Writing x-component orbital angular momentum current
          write(unit=iw+1005,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rcurrents(5,neighbor,i)),real(rcurrents(5,neighbor,i)),atan2(aimag(rcurrents(5,neighbor,i)),real(rcurrents(5,neighbor,i))),real(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i)),aimag(rcurrents(5,neighbor,i))/abs(rcurrents(5,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
          ! Writing y-component orbital angular momentum current
          write(unit=iw+1006,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rcurrents(6,neighbor,i)),real(rcurrents(6,neighbor,i)),atan2(aimag(rcurrents(6,neighbor,i)),real(rcurrents(6,neighbor,i))),real(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i)),aimag(rcurrents(6,neighbor,i))/abs(rcurrents(6,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
          ! Writing z-component orbital angular momentum current
          write(unit=iw+1007,fmt="(8(es16.9,2x))") dcfield_variable,aimag(rcurrents(7,neighbor,i)),real(rcurrents(7,neighbor,i)),atan2(aimag(rcurrents(7,neighbor,i)),real(rcurrents(7,neighbor,i))),real(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i)),aimag(rcurrents(7,neighbor,i))/abs(rcurrents(7,neighbor,i)),mtheta(i)*pi,mphi(i)*pi
        end if
      end do
    end do

    return
  end subroutine write_dclimit_currents
end module mod_currents