! This subroutine reads previous band-shifting results and writes new ones in a file
subroutine readwritesc(iflag,err)
  use mod_constants
  use mod_parameters
  use mod_magnet
  use mod_mpi_pars, only: myrank
  implicit none
  character(len=300)     :: file
  character(len=100)     :: fieldpart,socpart,folder,prefix
  integer,intent(inout)  :: iflag
  integer,intent(out)    :: err
  integer                :: i
  real(double),dimension(Npl)   :: mx,my

  if((trim(scfile).ne."").and.(iflag.eq.0)) then
    open(unit=99,file=scfile,status="old",iostat=err)
    if(err.ne.0) then
      if(myrank.eq.0) write(*,"('*** WARNING: Self-consistency file given on input file does not exist! Using default... ***')")
      scfile = " "
    end if
    close(99)
  end if

  folder    = "./results/selfconsistency/"
  prefix    = "selfconsistency_"
  fieldpart = ""
  socpart   = ""
  if(FIELD) then
    write(fieldpart,"('_hwx=',E8.1,'_hwy=',E8.1,'_hwz=',E8.1)") hwx,hwy,hwz
    if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
  end if
  if(SOC) then
    write(socpart,"('_magaxis=',A,'_socscale=',f5.2)") magaxis,socscale
    if((llineargfsoc).or.(llinearsoc)) socpart = trim(socpart) // "_linearsoc"
  end if

!   Reading previous results (mx, my, mz and eps1) from files (if available)
  if(iflag.eq.0) then
    if(trim(scfile).eq."") then ! If a filename is not given in inputcard, use the default one
      write(file,"(a,a,'Npl=',i0,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_ncp=',i0,'_eta=',e8.1,a,'.dat')") trim(folder),trim(prefix),Npl,dfttype,parts,Utype,trim(fieldpart),ncp,eta,trim(socpart)
      open(unit=99,file=file,status="old",iostat=err)
      if((err.eq.0).and.(myrank.eq.0)) write(*,"('[readwritesc] Self-consistency file already exists. Reading it now...')")
    else
      open(unit=99,file=scfile,status="old",iostat=err)
      if((err.eq.0).and.(myrank.eq.0)) then
        write(*,"('[readwritesc] Using another file as input for self-consistency:')")
        write(*,"(a)") scfile
      end if
    end if
    if(err.eq.0) then
      do i=1,Npl
        read(99,"(e21.11)") eps1(i)
        read(99,"(e21.11)") mx(i)
        read(99,"(e21.11)") my(i)
        read(99,"(e21.11)") mz(i)
      end do
      mp  = mx + zi*my
      do i=1,Npl
        hdel(i)   = 0.5d0*U(i+1)*mz(i)
        hdelp(i)  = 0.5d0*U(i+1)*mp(i)
      end do
      hdelm = conjg(hdelp)
      iflag   = 1 ! something was read
    else
      ! If file does not exist, try to read for parts-1
      close(99)
      write(file,"(a,a,'Npl=',i0,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_ncp=',i0,'_eta=',e8.1,a,'.dat')") trim(folder),trim(prefix),Npl,dfttype,parts-1,Utype,trim(fieldpart),ncp,eta,trim(socpart)
      open(unit=99,file=file,status="old",iostat=err)
      if(err.eq.0) then
        if(myrank.eq.0) then
          write(*,"('[readwritesc] Self-consistency file does not exist. Reading results for parts-1 now...')")
          write(*,"('[readwritesc] Updating values obtained for parts-1...')")
        end if
        do i=1,Npl
          read(99,"(e21.11)") eps1(i)
          read(99,"(e21.11)") mx(i)
          read(99,"(e21.11)") my(i)
          read(99,"(e21.11)") mz(i)
        end do
        mp  = mx + zi*my
        do i=1,Npl
          hdel(i)   = 0.5d0*U(i+1)*mz(i)
          hdelp(i)  = 0.5d0*U(i+1)*mp(i)
        end do
        hdelm = conjg(hdelp)
        iflag = 1 ! something was read
        err = 1   ! but not the same parameters
      end if
    end if
    close(99)
!   Writing new results (mx, my, mz and eps1) and mz to file
  else
    if(myrank.eq.0) write(*,"('[readwritesc] Writing new eps1, mx, my and mz to file...')")
    write(file,"(a,a,'Npl=',i0,'_dfttype=',a,'_parts=',i0,'_Utype=',i0,a,'_ncp=',i0,'_eta=',e8.1,a,'.dat')") trim(folder),trim(prefix),Npl,dfttype,parts,Utype,trim(fieldpart),ncp,eta,trim(socpart)
    open (unit=99,status="unknown",file=file)
    do i=1,Npl
      write(99,"(e21.11,2x,'! eps1')") eps1(i)
      write(99,"(e21.11,2x,'! mx')") real(mp(i))
      write(99,"(e21.11,2x,'! my')") aimag(mp(i))
      write(99,"(e21.11,2x,'! mz')") mz(i)
    end do
    close(99)
  end if
  return
end subroutine readwritesc