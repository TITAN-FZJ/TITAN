module mod_fermi_surface
  use mod_kind, only: dp
  implicit none
  logical  :: lfs_loop = .false.
  integer  :: fs_energy_npts = 0, fs_energy_npt1 = 1
  real(dp) :: fs_energy_i, fs_energy_f, fs_energy_s
  real(dp), allocatable :: fs_energies(:)

  character(len=4), dimension(2), private :: filename = ["orb","tot"]
  character(len=50)  :: epart

contains

  ! This subroutine sets up the loop over iso-energies surfaces
  subroutine setFSloop()
    use mod_system,        only: s => sys
    implicit none
    integer :: i

    ! Setting up iso-energy surface loop
    if(lfs_loop) then
      allocate(fs_energies(fs_energy_npt1))
      ! If 0 points is selected, use single energy
      if (fs_energy_npts==0) then
        fs_energy_f = fs_energy_i
        fs_energy_npts = 1
      end if
      ! Calculating energy-step size
      fs_energy_s = (fs_energy_f - fs_energy_i)/fs_energy_npts
      ! If 0.0, then final and initial are the same; use one point only
      if(abs(fs_energy_s) <= 1.e-15_dp) fs_energy_npt1 = 1

      ! Creating list of energies
      do i = 1, fs_energy_npt1
        fs_energies(i) = fs_energy_i + (i-1)*fs_energy_s
        if( abs(fs_energies(i))<1.e-12_dp ) fs_energies(i) = 0._dp
      end do
    else
      ! If no loop is given, calculate Fermi surface
      allocate(fs_energies(1))
      fs_energies(1) = s%Ef
    end if

  end subroutine setFSloop

  !   Calculates iso-energy surface (e=Ef for Fermi surface)
  subroutine fermi_surface()
    use mod_parameters,    only: output, laddresults
    use mod_System,        only: s => sys
    use mod_tools,         only: rtos
    use mod_mpi_pars,      only: rField
    implicit none
    integer :: i

    if(rField == 0) &
      write(output%unit_loop,"('CALCULATING CHARGE, SPIN AND ORBITAL DENSITY AT ENERGY SURFACE')")

    ! Set energy values for the loop
    call setFSloop()

    do i = 1,fs_energy_npt1
      if(rField == 0) &
        write(output%unit_loop,"(i0,' of ',i0,' iso-energy surfaces. Energy e = ',es8.1)") i,fs_energy_npt1,fs_energies(i)
      ! Setting filename
      if((fs_energies(i)-s%Ef>1.e-6_dp).or.(lfs_loop)) then
        write(epart,fmt="('isoe_e=',a,'_')") trim(rtos(fs_energies(i),"(es8.1)"))
      else
        write(epart,fmt="('fs_')")
      end if

      ! Opening files
      if(rField == 0) then
        if(.not.laddresults) call createFSfiles()
        call openFSfiles()
      end if

      call calculate_fermi_surface(fs_energies(i))

      ! Closing files
      if(rField == 0) then
        call sortFermiSurface()
        call closeFSfiles()
      end if

    end do

  end subroutine fermi_surface


  !   Calculates iso-energy surface (e=Ef for Fermi surface)
  subroutine calculate_fermi_surface(e)
    use mod_kind,          only: dp,int32,int64
    use AtomTypes,         only: default_nOrb
    use mod_constants,     only: pi,cZero,cOne
    use mod_parameters,    only: eta
    use mod_SOC,           only: llinearsoc,llineargfsoc
    use mod_system,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    use mod_hamiltonian,   only: hamilt_local
    use mod_greenfunction, only: calc_green
    use mod_mpi_pars,      only: rField,rFreq,sFreq,FreqComm,MPI_IN_PLACE,MPI_DOUBLE_PRECISION,MPI_SUM,ierr
    implicit none
    real(dp),intent(in) :: e
    integer(int64)      :: iz
    integer(int32)      :: i,mu,nu,mup,nup,orb,sigma,nOrb2
    real(dp)            :: fs_atom(s%nAtoms,realBZ%nkpt,7),fs_orb(default_nOrb,realBZ%nkpt,7),fs_total(realBZ%nkpt,7)
    real(dp)            :: kp(3)
    real(dp)            :: temp
    complex(dp)         :: templ
    complex(dp),dimension(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms)    :: gf
    complex(dp),dimension(s%nOrb2,s%nOrb2)    :: temp1,temp2,pauli_gf

    external :: zgemm,MPI_Reduce

    !  Distributing Brillouin Zone
    call realBZ%setup_fraction(s,rFreq(1),sFreq(1),FreqComm(1))

    ! Build local hamiltonian
    if((.not.llineargfsoc) .and. (.not.llinearsoc)) call hamilt_local(s)

    fs_atom  = 0._dp
    fs_orb   = 0._dp
    fs_total = 0._dp

    !$omp parallel do &
    !$omp& default(none) &
    !$omp& private(iz,kp,gf,i,mu,nu,mup,nup,orb,nOrb2,sigma,temp,temp1,temp2,templ,pauli_gf) &
    !$omp& shared(s,calc_green,realBZ,e,eta,fs_atom,fs_orb,fs_total)
    do iz = 1,realBZ%workload
      kp = realBZ%kp(1:3,iz)

      ! Green function on energy Ef + ieta, and wave vector kp
      call calc_green(e,eta,s,kp,gf)

      site_i: do i=1,s%nAtoms
        nOrb2 = s%Types(s%Basis(i)%Material)%nOrb2

        ! Spin and charge densities
        do sigma=1,4
          if(sigma==1) then
            pauli_gf(1:nOrb2,1:nOrb2) = gf(1:nOrb2,1:nOrb2,i,i)
          else
            temp1(1:nOrb2,1:nOrb2) = s%Types(s%Basis(i)%Material)%pauli_orb(sigma-1,1:nOrb2,1:nOrb2)
            temp2(1:nOrb2,1:nOrb2) = gf(1:nOrb2,1:nOrb2,i,i)

            call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp1,s%nOrb2,temp2,s%nOrb2,cZero,pauli_gf,s%nOrb2)
          end if
          do mu=1,s%Types(s%Basis(i)%Material)%nOrb
            nu = mu + s%Types(s%Basis(i)%Material)%nOrb
            temp = - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
            ! Storing orbital accumulation for this particular orbital type
            orb = s%Types(s%Basis(i)%Material)%Orbs(mu)
            fs_orb(orb,realBZ%first+iz-1,sigma) = fs_orb(orb,realBZ%first+iz-1,sigma) + temp
            ! Storing site accumulation
            fs_atom(i,realBZ%first+iz-1,sigma) = fs_atom(i,realBZ%first+iz-1,sigma) + temp
          end do
        end do

        ! OAM densities
        orb_mu: do mu=1,s%Types(s%Basis(i)%Material)%nOrb
          mup = mu+s%Types(s%Basis(i)%Material)%nOrb
          ! orbital type:
          orb = s%Types(s%Basis(i)%Material)%Orbs(mu)

          orb_nu: do nu=1,s%Types(s%Basis(i)%Material)%nOrb
            nup = nu+s%Types(s%Basis(i)%Material)%nOrb
            templ = (gf(nu,mu,i,i) + gf(nup,mup,i,i))/pi
            temp  = real( s%Types(s%Basis(i)%Material)%lvec(mu,nu,1)*templ )
            fs_orb(orb,realBZ%first+iz-1,5) = fs_orb(orb,realBZ%first+iz-1,5) + temp
            fs_atom(i,realBZ%first+iz-1,5) = fs_atom(i,realBZ%first+iz-1,5) + temp
            temp  = real( s%Types(s%Basis(i)%Material)%lvec(mu,nu,2)*templ )
            fs_orb(orb,realBZ%first+iz-1,6) = fs_orb(orb,realBZ%first+iz-1,6) + temp
            fs_atom(i,realBZ%first+iz-1,6) = fs_atom(i,realBZ%first+iz-1,6) + temp
            temp  = real( s%Types(s%Basis(i)%Material)%lvec(mu,nu,3)*templ )
            fs_orb(orb,realBZ%first+iz-1,7) = fs_orb(orb,realBZ%first+iz-1,7) + temp
            fs_atom(i,realBZ%first+iz-1,7) = fs_atom(i,realBZ%first+iz-1,7) + temp
          end do orb_nu
        end do orb_mu

        fs_total(realBZ%first+iz-1,:) = fs_total(realBZ%first+iz-1,:) + fs_atom(i,realBZ%first+iz-1,:)
      end do site_i

    end do
    !$omp end parallel do

    if(rField == 0) then
      call MPI_Reduce(MPI_IN_PLACE, fs_orb   , default_nOrb*realBZ%nkpt*7  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
      call MPI_Reduce(MPI_IN_PLACE, fs_atom  , s%nAtoms*realBZ%nkpt*7      , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
      call MPI_Reduce(MPI_IN_PLACE, fs_total , realBZ%nkpt*7               , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    else
      call MPI_Reduce(fs_orb      , fs_orb   , default_nOrb*realBZ%nkpt*7  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
      call MPI_Reduce(fs_atom     , fs_atom  , s%nAtoms*realBZ%nkpt*7      , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
      call MPI_Reduce(fs_total    , fs_total , realBZ%nkpt*7               , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    end if

    if(rField == 0) call writeFS(fs_atom,fs_orb,fs_total)

  end subroutine calculate_fermi_surface


  subroutine createFSfiles()
    use mod_parameters,    only: output
    use mod_System,        only: s => sys
    use AtomTypes,         only: default_Orbs,default_nOrb
    implicit none
    character(len=400) :: varm
    character(len=1200):: formatvar
    character(len=16)  :: temp
    character(len=6)   :: quantities(7)
    integer            :: i,j

    quantities(1) = 'charge'
    quantities(2:7) = ['ms_x','ms_y','ms_z','mo_x','mo_y','mo_z']
    do i=1,s%nAtoms
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'atom',i0,a,a,a,a,'.dat')") &
        output%SOCchar,trim(output%Sites),trim(epart),i,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=17+i, file=varm,status='replace', form='formatted')
      write(unit=17+i, fmt="('#      kx       ,        ky       ,        kz       ,      charge     ,        ms_x     ,        ms_y     ,        ms_z     ,        mo_x     ,        mo_y     ,        mo_z     ')")
      close(unit=17+i)
    end do
    ! Orbital-depedent
    write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,a,a,a,a,a,'.dat')") &
      output%SOCchar,trim(output%Sites),trim(epart),trim(filename(1)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=96, file=varm,status='replace', form='formatted')

    formatvar = "('#      kx       ,        ky       ,        kz       ,"
    do i = 1,size(quantities)
      do j = 1,default_nOrb
        temp = trim(quantities(i)) // "(" // trim(default_Orbs(j)) // "),"
        write(formatvar,fmt="(a,a18)") trim(formatvar), trim(temp)
      end do
    end do
    write(formatvar,fmt="(a,a)") trim( formatvar( :len(trim(formatvar))-1 ) ), "')"
    write(unit=96, fmt=formatvar)
    close(unit=96)
    ! Total
    write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,a,a,a,a,a,'.dat')") &
      output%SOCchar,trim(output%Sites),trim(epart),trim(filename(2)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
    open (unit=97, file=varm,status='replace', form='formatted')
    write(unit=97, fmt="('#      kx       ,        ky       ,        kz       ,      charge     ,        ms_x     ,        ms_y     ,        ms_z     ,        mo_x     ,        mo_y     ,        mo_z     ')")
    close(unit=97)

  end subroutine createFSfiles


  subroutine openFSfiles()
    use mod_parameters, only: output,missing_files
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=400) :: varm
    integer            :: i, j, err, errt=0

    ! Opening files for writing
    do i=1,s%nAtoms
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'atom',i0,a,a,a,a,'.dat')") &
        output%SOCchar,trim(output%Sites),trim(epart),i,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=17+i, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
    end do
    do j = 1, size(filename)
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,a,a,a,a,a,'.dat')") &
        output%SOCchar,trim(output%Sites),trim(epart),trim(filename(j)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      open (unit=95+j, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
    end do

    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openFSfiles] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

  end subroutine openFSfiles


  subroutine sortFermiSurface()
    use mod_parameters, only: output
    use mod_tools,      only: sort_command
    use mod_System,     only: s => sys
    implicit none
    integer :: i,j
    character(len=400) :: varm

    varm = ""
    do i=1,s%nAtoms
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'atom',i0,a,a,a,a,'.dat')") &
        output%SOCchar,trim(output%Sites),trim(epart),i,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      call sort_command(varm,3,[1,2,3])
    end do
    do j = 1, size(filename)
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,a,a,a,a,a,'.dat')") &
        output%SOCchar,trim(output%Sites),trim(epart),trim(filename(j)),trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
      call sort_command(varm,3,[1,2,3])
    end do
  end subroutine sortFermiSurface


  subroutine closeFSfiles()
    use mod_System, only: s => sys
    implicit none
    integer i,j
    do i=1,s%nAtoms
      close(unit=17+i)
    end do
    do j = 1, size(filename)
      close(unit=95+j)
    end do
  end subroutine closeFSfiles


  ! This subtoutine writes the iso-surfaces to the files
  subroutine writeFS(fs_atom,fs_orb,fs_total)
    use mod_kind,          only: int32, int64
    use AtomTypes,         only: default_nOrb
    use mod_System,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    use mod_mpi_pars,      only: FreqComm
    implicit none
    real(dp), intent(in) :: fs_atom(s%nAtoms,realBZ%nkpt,7),fs_orb(default_nOrb,realBZ%nkpt,7),fs_total(realBZ%nkpt,7)
    character(len=30)    :: formatvar
    real(dp)             :: kp(3)
    integer(int64)       :: iz
    integer(int32)       :: i, mu, sigma

    ! Generating all points at rank=0
    call realBZ % setup_fraction(s,0, 1, FreqComm(1))

    write(formatvar,fmt="(a,i0,a)") '(',default_nOrb*7+3,'(es16.9,2x))'

    do iz = 1, realBZ%nkpt
      kp = realBZ%kp(1:3,iz)
      do i = 1, s%nAtoms
        write(unit=17+i,fmt="(10(es16.9,2x))") kp(1), kp(2), kp(3), (fs_atom(i,iz,sigma),sigma=1,7)
      end do

      write(unit=96,fmt=formatvar) kp(1), kp(2), kp(3),( (fs_orb(mu,iz,sigma),mu=1,default_nOrb),sigma=1,7)
      write(unit=97,fmt="(10(es16.9,2x))") kp(1), kp(2), kp(3),(fs_total(iz,sigma),sigma=1,7)
    end do
  end subroutine writeFS



end module mod_fermi_surface
