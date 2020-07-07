module mod_fermi_surface
  use mod_f90_kind, only: double
  implicit none
  logical      :: lfs_loop = .false.
  integer      :: fs_energy_npts = 0, fs_energy_npt1 = 1
  real(double) :: fs_energy_i, fs_energy_f, fs_energy_s
  real(double), allocatable :: fs_energies(:)

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
      if(abs(fs_energy_s) <= 1.d-15) fs_energy_npt1 = 1

      ! Creating list of energies
      do i = 1, fs_energy_npt1
        fs_energies(i) = fs_energy_i + (i-1)*fs_energy_s
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
    use mod_system,        only: s => sys
    use mod_tools,         only: rtos
    use mod_mpi_pars,      only: rField
    implicit none
    integer :: i

    if(rField == 0) &
      write(output%unit_loop,"('CALCULATING CHARGE AND SPIN DENSITY AT ENERGY SURFACE')")

    ! Set energy values for the loop
    call setFSloop()

    do i = 1,fs_energy_npt1
      if(rField == 0) &
        write(output%unit_loop,"(i0,' of ',i0,' iso-energy surfaces. Energy e = ',es8.1)") i,fs_energy_npt1,fs_energies(i)
      ! Setting filename
      if((fs_energies(i)-s%Ef>1.d-6).or.(lfs_loop)) then
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
      if(rField == 0) call closeFSfiles()

    end do

  end subroutine fermi_surface


  !   Calculates iso-energy surface (e=Ef for Fermi surface)
  subroutine calculate_fermi_surface(e)
    use mod_f90_kind,      only: double
    use mod_constants,     only: pi, pauli_orb, cZero, cOne
    use mod_parameters,    only: nOrb, nOrb2, eta
    use mod_SOC,           only: llinearsoc, llineargfsoc
    use mod_system,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    use mod_mpi_pars
    implicit none
    real(double),intent(in)    :: e
    integer*8          :: iz
    integer*4          :: i,mu,nu,sigma
    real(double)       :: fs_atom(s%nAtoms,realBZ%nkpt,4),fs_orb(nOrb,realBZ%nkpt,4),fs_total(realBZ%nkpt,4)
    real(double)       :: kp(3)
    complex(double),dimension(nOrb2,nOrb2,s%nAtoms,s%nAtoms)    :: gf
    complex(double),dimension(nOrb2,nOrb2)    :: temp1,temp2,pauli_gf

    ! Allocating Brillouin Zone without parallelization
    call realBZ % setup_fraction(s,0, 1, FreqComm(1))

    !$omp parallel default(none) &
    !$omp& private(iz,kp,gf,i,mu,nu,sigma,temp1,temp2) &
    !$omp& shared(llineargfsoc,llinearsoc,s,nOrb,nOrb2,realBZ,e,eta,pauli_orb,pauli_gf,fs_atom,fs_orb,fs_total)

    fs_atom  = 0.d0
    fs_orb   = 0.d0
    fs_total = 0.d0

    !$omp do
    do iz = 1, realBZ%nkpt
      kp = realBZ%kp(1:3,iz)
      ! Green function on energy Ef + ieta, and wave vector kp
      if((llinearsoc).or.(llineargfsoc)) then
        call greenlineargfsoc(e,eta,s,kp,gf)
      else
        call green(e,eta,s,kp,gf)
      end if

      do sigma=1,4
        do i=1,s%nAtoms
           if(sigma==1) then
             pauli_gf = gf(:,:,i,i)
           else
             temp1 = pauli_orb(sigma-1,:,:)
             temp2 = gf(:,:,i,i)
             call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp1,nOrb2,temp2,nOrb2,cZero,pauli_gf,nOrb2)
           end if
           do mu=1,nOrb
             nu = mu + nOrb
             fs_orb(mu,iz,sigma) = fs_orb(mu,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
             fs_atom(i,iz,sigma) = fs_atom(i,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
           end do
         end do
      end do

      do mu=1,nOrb      
        fs_total(iz,:) = fs_total(iz,:) + fs_orb(mu,iz,:)
      end do

    end do
    !$omp end do
    !$omp end parallel

    if(rField == 0) call writeFS(fs_atom,fs_orb,fs_total)

  end subroutine calculate_fermi_surface

  subroutine createFSfiles()
    use mod_parameters,    only: output
    use mod_System,        only: s => sys
    implicit none
    character(len=400) :: varm
    integer            :: i

    do i=1,s%nAtoms
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'atom',i0,a,a,a,'.dat')") &
       output%SOCchar,trim(output%Sites),trim(epart),i,trim(output%info),trim(output%BField),trim(output%SOC)
      open (unit=17+i, file=varm,status='replace', form='formatted')
      write(unit=17+i, fmt="('#      kx       ,        ky       ,        kz       ,      charge     ,      spin x     ,      spin y     ,      spin z     ')")
      close(unit=17+i)
    end do
    ! Orbital-depedent
    write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,a,a,a,a,'.dat')") &
     output%SOCchar,trim(output%Sites),trim(epart),trim(filename(1)),trim(output%info),trim(output%BField),trim(output%SOC)
    open (unit=96, file=varm,status='replace', form='formatted')
    write(unit=96, fmt="('#      kx       ,        ky       ,        kz       ,  charge(*nOrb)  ,  spin x(*nOrb)  ,  spin y(*nOrb)  ,  spin z(*nOrb)  ')")
    close(unit=96)
    ! Total
    write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,a,a,a,a,'.dat')") &
     output%SOCchar,trim(output%Sites),trim(epart),trim(filename(2)),trim(output%info),trim(output%BField),trim(output%SOC)
    open (unit=97, file=varm,status='replace', form='formatted')
    write(unit=97, fmt="('#      kx       ,        ky       ,        kz       ,      charge     ,      spin x     ,      spin y     ,      spin z     ')")
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
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'atom',i0,a,a,a,'.dat')") &
       output%SOCchar,trim(output%Sites),trim(epart),i,trim(output%info),trim(output%BField),trim(output%SOC)
      open (unit=17+i, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
    end do
    do j = 1, size(filename)
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,a,a,a,a,'.dat')") &
       output%SOCchar,trim(output%Sites),trim(epart),trim(filename(j)),trim(output%info),trim(output%BField),trim(output%SOC)
      open (unit=95+j, file=varm, status='old', position='append', form='formatted', iostat=err)
      errt = errt + err
      if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
    end do

    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openFSfiles] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

  end subroutine openFSfiles

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
    use mod_parameters,    only: nOrb
    use mod_System,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    implicit none
    real(double), intent(in) :: fs_atom(s%nAtoms,realBZ%nkpt,4),fs_orb(nOrb,realBZ%nkpt,4),fs_total(realBZ%nkpt,4)
    character(len=30) :: formatvar
    real(double)      :: kp(3)
    integer*8         :: iz
    integer*4         :: i, mu, sigma

    write(formatvar,fmt="(a,i0,a)") '(',nOrb*4+3,'(es16.9,2x))'

    do iz = 1, realBZ%nkpt
      kp = realBZ%kp(1:3,iz)
      do i = 1, s%nAtoms
        write(unit=17+i,fmt="(7(es16.9,2x))") kp(1), kp(2), kp(3), (fs_atom(i,iz,sigma),sigma=1,4)
      end do

      write(unit=96,fmt=formatvar) kp(1), kp(2), kp(3),( (fs_orb(mu,iz,sigma),mu=1,nOrb),sigma=1,4)
      write(unit=97,fmt="(7(es16.9,2x))") kp(1), kp(2), kp(3),(fs_total(iz,sigma),sigma=1,4)
    end do
  end subroutine writeFS



end module mod_fermi_surface