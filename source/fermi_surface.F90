!   Calculates iso-energy surface (e=Ef for Fermi surface)
subroutine fermi_surface(e)
  use mod_f90_kind,      only: double
  use mod_constants,     only: pi, pauli_orb, cZero, cOne
  use mod_parameters,    only: output, eta
  use mod_SOC,           only: llinearsoc, llineargfsoc
  use mod_system,        only: s => sys
  use mod_BrillouinZone, only: realBZ
  use TightBinding,      only: nOrb,nOrb2
  use mod_mpi_pars
  implicit none
  character(len=400) :: varm
  character(len=50)  :: epart
  integer            :: i,mu,nu,iz,sigma
  real(double)       :: fs_layer(s%nAtoms,realBZ%workload,4),fs_orb(3,realBZ%workload,4),fs_total(realBZ%workload,4)
  real(double)       :: kp(3)
  real(double),intent(in)    :: e
  complex(double),dimension(nOrb2,nOrb2,s%nAtoms,s%nAtoms)    :: gf
  complex(double),dimension(nOrb2,nOrb2)    :: temp1,temp2,pauli_gf
  integer :: proc
  write(output%unit_loop,"('CALCULATING CHARGE AND SPIN DENSITY AT FERMI SURFACE')")

  ! Allocating Brillouin Zone
  call realBZ % setup_fraction(s,rField, sField, FieldComm)

  ! Opening files for writing
  if(abs(e-s%ef)>1.d-6) then
    write(epart,fmt="('e=',es8.1,'_')") e
  else
    write(epart,fmt="('fs_')")
  end if

  fs_layer = 0.d0
  fs_orb   = 0.d0
  fs_total = 0.d0

!$omp parallel default(none) &
!$omp& private(iz,kp,gf,i,mu,nu,sigma,temp1,temp2) &
!$omp& shared(llineargfsoc,llinearsoc,s,realBZ,e,eta,pauli_orb,pauli_gf,fs_layer,fs_orb,fs_total)

!$omp do
  do iz = 1, realBZ%workload
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
           if(mu==1) then
             fs_orb  (1,iz,sigma) = fs_orb  (1,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
             fs_layer(i,iz,sigma) = fs_layer(i,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
           else if ((mu>=2).and.(mu<=4)) then
             fs_orb  (2,iz,sigma) = fs_orb  (2,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
             fs_layer(i,iz,sigma) = fs_layer(i,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
           else if ((mu>=5).and.(mu<=9)) then
             fs_orb  (3,iz,sigma) = fs_orb  (3,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
             fs_layer(i,iz,sigma) = fs_layer(i,iz,sigma) - aimag(pauli_gf(mu,mu)+pauli_gf(nu,nu))/pi
           end if
         end do
       end do
    end do

    fs_total(iz,:) = fs_orb(1,iz,:) + fs_orb(2,iz,:) + fs_orb(3,iz,:)

  end do
!$omp end do
!$omp end parallel


  ! Writing on files
  do proc = 0, sFreq(1)-1
    if(proc == rFreq(1)) then
      do i=1,s%nAtoms
        write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'layer',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") output%SOCchar,trim( output%Sites),trim(epart),i,trim(output%info),trim( output%BField),trim(output%SOC)
        open (unit=17+i, file=varm,status='replace')
      end do
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'s_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')")  output%SOCchar,trim( output%Sites),trim(epart),trim(output%info),trim( output%BField),trim(output%SOC)
      open (unit=96, file=varm,status='replace')
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'p_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')")  output%SOCchar,trim( output%Sites),trim(epart),trim(output%info),trim( output%BField),trim(output%SOC)
      open (unit=97, file=varm,status='replace')
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'d_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')")  output%SOCchar,trim( output%Sites),trim(epart),trim(output%info),trim( output%BField),trim(output%SOC)
      open (unit=98, file=varm,status='replace')
      write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'total_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')")  output%SOCchar,trim( output%Sites),trim(epart),trim(output%info),trim( output%BField),trim(output%SOC)
      open (unit=99, file=varm,status='replace')

      do iz = 1, realBZ%workload
        kp = realBZ%kp(1:3,iz)
        do i = 1, s%nAtoms
          write(unit=17+i,fmt="(6(es16.9,2x))") kp(1), kp(2), kp(3), (fs_layer(i,iz,sigma),sigma=1,4)
        end do

        write(unit=96,fmt="(6(es16.9,2x))") kp(1), kp(2), kp(3),(fs_orb(1,iz,sigma),sigma=1,4)
        write(unit=97,fmt="(6(es16.9,2x))") kp(1), kp(2), kp(3),(fs_orb(2,iz,sigma),sigma=1,4)
        write(unit=98,fmt="(6(es16.9,2x))") kp(1), kp(2), kp(3),(fs_orb(3,iz,sigma),sigma=1,4)
        write(unit=99,fmt="(6(es16.9,2x))") kp(1), kp(2), kp(3),(fs_total(iz,sigma),sigma=1,4)
      end do
      ! Closing files
      do i=1,s%nAtoms
        close (17+i)
      end do
      close(96)
      close(97)
      close(98)
      close(99)
    end if
    call MPI_Barrier(FreqComm(1), ierr)
  end do
  return
end subroutine fermi_surface
