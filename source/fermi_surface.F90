!   Calculates iso-energy surface (e=Ef for Fermi surface)
subroutine fermi_surface(e)
  use mod_f90_kind, only: double
  use mod_constants, only: pi, pauli_orb, cZero, cOne
  use mod_parameters, only: outputunit_loop, lverbose, eta, Ef, strSites, Utype, fieldpart
  use mod_SOC, only: SOCc, socpart, llinearsoc, llineargfsoc
  use mod_system, only: s => sys
  use TightBinding, only: nOrb,nOrb2
  use mod_progress
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer          :: nthreads,mythread
  character(len=400) :: varm
  character(len=50)  :: epart
  integer            :: i,mu,nu,iz,sigma
  real(double)       :: fs_layer(s%nAtoms,s%nkpt,4),fs_orb(3,s%nkpt,4),fs_total(s%nkpt,4)
  real(double)       :: kp(3)
  real(double),intent(in)    :: e
  complex(double),dimension(nOrb2,nOrb2,s%nAtoms,s%nAtoms)    :: gf
  complex(double),dimension(nOrb2,nOrb2)    :: temp1,temp2,pauli_gf

  write(outputunit_loop,"('CALCULATING CHARGE AND SPIN DENSITY AT FERMI SURFACE')")

  ! Opening files for writing
  if(abs(e-ef)>1.d-6) then
    write(epart,fmt="('e=',es8.1,'_')") e
  else
    write(epart,fmt="('fs_')")
  end if
  do i=1,s%nAtoms
    write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'layer',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(strSites),trim(epart),i,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=17+(mpitag-1)*s%nAtoms+i, file=varm,status='replace')
  end do
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'s_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(strSites),trim(epart),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=96+mpitag, file=varm,status='replace')
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'p_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(strSites),trim(epart),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=97+mpitag, file=varm,status='replace')
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'d_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(strSites),trim(epart),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=98+mpitag, file=varm,status='replace')
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'total_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(strSites),trim(epart),s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=99+mpitag, file=varm,status='replace')

  fs_layer = 0.d0
  fs_orb   = 0.d0
  fs_total = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,mu,nu,sigma,temp1,temp2) &
!$omp& shared(lverbose,llineargfsoc,llinearsoc,s,e,Ef,eta,nthreads,pauli_orb,pauli_gf,fs_layer,fs_orb,fs_total,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if(mythread==0) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[fermi_surface] Number of threads: ',i0)") nthreads
!$  end if

!$omp do
  fermi_surface_kpoints: do iz=1,s%nkpt
!$  if((mythread==0)) then
    if(lverbose) call progress_bar(outputunit_loop,"kpoints",iz,s%nkpt)
!$   end if

    kp = s%kbz(:,iz)

    ! Green function on energy Ef + ieta, and wave vector kp
    if((llinearsoc).or.(llineargfsoc)) then
      call greenlineargfsoc(e,eta,kp,gf)
    else
      call green(e,eta,kp,gf)
    end if

    do sigma=1,4 ; do i=1,s%nAtoms
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
    end do ; end do

    fs_total(iz,:) = fs_orb(1,iz,:) + fs_orb(2,iz,:) + fs_orb(3,iz,:)

  end do fermi_surface_kpoints
!$omp end do
!$omp end parallel

  ! Writing on files
  writing_fermi_surface: do iz=1,s%nkpt
    write_plane_loop_fs: do i=1,s%nAtoms
      write(unit=17+(mpitag-1)*s%nAtoms+i,fmt="(6(es16.9,2x))") s%kbz(1,iz),s%kbz(2,iz), s%kbz(3,iz), (fs_layer(i,iz,sigma),sigma=1,4)
    end do write_plane_loop_fs

    write(unit=96+mpitag,fmt="(6(es16.9,2x))") s%kbz(1,iz),s%kbz(2,iz),s%kbz(3,iz),(fs_orb(1,iz,sigma),sigma=1,4)
    write(unit=97+mpitag,fmt="(6(es16.9,2x))") s%kbz(1,iz),s%kbz(2,iz),s%kbz(3,iz),(fs_orb(2,iz,sigma),sigma=1,4)
    write(unit=98+mpitag,fmt="(6(es16.9,2x))") s%kbz(1,iz),s%kbz(2,iz),s%kbz(3,iz),(fs_orb(3,iz,sigma),sigma=1,4)
    write(unit=99+mpitag,fmt="(6(es16.9,2x))") s%kbz(1,iz),s%kbz(2,iz),s%kbz(3,iz),(fs_total(iz,sigma),sigma=1,4)
  end do writing_fermi_surface
  ! Closing files
  do i=1,s%nAtoms
    close (17+(mpitag-1)*s%nAtoms+i)
  end do
  close(96+mpitag)
  close(97+mpitag)
  close(98+mpitag)
  close(99+mpitag)

  return
end subroutine fermi_surface
