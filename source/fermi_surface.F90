!   Calculates iso-energy surface (e=Ef for Fermi surface)
subroutine fermi_surface(e)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_system, only: nkpt, kbz
  use mod_progress
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer          :: nthreads,mythread
  character(len=400) :: varm
  character(len=50)  :: fieldpart,socpart,epart
  character(len=1)   :: SOCc
  integer            :: i,mu,nu,iz,sigma
  real(double)       :: fs_layer(Npl,nkpt,4),fs_orb(3,nkpt,4),fs_total(nkpt,4)
  real(double)       :: kp(3)
  real(double),intent(in)    :: e
  complex(double),dimension(Npl,Npl,18,18)    :: gf
  complex(double),dimension(18,18)    :: temp1,temp2,pauli_gf

  write(outputunit_loop,"('CALCULATING CHARGE AND SPIN DENSITY AT FERMI SURFACE')")

  ! Opening files for writing
  fieldpart = ""
  socpart   = ""
  if(SOC) then
    if((llinearsoc).or.(llineargfsoc)) then
      SOCc = "L"
    else
      SOCc = "T"
    end if
    if(abs(socscale-1.d0)>1.d-6) write(socpart,"('_socscale=',f5.2)") socscale
  else
    SOCc = "F"
  end if
  if(lfield) then
    write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
    if(ltesla)    fieldpart = trim(fieldpart) // "_tesla"
    if(lnolb)     fieldpart = trim(fieldpart) // "_nolb"
    if(lhwscale)  fieldpart = trim(fieldpart) // "_hwscale"
    if(lhwrotate) fieldpart = trim(fieldpart) // "_hwrotate"
  end if
  if(abs(e-ef)>1.d-6) then
    write(epart,fmt="('e=',es8.1,'_')") e
  else
    write(epart,fmt="('fs_')")
  end if
  do i=1,Npl
    write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'layer',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),trim(epart),i,nkpt,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=17+(mpitag-1)*Npl+i, file=varm,status='replace')
  end do
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'s_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),trim(epart),nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=96+mpitag, file=varm,status='replace')
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'p_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),trim(epart),nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=97+mpitag, file=varm,status='replace')
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'d_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),trim(epart),nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=98+mpitag, file=varm,status='replace')
  write(varm,"('./results/',a1,'SOC/',a,'/FS/',a,'total_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),trim(epart),nkpt,eta,Utype,trim(fieldpart),trim(socpart)
  open (unit=99+mpitag, file=varm,status='replace')

  fs_layer = 0.d0
  fs_orb   = 0.d0
  fs_total = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,i,mu,nu,sigma,temp1,temp2) &
!$omp& shared(lverbose,llineargfsoc,llinearsoc,e,kbz,nkpt,Ef,eta,Npl,nthreads,pi,pauli_orb,pauli_gf,fs_layer,fs_orb,fs_total,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if(mythread==0) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[fermi_surface] Number of threads: ',i0)") nthreads
!$  end if

!$omp do
  fermi_surface_kpoints: do iz=1,nkpt
!$  if((mythread==0)) then
    if(lverbose) call progress_bar(outputunit_loop,"kpoints",iz,nkpt)
!$   end if

    kp = kbz(:,iz)

    ! Green function on energy Ef + ieta, and wave vector kp
    if((llinearsoc).or.(llineargfsoc)) then
      call greenlineargfsoc(e,eta,kp,gf)
    else
      call green(e,eta,kp,gf)
    end if

    do sigma=1,4 ; do i=1,Npl
      if(sigma==1) then
        pauli_gf = gf(i,i,:,:)
      else
        temp1 = pauli_orb(sigma-1,:,:)
        temp2 = gf(i,i,:,:)
        call zgemm('n','n',18,18,18,zum,temp1,18,temp2,18,zero,pauli_gf,18)
      end if
      do mu=1,9
        nu = mu + 9
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
  writing_fermi_surface: do iz=1,nkpt
    write_plane_loop_fs: do i=1,Npl
      write(unit=17+(mpitag-1)*Npl+i,fmt="(6(es16.9,2x))") kbz(1,iz),kbz(2,iz), kbz(3,iz), (fs_layer(i,iz,sigma),sigma=1,4)
    end do write_plane_loop_fs

    write(unit=96+mpitag,fmt="(6(es16.9,2x))") kbz(1,iz),kbz(2,iz),kbz(3,iz),(fs_orb(1,iz,sigma),sigma=1,4)
    write(unit=97+mpitag,fmt="(6(es16.9,2x))") kbz(1,iz),kbz(2,iz),kbz(3,iz),(fs_orb(2,iz,sigma),sigma=1,4)
    write(unit=98+mpitag,fmt="(6(es16.9,2x))") kbz(1,iz),kbz(2,iz),kbz(3,iz),(fs_orb(3,iz,sigma),sigma=1,4)
    write(unit=99+mpitag,fmt="(6(es16.9,2x))") kbz(1,iz),kbz(2,iz),kbz(3,iz),(fs_total(iz,sigma),sigma=1,4)
  end do writing_fermi_surface
  ! Closing files
  do i=1,Npl
    close (17+(mpitag-1)*Npl+i)
  end do
  close(96+mpitag)
  close(97+mpitag)
  close(98+mpitag)
  close(99+mpitag)

  return
end subroutine fermi_surface
