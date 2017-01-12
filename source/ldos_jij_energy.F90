!   Calculates spin-resolved LDOS and energy-dependence of exchange interactions
subroutine ldos_jij_energy(e,ldosu,ldosd,Jijint)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_generate_kpoints
  use mod_progress
  use mod_magnet, only: hdel,mz
  use mod_mpi_pars
!$  use omp_lib
  implicit none
!$  integer           :: nthreads,mythread
  integer             :: i,j,mu,nu,iz,alpha
  real(double)        :: kp(3),Jijkan(nmaglayers,3,3),Jijk(nmaglayers,nmaglayers,3,3)
  real(double),intent(in)     :: e
  real(double),intent(out)    :: ldosu(Npl,9),ldosd(Npl,9)
  real(double),intent(out)    :: Jijint(nmaglayers,nmaglayers,3,3)
  complex(double)     :: paulimatan(3,3,18,18)
  complex(double),dimension(Npl,Npl,18,18)    :: gf
  complex(double),dimension(Npl,9)            :: gfdiagu,gfdiagd
  complex(double),dimension(18,18)            :: gij,gji,temp1,temp2,paulia,paulib

! (x,y,z)-tensor formed by Pauli matrices to calculate anisotropy term (when i=j)
  paulimatan = zero
  paulimatan(1,1,:,:) = -pauli_dorb(3,:,:)
  paulimatan(2,2,:,:) = -pauli_dorb(3,:,:)
  paulimatan(1,3,:,:) = -pauli_dorb(1,:,:)
  paulimatan(3,1,:,:) = -pauli_dorb(1,:,:)
  paulimatan(2,3,:,:) = -pauli_dorb(2,:,:)
  paulimatan(3,2,:,:) = -pauli_dorb(2,:,:)

  ldosu = 0.d0
  ldosd = 0.d0
  Jijint = 0.d0

!$omp parallel default(none) &
!$omp& private(mythread,iz,kp,gf,gij,gji,paulia,paulib,i,j,mu,nu,alpha,gfdiagu,gfdiagd,Jijk,Jijkan,temp1,temp2) &
!$omp& shared(lverbose,kbz,nkpoints,wkbz,e,eta,Npl,hdel,mz,nmaglayers,mmlayermag,pauli_dorb,paulimatan,ldosu,ldosd,Jijint,nthreads,outputunit_loop)
!$  mythread = omp_get_thread_num()
!$  if(mythread.eq.0) then
!$    nthreads = omp_get_num_threads()
!$    write(outputunit_loop,"('[ldos_jij_energy] Number of threads: ',i0)") nthreads
!$  end if

!$omp do reduction(+:ldosu,ldosd,Jijint)
  kpoints: do iz=1,nkpoints
!$  if((mythread.eq.0)) then
      if(lverbose) call progress_bar(outputunit_loop,"kpoints",iz,nkpoints)
!$  end if
    kp = kbz(iz,:)

    ! Green function on energy E + ieta, and wave vector kp
    call green(e,eta,kp,gf)

    ! Exchange interaction tensor
    Jijk   = 0.d0
    Jijkan = 0.d0
    do nu = 1,3 ; do mu = 1,3 ; do j = 1,nmaglayers ; do i = 1,nmaglayers
      paulia = pauli_dorb(mu,:,:)
      gij = gf(mmlayermag(i)-1,mmlayermag(j)-1,:,:)
      paulib = pauli_dorb(nu,:,:)
      gji = gf(mmlayermag(j)-1,mmlayermag(i)-1,:,:)
      call zgemm('n','n',18,18,18,zum,paulia,18,gij,18,zero,temp1,18)
      call zgemm('n','n',18,18,18,zum,temp1,18,paulib,18,zero,temp2,18)
      call zgemm('n','n',18,18,18,zum,temp2,18,gji,18,zero,temp1,18)
      do alpha = 1,18
        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
      end do

      Jijk(i,j,mu,nu) = -hdel(mmlayermag(i)-1)*hdel(mmlayermag(j)-1)*Jijk(i,j,mu,nu)*wkbz(iz)/(mz(mmlayermag(i)-1)*mz(mmlayermag(j)-1))

      ! Anisotropy (on-site) term
      if(i.eq.j) then
        gij = gf(mmlayermag(i)-1,mmlayermag(i)-1,:,:)
        paulia = paulimatan(mu,nu,:,:)
        call zgemm('n','n',18,18,18,zum,gij,18,paulia,18,zero,temp1,18)

        do alpha = 1,18
          Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
        end do

        Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + hdel(mmlayermag(i)-1)*Jijkan(i,mu,nu)*wkbz(iz)/(mz(mmlayermag(i)-1)**2)
      end if
    end do ; end do ; end do ; end do

    ! Density of states
    do mu=1,9; do i=1,Npl
      nu=mu+9
      gfdiagu(i,mu) = - aimag(gf(i,i,mu,mu))*wkbz(iz)
      gfdiagd(i,mu) = - aimag(gf(i,i,nu,nu))*wkbz(iz)
    end do ; end do

    ldosu = ldosu + gfdiagu
    ldosd = ldosd + gfdiagd
    Jijint = Jijint + Jijk

  end do kpoints
!$omp end do
!$omp end parallel

  ldosu  = ldosu/pi
  ldosd  = ldosd/pi
  Jijint = Jijint/pi

  return
end subroutine ldos_jij_energy

! !   Calculates spin-resolved LDOS (including empty spheres)
! subroutine ldos_jij_energy_es(e)
!   use mod_f90_kind
!   use mod_constants
!   use mod_parameters
!   use mod_generate_kpoints
!   use mod_progress
!   use mod_mpi_pars
! !$  use omp_lib
!   implicit none
! !$  integer       :: nthreads,mythread
!   character(len=400) :: varm
!   integer         :: i,j,mu,nu,iz
!   real(double)    :: ldosu(Npl+2,9),ldosd(Npl+2,9)
!   real(double)    :: e,kp(3)
!   complex(double),dimension(Npl,9)                :: gfdiagu,gfdiagd
!   complex(double),dimension(Npl+2,Npl+2,18,18)    :: gf

!   ldosu = 0.d0
!   ldosd = 0.d0

! !$omp parallel default(none) &
! !$omp& private(mythread,iz,kp,gf,i,j,mu,nu,nthreads) &
! !$omp& shared(prog,spiner,elapsed_time,start_program,progbar,kbz,nkpoints,wkbz,e,eta,Npl,gfdiagu,gfdiagd,ldosu,ldosd,myrank_row_hw)
! !$  mythread = omp_get_thread_num()
! !$  if((mythread.eq.0).and.(myrank_row_hw.eq.0)) then
! !$    nthreads = omp_get_num_threads()
! !$    write(outputunit_loop,"('[ldos_jij_energy_es] Number of threads: ',i0)") nthreads
! !$  end if

! !$omp do reduction(+:ldosu,ldosd)
!   kpoints: do iz=1,nkpoints
! !$  if((mythread.eq.0)) then
!       if(lverbose) call progress_bar(outputunit_loop,"kpoints",iz,nkpoints)
! !$  end if

!     kp = kbz(iz,:)

!     ! Green function on energy E + ieta, and wave vector kp
!     call green_es(e,eta,kp,gf)

!     do mu=1,9; do i=1,Npl+2
!       nu=mu+9
!       gfdiagu(i,mu) = - aimag(gf(i,i,mu,mu))*wkbz(iz)
!       gfdiagd(i,mu) = - aimag(gf(i,i,nu,nu))*wkbz(iz)
!     end do ; end do

!     ldosu = ldosu + gfdiagu
!     ldosd = ldosd + gfdiagd

!   end do kpoints
! !$omp end do
! !$omp end parallel

!   ldosu  = ldosu/pi
!   ldosd  = ldosd/pi

!   ! Transform energy to eV if runoption is on
!   e = e

!   !Writing on files
!   do i=1,Npl+2
!     write(varm,"('./results/',l1,'SOC/',a,'/LDOS/ldosu_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,trim(Npl_folder),i,magaxis,socscale,ncp,eta,Utype,hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
!     open (unit=117+i, file=varm,status='unknown')
!     write(varm,"('./results/',l1,'SOC/',a,'/LDOS/ldosd_layer',I0,'_magaxis=',A,'_socscale=',f5.2,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,'_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2,'.dat')") SOC,trim(Npl_folder),i,magaxis,socscale,ncp,eta,Utype,hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
!     open (unit=517+i, file=varm,status='unknown')
!   end do
!   ldos_writing_plane_loop: do i=1,Npl+2
!       write(unit=117+i,fmt="(5(es16.9,2x))") e,sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
!       write(unit=517+i,fmt="(5(es16.9,2x))") e,sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
!   end do ldos_writing_plane_loop

!   do i=1,Npl+2
!     close (117+i)
!     close (517+i)
!   end do

!   return
! end subroutine ldos_jij_energy_es