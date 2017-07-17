! This subroutine calculates LDOS and coupling as a function of energy
subroutine ldos_and_coupling()
  use mod_f90_kind, only: double
  use mod_parameters, only: outputunit_loop,Npl_folder,mmlayermag,nmaglayers, emin, deltae, npt1, Utype, eta, fieldpart
  use mod_SOC, only: SOCc, socpart
  use mod_mpi_pars, only: mpitag
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  implicit none
  character(len=400) :: varm
  integer            :: i,iw,j,mu, count
  real(double)       :: e
  real(double),dimension(:,:),allocatable :: ldosu,ldosd
  ! Exchange interaction
  real(double),dimension(:,:),allocatable     :: trJij
  real(double),dimension(:,:,:,:),allocatable :: Jij,Jijs,Jija

  allocate(ldosu(s%nAtoms,nOrb),ldosd(s%nAtoms,nOrb))
  allocate(trJij(nmaglayers,nmaglayers),Jij(nmaglayers,nmaglayers,3,3),Jijs(nmaglayers,nmaglayers,3,3),Jija(nmaglayers,nmaglayers,3,3))

  write(outputunit_loop,"('CALCULATING LDOS AND EXCHANGE INTERACTIONS AS A FUNCTION OF ENERGY')")

  ! Opening files
  ! LDOS
  do i=1,s%nAtoms
    iw = 1000+(mpitag-1)*s%nAtoms*2 + (i-1)*2 + 1
    write(varm,"('./results/',a1,'SOC/',a,'/LDOS/ldosu_layer',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),i,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=iw, file=varm,status='replace')
    iw = iw+1
    write(varm,"('./results/',a1,'SOC/',a,'/LDOS/ldosd_layer',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),i,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
    open (unit=iw, file=varm,status='replace')
  end do
  ! Exchange interactions
  do j=1,nmaglayers ; do i=1,nmaglayers
    iw = 2000+(mpitag-1)*nmaglayers*nmaglayers*2+(j-1)*nmaglayers*2+(i-1)*2
    if(i==j) then
      iw = iw + 1
      write(varm,"('./results/',a1,'SOC/',a,'/Jij/Jii_',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),mmlayermag(i)-1,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
      open (unit=iw, file=varm,status='replace')
      write(unit=iw, fmt="('#   energy      ,  Jii_xx           ,   Jii_yy  ')")
      iw = iw + 1
      ! Anisotropy energy is given by K^a = 2*J_ii^aa
      ! omega_res ~ gamma*m_i*J_ii (*2?) ,
      ! where J_ii is the one calculated here
    else
      iw = iw + 1
      write(varm,"('./results/',a1,'SOC/',a,'/Jij/J_',i0,'_',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),mmlayermag(i)-1,mmlayermag(j)-1,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
      open (unit=iw, file=varm,status='replace')
      write(unit=iw, fmt="('#   energy      ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
      iw = iw + 1
      write(varm,"('./results/',a1,'SOC/',a,'/Jij/Dz_',i0,'_',i0,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,trim(Npl_folder),mmlayermag(i)-1,mmlayermag(j)-1,s%nkpt,eta,Utype,trim(fieldpart),trim(socpart)
      open (unit=iw, file=varm,status='replace')
      write(unit=iw, fmt="('#   energy      , Dz = (Jxy - Jyx)/2       ')")
    end if
  end do ; end do

  ldos_energy_loop: do count=1,npt1
    e = emin + (count-1)*deltae
    write(outputunit_loop,"('[ldos_and_coupling] ',i0,' of ',i0,' points',', e = ',es10.3)") count,npt1,e

    call ldos_jij_energy(e,ldosu,ldosd,Jij)

    do i=1,nmaglayers ; do j=1,nmaglayers
      trJij(i,j)    = 0.5d0*(Jij(i,j,1,1)+Jij(i,j,2,2))
      Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
      Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
      do mu=1,3
        Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
      end do
    end do ; end do

    ! Writing into files
    ! LDOS
    ldos_writing_plane_loop: do i=1,s%nAtoms
        iw = 1000+(mpitag-1)*s%nAtoms*2 + (i-1)*2 + 1
        write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
        iw = iw+1
        write(unit=iw,fmt="(5(es16.9,2x))") e,sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
    end do ldos_writing_plane_loop

    ! Exchange interactions
    jij_writing_loop: do j=1,nmaglayers ; do i=1,nmaglayers
      iw = 2000+(mpitag-1)*nmaglayers*nmaglayers*2+(j-1)*nmaglayers*2+(i-1)*2
      if(i==j) then
        iw = iw + 1
        write(unit=iw,fmt="(3(es16.9,2x))") e,Jij(i,j,1,1),Jij(i,j,2,2)
        iw = iw + 1
      else
        iw = iw + 1
        write(unit=iw,fmt="(4(es16.9,2x))") e,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
        iw = iw + 1
        write(unit=iw,fmt="(2(es16.9,2x))") e,Jija(i,j,1,2)
      end if
    end do ; end do jij_writing_loop

  end do ldos_energy_loop

  ! Closing files
  do i=1,s%nAtoms
    iw = 1000+(mpitag-1)*s%nAtoms*2 + (i-1)*2 + 1
    close (iw)
    iw = iw+1
    close (iw)
  end do
  do j=1,nmaglayers ; do i=1,nmaglayers
    iw = 2000+(mpitag-1)*nmaglayers*nmaglayers*2+(j-1)*nmaglayers*2+(i-1)*2
    if(i==j) then
      iw = iw + 1
      close (iw)
      iw = iw + 1
    else
      iw = iw + 1
      close (iw)
      iw = iw + 1
      close (iw)
    end if
  end do ; end do

  deallocate(ldosu,ldosd)
  deallocate(trJij,Jij,Jijs,Jija)

  return
end subroutine ldos_and_coupling
