module mod_Coupling
  use mod_f90_kind, only: double
  implicit none
  character(len=3), private :: folder = "Jij"
  character(len=3), dimension(3), private :: filename = ["Jij", "J  ", "Dz "]

  ! Exchange interaction
  real(double), dimension(:,:), allocatable :: trJij
  real(double), dimension(:,:,:,:), allocatable :: Jij,Jijs,Jija

contains

  subroutine allocateCoupling()
    use mod_parameters, only: nmaglayers
    implicit none

    if(allocated(trJij)) deallocate(trJij)
    if(allocated(Jij)) deallocate(Jij)
    if(allocated(Jijs)) deallocate(Jijs)
    if(allocated(Jija)) deallocate(Jija)
    allocate(trJij(nmaglayers,nmaglayers))
    allocate(Jij(nmaglayers,nmaglayers,3,3))
    allocate(Jijs(nmaglayers,nmaglayers,3,3))
    allocate(Jija(nmaglayers,nmaglayers,3,3))

  end subroutine allocateCoupling

  subroutine deallocateCoupling()
    implicit none
    if(allocated(trJij)) deallocate(trJij)
    if(allocated(Jij)) deallocate(Jij)
    if(allocated(Jijs)) deallocate(Jijs)
    if(allocated(Jija)) deallocate(Jija)
  end subroutine deallocateCoupling

  subroutine openCouplingFiles()
    use mod_parameters, only: nmaglayers, mmlayermag, output
    implicit none
    character(len=400) :: varm
    integer :: i, j, iw

    do j=1,nmaglayers
      do i=1,nmaglayers
        iw = 2000 + (j-1) * nmaglayers * 2 + (i-1) * 2
        if(i==j) then
         iw = iw + 1
         write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(1)),mmlayermag(i)-1,trim(output%info),trim(output%BField),trim(output%SOC)
         open (unit=iw, file=varm,status='replace')
         write(unit=iw, fmt="('#   energy      ,  Jii_xx           ,   Jii_yy  ')")
         ! Anisotropy energy is given by K^a = 2*J_ii^aa
         ! omega_res ~ gamma*m_i*J_ii (*2?) ,
         ! where J_ii is the one calculated here
        else
         iw = iw + 1
         write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,'_',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(2)),mmlayermag(i)-1,mmlayermag(j)-1,trim(output%info),trim(output%BField),trim(output%SOC)
         open (unit=iw, file=varm,status='replace')
         write(unit=iw, fmt="('#   energy      ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
         iw = iw + 1
         write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_',i0,'_',i0,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(3)),mmlayermag(i)-1,mmlayermag(j)-1,trim(output%info),trim(output%BField),trim(output%SOC)
         open (unit=iw, file=varm,status='replace')
         write(unit=iw, fmt="('#   energy      , Dz = (Jxy - Jyx)/2       ')")
        end if
      end do
    end do
  end subroutine openCouplingFiles

  subroutine closeCouplingFiles()
    use mod_parameters, only: nmaglayers
    implicit none
    integer :: i, j, iw
    do j=1,nmaglayers
      do i=1,nmaglayers
        iw = 2000 + (j-1) * nmaglayers * 2 + (i-1)*2 + 1
        close (iw)
        if(i/=j) close(iw+1)
      end do
    end do
  end subroutine closeCouplingFiles

  subroutine writeCoupling(e)
    use mod_f90_kind, only: double
    use mod_parameters, only: nmaglayers
    implicit none
    real(double), intent(in) :: e
    integer :: i, j, iw

    do j=1,nmaglayers
      do i=1,nmaglayers
        iw = 2000 + (j-1) * nmaglayers * 2 + (i-1) * 2
        if(i==j) then
         iw = iw + 1
         write(unit=iw,fmt="(3(es16.9,2x))") e,Jij(i,j,1,1),Jij(i,j,2,2)
        else
         iw = iw + 1
         write(unit=iw,fmt="(4(es16.9,2x))") e,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
         iw = iw + 1
         write(unit=iw,fmt="(2(es16.9,2x))") e,Jija(i,j,1,2)
        end if
      end do
    end do
  end subroutine writeCoupling

end module mod_Coupling
