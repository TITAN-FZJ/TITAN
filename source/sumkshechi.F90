! ----------- Sum over wave vectors to calculate spin disturbance -----------
subroutine sumkshechi(e,ep,Fint,iflag)
  use mod_f90_kind, only: double
  use mod_parameters, only: eta, ef, dim, sigmaijmunu2i, sigmaimunu2i
  use mod_constants, only: zero, zum, zi, tpi
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use mod_SOC, only: llineargfsoc
  use mod_mpi_pars, only: abortProgram
!$  use omp_lib
  implicit none
  integer         :: AllocateStatus
  integer         :: i,j,mu,nu,gamma,xi,iz
  integer, intent(in) :: iflag
  real(double),intent(in)     :: e,ep
  real(double)                :: kp(3)
  complex(double),dimension(:,:,:,:),allocatable    :: gf
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd
  complex(double),dimension(dim,dim),intent(out)    :: Fint
  complex(double),dimension(:,:),allocatable        :: df1, Fint_loc

  Fint      = zero

  if(iflag==0) then
    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,iz,kp,gf,gfuu,gfud,gfdu,gfdd,i,j,mu,nu,gamma,xi,df1, Fint_loc) &
    !$omp& shared(llineargfsoc,s,e,ep,iflag,Fint,Ef,eta,dim,sigmaimunu2i,sigmaijmunu2i)
    allocate(df1(dim,dim), Fint_loc(dim,dim), &
             gf  (s%nAtoms,s%nAtoms,2*nOrb,2*nOrb), &
             gfuu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
             gfud(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
             gfdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
             gfdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumkshechi] Not enough memory for: df1,gf,gfuu,gfud,gfdu,gfdd")
    Fint_loc = zero

    !$omp do schedule(static)
    do iz = 1, s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E_F+E+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef+e,ep,kp,gf)
      else
        call green(Ef+e,ep,kp,gf)
      end if
      gfuu(:,:,:,:,1) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,1) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E_F+iy)
      if(llineargfsoc) then
        call greenlineargfsoc(Ef,ep,kp,gf)
      else
        call green(Ef,ep,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,2) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

!dir$ simd
      do xi=1,nOrb
        do gamma=1,nOrb
          do j=1,s%nAtoms
            do nu=1,nOrb
              do mu=1,nOrb
                do i=1,s%nAtoms
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)

                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)

                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)

                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
                  df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)
                end do
              end do
            end do
          end do
        end do
      end do

      ! Locally add up df1
#ifdef _JUQUEEN
      call zgeadd(Fint_loc,dim,'N',df1,dim,'N',Fint_loc,dim,dim,dim)
#else
      call ZAXPY(dim*dim,zum,df1,1,Fint_loc,1)
#endif
    end do
    !$omp end do

    ! Add up df1 on MPI thread
    !$omp critical
      Fint = Fint + Fint_loc
    !$omp end critical
    deallocate(Fint_loc)
    deallocate(df1)
    deallocate(gf,gfuu,gfud,gfdu,gfdd)
    !$omp end parallel


  else
    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,iz,kp,gf,gfuu,gfud,gfdu,gfdd,i,j,mu,nu,gamma,xi,df1, Fint_loc) &
    !$omp& shared(llineargfsoc,s,e,ep,iflag,Fint,Ef,eta,dim,sigmaimunu2i,sigmaijmunu2i)
    allocate(df1(dim,dim), Fint_loc(dim,dim), &
             gf  (s%nAtoms,s%nAtoms,2*nOrb,2*nOrb), &
             gfuu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
             gfud(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
             gfdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2), &
             gfdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumkshechi] Not enough memory for: df1,gf,gfuu,gfud,gfdu,gfdd")
    Fint_loc = zero

    !$omp do schedule(static)
    do iz = 1, s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E'+E+i.eta)
      if(llineargfsoc) then
        call greenlineargfsoc(ep+e,eta,kp,gf)
      else
        call green(ep+e,eta,kp,gf)
      end if
      gfuu(:,:,:,:,1) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,1) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E'+i.eta)
      if(llineargfsoc) then
        call greenlineargfsoc(ep,eta,kp,gf)
      else
        call green(ep,eta,kp,gf)
      end if
      gfuu(:,:,:,:,2) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,2) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

  !dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1) - conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1) - conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1) - conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1) - conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do

      ! Locally add up df1
#ifdef _JUQUEEN
      call zgeadd(Fint_loc,dim,'N',df1,dim,'N',Fint_loc,dim,dim,dim)
#else
      call ZAXPY(dim*dim,zum,df1,1,Fint_loc,1)
#endif
    end do
    !$omp end do

    ! Add up df1 on MPI thread
    !$omp critical
      Fint = Fint + Fint_loc
    !$omp end critical

    deallocate(df1)
    deallocate(gf,gfuu,gfud,gfdu,gfdd)
    !$omp end parallel
  end if

  Fint = Fint/tpi

  return
end subroutine sumkshechi


! ----------- Sum over wave vectors to calculate spin susceptibility -----------
! -------------- to be used in the calculation of linear SOC chi ---------------
subroutine sumkshechilinearsoc(e,ep,Fint,Fintlsoc,iflag)
  use mod_f90_kind, only: double
  use mod_parameters, only: eta, ef, dim, sigmaijmunu2i, sigmaimunu2i
  use mod_constants, only: zero, zum, zi, tpi
  use mod_mpi_pars, only: abortProgram
  use mod_system, only: s => sys
  use TightBinding, only: nOrb
  use mod_SOC, only: llineargfsoc
!$  use omp_lib
  implicit none
  integer         :: AllocateStatus
  integer         :: i,j,mu,nu,gamma,xi,iz
  integer, intent(in) :: iflag
  real(double),intent(in)     :: e,ep
  real(double)                :: kp(3)
  complex(double),dimension(:,:,:,:),allocatable    :: gf,gvg
  complex(double),dimension(:,:,:,:,:),allocatable  :: gfuu,gfud,gfdu,gfdd,gvguu,gvgud,gvgdu,gvgdd
  complex(double),dimension(dim,dim),intent(out)    :: Fint,Fintlsoc
  complex(double),dimension(:,:),allocatable        :: df1,df1lsoc

  Fint      = zero
  Fintlsoc  = zero


  if(iflag==0)then
    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,iz,kp,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,i,j,mu,nu,gamma,xi,df1,df1lsoc) &
    !$omp& shared(llineargfsoc,s,e,ep,iflag,Fint,Fintlsoc,Ef,eta,dim,sigmaimunu2i,sigmaijmunu2i)
    allocate(df1(dim,dim),gf(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb),gfuu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gfud(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gfdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gfdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumkshechilinearsoc] Not enough memory for: df1,gf,gfuu,gfud,gfdu,gfdd")

    allocate( df1lsoc(dim,dim),gvg(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb),gvguu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gvgud(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gvgdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gvgdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumkshechilinearsoc] Not enough memory for: df1lsoc,gvg,gvguu,gvgud,gvgdu,gvgdd")

    !$omp do schedule(static), reduction(+:Fint), reduction(+:Fintlsoc)
    do iz=1,s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E_F+E+iy)
      call greenlinearsoc(Ef+e,ep,kp,gf,gvg)
      gfuu(:,:,:,:,1) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,1) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,1) = gvg(:,:, 1: nOrb, 1: nOrb)
      gvgud(:,:,:,:,1) = gvg(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gvgdd(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(Ef,ep,kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,2) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,2) = gvg(:,:, 1: nOrb, 1: nOrb)
      gvgud(:,:,:,:,2) = gvg(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gvgdd(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

!dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1)))*s%wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gfud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gfuu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1)))*s%wkbz(iz)


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*s%wkbz(iz)

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfuu(j,i,xi,mu,2) + conjg(gvguu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvguu(j,i,xi,mu,2) + conjg(gfuu(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfdu(j,i,xi,mu,2) + conjg(gvgud(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdu(j,i,xi,mu,2) + conjg(gfud(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*s%wkbz(iz)

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgdd(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfdd(j,i,gamma,nu,1))  +  gfdd(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgdd(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvgdu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfud(j,i,gamma,nu,1))  +  gfdu(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgud(j,i,gamma,nu,1)))*s%wkbz(iz)

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfud(j,i,xi,mu,2) + conjg(gvgdu(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgud(j,i,xi,mu,2) + conjg(gfdu(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = (gvgud(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfdu(j,i,gamma,nu,1))  +  gfud(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvgdu(j,i,gamma,nu,1)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = (gvguu(i,j,nu,gamma,1)*gfdd(j,i,xi,mu,2) + conjg(gvgdd(i,j,mu,xi,2)*gfuu(j,i,gamma,nu,1))  +  gfuu(i,j,nu,gamma,1)*gvgdd(j,i,xi,mu,2) + conjg(gfdd(i,j,mu,xi,2)*gvguu(j,i,gamma,nu,1)))*s%wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do

#ifdef _JUQUEEN
      call zgeadd(Fint,dim,'N',df1,dim,'N',Fint,dim,dim,dim)
      call zgeadd(Fintlsoc,dim,'N',df1lsoc,dim,'N',Fintlsoc,dim,dim,dim)
#else
      Fint = Fint + df1
      Fintlsoc = Fintlsoc + df1lsoc
#endif
    end do
    !$omp end do
    deallocate(df1)
    deallocate(gvg,gvguu,gvgud,gvgdu,gvgdd)
    deallocate(gf,gfuu,gfud,gfdu,gfdd)
    !$omp end parallel

  else
    !$omp parallel default(none) &
    !$omp& private(AllocateStatus,iz,kp,gf,gfuu,gfud,gfdu,gfdd,gvg,gvguu,gvgud,gvgdu,gvgdd,i,j,mu,nu,gamma,xi,df1,df1lsoc) &
    !$omp& shared(llineargfsoc,s,e,ep,iflag,Fint,Fintlsoc,Ef,eta,dim,sigmaimunu2i,sigmaijmunu2i)
    allocate(df1(dim,dim),gf(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb),gfuu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gfud(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gfdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gfdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumkshechilinearsoc] Not enough memory for: df1,gf,gfuu,gfud,gfdu,gfdd")

    allocate( df1lsoc(dim,dim),gvg(s%nAtoms,s%nAtoms,2*nOrb,2*nOrb),gvguu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gvgud(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gvgdu(s%nAtoms,s%nAtoms,nOrb,nOrb,2),gvgdd(s%nAtoms,s%nAtoms,nOrb,nOrb,2), STAT = AllocateStatus  )
    if (AllocateStatus/=0) call abortProgram("[sumkshechilinearsoc] Not enough memory for: df1lsoc,gvg,gvguu,gvgud,gvgdu,gvgdd")

    !$omp do schedule(static), reduction(+:Fint), reduction(+:Fintlsoc)
    do iz=1,s%nkpt
      kp = s%kbz(:,iz)

      ! Green function at (k+q,E_F+E+iy)
      call greenlinearsoc(ep+e,eta,kp,gf,gvg)
      gfuu(:,:,:,:,1) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,1) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,1) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,1) = gvg(:,:, 1: nOrb, 1: nOrb)
      gvgud(:,:,:,:,1) = gvg(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gvgdd(:,:,:,:,1) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

      ! Green function at (k,E_F+iy)
      call greenlinearsoc(ep,eta,kp,gf,gvg)
      gfuu(:,:,:,:,2) = gf(:,:, 1: nOrb, 1: nOrb)
      gfud(:,:,:,:,2) = gf(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gfdu(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gfdd(:,:,:,:,2) = gf(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)
      gvguu(:,:,:,:,2) = gvg(:,:, 1: nOrb, 1: nOrb)
      gvgud(:,:,:,:,2) = gvg(:,:, 1: nOrb,nOrb+1:2*nOrb)
      gvgdu(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb, 1: nOrb)
      gvgdd(:,:,:,:,2) = gvg(:,:,nOrb+1:2*nOrb,nOrb+1:2*nOrb)

!dir$ simd
      do xi=1,nOrb ; do gamma=1,nOrb ; do j=1,s%nAtoms ; do nu=1,nOrb ; do mu=1,nOrb ; do i=1,s%nAtoms
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)

        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))*s%wkbz(iz)

        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)

        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*(gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)
        df1(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*(gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))*s%wkbz(iz)


        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(1,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*s%wkbz(iz)

        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfuu(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvguu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(2,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfud(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgud(i,j,mu,xi,2)))*s%wkbz(iz)

        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgdd(i,j,nu,gamma,1)-conjg(gvgdd(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfdd(i,j,nu,gamma,1)-conjg(gfdd(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(3,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvgdu(i,j,nu,gamma,1)-conjg(gvgud(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfdu(i,j,nu,gamma,1)-conjg(gfud(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*s%wkbz(iz)

        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(1,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(2,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdu(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdu(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(3,j,gamma,xi)) = -zi*((gvgud(i,j,nu,gamma,1)-conjg(gvgdu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfud(i,j,nu,gamma,1)-conjg(gfdu(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*s%wkbz(iz)
        df1lsoc(sigmaimunu2i(4,i,mu,nu),sigmaimunu2i(4,j,gamma,xi)) = -zi*((gvguu(i,j,nu,gamma,1)-conjg(gvguu(j,i,gamma,nu,1)))*conjg(gfdd(i,j,mu,xi,2))  +  (gfuu(i,j,nu,gamma,1)-conjg(gfuu(j,i,gamma,nu,1)))*conjg(gvgdd(i,j,mu,xi,2)))*s%wkbz(iz)
      end do ; end do ; end do ; end do ; end do ; end do

#ifdef _JUQUEEN
      call zgeadd(Fint,dim,'N',df1,dim,'N',Fint,dim,dim,dim)
      call zgeadd(Fintlsoc,dim,'N',df1lsoc,dim,'N',Fintlsoc,dim,dim,dim)
#else
      Fint = Fint + df1
      Fintlsoc = Fintlsoc + df1lsoc
#endif
    end do
    !$omp end do
    deallocate(df1)
    deallocate(gvg,gvguu,gvgud,gvgdu,gvgdd)
    deallocate(gf,gfuu,gfud,gfdu,gfdd)
    !$omp end parallel
  end if

  Fint     = Fint/tpi
  Fintlsoc = Fintlsoc/tpi

  return
end subroutine sumkshechilinearsoc
