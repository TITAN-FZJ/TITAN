! Calculate hamiltonian of a slab containing
! Npl layers + 1 layer of Empty spheres on each side
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!  1     2     3     4      Npl-1   Npl  Npl+1  Npl+2
!         <-S-> <S-1>           <S-1> <-S->
subroutine hamiltk(kp,hk)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_lattice, only: plnn
  use mod_tight_binding, only: lambda, ls
  use mod_magnet
  implicit none
  integer     :: i,i0,i1,j0,j1
  real(double), intent(in)  :: kp(3)
  complex(double) :: hee(Npl+2,18,18)
  complex(double),dimension(18,18)    :: h00,h01,h10,h20,h02
  complex(double),dimension((Npl+2)*18,(Npl+2)*18),intent(out)  :: hk

  hk = zero

  call U_matrix(hee)

! Mouting slab hamiltonian
  do i=1,Npl+2
    i0 = (i-1)*18+1
    i1 = i0+17

    select case (plnn)
    case( 1 )
      call helphbccsoc(kp,h00,h01,h10,i)

      hk(i0:i1,i0:i1) = h00 + lb(i,:,:) + sb(i,:,:) + hee(i,:,:) + (socscale*lambda(i)*ls)
      if (i/=(Npl+2)) then
        j0 = i0+18
        j1 = i1+18
        hk(i0:i1,j0:j1) = h01(:,:)
        hk(j0:j1,i0:i1) = h10(:,:)
      end if
    case( 2 )
      call helphfccsoc(kp,h00,h01,h10,h02,h20,i)

      hk(i0:i1,i0:i1) = h00 + lb(i,:,:) + sb(i,:,:) + hee(i,:,:) + (socscale*lambda(i)*ls)
      if (i/=(Npl+2)) then
        j0 = i0+18
        j1 = i1+18
        hk(i0:i1,j0:j1) = h01(:,:)
        hk(j0:j1,i0:i1) = h10(:,:)
        if (i/=(Npl+1)) then
          j0 = j0+18
          j1 = j1+18
          hk(i0:i1,j0:j1) = h02(:,:)
          hk(j0:j1,i0:i1) = h20(:,:)
        end if
      end if
    end select
  end do

  return
end subroutine hamiltk

! Calculate hamiltonian of a slab containing
! Npl layers + 1 layer of Empty spheres on each side
!  ES    S    S-1   S-2       S-2   S-1    S     ES
!  o-----|-----|-----|---...---|-----|-----|-----o
!  1     2     3     4      Npl-1   Npl  Npl+1  Npl+2
!         <-S-> <S-1>           <S-1> <-S->
subroutine hamiltklinearsoc(kp,hk,vsoc)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_lattice, only: plnn
  use mod_tight_binding, only: lambda, ls
  use mod_magnet
  implicit none
  integer     :: i,i0,i1,j0,j1
  real(double), intent(in)  :: kp(3)
  complex(double) :: hee(Npl+2,18,18)
  complex(double),dimension(18,18)    :: h00,h01,h10,h20,h02
  complex(double),dimension((Npl+2)*18,(Npl+2)*18),intent(out)  :: hk,vsoc

  hk = zero
  vsoc = zero

  call U_matrix(hee)

! Mouting slab hamiltonian
  do i=1,Npl+2
    i0 = (i-1)*18+1
    i1 = i0+17

    select case (plnn)
    case( 1 )
      call helphbccsoc(kp,h00,h01,h10,i)

      hk(i0:i1,i0:i1)   = h00 + lb(i,:,:) + sb(i,:,:) + hee(i,:,:)
      vsoc(i0:i1,i0:i1) = socscale*lambda(i)*ls
      if (i/=(Npl+2)) then
        j0 = i0+18
        j1 = i1+18
        hk(i0:i1,j0:j1) = h01(:,:)
        hk(j0:j1,i0:i1) = h10(:,:)
      end if
    case( 2 )
      call helphfccsoc(kp,h00,h01,h10,h02,h20,i)

      hk(i0:i1,i0:i1)   = h00 + lb(i,:,:) + sb(i,:,:) + hee(i,:,:)
      vsoc(i0:i1,i0:i1) = socscale*lambda(i)*ls
      if (i/=(Npl+2)) then
        j0 = i0+18
        j1 = i1+18
        hk(i0:i1,j0:j1) = h01(:,:)
        hk(j0:j1,i0:i1) = h10(:,:)
        if (i/=(Npl+1)) then
          j0 = j0+18
          j1 = j1+18
          hk(i0:i1,j0:j1) = h02(:,:)
          hk(j0:j1,i0:i1) = h20(:,:)
        end if
      end if
    end select
  end do

  return
end subroutine hamiltklinearsoc