module mod_SOC
  use mod_kind, only: dp
  implicit none
  logical :: SOC
  !! Turn on/off SOC
  logical :: llineargfsoc = .false.
  logical :: llinearsoc = .false.
  !! Linear SOC
  real(dp) :: socscale = 1._dp
  !! Rescale of SOC parameter
  complex(dp), dimension(:,:,:), allocatable :: ls
  !! L.S matrix
contains

  subroutine updateLS(s,theta,phi)
    use mod_kind,              only: dp
    use mod_System,            only: System_type
    use mod_constants,         only: pauli_mat
    use mod_magnet,            only: lvec
    use mod_rotation_matrices, only: rotation_matrix_ry,rotation_matrix_rz
    implicit none
    type(System_type), intent(in) :: s
    real(dp),          intent(in) :: theta, phi
    real(dp),      dimension(3,3) :: ry,rz,rzy
    integer     :: i,m,n,mu,nu,mup,nup,mupp,nupp

    external :: dgemm

    call rotation_matrix_ry(theta,ry)
    call rotation_matrix_rz(phi,rz)

    call dgemm('n','n',3,3,3,1._dp,rz,3,ry,3,0._dp,rzy,3)

    sites: do i=1,s%nAtoms
      cart1: do n=1,3
        cart2: do m=1,3
          ! p block
          orbp1: do nupp=1,s%npOrb
            nu = s%pOrbs(nupp)
            nup = nu+s%nOrb
            orbp2: do mupp=1,s%npOrb
              mu  = s%pOrbs(mupp)
              mup = mu+s%nOrb
              ls(mu ,nu ,i) = ls(mu ,nu ,i) + s%Types(s%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,1,n)
              ls(mu ,nup,i) = ls(mu ,nup,i) + s%Types(s%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,2,n)
              ls(mup,nu ,i) = ls(mup,nu ,i) + s%Types(s%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,1,n)
              ls(mup,nup,i) = ls(mup,nup,i) + s%Types(s%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,2,n)
            end do orbp2
          end do orbp1
          ! d-block
          orbd1: do nupp=1,s%ndOrb
            nu = s%dOrbs(nupp)
            nup = nu+s%nOrb
            orbd2: do mupp=1,s%ndOrb
              mu  = s%dOrbs(mupp)
              mup = mu+s%nOrb
              ls(mu ,nu ,i) = ls(mu ,nu ,i) + s%Types(s%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,1,n)
              ls(mu ,nup,i) = ls(mu ,nup,i) + s%Types(s%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,2,n)
              ls(mup,nu ,i) = ls(mup,nu ,i) + s%Types(s%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,1,n)
              ls(mup,nup,i) = ls(mup,nup,i) + s%Types(s%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,2,n)
            end do orbd2
          end do orbd1
        end do cart2
      end do cart1
    end do sites

    ! rescale
    ls = 0.5_dp * ls

  end subroutine updateLS

  subroutine allocateLS(nAtoms,nOrb)
    use mod_constants, only: cZero
    implicit none
    integer, intent(in) :: nOrb,nAtoms

    if(allocated(ls)) deallocate(ls)
    allocate(ls(2*nOrb, 2*nOrb,nAtoms))

    ! the spin-orbit matrix
    ls = cZero

  end subroutine allocateLS
  

  subroutine deallocateLS()
    implicit none

    deallocate(ls)

  end subroutine deallocateLS

end module mod_SOC
