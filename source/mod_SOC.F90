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
contains

  subroutine updateLS(s,theta,phi)
    use mod_kind,              only: dp
    use mod_System,            only: System_type
    use mod_constants,         only: pauli_mat
    use mod_rotation_matrices, only: rotation_matrix_ry,rotation_matrix_rz
    implicit none
    type(System_type), intent(inout) :: s
    real(dp),          intent(in)    :: theta, phi
    real(dp),         dimension(3,3) :: ry,rz,rzy
    integer     :: i,m,n,mu,nu,mup,nup,mupp,nupp
    complex(dp) :: l_el

    external :: dgemm

    call rotation_matrix_ry(theta,ry)
    call rotation_matrix_rz(phi,rz)

    call dgemm('n','n',3,3,3,1._dp,rz,3,ry,3,0._dp,rzy,3)

    sites: do i=1,s%nAtoms
      cart1: do n=1,3
        cart2: do m=1,3
          ! p block
          orbp1: do nupp=1,s%Types(s%Basis(i)%Material)%npOrb
            nu = s%Types(s%Basis(i)%Material)%pOrbs(nupp)
            nup = nu+s%Types(s%Basis(i)%Material)%nOrb
            orbp2: do mupp=1,s%Types(s%Basis(i)%Material)%npOrb
              mu  = s%Types(s%Basis(i)%Material)%pOrbs(mupp)
              mup = mu+s%Types(s%Basis(i)%Material)%nOrb
              ! Matrix element of L
              l_el = s%Types(s%Basis(i)%Material)%lvec(mu,nu,m)
              
              s%Basis(i)%ls(mu ,nu ) = s%Basis(i)%ls(mu ,nu ) + s%Types(s%Basis(i)%Material)%LambdaP * l_el * rzy(m,n) * pauli_mat(1,1,n)
              s%Basis(i)%ls(mu ,nup) = s%Basis(i)%ls(mu ,nup) + s%Types(s%Basis(i)%Material)%LambdaP * l_el * rzy(m,n) * pauli_mat(1,2,n)
              s%Basis(i)%ls(mup,nu ) = s%Basis(i)%ls(mup,nu ) + s%Types(s%Basis(i)%Material)%LambdaP * l_el * rzy(m,n) * pauli_mat(2,1,n)
              s%Basis(i)%ls(mup,nup) = s%Basis(i)%ls(mup,nup) + s%Types(s%Basis(i)%Material)%LambdaP * l_el * rzy(m,n) * pauli_mat(2,2,n)
            end do orbp2
          end do orbp1
          ! d-block
          orbd1: do nupp=1,s%Types(s%Basis(i)%Material)%ndOrb
            nu = s%Types(s%Basis(i)%Material)%dOrbs(nupp)
            nup = nu+s%Types(s%Basis(i)%Material)%nOrb
            orbd2: do mupp=1,s%Types(s%Basis(i)%Material)%ndOrb
              mu  = s%Types(s%Basis(i)%Material)%dOrbs(mupp)
              mup = mu+s%Types(s%Basis(i)%Material)%nOrb
              ! Matrix element of L
              l_el = s%Types(s%Basis(i)%Material)%lvec(mu,nu,m)

              s%Basis(i)%ls(mu ,nu ) = s%Basis(i)%ls(mu ,nu ) + s%Types(s%Basis(i)%Material)%LambdaD * l_el * rzy(m,n) * pauli_mat(1,1,n)
              s%Basis(i)%ls(mu ,nup) = s%Basis(i)%ls(mu ,nup) + s%Types(s%Basis(i)%Material)%LambdaD * l_el * rzy(m,n) * pauli_mat(1,2,n)
              s%Basis(i)%ls(mup,nu ) = s%Basis(i)%ls(mup,nu ) + s%Types(s%Basis(i)%Material)%LambdaD * l_el * rzy(m,n) * pauli_mat(2,1,n)
              s%Basis(i)%ls(mup,nup) = s%Basis(i)%ls(mup,nup) + s%Types(s%Basis(i)%Material)%LambdaD * l_el * rzy(m,n) * pauli_mat(2,2,n)
            end do orbd2
          end do orbd1
        end do cart2
      end do cart1
      ! rescale
      s%Basis(i)%ls(:,:) = 0.5_dp * s%Basis(i)%ls(:,:)
    end do sites


  end subroutine updateLS

end module mod_SOC
