module mod_SOC
  use mod_f90_kind, only: double
  implicit none
  logical :: SOC
  !! Turn on/off SOC
  logical :: llineargfsoc = .false.
  logical :: llinearsoc = .false.
  !! Linear SOC
  real(double) :: socscale = 1.d0
  !! Rescale of SOC parameter
  complex(double), dimension(:,:,:), allocatable :: ls
  !! L.S matrix
contains

  subroutine updateLS(sys,theta, phi)
    use mod_System,            only: System
    use mod_f90_kind,          only: double
    use mod_constants,         only: pauli_mat, cZero
    use mod_magnet,            only: lvec
    use mod_rotation_matrices, only: rotation_matrix_ry, rotation_matrix_rz
    use mod_System,            only: System
    implicit none
    type(System), intent(in) :: sys
    real(double), intent(in) :: theta, phi
    real(double),dimension(3,3) :: ry,rz,rzy
    integer :: i,m,n,mu,nu,mup,nup

    call rotation_matrix_ry(theta,ry)
    call rotation_matrix_rz(phi,rz)

    call dgemm('n','n',3,3,3,1.d0,rz,3,ry,3,0.d0,rzy,3)

    ls = cZero
    do i=1,sys%nAtoms
      do n=1,3
        do m=1,3
          ! p-block
          do nu=2,4
            nup = nu+9
            do mu=2,4
              mup = mu+9
              ls(mu ,nu ,i) = ls(mu ,nu ,i) + sys%Types(sys%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,1,n)
              ls(mu ,nup,i) = ls(mu ,nup,i) + sys%Types(sys%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,2,n)
              ls(mup,nu ,i) = ls(mup,nu ,i) + sys%Types(sys%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,1,n)
              ls(mup,nup,i) = ls(mup,nup,i) + sys%Types(sys%Basis(i)%Material)%LambdaP*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,2,n)
            end do
          end do
          ! d-block
          do nu=5,9
            nup = nu+9
            do mu=5,9
              mup = mu+9
              ls(mu ,nu ,i) = ls(mu ,nu ,i) + sys%Types(sys%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,1,n)
              ls(mu ,nup,i) = ls(mu ,nup,i) + sys%Types(sys%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(1,2,n)
              ls(mup,nu ,i) = ls(mup,nu ,i) + sys%Types(sys%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,1,n)
              ls(mup,nup,i) = ls(mup,nup,i) + sys%Types(sys%Basis(i)%Material)%LambdaD*lvec(mu,nu,m)*rzy(m,n)*pauli_mat(2,2,n)
            end do
          end do

        end do
      end do
    end do

    ! rescale
    ls = 0.5d0 * socscale * ls

  end subroutine updateLS

  subroutine allocateLS(nAtoms,nOrb)
    use mod_constants, only: cZero
    use mod_mpi_pars,  only: abortProgram
    use mod_mpi_pars,  only: myrank
    implicit none
    integer, intent(in) :: nOrb,nAtoms
    if((myrank==0).and.(nOrb /= 9)) call abortProgram("[allocateLS] LS Matrix only implemented for nOrb = 9.")

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
