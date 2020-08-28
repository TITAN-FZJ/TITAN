!! This subroutine mounts the conversion matrices from 4 to 2 ranks
subroutine initConversionMatrices(nAtoms, nOrbs)
  use mod_parameters, only: sigmai2i, sigmaimunu2i, sigmaijmunu2i, isigmamu2n, n2isigmamu
  implicit none
  integer, intent(in) :: nAtoms, nOrbs
  integer :: nu, mu, i, sigma, j, kount

  !------------------------- Conversion arrays  --------------------------
  do nu = 1, nOrbs
    do mu = 1, nOrbs
      do i = 1, nAtoms
        do sigma = 1, 4
          sigmaimunu2i(sigma,i,mu,nu) = (sigma-1)*nAtoms*nOrbs*nOrbs + (i-1)*nOrbs*nOrbs + (mu-1)*nOrbs + nu
          do j = 1, nAtoms
            sigmaijmunu2i(sigma,i,j,mu,nu) = (sigma-1)*nAtoms*nAtoms*nOrbs*nOrbs + (i-1)*nAtoms*nOrbs*nOrbs + (j-1)*nOrbs*nOrbs + (mu-1)*nOrbs + nu
          end do
        end do
      end do
    end do
  end do

  do i = 1, nAtoms
    do sigma = 1, 4
      sigmai2i(sigma,i) = (sigma-1)*nAtoms + i
    end do
  end do

  ! Conversion array from local atomic orbitals to (i,sigma,mu) to eigenvectors n
  do i = 1, nAtoms
    do sigma = 1, 2
      do mu = 1, nOrbs
        kount = (i-1)*2*nOrbs + (sigma-1)*nOrbs + mu
        isigmamu2n(i,sigma,mu) = kount
        n2isigmamu(kount,1) = i
        n2isigmamu(kount,2) = sigma
        n2isigmamu(kount,3) = mu
      end do
    end do
  end do

end subroutine initConversionMatrices
