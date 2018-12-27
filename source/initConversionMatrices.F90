!! This subroutine mounts the conversion matrices from 4 to 2 ranks
subroutine initConversionMatrices(nAtoms, nOrbs)
  use mod_parameters, only: sigmai2i, sigmaimunu2i, sigmaijmunu2i
  implicit none
  integer, intent(in) :: nAtoms, nOrbs
  integer :: nu, mu, i, sigma, j

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
end subroutine initConversionMatrices
