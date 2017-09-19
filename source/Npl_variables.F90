subroutine initConversionMatrices(nAtoms, nOrbs)
  use mod_parameters, only: sigmai2i, sigmaimunu2i, sigmaijmunu2i
  use mod_mpi_pars, only: abortProgram
  implicit none
  integer, intent(in) :: nAtoms, nOrbs
  integer :: AllocateStatus
  integer :: nu, mu, i, sigma, j

  allocate( sigmai2i(4,nAtoms),sigmaimunu2i(4,nAtoms,nOrbs,nOrbs),sigmaijmunu2i(4,nAtoms,nAtoms,nOrbs,nOrbs),STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: sigmai2i,sigmaimunu2i,sigmaijmunu2i")
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
  return
end subroutine initConversionMatrices

! This subroutine allocates variables that depend on Npl
subroutine allocate_Npl_variables(nAtoms)
  use mod_parameters, only: mmlayermag, U, layertype, mmlayer
  use mod_magnet, only: allocate_magnet_variables
  use TightBinding, only: nOrb
  use mod_mpi_pars, only: abortProgram
  implicit none
  integer, intent(in) :: nAtoms
  integer :: AllocateStatus

  call allocate_magnet_variables(nAtoms)

  allocate( mmlayer(nAtoms),layertype(nAtoms),U(nAtoms),mmlayermag(nAtoms), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: mmlayer,layertype,U,mmlayermag,lambda,npart0")

  ! pln_cnt(1) contains # of in-plane neighbors TODO: Re-Include
  ! allocate(sha_longitudinal(pln_cnt(1)),sha_transverse(pln_cnt(1)),long_cos(pln_cnt(1)),transv_cos(pln_cnt(1)))
  ! if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: sha_longitudinal,sha_transverse,long_cos,transv_cos")

  return
end subroutine allocate_Npl_variables

! This subroutine allocates variables that depend on Npl
subroutine deallocate_Npl_variables()
  use mod_parameters, only: sigmai2i, sigmaimunu2i, sigmaijmunu2i, mmlayer, layertype, U, mmlayermag
  use mod_magnet, only: deallocate_magnet_variables
  implicit none

  call deallocate_magnet_variables()

  deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i)
  deallocate(mmlayer,layertype,U,mmlayermag)
  !deallocate(t0, t0i)
  !deallocate(sha_longitudinal,sha_transverse,long_cos,transv_cos)

  return
end subroutine deallocate_Npl_variables
