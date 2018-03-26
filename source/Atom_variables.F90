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
end subroutine initConversionMatrices

! This subroutine allocates variables that depend on nAtom
subroutine allocate_Atom_variables(nAtoms)
  use mod_f90_kind, only: double
  use mod_parameters, only: mmlayermag, U, layertype, mmlayer, nmaglayers
  use TightBinding, only: nOrb
  use mod_mpi_pars, only: abortProgram
  implicit none
  integer, intent(in) :: nAtoms
  real(double) :: U_tmp
  integer :: AllocateStatus
  integer :: i

  if(size(U) == 1) then
    U_tmp = U(1)
    deallocate(U)
    allocate(U(nAtoms))
    U = U_tmp
  else if(size(U) /= nAtoms) then
    call abortProgram("[allocate_Atom_variables] U has wrong size")
  end if


  allocate( mmlayer(nAtoms),layertype(nAtoms),mmlayermag(nAtoms), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: mmlayer,layertype,mmlayermag")

  do i = 1, nAtoms
    if(abs(U(i)) > 1.d-9) then
      layertype(i) = 2
      nmaglayers = nmaglayers + 1
      mmlayermag(nmaglayers) = i
    else
      layertype(i) = 0
    end if
  end do
  ! pln_cnt(1) contains # of in-plane neighbors TODO: Re-Include
  ! allocate(sha_longitudinal(pln_cnt(1)),sha_transverse(pln_cnt(1)),long_cos(pln_cnt(1)),transv_cos(pln_cnt(1)))
  ! if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: sha_longitudinal,sha_transverse,long_cos,transv_cos")

end subroutine allocate_Atom_variables

! This subroutine allocates variables that depend on nAtom
subroutine deallocate_Atom_variables()
  use mod_parameters, only: sigmai2i, sigmaimunu2i, sigmaijmunu2i, mmlayer, layertype, U, mmlayermag
  use mod_magnet, only: deallocate_magnet_variables,rho0,rhod0
  implicit none

  call deallocate_magnet_variables()

  deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i)
  deallocate(mmlayer,layertype,U,mmlayermag)
  deallocate(rho0,rhod0)

  !deallocate(t0, t0i)
  !deallocate(sha_longitudinal,sha_transverse,long_cos,transv_cos)

end subroutine deallocate_Atom_variables
