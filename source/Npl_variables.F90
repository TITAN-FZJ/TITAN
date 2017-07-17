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
  use mod_magnet, only: eps1, mx, my, mz, mvec_cartesian, mvec_spherical, hdel, mp, &
                        lptheta, lpabs, lphi, ltheta, labs, mm, hdelm, mabs, mtheta, &
                        mphi, lpphi, hdelp
  use TightBinding, only: nOrb
  use mod_mpi_pars, only: abortProgram
  implicit none
  integer, intent(in) :: nAtoms
  integer :: AllocateStatus

  allocate(eps1(nAtoms), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: eps1")

  allocate( mx(nAtoms),my(nAtoms),mz(nAtoms),mvec_cartesian(nAtoms,3),mvec_spherical(nAtoms,3),hdel(nAtoms),mp(nAtoms),hdelp(nAtoms),mm(nAtoms),hdelm(nAtoms), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm")

  allocate( mabs(nAtoms),mtheta(nAtoms),mphi(nAtoms),labs(nAtoms),ltheta(nAtoms),lphi(nAtoms),lpabs(nAtoms),lptheta(nAtoms),lpphi(nAtoms), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi")

  allocate( mmlayer(nAtoms),layertype(nAtoms),U(nAtoms),mmlayermag(nAtoms), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: mmlayer,layertype,U,mmlayermag,lambda,npart0")

  ! pln_cnt(1) contains # of in-plane neighbors TODO: Re-Include
  ! allocate(sha_longitudinal(pln_cnt(1)),sha_transverse(pln_cnt(1)),long_cos(pln_cnt(1)),transv_cos(pln_cnt(1)))
  ! if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: sha_longitudinal,sha_transverse,long_cos,transv_cos")

  return
end subroutine allocate_Npl_variables

! This subroutine allocates variables that depend on Npl
subroutine deallocate_Npl_variables()
  use mod_parameters
  use mod_magnet
  implicit none

  deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1)
  deallocate(mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm)
  deallocate(mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi)
  if(lGSL) deallocate(lxm,lym,lzm,lxpm,lypm,lzpm)
  deallocate(mmlayer,layertype,U,mmlayermag)
  deallocate(hhwx,hhwy,hhwz,sb,lb)
  !deallocate(t0, t0i)
  deallocate(sha_longitudinal,sha_transverse,long_cos,transv_cos)

  return
end subroutine deallocate_Npl_variables
