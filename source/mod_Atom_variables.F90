module mod_Atom_variables
  implicit none
contains

  subroutine allocate_Atom_variables(nAtoms,nOrbs)
    !! This subroutine allocates variables that depend on nAtom
    use mod_mpi_pars,   only: abortProgram
    use mod_System,     only: System_type
#ifdef _GPU
    use mod_parameters, only: sigmai2i,sigmaimunu2i,sigmaijmunu2i,isigmamu2n,isigmamu2n_d,n2isigmamu
#else
    use mod_parameters, only: sigmai2i,sigmaimunu2i,sigmaijmunu2i,isigmamu2n,n2isigmamu
#endif
    implicit none
    integer, intent(in) :: nAtoms,nOrbs

    allocate( sigmai2i(4,nAtoms),sigmaimunu2i(4,nAtoms,nOrbs,nOrbs),sigmaijmunu2i(4,nAtoms,nAtoms,nOrbs,nOrbs),isigmamu2n(nAtoms,2,nOrbs), n2isigmamu(nAtoms*2*nOrbs,3))
    
#ifdef _GPU
    allocate( isigmamu2n_d(nAtoms,2,nOrbs) )
#endif

    ! pln_cnt(1) contains # of in-plane neighbors TODO: Re-Include
    ! allocate(sha_longitudinal(pln_cnt(1)),sha_transverse(pln_cnt(1)),long_cos(pln_cnt(1)),transv_cos(pln_cnt(1)))
    ! if (AllocateStatus/=0) call abortProgram("[main] Not enough memory for: sha_longitudinal,sha_transverse,long_cos,transv_cos")

  end subroutine allocate_Atom_variables

  subroutine deallocate_Atom_variables()
    !! This subroutine deallocates variables that depend on nAtom
#ifdef _GPU
    use mod_parameters, only: sigmai2i,sigmaimunu2i,sigmaijmunu2i,isigmamu2n,isigmamu2n_d,n2isigmamu
#else
    use mod_parameters, only: sigmai2i,sigmaimunu2i,sigmaijmunu2i,isigmamu2n,n2isigmamu
#endif
    use mod_magnet,     only: mzd0,mpd0,rho0,rhod0
    implicit none

    if(allocated(sigmai2i))      deallocate( sigmai2i )
    if(allocated(sigmaimunu2i))  deallocate( sigmaimunu2i )
    if(allocated(sigmaijmunu2i)) deallocate( sigmaijmunu2i )
    if(allocated(isigmamu2n))    deallocate( isigmamu2n )
    if(allocated(n2isigmamu))    deallocate( n2isigmamu )
    if(allocated(rho0))          deallocate( rho0 )
    if(allocated(rhod0))         deallocate( rhod0 )
    if(allocated(mzd0))          deallocate( mzd0 )
    if(allocated(mpd0))          deallocate( mpd0 )
#ifdef _GPU
    if(allocated(isigmamu2n_d))    deallocate( isigmamu2n_d )
#endif

    !deallocate(t0, t0i)
    !deallocate(sha_longitudinal,sha_transverse,long_cos,transv_cos)

  end subroutine deallocate_Atom_variables
end module mod_Atom_variables
