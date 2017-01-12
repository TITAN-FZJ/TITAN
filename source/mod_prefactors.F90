module mod_prefactors
  use mod_f90_kind
  implicit none
  complex(double), dimension(:,:), allocatable  :: prefactor,prefactorlsoc
  complex(double),dimension(:,:,:,:), allocatable :: lxpt,lypt,lzpt,tlxp,tlyp,tlzp

contains

  ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
  subroutine OAM_curr_hopping_times_L()
    use mod_magnet, only: lxp, lyp, lzp
    use mod_tight_binding
    use mod_parameters, only: Npl, n0sc1, n0sc2
    use mod_constants, only: zero
    integer         :: i,mu,nu,neighbor,alpha

    lxpt = zero
    lypt = zero
    lzpt = zero
    tlxp = zero
    tlyp = zero
    tlzp = zero
    do nu=1,9 ; do mu=1,9 ; do neighbor=n0sc1,n0sc2 ; do i=1,Npl ; do alpha=1,9
        lxpt(i,neighbor,mu,nu) = lxpt(i,neighbor,mu,nu) + lxp(mu,alpha)*t00(i+1,neighbor,alpha,nu)
        lypt(i,neighbor,mu,nu) = lypt(i,neighbor,mu,nu) + lyp(mu,alpha)*t00(i+1,neighbor,alpha,nu)
        lzpt(i,neighbor,mu,nu) = lzpt(i,neighbor,mu,nu) + lzp(mu,alpha)*t00(i+1,neighbor,alpha,nu)
        tlxp(i,neighbor,mu,nu) = tlxp(i,neighbor,mu,nu) + t00(i+1,neighbor,alpha,mu)*lxp(alpha,nu)
        tlyp(i,neighbor,mu,nu) = tlyp(i,neighbor,mu,nu) + t00(i+1,neighbor,alpha,mu)*lyp(alpha,nu)
        tlzp(i,neighbor,mu,nu) = tlzp(i,neighbor,mu,nu) + t00(i+1,neighbor,alpha,mu)*lzp(alpha,nu)
    end do ; end do ; end do ; end do ; end do
    return
  end subroutine OAM_curr_hopping_times_L

  subroutine allocate_prefactors()
    use mod_parameters, only: Npl, n0sc1, n0sc2
    use mod_constants, only: zero

    allocate(lxpt(Npl,n0sc1:n0sc2,9,9),lypt(Npl,n0sc1:n0sc2,9,9),lzpt(Npl,n0sc1:n0sc2,9,9),tlxp(Npl,n0sc1:n0sc2,9,9),tlyp(Npl,n0sc1:n0sc2,9,9),tlzp(Npl,n0sc1:n0sc2,9,9))
    return
  end subroutine allocate_prefactors

  subroutine deallocate_prefactors()

    deallocate(lxpt,lypt,lzpt,tlxp,tlyp,tlzp)
    return
  end subroutine deallocate_prefactors

end module mod_prefactors
