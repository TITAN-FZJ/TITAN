module mod_prefactors
  use mod_kind, only: dp
  implicit none
  complex(dp), dimension(:,:), allocatable  :: prefactor,prefactorlsoc
  complex(dp), dimension(:,:,:,:), allocatable :: lxpt,lypt,lzpt,tlxp,tlyp,tlzp

contains

  ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
  ! subroutine OAM_curr_hopping_times_L() !TODO: Re-Include
  !   use mod_magnet,        only: lxp, lyp, lzp
  !   use mod_parameters,    only: Npl, offset
  !   use mod_system,        only: n0sc1, n0sc2
  !   use mod_constants,     only: cZero
  !   implicit none
  !   integer         :: i,mu,nu,neighbor,alpha
  !
  !   lxpt = cZero
  !   lypt = cZero
  !   lzpt = cZero
  !   tlxp = cZero
  !   tlyp = cZero
  !   tlzp = cZero
  !   do nu=1,9
  !     do mu=1,9
  !       do neighbor=n0sc1,n0sc2
  !         do i=1,Npl
  !           do alpha=1,9
  !             lxpt(i,neighbor,mu,nu) = lxpt(i,neighbor,mu,nu) + lxp(mu,alpha,i)*t0i(nu,alpha,neighbor,i+offset)
  !             lypt(i,neighbor,mu,nu) = lypt(i,neighbor,mu,nu) + lyp(mu,alpha,i)*t0i(nu,alpha,neighbor,i+offset)
  !             lzpt(i,neighbor,mu,nu) = lzpt(i,neighbor,mu,nu) + lzp(mu,alpha,i)*t0i(nu,alpha,neighbor,i+offset)
  !             tlxp(i,neighbor,mu,nu) = tlxp(i,neighbor,mu,nu) + t0i(mu,alpha,neighbor,i+offset)*lxp(alpha,nu,i)
  !             tlyp(i,neighbor,mu,nu) = tlyp(i,neighbor,mu,nu) + t0i(mu,alpha,neighbor,i+offset)*lyp(alpha,nu,i)
  !             tlzp(i,neighbor,mu,nu) = tlzp(i,neighbor,mu,nu) + t0i(mu,alpha,neighbor,i+offset)*lzp(alpha,nu,i)
  !           end do
  !         end do
  !       end do
  !     end do
  !   end do
  !     ! end subroutine OAM_curr_hopping_times_L

  subroutine allocate_prefactors() !TODO:Re-Include
    use mod_parameters, only: dimens
    use mod_SOC, only: llinearsoc
    use mod_mpi_pars, only: abortProgram
    !use mod_system, only: n0sc1, n0sc2
    implicit none
    integer :: AllocateStatus

    if (.not. allocated(prefactor)) then
      allocate(prefactor(dimens,dimens), STAT = AllocateStatus  )
      if (AllocateStatus/=0) call abortProgram("[allocate_disturbances] Not enough memory for: prefactor")
    end if

    if (.not. allocated(prefactorlsoc)) then
      if(llinearsoc) then
        allocate(prefactorlsoc(dimens,dimens), STAT = AllocateStatus  )
        if (AllocateStatus/=0) call abortProgram("[allocate_disturbances] Not enough memory for: prefactorlsoc")
      end if
    end if

    ! TODO: Re-Include
    ! allocate(lxpt(Npl,n0sc1:n0sc2,9,9), &
    !          lypt(Npl,n0sc1:n0sc2,9,9), &
    !          lzpt(Npl,n0sc1:n0sc2,9,9), &
    !          tlxp(Npl,n0sc1:n0sc2,9,9), &
    !          tlyp(Npl,n0sc1:n0sc2,9,9), &
    !          tlzp(Npl,n0sc1:n0sc2,9,9), stat = AllocateStatus)
    ! if (AllocateStatus/=0) call abortProgram("[allocate_prefactors] Not enough memory for: lxpt,lypt,lzpt,tlxp,tlyp,tlzp")

  end subroutine allocate_prefactors

  subroutine deallocate_prefactors()
    implicit none
    if(allocated(prefactor)) deallocate(prefactor)
    if(allocated(prefactorlsoc)) deallocate(prefactorlsoc)

    if(allocated(lxpt)) deallocate(lxpt)
    if(allocated(lypt)) deallocate(lypt)
    if(allocated(lypt)) deallocate(lzpt)
    if(allocated(lypt)) deallocate(tlxp)
    if(allocated(lypt)) deallocate(tlyp)
    if(allocated(lypt)) deallocate(tlzp)

  end subroutine deallocate_prefactors

end module mod_prefactors
