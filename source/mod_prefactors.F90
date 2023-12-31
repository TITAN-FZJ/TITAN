module mod_prefactors
  !! This module stores the variables and procedures to calculate prefactors of current calculations
  !! ATTENTION: Since the currents are currently not working, this module is not used nor adapted to the new formats
  use mod_kind, only: dp
  implicit none
  complex(dp), dimension(:,:), allocatable  :: prefactor,prefactorlsoc
  complex(dp), dimension(:,:,:,:), allocatable :: lxpt,lypt,lzpt,tlxp,tlyp,tlzp

contains

  ! Calculates matrices hopping x angular momentum matrix for orbital angular momentum current calculation
  ! subroutine OAM_curr_hopping_times_L(s) !TODO: Re-Include
  !   
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
  !   do nu=1,s%nOrb
  !     do mu=1,s%nOrb
  !       do neighbor=n0sc1,n0sc2
  !         do i=1,s%nAtoms
  !           do alpha=1,s%nOrb
  !             lxpt(i,neighbor,mu,nu) = lxpt(i,neighbor,mu,nu) + s%Basis(i)%lpvec(mu,alpha,1)*t0i(nu,alpha,neighbor,i+offset)
  !             lypt(i,neighbor,mu,nu) = lypt(i,neighbor,mu,nu) + s%Basis(i)%lpvec(mu,alpha,2)*t0i(nu,alpha,neighbor,i+offset)
  !             lzpt(i,neighbor,mu,nu) = lzpt(i,neighbor,mu,nu) + s%Basis(i)%lpvec(mu,alpha,3)*t0i(nu,alpha,neighbor,i+offset)
  !             tlxp(i,neighbor,mu,nu) = tlxp(i,neighbor,mu,nu) + t0i(mu,alpha,neighbor,i+offset)*s%Basis(i)%lpvec(alpha,nu,1)
  !             tlyp(i,neighbor,mu,nu) = tlyp(i,neighbor,mu,nu) + t0i(mu,alpha,neighbor,i+offset)*s%Basis(i)%lpvec(alpha,nu,2)
  !             tlzp(i,neighbor,mu,nu) = tlzp(i,neighbor,mu,nu) + t0i(mu,alpha,neighbor,i+offset)*s%Basis(i)%lpvec(alpha,nu,3)
  !           end do
  !         end do
  !       end do
  !     end do
  !   end do
  !     ! end subroutine OAM_curr_hopping_times_L

  subroutine allocate_prefactors() !TODO:Re-Include
    use mod_parameters, only: dimens
    use mod_SOC,        only: llinearsoc
    use mod_mpi_pars,   only: abortProgram
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
    ! allocate(lxpt(Npl,n0sc1:n0sc2,s%nOrb,s%nOrb), &
    !          lypt(Npl,n0sc1:n0sc2,s%nOrb,s%nOrb), &
    !          lzpt(Npl,n0sc1:n0sc2,s%nOrb,s%nOrb), &
    !          tlxp(Npl,n0sc1:n0sc2,s%nOrb,s%nOrb), &
    !          tlyp(Npl,n0sc1:n0sc2,s%nOrb,s%nOrb), &
    !          tlzp(Npl,n0sc1:n0sc2,s%nOrb,s%nOrb), stat = AllocateStatus)
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
