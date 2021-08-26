module mod_greenfunction

  procedure(calc_green_sub), pointer :: calc_green => green

  abstract interface
    subroutine calc_green_sub(er,ei,s,kp,gf)
      use mod_kind,              only: dp
      use mod_System,            only: System_type
      implicit none
      real(dp),          intent(in) :: er,ei,kp(3)
      type(System_type), intent(in) :: s
      complex(dp), dimension(s%nOrb2sc,s%nOrb2sc,s%nAtoms,s%nAtoms), intent(out)  :: gf
    end subroutine calc_green_sub
  end interface

contains

  ! Calculate green function (E - H)^-1
  subroutine green(er,ei,s,kp,gf)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,cOne
    use mod_System,            only: ia,ia_sc,System_type
    use mod_parameters,        only: dimHsc
    use mod_hamiltonian,       only: calchk,h0
    use mod_superconductivity, only: lsupercond,superCond
    use mod_tools,             only: invers
    implicit none
    real(dp),          intent(in) :: er,ei,kp(3)
    type(System_type), intent(in) :: s
    complex(dp), dimension(s%nOrb2*superCond,s%nOrb2*superCond,s%nAtoms,s%nAtoms), intent(out)  :: gf

    integer     :: i,j
    complex(dp) :: ec
    complex(dp),dimension(dimHsc,dimHsc) :: gslab,hk

    external :: zaxpy

    ec    = cmplx(er,ei,dp)

    gslab = cZero
    do i=1,dimHsc
      gslab(i,i) = ec
    end do

    hk = h0 + calchk(s,kp)

    !zaxpy performs gslab = gslab + (-cOne)*hk
    !dimHsc*dimHsc si the dimension of the array
    !both 1's are for the strides to use on each array
    !for example we can add each second or third entry of the array
    !on this case we just use every component
    call zaxpy(dimHsc*dimHsc,-cOne,hk,1,gslab,1)
    !gslab(i,j) = gslab(i,j) - hk(i,j)
    call invers(gslab, dimHsc)

    if(.not.lsupercond) then
      ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          gf(:,:,i,j) = gslab(ia(1,i):ia(4,i),ia(1,j):ia(4,j))
        end do
      end do
    else
      ! Put the slab Green's function [A(nAtoms*36,nAtoms*36)] in the A(i,j,mu,nu) form
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          gf(        1:  s%nOrb2,        1:  s%nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(1,j):ia_sc(2,j))
          gf(        1:  s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(3,j):ia_sc(4,j))
          gf(s%nOrb2+1:2*s%nOrb2,        1:  s%nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(1,j):ia_sc(2,j))
          gf(s%nOrb2+1:2*s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(3,j):ia_sc(4,j))
        end do
      end do
    end if

  end subroutine green


  ! ! Calculate green function (E - H)^-1
  ! subroutine green_fullhk(er,ei,s,kp,gf)
  !   use mod_kind,              only: dp
  !   use mod_constants,         only: cZero,cOne
  !   use mod_System,            only: ia,ia_sc,System_type
  !   use mod_parameters,        only: dimHsc
  !   use mod_hamiltonian,       only: calchk,h0
  !   use mod_superconductivity, only: lsupercond,superCond
  !   use mod_tools,             only: invers
  !   implicit none
  !   real(dp),          intent(in) :: er,ei,kp(3)
  !   type(System_type), intent(in) :: s
  !   complex(dp), dimension(s%nOrb2*superCond,s%nOrb2*superCond,s%nAtoms,s%nAtoms), intent(out)  :: gf

  !   integer     :: i,j
  !   complex(dp) :: ec
  !   complex(dp),dimension(dimHsc,dimHsc) :: gslab,hk

  !   external :: zaxpy

  !   ec    = cmplx(er,ei,dp)

  !   gslab = cZero
  !   do i=1,dimHsc
  !     gslab(i,i) = ec
  !   end do

  !   hk = h0 + fullhk(:,:,iz)

  !   !zaxpy performs gslab = gslab + (-cOne)*hk
  !   !dimHsc*dimHsc si the dimension of the array
  !   !both 1's are for the strides to use on each array
  !   !for example we can add each second or third entry of the array
  !   !on this case we just use every component
  !   call zaxpy(dimHsc*dimHsc,-cOne,hk,1,gslab,1)
  !   !gslab(i,j) = gslab(i,j) - hk(i,j)
  !   call invers(gslab, dimHsc)

  !   if(.not.lsupercond) then
  !     ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
  !     do j = 1,s%nAtoms
  !       do i = 1,s%nAtoms
  !         gf(:,:,i,j) = gslab(ia(1,i):ia(4,i),ia(1,j):ia(4,j))
  !       end do
  !     end do
  !   else
  !     ! Put the slab Green's function [A(nAtoms*36,nAtoms*36)] in the A(i,j,mu,nu) form
  !     do j = 1,s%nAtoms
  !       do i = 1,s%nAtoms
  !         gf(        1:  s%nOrb2,        1:  s%nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(1,j):ia_sc(2,j))
  !         gf(        1:  s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(3,j):ia_sc(4,j))
  !         gf(s%nOrb2+1:2*s%nOrb2,        1:  s%nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(1,j):ia_sc(2,j))
  !         gf(s%nOrb2+1:2*s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(3,j):ia_sc(4,j))
  !       end do
  !     end do
  !   end if

  ! end subroutine green_fullhk


  ! Calculate green function (E - H)^-1 with linear SOC: G = G0+G0.H_so.G0
  subroutine greenlineargfsoc(er,ei,s,kp,gf)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,cOne
    use mod_parameters,        only: dimHsc
    use mod_System,            only: ia,ia_sc,System_type
    use mod_hamiltonian,       only: hamiltklinearsoc
    use mod_tools,             only: invers
    use mod_superconductivity, only: lsupercond,superCond
    implicit none
    real(dp),          intent(in) :: er,ei,kp(3)
    type(System_type), intent(in) :: s
    complex(dp), dimension(s%nOrb2*superCond,s%nOrb2*superCond,s%nAtoms,s%nAtoms),intent(out)  :: gf

    integer     :: i,j
    complex(dp) :: ec
    complex(dp), dimension(dimHsc,dimHsc) :: gslab,gslab0,hk,vsoc,temp

    external :: zgemm

    ec    = cmplx(er,ei,dp)

    temp   = cZero
    gslab0 = cZero
    do i=1,dimHsc
     temp(i,i)   = cOne
     gslab0(i,i) = ec
    end do

    call hamiltklinearsoc(s,kp,hk,vsoc)

    gslab0 = gslab0 - hk

    call invers(gslab0,dimHsc)

    ! Adding linear SOC
    call zgemm('n','n',dimHsc,dimHsc,dimHsc,cOne,vsoc,dimHsc,gslab0,dimHsc,cOne,temp,dimHsc)   ! temp = 1+vsoc*gslab0
    call zgemm('n','n',dimHsc,dimHsc,dimHsc,cOne,gslab0,dimHsc,temp,dimHsc,cZero,gslab,dimHsc) ! gslab = gslab0(1+vsoc*gslab0)

    if(.not.lsupercond) then
      ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          gf(:,:,i,j) = gslab(ia(1,i):ia(4,i),ia(1,j):ia(4,j))
        end do
      end do
    else
      ! Put the slab Green's function [A(nAtoms*36,nAtoms*36)] in the A(i,j,mu,nu) form
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          gf(        1:  s%nOrb2,        1:  s%nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(1,j):ia_sc(2,j))
          gf(        1:  s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(3,j):ia_sc(4,j))
          gf(s%nOrb2+1:2*s%nOrb2,        1:  s%nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(1,j):ia_sc(2,j))
          gf(s%nOrb2+1:2*s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(3,j):ia_sc(4,j))
        end do
      end do
    end if
  end subroutine greenlineargfsoc



  ! Calculate green function (E - H)^-1 without SOC and the linear term G0.H_so.G0
  subroutine greenlinearsoc(er,ei,s,kp,g0,g0vsocg0)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,cOne
    use mod_parameters,        only: dimHsc
    use mod_System,            only: ia,ia_sc,System_type
    use mod_hamiltonian,       only: hamiltklinearsoc
    use mod_tools,             only: invers
    use mod_superconductivity, only: lsupercond,superCond
    implicit none
    real(dp),          intent(in) :: er,ei,kp(3)
    type(System_type), intent(in) :: s
    complex(dp), dimension(s%nOrb2*superCond,s%nOrb2*superCond,s%nAtoms,s%nAtoms), intent(out)  :: g0,g0vsocg0

    integer     :: i,j
    complex(dp) :: ec
    complex(dp), dimension(dimHsc,dimHsc)  :: gslab0,hk,vsoc,temp,temp2

    external :: zgemm

    ec = cmplx(er,ei,dp)

    gslab0 = cZero
    do i=1,dimHsc
     gslab0(i,i) = ec
    end do

    call hamiltklinearsoc(s,kp,hk,vsoc)

    gslab0 = gslab0 - hk

    call invers(gslab0,dimHsc)

    ! Adding linear SOC
    call zgemm('n','n',dimHsc,dimHsc,dimHsc,cOne,vsoc,dimHsc,gslab0,dimHsc,cZero,temp,dimHsc)  ! temp = vsoc*gslab0
    call zgemm('n','n',dimHsc,dimHsc,dimHsc,cOne,gslab0,dimHsc,temp,dimHsc,cZero,temp2,dimHsc) ! temp2 = gslab0*vsoc*gslab0


    if(.not.lsupercond) then
      ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          g0(:,:,i,j)      = gslab0(ia(1,i):ia(4,i),ia(1,j):ia(4,j))
          g0vsocg0(:,:,i,j) = temp2(ia(1,i):ia(4,i),ia(1,j):ia(4,j))
        end do
      end do
    else
      ! Put the slab Green's function [A(nAtoms*36,nAtoms*36)] in the A(i,j,mu,nu) form
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          g0(        1:  s%nOrb2,        1:  s%nOrb2,i,j) = gslab0(ia_sc(1,i):ia_sc(2,i),ia_sc(1,j):ia_sc(2,j))
          g0(        1:  s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab0(ia_sc(1,i):ia_sc(2,i),ia_sc(3,j):ia_sc(4,j))
          g0(s%nOrb2+1:2*s%nOrb2,        1:  s%nOrb2,i,j) = gslab0(ia_sc(3,i):ia_sc(4,i),ia_sc(1,j):ia_sc(2,j))
          g0(s%nOrb2+1:2*s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = gslab0(ia_sc(3,i):ia_sc(4,i),ia_sc(3,j):ia_sc(4,j))

          g0vsocg0(        1:  s%nOrb2,        1:  s%nOrb2,i,j) = temp2(ia_sc(1,i):ia_sc(2,i),ia_sc(1,j):ia_sc(2,j))
          g0vsocg0(        1:  s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = temp2(ia_sc(1,i):ia_sc(2,i),ia_sc(3,j):ia_sc(4,j))
          g0vsocg0(s%nOrb2+1:2*s%nOrb2,        1:  s%nOrb2,i,j) = temp2(ia_sc(3,i):ia_sc(4,i),ia_sc(1,j):ia_sc(2,j))
          g0vsocg0(s%nOrb2+1:2*s%nOrb2,s%nOrb2+1:2*s%nOrb2,i,j) = temp2(ia_sc(3,i):ia_sc(4,i),ia_sc(3,j):ia_sc(4,j))
        end do
      end do
    end if

  end subroutine greenlinearsoc

end module mod_greenfunction
