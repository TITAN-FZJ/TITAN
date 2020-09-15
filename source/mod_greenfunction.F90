module mod_greenfunction

contains
  ! Calculate green function (E - H)^-1
  subroutine green(er,ei,sys,kp,gf)
    use mod_kind,              only: dp
    use mod_constants,         only: cZero,cOne
    use mod_System,            only: ia,ia_sc,System_type
    use mod_parameters,        only: nOrb2,offset,dimH
    use mod_hamiltonian,       only: hamiltk
    use mod_superconductivity, only: lsupercond,superCond
    use mod_tools,             only: invers
    implicit none
    real(dp),          intent(in) :: er,ei,kp(3)
    type(System_type), intent(in) :: sys
    complex(dp), dimension(nOrb2*superCond,nOrb2*superCond,sys%nAtoms,sys%nAtoms), intent(out)  :: gf

    integer     :: i,j,d
    complex(dp) :: ec
    complex(dp),dimension(dimH*superCond, dimH*superCond) :: gslab,hk

    external :: zaxpy

    d = dimH * superCond

    ec    = cmplx(er,ei,dp)

    gslab = cZero
    do i=1,d
      gslab(i,i) = ec
    end do

    call hamiltk(sys,kp,hk)

    !zaxpy performs gslab = gslab + (-cOne)*hk
    !d*d si the dimension of the array
    !both 1's are for the strides to use on each array
    !for example we can add each second or third entry of the array
    !on this case we just use every component
    call zaxpy(d*d,-cOne,hk,1,gslab,1)
    !gslab(i,j) = gslab(i,j) - hk(i,j)
    call invers(gslab, d)


    if(.not.lsupercond) then
      ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
      do j = 1,sys%nAtoms
        do i = 1,sys%nAtoms
          gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
        end do
      end do
    else
      ! Put the slab Green's function [A(nAtoms*36,nAtoms*36)] in the A(i,j,mu,nu) form
      do j = 1,sys%nAtoms
        do i = 1,sys%nAtoms
          gf(      1:  nOrb2,      1:  nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(1,j):ia_sc(2,j))
          gf(      1:  nOrb2,nOrb2+1:2*nOrb2,i,j) = gslab(ia_sc(1,i):ia_sc(2,i),ia_sc(3,j):ia_sc(4,j))
          gf(nOrb2+1:2*nOrb2,      1:  nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(1,j):ia_sc(2,j))
          gf(nOrb2+1:2*nOrb2,nOrb2+1:2*nOrb2,i,j) = gslab(ia_sc(3,i):ia_sc(4,i),ia_sc(3,j):ia_sc(4,j))
        end do
      end do
    end if

  end subroutine green


  ! Calculate green function (E - H)^-1 without SOC and the linear term G0.H_so.G0
  subroutine greenlinearsoc(er,ei,sys,kp,g0,g0vsocg0)
    use mod_kind,        only: dp
    use mod_constants,   only: cZero, cOne
    use mod_parameters,  only: nOrb2, offset
    use mod_System,      only: ia, System_type
    use mod_hamiltonian, only: hamiltklinearsoc
    use mod_tools,       only: invers
    implicit none
    real(dp),          intent(in) :: er,ei,kp(3)
    type(System_type), intent(in) :: sys
    complex(dp), dimension(nOrb2,nOrb2,sys%nAtoms,sys%nAtoms), intent(out)  :: g0,g0vsocg0

    integer     :: i,j,d
    complex(dp) :: ec
    complex(dp), dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2)  :: gslab0,hk,vsoc,temp,temp2

    external :: zgemm

    d = sys%nAtoms*nOrb2

    ec = cmplx(er,ei,dp)

    gslab0 = cZero
    do i=1,d
     gslab0(i,i) = ec
    end do

    call hamiltklinearsoc(sys,kp,hk,vsoc)

    gslab0 = gslab0 - hk

    call invers(gslab0,d)

    ! Adding linear SOC
    call zgemm('n','n',d,d,d,cOne,vsoc,d,gslab0,d,cZero,temp,d)  ! temp = vsoc*gslab0
    call zgemm('n','n',d,d,d,cOne,gslab0,d,temp,d,cZero,temp2,d) ! temp2 = gslab0*vsoc*gslab0

    ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
    do j = 1,sys%nAtoms
      do i = 1,sys%nAtoms
        g0(:,:,i,j) = gslab0(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
        g0vsocg0(:,:,i,j) = temp2(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
      end do
    end do
  end subroutine greenlinearsoc


  ! Calculate green function (E - H)^-1 with linear SOC: G = G0+G0.H_so.G0
  subroutine greenlineargfsoc(er,ei,sys,kp,gf)
    use mod_kind,        only: dp
    use mod_constants,   only: cZero, cOne
    use mod_parameters,  only: nOrb2, offset
    use mod_System,      only: ia, System_type
    use mod_hamiltonian, only: hamiltklinearsoc
    use mod_tools,       only: invers
    implicit none
    real(dp),          intent(in) :: er,ei,kp(3)
    type(System_type), intent(in) :: sys
    complex(dp),dimension(nOrb2, nOrb2,sys%nAtoms,sys%nAtoms),intent(out)  :: gf

    integer     :: i,j,d
    complex(dp) :: ec
    complex(dp),dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2)  :: gslab,gslab0,hk,vsoc,temp

    external :: zgemm

    d = sys%nAtoms*nOrb2

    ec    = cmplx(er,ei,dp)

    temp   = cZero
    gslab0 = cZero
    do i=1,d
     temp(i,i)   = cOne
     gslab0(i,i) = ec
    end do

    call hamiltklinearsoc(sys,kp,hk,vsoc)

    gslab0 = gslab0 - hk

    call invers(gslab0,d)

    ! Adding linear SOC
    call zgemm('n','n',d,d,d,cOne,vsoc,d,gslab0,d,cOne,temp,d)   ! temp = 1+vsoc*gslab0
    call zgemm('n','n',d,d,d,cOne,gslab0,d,temp,d,cZero,gslab,d) ! gslab = gslab0(1+vsoc*gslab0)

    ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
    do j = 1,sys%nAtoms
      do i = 1,sys%nAtoms
        gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
      end do
    end do
  end subroutine greenlineargfsoc
end module mod_greenfunction