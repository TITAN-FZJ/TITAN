! Calculate green function (E - H)^-1
subroutine green(er,ei,sys,kp,gf)
  use mod_f90_kind,   only: double
  use mod_constants,  only: cZero,cOne
  use mod_System,     only: ia, System
  use mod_parameters, only: nOrb2, offset
  implicit none
  integer     :: i,j,d
  real(double), intent(in) :: er,ei,kp(3)
  type(System), intent(in) :: sys
  complex(double) :: ec
  complex(double),dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2) :: gslab,hk
  complex(double),dimension(nOrb2, nOrb2, sys%nAtoms, sys%nAtoms), intent(out)  :: gf

  d = sys%nAtoms * nOrb2

  ec    = cmplx(er,ei,double)

  gslab = cZero
  do i = 1, d
    gslab(i,i) = ec
  end do

  call hamiltk(sys,kp,hk)

  call zaxpy(d*d,-cOne,hk,1,gslab,1)
  !gslab(i,j) = gslab(i,j) - hk(i,j)
  call invers(gslab, d)

  ! Put the slab Green's function [A(nAtoms*18,nAtoms*18)] in the A(i,j,mu,nu) form
  !dir$ ivdep:loop
  do j = 1, sys%nAtoms
    !dir$ ivdep:loop
    do i = 1, sys%nAtoms
      gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    end do
  end do

end subroutine green


! Calculate green function (E - H)^-1 without SOC and the linear term G0.H_so.G0
subroutine greenlinearsoc(er,ei,sys,kp,g0,g0vsocg0)
  use mod_f90_kind,   only: double
  use mod_constants,  only: cZero, cOne
  use mod_parameters, only: nOrb2, offset
  use mod_System,     only: ia, System
  !use mod_magnet, only:
  implicit none
  integer     :: i,j,d
  real(double), intent(in) :: er,ei,kp(3)
  type(System), intent(in) :: sys
  complex(double) :: ec
  complex(double), dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2)  :: gslab0,hk,vsoc,temp,temp2
  complex(double), dimension(nOrb2,nOrb2,sys%nAtoms,sys%nAtoms), intent(out)  :: g0,g0vsocg0

  d = sys%nAtoms*nOrb2

  ec = cmplx(er,ei,double)

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
  do j=1,sys%nAtoms
    do i=1,sys%nAtoms
      g0(:,:,i,j) = gslab0(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
      g0vsocg0(:,:,i,j) = temp2(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    end do
  end do
end subroutine greenlinearsoc


! Calculate green function (E - H)^-1 with linear SOC: G = G0+G0.H_so.G0
subroutine greenlineargfsoc(er,ei,sys,kp,gf)
  use mod_f90_kind,   only: double
  use mod_constants,  only: cZero, cOne
  use mod_parameters, only: nOrb2, offset
  use mod_System,     only: ia, System
  implicit none
  integer     :: i,j,d
  real(double), intent(in) :: er,ei,kp(3)
  type(System), intent(in) :: sys
  complex(double) :: ec
  complex(double),dimension(sys%nAtoms*nOrb2, sys%nAtoms*nOrb2)  :: gslab,gslab0,hk,vsoc,temp
  complex(double),dimension(nOrb2, nOrb2,sys%nAtoms,sys%nAtoms),intent(out)  :: gf

  d = sys%nAtoms*nOrb2

  ec    = cmplx(er,ei,double)

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
  do j=1,sys%nAtoms
    do i=1,sys%nAtoms
      gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    end do
  end do
end subroutine greenlineargfsoc
