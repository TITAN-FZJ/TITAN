! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine green(er,ei,kp,gf)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero,cOne
  use mod_System, only: ia, s => sys
  use mod_parameters, only: offset
  use TightBinding, only: nOrb,nOrb2
  implicit none
  integer     :: i,j,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension(s%nAtoms*nOrb2, s%nAtoms*nOrb2) :: gslab,hk
  complex(double),dimension(nOrb2, nOrb2, s%nAtoms, s%nAtoms), intent(out)  :: gf

  d = s%nAtoms * 2 * nOrb

  ec    = cmplx(er,ei,double)

  gslab = cZero
  do i = 1, d
    gslab(i,i) = ec
  end do

  call hamiltk(kp,hk)

  call zaxpy(d*d,-cOne,hk,1,gslab,1)
  !gslab(i,j) = gslab(i,j) - hk(i,j)
  call invers(gslab, d)

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  !dir$ ivdep:loop
  do j = 1, s%nAtoms
    !dir$ ivdep:loop
    do i = 1, s%nAtoms
      gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    end do
  end do

  return
end subroutine green


! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine greenlinearsoc(er,ei,kp,g0,g0vsocg0)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cOne
  use mod_parameters, only: offset
  use TightBinding, only: nOrb,nOrb2
  use mod_System, only: ia, s => sys
  !use mod_magnet, only:
  implicit none
  integer     :: i,j,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double), dimension(s%nAtoms*nOrb2, s%nAtoms*nOrb2)  :: gslab0,hk,vsoc,temp,temp2
  complex(double), dimension(nOrb2,nOrb2,s%nAtoms,s%nAtoms), intent(out)  :: g0,g0vsocg0

  d = s%nAtoms*nOrb2

  ec = cmplx(er,ei,double)

  gslab0 = cZero
  do i=1,d
   gslab0(i,i) = ec
  end do

  call hamiltklinearsoc(kp,hk,vsoc)

  gslab0 = gslab0 - hk

  call invers(gslab0,d)

  ! Adding linear SOC
  call zgemm('n','n',d,d,d,cOne,vsoc,d,gslab0,d,cZero,temp,d)  ! temp = vsoc*gslab0
  call zgemm('n','n',d,d,d,cOne,gslab0,d,temp,d,cZero,temp2,d) ! temp2 = gslab0*vsoc*gslab0

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  do j=1,s%nAtoms
    do i=1,s%nAtoms
      g0(:,:,i,j) = gslab0(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
      g0vsocg0(:,:,i,j) = temp2(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    end do
  end do
  return
end subroutine greenlinearsoc


! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine greenlineargfsoc(er,ei,kp,gf)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cOne
  use mod_parameters, only: offset
  use mod_System, only: ia, s => sys
  use TightBinding, only: nOrb,nOrb2
  implicit none
  integer     :: i,j,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension(s%nAtoms*nOrb2, s%nAtoms*nOrb2)  :: gslab,gslab0,hk,vsoc,temp
  complex(double),dimension(nOrb2, nOrb2,s%nAtoms,s%nAtoms),intent(out)  :: gf

  d = s%nAtoms*nOrb2

  ec    = cmplx(er,ei,double)

  temp   = cZero
  gslab0 = cZero
  do i=1,d
   temp(i,i)   = cOne
   gslab0(i,i) = ec
  end do

  call hamiltklinearsoc(kp,hk,vsoc)

  gslab0 = gslab0 - hk

  call invers(gslab0,d)

  ! Adding linear SOC
  call zgemm('n','n',d,d,d,cOne,vsoc,d,gslab0,d,cOne,temp,d)   ! temp = 1+vsoc*gslab0
  call zgemm('n','n',d,d,d,cOne,gslab0,d,temp,d,cZero,gslab,d) ! gslab = gslab0(1+vsoc*gslab0)

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  do j=1,s%nAtoms
    do i=1,s%nAtoms
      gf(:,:,i,j) = gslab(ia(1,i+offset):ia(4,i+offset),ia(1,j+offset):ia(4,j+offset))
    end do
  end do
  return
end subroutine greenlineargfsoc

! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine green_es(er,ei,kp,gf)
  use mod_f90_kind, only: double
  use mod_constants, only: cZero, cOne
  use mod_System, only: ia, s => sys
  use TightBinding, only: nOrb,nOrb2
  implicit none
  integer     :: i,j
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension(s%nAtoms*nOrb2, s%nAtoms*nOrb2)  :: gslab,identes,hk
  complex(double),dimension(nOrb2, nOrb2, s%nAtoms, s%nAtoms),intent(out)  :: gf

  ec    = cmplx(er,ei,double)
  identes     = cZero
  do i=1,s%nAtoms*nOrb2
   identes(i,i) = cOne
  end do

  gslab = cZero

  call hamiltk(kp,hk)

  gslab = ec*identes - hk

  call invers(gslab,s%nAtoms*nOrb2)

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  do j=1,s%nAtoms
    do i=1,s%nAtoms
      gf(:,:,i,j) = gslab(ia(1,i):ia(4,i),ia(1,j):ia(4,j))
    end do
  end do

  return
end subroutine green_es
