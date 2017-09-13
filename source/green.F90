! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine green(er,ei,kp,gf)
  use mod_f90_kind, only: double
  use mod_constants, only: zero
  use mod_System, only: ia, s => sys
  use mod_parameters, only: offset
  use TightBinding, only: nOrb
  implicit none
  integer     :: i,j,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension(s%nAtoms*2*nOrb, s%nAtoms*2*nOrb) :: gslab,hk
  complex(double),dimension(2*nOrb, 2*nOrb, s%nAtoms, s%nAtoms), intent(out)  :: gf

  d = s%nAtoms * 2 * nOrb

  ec    = cmplx(er,ei,double)

  gslab = zero
  do i = 1, d
    gslab(i,i) = ec
  end do

  call hamiltk(kp,hk)

  gslab = gslab - hk

  call invers(gslab, d)

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  do j = 1, s%nAtoms
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
  use mod_constants, only: zero, zum
  use mod_parameters, only: offset
  use TightBinding, only: nOrb
  use mod_System, only: ia, s => sys
  !use mod_magnet, only:
  implicit none
  integer     :: i,j,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double), dimension(s%nAtoms*2*nOrb, s%nAtoms*2*nOrb)  :: gslab0,hk,vsoc,temp,temp2
  complex(double), dimension(2*nOrb,2*nOrb,s%nAtoms,s%nAtoms), intent(out)  :: g0,g0vsocg0

  d = s%nAtoms*2*nOrb

  ec = cmplx(er,ei,double)

  gslab0 = zero
  do i=1,d
   gslab0(i,i) = ec
  end do

  call hamiltklinearsoc(kp,hk,vsoc)

  gslab0 = gslab0 - hk

  call invers(gslab0,d)

  ! Adding linear SOC
  call zgemm('n','n',d,d,d,zum,vsoc,d,gslab0,d,zero,temp,d)  ! temp = vsoc*gslab0
  call zgemm('n','n',d,d,d,zum,gslab0,d,temp,d,zero,temp2,d) ! temp2 = gslab0*vsoc*gslab0

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
  use mod_constants, only: zero, zum
  use mod_parameters, only: offset
  use mod_System, only: ia, s => sys
  use TightBinding, only: nOrb
  implicit none
  integer     :: i,j,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension(s%nAtoms*2*nOrb, s%nAtoms*2*nOrb)  :: gslab,gslab0,hk,vsoc,temp
  complex(double),dimension(2*nOrb, 2*nOrb,s%nAtoms,s%nAtoms),intent(out)  :: gf

  d = s%nAtoms*2*nOrb

  ec    = cmplx(er,ei,double)

  temp   = zero
  gslab0 = zero
  do i=1,d
   temp(i,i)   = zum
   gslab0(i,i) = ec
  end do

  call hamiltklinearsoc(kp,hk,vsoc)

  gslab0 = gslab0 - hk

  call invers(gslab0,d)

  ! Adding linear SOC
  call zgemm('n','n',d,d,d,zum,vsoc,d,gslab0,d,zum,temp,d)   ! temp = 1+vsoc*gslab0
  call zgemm('n','n',d,d,d,zum,gslab0,d,temp,d,zero,gslab,d) ! gslab = gslab0(1+vsoc*gslab0)

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
  use mod_constants, only: zero, zum
  use mod_System, only: ia, s => sys
  use TightBinding, only: nOrb
  implicit none
  integer     :: i,j
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension(s%nAtoms*2*nOrb, s%nAtoms*2*nOrb)  :: gslab,identes,hk
  complex(double),dimension(2*nOrb, 2*nOrb, s%nAtoms, s%nAtoms),intent(out)  :: gf

  ec    = cmplx(er,ei,double)
  identes     = zero
  do i=1,s%nAtoms*2*nOrb
   identes(i,i) = zum
  end do

  gslab = zero

  call hamiltk(kp,hk)

  gslab = ec*identes - hk

  call invers(gslab,s%nAtoms*2*nOrb)

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  do j=1,s%nAtoms
    do i=1,s%nAtoms
      gf(:,:,i,j) = gslab(ia(1,i):ia(4,i),ia(1,j):ia(4,j))
    end do
  end do

  return
end subroutine green_es
