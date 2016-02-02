! Calculate green function of a slab containing Npl W layers
!  S    S-1   S-2       S-2   S-1    S
!  |-----|-----|---...---|-----|-----|
!  1     2     3       Npl-2 Npl-1  Npl
!   <-1-> <-2->           <-2-> <-1->
subroutine green(er,ei,kp,gf)
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_magnet
  implicit none
  integer     :: i,j,i0,i1,j0,j1,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension((Npl+2)*18,(Npl+2)*18) :: gslab,hk
  complex(double),dimension(Npl,Npl,18,18),intent(out)  :: gf

  d=(Npl+2)*18

  ec    = cmplx(er,ei,double)

  gslab = zero
  do i=1,d
   gslab(i,i) = ec
  end do

  call hamiltk(kp,hk)

  gslab = gslab - hk

  call invers(gslab,d)

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  do j=1,Npl
    do i=1,Npl
      i0 = i*18+1
      i1 = i0+17
      j0 = j*18+1
      j1 = j0+17
      gf(i,j,:,:) = gslab(i0:i1,j0:j1)
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
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_magnet
  implicit none
  integer     :: i,j,i0,i1,j0,j1,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension((Npl+2)*18,(Npl+2)*18)  :: gslab0,hk,vsoc,temp,temp2
  complex(double),dimension(Npl,Npl,18,18),intent(out)  :: g0,g0vsocg0

  d = (Npl+2)*18

  ec    = cmplx(er,ei,double)

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
  do j=1,Npl
    do i=1,Npl
      i0 = i*18+1
      i1 = i0+17
      j0 = j*18+1
      j1 = j0+17
      g0(i,j,:,:) = gslab0(i0:i1,j0:j1)
      g0vsocg0(i,j,:,:) = temp2(i0:i1,j0:j1)
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
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_magnet
  implicit none
  integer     :: i,j,i0,i1,j0,j1,d
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension((Npl+2)*18,(Npl+2)*18)  :: gslab,gslab0,hk,vsoc,temp
  complex(double),dimension(Npl,Npl,18,18),intent(out)  :: gf

  d = (Npl+2)*18

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
  do j=1,Npl
    do i=1,Npl
      i0 = i*18+1
      i1 = i0+17
      j0 = j*18+1
      j1 = j0+17
      gf(i,j,:,:) = gslab(i0:i1,j0:j1)
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
  use mod_f90_kind
  use mod_constants
  use mod_parameters
  use mod_magnet
  implicit none
  integer     :: i,j,i0,i1,j0,j1
  real(double), intent(in)  :: er,ei,kp(3)
  complex(double) :: ec
  complex(double),dimension((Npl+2)*18,(Npl+2)*18)  :: gslab,identes,hk
  complex(double),dimension(Npl+2,Npl+2,18,18),intent(out)  :: gf

  ec    = cmplx(er,ei,double)
  identes     = zero
  do i=1,(Npl+2)*18
   identes(i,i) = zum
  end do

  gslab = zero

  call hamiltk(kp,hk)

  gslab = ec*identes - hk

  call invers(gslab,(Npl+2)*18)

  ! Put the slab Green's function [A(Npl*18,Npl*18)] in the A(i,j,mu,nu) form
  do j=1,Npl+2
    do i=1,Npl+2
      i0 = (i-1)*18+1
      i1 = i0+17
      j0 = (j-1)*18+1
      j1 = j0+17
      gf(i,j,:,:) = gslab(i0:i1,j0:j1)
    end do
  end do

  return
end subroutine green_es
