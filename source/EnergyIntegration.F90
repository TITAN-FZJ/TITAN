module EnergyIntegration
  use mod_kind, only: dp
  use mod_System, only: s => sys
  implicit none

  integer :: pn1, pn2, pnt
  integer :: parts, parts3
  integer :: n1gl, n3gl

!========================================================================================!
! Energy integration variables
! Number of parts to divide energy integral I1+I2 and I3
  integer :: nepoints
  real(dp), allocatable :: y(:)
  !! Energy integration points along imaginary axis
  real(dp), allocatable :: wght(:)
  !! Weights of energy integration points along imaginary axis

  real(dp), allocatable :: x2(:)
  !! Energy integration points along real axis
  real(dp), allocatable :: p2(:)
  !! Weights of energy integration points along real axis

contains

  subroutine allocate_energy_points()
  !! Allocating arrays for real and imaginary energy integrations
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer :: AllocateStatus
    allocate( y(pn1), &
              wght(pn1), &
              x2(pn2), &
              p2(pn2), STAT = AllocateStatus )
    if(AllocateStatus/=0) call abortProgram("[allocate_energy_points] Not enough memory for: x1, p1, y, wght, x2, p2")
  end subroutine allocate_energy_points

  subroutine deallocate_energy_points()
  !! Allocating arrays for real and imaginary energy integrations
    implicit none
    deallocate( y,wght,x2,p2 )
  end subroutine deallocate_energy_points


  subroutine generate_imag_epoints()
  !! Generating integration points of the complex energy integral
  !! Generates points from zero until infinity
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer                       :: i,j,k,AllocateStatus
    real(dp)                  :: e1,e2,et1,xx, tmp1, tmp2
    real(dp),allocatable      :: x1(:),p1(:)

    allocate(x1(pn1),p1(pn1), STAT=AllocateStatus)
    if(AllocateStatus/=0) call abortProgram("[generate_imag_epoints] Not enough memory: x1,p1")

    do k=parts,1,-1
      if(parts/=1)then
        if(k==1) then
          e1 = 10._dp**(-2)
          e2 = 1._dp
        else if(k/=parts) then
          e1 = 10._dp**(-k-1)
          e2 = 10._dp**(-k)
        else
          e1 = 0._dp
          e2 = 10._dp**(-k)
        end if
      else
        e1 = 0._dp
        e2 = 1._dp
      end if
      i = (k-1)*n1gl+1
      j = k*n1gl
      call gauleg(e1,e2,x1(i:j),p1(i:j),n1gl)
    end do
    et1 = 1._dp
    do k = 1, pn1
      y(k)  = x1(k)/(1._dp-x1(k))
      xx    = (1._dp-x1(k))*(1._dp-x1(k))
      wght(k) = et1*p1(k)/xx
    end do
    deallocate(x1,p1)

    do i = 1, pn1
      do j = i+1, pn1
         if(y(i) > y(j)) then
            tmp1 = y(i)
            y(i) = y(j)
            y(j) = tmp1
            tmp2 = wght(i)
            wght(i) = wght(j)
            wght(j) = tmp2
         end if
      end do
   end do
  end subroutine generate_imag_epoints

  subroutine generate_real_epoints(e)
  !! Generating energy points in the real axis for third integration
  !! Generates points from Ef-e to Ef
    implicit none
    real(dp),intent(in) :: e
    !! Energy defining the lower boundary [Ef-e,Ef]
    integer :: i,j,k
    real(dp) :: e1,e2

    if (abs(e)>=1.e-15_dp) then
      nepoints = pnt
      do k=1,parts3
        e1 = s%Ef-((parts3-(k-1))*e/parts3)
        e2 = s%Ef-((parts3-k)*e/parts3)

        i = (k-1)*n3gl+1
        j = k*n3gl

        call gauleg(e1,e2,x2(i:j),p2(i:j),n3gl)
      end do
    else
      nepoints = pn1
    end if
  end subroutine generate_real_epoints

  ! This subroutine returns the N points in the abcissae x and its weight w using Gauss-Legendre
  subroutine gauleg(x1,x2,x,w,n)
    implicit none
    integer,      intent(in) :: n
    real(dp), intent(in) :: x1, x2
    real(dp), dimension(n), intent(out) :: x, w
    integer :: i, j, m
    real(dp) :: p1, p2, p3, pp, xl, xm, z, z1
    real(dp), parameter :: eps=3.e-14_dp

    m = (n+1)/2
    xm = 0.5_dp*(x2+x1)
    xl = 0.5_dp*(x2-x1)
    do i=1,m
      z = cos(3.141592654_dp*(i-0.25_dp)/(n+0.5_dp))
      z1 = 0.0
      do while(abs(z-z1)> eps)
        p1 = 1.0_dp
        p2 = 0.0_dp
        do j=1,n
          p3 = p2
          p2 = p1
          p1 = ((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1.0_dp)
        z1 = z
        z = z1 - p1/pp
      end do
      x(i) = xm - xl*z
      x(n+1-i) = xm + xl*z
      w(i) = (2.0_dp*xl)/((1.0_dp-z*z)*pp*pp)
      w(n+1-i) = w(i)
    end do
  end subroutine gauleg

end module EnergyIntegration
