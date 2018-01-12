module SK_TightBinding
  use mod_f90_kind, only: double
  implicit none

contains

  subroutine get_parameter(s, fermi_layer, nOrb, Orbitals)
    use mod_f90_kind, only: double
    use AtomTypes, only: NeighborIndex
    use mod_system, only: System
    use mod_parameters, only: Ef
    implicit none
    type(System), intent(inout) :: s
    integer, intent(in) :: fermi_layer
    integer, intent(in) :: nOrb
    logical, dimension(9) :: Orbitals
    integer :: i,j,k, mu, nu
    real(double), dimension(:,:), allocatable :: bp
    real(double) :: scale_factor
    type(NeighborIndex), pointer :: current
    nullify(current)
    allocate(bp(nOrb,nOrb))

    do i = 1, s%nTypes
      call read_Papa_2C_param(s%Types(i), s%nStages, nOrb)
    end do

    Ef = s%Types(fermi_layer)%FermiLevel
    do i = 1, s%nTypes
      do j = 1, nOrb
        s%Types(i)%onSite(j,j) = s%Types(i)%onSite(j,j) - s%Types(i)%FermiLevel + Ef
      end do
    end do

    ! Allocate & initialize Hopping variables
    do i = 1, s%nNeighbors
      allocate(s%Neighbors(i)%t0i(nOrb,nOrb,s%nAtoms))
      allocate(s%Neighbors(i)%isHopping(s%nAtoms))
      s%Neighbors(i)%t0i = 0.d0
      s%Neighbors(i)%isHopping = .false.
    end do

    do i = 1, s%nAtoms
      do j = 1, s%nAtoms
        do k = 1, s%nStages
          current => s%Basis(i)%NeighborList(k,j)%head
          do while(associated(current))
            call set_hopping_matrix(s%Neighbors(current%index)%t0i(:,:,i), &
                                    s%Neighbors(current%index)%dirCos(:,i), &
                                    s%Types(s%Basis(i)%Material)%Hopping(:,k), &
                                    s%Types(s%Basis(j)%Material)%Hopping(:,k), nOrb)
            s%Neighbors(current%index)%isHopping(i) = .true.

            ! Scaling law by Andersen et al. O.K. Andersen, O. Jepsen, Physica 91B, 317 (1977); O.K. Andersen, W. Close. H. Nohl, Phys. Rev. B17, 1209 (1978)
            ! Distance dependence of tight binding matrix elements is given by V = C * d^(-[l+l'+1])
            ! e.g for ss hopping distance dependence is d^-1, for sp hopping d^-2
            scale_factor = s%Types(s%Basis(i)%Material)%stage(k) / s%Neighbors(current%index)%Distance(i)
            do mu = 1, nOrb
              do nu = 1, nOrb
                s%Neighbors(current%index)%t0i(mu,nu,i) = s%Neighbors(current%index)%t0i(mu,nu,i) * scale_factor ** (mu + nu + 1) 
              end do
            end do
            current => current%next
          end do
        end do
      end do
    end do
  end subroutine get_parameter

  subroutine set_hopping_matrix(t0i, dirCos, t1, t2, nOrb)
    use mod_f90_kind, only: double
    implicit none
    real(double), dimension(nOrb,nOrb), intent(inout) :: t0i
    real(double), dimension(3), intent(in) :: dirCos
    real(double), dimension(10), intent(in) :: t1, t2
    integer, intent(in) :: nOrb

    integer :: i
    real(double), dimension(10) :: mix

    do i = 1, 10
      mix(i) = sign(sqrt(abs(t1(i)) * abs(t2(i))), t1(i) + t2(i))
    end do
    call intd(mix(1), mix(2), mix(3), mix(4), mix(5), mix(6), mix(7), mix(8), mix(9), mix(10), dirCos, t0i)

    return
  end subroutine set_hopping_matrix

  subroutine read_Papa_2C_param(material, nStages, nOrb)
    use mod_f90_kind, only: double
    use AtomTypes, only: AtomType
    use mod_mpi_pars, only: abortProgram
    implicit none
    type(AtomType), intent(inout) :: material
    integer, intent(in) :: nStages
    integer, intent(in) :: nOrb
    real(double), dimension(3) :: dens
    integer :: i, j, k, l, ios, line_count = 0
    character(200) :: line
    character(50) :: words(10)
    real(double), dimension(3) :: a1, a2, a3, vec
    real(double), dimension(4) :: on_site
    real(double) :: dist

    allocate(material%onSite(nOrb,nOrb))
    allocate(material%Hopping(10,nStages))

    open(unit=995594, file=trim(material%Name), status='old', iostat=ios)
    if(ios /= 0) call abortProgram("[read_Papa_2C_param] Error occured")
    line_count = 0

    ! Counting lines to determine whether the parameter are available for nn_stages next neighbours
    read(unit=995594, fmt='(A)', iostat=ios) line
    do while (ios == 0)
      if(len_trim(line) > 0 .and. len_trim(line) < 200) line_count = line_count + 1
      read(unit=995594, fmt='(A)', iostat=ios) line
    end do
    rewind(995594)

    if(line_count < 10*nStages) call abortProgram("[read_Papa_2C_param] Parameter File corrupt")

    read(unit=995594, fmt='(A)', iostat = ios) line
    read(unit=line, fmt=*, iostat=ios) material%LatticeConstant

    read(unit=995594, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) material%Lambda

    read(unit=995594, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) material%FermiLevel

    read(unit=995594, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) dens(1), dens(2), dens(3)
    material%Occupation = dens(1)+dens(2)+dens(3)


    read(unit=995594, fmt='(A)', iostat=ios) line
    read(unit=line, fmt=*, iostat=ios) (a1(i), i=1,3)
    if(dot_product(a1,a1) <= 1.d-9) call abortProgram("[read_Papa_2C_param] a1 not given properly")
    a1 = a1 * material%LatticeConstant

    read(unit=995594, fmt='(A)', iostat=ios) line
    read(unit=line, fmt=*, iostat=ios) (a2(i), i=1,3)
    if(dot_product(a2,a2) <= 1.d-9) call abortProgram("[read_Papa_2C_param] a1 not given properly")
    a2 = a2 * material%LatticeConstant

    read(unit=995594, fmt='(A)', iostat=ios) line
    read(unit=line, fmt=*, iostat=ios) (a3(i), i=1,3)
    if(dot_product(a3,a3) <= 1.d-9) call abortProgram("[read_Papa_2C_param] a1 not given properly")
    a3 = a3 * material%LatticeConstant

    do i = 1, 4
      read(unit=995594, fmt='(A)', iostat = ios) line
      read(unit=line, fmt=*, iostat=ios) (words(j), j=1,10)
      read(unit=words(3), fmt=*, iostat=ios) on_site(i)
    end do

    material%onSite = 0.d0
    material%onSite(1,1) = on_site(1)
    do j=2,4
       material%onSite(j,j) = on_site(2)
    end do
    do j=5,7
       material%onSite(j,j) = on_site(3)
    end do
    do j=8,9
       material%onSite(j,j) = on_site(4)
    end do

    do i = 1, nstages
      do j = 1, 10
        read(unit=995594, fmt='(A)', iostat = ios) line
        read(unit=line, fmt=*, iostat=ios) (words(k), k=1,10)
        read(unit=words(4), fmt=*, iostat=ios) material%Hopping(j,i)
        !material%Hopping(j,i) = material%Hopping(j,i) * (a0_corr ** exponent(j)) ! Correction of hopping parameter by scaling law.
      end do
    end do

    allocate(material%stage(nStages))
    material % stage = 10.d0 * material % LatticeConstant
    do i = 0, 3*nStages
      do j = 0, 3*nStages
        do k = 0, 3*nStages
          if(i == 0 .and. j == 0 .and. k == 0) cycle
          vec = i * a1 + j * a2 + k * a3
          dist = sqrt(dot_product(vec,vec))
          l = nStages
          do while(1 <= l)
            if(material % stage(l) - dist < 1.d-9) exit
            if(l < nStages) material%stage(l+1) = material%stage(l)
            l = l - 1
          end do
          if(l == 0) then
            material%stage(l+1) = dist
          elseif(abs(material % stage(l) - dist) >= 1.d-9 .and. l < nStages) then
            material%stage(l+1) = dist
          end if

        end do
      end do
    end do
    close(995594)
  end subroutine

  pure subroutine intd(sss,pps,ppp,dds,ddp,ddd,sps,sds,pds,pdp,w,b)
    use mod_f90_kind, only: double
    use mod_constants, only: sq3
    implicit none
    real(double), intent(in)  :: sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd,w(3)
    real(double), intent(out) :: b(9,9)
    real(double)  :: x,y,z,xx,xy,yy,yz,zz,zx,xxyy,yyzz,zzxx
    real(double)  :: aux,aux1,aux2,aux3,aux4,r3,f1,f2,f3,f4,f5,f8,g1,g2
    x=w(1)
    y=w(2)
    z=w(3)
    xx=x*x
    xy=x*y
    yy=y*y
    yz=y*z
    zz=z*z
    zx=z*x
    xxyy=xx*yy
    yyzz=yy*zz
    zzxx=zz*xx

    r3=sq3

    f1=xx+yy
    f2=xx-yy
    f3=zz-.5d0*f1
    f8=3.d0*zz-1.d0
    f4=.5d0*r3*f2*pds
    f5=.5d0*f8*pds

    g1=1.5d0*f2*dds
    g2=r3*f3*dds

    aux=pps-ppp
    aux1=r3*sds
    aux2=r3*xx*pds+(1.d0-2.d0*xx)*pdp
    aux3=(r3*yy*pds+(1.d0-2.d0*yy)*pdp)
    aux4=r3*zz*pds+(1.d0-2.d0*zz)*pdp

    b(1,1)=sss
    b(1,2)=x*sps
    b(1,3)=y*sps
    b(1,4)=z*sps
    b(1,5)=xy*aux1
    b(1,6)=yz*aux1
    b(1,7)=zx*aux1
    b(1,8)=.5d0*f2*aux1
    b(1,9)=.5d0*f8*sds

    b(2,1)=-b(1,2)
    b(2,2)=xx*pps+(1.d0-xx)*ppp
    b(2,3)=xy*aux
    b(2,4)=zx*aux
    b(2,5)=aux2*y
    b(2,6)=(r3*pds-2.d0*pdp)*xy*z
    b(2,7)=aux2*z
    b(2,8)=(f4+(1.d0-f2)*pdp)*x
    b(2,9)=(f5-r3*zz*pdp)*x

    b(3,1)=-b(1,3)
    b(3,2)= b(2,3)
    b(3,3)=yy*pps+(1.d0-yy)*ppp
    b(3,4)=yz*aux
    b(3,5)=aux3*x
    b(3,6)=aux3*z
    b(3,7)= b(2,6)
    b(3,8)=(f4-(1.d0+f2)*pdp)*y
    b(3,9)=(f5-r3*zz*pdp)*y

    b(4,1)=-b(1,4)
    b(4,2)= b(2,4)
    b(4,3)= b(3,4)
    b(4,4)=zz*pps+(1.d0-zz)*ppp
    b(4,5)= b(2,6)
    b(4,6)=aux4*y
    b(4,7)=aux4*x
    b(4,8)=(f4-f2*pdp)*z
    b(4,9)=(f5+r3*f1*pdp)*z

    b(5,1)= b(1,5)
    b(5,2)=-b(2,5)
    b(5,3)=-b(3,5)
    b(5,4)=-b(4,5)
    b(5,5)=3.d0*xxyy*dds+(f1-4.d0*xxyy)*ddp+(zz+xxyy)*ddd
    b(5,6)=(3.d0*yy*dds+(1.d0-4.d0*yy)*ddp+(yy-1.d0)*ddd)*zx
    b(5,7)=(3.d0*xx*dds+(1.d0-4.d0*xx)*ddp+(xx-1.d0)*ddd)*yz
    b(5,8)=(g1-2.d0*f2*ddp+.5d0*f2*ddd)*xy
    b(5,9)=(g2-2.d0*r3*zz*ddp+.5d0*r3*(1.d0+zz)*ddd)*xy

    b(6,1)= b(1,6)
    b(6,2)=-b(2,6)
    b(6,3)=-b(3,6)
    b(6,4)=-b(4,6)
    b(6,5)= b(5,6)
    b(6,6)=3.d0*yyzz*dds+(yy+zz-4.d0*yyzz)*ddp+(xx+yyzz)*ddd
    b(6,7)=(3.d0*zz*dds+(1.d0-4.d0*zz)*ddp+(zz-1.d0)*ddd)*xy
    b(6,8)=(g1-(1.d0+2.d0*f2)*ddp+(1.d0+.5d0*f2)*ddd)*yz
    b(6,9)=(g2+r3*(f1-zz)*ddp-.5d0*r3*f1*ddd)*yz

    b(7,1)= b(1,7)
    b(7,2)=-b(2,7)
    b(7,3)=-b(3,7)
    b(7,4)=-b(4,7)
    b(7,5)= b(5,7)
    b(7,6)= b(6,7)
    b(7,7)=3.d0*zzxx*dds+(zz+xx-4.d0*zzxx)*ddp+(yy+zzxx)*ddd
    b(7,8)=(g1+(1.d0-2.d0*f2)*ddp-(1.d0-.5d0*f2)*ddd)*zx
    b(7,9)=(g2+r3*(f1-zz)*ddp-.5d0*r3*f1*ddd)*zx

    b(8,1)= b(1,8)
    b(8,2)=-b(2,8)
    b(8,3)=-b(3,8)
    b(8,4)=-b(4,8)
    b(8,5)= b(5,8)
    b(8,6)= b(6,8)
    b(8,7)= b(7,8)
    b(8,8)=.75d0*f2*f2*dds+(f1-f2*f2)*ddp+(zz+.25d0*f2*f2)*ddd
    b(8,9)=.5d0*f2*g2-r3*zz*f2*ddp+.25d0*r3*(1.d0+zz)*f2*ddd

    b(9,1)= b(1,9)
    b(9,2)=-b(2,9)
    b(9,3)=-b(3,9)
    b(9,4)=-b(4,9)
    b(9,5)= b(5,9)
    b(9,6)= b(6,9)
    b(9,7)= b(7,9)
    b(9,8)= b(8,9)
    b(9,9)=f3*f3*dds+3.d0*zz*f1*ddp+.75d0*f1*f1*ddd

    return
  end subroutine intd

end module SK_TightBinding
