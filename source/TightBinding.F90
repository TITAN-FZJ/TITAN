module TightBinding
  implicit none

contains

  subroutine initTightBinding(s)
    use mod_system,     only: System
    use mod_mpi_pars,   only: abortProgram
    use mod_parameters, only: tbmode, fermi_layer, nOrb
    implicit none
    type(System), intent(inout) :: s
    if(tbmode == 1) then
      call get_SK_parameter(s, fermi_layer, nOrb)
    else if(tbmode == 2) then
      call abortProgram("[initTightBinding] tbmode == 2 not implemented yet.")
      !call get_DFT_hopping()
    end if
  end subroutine initTightBinding


  subroutine get_SK_parameter(s, fermi_layer, nOrb)
    use mod_f90_kind, only: double
    use AtomTypes,    only: NeighborIndex
    use mod_system,   only: System
    implicit none
    type(System), intent(inout) :: s
    integer,      intent(in)    :: fermi_layer
    integer,      intent(in)    :: nOrb
    integer                     :: i,j,k,l
    real(double), dimension(:,:), allocatable :: bp
    real(double)                              :: scale_factor(2), mix(10,2)
    type(NeighborIndex), pointer              :: current
    real(double), dimension(10), parameter    :: expon = [1.0d0,3.0d0,3.0d0,5.0d0,5.0d0,5.0d0,2.0d0,3.0d0,4.0d0,4.0d0]
    nullify(current)
    allocate(bp(nOrb,nOrb))

    do i = 1, s%nTypes
      call readElementFile(s%Types(i), s%nStages, nOrb)
    end do

    s%Ef = s%Types(fermi_layer)%FermiLevel
    do i = 1, s%nTypes
      do j = 1, nOrb
        s%Types(i)%onSite(j,j) = s%Types(i)%onSite(j,j) - s%Types(i)%FermiLevel + s%Ef
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
      s%totalOccupation = s%totalOccupation + s%Types(s%Basis(i)%Material)%Occupation
      do j = 1, s%nAtoms
        do k = 1, s%nStages
          current => s%Basis(i)%NeighborList(k,j)%head
          do while(associated(current))
            scale_factor(1) = s%Types(s%Basis(i)%Material)%stage(k) / s%Neighbors(current%index)%Distance(i)
            scale_factor(2) = s%Types(s%Basis(j)%Material)%stage(k) / s%Neighbors(current%index)%Distance(i)
            do l = 1, 10
              mix(l,1) = s%Types(s%Basis(i)%Material)%Hopping(l,k) * scale_factor(1) ** expon(l)
              mix(l,2) = s%Types(s%Basis(j)%Material)%Hopping(l,k) * scale_factor(2) ** expon(l)
            end do
            call set_hopping_matrix(s%Neighbors(current%index)%t0i(:,:,i), &
                                    s%Neighbors(current%index)%dirCos(:,i), &
                                    mix(:,1), mix(:,2), nOrb)
            s%Neighbors(current%index)%isHopping(i) = .true.
            ! Scaling law by Andersen et al. O.K. Andersen, O. Jepsen, Physica 91B, 317 (1977); O.K. Andersen, W. Close. H. Nohl, Phys. Rev. B17, 1209 (1978)
            ! Distance dependence of tight binding matrix elements is given by V = C * d^(-[l+l'+1])
            ! e.g for ss hopping distance dependence is d^-1, for sp hopping d^-2
            ! do mu = 1, nOrb
            !   do nu = 1, nOrb
            !     s%Neighbors(current%index)%t0i(mu,nu,i) = s%Neighbors(current%index)%t0i(mu,nu,i) * scale_factor ** (mu + nu + 1)
            !   end do
            ! end do
            current => current%next
          end do
        end do
      end do
    end do
  end subroutine get_SK_parameter

  subroutine set_hopping_matrix(t0i, dirCos, t1, t2, nOrb)
    use mod_f90_kind,   only: double
    use mod_parameters, only: lsimplemix
    implicit none
    real(double), dimension(nOrb,nOrb), intent(inout) :: t0i
    real(double), dimension(3),         intent(in)    :: dirCos
    real(double), dimension(10),        intent(in)    :: t1, t2
    integer,                            intent(in)    :: nOrb
    integer                     :: i
    real(double), dimension(10) :: mix

    do i = 1, 10
      if(lsimplemix) then
        mix(i) = 0.5d0*(t1(i) + t2(i))
      else
        mix(i) = sign(sqrt(abs(t1(i)) * abs(t2(i))), t1(i) + t2(i))
      end if
    end do
    call intd(mix(1), mix(2), mix(3), mix(4), mix(5), mix(6), mix(7), mix(8), mix(9), mix(10), dirCos, t0i)

  end subroutine set_hopping_matrix

  !! Reading element file, including all the parameters
  subroutine readElementFile(material, nStages, nOrb)
    use mod_f90_kind, only: double
    use AtomTypes,    only: AtomType
    use mod_mpi_pars, only: abortProgram
    use mod_tools,    only: itos, vec_norm
    use mod_io,       only: log_warning
    implicit none
    type(AtomType), intent(inout) :: material
    integer,        intent(in)    :: nStages
    integer,        intent(in)    :: nOrb
    integer :: nTypes
    integer :: nn_stages
    integer :: nAtoms
    integer,      dimension(:), allocatable :: type_count
    real(double), dimension(3) :: dens
    integer            :: i, j, k, l, m, ios, line_count = 0
    integer, parameter :: line_length = 300, word_length = 50, max_elements = 50
    character(len=word_length), dimension(max_elements) :: str_arr
    character(200) :: line
    character(50)  :: words(10)
    real(double)                 :: dist
    real(double), dimension(3)   :: vec
    real(double), dimension(3,3) :: Bravais
    real(double), dimension(4)   :: on_site
    real(double), dimension(:,:), allocatable :: position
    real(double), dimension(:),   allocatable :: localDistances
    character(len=1)   :: coord_type
    character(len=100) :: Name
    integer :: f_unit = 995594

    allocate(material%onSite(nOrb,nOrb))
    allocate(material%Hopping(10,nStages))

    ! Opening file
    open(f_unit, file=trim(material%Name), status='old', iostat=ios)
    if(ios /= 0) call abortProgram("[readElementFile] Error occured")
    line_count = 0

    ! Read parameter name
    read(f_unit, fmt='(A)', iostat=ios) Name
    Name = trim(adjustl(Name))

    ! Reading lattice parameter a0
    read(f_unit, fmt='(A)', iostat=ios) line
    read(unit=line, fmt=*, iostat=ios) material%LatticeConstant

    ! Bravais lattice
    do j = 1, 3
      read(f_unit, fmt='(A)', iostat=ios) line
      read(unit=line, fmt=*, iostat=ios) (Bravais(i,j), i=1,3)
      if(vec_norm(Bravais(:,j),3) <= 1.d-9) &
        call log_warning("readElementFile", "Bravais vector a" // trim(itos(j)) // "'not given.")
    end do
    Bravais = Bravais * material%LatticeConstant
    material%a1 = Bravais(:,1)
    material%a2 = Bravais(:,2)
    material%a3 = Bravais(:,3)

    ! Read Different Elements in File
    str_arr = ""
    read(f_unit, fmt='(A)', iostat=ios) line
    read(unit=line, fmt=*, iostat=ios) (str_arr(i), i = 1, max_elements)
    nTypes = 0
    do i = 1, max_elements
      if(str_arr(i)(1:1) == "!") exit
      if(len_trim(str_arr(i)) == 0 .or. len_trim(str_arr(i)) == word_length) cycle
      nTypes = nTypes + 1
    end do

    ! Read amount of atoms per element
    allocate(type_count(nTypes))
    read(f_unit, fmt='(A)', iostat=ios) line
    read(unit=line, fmt=*, iostat=ios) (type_count(i), i = 1, nTypes)

    ! Count number of atoms
    nAtoms = sum(type_count(1:nTypes))
    if(nAtoms <= 0) call abortProgram("[readElementFile] No basis atoms given!")

    ! Read coordinate type
    read(f_unit, fmt='(A)', iostat=ios) line
    read(unit=line, fmt=*, iostat=ios) coord_type

    ! Read atom positions
    allocate(position(3,nAtoms))
    k = 0
    do i = 1, nTypes
      do j = 1, type_count(i)
        k = k + 1
        read(f_unit, fmt='(A)', iostat=ios) line
        if(ios /= 0) call abortProgram("[readElementFile] Not enough basis atoms given!")
        read(unit=line, fmt=*, iostat=ios) (position(l,k), l=1,3)
        if(coord_type == 'C' .or. coord_type == 'c' .or. coord_type == 'K' .or. coord_type == 'k') then
          ! Position of atoms given in Cartesian coordinates
          position(:,k) = position(:,k) * material%LatticeConstant
        else
          ! Position of atoms given in Bravais (or Direct, Internal, Lattice) coordinates
          position(:,k) = position(1,k) * material%a1 + position(2,k) * material%a2 + position(3,k) * material%a3
        end if
      end do
    end do

    ! Read dimension of the system
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) material%isysdim

    ! Read Fermi level
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) material%FermiLevel

    ! Read charge densitites for s p d
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) dens(1), dens(2), dens(3)
    material%OccupationS = dens(1)
    material%OccupationP = dens(2)
    material%OccupationD = dens(3)
    material%Occupation  = material%OccupationS+material%OccupationP+material%OccupationD

    ! Read Hubbard effective Coulomb interaction strength
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) material%U

    ! Read Spin-Orbit interaction strength for p and d
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) material%LambdaP, material%LambdaD

    ! Read next nearest neighbor stages
    read(f_unit, fmt='(A)', iostat = ios) line
    read(unit= line, fmt=*, iostat=ios) nn_stages

    ! Read Hopping Parameter
    do i = 1, nTypes
      ! Read on-site terms
      do j = 1, 4
        read(f_unit, fmt='(A)', iostat = ios) line
        read(unit=line, fmt=*, iostat=ios) (words(k), k=1,10)
        read(unit=words(3), fmt=*, iostat=ios) on_site(j)
      end do

      ! Setting up on-site terms
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

      ! Reading two-center integrals
      do j = 1, nn_stages
        do k = 1, 10
          read(f_unit, fmt='(A)', iostat = ios) line
          read(unit=line, fmt=*, iostat=ios) (words(l), l=1,10)
          if(j<=nStages) read(unit=words(4), fmt=*, iostat=ios) material%Hopping(k,j)
          !material%Hopping(j,i) = material%Hopping(j,i) * (a0_corr ** expon(j)) ! Correction of hopping parameter by scaling law.
        end do
      end do
    end do

    close(f_unit)

    ! Determine neighbor distances
    allocate(material%stage(nStages))
    allocate(localDistances((6*nStages+1)**3))
    m = 0
    material % stage = 10.d0 * material % LatticeConstant
    do i = -3*nStages, 3*nStages
      if( (material%isysdim==1).and.(i == 0) ) cycle ! 1D case
      do j = -3*nStages, 3*nStages
        if( (material%isysdim==2).and.(i == 0 .and. j == 0) ) cycle ! 2D case
        do k = -3*nStages, 3*nStages
          if(i == 0 .and. j == 0 .and. k == 0)  cycle
          m = m + 1
          vec = i * material%a1 + j * material%a2 + k * material%a3
          dist = vec_norm(vec,3)
          localDistances(m) = dist
          l = m - 1
          do while(1 <= l)
            if(localDistances(l) - dist < 1.d-9) exit
            localDistances(l+1) = localDistances(l)
            l = l - 1
          end do
          localDistances(l+1) = dist
        end do
      end do
    end do

    l = 1
    material%stage(1) = localDistances(1)
    do i = 1, m
      if(abs(localDistances(i) - material%stage(l)) < material%stage(1) * 0.01) cycle
      l = l + 1
      if(l > nStages) exit
      material%stage(l) = localDistances(i)
    end do
    deallocate(localDistances)
  end subroutine readElementFile

  pure subroutine intd(sss,pps,ppp,dds,ddp,ddd,sps,sds,pds,pdp,w,b)
    use mod_f90_kind,  only: double
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
  end subroutine intd

end module TightBinding
