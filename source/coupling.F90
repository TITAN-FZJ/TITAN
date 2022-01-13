! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
subroutine coupling()
  use mod_kind,       only: dp
  use mod_parameters, only: output,nQvec1,nQvec,deltak,kpoints
  ! use mod_parameters, only: kdirection,bsfile,wsfile
  use mod_magnet,     only: mvec_cartesian,mabs
  use mod_system,     only: s => sys
  use mod_tools,      only: vec_norm
  use mod_mpi_pars,   only: abortProgram,rField,sField,FieldComm
  use adaptiveMesh,   only: genLocalEKMesh,freeLocalEKMesh
  use mod_Coupling,   only: Jij,trJij,Jija,Jijs,allocateCoupling,deallocateCoupling,openCouplingFiles,closeCouplingFiles,writeCoupling
  use adaptiveMesh,   only: bzs
  implicit none
  integer  :: i,j,mu,kount
  logical  :: lprint = .true.
  real(dp) :: q(3)

  external :: jij_energy

  if(rField == 0) write(output%unit_loop,"('CALCULATING FULL TENSOR OF EXHANGE INTERACTIONS AND ANISOTROPIES')")

  if(sum(mabs(:))<1.e-8_dp) &
    call abortProgram("[coupling] No magnetic layers for coupling calculation!")

  if(rField == 0) call openCouplingFiles()

  call allocateCoupling()

  call genLocalEKMesh(s,rField,sField, FieldComm,bzs)

  do kount=1,nQvec1
    if(rField == 0) write(output%unit_loop,"('[coupling] ',i0,' of ',i0,' points',', i = ',es10.3)") kount,nQvec1,dble((kount-1._dp)/nQvec)
    q = kpoints(:,kount)
    call jij_energy(q,Jij)

    if(rField == 0) then

      do i=1,s%nAtoms
        do j=1,s%nAtoms
          trJij(i,j)    = 0.5_dp*(Jij(i,j,1,1)+Jij(i,j,2,2))
          Jija(i,j,:,:) = 0.5_dp*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
          Jijs(i,j,:,:) = 0.5_dp*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
          do mu = 1, 3
            Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
          end do
        end do
      end do

      ! Print only Gamma point (and only once)
      if((vec_norm(q,3)<1.e-12_dp).and.(lprint)) then
        ! Writing exchange couplings and anisotropies
        write(output%unit_loop,"('  ************************* Full tensor Jij:  *************************')")
        do i=1,s%nAtoms
          do j=1,s%nAtoms
          ! Writing on screen
          ! Writing original full tensor Jij (in units of Ry)
          ! Only the transverse components are supposed to be non-zero (e.g., for m //z, only Jxx,Jxy,Jyx,Jyy)
          ! Relation between J_ii calculated and the position of the peak in the susceptibility:
          ! w_res = 2*gamma*sqrt( (K_z-K_x)*(K_z-K_y) )/mz  - for m // z (local frame of reference)
          ! where K_x = J_ii^xx/2 ; K_y = J_ii^yy/2 ; K_z = J_ii^zz/2
          ! K > 0 - easy axis ; K < 0 - hard axis
            if(i==j) then
              write(output%unit_loop,"(3x,' ******************* Magnetization components:  ******************')")
              write(output%unit_loop,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mvec_cartesian(1,i),i,mvec_cartesian(2,i),i,mvec_cartesian(3,i)
              write(output%unit_loop,"(' |--------------- i = ',i0,'   j = ',i0,': anisotropies ---------------|')") i,j
            else
              write(output%unit_loop,"(' |----------- i = ',i0,'   j = ',i0,': exchange couplings -------------|')") i,j
            end if
            write(output%unit_loop,"('             x                  y                  z')")
            write(output%unit_loop,"('  x  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,1,1),Jij(i,j,1,2),Jij(i,j,1,3)
            write(output%unit_loop,"('  y  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,2,1),Jij(i,j,2,2),Jij(i,j,2,3)
            write(output%unit_loop,"('  z  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,3,1),Jij(i,j,3,2),Jij(i,j,3,3)
          end do
        end do
        if(s%nAtoms>1) write(output%unit_loop,"('  *** Symmetric and antisymmetric exchange interactions:  ***')")
        do i=1,s%nAtoms
          do j=1,s%nAtoms
            if(i==j) cycle
            write(output%unit_loop,"(' |--------------------- i = ',i0,'   j = ',i0,' -----------------------|')") i,j
          ! Writing Heisenberg exchange interactions
            write(output%unit_loop,"('     Isotropic:     J     = ',es16.9)") trJij(i,j)
            write(output%unit_loop,"('   Anisotropic:     Js_xx = ',es16.9)") Jijs(i,j,1,1)
            write(output%unit_loop,"('                    Js_yy = ',es16.9)") Jijs(i,j,2,2)
            write(output%unit_loop,"('  DMI: Dz = (Jxy - Jyx)/2 = ',es16.9)") Jija(i,j,1,2)
            write(output%unit_loop,"(' --------- z components of Jij  ---------')")
            write(output%unit_loop,"('  Anisotropic:  Js_zz = ',es16.9)") Jijs(i,j,3,3)
            write(output%unit_loop,"('  DMI: Dy = (Jzx - Jxz)/2 = ',es16.9)") -Jija(i,j,1,3)
            write(output%unit_loop,"('  DMI: Dx = (Jyz - Jzy)/2 = ',es16.9)") Jija(i,j,2,3)
          end do
        end do
        lprint = .false.
      end if ! Print only Gamma point (and only once)

      ! Writing into files
      ! Exchange interactions
      call writeCoupling(0._dp,dble((kount-1._dp)*deltak))
    end if

  end do

  ! Closing files
  call closeCouplingFiles()

  call freeLocalEKMesh()

  call deallocateCoupling()

end subroutine coupling

subroutine real_coupling()
  use mod_kind,          only: dp, int64, int32
  use mod_parameters,    only: output, cluster_layers, nqpt, total_nkpt => kptotal_in, total_nqpt => qptotal_in, kp_in, qp_in

  ! use mod_parameters, only: kdirection,bsfile,wsfile
  use mod_magnet,        only: mabs
  use mod_system,        only: s => sys, System_type
  use mod_tools,         only: vec_norm
  use mod_mpi_pars,      only: abortProgram,rField,sField,FieldComm,FreqComm,ierr
  use adaptiveMesh,      only: genLocalEKMesh,freeLocalEKMesh
  use mod_Coupling,      only: Jij,Jij_q,allocateCoupling,deallocateCoupling,openRealCouplingFiles,closeCouplingFiles,writeCoupling
  use adaptiveMesh,      only: bzs
  use mod_BrillouinZone, only: realBZ, q_realBZ
  use mod_constants,     only: cI
  use Lattice,           only: initLattice
  use mod_progress,      only: write_time
  implicit none
  integer  :: i,j,k, size, iw, counter
  integer(int64) :: iz
  integer(int32) :: stages
  integer  :: cell_index(3)
  real(dp) :: q(3)
  real(dp), dimension(3) :: cell_vector, atom_vector ,norms_vec
  real(dp) :: w, rx, ry, rz, rx_0, ry_0, rz_0, r_norm, sphere_radius
  complex(dp) :: kpExp
  real(dp), dimension(:,:,:,:,:), allocatable :: Jij_real
  integer :: cells
  real(dp), dimension(:,:), allocatable :: cluster

  external :: jij_energy

  if(rField == 0) write(output%unit_loop,"('CALCULATING REAL SPACE FULL TENSOR OF EXHANGE INTERACTIONS AND ANISOTROPIES')")

  if(sum(mabs(:))<1.e-8_dp) &
    call abortProgram("[coupling] No magnetic layers for coupling calculation!")

  if(rField == 0) call openRealCouplingFiles()

  ! In this part of the code we use two meshes. One is used to run over the points in which we will
  ! calculate Jij(q), this mesh will be at q_realBZ. The second mesh will be used to calculate each 
  ! Jij(q), inside this calculation there is a sum over k-points which typically needs a finer mesh
  ! this one will be stored with the standard names.
  ! To generate the q_realBZ we need to establish the number of kpoints it will consist of, this is done
  ! below
  select case(s%isysdim)
  case(3)
    q_realBZ % nkpt_x = ceiling((dble(total_nqpt))**(0.333333333333333_dp),kind(qp_in(1)) )
    q_realBZ % nkpt_y = q_realBZ % nkpt_x
    q_realBZ % nkpt_z = q_realBZ % nkpt_x
  case(2)
    q_realBZ % nkpt_x = ceiling((dble(total_nqpt))**(0.5_dp),kind(qp_in(1)) )
    q_realBZ % nkpt_y = q_realBZ % nkpt_x
    q_realBZ % nkpt_z = 1
  case default
    q_realBZ % nkpt_x = ceiling((dble(total_nqpt)), kind(qp_in(1)) )
    q_realBZ % nkpt_y = 1
    q_realBZ % nkpt_z = 1
  end select
  call q_realBZ % countBZ(s)  

  ! Generate the q mesh. All processors will have a copy of ALL the points in the BZ
  call q_realBZ % setup_fraction(s,0, 1, FreqComm(1))

  ! [To delete]
  ! call realBZ % setup_fraction(s,0, 1, FreqComm(1))
  ! write(*,*) q_realBZ % nkpt_x, q_realBZ % nkpt_y, q_realBZ % nkpt_z
  ! call MPI_Barrier(FieldComm, ierr)
  ! call q_realBZ % setup_fraction(s,0, 1, FreqComm(1))

  if(rField == 0) then
    write(*,*) "Qgrid"
    write(*,*) q_realBZ % nkpt_x, q_realBZ % nkpt_y, q_realBZ % nkpt_z
  end if
  call MPI_Barrier(FieldComm, ierr)
  ! [To delete]

  ! [To delete]
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! select case(s%isysdim)
  ! case(3)
  !   realBZ % nkpt_x = ceiling((dble(total_nkpt))**(0.333333333333333_dp),kind(kp_in(1)) )
  !   realBZ % nkpt_y = realBZ % nkpt_x
  !   realBZ % nkpt_z = realBZ % nkpt_x
  ! case(2)
  !   realBZ % nkpt_x = ceiling((dble(total_nkpt))**(0.5_dp),kind(kp_in(1)) )
  !   realBZ % nkpt_y = realBZ % nkpt_x
  !   realBZ % nkpt_z = 1
  ! case default
  !   realBZ % nkpt_x = ceiling((dble(total_nkpt)), kind(kp_in(1)) )
  !   realBZ % nkpt_y = 1
  !   realBZ % nkpt_z = 1
  ! end select
  ! call realBZ % countBZ(s)  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(rField == 0) then
    write(*,*) "Kgrid"
    ! write(*,*) realBZ % nkpt_x, realBZ % nkpt_y, realBZ % nkpt_z
  end if
  call MPI_Barrier(FieldComm, ierr)
  ! [To delete]

  ! Generating the mesh of local points on the BZ zone
  call genLocalEKMesh(s,rField,sField, FieldComm,bzs)
  ! All processors will have a copy of all the points in the BZ after we call the function below
  ! call realBZ % setup_fraction(s,0, 1, FreqComm(1))

  call allocateCoupling()

  ! [To delete] temporary prints to terminal
  if(rField == 0) then
    write(*,*) "realBZ%workload", q_realBZ%workload
    write(*,*) "nqpt = ", nqpt, " total_nqpt = ", total_nqpt, total_nkpt
    write(*,*) "Meshes generated"
  end if
  ! [To delete]

! [To delete]
  ! stop
! [To delete]
  
  if(rField == 0) then
      call write_time('[real_coupling] Started Jij(q) calculation: ',output%unit_loop)
      open(unit=13131, file='jij_q.dat', status = 'replace')
      allocate(Jij_q(q_realBZ%workload,s%nAtoms,s%nAtoms,3,3))
  end if

  do iz=1,q_realBZ%workload
    ! if(rField == 0) write(*,*) iz, " out of ", realBZ%workload
    if(rField == 0) write(output%unit,* ) iz, " out of ", q_realBZ%workload
    if(rField == 0) call write_time('[real_coupling] Step: ',output%unit_loop)
    q = q_realBZ%kp(1:3,iz)
    ! w = realBZ%w(iz)
    ! Here we call the subroutine to calculate jij at point q. Inside the function
    ! it takes the local points k and calculates their share for the total Jij(k)
    call jij_energy(q,Jij)
    if(rField == 0) then
      Jij_q(iz,:,:,:,:) = Jij
      write(13131,*) q_realBZ%w(iz), q, Jij(:,:,:,:) ! This line is targeting the 6th atom
    end if
  end do

  if(rField == 0) then
    close(13131)
    call write_time('[real_coupling] Finished Jij(q) calculation: ',output%unit_loop)
  end if

  stages = cluster_layers
  size = 0

  if(rField ==0) then
    cells = (2*stages+1)**(s%isysdim)

    ! Getting the radius of the cluster's sphere
    select case(s%isysdim)
    case(3)
        norms_vec(1) = vec_norm(s%a1*cluster_layers, 3)
        norms_vec(2) = vec_norm(s%a2*cluster_layers, 3)
        norms_vec(3) = vec_norm(s%a3*cluster_layers, 3)
        ! The small value is to account for small numerical differences
        sphere_radius = 0.0000000001_dp + MINVAL(norms_vec)
        write(*,*) "Radius = ", MINVAL(norms_vec), s%a1*cluster_layers
    case(2)
        norms_vec(1) = vec_norm(s%a1*cluster_layers, 3)
        norms_vec(2) = vec_norm(s%a2*cluster_layers, 3)
        ! The small value is to account for small numerical differences
        sphere_radius = 0.0000000001_dp + MINVAL(norms_vec(1:2))
        write(*,*) "Radius = ", sphere_radius, s%a1*cluster_layers, s%a2*cluster_layers
    case default
        ! The small value is to account for small numerical differences
        sphere_radius = 0.0000000001_dp + vec_norm(s%a1*cluster_layers, 3)
    end select

    ! Checking which atoms are inside the sphere
    counter = 0
    do i = 1, cells
        ! ...and atoms in the unit cell
        ! do j = 1, s%nAtoms
        ! "size" is the current atom
        size = size + 1

        select case(s%isysdim)
        case(3)
            cell_index(1) = mod( (i-1),(2*stages+1) ) - stages
            cell_index(2) = mod( (i-1)/(2*stages+1),(2*stages+1) ) - stages
            cell_index(3) = mod( (i-1)/((2*stages+1)*(2*stages+1)),(2*stages+1) ) - stages
        case(2)
            cell_index(1) = mod( (i-1),(2*stages+1) ) - stages
            cell_index(2) = mod( (i-1)/(2*stages+1),(2*stages+1) ) - stages
            cell_index(3) = 0
        case default
            cell_index(1) = mod( (i-1),(2*stages+1) ) - stages
            cell_index(2) = 0
            cell_index(3) = 0
        end select

        cell_vector = cell_index(1) * s%a1 + cell_index(2) * s%a2 + cell_index(3) * s%a3
        ! Atom position is r = R_i + r_j
        ! atom_vector = s%Basis(j)%Position + cell_vector

        ! write(*,*) cluster(size,:), "norm = ", vec_norm(atom_vector, 3), "Rad = ", sphere_radius, "Bool = ", vec_norm(atom_vector, 3) <= sphere_radius
        if(vec_norm(cell_vector, 3) <= sphere_radius) then
            counter = counter + 1
            write(*,*) cell_vector
        end if

        ! [To delete]
        ! write(*,*) "cell_vector", cell_vector
        ! write(*,*) "cluster(size,:)", cluster(size,:)
        ! write(*,*) cluster(size,:), "norm = ", vec_norm(atom_vector, 3)
        ! write(*,*) " "
        ! [To delete]

        ! end do
    end do

    !TO DELETE
    !write(*,*) "cluster", counter
    !
    if(allocated(cluster)) deallocate(cluster)
    allocate(cluster(counter,3))
    cluster = 0._dp

    ! Pupulating cluster array
    counter = 0

    do i = 1, cells
        ! ...and atoms in the unit cell
        ! do j = 1, s%nAtoms
        ! "size" is the current atom
        size = size + 1

        select case(s%isysdim)
        case(3)
            cell_index(1) = mod( (i-1),(2*stages+1) ) - stages
            cell_index(2) = mod( (i-1)/(2*stages+1),(2*stages+1) ) - stages
            cell_index(3) = mod( (i-1)/((2*stages+1)*(2*stages+1)),(2*stages+1) ) - stages
        case(2)
            cell_index(1) = mod( (i-1),(2*stages+1) ) - stages
            cell_index(2) = mod( (i-1)/(2*stages+1),(2*stages+1) ) - stages
            cell_index(3) = 0
        case default
            cell_index(1) = mod( (i-1),(2*stages+1) ) - stages
            cell_index(2) = 0
            cell_index(3) = 0
        end select

        cell_vector = cell_index(1) * s%a1 + cell_index(2) * s%a2 + cell_index(3) * s%a3
        ! Atom position is r = R_i + r_j
        ! atom_vector = s%Basis(j)%Position + cell_vector

        ! write(*,*) cluster(size,:), "norm = ", vec_norm(atom_vector, 3), "Rad = ", sphere_radius, "Bool = ", vec_norm(atom_vector, 3) <= sphere_radius
        if(vec_norm(cell_vector, 3) <= sphere_radius) then
            counter = counter + 1
            cluster(counter,:) = cell_vector
            ! cluster(counter,:) = atom_vector
            ! write(*,*) cluster(counter,:)
        end if

        ! write(*,*) " "

        ! end do
    end do

    call write_time('[real_coupling] Finished cluster generation: ',output%unit_loop)

! ##############################################################################################################################


    if(allocated(Jij_real)) deallocate(Jij_real)
    allocate(Jij_real(counter,s%nAtoms,s%nAtoms,3,3))

    Jij_real = 0._dp

    ! [To delete]
    ! write(*,*) "Positions"
    ! [To delete]
    ! pragma omp for
    do k = 1, counter
        do iz=1,q_realBZ%workload
            kpExp = exp(-1._dp * cI * dot_product(q_realBZ%kp(1:3,iz), cluster(k,:)))
            w = q_realBZ%w(iz)
            Jij_real(k,:,:,:,:) = Jij_real(k,:,:,:,:) + real(kpExp*w*Jij_q(iz,:,:,:,:))
        end do

        ! [To delete]
        ! This was a trial when we wanted to see if the problem was on the fourier transform
        ! do i=1,s%nAtoms
        !   do j=1,s%nAtoms
        !     do iz=1,q_realBZ%workload
        !       kpExp = exp(-1._dp * cI * dot_product(q_realBZ%kp(1:3,iz), cluster(k,:)+s%Basis(j)%Position -  s%Basis(i)%Position))
        !       w = q_realBZ%w(iz)
        !       Jij_real(k,i,j,:,:) = Jij_real(k,i,j,:,:) + real(kpExp*w*Jij_q(iz,i,j,:,:))
        !     end do
        !   end do
        ! end do
        ! [To delete]

        ! [To delete]
        ! write(*,*) "s%Neighbors(k)%CellVector", cluster(k,:)
        ! write(*,*) " "
        ! write(*,*) "Jij_real(r,:,:,:,:)", Jij_real(k,:,:,:,:)
        ! write(*,*) " "
        ! [To delete]
    end do

    !!!!!! TO DELETE
    !write(*,*) "r0 = ", s%Basis(12)%Position(1), s%Basis(12)%Position(2) , s%Basis(12)%Position(3)
    !!!!!!!!!!!!!!!!
    ! Writing files
    do k = 1, counter
        do j=1,s%nAtoms
            rx_0 = s%Basis(j)%Position(1)
            ry_0 = s%Basis(j)%Position(2)
            rz_0 = s%Basis(j)%Position(3)
            do i=1,s%nAtoms
                iw = 2000 + (j-1) * s%nAtoms * 2 + (i-1) * 2
                if(i==j) then
                    iw = iw + 1
                    !!!!!! TO DELETE
                    ! write(*,*) "i, j, iw ", i , j, iw
                    !!!!!! 
                    rx = cluster(k,1) + s%Basis(i)%Position(1)
                    ry = cluster(k,2) + s%Basis(i)%Position(2)
                    rz = cluster(k,3) + s%Basis(i)%Position(3)
                    r_norm = SQRT((rx-rx_0)*(rx-rx_0) + (ry-ry_0)*(ry-ry_0) + (rz-rz_0)*(rz-rz_0))
                    write(unit=iw,fmt="(13(es16.9,2x))") rx, ry, rz ,r_norm,&
                    Jij_real(k,i,j,1,1),Jij_real(k,i,j,1,2), Jij_real(k,i,j,1,3), &
                    Jij_real(k,i,j,2,1), Jij_real(k,i,j,2,2), Jij_real(k,i,j,2,3), &
                    Jij_real(k,i,j,3,1), Jij_real(k,i,j,3,1), Jij_real(k,i,j,3,1)
                else
                    iw = iw + 1
                    !!!!!! TO DELETE
                    ! write(*,*) "i, j, iw ", i , j, iw
                    !!!!!!
                    !!!!!!!!!! Testing this block
                    rx = cluster(k,1) + s%Basis(i)%Position(1)
                    ry = cluster(k,2) + s%Basis(i)%Position(2)
                    rz = cluster(k,3) + s%Basis(i)%Position(3)
                    r_norm = SQRT((rx-rx_0)*(rx-rx_0) + (ry-ry_0)*(ry-ry_0) + (rz-rz_0)*(rz-rz_0))
                    !!!!!!!!!!!
                    write(unit=iw,fmt="(13(es16.9,2x))") rx, ry, rz ,r_norm,&
                    Jij_real(k,i,j,1,1),Jij_real(k,i,j,1,2), Jij_real(k,i,j,1,3), &
                    Jij_real(k,i,j,2,1), Jij_real(k,i,j,2,2), Jij_real(k,i,j,2,3), &
                    Jij_real(k,i,j,3,1), Jij_real(k,i,j,3,1), Jij_real(k,i,j,3,1)
                    ! iw = iw + 1
                    ! write(unit=iw,fmt="(2(es16.9,2x))") q,Jija(i,j,1,2)
                end if
            end do
        end do
    end do

    call write_time('[real_coupling] Finished Fourier transform: ',output%unit_loop)

  end if

  ! Closing files
  call closeCouplingFiles()

  call freeLocalEKMesh()

  call deallocateCoupling()

  if(allocated(Jij_real)) deallocate(Jij_real)
  if(allocated(cluster)) deallocate(cluster)

end subroutine real_coupling
