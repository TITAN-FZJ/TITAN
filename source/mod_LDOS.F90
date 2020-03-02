module mod_LDOS
  use mod_f90_kind, only: double
  implicit none
  character(len=4), private :: folder = "LDOS"
  character(len=6), dimension(4), private :: filename = ["ldosu ","ldosd ","lmdosu", "lmdosd"]

  real(double),dimension(:,:),allocatable :: ldosu,ldosd

contains
  subroutine allocateLDOS()
    use mod_System,     only: s => sys
    use mod_parameters, only: nOrb
    use mod_superconductivity, only: superCond
    implicit none

    if(allocated(ldosu)) deallocate(ldosu)
    if(allocated(ldosd)) deallocate(ldosd)
    allocate(ldosu(s%nAtoms,nOrb*superCond))
    allocate(ldosd(s%nAtoms,nOrb*superCond))

  end subroutine allocateLDOS

  subroutine deallocateLDOS()
    implicit none
    if(allocated(ldosu)) deallocate(ldosu)
    if(allocated(ldosd)) deallocate(ldosd)
  end subroutine deallocateLDOS


  ! This subroutine calculates LDOS
  subroutine ldos()
    use mod_f90_kind,      only: double
    use mod_parameters,    only: nOrb, output, nEner1, emin, deltae,laddresults
    use mod_system,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    use mod_superconductivity, only: superCond
    use mod_mpi_pars
    implicit none
    integer :: i, j
    real(double) :: e

    call allocateLDOS()
    call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))

    ! Opening files
    if(rField == 0) then
      write(output%unit_loop,"('CALCULATING LDOS')")
      if(.not.laddresults) call createLDOSFiles()
      call openLDOSFiles()
    end if

    do i = startFreq, endFreq
      e = emin + (i-1)*deltae
      if(rFreq(1) == 0) write(output%unit_loop,"('[ldos] ',i0,' of ',i0,' points',', e = ',es10.3)") i,nEner1,e
      call ldos_energy(e,ldosu,ldosd)

      ! These IFs may be wrong!
      if(rFreq(1) == 0) then
        if(rFreq(2) == 0) then
          do j = 1, sFreq(2)
            if (j /= 1) then
              call MPI_Recv(e,     1            ,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE  ,1000,FreqComm(2),stat,ierr)
              call MPI_Recv(ldosd, s%nAtoms*nOrb*superCond,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1100,FreqComm(2),stat,ierr)
              call MPI_Recv(ldosu, s%nAtoms*nOrb*superCond,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1200,FreqComm(2),stat,ierr)
            end if

            ! Writing into files
            call writeLDOS(e)
          end do
        else
          call MPI_Send(e,     1            ,MPI_DOUBLE_PRECISION,0,1000,FreqComm(2),stat,ierr)
          call MPI_Send(ldosd, s%nAtoms*nOrb*superCond,MPI_DOUBLE_PRECISION,0,1100,FreqComm(2),stat,ierr)
          call MPI_Send(ldosu, s%nAtoms*nOrb*superCond,MPI_DOUBLE_PRECISION,0,1200,FreqComm(2),stat,ierr)
        end if
      end if
      call MPI_Barrier(FieldComm, ierr)
     end do

    call deallocateLDOS()

    ! Closing files
    if(rField == 0) call closeLDOSFiles()

  end subroutine ldos


  ! Calculates spin-resolved LDOS and energy-dependence of exchange interactions
  subroutine ldos_energy(e,ldosu,ldosd)
    use mod_f90_kind,      only: double
    use mod_constants,     only: pi
    use mod_parameters,    only: nOrb, nOrb2, eta
    use mod_system,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    use mod_superconductivity, only: lsupercond, green_sc, superCond
    use mod_mpi_pars
    implicit none
    real(double), intent(in) :: e
    real(double), dimension(s%nAtoms, nOrb*superCond), intent(out) :: ldosu, ldosd
    complex(double), dimension(nOrb2*superCond, nOrb2*superCond, s%nAtoms, s%nAtoms) :: gf
    complex(double), dimension(s%nAtoms, nOrb*superCond) :: gfdiagu,gfdiagd
    real(double), dimension(3) :: kp
    real(double) :: weight
    integer :: i,mu,nu
    integer*8 :: iz

    ldosu = 0.d0
    ldosd = 0.d0

    !$omp parallel default(none) &
    !$omp& private(iz,kp,weight,gf,i,mu,nu,gfdiagu,gfdiagd) &
    !$omp& shared(s,realBZ,e,nOrb,eta,ldosu,ldosd,lsupercond,nOrb2)
    !$omp do reduction(+:ldosu,ldosd)
    do iz = 1,realBZ%workload
      kp = realBZ%kp(1:3,iz)
      weight = realBZ%w(iz)
      ! Green function on energy E + ieta, and wave vector kp
      if(lsupercond == .true.) then
          call green_sc(e,eta,s,kp,gf)
      else
          call green(e,eta,s,kp,gf)
      end if

      ! Density of states
      do mu=1,nOrb
        do i=1,s%nAtoms
          nu = mu + nOrb
          gfdiagu(i,mu) = - aimag(gf(mu,mu,i,i)) * weight
          gfdiagd(i,mu) = - aimag(gf(nu,nu,i,i)) * weight
          if(lsupercond) then
              gfdiagu(i,mu+nOrb) = - aimag(gf(mu+nOrb2,mu+nOrb2,i,i)) * weight
              gfdiagd(i,mu+nOrb) = - aimag(gf(nu+nOrb2,nu+nOrb2,i,i)) * weight
          end if
        end do
     end do

      ldosu = ldosu + gfdiagu
      ldosd = ldosd + gfdiagd

    end do
    !$omp end do
    !$omp end parallel

    ldosu  = ldosu/pi
    ldosd  = ldosd/pi

    if(rFreq(1) == 0) then
      call MPI_Reduce(MPI_IN_PLACE, ldosu , s%nAtoms*nOrb*superCond, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
      call MPI_Reduce(MPI_IN_PLACE, ldosd , s%nAtoms*nOrb*superCond, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    else
      call MPI_Reduce(ldosu , ldosu , s%nAtoms*nOrb*superCond, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
      call MPI_Reduce(ldosd , ldosd , s%nAtoms*nOrb*superCond, MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    end if

  end subroutine ldos_energy


  ! This subroutine calculates LDOS and coupling as a function of energy
  subroutine ldos_and_coupling()
    use mod_f90_kind,      only: double
    use mod_parameters,    only: nOrb, output, emin, deltae, nEner1, skip_steps
    use mod_system,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    use mod_Coupling
    use mod_mpi_pars
    implicit none
    integer :: i, j, mu, count, ncount,ncount2,ncount3
    real(double) :: e

    ncount  = s%nAtoms*nOrb
    ncount2 = s%nAtoms*s%nAtoms
    ncount3 = ncount*9

    if(rField == 0) write(output%unit_loop,"('CALCULATING LDOS AND EXCHANGE INTERACTIONS AS A FUNCTION OF ENERGY')")

    call realBZ % setup_fraction(s,rFreq(1), sFreq(1), FreqComm(1))

    ! Opening files
    if(rField == 0) then
      ! LDOS
      call openLDOSFiles()
      ! Exchange interactions
      call openCouplingFiles()
    end if

    call allocateLDOS()
    call allocateCoupling()

    do count = startFreq + skip_steps, endFreq + skip_steps
      e = emin + (count-1) * deltae
      if(rFreq(1) == 0) write(output%unit_loop,"('[ldos_and_coupling] ',i0,' of ',i0,' points',', e = ',es10.3)") count,nEner1,e

      call ldos_jij_energy(e,ldosu,ldosd,Jij)

      if(rFreq(1) == 0) then
        do i = 1, s%nAtoms
          do j = 1, s%nAtoms
            trJij(i,j)    = 0.5d0*(Jij(i,j,1,1) + Jij(i,j,2,2))
            Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
            Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
            do mu = 1, 3
              Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
            end do
          end do
        end do


        if(rFreq(2) == 0) then
          do i = 1, sFreq(2)
            if (i /= 1) then
              call MPI_Recv(e,     1       ,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE  ,1000,FreqComm(2),stat,ierr)
              call MPI_Recv(ldosd, ncount  ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1100,FreqComm(2),stat,ierr)
              call MPI_Recv(ldosu, ncount  ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1200,FreqComm(2),stat,ierr)
              call MPI_Recv(trJij, ncount2 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1300,FreqComm(2),stat,ierr)
              call MPI_Recv(Jij,   ncount3 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1400,FreqComm(2),stat,ierr)
              call MPI_Recv(Jijs,  ncount3 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1500,FreqComm(2),stat,ierr)
              call MPI_Recv(Jija,  ncount3 ,MPI_DOUBLE_PRECISION,stat(MPI_SOURCE),1600,FreqComm(2),stat,ierr)
            end if

            ! Writing into files
            call writeLDOS(e)
            ! Exchange interactions
            call writeCoupling(e)
          end do
        else
          call MPI_Recv(e,     1       ,MPI_DOUBLE_PRECISION,0,1000,FreqComm(2),stat,ierr)
          call MPI_Recv(ldosd, ncount  ,MPI_DOUBLE_PRECISION,0,1100,FreqComm(2),stat,ierr)
          call MPI_Recv(ldosu, ncount  ,MPI_DOUBLE_PRECISION,0,1200,FreqComm(2),stat,ierr)
          call MPI_Recv(trJij, ncount2 ,MPI_DOUBLE_PRECISION,0,1300,FreqComm(2),stat,ierr)
          call MPI_Recv(Jij,   ncount3 ,MPI_DOUBLE_PRECISION,0,1400,FreqComm(2),stat,ierr)
          call MPI_Recv(Jijs,  ncount3 ,MPI_DOUBLE_PRECISION,0,1500,FreqComm(2),stat,ierr)
          call MPI_Recv(Jija,  ncount3 ,MPI_DOUBLE_PRECISION,0,1600,FreqComm(2),stat,ierr)
        end if
      end if
      call MPI_Barrier(FieldComm, ierr)

    end do

    call deallocateLDOS()
    call deallocateCoupling()

    ! Closing files
    if(rField == 0) then
      call closeLDOSFiles()
      call closeCouplingFiles()
    end if

  end subroutine ldos_and_coupling


  ! Calculates spin-resolved LDOS and energy-dependence of exchange interactions
  subroutine ldos_jij_energy(e,ldosu,ldosd,Jijint)
    use mod_f90_kind,      only: double
    use mod_constants,     only: pi, cZero, cOne, pauli_dorb
    use mod_parameters,    only: nOrb, nOrb2, eta, U
    use mod_system,        only: s => sys
    use mod_BrillouinZone, only: realBZ
    use mod_magnet,        only: mvec_cartesian, mabs
    use mod_progress
    use mod_mpi_pars
    implicit none
    real(double),intent(in)     :: e
    real(double),intent(out)    :: ldosu(s%nAtoms, nOrb),ldosd(s%nAtoms, nOrb)
    real(double),intent(out)    :: Jijint(s%nAtoms,s%nAtoms,3,3)
    complex(double), dimension(nOrb2, nOrb2, s%nAtoms, s%nAtoms) :: gf
    complex(double), dimension(s%nAtoms, nOrb)   :: gfdiagu,gfdiagd
    complex(double), dimension(nOrb2, nOrb2)     :: gij,gji,temp1,temp2,paulia,paulib
    real(double),    dimension(:,:),allocatable     :: ldosu_loc,ldosd_loc
    real(double),    dimension(:,:,:,:),allocatable :: Jijint_loc
    real(double),    dimension(3) :: kp
    complex(double) :: paulimatan(3,3,nOrb2, nOrb2)
    real(double)    :: Jijkan(s%nAtoms,3,3), Jijk(s%nAtoms,s%nAtoms,3,3)
    real(double)    :: weight
    integer         :: i,j,mu,nu,alpha,ncount,ncount2
    integer*8       :: iz

    real(double),    dimension(3,s%nAtoms)               :: evec
    complex(double), dimension(s%nAtoms,3,nOrb2,nOrb2)   :: dbxcdm
    complex(double), dimension(s%nAtoms,3,3,nOrb2,nOrb2) :: d2bxcdm2
    complex(double), dimension(s%nAtoms,nOrb2,nOrb2)     :: paulievec

    ncount  = s%nAtoms*nOrb
    ncount2 = s%nAtoms*s%nAtoms*9

  ! (x,y,z)-tensor formed by Pauli matrices to calculate anisotropy term (when i=j)
    paulimatan = cZero
    paulimatan(1,1,:,:) = -pauli_dorb(3,:,:)
    paulimatan(2,2,:,:) = -pauli_dorb(3,:,:)
    paulimatan(1,3,:,:) = -pauli_dorb(1,:,:)
    paulimatan(3,1,:,:) = -pauli_dorb(1,:,:)
    paulimatan(2,3,:,:) = -pauli_dorb(2,:,:)
    paulimatan(3,2,:,:) = -pauli_dorb(2,:,:)

    ldosu = 0.d0
    ldosd = 0.d0
    Jijint = 0.d0

    do iz = 1, s%nAtoms
      ! Unit vector along the direction of the magnetization of each magnetic plane
      evec(:,iz) = [ mvec_cartesian(1,iz), mvec_cartesian(2,iz), mvec_cartesian(3,iz) ]/mabs(iz)

      ! Inner product of pauli matrix in spin and orbital space and unit vector evec
      paulievec(iz,:,:) = pauli_dorb(1,:,:) * evec(1,iz) + pauli_dorb(2,:,:) * evec(2,iz) + pauli_dorb(3,:,:) * evec(3,iz)

      do i = 1, 3
        ! Derivative of Bxc*sigma*evec w.r.t. m_i (Bxc = -U.m/2)
        dbxcdm(iz,i,:,:) = -0.5d0 * U(iz) * (pauli_dorb(i,:,:) - (paulievec(iz,:,:)) * evec(i,iz))

        ! Second derivative of Bxc w.r.t. m_i (Bxc = -U.m/2)
        do j=1,3
          d2bxcdm2(iz,i,j,:,:) = evec(i,iz)*pauli_dorb(j,:,:) + pauli_dorb(i,:,:)*evec(j,iz) - 3*paulievec(iz,:,:)*evec(i,iz)*evec(j,iz)
          if(i==j) d2bxcdm2(iz,i,j,:,:) = d2bxcdm2(iz,i,j,:,:) + paulievec(iz,:,:)
          d2bxcdm2(iz,i,j,:,:) = 0.5d0*U(iz)*d2bxcdm2(iz,i,j,:,:)/(mabs(iz))
        end do
      end do
    end do



    !$omp parallel default(none) &
    !$omp& private(iz,kp,weight,gf,gij,gji,paulia,paulib,i,j,mu,nu,alpha,gfdiagu,gfdiagd,Jijk,Jijkan,temp1,temp2,ldosu_loc,ldosd_loc,Jijint_loc) &
    !$omp& shared(s,realBZ,nOrb,nOrb2,e,eta,U,dbxcdm,d2bxcdm2,pauli_dorb,paulimatan,ldosu,ldosd,Jijint)
    allocate(ldosu_loc(s%nAtoms, nOrb), ldosd_loc(s%nAtoms, nOrb), Jijint_loc(s%nAtoms,s%nAtoms,3,3))
    ldosu_loc = 0.d0
    ldosd_loc = 0.d0
    Jijint_loc = 0.d0

    !$omp do schedule(static)
    do iz = 1, realBZ%workload
      kp = realBZ%kp(:,iz)
      weight = realBZ%w(iz)
      ! Green function on energy E + ieta, and wave vector kp
      call green(e,eta,s,kp,gf)

      ! Exchange interaction tensor
      Jijk   = 0.d0
      Jijkan = 0.d0
      do j = 1,s%nAtoms
        do i = 1,s%nAtoms
          gij = gf(:,:,i,j)
          gji = gf(:,:,j,i)
          do nu = 1,3
            do mu = 1,3
              paulia = dbxcdm(i,mu,:,:)
              paulib = dbxcdm(j,nu,:,:)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,paulia,nOrb2,gij,   nOrb2,cZero,temp1,nOrb2)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp1, nOrb2,paulib,nOrb2,cZero,temp2,nOrb2)
              call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,temp2, nOrb2,gji,   nOrb2,cZero,temp1,nOrb2)
              ! Trace over orbitals and spins
              do alpha = 1,nOrb2
                Jijk(i,j,mu,nu) = Jijk(i,j,mu,nu) + real(temp1(alpha,alpha))
              end do
            end do
          end do

          ! Anisotropy (on-site) term
          if(i==j) then
            gij = gf(:,:,i,i)
            do nu = 1,3
              do mu = 1,3
                paulia = d2bxcdm2(i,mu,nu,:,:)
                call zgemm('n','n',nOrb2,nOrb2,nOrb2,cOne,gij,nOrb2,paulia,nOrb2,cZero,temp1,nOrb2)
                ! Trace over orbitals and spins
                do alpha = 1,nOrb2
                  Jijkan(i,mu,nu) = Jijkan(i,mu,nu) + real(temp1(alpha,alpha))
                end do
                Jijk(i,i,mu,nu) = Jijk(i,i,mu,nu) + Jijkan(i,mu,nu)
              end do
            end do
          end if
        end do
      end do
      Jijint = Jijint + Jijk * weight

      ! Density of states
      do mu=1,nOrb
        do i=1,s%nAtoms
           nu=mu+nOrb
           gfdiagu(i,mu) = - aimag(gf(mu,mu,i,i))*weight
           gfdiagd(i,mu) = - aimag(gf(nu,nu,i,i))*weight
         end do
      end do

      ldosu_loc = ldosu_loc + gfdiagu
      ldosd_loc = ldosd_loc + gfdiagd
      Jijint_loc = Jijint_loc + Jijk

    end do
    !$omp end do nowait

    !$omp critical
      ldosu = ldosu + ldosu_loc
      ldosd = ldosd + ldosd_loc
      Jijint = Jijint + Jijint_loc
    !$omp end critical
    !$omp end parallel

    ldosu  = ldosu/pi
    ldosd  = ldosd/pi
    Jijint = Jijint/pi

    if(rFreq(1) == 0) then
       call MPI_Reduce(MPI_IN_PLACE, ldosu , ncount  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
       call MPI_Reduce(MPI_IN_PLACE, ldosd , ncount  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
       call MPI_Reduce(MPI_IN_PLACE, Jijint, ncount2 , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    else
       call MPI_Reduce(ldosu , ldosu , ncount  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
       call MPI_Reduce(ldosd , ldosd , ncount  , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
       call MPI_Reduce(Jijint, Jijint, ncount2 , MPI_DOUBLE_PRECISION, MPI_SUM, 0, FreqComm(1), ierr)
    end if
  end subroutine ldos_jij_energy


  subroutine createLDOSFiles()
    use mod_parameters,        only: output
    use mod_System,            only: s => sys
    use mod_io,                only: write_header
    use mod_superconductivity, only: lsuperCond
    implicit none
    character(len=400) :: varm,title(size(filename))
    integer            :: i, iw, j

    title(1:2) = "#    energy      ,  LDOS SUM       ,  LDOS S         ,  LDOS P         ,  LDOS D         "
    title(3:4) = "#    energy      ,  LDOS S         ,  LDOS PX        ,  LDOS PY        ,  LDOS PZ        ,  LDOS DXY       ,  LDOS DYZ       ,  LDOS DZX       ,  LDOS DX2       ,  LDOS DZ2       "

    do i = 1, s%nAtoms
      do j = 1, size(filename)
        iw = 1000 + (i-1) * size(filename) + j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_site=',i0,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(j)),i,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open (unit=iw, file=varm,status='replace', form='formatted')
        !                             if lsuperCond is true then it writes 0.0, otherwise it writes s%Ef
        call write_header(iw,title(j),merge(0.0,s%Ef,lsuperCond))
        ! write(unit=iw, fmt="('#   energy      ,  LDOS SUM        ,  LDOS S          ,  LDOS P          ,  LDOS T2G        ,  LDOS EG         ')")
        close(unit=iw)
      end do
    end do
  end subroutine createLDOSFiles

  subroutine openLDOSFiles()
    use mod_parameters, only: output,missing_files
    use mod_System,     only: s => sys
    use mod_mpi_pars,   only: abortProgram
    implicit none
    character(len=400) :: varm
    integer            :: i, iw, j, err, errt=0

    do i = 1, s%nAtoms
      do j = 1, size(filename)
        iw = 1000 + (i-1) * size(filename) + j
        write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_site=',i0,a,a,a,a,'.dat')") output%SOCchar,trim(output%Sites),trim(folder),trim(filename(j)),i,trim(output%info),trim(output%BField),trim(output%SOC),trim(output%suffix)
        open (unit=iw, file=varm, status='old', position='append', form='formatted', iostat=err)
        errt = errt + err
        if(err/=0) missing_files = trim(missing_files) // " " // trim(varm)
      end do
    end do
    ! Stop if some file does not exist
    if(errt/=0) call abortProgram("[openLDOSFiles] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

  end subroutine openLDOSFiles

  subroutine closeLDOSFiles()
    use mod_System, only: s => sys
    implicit none
    integer i,j, iw
    do i = 1, s%nAtoms
      do j = 1, size(filename)
        iw = 1000 + (i-1) * size(filename) + j
        close(iw)
      end do
    end do
  end subroutine closeLDOSFiles

  subroutine writeLDOS(e)
    use mod_f90_kind, only: double
    use mod_System,   only: s => sys
    use mod_superconductivity, only: lsuperCond
    implicit none
    real(double), intent(in) :: e
    integer :: i, iw, j

    if(lsuperCond) then
        do i = 1, s%nAtoms
           iw = 1000 + (i-1) * size(filename) + 1
           write(unit=iw,fmt="( 8(es16.9,2x))") e, sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9)), ldosu(i,10),sum(ldosu(i,11:13)),sum(ldosu(i,14:18))
           iw = iw + 1
           write(unit=iw,fmt="( 8(es16.9,2x))") e, sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9)), ldosd(i,10),sum(ldosd(i,11:13)),sum(ldosd(i,14:18))
           iw = iw + 1
           write(unit=iw,fmt="(19(es16.9,2x))") e, (ldosu(i,j), j = 1,18) !sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="(19(es16.9,2x))") e, (ldosd(i,j), j = 1,18) !sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
        end do
    else
        do i = 1, s%nAtoms
           iw = 1000 + (i-1) * size(filename) + 1
           write(unit=iw,fmt="( 5(es16.9,2x))") e, sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="( 5(es16.9,2x))") e, sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="(10(es16.9,2x))") e, (ldosu(i,j), j = 1,9) !sum(ldosu(i,:)),ldosu(i,1),sum(ldosu(i,2:4)),sum(ldosu(i,5:9))
           iw = iw + 1
           write(unit=iw,fmt="(10(es16.9,2x))") e, (ldosd(i,j), j = 1,9) !sum(ldosd(i,:)),ldosd(i,1),sum(ldosd(i,2:4)),sum(ldosd(i,5:9))
        end do
    end if

  end subroutine writeLDOS

  ! This subroutine sorts LDOS files
  subroutine sortLDOS()
    use mod_f90_kind, only: double
    use mod_tools,    only: sort_file
    use mod_system,   only: s => sys
    implicit none
    integer :: i,j,iw

    ! Opening LDOS files
    call openLDOSFiles()

    do i = 1, s%nAtoms
      do j = 1, size(filename)
        iw = 1000 + (i-1) * size(filename) + j
        call sort_file(iw)
      end do
    end do

    ! Closing LDOS files
    call closeLDOSFiles()

  end subroutine sortLDOS

end module mod_LDOS
