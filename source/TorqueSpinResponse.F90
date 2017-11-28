module TorqueSpinResponse
  use mod_f90_kind, only: double
  implicit none

  complex(double), dimension(:,:,:,:), allocatable :: TSResponse, TSResponseHF
  character(len=7), dimension(2), parameter, private :: folder = ["A/Slope", "A/Slope"]
  character(len=5), dimension(2), parameter, private :: filename = ["TSR  ", "TSAHF"]
  integer, dimension(2) :: unitBase = [12000,13000]
  character(len=50), dimension(2) :: FileHeader = [ "#     energy    ,  Torque Spin Response ((i,j, i=1,3), j=1,3)", &
                                                    "#     energy    ,  Torque Spin Response ((i,j, i=1,3), j=1,3)" ]
contains
  subroutine allocTSResponse(nAtoms)
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms
    integer :: AllocateStatus

    if(allocated(TSResponse)) deallocate(TSResponse)

    allocate(TSResponse(4,4,nAtoms,nAtoms), TSResponseHF(4,4,nAtoms,nAtoms), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocateTSResponse] Not enough memory for: TSResponse, TSResponseHF")

    return
  end subroutine allocTSResponse

  subroutine create_TSR_files()
    use mod_System, only: s => sys
    use mod_BrillouinZone, only: BZ
    use mod_parameters, only: strSites, eta, Utype, suffix, fieldpart
    use mod_SOC, only: SOCc, socpart
    use EnergyIntegration, only: strEnergyParts
    implicit none
    character(len=500)  :: varm
    integer :: i,j,k

    do k = 1, size(filename)
      do i=1, s%nAtoms
        do j = 1, s%nAtoms
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(k)),trim(filename(k)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=unitBase(k)+s%nAtoms*i+j, file=varm, status='replace', form='formatted')
          write(unit=unitBase(k)+s%nAtoms*i+j, fmt="(a)") FileHeader(k)
          close(unit=unitBase(k)+s%nAtoms*i+j)
        end do
       end do
    end do

    return

  end subroutine create_TSR_files

  subroutine open_TSR_files()
    use mod_parameters, only: fieldpart, suffix, eta, Utype, strSites
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use mod_BrillouinZone, only: BZ
    use EnergyIntegration, only: strEnergyParts
    use mod_mpi_pars, only: abortProgram
    implicit none

    character(len=500) :: varm
    character(len=500) :: missing_files
    integer :: i, j, k, err, errt
    errt = 0

    do i=1, s%nAtoms
       do j = 1, s%nAtoms
         do k = 1, size(filename)
           write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(k)),trim(filename(k)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
           open (unit=unitBase(k)+s%nAtoms*i+j, file=varm, status='old', position='append', form='formatted', iostat=err)
           errt = errt + err
           if(err .ne. 0) missing_files = trim(missing_files) // " " // trim(varm)
         end do
      end do
    end do
    if(errt/=0) call abortProgram("[open_TSR_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

    return
  end subroutine open_TSR_files

  subroutine close_TSR_files()
    use mod_system, only: s => sys
    implicit none
    integer :: i,j,k

    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
         do k = 1, size(filename)
           close(unit=unitBase(k)+s%nAtoms*i+j)
         end do
       end do
    end do

    return
  end subroutine close_TSR_files


  subroutine calcTSResponse(e)
    use mod_constants, only: levi_civita, StoC, CtoS, cZero
    use mod_System, only: s => sys
    use mod_magnet, only: Lxp, Lyp, Lzp
    use mod_parameters, only: sigmai2i
    use TightBinding, only: nOrb
    use mod_parameters, only: sigmaimunu2i
    use mod_susceptibilities, only: schi, chiorb, chiorb_hf
    use mod_mpi_pars, only: abortProgram
    use mod_magnet, only: mabs
    implicit none
    complex(double), dimension(9,9,3,s%nAtoms) :: L
    integer :: i,j, m,n,k, mp,np,kp, mu,nu, gamma,zeta, p,q
    real(double), intent(in) :: e

    do i = 1, s%nAtoms
       L(1:9,1:9,1,i) = lxp(1:9,1:9,i)
       L(1:9,1:9,2,i) = lyp(1:9,1:9,i)
       L(1:9,1:9,3,i) = lzp(1:9,1:9,i)
    end do

    TSResponse = cmplx(0.d0, 0.d0)
    TSResponseHF = cmplx(0.d0, 0.d0)
    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
          do m = 1, 3
             do n = 1, 3
                do k = 1, 3
                   if(levi_civita(m,n,k) == 0.d0) cycle
                   do mp = 1,3
                      do mu = 1, nOrb
                         do nu = 1, nOrb
                            if(L(nu,mu,n,i) == cZero) cycle
                            do gamma = 1, nOrb
                               do p = 1, 4
                                  do q = 1, 4
                                     if(StoC(k+1,p) == cZero .or. CtoS(q,mp+1) == 0) cycle
                                     TSResponse(mp, m, j, i) = TSResponse(mp, m, j, i) &
                                          - 2.0d0 / sqrt(mabs(i)*mabs(j)) * s%Types(s%Basis(i)%Material)%Lambda * levi_civita(m,n,k) * L(mu, nu, n, i) &
                                          * StoC(k+1,p) * chiorb(sigmaimunu2i(p,i,mu,nu), sigmaimunu2i(q,j,gamma, gamma)) * CtoS(q,mp+1)
                                     TSResponseHF(mp, m, j, i) = TSResponse(mp, m, j, i) &
                                          - 2.0d0 / sqrt(mabs(i)*mabs(j)) * s%Types(s%Basis(i)%Material)%Lambda * levi_civita(m,n,k) * L(mu, nu, n, i) &
                                          * StoC(k+1,p) * chiorb_hf(sigmaimunu2i(p,i,mu,nu), sigmaimunu2i(q,j,gamma, gamma)) * CtoS(q,mp+1)
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    call open_TSR_files()
    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
          write(unitBase(1)+s%nAtoms*i+j, "(19(es16.9,2x))") e, ((real(TSResponse(q,p,j,i)), aimag(TSResponse(q,p,j,i)), q = 1, 3), p = 1, 3)
          write(unitBase(2)+s%nAtoms*i+j, "(19(es16.9,2x))") e, ((real(TSResponseHF(q,p,j,i)), aimag(TSResponseHF(q,p,j,i)), q = 1, 3), p = 1, 3)
       end do
    end do
    call close_TSR_files()
    return
  end subroutine calcTSResponse

end module TorqueSpinResponse
