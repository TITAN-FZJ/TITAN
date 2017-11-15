module TorqueTorqueResponse
  use mod_f90_kind, only: double
  implicit none

  complex(double), dimension(:,:,:,:), allocatable :: TTResponse, TTResponseHF
  character(len=7), dimension(2), parameter, private :: folder = ["RPA    ", "A/Slope"]
  character(len=6), dimension(3), parameter, private :: filename = ["TTR   ", "TTAinv", "TTA   "]

contains
  subroutine allocTTResponse(nAtoms)
    use mod_mpi_pars, only: abortProgram
    implicit none
    integer, intent(in) :: nAtoms
    integer :: AllocateStatus

    if(allocated(TTResponse)) deallocate(TTResponse)

    allocate(TTResponse(4,4,nAtoms,nAtoms), TTResponseHF(4,4,nAtoms,nAtoms), stat = AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[allocateTTResponse] Not enough memory for: TTResponse, TTResponseHF")

    return
  end subroutine allocTTResponse

  subroutine create_TTR_files()
    use mod_System, only: s => sys
    use mod_BrillouinZone, only: BZ
    use mod_parameters, only: strSites, eta, Utype, suffix, fieldpart
    use mod_SOC, only: SOCc, socpart
    use EnergyIntegration, only: strEnergyParts
    implicit none
    character(len=500)  :: varm
    integer :: i,j

    do i=1, s%nAtoms
       do j = 1, s%nAtoms
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(1)),trim(filename(1)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=555+s%nAtoms*i+j, file=varm, status='replace', form='formatted')
          write(unit=555+s%nAtoms*i+j, fmt="('#     energy    ,  Torque Torque Response ((i,j, i=1,3), j=1,3)')")

          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_','asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(2)),trim(filename(2)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=666+s%nAtoms*i+j, file=varm, status='replace', position='append', form='formatted')
          write(unit=666+s%nAtoms*i+j, fmt="('#     energy    ,  Torque Torque Response ((i,j, i=1,3), j=1,3)')")

          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_','asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(2)),trim(filename(3)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=777+s%nAtoms*i+j, file=varm, status='replace', position='append', form='formatted')
          write(unit=777+s%nAtoms*i+j, fmt="('#     energy    ,  Inverse Torque Torque Response ((i,j, i=1,3), j=1,3)')")

          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'HF_','asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(2)),trim(filename(2)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=888+s%nAtoms*i+j, file=varm, status='replace', position='append', form='formatted')
          write(unit=888+s%nAtoms*i+j, fmt="('#     energy    ,  Torque Torque Response ((i,j, i=1,3), j=1,3)')")


       end do
    end do

    return

  end subroutine create_TTR_files

  subroutine open_TTR_files()
    use mod_parameters, only: fieldpart, suffix, eta, Utype, strSites
    use mod_SOC, only: SOCc, socpart
    use mod_system, only: s => sys
    use mod_BrillouinZone, only: BZ
    use EnergyIntegration, only: strEnergyParts
    use mod_mpi_pars, only: abortProgram
    implicit none

    character(len=500) :: varm
    character(len=500) :: missing_files
    integer :: i, j,err, errt
    errt = 0

    do i=1, s%nAtoms
       do j = 1, s%nAtoms
          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(1)),trim(filename(1)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=555+s%nAtoms*i+j, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err .ne. 0) missing_files = trim(missing_files) // " " // trim(varm)

          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_','asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(2)),trim(filename(2)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=666+s%nAtoms*i+j, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err .ne. 0) missing_files = trim(missing_files) // " " // trim(varm)

          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'_','asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(2)),trim(filename(3)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=777+s%nAtoms*i+j, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err .ne. 0) missing_files = trim(missing_files) // " " // trim(varm)

          write(varm,"('./results/',a1,'SOC/',a,'/',a,'/',a,'HF_','asite=',i0,'_bsite=',i0,a,'_nkpt=',i0,'_eta=',es8.1,'_Utype=',i0,a,a,a,'.dat')") SOCc,trim(strSites),trim(folder(2)),trim(filename(2)),i,j,trim(strEnergyParts),BZ%nkpt,eta,Utype,trim(fieldpart),trim(socpart),trim(suffix)
          open (unit=888+s%nAtoms*i+j, file=varm, status='old', position='append', form='formatted', iostat=err)
          errt = errt + err
          if(err .ne. 0) missing_files = trim(missing_files) // " " // trim(varm)

       end do
    end do
    if(errt/=0) call abortProgram("[open_alpha_files] Some file(s) do(es) not exist! Stopping before starting calculations..." // NEW_line('A') // trim(missing_files))

    return
  end subroutine open_TTR_files

  subroutine close_TTR_files()
    use mod_system, only: s => sys
    implicit none
    integer :: i,j

    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
          close(unit=555+s%nAtoms*i+j)
          close(unit=666+s%nAtoms*i+j)
          close(unit=777+s%nAtoms*i+j)
          close(unit=888+s%nAtoms*i+j)
       end do
    end do

    return
  end subroutine close_TTR_files


  subroutine calcTTResponse(e)
    use mod_constants, only: levi_civita, StoC, CtoS, cZero
    use mod_System, only: s => sys
    use mod_magnet, only: Lxp, Lyp, Lzp
    use mod_parameters, only: sigmai2i
    use TightBinding, only: nOrb
    use mod_parameters, only: sigmaimunu2i
    use mod_susceptibilities, only: schi, chiorb, chiorb_hf
    use mod_mpi_pars, only: abortProgram
    implicit none
    complex(double), dimension(9,9,3,s%nAtoms) :: L
    integer :: AllocateStatus
    integer :: i,j, m,n,k, mp,np,kp, mu,nu, gamma,zeta, p,q
    real(double), intent(in) :: e
    complex(double), dimension(:,:), allocatable :: TTInverse
    allocate(TTInverse(s%nAtoms*3, s%nAtoms*3), stat=AllocateStatus)
    if(AllocateStatus /= 0) call abortProgram("[calcTTResponse] Couldn't allocate memory for: TTInverse")

    do i = 1, s%nAtoms
       L(1:9,1:9,1,i) = lxp(1:9,1:9,i)
       L(1:9,1:9,2,i) = lyp(1:9,1:9,i)
       L(1:9,1:9,3,i) = lzp(1:9,1:9,i)
    end do

    TTResponse = cmplx(0.d0, 0.d0)

    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
          do m = 1, 3
             do n = 1, 3
                do k = 1, 3
                   if(levi_civita(m,n,k) == 0.d0) cycle
                   do mp = 1,3
                      do np = 1, 3
                         do kp = 1, 3
                            if(levi_civita(mp,np,kp) == 0.d0) cycle
                            do mu = 1, nOrb
                               do nu = 1, nOrb
                                 if(L(nu,mu,n,i) == cZero) cycle
                                  do gamma = 1, nOrb
                                     do zeta = 1, nOrb
                                        if(L(gamma,zeta,np,j) == cZero) cycle
                                        do p = 1, 4
                                           do q = 1, 4
                                              if(StoC(k+1,p) == cZero .or. CtoS(q,kp+1) == 0) cycle
                                              TTResponse(mp, m, j, i) = TTResponse(mp, m, j, i) &
                                                   - s%Types(s%Basis(i)%Material)%Lambda * s%Types(s%Basis(j)%Material)%Lambda &
                                                   * levi_civita(m,n,k) * levi_civita(mp, np, kp) &
                                                   * L(mu, nu, n, i) * L(gamma, zeta, np, j) &
                                                   * StoC(k+1,p) * chiorb(sigmaimunu2i(p,i,mu,nu), sigmaimunu2i(q,j,gamma,zeta)) * CtoS(q,kp+1)
                                              TTResponseHF(mp, m, j, i) = TTResponseHF(mp, m, j, i) &
                                                   - s%Types(s%Basis(i)%Material)%Lambda * s%Types(s%Basis(j)%Material)%Lambda &
                                                   * levi_civita(m,n,k) * levi_civita(mp, np, kp) &
                                                   * L(mu, nu, n, i) * L(gamma, zeta, np, j) &
                                                   * StoC(k+1,p) * chiorb_hf(sigmaimunu2i(p,i,mu,nu), sigmaimunu2i(q,j,gamma,zeta)) * CtoS(q,kp+1)

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
          end do
       end do
    end do

    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
          write(unit=555+s%nAtoms*i+j, fmt="(es16.9,2x,i0,2x,i0,2x,18(es16.9,2x))") e, i,j,(((real(TTResponse(n,m,j,i)), aimag(TTResponse(n,m,j,i))), n=1,3), m = 1,3)
       end do
    end do

    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
          do m = 1,3
             do n = 1,3
                TTInverse(3*(j-1) + n, 3*(i-1)+m) = TTResponse(n,m,j,i)
             end do
          end do
       end do
    end do
    call invers(TTInverse, 3*s%nAtoms)

    do i = 1, s%nAtoms
       do j = 1, s%nAtoms
          write(666+s%nAtoms*i+j, "(19(es16.9,2x))") e, ((real(TTInverse(3*(j-1)+q, 3*(i-1)+p)), aimag(TTInverse(3*(j-1)+q, 3*(i-1)+p)), q = 1, 3), p = 1, 3)
          write(777+s%nAtoms*i+j, "(19(es16.9,2x))") e, ((real(TTResponse(q,p,j,i)), aimag(TTResponse(q,p,j,i)), q = 1, 3), p = 1, 3)
          write(888+s%nAtoms*i+j, "(19(es16.9,2x))") e, ((real(TTResponseHF(q,p,j,i)), aimag(TTResponseHF(q,p,j,i)), q = 1, 3), p = 1, 3)
       end do
    end do

    return
  end subroutine calcTTResponse

end module TorqueTorqueResponse
