! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
subroutine coupling()
  use mod_f90_kind,   only: double
  use mod_parameters, only: output, mmlayermag, nmaglayers, q
  use mod_magnet,     only: mvec_cartesian
  use mod_system,     only: s => sys
  use adaptiveMesh
  use mod_mpi_pars
  use mod_Coupling
  implicit none
  integer            :: i,j,mu

  if(rField == 0) write(output%unit_loop,"('CALCULATING FULL TENSOR OF EXHANGE INTERACTIONS AND ANISOTROPIES')")

  if(rField == 0) call openCouplingFiles()

  call allocateCoupling()

  call genLocalEKMesh(s,rField,sField, FieldComm)
  q = [ 0.d0 , 0.d0 , 0.d0 ]
  call jij_energy(Jij)
  call freeLocalEKMesh()

  if(rField == 0) then
    do i=1,nmaglayers
      do j=1,nmaglayers
         trJij(i,j)    = 0.5d0*(Jij(i,j,1,1)+Jij(i,j,2,2))
         Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
         Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
         do mu = 1, 3
            Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
         end do
      end do
    end do

    ! Writing exchange couplings and anisotropies
    write(output%unit_loop,"('  ************************* Full tensor Jij:  *************************')")
    do i=1,nmaglayers
      do j=1,nmaglayers
      ! Writing on screen
      ! Writing original full tensor Jij
      ! Only the transverse components are supposed to be non-cZero (e.g., for m //z, only Jxx,Jxy,Jyx,Jyy)
      ! Relation between J_ii calculated and the position of the peak in the susceptibility:
      ! w_res = 2*gamma*mz*sqrt( (K_z-K_x)*(K_z-K_y) )  - for m // z
      ! where K_x = J_ii^xx/2 ; K_y = J_ii^yy/2 ; K_z = J_ii^zz/2
      ! K > 0 - easy axis ; K < 0 - hard axis
        if(i==j) then
          write(output%unit_loop,"(3x,' ******************* Magnetization components:  ******************')")
          write(output%unit_loop,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mvec_cartesian(1,i),i,mvec_cartesian(2,i),i,mvec_cartesian(3,i)
          write(output%unit_loop,"(' |--------------- i = ',i0,'   j = ',i0,': anisotropies ---------------|')") mmlayermag(i),mmlayermag(j)
        else
          write(output%unit_loop,"(' |----------- i = ',i0,'   j = ',i0,': exchange couplings -------------|')") mmlayermag(i),mmlayermag(j)
        end if
        write(output%unit_loop,"('             x                  y                  z')")
        write(output%unit_loop,"('  x  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,1,1),Jij(i,j,1,2),Jij(i,j,1,3)
        write(output%unit_loop,"('  y  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,2,1),Jij(i,j,2,2),Jij(i,j,2,3)
        write(output%unit_loop,"('  z  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,3,1),Jij(i,j,3,2),Jij(i,j,3,3)
      end do
    end do
    if(nmaglayers>1) write(output%unit_loop,"('  *** Symmetric and antisymmetric exchange interactions:  ***')")
    do i=1,nmaglayers
      do j=1,nmaglayers
        if(i==j) cycle
        write(output%unit_loop,"(' |--------------------- i = ',i0,'   j = ',i0,' -----------------------|')") mmlayermag(i),mmlayermag(j)
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

    ! Writing into files
    ! Exchange interactions
    call writeCoupling(0.d0)

    ! Closing files
    call closeCouplingFiles()
  end if
  ! end if

  call deallocateCoupling()

  return
end subroutine coupling
