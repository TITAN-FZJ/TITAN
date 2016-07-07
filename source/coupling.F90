! Calculates the full 3x3 J tensor (including coupling, DMI and anisotropic pair interactions)
	subroutine coupling()
		use mod_f90_kind
		use mod_parameters
		use mod_magnet
		use mod_mpi_pars,only: myrank
  character(len=400) :: varm
  character(len=50)  :: fieldpart,socpart
  character(len=1)   :: SOCc
	integer            :: i,j,mu,iw
  real(double),dimension(:,:),allocatable     :: trJij
  real(double),dimension(:,:,:,:),allocatable :: Jij,Jijs,Jija

  allocate(Jij(nmaglayers,nmaglayers,3,3),trJij(nmaglayers,nmaglayers),Jijs(nmaglayers,nmaglayers,3,3),Jija(nmaglayers,nmaglayers,3,3))

	if(myrank.eq.0) write(outputunit,"('CALCULATING FULL TENSOR OF EXHANGE INTERACTIONS AND ANISOTROPIES AS A FUNCTION OF THICKNESS')")

  fieldpart = ""
  socpart   = ""
  if(SOC) then
    if(llinearsoc) then
      SOCc = "L"
    else
      SOCc = "T"
    end if
    write(socpart,"('_magaxis=',a,'_socscale=',f5.2)") magaxis,socscale
  else
    SOCc = "F"
  end if
  if(lfield) then
    write(fieldpart,"('_hwa=',es9.2,'_hwt=',f5.2,'_hwp=',f5.2)") hw_list(hw_count,1),hw_list(hw_count,2),hw_list(hw_count,3)
    if(ltesla) fieldpart = trim(fieldpart) // "_tesla"
  end if

  ! Opening files for position dependence
  if((myrank.eq.0).and.(Npl.eq.Npl_i)) then
    ! Exchange interactions
    do j=1,nmaglayers ; do i=1,nmaglayers
      iw = 199+(j-1)*nmaglayers*2+(i-1)*2
      if(i.eq.j) then
        iw = iw + 1
        write(varm,"('./results/',a1,'SOC/Jij/Jii_',i0,'_parts=',I0,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,i,parts,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=iw, file=varm,status='unknown')
        write(unit=iw, fmt="('#  Npl ,      Jii_xx       ,       Jii_yy  ')")
        iw = iw + 1
        ! TODO : Check how to write the anisotropy term here
      else
        iw = iw + 1
        write(varm,"('./results/',a1,'SOC/Jij/J_',i0,'_',i0,'_parts=',I0,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,i,j,parts,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=iw, file=varm,status='unknown')
        write(unit=iw, fmt="('#  Npl ,   isotropic Jij    ,   anisotropic Jij_xx    ,   anisotropic Jij_yy     ')")
        iw = iw + 1
        write(varm,"('./results/',a1,'SOC/Jij/Dz_',i0,'_',i0,'_parts=',I0,'_ncp=',I0,'_eta=',es8.1,'_Utype=',i0,a,a,'.dat')") SOCc,i,j,parts,ncp,eta,Utype,trim(fieldpart),trim(socpart)
        open (unit=iw, file=varm,status='unknown')
        write(unit=iw, fmt="('#  Npl , Dz = (Jxy - Jyx)/2       ')")
      end if
    end do ; end do
  end if

  q = [ 0.d0 , 0.d0 , 0.d0 ]
  call jij_energy(Jij)

  if(myrank.eq.0) then
    do i=1,nmaglayers ; do j=1,nmaglayers
      trJij(i,j)    = 0.5d0*(Jij(i,j,1,1)+Jij(i,j,2,2))
      Jija(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) - transpose(Jij(i,j,:,:)))
      Jijs(i,j,:,:) = 0.5d0*(Jij(i,j,:,:) + transpose(Jij(i,j,:,:)))
      do mu=1,3
        Jijs(i,j,mu,mu) = Jijs(i,j,mu,mu) - trJij(i,j)
      end do
    end do ; end do

    ! Writing exchange couplings and anisotropies
    write(outputunit,"('  ************************* Full tensor Jij:  *************************')")
    do i=1,nmaglayers ; do j=1,nmaglayers
    ! Writing on screen
    ! Writing original full tensor Jij
    ! Only the transverse components are "reliable" (e.g., for m //z, only Jxx,Jxy,Jyx,Jyy are correct)
    ! Relation between J_ii calculated and the position of the peak in the susceptibility:
    ! w_res = 2*gamma*mz*sqrt( (K_z-K_x)*(K_z-K_y) )  - for m // z
    ! where K_x = J_ii^xx ; K_y = J_ii^yy ; K_z = J_ii^zz
      if(i.eq.j) then
        write(outputunit,"(3x,' ******** Magnetization components: (magaxis = ',a,') *******')") magaxis
        write(outputunit,"(4x,'Mx (',i2.0,')=',f11.8,4x,'My (',i2.0,')=',f11.8,4x,'Mz (',i2.0,')=',f11.8)") i,mx(i),i,my(i),i,mz(i)
        write(outputunit,"(' |--------------- i = ',i0,'   j = ',i0,': anisotropies ---------------|')") mmlayermag(i),mmlayermag(j)
      else
        write(outputunit,"(' |----------- i = ',i0,'   j = ',i0,': exchange couplings -------------|')") mmlayermag(i),mmlayermag(j)
      end if
      write(outputunit,"('             x                  y                  z')")
      write(outputunit,"('  x  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,1,1),Jij(i,j,1,2),Jij(i,j,1,3)
      write(outputunit,"('  y  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,2,1),Jij(i,j,2,2),Jij(i,j,2,3)
      write(outputunit,"('  z  (',es16.9,') (',es16.9,') (',es16.9,')')") Jij(i,j,3,1),Jij(i,j,3,2),Jij(i,j,3,3)
    end do ; end do
    if(nmaglayers.gt.1) write(outputunit,"('  *** Symmetric and antisymmetric exchange interactions:  ***')")
    do i=1,nmaglayers ; do j=1,nmaglayers
      if(i.eq.j) cycle
      write(outputunit,"(' |--------------------- i = ',i0,'   j = ',i0,' -----------------------|')") mmlayermag(i),mmlayermag(j)
    ! Writing Heisenberg exchange interactions
      write(outputunit,"('     Isotropic:     J     = ',es16.9)") trJij(i,j)
      write(outputunit,"('   Anisotropic:     Js_xx = ',es16.9)") Jijs(i,j,1,1)
      write(outputunit,"('                    Js_yy = ',es16.9)") Jijs(i,j,2,2)
      write(outputunit,"('  DMI: Dz = (Jxy - Jyx)/2 = ',es16.9)") Jija(i,j,1,2)
      write(outputunit,"(' --- z components of Jij (not physically correct) ---')")
      write(outputunit,"('  Anisotropic:  Js_zz = ',es16.9)") Jijs(i,j,3,3)
      write(outputunit,"('  DMI: Dy = (Jzx - Jxz)/2 = ',es16.9)") -Jija(i,j,1,3)
      write(outputunit,"('  DMI: Dx = (Jyz - Jzy)/2 = ',es16.9)") Jija(i,j,2,3)
    end do ; end do

    ! Writing into files
    ! Exchange interactions
    exchange_writing_loop: do j=1,nmaglayers ; do i=1,nmaglayers
      iw = 199+(j-1)*nmaglayers*2+(i-1)*2
      if(i.eq.j) then
        iw = iw + 1
        write(unit=iw,fmt="(4x,i3,3x,2(es16.9,2x))") Npl,Jij(i,j,1,1),Jij(i,j,2,2)
        iw = iw + 1
      else
        iw = iw + 1
        write(unit=iw,fmt="(4x,i3,3x,3(es16.9,2x))") Npl,trJij(i,j),Jijs(i,j,1,1),Jijs(i,j,2,2)
        iw = iw + 1
        write(unit=iw,fmt="(4x,i3,3x,es16.9,2x)") Npl,Jija(i,j,1,2)
      end if
    end do ; end do exchange_writing_loop

    ! Closing files
    if(Npl.eq.Npl_f) then
      ! Closing files
      do j=1,nmaglayers ; do i=1,nmaglayers
        iw = 199+(j-1)*nmaglayers*2+(i-1)*2
        if(i.eq.j) then
          iw = iw + 1
          close (iw)
          iw = iw + 1
        else
          iw = iw + 1
          close (iw)
          iw = iw + 1
          close (iw)
        end if
      end do ; end do
    end if
  end if

  deallocate(trJij,Jij,Jijs,Jija)

	return
	end subroutine coupling