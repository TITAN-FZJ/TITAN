! This subroutine allocates variables that depend on Npl
subroutine allocate_Npl_variables()
  use mod_parameters
  use mod_magnet
  use mod_tight_binding
  use mod_mpi_pars
  use mod_lattice, only: n0
  implicit none
  integer           :: AllocateStatus

  allocate( sigmai2i(4,Npl),sigmaimunu2i(4,Npl,9,9),sigmaijmunu2i(4,Npl,Npl,9,9),eps1(Npl), STAT = AllocateStatus )
  if (AllocateStatus.ne.0) then
    write(outputunit,"('[main] Not enough memory for: sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  allocate( mx(Npl),my(Npl),mz(Npl),mvec_cartesian(Npl,3),mvec_spherical(Npl,3),hdel(Npl),mp(Npl),hdelp(Npl),mm(Npl),hdelm(Npl), STAT = AllocateStatus )
  if (AllocateStatus.ne.0) then
    write(outputunit,"('[main] Not enough memory for: mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  allocate( mabs(Npl),mtheta(Npl),mphi(Npl),labs(Npl),ltheta(Npl),lphi(Npl),lpabs(Npl),lptheta(Npl),lpphi(Npl), STAT = AllocateStatus )
  if (AllocateStatus.ne.0) then
    write(outputunit,"('[main] Not enough memory for: mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  allocate( mmlayer(Npl+2),layertype(Npl+2),U(Npl+2),mmlayermag(Npl+2),lambda(Npl+2),npart0(Npl+2), STAT = AllocateStatus )
  if (AllocateStatus.ne.0) then
    write(outputunit,"('[main] Not enough memory for: mmlayer,layertype,U,mmlayermag,lambda,npart0')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  allocate( hhwx(Npl+2),hhwy(Npl+2),hhwz(Npl+2),sb(Npl+2,18,18),lb(Npl+2,18,18), STAT = AllocateStatus )
  if (AllocateStatus.ne.0) then
    write(outputunit,"('[main] Not enough memory for: hhwx,hhwy,hhwz,sb,lb')")
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ierr)
  end if
  allocate(sha_longitudinal(n0),sha_transverse(n0),long_cos(n0),transv_cos(n0))

  return
end subroutine allocate_Npl_variables

! This subroutine allocates variables that depend on Npl
subroutine deallocate_Npl_variables()
  use mod_parameters
  use mod_magnet
  use mod_tight_binding
  use mod_lattice
  implicit none

  deallocate(sigmai2i,sigmaimunu2i,sigmaijmunu2i,eps1)
  deallocate(mx,my,mz,mvec_cartesian,mvec_spherical,hdel,mp,hdelp,mm,hdelm)
  deallocate(mabs,mtheta,mphi,labs,ltheta,lphi,lpabs,lptheta,lpphi)
  if(lGSL) deallocate(lxm,lym,lzm,lxpm,lypm,lzpm)
  deallocate(mmlayer,layertype,U,mmlayermag,lambda,npart0)
  deallocate(hhwx,hhwy,hhwz,sb,lb)
  select case (plnn)
  case(1)
    deallocate(t00,t01)
  case(2)
    deallocate(t00,t01,t02)
  end select
  deallocate(sha_longitudinal,sha_transverse,long_cos,transv_cos)

  return
end subroutine deallocate_Npl_variables