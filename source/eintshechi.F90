! ---------- Spin disturbance: Energy integration ---------
subroutine eintshechi(e, count)
  use mod_f90_kind, only: double
  use mod_constants, only: zero
  use mod_parameters, only: dim, outputunit, outputunit_loop, lverbose, host
  use EnergyIntegration, only: generate_real_epoints, y, wght, x2, p2, nepoints, pn1
  use mod_susceptibilities, only: chiorb_hf
  use mod_mpi_pars
  implicit none
  real(double), intent(in)    :: e
  integer, intent(in) :: count

  integer :: AllocateStatus
  integer :: i
  real(double) :: start_time,elapsed_time,sizemat,speed
  complex(double), dimension(:,:),allocatable :: Fint
!--------------------- begin MPI vars --------------------
  integer :: ix,ix2,itask
  integer :: ncount
  ncount=dim*dim
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  allocate( Fint(dim,dim), STAT = AllocateStatus )
  if (AllocateStatus/=0) then
    write(outputunit,"('[eintshechi] Not enough memory for: Fint')")
    call MPI_Abort(MPI_Comm_Row,errorcode,ierr)
  end if

! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  ix = myrank_row+1
  itask = numprocs_row ! Number of tasks done initially

  ! Starting to calculate energy integral
  chiorb_hf = zero
  do while(ix <= nepoints)
    if (ix <= pn1) then ! First and second integrations (in the complex plane)
      call sumkshechi(e,y(ix),Fint,0)
      Fint     = Fint*wght(ix)
    else ! Third integration (on the real axis)
      ix2 = ix-pn1
      call sumkshechi(e,x2(ix2),Fint,1)
      Fint   = Fint*p2(ix2)
    end if

    chiorb_hf = chiorb_hf + Fint
    ix = ix + numprocs_row
  end do

  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)
  deallocate(Fint)

  return
end subroutine eintshechi

! -------------------- Spin disturbance: Energy integration --------------------
! -------------- to be used in the calculation of linear SOC chi ---------------
subroutine eintshechilinearsoc(e, count)
  use mod_f90_kind, only: double
  use mod_constants, only: zero
  use mod_parameters, only: dim, outputunit_loop, host, lverbose
  use EnergyIntegration, only: generate_real_epoints,y, wght, x2, p2, nepoints, pn1
  use mod_susceptibilities, only: chiorb_hf,chiorb_hflsoc
  use mod_mpi_pars
  implicit none
  real(double), intent(in)    :: e
  integer, intent(in) :: count
  integer           :: AllocateStatus
  integer           :: i
  real(double)                :: start_time,elapsed_time,sizemat,speed
  complex(double), dimension(:,:),allocatable         :: Fint,Fintlsoc
!--------------------- begin MPI vars --------------------
  integer :: ix,ix2,itask
  integer :: ncount
  ncount=dim*dim
!^^^^^^^^^^^^^^^^^^^^^ end MPI vars ^^^^^^^^^^^^^^^^^^^^^^

  allocate( Fint(dim,dim),Fintlsoc(dim,dim), STAT = AllocateStatus )
  if (AllocateStatus/=0) call abortProgram("[eintshechilinearsoc] Not enough memory for: Fint,Fintlsoc")

! Generating energy points in the real axis for third integration
  call generate_real_epoints(e)

  ix = myrank_row+1
  itask = numprocs ! Number of tasks done initially

  chiorb_hf     = zero
  chiorb_hflsoc = zero
  ! Starting to calculate energy integral
  do while(ix < nepoints)
    if (ix<=pn1) then ! First and second integrations (in the complex plane)
      call sumkshechilinearsoc(e,y(ix),Fint,Fintlsoc,0)
      Fint     = Fint*wght(ix)
      Fintlsoc = Fintlsoc*wght(ix)
    else ! Third integration (on the real axis)
      ix2 = ix-pn1
      call sumkshechilinearsoc(e,x2(ix2),Fint,Fintlsoc,1)
      Fint     = Fint*p2(ix2)
      Fintlsoc = Fintlsoc*p2(ix2)
    end if

    chiorb_hf     = chiorb_hf + Fint
    chiorb_hflsoc = chiorb_hflsoc + Fintlsoc

    ix = ix + numprocs_row
  end do
  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hf, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)
  call MPI_Allreduce(MPI_IN_PLACE, chiorb_hflsoc, ncount, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_Comm_Row_hw, ierr)

  deallocate(Fint,Fintlsoc)

  return
end subroutine eintshechilinearsoc
