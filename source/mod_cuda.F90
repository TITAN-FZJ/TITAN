#include "macros.fi"

module mod_cuda
#ifdef _GPU
  use cudafor,    only: cudaGetErrorString,cudaGetDeviceCount,cudaSetDevice
  use cusolverDn, only: cusolverDnHandle
  implicit none
  integer :: num_gpus
  !! Number of GPUS
  integer :: result
  !! Result from calls
  type(cusolverDnHandle) :: h
  !! Cuda handle

  !> Interface for the C++ diagonalization subroutine
  interface diagonalize_gpu
    subroutine diagonalize_gpu(n,A_d,eval_d,handle) bind(C,name='diagonalize_gpu')
    use mod_kind,   only: dp
    use cusolverDn, only: cusolverDnHandle
    implicit none
    integer, value, intent(in) :: n
    complex(dp), device, dimension(*) :: A_d         !, intent(inout)
    real(dp),    device, dimension(*) :: eval_d      !, intent(out)  
    type(cusolverDnHandle), value     :: handle      !  intent(in)   
    end subroutine diagonalize_gpu
  end interface diagonalize_gpu
  
contains

  ! --------------------------------------------------------------------
  ! subroutine create_acc_handle():
  !    This subruotine creates a CUDA handle for the cusolver zheev subroutine
  ! --------------------------------------------------------------------
  subroutine create_handle()
    use cusolverDn, only: cusolverDnCreate
    implicit none
    integer :: stat

    CUDA_CALL( cusolverDnCreate(h) ) 

  end subroutine create_handle


  ! --------------------------------------------------------------------
  ! subroutine diagonalize_gpu():
  !    This subruotine is an interface for the cusolverDnZheevd subroutine
  ! to diagonalize an hermitian matrix on the GPUs.
  ! --------------------------------------------------------------------
  ! subroutine diagonalize_gpu(n,A_d,eval_d,handle) 
  !   use mod_kind,     only: dp
  !   use cudafor,      only: cudadevicesynchronize
  !   use cusolverDn,   only: CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,cusolverDnZheevd_bufferSize,cusolverDnZheevd
  !   implicit none
  !   integer,                intent(in)    :: n
  !   complex(dp), device   , intent(inout) :: A_d(n,n)
  !   real(dp),    device   , intent(out)   :: eval_d(n)
  !   type(cusolverDnHandle), intent(in)    :: handle
  !   ! Workspace variables
  !   integer,     device :: info
  !   integer             :: lwork
  !   complex(dp), device, allocatable :: work(:) ! dimension of lwork below

  !   CUDA_CALL( cusolverDnZheevd_bufferSize(handle,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,n,A_d,n,eval_d,lwork) )

  !   allocate(work(lwork))

  !   CUDA_CALL( cusolverDnZheevd(handle,CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,n,A_d,n,eval_d,work,lwork,info) )

  !   CUDA_CALL( cudaDeviceSynchronize() )

  !   deallocate(work)

  ! end subroutine diagonalize_gpu


  ! --------------------------------------------------------------------
  ! subroutine create_acc_handle():
  !    This subruotine creates an openacc handle for the cuda LAPACK zheev subroutine
  ! --------------------------------------------------------------------
  subroutine destroy_handle()
    use cusolverDn,   only: cusolverDnDestroy
    use mod_mpi_pars, only: abortProgram
    implicit none

    CUDA_CALL( cusolverDnDestroy(h) )

  end subroutine destroy_handle

#endif
end module mod_cuda