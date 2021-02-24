// #include <stdio.h>
// #include <stdlib.h>


// #define CUDA_CALL( call )               \
// {                                       \
// cudaError_t result = call;              \
// assert (cudaSuccess == result);         \
// }

#include <iostream>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>

extern "C" void diagonalize_gpu(const int n, cuDoubleComplex* const A_d, double*const eval_d, cusolverDnHandle_t handle)
{
  // Initializing variables
  int  lwork = 0;
  int *devInfo = NULL;
  cuDoubleComplex *work_d = NULL;
  
  // Allocating devInfo
  cudaError_t cudaStat = cudaMalloc (&devInfo, sizeof(int)) ;
  assert(cudaSuccess == cudaStat);

  // Selecting modes of diagonalization
  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

  // Querying working space for diagonalization
  cusolverStatus_t cusolver_status = cusolverDnZheevd_bufferSize(
    handle,
    jobz,
    uplo,
    n,
    A_d,
    n,
    eval_d,
    &lwork) ;
  assert (cusolver_status == CUSOLVER_STATUS_SUCCESS);

  // Allocating working space
  cudaStat = cudaMalloc(&work_d, sizeof(cuDoubleComplex)*lwork);
  assert(cudaSuccess == cudaStat);

  // Computing spectrum
  cusolver_status = cusolverDnZheevd(
    handle,
    jobz,
    uplo,
    n,
    A_d,
    n,
    eval_d,
    work_d,
    lwork,
    devInfo);
  assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

  // Synchronizing
  cudaStat = cudaDeviceSynchronize();
  assert(cudaSuccess == cudaStat);

  // Freeing memory
  cudaFree(devInfo);
  cudaFree(work_d);

  return;
}