// #include <stdio.h>
// #include <stdlib.h>


// #define CUDA_CALL( call )               \
// {                                       \
// cudaError_t result = call;              \
// assert (cudaSuccess == result);         \
// }

#include <assert.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>


extern "C" void diagonalize_gpu(int n, cuDoubleComplex *A_d, double *eval_d, cusolverDnHandle_t handle)
{
  int *devInfo = NULL;
  cuDoubleComplex *work_d = NULL;
  int  lwork = 0;
  cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
  cudaError_t cudaStat = cudaSuccess;

  cudaStat = cudaMalloc ((void**)&devInfo, sizeof(int)) ;
  assert(cudaSuccess == cudaStat);

// step 3: query working space of syevd
  cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
  cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
  cusolver_status = cusolverDnZheevd_bufferSize(
    handle,
    jobz,
    uplo,
    n,
    A_d,
    n,
    eval_d,
    &lwork) ;
  assert (cusolver_status == CUSOLVER_STATUS_SUCCESS);

  cudaStat = cudaMalloc((void**)&work_d, sizeof(double)*lwork) ;
  assert(cudaSuccess == cudaStat);

// step 4: compute spectrum
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

  cudaStat = cudaDeviceSynchronize();
  assert(cudaSuccess == cudaStat);

  if (devInfo) cudaFree(devInfo);
  if (work_d ) cudaFree(work_d);

}



