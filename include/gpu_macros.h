/*! \file
  A bunch of macros to choose either CUDA or HIP for GPU operations. Needs the compiler to define `__NVCC__` or `__HIPCC__` to choose the correct functionality

  @author Marcin Sokolowski and Cristian Di Pietrantonio, edited by Jack Line
*/

#ifndef __GPU_MACROS_H__
#define __GPU_MACROS_H__

#include <stdio.h>
#include <stdbool.h>

// #ifndef __NVCC__
// #define __NVCC__ // should be set by compiler !!!
// #endif

// TODO: Cristian should know how to avoid being explicit here:
// #define __HIPCC__

#if defined (__NVCC__) || defined (__HIPCC__)

#define __GPU__

// bool gpu_support() { return true;}

// I first define the error handling macro and related definitions. I will
// then use those to wrap all other macros, so that error handling is done
// automatically when using "gpu*" calls.

#ifdef __NVCC__
#define gpuError_t cudaError_t
#define gpuSuccess cudaSuccess
#define gpuGetErrorString cudaGetErrorString
#include <cuda_runtime.h>
#else
#include <hip/hip_runtime.h>
#define gpuError_t hipError_t
#define gpuSuccess hipSuccess
#define gpuGetErrorString hipGetErrorString
#endif

#define GPU_CHECK_ERROR(X)({\
    if(X != gpuSuccess){\
        fprintf(stderr, "GPU error (%s:%d): %s\n", __FILE__ , __LINE__ , gpuGetErrorString(X));\
        exit(1);\
    }\
})

/*! `true` \n
Used within `GPUErrorCheck`. If true, exit if an error is found */
#define EXITERROR true


/**
@brief Take a GPU error message (code), and checks whether an error
occurred.

@details If an error happened, uses `message` to give more information to the
user, along with the decoded CUDA error message. Uses `file` and `line` to
report where the error happened. Optional bool `abort` means you can switch off
exiting if an error is found (default true)

@param[in] message User supplied error message
@param[in] code Error message out of CUDA call (e.g. cudaMalloc)
@param[in] file Name of file call was made in
@param[in] line Line of file call was made in
@param[in] abort If true, exit the CUDA code when an error is found

*/
inline void GPUErrorCheck(const char *message, gpuError_t code, const char *file, int line, bool abort=EXITERROR){
  if (code != gpuSuccess) {
    fprintf(stderr,"GPU ERROR %s: %s\n %s:%d\n",
                    message, gpuGetErrorString(code), file, line);
    if (abort) {
      printf("GPU IS EXITING\n");
      exit(code);
    }
  }
}


#ifdef __NVCC__

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)


#include <cuComplex.h>

#define gpuMalloc(...) GPUErrorCheck("cudaMalloc", cudaMalloc(__VA_ARGS__),__FILE__, __LINE__)
#define gpuHostAlloc(...) GPUErrorCheck("cudaHostAlloc", cudaHostAlloc(__VA_ARGS__, 0),__FILE__, __LINE__)
#define gpuHostAllocDefault cudaHostAllocDefault
#define gpuMemcpy(...) GPUErrorCheck("cudaMemcpy", cudaMemcpy(__VA_ARGS__),__FILE__, __LINE__)
#define gpuMemcpyAsync(...) GPUErrorCheck("cudaMemcpyAsync", cudaMemcpyAsync(__VA_ARGS__),__FILE__, __LINE__)
#define gpuMemset(...) GPUErrorCheck("cudaMemset", cudaMemset(__VA_ARGS__),__FILE__, __LINE__)
// #define gpuDeviceSynchronize(...) GPUErrorCheck("cudaDeviceSynchronize", cudaDeviceSynchronize(__VA_ARGS__),__FILE__, __LINE__)
#define gpuDeviceSynchronize cudaDeviceSynchronize
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
#define gpuFree(...) GPUErrorCheck("cudaFree", cudaFree(__VA_ARGS__),__FILE__, __LINE__)
#define gpuHostFree(...) GPUErrorCheck("cudaFreeHost", cudaFreeHost(__VA_ARGS__),__FILE__, __LINE__)
#define gpuStream_t cudaStream_t
#define gpuStreamCreate(...) GPUErrorCheck("cudaStreamCreate", cudaStreamCreate(__VA_ARGS__),__FILE__, __LINE__)
#define gpuStreamDestroy(...) GPUErrorCheck("cudaStreamDestroy", cudaStreamDestroy(__VA_ARGS__),__FILE__, __LINE__)
#define gpuEventCreate(...) GPUErrorCheck("cudaEventCreate", cudaEventCreate(__VA_ARGS__),__FILE__, __LINE__)
#define gpuGetDeviceCount(...) GPUErrorCheck("cudaGetDeviceCount", cudaGetDeviceCount(__VA_ARGS__),__FILE__, __LINE__)
#define gpuGetLastError cudaGetLastError
#define gpuMemGetInfo(...) GPUErrorCheck("cudaMemGetInfo", cudaMemGetInfo(__VA_ARGS__),__FILE__, __LINE__)
#define gpuMallocHost(...) GPUErrorCheck("cudaMallocHost", cudaMallocHost(__VA_ARGS__),__FILE__, __LINE__)
#define gpuCheckErrors(...) cudaCheckErrors(__VA_ARGS__)
#define gpuFreeHost(...) GPUErrorCheck(" cudaFreeHost",  cudaFreeHost(__VA_ARGS__),__FILE__, __LINE__ )
#define gpuGetDeviceProperties(...) cudaGetDeviceProperties(__VA_ARGS__)
#define gpuDeviceProp cudaDeviceProp
#define gpuPeekAtLastError cudaPeekAtLastError

// Complex number operations:
#define gpuCreal cuCreal
#define gpuCrealf cuCrealf
#define gpuCimag cuCimag
#define gpuCimagf cuCimagf
#define gpuCadd  cuCadd
#define gpuCmul  cuCmul
#define gpuCdiv  cuCdiv
#define gpuConj  cuConj
#define gpuCsub  cuCsub
#define gpuCabs  cuCabs
#define gpuCaddf cuCaddf
#define gpuCsubf cuCsubf
#define gpuCmulf cuCmulf
#define gpuCdivf cuCdivf
#define gpuDoubleComplex cuDoubleComplex
#define gpuFloatComplex cuFloatComplex
#define make_gpuDoubleComplex make_cuDoubleComplex
#define make_gpuFloatComplex make_cuFloatComplex

/*inline int num_available_gpus()
{
    int num_gpus;
    gpuGetDeviceCount(&num_gpus);
    return num_gpus;
} */   


#else

// no need in HIP
// #include <hipComplex.h>
#include <hip/hip_complex.h>

#define hipCheckErrors(msg) \
    do { \
        hipError_t __err = hipGetLastError(); \
        if (__err != hipSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, hipGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)


#define gpuMalloc(...) GPUErrorCheck("hipMalloc", hipMalloc(__VA_ARGS__),__FILE__, __LINE__)
#define gpuHostAlloc(...) GPUErrorCheck("hipHostMalloc", hipHostMalloc(__VA_ARGS__, 0),__FILE__, __LINE__)
#define gpuHostAllocDefault 0
#define gpuMemcpy(...) GPUErrorCheck("hipMemcpy", hipMemcpy(__VA_ARGS__),__FILE__, __LINE__)
#define gpuMemcpyAsync(...) GPUErrorCheck("hipMemcpyAsync", hipMemcpyAsync(__VA_ARGS__),__FILE__, __LINE__)
#define gpuMemset(...) GPUErrorCheck("hipMemset", hipMemset(__VA_ARGS__),__FILE__, __LINE__)
// #define gpuDeviceSynchronize(...) GPUErrorCheck("hipDeviceSynchronize",hipDeviceSynchronize(__VA_ARGS__),__FILE__, __LINE__)
#define gpuDeviceSynchronize hipDeviceSynchronize
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define gpuFree(...) GPUErrorCheck("hipFree", hipFree(__VA_ARGS__),__FILE__, __LINE__)
#define gpuHostFree(...) GPUErrorCheck("hipHostFree", hipHostFree(__VA_ARGS__),__FILE__, __LINE__)
#define gpuStream_t hipStream_t
#define gpuStreamCreate(...) GPUErrorCheck("hipStreamCreate", hipStreamCreate(__VA_ARGS__),__FILE__, __LINE__)
#define gpuStreamDestroy(...) GPUErrorCheck("hipStreamDestroy", hipStreamDestroy(__VA_ARGS__),__FILE__, __LINE__)
#define gpuEventCreate(...) GPUErrorCheck("hipEventCreate", hipEventCreate(__VA_ARGS__),__FILE__, __LINE__)
#define gpuGetDeviceCount(...) GPUErrorCheck("hipGetDeviceCount", hipGetDeviceCount(__VA_ARGS__),__FILE__, __LINE__)
#define gpuGetLastError hipGetLastError
#define gpuMemGetInfo(...) GPUErrorCheck("hipMemGetInfo", hipMemGetInfo(__VA_ARGS__),__FILE__, __LINE__)
#define gpuMallocHost(...) GPUErrorCheck("hipHostMalloc", hipHostMalloc(__VA_ARGS__, 0),__FILE__, __LINE__) // TODO : double check this may be temporary only
#define gpuCheckErrors(...) hipCheckErrors(__VA_ARGS__)
#define gpuFreeHost(...)  GPUErrorCheck( "hipFreeHost", hipFreeHost(__VA_ARGS__),__FILE__, __LINE__ )
#define gpuGetDeviceProperties(...) GPUErrorCheck( "hipGetDeviceProperties", hipGetDeviceProperties(__VA_ARGS__),__FILE__, __LINE__ )
#define gpuDeviceProp hipDeviceProp_t
#define gpuPeekAtLastError hipPeekAtLastError


// Complex number operations:
#define gpuCreal hipCreal
#define gpuCrealf hipCrealf
#define gpuCimag hipCimag
#define gpuCimagf hipCimagf
#define gpuCadd  hipCadd
#define gpuCmul  hipCmul
#define gpuCdiv  hipCdiv
#define gpuConj  hipConj
#define gpuCsub  hipCsub
#define gpuCabs  hipCabs
#define gpuCaddf hipCaddf
#define gpuCsubf hipCsubf
#define gpuCmulf hipCmulf
#define gpuCdivf hipCdivf
#define gpuDoubleComplex hipDoubleComplex
#define gpuFloatComplex  hipFloatComplex
#define make_gpuDoubleComplex make_hipDoubleComplex
#define make_gpuFloatComplex make_hipFloatComplex

#endif
#define gpuCheckLastError(...) GPUErrorCheck(gpuGetLastError())
#else
// bool gpu_support() { return false;}
// inline int num_available_gpus(){ return 0; } 


#endif
#endif



/**
@brief Takes a GPU kernel, runs it with given arguments, and passes results
onto `ErrorCheck` to check for errors.

@details All arguements need to run `kernel` must be included after the listed
arguments here. `kernel` is then run via

        kernel <<< grid , threads, 0 >>>(__VA_ARGS__)

where `__VA_ARGS__` passes on the arguments located at `...`

`gpuErrorCheckKernel` then passes the string `message` on to `ErrorCheck`,
along with the file name and line number via `__FILE__` and `__LINE__`, and
checks the errors from both `gpuGetLastError()` and `gpuDeviceSynchronize()`
after running the kernel.

For example, if `fancy_kernel` takes the arguments `arg1` and `arg2`, to run it
with 10 grids of 64 threads, run the following:

          dim3 grid, threads;
          grid.x = 10
          threads.x = 64
          gpuErrorCheckKernel("Call to fancy_kernel",
                              fancy_kernel, grid, threads,
                              arg1, arg2);



@param[in] message Message to report when an error occurs
@param[in] kernel Name of kernel to be run
@param[in] grid A `dim3` containing grid specifications to run `kernel` with
@param[in] threads A `dim3` containing thread specifications to run `kernel` with
@param[in] ... All arguments to be passed into `kernel`
*/
#define gpuErrorCheckKernel(message, kernel, grid, threads, ...) \
  kernel <<< grid , threads, 0 >>>(__VA_ARGS__); \
  GPUErrorCheck(message, gpuGetLastError(), __FILE__, __LINE__); \
  GPUErrorCheck(message, gpuDeviceSynchronize(), __FILE__, __LINE__);
//   gpuDeviceSynchronize();