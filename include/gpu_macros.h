#ifndef __GPU_MACROS_H__
#define __GPU_MACROS_H__

#include <stdio.h>

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

// void __gpu_check_error(gpuError_t x, const char *file, int line);
/*inline void __gpu_check_error(gpuError_t x, const char *file, int line){
    if(x != gpuSuccess){
        fprintf(stderr, "GPU error (%s:%d): %s\n", file, line, gpuGetErrorString(x));
        exit(1);
    }
}*/


/*#define GPU_CHECK_ERROR(X)({\
    __gpu_check_error((X), __FILE__, __LINE__);\
})*/

#define GPU_CHECK_ERROR(X)({\
    if(X != gpuSuccess){\
        fprintf(stderr, "GPU error (%s:%d): %s\n", __FILE__ , __LINE__ , gpuGetErrorString(X));\
        exit(1);\
    }\
})


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

#define gpuMalloc(...) GPU_CHECK_ERROR(cudaMalloc(__VA_ARGS__))
#define gpuHostAlloc(...) GPU_CHECK_ERROR(cudaHostAlloc(__VA_ARGS__, 0))
#define gpuHostAllocDefault cudaHostAllocDefault
#define gpuMemcpy(...) GPU_CHECK_ERROR(cudaMemcpy(__VA_ARGS__))
#define gpuMemcpyAsync(...) GPU_CHECK_ERROR(cudaMemcpyAsync(__VA_ARGS__))
#define gpuMemset(...) GPU_CHECK_ERROR(cudaMemset(__VA_ARGS__))
#define gpuDeviceSynchronize(...) GPU_CHECK_ERROR(cudaDeviceSynchronize(__VA_ARGS__))
#define gpuMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define gpuMemcpyHostToDevice cudaMemcpyHostToDevice
#define gpuMemcpyDeviceToDevice cudaMemcpyDeviceToDevice
#define gpuFree(...) GPU_CHECK_ERROR(cudaFree(__VA_ARGS__))
#define gpuHostFree(...) GPU_CHECK_ERROR(cudaFreeHost(__VA_ARGS__))
#define gpuStream_t cudaStream_t
#define gpuStreamCreate(...) GPU_CHECK_ERROR(cudaStreamCreate(__VA_ARGS__))
#define gpuStreamDestroy(...) GPU_CHECK_ERROR(cudaStreamDestroy(__VA_ARGS__))
#define gpuEventCreate(...) GPU_CHECK_ERROR(cudaEventCreate(__VA_ARGS__))
#define gpuGetDeviceCount(...) GPU_CHECK_ERROR(cudaGetDeviceCount(__VA_ARGS__))
#define gpuGetLastError cudaGetLastError
#define gpuMemGetInfo(...) GPU_CHECK_ERROR(cudaMemGetInfo(__VA_ARGS__))
#define gpuMallocHost(...) GPU_CHECK_ERROR(cudaMallocHost(__VA_ARGS__))
#define gpuCheckErrors(...) cudaCheckErrors(__VA_ARGS__)
#define gpuFreeHost(...) GPU_CHECK_ERROR( cudaFreeHost(__VA_ARGS__) )
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


#define gpuMalloc(...) GPU_CHECK_ERROR(hipMalloc(__VA_ARGS__))
#define gpuHostAlloc(...) GPU_CHECK_ERROR(hipHostMalloc(__VA_ARGS__, 0))
#define gpuHostAllocDefault 0
#define gpuMemcpy(...) GPU_CHECK_ERROR(hipMemcpy(__VA_ARGS__))
#define gpuMemcpyAsync(...) GPU_CHECK_ERROR(hipMemcpyAsync(__VA_ARGS__))
#define gpuMemset(...) GPU_CHECK_ERROR(hipMemset(__VA_ARGS__))
#define gpuDeviceSynchronize(...) GPU_CHECK_ERROR(hipDeviceSynchronize(__VA_ARGS__))
#define gpuMemcpyDeviceToHost hipMemcpyDeviceToHost
#define gpuMemcpyHostToDevice hipMemcpyHostToDevice
#define gpuMemcpyDeviceToDevice hipMemcpyDeviceToDevice
#define gpuFree(...) GPU_CHECK_ERROR(hipFree(__VA_ARGS__))
#define gpuHostFree(...) GPU_CHECK_ERROR(hipHostFree(__VA_ARGS__))
#define gpuStream_t hipStream_t
#define gpuStreamCreate(...) GPU_CHECK_ERROR(hipStreamCreate(__VA_ARGS__))
#define gpuStreamDestroy(...) GPU_CHECK_ERROR(hipStreamDestroy(__VA_ARGS__))
#define gpuEventCreate(...) GPU_CHECK_ERROR(hipEventCreate(__VA_ARGS__))
#define gpuGetDeviceCount(...) GPU_CHECK_ERROR(hipGetDeviceCount(__VA_ARGS__))
#define gpuGetLastError hipGetLastError
#define gpuMemGetInfo(...) GPU_CHECK_ERROR(hipMemGetInfo(__VA_ARGS__))
#define gpuMallocHost(...) GPU_CHECK_ERROR(hipHostMalloc(__VA_ARGS__, 0)) // TODO : double check this may be temporary only
#define gpuCheckErrors(...) hipCheckErrors(__VA_ARGS__)
#define gpuFreeHost(...)  GPU_CHECK_ERROR( hipFreeHost(__VA_ARGS__) )
#define gpuGetDeviceProperties(...) GPU_CHECK_ERROR( hipGetDeviceProperties(__VA_ARGS__) )
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
#define gpuCheckLastError(...) GPU_CHECK_ERROR(gpuGetLastError())
#else
// bool gpu_support() { return false;}
// inline int num_available_gpus(){ return 0; } 


#endif
#endif
