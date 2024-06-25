#ifndef __GPUFFT_H__
#define __GPUFFT_H__

#ifdef __NVCC__
  // NVIDIA / CUFFT :

  #include <cufftw.h>
  #include <cufft.h>
  
  #define gpufftResult      cufftResult
  #define gpufftComplex  cufftComplex
  #define gpufftPlanMany cufftPlanMany
  #define gpufftHandle   cufftHandle
  #define gpufftPlan2d   cufftPlan2d
  #define gpufftExecC2C  cufftExecC2C
  #define gpufftDestroy  cufftDestroy
  #define GPUFFT_C2C     CUFFT_C2C
  #define GPUFFT_FORWARD CUFFT_FORWARD
  #define GPUFFT_SUCCESS CUFFT_SUCCESS
#else
  // AMD / HIP :
  
  #include <hipfft.h>
  
  #define gpufftResult      hipfftResult
  #define gpufftComplex  hipfftComplex
  #define gpufftPlanMany hipfftPlanMany
  #define gpufftHandle   hipfftHandle
  #define gpufftPlan2d   hipfftPlan2d
  #define gpufftExecC2C  hipfftExecC2C
  #define gpufftDestroy  hipfftDestroy
  #define GPUFFT_C2C     HIPFFT_C2C
  #define GPUFFT_FORWARD HIPFFT_FORWARD
  #define GPUFFT_SUCCESS HIPFFT_SUCCESS
#endif


#endif
