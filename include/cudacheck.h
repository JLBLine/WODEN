/*! \file
  File to hold some useful CUDA macros
  Largely taken from the SDK
*/

#ifndef CUDACHECK_H
#define CUDACHECK_H

#if defined(__cplusplus)
#include <cstdio>
#include <cstdlib>
#else
#include <stdio.h>
#include <stdlib.h>
#endif


#define CUDA_SAFE_CALL( call ) do {                                     \
    cudaError_t err = call;						\
    if( cudaSuccess != err ) {						\
      fprintf( stderr, "CUDA Error in file '%s' on line %i : %s.\n",	\
	       __FILE__, __LINE__, cudaGetErrorString( err ) );		\
      exit( EXIT_FAILURE );						\
    } } while( 0 );

#define CUDA_CHECK_KERNEL( errorMessage ) do {	\
    cudaError_t err = cudaGetLastError();	\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",	\
	      errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) ); \
      exit(EXIT_FAILURE);						\
    }									\
    err = cudaDeviceSynchronize();					\
    if( cudaSuccess != err) {						\
      fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",	\
	      errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) ); \
      exit(EXIT_FAILURE);						\
    }									\
  } while( 0 );


#endif
