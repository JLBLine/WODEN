#pragma once

#include <stdbool.h>

#define EXITERROR true

#define cudaErrorCheckCall(code) { ErrorCheck("Call", code, __FILE__, __LINE__); }

inline void ErrorCheck(const char *message, cudaError_t code, const char *file,
                       int line, bool abort=EXITERROR){
  if (code != cudaSuccess) {
    fprintf(stderr,"CUDA ERROR CHECK %s: %s\n %s:%d\n",
                    message, cudaGetErrorString(code), file, line);
    if (abort) {
      printf("CUDA IS EXITING\n");
      exit(code);
    }
  }
}

//Runs a given kernel with the appropirate number of grids/threads
//Checks for errors and includes "message" in the error message

#define cudaErrorCheckKernel(message, kernel, grid, threads, ...) \
  kernel <<< grid , threads, 0 >>>(__VA_ARGS__); \
  ErrorCheck(message, cudaGetLastError(), __FILE__, __LINE__); \
  ErrorCheck(message, cudaDeviceSynchronize(), __FILE__, __LINE__);

// #define cudaErrorCheckKernel(kernel, grid, threads, ...){
//   kernel<<< grid, threads >>>( __VA_ARGS__ );
//   cudaError_t err = cudaGetLastError();
//   ErrorCheck(err, __FILE__, __LINE__);
//   err = cudaDeviceSynchronize();
//   ErrorCheck(err, __FILE__, __LINE__);
// }


// #define CUDA_SAFE_CALL( call ) do {                                     \
//     cudaError_t err = call;						\
//     if( cudaSuccess != err ) {						\
//       fprintf( stderr, "CUDA Error in file '%s' on line %i : %s.\n",	\
// 	       __FILE__, __LINE__, cudaGetErrorString( err ) );		\
//       exit( EXIT_FAILURE );						\
//     } } while( 0 );
