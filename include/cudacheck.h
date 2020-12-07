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
