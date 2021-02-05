#pragma once
#include <stdbool.h>

//A Bool to say whether we exit if an error is found. Might be useful as an
//option somewhere down the line
#define EXITERROR true

//Take a CUDA error (code), and checks whether an error occurred. If an error
//happens, uses (message) to give more information to the user. Uses (file)
//and (line) to report where the error happened. Optional bool (abort) means
//you can switch off exiting if an error is found.
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

//Function to wrap any CUDA call in and error check
#define cudaErrorCheckCall(code) { ErrorCheck("Call", code, __FILE__, __LINE__); }

//Runs a given kernel with the appropirate number of grids/threads
//Checks for errors and includes "message" in the error message
#define cudaErrorCheckKernel(message, kernel, grid, threads, ...) \
  kernel <<< grid , threads, 0 >>>(__VA_ARGS__); \
  ErrorCheck(message, cudaGetLastError(), __FILE__, __LINE__); \
  ErrorCheck(message, cudaDeviceSynchronize(), __FILE__, __LINE__);
