/*! \file cudacheck.h
  Error checking code for CUDA calls and kernels.
*/

#pragma once
#include <stdbool.h>
#include "gpu_macros.h"

//A Bool to say whether we exit if an error is found. Might be useful as an
//option somewhere down the line

/*! `true` \n
Used within `ErrorCheck`. If true, exit if an error is found */
#define EXITERROR true


/**
@brief Take a CUDA error message (code), and checks whether an error
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
inline void ErrorCheck(const char *message, gpuError_t code, const char *file, int line, bool abort=EXITERROR){
  if (code != gpuSuccess) {
    fprintf(stderr,"CUDA ERROR CHECK %s: %s\n %s:%d\n",
                    message, gpuGetErrorString(code), file, line);
    if (abort) {
      printf("CUDA IS EXITING\n");
      exit(code);
    }
  }
}

/**
@brief Take a CUDA error message `code`, and passes it onto `ErrorCheck` to
check for errors. This can be wrapped around any call to `gpuMalloc`,
`gpuMemcpy`, or `gpuFree`.

@details Example usage:

          float *d_array=NULL;
          int num_values=1e6;
          gpuErrorCheckCall( gpuMalloc( (void**)&(d_array), num_values*sizeof(float)) );

Uses `__FILE__` and `__LINE__` to get file name and line
number to pass on to `ErrorCheck`

@param[in] code A `gpuError_t` code
*/
#define gpuErrorCheckCall( code) { ErrorCheck("Call", code, __FILE__, __LINE__); }

//Runs a given kernel with the appropirate number of grids/threads
//Checks for errors and includes "message" in the error message

/**
@brief Takes a CUDA kernel, runs it with given arguments, and passes results
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
//  ErrorCheck(message, gpuGetLastError(), __FILE__, __LINE__); \
  gpuGetLastError(); \\
  ErrorCheck(message, gpuDeviceSynchronize(), __FILE__, __LINE__);
