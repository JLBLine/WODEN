/*! \file
  Handle the error messages that come out of hyperbeam version >= 0.5.0
*/
// #pragma once
// #include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mwa_hyperbeam.h>
#include "logger.h"

// #ifndef __hyperbeam_error_h__
// #define __hyperbeam_error_h__

//Going to be calling this code in both C, and CUDA, so stick the
//conditional extern C around so linkage by compilers survives all
//the mangling fun
// #ifndef __HIPCC__
#ifdef __cplusplus
  extern "C" {
#endif
// #endif

/**
@brief Use functions out of `mwa_hyperbeam` to handle errors out of said package

@details Given the name of the source code `file[]`, which can be passed
via `__FILE__`, and the line of the file `line_num`, which can be passed as
`__LINE__`, and the name of the `mwa_hyperbeam` function being called
`function_name[]`, make sense of the error that occured and print out something
useful, the exit.

@param[in] file[] Name of source code file (passed normally via `__FILE__`)
@param[in] line_num Line number of function call inside source code (passed
normally via `__LINE__`)
@param[in] function_name[] Name of `mwa_hyperbeam` function being called

*/
void handle_hyperbeam_error(const char file[], int line_num, const char function_name[]);

// #ifndef __HIPCC__
#ifdef __cplusplus
}
#endif

// #endif // __hyperbeam_error_h__
