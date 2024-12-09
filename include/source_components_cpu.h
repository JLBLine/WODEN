/*! \file
  CPU methods to extrapolate flux densities, assign and apply primary beam
  gains, and visibility responses for different sky model COMPONENT types.
*/

#include "woden_precision_defs.h"

// #if defined (__NVCC__) || defined (__HIPCC__)
// #include "gpu_macros.h"
// #endif


// /*!
// A struct to contain primary beam values for a give COMPONENT. `d_gxs,d_Dxs,d_Dys,d_gys`
// should be used when all antennas have the same primary beam, and `d_gxs,d_Dxs,d_Dys,d_gys` used when all primary beams are different.
// */
// typedef struct _beam_gains_t {

//   //OK, so when this is compiled for CPU, the compiler doesn't have access
//   //to GPU libraries. So wrap the GPU specific stuff in an ifdef
//   //that only triggers when we're compiling for GPU
//   #if defined (__NVCC__) || defined (__HIPCC__)
//   gpuUserComplex *d_gxs = NULL; /*!< Device copy of North-South Beam gain values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/
//   gpuUserComplex *d_Dxs = NULL; /*!< Device copy of North-South Beam leakage values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/
//   gpuUserComplex *d_Dys = NULL; /*!< Device copy of East-West Beam leakage values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/
//   gpuUserComplex *d_gys = NULL; /*!< Device copy of East-West Beam gain values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/
//   #endif

//   user_precision_complex_t *gxs = NULL; /*!< Host copy of North-South Beam gain values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/
//   user_precision_complex_t *Dxs = NULL; /*!< Host copy of North-South Beam leakage values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/
//   user_precision_complex_t *Dys = NULL; /*!< Host copy of East-West Beam leakage values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/
//   user_precision_complex_t *gys = NULL; /*!< Host copy of East-West Beam gain values
//   for all beams, directions, frequencies, and times for these COMPONENTS*/


//   int *d_ant1_to_baseline_map = NULL; /*!< The index of antenna 1 in all unique pairs of
// antennas. Used to map iBaseline to the correct antenna 1 */
//   int *d_ant2_to_baseline_map = NULL; /*!< The index of antenna 2 in all unique pairs of
// antennas. Used to map iBaseline to the correct antenna 2 */
//   int use_twobeams; /*!< The beam gains were made with unique primary beams so
//   should use two antenna patterns per visibility */

// } beam_gains_t;