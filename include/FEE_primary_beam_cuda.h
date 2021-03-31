/*! \file
  Device methods to calculate the MWA FEE primary beam response
  @author J.L.B. Line
*/
// #define kP1BlockSize 64
#include "cudacheck.h"
#include "woden_struct_defs.h"

/**
@brief Copy the spherical harmonic MWA FEE coeffs from host to device

@details Copies `FEE_beam->M`, `FEE_beam->N`, `FEE_beam->Q1`, `FEE_beam->Q2`
into `FEE_beam->d_M`, `FEE_beam->d_N`, `FEE_beam->d_Q1`, `FEE_beam->d_Q2`,
reshaping from 2D arrays into 1D arrays.

@param[in,out] FEE_beam `RTS_MWA_FEE_beam_t` which has been initialised with
`FEE_primary_beam.RTS_MWAFEEInit`

*/
extern "C" void copy_FEE_primary_beam_to_GPU(RTS_MWA_FEE_beam_t *FEE_beam);

/**
@brief RTS GPU implemtation of John Burkardt's
`legendre_polynomial.pm_polynomial_value`. Evaluates all legendre polynomials
up to the specified orders, for a single evaluation point `x`.
Original documentation from `pm_polynomial_value` included below.

@details Note the original function evaluated for mulitple `x` vales, so needed
a number of evalutation points `mm`, which is set to 1.0 in this RTS GPU code.

Purpose:

  PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).

Differential equation:

  (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0

First terms:

  M = 0  ( = Legendre polynomials of first kind P(N,X) )

  Pm(0,0,x) =    1
  Pm(1,0,x) =    1 X
  Pm(2,0,x) = (  3 X^2 -   1)/2
  Pm(3,0,x) = (  5 X^3 -   3 X)/2
  Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
  Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
  Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
  Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16

  M = 1

  Pm(0,1,x) =   0
  Pm(1,1,x) =   1 * SQRT(1-X^2)
  Pm(2,1,x) =   3 * SQRT(1-X^2) * X
  Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
  Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)

  M = 2

  Pm(0,2,x) =   0
  Pm(1,2,x) =   0
  Pm(2,2,x) =   3 * (1-X^2)
  Pm(3,2,x) =  15 * (1-X^2) * X
  Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)

  M = 3

  Pm(0,3,x) =   0
  Pm(1,3,x) =   0
  Pm(2,3,x) =   0
  Pm(3,3,x) =  15 * (1-X^2)^1.5
  Pm(4,3,x) = 105 * (1-X^2)^1.5 * X

  M = 4

  Pm(0,4,x) =   0
  Pm(1,4,x) =   0
  Pm(2,4,x) =   0
  Pm(3,4,x) =   0
  Pm(4,4,x) = 105 * (1-X^2)^2

Recursion:

  if N < M:
    Pm(N,M,x) = 0
  if N = M:
    Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
    all the odd integers less than or equal to N.
  if N = M+1:
    Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
  if M+1 < N:
    Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)

Licensing:

  This code is distributed under the GNU LGPL license.

Modified:

  08 August 2013

Author:

  John Burkardt

Reference:

  Milton Abramowitz, Irene Stegun,
  Handbook of Mathematical Functions,
  National Bureau of Standards, 1964,
  ISBN: 0-486-61272-4,
  LC: QA47.A34.

Parameters:

  Input, int MM, the number of evaluation points.

  Input, int N, the maximum first index of the Legendre
  function, which must be at least 0.

  Input, int M, the second index of the Legendre function,
  which must be at least 0, and no greater than N.

  Input, double X[MM], the point at which the function is to be
  evaluated.

  Output, double PM_POLYNOMIAL_VALUE[MM*(N+1)], the function values.


@param[in] n The maximum first index of the Legendre function, which must
be at least 0
@param[in] m The second index of the Legendre function,which must be at least 0,
and no greater than N
@param[in] x The point at which the function is to be
evaluated
@param[in,out] values Array to store calculated function values within
*/


__device__ void RTS_CUDA_pm_polynomial_value_singlef(int n, int m, float x,
                                                     float *values );

/**
@brief This this grabs legendre polynomial values for all theta values

@details

@param[in]
*/
__global__ void RTS_P1SINfKernel( float *d_theta, cuFloatComplex *rts_P_sin,
           cuFloatComplex *rts_p1, int nmax, int num_coords);

/**
@brief I think this grabs all the FEE beam coeffs and generated legendre
polynomial values and combines them to get the spherical harmonic values for all
avaible orders, and converts them into complex field values

@details

@param[in]
*/
__global__ void RTS_getTileGainsKernel( float *d_phi, float *d_theta,
           int nMN, int num_coords,
           float *pb_M, float *pb_N,
           cuFloatComplex *pb_Q1, cuFloatComplex *pb_Q2,
           cuFloatComplex *rts_P_sin, cuFloatComplex *rts_P1,
           cuFloatComplex *emn_T, cuFloatComplex *emn_P);

/**
@brief I think this kernel takes the individual spherical harmonic values,
which are still separated by their \f$n,m\f$ values, and sums over them to
get the beam response for a given direction on the sky.

@details I changed this to be an `atomicAdd` to save memory

@param[in]
*/
__global__ void kern_sum_emn_PT_by_M(cuFloatComplex *emn_T, cuFloatComplex *emn_P,
           float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
           float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
           float *d_m_range, float *d_M, int nMN, int nmax, int num_coords);

/**
@brief This basically just maps the outputs of the beam codes into a single
complex array, ordered by polarisation and then sky coordinate

@details

@param[in]
*/
__global__ void kern_calc_sigmaTP(cuFloatComplex *TileGainMatrices,
                float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
                float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
                int nmax, int num_coords);

/**
@brief

@details

@param[in]
*/
__global__ void kern_apply_FEE_norm(cuFloatComplex *TileGainMatrices,
           cuFloatComplex *d_norm_fac, int num_coords );

/**
@brief

@details

@param[in]
*/
__global__ void kern_rotate_FEE_beam(cuFloatComplex *d_FEE_beam_gain_matrices,
                                float *d_sin_para_angs, float *d_cos_para_angs,
                                int num_components, int num_time_steps);

/**
@brief

@details

@param[in]
*/
extern "C" void RTS_CUDA_get_TileGains(float *phi, float *theta,
           float *sin_para_angs, float *cos_para_angs,
           int num_time_steps, int num_components,
           float rotation, RTS_MWA_FEE_beam_t *primary_beam,
           float _Complex *TileGainMatrices, int scaling);

extern "C" void calc_CUDA_FEE_beam(float *azs, float *zas,
                                   float *sin_para_angs, float *cos_para_angs,
                                   int num_components, int num_time_steps,
                                   RTS_MWA_FEE_beam_t *FEE_beam,
                                   int rotation, int scaling);
/**
@brief

@details

@param[in]
*/
extern "C" void get_HDFBeam_normalisation(RTS_MWA_FEE_beam_t *FEE_beam_zenith,
                RTS_MWA_FEE_beam_t *FEE_beam);

/**
@brief

@details

@param[in]
*/
__global__ void kern_map_FEE_beam_gains(cuFloatComplex *d_FEE_beam_gain_matrices,
    cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
    cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
    int num_freqs, int num_components, int num_visis, int num_baselines,
    int num_times);

/**
@brief

@details

@param[in]
*/
extern "C" void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

/**
@brief

@details

@param[in]
*/
extern "C" void test_RTS_CUDA_FEE_beam(int num_components,
           float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           RTS_MWA_FEE_beam_t *FEE_beam_zenith,
           RTS_MWA_FEE_beam_t *FEE_beam,
           int rotation, int scaling,
           float _Complex *FEE_beam_gains);
