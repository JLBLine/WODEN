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
@param[in] x The point at which the function is to be evaluated
@param[in,out] values Array to store calculated function values within
*/
__device__ void RTS_CUDA_pm_polynomial_value_singlef(int n, int m, float x,
                                                     float *values );

/**
@brief This grabs a legendre polynomial value for a `theta` value in `d_theta`
and a given order, and inserts the output into `rts_p1`, and the
output / `sin(theta)` into `rts_P_sin`.

@details The outputs are scaled correctly for the FEE model resolution, and
dumped into `rts_p1`, `rts_P_sin` in order of `theta` value first, order of
legendre polynomial second. `rts_p1` and `rts_P_sin` should be
`num_coords*(nmax*nmax + 2*nmax)` long.

When called with `dim3 grid, threads`, kernel should be called with both
`grid.x` and `grid.y` defined, where:
 - grid.x * threads.x >= `num_coords`
 - grid.y * threads.y >= `nmax*nmax + 2*nmax`

@param[in] d_theta Array of theta values (radians)
@param[in,out] rts_P_sin Output array to store l. poly. value / sin(theta) in
@param[in,out] rts_p1 Output array to store l. poly. value in
@param[in] nmax Maximum order of legendre polynomial to calculate to
@param[in] num_coords Number of theta coords


*/
__global__ void RTS_P1SINfKernel( float *d_theta, cuFloatComplex *rts_P_sin,
           cuFloatComplex *rts_p1, int nmax, int num_coords);

/**
@brief This kernel grabs all the FEE beam coeffs and generated legendre
polynomial values, combines them to get the spherical harmonic values for all
available orders, and converts them into complex field values.

@details `pb_M`, `pb_N`, `pb_Q1`, `pb_Q2` should have been calculated using
`FEE_primary_beam.RTS_MWAFEEInit` and transferred to the device using
`FEE_primary_beam_cuda.copy_FEE_primary_beam_to_GPU`. `rts_P_sin`, `rts_p1` and
should have been calculated using `FEE_primary_beam_cuda.RTS_P1SINfKernel`.
The output arrays `emn_T`,`emn_P` contain the theta and phi polarisation outputs
for both the north-south and east-west dipoles, and contain separate values
for all spherical harmonic values, and so contain `2*num_coords*nMN` values.

When called with `dim3 grid, threads`, kernel should be called with all
`grid.x`, `grid.y`, and `grid.z` defined, where:
 - grid.x * threads.x >= `num_coords`
 - grid.y * threads.y >= `nMN`
 - grid.z * threads.z >= 2 (number of polarisations)

@param[in] d_phi Array of phi values (radians) (these are azimuth values)
@param[in] d_theta Array of theta values (radians) (these are zenith angles)
@param[in] nMN Total number of legendre polynomial values to be used (if `nmax`
is the maximum order to use, `nMN = nmax*nmax + 2*nmax`)
@param[in] num_coords Number of theta/phi coords
@param[in] pb_M Precalculated FEE beam params
@param[in] pb_N Precalculated FEE beam params
@param[in] pb_Q1 Precalculated FEE beam params
@param[in] pb_Q2 Precalculated FEE beam params
@param[in] rts_P_sin Calculated l. poly. value / sin(theta)
@param[in] rts_P1 Calculated l. poly. value in
@param[in,out] emn_T Complex theta polarisation output for all sky coords,
still separated by order of spherical harmonic
@param[in,out] emn_P Complex phi polarisation output for all sky coords,
still separated by order of spherical harmonic

*/
__global__ void RTS_getTileGainsKernel( float *d_phi, float *d_theta,
           int nMN, int num_coords,
           float *pb_M, float *pb_N,
           cuFloatComplex *pb_Q1, cuFloatComplex *pb_Q2,
           cuFloatComplex *rts_P_sin, cuFloatComplex *rts_P1,
           cuFloatComplex *emn_T, cuFloatComplex *emn_P);

/**
@brief This kernel takes comlex polarisation outputs of `RTS_getTileGainsKernel`,
which are separated by spherical harmonic order, and sums over the spherical
harmonics, to return a single complex beam gain direction on sky.

@details The summation happens over the `m` index of the spherical harmonic,
and so we require an array containing the range of possible `m` values
(`d_m_range`) and an array that contains the `m` value of all outputs in
`emn_T`, `emn_P` (`d_M`), to ensure we are summing the correct values and not
double counting anything.

When called with `dim3 grid, threads`, kernel should be called with all
`grid.x` and `grid.y` defined, where:
 - grid.x * threads.x >= `num_coords`
 - grid.y * threads.y >= `2*nmax + 1`

@param[in] emn_T theta polarisation outputs of `RTS_getTileGainsKernel`
@param[in] emn_P phi polarisation outputs of `RTS_getTileGainsKernel`
@param[in,out] d_emn_T_sum_real Real part of summed theta polarisation
@param[in,out] d_emn_T_sum_imag Imaginary part of summed theta polarisation
@param[in,out] d_emn_P_sum_real Real part of summed phi polarisation
@param[in,out] d_emn_P_sum_imag Imaginary part of summed phi polarisation
@param[in] d_m_range All possible `m` values (from -`nmax` to +`nmax`)
@param[in] d_M The `m` value corresponding to values in `emn_T`, `emn_P`
@param[in] nMN Total number of spherical harmonics, `nMN = nmax*nmax + 2*nmax`
@param[in] nmax Maximum order of spherical harmonic
@param[in] num_coords Number of theta/phi coords
*/
__global__ void kern_sum_emn_PT_by_M(cuFloatComplex *emn_T, cuFloatComplex *emn_P,
           float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
           float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
           float *d_m_range, float *d_M, int nMN, int nmax, int num_coords);

/**
@brief This basically just maps the outputs of `kern_sum_emn_PT_by_M` into a
single complex array, ordered by polarisation and dipole orientation,
and then sky coordinate.

@details The summation happens over the `m` index of the spherical harmonic,
and so we require an array containing the range of possible `m` values
(`d_m_range`) and an array that contains the `m` value of all outputs in
`emn_T`, `emn_P` (`d_M`), to ensure we are summing the correct values and not
double counting anything.

Output ordering in `TileGainMatrices` is ordered by
`J_theta_x,J_phi_x,J_theta_y,J_phi_y` (where `theta,phi` is polarisation,
`x,y` are dipole orientation), and these four values are then repeated for each
direction on the sky.

When called with `dim3 grid, threads`, kernel should be called with
`grid.x` and `grid.y` defined, where:
 - grid.x * threads.x >= `num_coords`
 - grid.y * threads.y >= `2`

@param[in] TileGainMatrices A single 1D array that contains all complex beam
gains for all polaristations, dipole orientations, and sky coords
@param[in,out] d_emn_T_sum_real Real part of summed theta polarisation
@param[in,out] d_emn_T_sum_imag Imaginary part of summed theta polarisation
@param[in,out] d_emn_P_sum_real Real part of summed phi polarisation
@param[in,out] d_emn_P_sum_imag Imaginary part of summed phi polarisation
@param[in] nmax Maximum order of spherical harmonic
@param[in] num_coords Number of sky direction coords
*/
__global__ void kern_map_emn(cuFloatComplex *TileGainMatrices,
                float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
                float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
                int nmax, int num_coords);

/**
@brief Normalises the complex beam gains in `TileGainMatrices` to the zenith
reponse of a zenith beam pointing, `d_norm_fac`. The latter array should have
been calcuated using `FEE_primary_beam_cuda.get_HDFBeam_normalisation`.

@details This normalisation should be applied while the polarisations values
in `TileGainMatrices` are still tied to the `phi,theta` coord system (i.e.)
before any parallactic angle rotation). `d_norm_fac` should be an array of 4
values, appropriate to each polarisation in `TileGainMatrices`. The values in
`d_norm_fac` should be the absolute value of the zenith gain, but input as a
`float _Complex` type (so just have imaginary set to zero).

@param[in] TileGainMatrices A single 1D array that contains all complex beam
gains for all polaristations, dipole orientations, and sky coords
@param[in] d_norm_fac Normalisation factor for each polarisation
@param[in] num_coords Number of sky direction coords

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
defined, where:
 - grid.x * threads.x >= `num_coords`
*/
__global__ void kern_apply_FEE_norm(cuFloatComplex *TileGainMatrices,
           cuFloatComplex *d_norm_fac, int num_coords );

/**
@brief Takes the beam gains in `d_FEE_beam_gain_matrices`, which are in `theta`,
`phi` polarisations, and rotates them by the parallactic angle, to align them
into north-south and east-west gain and leakage terms, which can be used
to create XX,XY,YX,YY polarisations.

@details `d_FEE_beam_gain_matrices` should be as output by
`FEE_primary_beam_cuda.kern_map_emn`, ordered by
`J_theta_x,J_phi_x,J_theta_y,J_phi_y` (where `theta,phi` is polarisation,
`x,y` are dipole orientation). `d_sin_para_angs` and `d_cos_para_angs` should
be the sine and cosine of the parallactic angle that corresponds to each sky
direction in `d_FEE_beam_gain_matrices`, repectively. Once the rotation has
been applied, the complex beam values are ordered as `g_x, D_x, D_y, g_y`
where `g` means gain, `D` means leakage, `x` means dipole aligned north-south,
`y` aligned east-west.

The original FEE beam code has a different east-west / north-south convention,
and the sense of azimuth rotation, and so a reordering / sign flip is applied
here to match up with Stokes parameters.

@param[in,out] d_FEE_beam_gain_matrices A single 1D array that contains all complex beam
gains for all polaristations, dipole orientations, and sky coords
@param[in] d_sin_para_angs Sine of the parallactic angle for all sky directions
@param[in] d_cos_para_angs Cosine of the parallactic angle for all sky directions
@param[in] num_coords Number of sky directions in `d_FEE_beam_gain_matrices`

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
defined, where:
 - grid.x * threads.x >= `num_coords`
*/
__global__ void kern_rotate_FEE_beam(cuFloatComplex *d_FEE_beam_gain_matrices,
                                float *d_sin_para_angs, float *d_cos_para_angs,
                                int num_coords);

/**
@brief For the given `phi` (azimuth) and `theta` (zenith angle) sky coords,
calculate the MWA FEE beam model response for all polarisation and dipole
orientations, given the initialised beam model `primary_beam`. Optionally
normalise the beam to zenith (scaling=1) and rotate into a frame compatible
with Stokes parameters (rotation=1).

@details `primary_beam` should have been intialised using
`FEE_primary_beam.RTS_MWAFEEInit` and
`FEE_primary_beam_cuda.copy_FEE_primary_beam_to_GPU`.

If `scaling == 1`, normalises the beam to zenith. For this,
`FEE_primary_beam_cuda.get_HDFBeam_normalisation` should have already be run
to setup the correct attributes in `primary_beam`.

If `rotation == 1`, rotate the beam gains from a instrument sky-locked `theta`,
`phi` polarisation system into one aligned with Stokes parameters (for
creating `XX, XY, YX, YY` beams that can be combined with Stokes I,Q,U,V to
generate instrumental visibilities). For this, `sin_para_angs`,`cos_para_angs`
must be allocated.

Calls the following kernels from `FEE_primary_beam_cuda`. See their
documentation for more detail.

 - `RTS_P1SINfKernel`
 - `RTS_getTileGainsKernel`
 - `kern_sum_emn_PT_by_M`
 - `kern_map_emn`
 - `kern_apply_FEE_norm` (optional)
 - `kern_rotate_FEE_beam` (optional)

@param[in] phi Array of phi values (radians) (these are azimuth values)
@param[in] theta Array of theta values (radians) (these are zenith angles)
@param[in] sin_para_angs Sine of the parallactic angle for all sky directions
@param[in] cos_para_angs Cosine of the parallactic angle for all sky directions
@param[in] num_time_steps Number of time steps in simulation
@param[in] num_components Number of COMPONENTs being used here
@param[in] rotation 0=False, 1=True, rotate results by parallactic angle
@param[in] primary_beam An initialised `RTS_MWA_FEE_beam_t` containing MWA FEE
spherical harmonic coeffs for this pointing
@param[in, out] TileGainMatrices A single 1D array that contains all complex beam
gains for all polaristations, dipole orientations, and sky coords
@param[in] scaling 0=False, 1=True, Normlise results to zenith
*/
extern "C" void RTS_CUDA_get_TileGains(float *phi, float *theta,
           float *sin_para_angs, float *cos_para_angs,
           int num_time_steps, int num_components,
           float rotation, RTS_MWA_FEE_beam_t *primary_beam,
           float _Complex *TileGainMatrices, int scaling);

/**
@brief This function is basically a wrapper to
`FEE_primary_beam_cuda.RTS_CUDA_get_TileGains`, by cudaMalloc-ing an
appropriate array to store the outputs in. If you want to call the MWA FEE
beam, and only use the outputs in device memory, this is the function you want.

@details This function cudaMallocs `FEE_beam->d_FEE_beam_gain_matrices` and sets all array
entries to zero, before passing into `RTS_CUDA_get_TileGains`. See
`RTS_CUDA_get_TileGains` for more detail.

`FEE_beam` should have been intialised using `FEE_primary_beam.RTS_MWAFEEInit`
and `FEE_primary_beam_cuda.copy_FEE_primary_beam_to_GPU`.

If `scaling == 1`, normalises the beam to zenith. For this,
`FEE_primary_beam_cuda.get_HDFBeam_normalisation` should have already be run
to setup the correct attributes in `primary_beam`.

If `rotation == 1`, rotate the beam gains from a instrument sky-locked `theta`,
`phi` polarisation system into one aligned with Stokes parameters (for
creating `XX, XY, YX, YY` beams that can be combined with Stokes I,Q,U,V to
generate instrumental visibilities). For this, `sin_para_angs`,`cos_para_angs`
must be allocated.

@param[in] azs Array of azimuth values (radians)
@param[in] zas Array of zenith angle values (radians)
@param[in] sin_para_angs Sine of the parallactic angle for all az,za
@param[in] cos_para_angs Cosine of the parallactic angle for all az,za
@param[in] num_components Number of COMPONENTs being used here
@param[in] num_time_steps Number of time steps in simulation
@param[in, out] FEE_beam An initialised `RTS_MWA_FEE_beam_t` containing MWA FEE
spherical harmonic coeffs for this pointing
@param[in] rotation 0=False, 1=True, rotate results by parallactic angle
@param[in] scaling 0=False, 1=True, Normlise results to zenith
*/
extern "C" void calc_CUDA_FEE_beam(float *azs, float *zas,
                                   float *sin_para_angs, float *cos_para_angs,
                                   int num_components, int num_time_steps,
                                   RTS_MWA_FEE_beam_t *FEE_beam,
                                   int rotation, int scaling);

/**
@brief We want to normalise by the absolute value of the FEE primary beam,
so take the absolute of the cuFloatComplex values in `d_norm_fac`

@details Beam values to be normalised are stil cuFloatComplex, so easy to just
keep `d_norm_fac` as a complex, and put the absolute value in the real and set
imaginary to zero. There should be 16 values to normalise, as there are four
az/za directions and four polarisations for each

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
set, where:
 - grid.x * threads.x >= 16

@param[in,out] d_norm_fac Complex values to convert to absolute
*/
__global__ void kern_make_norm_abs(cuFloatComplex *d_norm_fac);


/**
@brief Using `FEE_beam_zenith`, calculate the zenith normlisation values for
both dipoles and both polarisations in the native FEE beam coodinate system

@details The FEE beam model is stored in theta/phi (instrument locked)
polarisations. To therefore calculate the zenith normalisation values, actually
need to calculate the beam at zenith with multiple azimuth values, and then
from those directions, select the correct east-west/north south dipole and
theta/phi polarisation output.

This function copies the zenith beam to device, calculates normalistaion values
by calling
 - `copy_FEE_primary_beam_to_GPU`
 - `calc_CUDA_FEE_beam`
 - `kern_make_norm_abs`

and then inserts the correct normalisation values into `FEE_beam->norm_fac`,
which can be used later on by `calc_CUDA_FEE_beam` with `scaling=1`.

Both `FEE_beam_zenith` and `FEE_beam` should have been initialised using
`FEE_primary_beam.RTS_MWAFEEInit`, with the former having all delays set to
zero.

@param[in] FEE_beam_zenith An initialised `RTS_MWA_FEE_beam_t` with zero delays
@param[in] FEE_beam An initialised `RTS_MWA_FEE_beam_t` of the same frequency
in which to store the normalisation outputs for later use
*/
extern "C" void get_HDFBeam_normalisation(RTS_MWA_FEE_beam_t *FEE_beam_zenith,
                RTS_MWA_FEE_beam_t *FEE_beam);

/**
@brief Map the RTS ordered FEE gain array in four separate arrays, each
corresponding to a different element in the beam Jones matrix

@details The visibility functions in `source_components.cu` expect the
primary beam Jones matrix in four 1D cuFloatComplex arrays. This mapping
should be done after rotation by parallatic angle, to be consistent with other
primary beam functions in WODEN. The arrays map to the physical dipoles as

- `d_primay_beam_J00`: north-south gain
- `d_primay_beam_J01`: north-south leakage
- `d_primay_beam_J10`: east-west leakage
- `d_primay_beam_J11`: east-west gain

Within each array, the outputs are stored by time index (slowest changing),
frequency index, then COMPONENT index (fastest changing). This is important
to ensure visibility functions in `source_components.cu` grab the correct
values later on.

When called with `dim3 grid, threads`, kernel should be called with both
`grid.x` and `grid.y` set, where:
 - grid.x * threads.x >= `num_visis`
 - grid.y * threads.y >= `num_components`

@param[in] d_FEE_beam_gain_matrices A single 1D array that contains all complex
beam gains for all polaristations, dipole orientations, and sky coords, as
output by `calc_CUDA_FEE_beam`
@param[in,out] d_primay_beam_J00 Array to store north-south gain on device
memory
@param[in,out] d_primay_beam_J01 Array to store north-south leakage on device
memory
@param[in,out] d_primay_beam_J10 Array to store east-west leakage on device
memory
@param[in,out] d_primay_beam_J11 Array to store east-west gain on device memory
@param[in] num_freqs Number of frequencies
@param[in] num_components Number of COMPONENTs
@param[in] num_visis Total number of visibilities
(`num_freqs*num_baselines*num_times`)
@param[in] num_baselines Number of baselines in the array
@param[in] num_times Number of times steps

*/
__global__ void kern_map_FEE_beam_gains(cuFloatComplex *d_FEE_beam_gain_matrices,
    cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
    cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
    int num_freqs, int num_components, int num_visis, int num_baselines,
    int num_times);

/**
@brief Frees device memory from `primary_beam`

@details Explicitly, cudaFrees:
 - `primary_beam->d_M`
 - `primary_beam->d_N`
 - `primary_beam->d_Q1`
 - `primary_beam->d_Q2`

@param[in] primary_beam A `RTS_MWA_FEE_beam_t` which has already been used on
the device to calculate MWA FEE beam values
*/
extern "C" void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam);

/**
@brief This function is basically a host wrapper to
`FEE_primary_beam_cuda.calc_CUDA_FEE_beam`, by copying the resultant
`FEE_beam->d_FEE_beam_gain_matrices` into host memory. It also calculates the
normlisation factors. If you want to run the MWA FEE beam code on the GPU and
copy outputs onto the CPU, this is the function you want.

@details Both `FEE_beam_zenith` and `FEE_beam` should have been initialised using
`FEE_primary_beam.RTS_MWAFEEInit`, with the former having all delays set to
zero.

If `scaling == 1`, normalises the beam to zenith. For this,
`FEE_primary_beam_cuda.get_HDFBeam_normalisation` should have already be run
to setup the correct attributes in `primary_beam`.

If `rotation == 1`, rotate the beam gains from a instrument sky-locked `theta`,
`phi` polarisation system into one aligned with Stokes parameters (for
creating `XX, XY, YX, YY` beams that can be combined with Stokes I,Q,U,V to
generate instrumental visibilities). For this, `sin_para_angs`,`cos_para_angs`
must be allocated.

@param[in] num_components Number of azimuth and zenith angles
@param[in] azs Array of azimuth values (radians)
@param[in] zas Array of zenith angle values (radians)
@param[in] sin_para_angs Sine of the parallactic angle for all az,za
@param[in] cos_para_angs Cosine of the parallactic angle for all az,za
@param[in, out] FEE_beam_zenith An initialised `RTS_MWA_FEE_beam_t` containing
MWA FEE spherical harmonic coeffs for a zenith pointing
@param[in, out] FEE_beam An initialised `RTS_MWA_FEE_beam_t` containing MWA FEE
spherical harmonic coeffs for desired pointing
@param[in] rotation 0=False, 1=True, rotate results by parallactic angle
@param[in] scaling 0=False, 1=True, Normlise results to zenith
@param[in] FEE_beam_gains Complex array to store outputs in (ordered by
`XX,XY,YX,YY`, and then by az/za)
*/
extern "C" void test_RTS_CUDA_FEE_beam(int num_components,
           float *azs, float *zas, float latitude,
           RTS_MWA_FEE_beam_t *FEE_beam_zenith,
           RTS_MWA_FEE_beam_t *FEE_beam,
           int rotation, int scaling,
           float _Complex *FEE_beam_gains);
