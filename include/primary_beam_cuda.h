/*! \file
  Device methods to calculate perfect Gaussian and Analytic Dipole primary beam
  responses. Currently, the analytic dipole is fixed to being an MWA dipole.
  Both models assume there is no leakage and beams are purely real.
*/
#pragma once
#include "woden_precision_defs.h"
#include <mwa_hyperbeam.h>
#include "cudacomplex.h"
#include "hyperbeam_error.h"

/**
@brief Calculate a two dimensional Gaussian

@details Returns the Gaussian as defined on Wikipedia
[here](https://en.wikipedia.org/wiki/Gaussian_function) via the equation:

\f[
G(x,y) = \exp \left( -\left( a(x-x_o)^2 + 2b(x-x_o)(y-y_o) + c(y-y_o)^2 \right)  \right)
\f]

where

\f{eqnarray*}{
a  &=&  \frac{\cos(\theta)^2}{2\sigma_x^2} + \frac{\sin(\theta)^2}{2\sigma_y^2} \\
b  &=&  -\frac{\sin(2\theta)}{4\sigma_x^2} + \frac{\sin(2\theta)}{4\sigma_y^2} \\
c  &=&  \frac{\sin(\theta)^2}{2\sigma_x^2} + \frac{\cos(\theta)^2}{2\sigma_y^2}
\f}

with \f$(x_0,y_0)\f$ the central coordinates, \f$(\sigma_x,\sigma_y)\f$ the
standard deviations.

@param[in] x x coordinate
@param[in] y y coordinate
@param[in] xo Central x coordinate
@param[in] yo Central y coordinate
@param[in] sigma_x Square root of variance in x
@param[in] sigma_y Square root of variance in x
@param[in] cos_theta Cosine of the rotation angle
@param[in] sin_theta Sine of the rotation angle
@param[in] sin_2theta Sine of two times the rotation angle
@param[in,out] *d_beam_real Real part of the Gaussian repsonse
@param[in,out] *d_beam_imag Imaginary part of the Gaussian reponse
*/
__device__ void twoD_Gaussian(user_precision_t x, user_precision_t y,
           user_precision_t xo, user_precision_t yo,
           user_precision_t sigma_x, user_precision_t sigma_y,
           user_precision_t cos_theta, user_precision_t sin_theta,
           user_precision_t sin_2theta,
           user_precision_t * d_beam_real, user_precision_t * d_beam_imag);

/**
@brief Kernel to calculate a Gaussian primary beam response at the given
interferometric sky \f$(l,m)\f$ coords and frequency

@details The primary beam is to be calculated for each sky direction, each time
step, and each frequency. The size of the beam repsonse on the sky changes with
frequency, so need a reference frequency `beam_ref_freq` and full-width
half-maximum (`fwhm_lm`, in \f$l,m\f$ units) to scale the beam width with frequency.
The Gaussian beam is held to point at a given az/za, so the sky positions to
calculate in `d_beam_ls, d_beam_ms` should contain `num_components*num_times`
values, as the COMPONENTs move through the beam with time. The outputs are
stored in `d_primay_beam_J00, d_primay_beam_J11`, where `00` refers to the
north-south polarisation, `11` the east-west polarisation, in order of time,
frequency, COMPONENT.

When called with `dim3 grid, threads`, kernel should be called with both
`grid.x` and `grid.y` defined, where:
 - grid.x * threads.x >= `num_components` * `num_time_steps`
 - grid.y * threads.y >= `num_freqs`

@param[in] *d_beam_ls Array of \f$l\f$ coords to calculate beam at
@param[in] *d_beam_ms Array pf \f$m\f$ coords to calculate beam at
@param[in] beam_ref_freq Reference frequency at which the FWHM is applicable (Hz)
@param[in] *d_freqs Array of frequencies to calculate beam at (Hz)
@param[in] fwhm_lm FWHM of the beam in \f$l,m\f$ coords
@param[in] cos_theta Cosine of the rotation angle
@param[in] sin_theta Sine of the rotation angle
@param[in] sin_2theta Sine of two times the rotation angle
@param[in] num_freqs Number of frequencies being calculated
@param[in] num_times Number of time steps being calculated
@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in,out] *d_primay_beam_J00 Device array to store the beam Jones
complex `J[0,0]` response in
@param[in,out] *d_primay_beam_J11 Device array to store the beam Jones
complex `J[1,1]` response in

*/
__global__ void kern_gaussian_beam(double *d_beam_ls, double *d_beam_ms,
           double beam_ref_freq, double *d_freqs,
           user_precision_t fwhm_lm, user_precision_t cos_theta,
           user_precision_t sin_theta, user_precision_t sin_2theta,
           int num_freqs, int num_times, int num_components,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11);

/**
@brief Calculate the Gaussian primary beam response at the given hour angle and
declinations `beam_point_has, beam_point_decs`. Note the XX and YY repsonses
are equal in this toy example.

@details The primary beam is to be calculated for each sky direction, each time
step, and each frequency. The size of the beam repsonse on the sky changes with
frequency, so need a reference frequency `beam_ref_freq` and full-width
half-maximum (`fwhm_lm`, in \f$l,m\f$ units) to scale the beam width with frequency.
The Gaussian beam is held to point at a given az/za, so the sky positions to
calculate in `d_beam_ls, d_beam_ms` should contain `num_components*num_times`
values, as the COMPONENTs move through the beam with time. The outputs are
stored in `d_primay_beam_J00, d_primay_beam_J11`, where `00` refers to the
north-south polarisation, `11` the east-west polarisation, in order of time,
frequency, COMPONENT. This function uses the beam centre pointing `ha0, dec0`
to calculate an \f$l,m\f$ coord system in which to calculate the Gaussian beam,
using `kern_gaussian_beam`.

@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in] num_time_steps Number of time steps being calculated
@param[in] num_freqs Number of frequencies being calculated
@param[in] ha0 Hour angle of beam pointing centre (radians)
@param[in] sdec0 Sine of Declination of the beam pointing centre (radians)
@param[in] cdec0 Cosine of Declination of the beam pointing centre (radians)
@param[in] fwhm_lm FWHM of the beam in \f$l,m\f$ coords
@param[in] cos_theta Cosine of the rotation angle
@param[in] sin_theta Sine of the rotation angle
@param[in] sin_2theta Sine of two times the rotation angle
@param[in] beam_ref_freq Reference frequency at which the FWHM is applicable (Hz)
@param[in] *d_freqs Array of frequencies to calculate beam at (Hz)
@param[in] *beam_has Array of Hour Angles to calculate the beam toward
@param[in] *beam_decs Array of Declinations to calculate the beam toward
@param[in,out] *d_primay_beam_J00 Device array to store the beam Jones
complex `J[0,0]` response in
@param[in,out] *d_primay_beam_J11 Device array to store the beam Jones
complex `J[1,1]` response in

*/
extern "C" void calculate_gaussian_beam(int num_components, int num_time_steps,
           int num_freqs, user_precision_t ha0,
           user_precision_t sdec0, user_precision_t cdec0,
           user_precision_t fwhm_lm, user_precision_t cos_theta,
           user_precision_t sin_theta, user_precision_t sin_2theta,
           double beam_ref_freq, double *d_freqs,
           double *beam_has, double *beam_decs,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11);

/**
@brief Calculate the beam response of a north-south (X) and east-west (Y)
analytic dipole on an infinite ground screen, for the given sky direciton
`az,za` and `wavelength`.

@details Dipoles are assumed to be MWA, and given a length of 0.3 metres. Beam
size on the sky scales with frequency hence the need for `wavelength`

@param[in] az Azimuth (radians)
@param[in] za Zenith Angle (radians)
@param[in] wavelength Wavelength (metres)
@param[in,out] d_beam_X Complex beam value for north-south dipole
@param[in,out] d_beam_Y Complex beam value for east-west dipole

*/
__device__ void analytic_dipole(user_precision_t az, user_precision_t za,
           user_precision_t wavelength,
           cuUserComplex * d_beam_X, cuUserComplex * d_beam_Y);

/**
@brief Kernel to calculate an Analytic MWA Dipole over an infinite ground screen
at the given Azimuth and Zenith Angles `d_azs, d_zas` and frequencies `d_freqs`.

@details The primary beam is to be calculated for each sky direction, each time
step, and each frequency. The Analytic dipole beam is stationary on the sky, so
the Azimuth and Zenith Angles in `azs,zas` should contain
`num_components*num_times` values, as the COMPONENTs move through the beam with
time. The outputs are stored in `d_primay_beam_J00, d_primay_beam_J11`,
where `00` refers to the north-south polarisation, `11` the east-west
polarisation, in order of time, frequency, COMPONENT. Beam outputs are
normalised to zenith

When called with `dim3 grid, threads`, kernel should be called with both
`grid.x` and `grid.y` defined, where:
 - grid.x * threads.x >= `num_components` * `num_time_steps`
 - grid.y * threads.y >= `num_freqs`

@todo Make the zenith normalisation an option

@param[in] *d_azs Array of Azimuth angles to calculate the beam towards (radians)
@param[in] *d_zas Array of Zenith Angles to calculate the beam towards (radians)
@param[in] *d_freqs Array of frequencies to calculate beam at (Hz)
@param[in] num_freqs Number of frequencies being calculated
@param[in] num_times Number of time steps being calculated
@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in,out] *d_primay_beam_J00 Device array to store the beam Jones
complex `J[0,0]` response in
@param[in,out] *d_primay_beam_J11 Device array to store the beam Jones
complex `J[1,1]` response in

*/
__global__ void kern_analytic_dipole_beam(user_precision_t *d_azs,
           user_precision_t *d_zas,  double *d_freqs, int num_freqs,
           int num_times, int num_components,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11);

/**
@brief Calculate the Analytic Dipole over an infinite ground screen primary beam
response at the given Azimuth and Zenith Angles `azs,zas`, and frequencies
`d_freqs`.

@details The primary beam is to be calculated for each sky direction, each time
step, and each frequency. The Analytic dipole beam is stationary on the sky, so
the Azimuth and Zenith Angles in `azs,zas` should contain
`num_components*num_times` values, as the COMPONENTs move through the beam with
time. The outputs are stored in `d_primay_beam_J00, d_primay_beam_J11`,
where `00` refers to the north-south polarisation, `11` the east-west
polarisation, in order of time, frequency, COMPONENT. Beam outputs are
normalised to zenith. Note eveything starting with `d_` should be in device
memory.

When called with `dim3 grid, threads`, kernel should be called with both
`grid.x` and `grid.y` defined, where:
 - grid.x * threads.x >= `num_components` * `num_time_steps`
 - grid.y * threads.y >= `num_freqs`

@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in] num_time_steps Number of time steps being calculated
@param[in] num_freqs Number of frequencies being calculated
@param[in] *d_freqs Array of frequencies to calculate beam at (Hz)
@param[in] *azs Array of Azimuth angles to calculate the beam towards (radians)
@param[in] *zas Array of Zenith Angles to calculate the beam towards (radians)
@param[in,out] *d_primay_beam_J00 Device array to store the beam Jones
complex `J[0,0]` response in
@param[in,out] *d_primay_beam_J11 Device array to store the beam Jones
complex `J[1,1]` response in

*/
extern "C" void calculate_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, double *d_freqs,
     cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J11);

/**
@brief Calculate the analytic MWA beam response to a single az/za direction
for a give wavelength and set of delays (`d_metre_delays`)

@details Based on the RTS code found in ``mwa_tile::local_FillMatrices``. I
think this code inherently does some kind of parallatic rotation,
hence needs the ha/dec along with the az/za. The delays added to the
paths of individual dipoles allow the beam to be 'pointed'. The physical
length of these paths should be given in `d_metre_delays` (this conversion
from the delays given in the MWA metafits is handled by
`primary_beam_cuda::calculate_RTS_MWA_analytic_beam`)


@param[in] az Azimuth to calculate the beam toward (radians)
@param[in] za Zenith angle to calculate the beam toward (radians)
@param[in] ha Corresponding hour angle to the given `az, za`
@param[in] dec Corresponding Declination to the given `az, za`
@param[in] wavelength Wavelength to calculate the beam at (metres)
@param[in] d_metre_delays The physical delays lengths added to each dipole to
steer the beam (metres)
@param[in] latitude Latitude of the array (radians)
@param[in] norm Whether to normalise to zenith or not (1 to normalise, 0 to not)
@param[in,out] *gx The gain for the north-south beam
@param[in,out] *Dx The leakage for the north-south beam
@param[in,out] *Dy The leakage for the east-west beam
@param[in,out] *gy The gain for the east-west beam

*/
__device__ void RTS_MWA_beam(user_precision_t az, user_precision_t za,
           double ha, double dec,
           double wavelength, double *d_metre_delays,
           double latitude, int norm,
           cuUserComplex * gx, cuUserComplex * Dx,
           cuUserComplex * Dy, cuUserComplex * gy);

/**
@brief Kernel to calculate the analytic MWA primary beam to a set of sky directions
`d_azs` and `d_zas` for a given set of delays `d_metre_delays` and frequencies
`d_freqs`.

@details Kernel calls `primary_beam_cuda::RTS_MWA_beam`. The MWA primary beam
is stationary on the sky for a given set of delays, so the Azimuth and Zenith
Angles in `azs,zas` should contain `num_components*num_times` values,
as the COMPONENTs move through the beam with time. The delays added to the
paths of individual dipoles allow the beam to be 'pointed'. The physical
length of these paths should be given in `d_metre_delays` (this conversion
from the delays given in the MWA metafits is handled by
`primary_beam_cuda::calculate_RTS_MWA_analytic_beam`)

When called with `dim3 grid, threads`, kernel should be called with both
`grid.x` and `grid.y` defined, where:
 - grid.x * threads.x >= `num_components` * `num_time_steps`
 - grid.y * threads.y >= `num_freqs`

@param[in] *d_azs Array of Azimuth angles to calculate the beam towards (radians)
@param[in] *d_zas Array of Zenith Angles to calculate the beam towards (radiany)
@param[in] *d_beam_has Hour angles corresponding to `azs,zas` (radians)
@param[in] *d_beam_decs Declinations corresponding to `azs,zas` (radians)
@param[in] *d_metre_delays Delays that specify the pointing (path length, metres)
@param[in] *d_freqs Array of frequencies to calculate beam at (Hz)
@param[in] latitude Latitude of the array (radians)
@param[in] norm Whether to normalise to zenith or not (1 to normalise, 0 to not)
@param[in] num_freqs Number of frequencies being calculated
@param[in] num_times Number of time steps being calculated
@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in,out] *d_gxs The gains for the north-south beam
@param[in,out] *d_Dxs The leakages for the north-south beam
@param[in,out] *d_Dys The leakages for the east-west beam
@param[in,out] *d_gys The gains for the east-west beam

*/
__global__ void kern_RTS_analytic_MWA_beam(user_precision_t *d_azs,
           user_precision_t *d_zas,
           user_precision_t *d_beam_has, user_precision_t *d_beam_decs,
           double *d_metre_delays,
           double *d_freqs, double latitude, int norm,
           int num_freqs, int num_times, int num_components,
           cuUserComplex *d_gxs, cuUserComplex *d_Dxs,
           cuUserComplex *d_Dys, cuUserComplex *d_gys);

/**
@brief Calculate the analytic MWA primary beam to a set of sky directions
`azs` and `zas` for a given set of delays `delays` as listed in an MWA metafits
file, and frequencies `d_freqs`.

@details Uses the kernel `primary_beam_cuda::kern_RTS_analytic_MWA_beam`.
The MWA primary beam is stationary on the sky for a given set of delays,
so the Azimuth and Zenith Angles in `azs,zas` should contain
`num_components*num_times` values, as the COMPONENTs move through the beam
with time. The delays added to the paths of individual dipoles allow the beam
to be 'pointed'. The delays as listed in the metafits (given as `delays`) are
listed in units of the delay time added internally to the tile (in seconds).
This function coverts them into a path length (metres), as needed by
`primary_beam_cuda::RTS_MWA_beam``.

@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in] num_time_steps Number of time steps being calculated
@param[in] num_freqs Number of frequencies being calculated
@param[in] *azs Array of Azimuth angles to calculate the beam towards (radians)
@param[in] *zas Array of Zenith Angles to calculate the beam towards (radiany)
@param[in] *delays Delays that specificy the pointing (as reported in the MWA metafits)
@param[in] latitude Latitude of the array (radians)
@param[in] norm Whether to normalise to zenith or not (1 to normalise, 0 to not)
@param[in] *beam_has Hour angles corresponding to `azs,zas` (radians)
@param[in] *beam_decs Declinations corresponding to `azs,zas` (radians)
@param[in] *d_freqs Array of frequencies to calculate beam at (Hz)
@param[in,out] *d_gxs The gains for the north-south beam
@param[in,out] *d_Dxs The leakages for the north-south beam
@param[in,out] *d_Dys The leakages for the east-west beam
@param[in,out] *d_gys The gains for the east-west beam

*/
extern "C" void calculate_RTS_MWA_analytic_beam(int num_components,
     int num_time_steps, int num_freqs,
     user_precision_t *azs, user_precision_t *zas, user_precision_t *delays,
     double latitude, int norm,
     double *beam_has, double *beam_decs, double *d_freqs,
     cuUserComplex *d_gxs, cuUserComplex *d_Dxs,
     cuUserComplex *d_Dys, cuUserComplex *d_gys);








/**
@brief Calculate the FEE MWA primary beam model to a set of sky directions
`azs` and `zas` for a given initialised `mwa_hyperbeam` device beam object
`*cuda_fee_beam`. NOTE that the azs, zas need to increment by component
(fastest changing), then time (slowest changing). This is the OPPOSITE
of what happens for all other functions, but is needed for pointer arithmatic
that must be done to feed things into `mwa_hyperbeam` efficiently. Soz boz.

@details Calls `mwa_hyperbeam::calc_jones_cuda_device` to calculate the beam
responses on the device. This function requires an initialised
`struct FEEBeamGpu *cuda_fee_beam` object (initialised using
`mwa_hyperbeam::new_cuda_fee_beam`), which in turn needs a
`struct FEEBeam *fee_beam` (initialised using `mwa_hyperbeam::new_fee_beam`).
Running these functions gathers the spherical harmnoic coefficients for the
requested frequencies, as well as the delays to point the beam. If these aren't
setup correctly, this will fall flat on it's face.

Once the beam repsonses have been calculated, split them up into the `WODEN`
d_primay_beam_J* arrays using the kernel `primary_beam_cuda::kern_map_hyperbeam_gains`.

@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in] num_time_steps Number of time steps being calculated
@param[in] num_freqs Number of frequencies being calculated
@param[in] num_ants How many antennas(tiles) are in the array. If making all
primary beams the same, set to 1
@param[in] parallactic Whether to rotate by parallactic angle or not
@param[in] *cuda_fee_beam An initialised `mwa_hyperbeam` `struct FEEBeamGpu`
@param[in] *azs Array of Azimuth angles to calculate the beam towards (radians)
@param[in] *zas Array of Zenith Angles to calculate the beam towards (radians)
@param[in] *latitudes The latitude of the array for each time step (radians); this
can be NULL is parallactic = 0
@param[in,out] *d_primay_beam_J00 The gains for the north-south beam
@param[in,out] *d_primay_beam_J01 The leakages for the north-south beam
@param[in,out] *d_primay_beam_J10 The leakages for the east-west beam
@param[in,out] *d_primay_beam_J11 The gains for the east-west beam

*/
extern "C" void run_hyperbeam_cuda(int num_components,
           int num_time_steps, int num_freqs,
           int num_ants, uint8_t parallactic,
           struct FEEBeamGpu *cuda_fee_beam,
           double *azs, double *zas,
           double *latitudes,
           cuUserComplex *d_primay_beam_J00,
           cuUserComplex *d_primay_beam_J01,
           cuUserComplex *d_primay_beam_J10,
           cuUserComplex *d_primay_beam_J11);
