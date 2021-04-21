/*! \file
  Device methods to calculate perfect Gaussian and Analytic Dipole primary beam
  responses. Currently, the analytic dipole is fixed to being an MWA dipole.
  Both models assume there is no leakage and beams are purely real.
*/

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
@param[in,out] *d_beam_real Real part of the Gaussian repsonse
@param[in,out] *d_beam_imag Imaginary part of the Gaussian reponse

*/
__device__ void twoD_Gaussian(float x, float y, float xo, float yo,
           float sigma_x, float sigma_y, float cos_theta, float sin_theta,
           float * d_beam_real, float * d_beam_imag);

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
__global__ void kern_gaussian_beam(float *d_beam_ls, float *d_beam_ms,
           float beam_ref_freq, float *d_freqs,
           float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
           int num_freqs, int num_times, int num_components,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J11);

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
     int num_freqs, float ha0, float sdec0, float cdec0,
     float fwhm_lm, float cos_theta, float sin_theta, float sin_2theta,
     float beam_ref_freq, float *d_freqs,
     float *beam_has, float *beam_decs,
     cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J11);

// extern "C" void testing_gaussian_beam( float *beam_has, float *beam_decs,
//            float *beam_angles_array, float *beam_freqs, float *ref_freq_array,
//            float *beam_ls, float *beam_ms,
//            int num_components, int num_times, int num_freqs,
//            float *beam_reals, float *beam_imags);

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

__device__ void analytic_dipole(float az, float za, float wavelength,
           cuFloatComplex * d_beam_X, cuFloatComplex * d_beam_Y);

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
__global__ void kern_analytic_dipole_beam(float *d_azs, float *d_zas,
           float *d_freqs, int num_freqs, int num_times, int num_components,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J11);

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
     float *azs, float *zas, float *d_freqs,
     cuFloatComplex * d_primay_beam_J00, cuFloatComplex * d_primay_beam_J11);

/**
@brief Run `calculate_analytic_dipole_beam` from the host, performing all the
necessary device memory management to retrieve the analytic dipole response
back to the host.

@details Calculates the beam response of a north-south (X) and east-west (Y)
analytic dipole on an infinite ground screen. See
`calculate_analytic_dipole_beam` for details - this function wraps that one
to run from the host, grabbing the device outputs and copying into the host
arrays `analy_beam_X, analy_beam_Y`.

@param[in] num_components Number of COMPONENTS the beam is calculated for
@param[in] num_time_steps Number of time steps being calculated
@param[in] num_freqs Number of frequencies being calculated
@param[in] azs Array of Azimuth angles to calculate the beam towards (radians)
@param[in] zas Array of Zenith Angles to calculate the beam towards (radians)
@param[in] *freqs Array of frequencies to calculate beam at (Hz)
@param[in,out] *analy_beam_X Array to store the north-south beam response in
@param[in,out] *analy_beam_Y Array to store the east-west beam response in

*/
extern "C" void test_analytic_dipole_beam(int num_components,
     int num_time_steps, int num_freqs,
     float *azs, float *zas, float *freqs,
     float _Complex *analy_beam_X, float _Complex *analy_beam_Y);
