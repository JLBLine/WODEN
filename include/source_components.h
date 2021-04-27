#include "calculate_visibilities.h"
/*! \file
  Device methods to extrapolate flux densities, assign and apply primary beam
  gains, and visibility responses for different sky model COMPONENT types.
*/

/**
@brief Assuming a simple power-law SED of \f$ S \propto \nu^\alpha \f$,
extrapolate Stokes flux density parameters at the requested allsteps_wavelengths

@details The 4 arrays of reference flux densities (`d_ref_stokesI`,
`d_ref_stokesQ`, `d_ref_stokesU`, `d_ref_stokesV`), reference frequencies
`d_ref_freqs`, and spectral indexes `d_SIs` should all be the same length.
`iComponent` indexes these 6 arrays. `iBaseline` indexes the array
`d_allsteps_wavelengths`, containging the wavelengths to extrapolate the
reference fluxes to. Returns the flux density values extrapolated to the
selected wavelength and reference data (`flux_I`, `flux_Q`, `flux_U`, `flux_V`).

@param[in] d_allsteps_wavelengths Wavelengths array to index with iBaseline (metres)
@param[in] d_ref_freqs Component reference frequencies to index with iBaseline (Hz)
@param[in] d_ref_stokesI Stokes I reference flux densities to reference with
iComponent (Jy)
@param[in] d_ref_stokesQ Stokes Q reference flux densities to reference with
iComponent (Jy)
@param[in] d_ref_stokesU Stokes U reference flux densities to reference with
iComponent (Jy)
@param[in] d_ref_stokesV Stokes V reference flux densities to reference with
iComponent (Jy)
@param[in] d_SIs Spectral indicies to index with iBaseline
@param[in] iComponent Component index
@param[in] iBaseline Baseline index
@param[in,out] *flux_I Pointer to extrapolated Stokes I (Jy)
@param[in,out] *flux_Q Pointer to extrapolated Stokes Q (Jy)
@param[in,out] *flux_U Pointer to extrapolated Stokes U (Jy)
@param[in,out] *flux_V Pointer to extrapolated Stokes V (Jy)

*/
__device__ void extrap_stokes(float *d_allsteps_wavelengths, float *d_ref_freqs,
           float *d_ref_stokesI, float *d_ref_stokesQ,
           float *d_ref_stokesU, float *d_ref_stokesV,
           float *d_SIs, int iComponent, int iBaseline,
           float * flux_I, float * flux_Q, float * flux_U, float * flux_V);

/**
@brief Calculate the complex exponential phase part of the visibility function
from the given u,v,w and l,m,n coord arrays

@details Calculates the following complex exponential:

\f[
\exp \left( 2\pi i\left( ul + vm + w(n-1) \right) \right)
\f]

where

    u  =  d_us[iBaseline]
    v  =  d_vs[iBaseline]
    w  =  d_ws[iBaseline]
    l  =  d_ns[iComponent]
    m  =  d_ms[iComponent]
    n  =  d_ns[iComponent]

Note there is no negative sign here - testing with the current implementation
through imaging software, and from comparisons to the OSKAR and RTS simulators
this yields the correct output visibilities.

@param[in] *d_us Pointer to \f$u\f$ coodinates (wavelengths)
@param[in] *d_vs Pointer to \f$v\f$ coodinates (wavelengths)
@param[in] *d_ws Pointer to \f$w\f$ coodinates (wavelengths)
@param[in] *d_ls Pointer to \f$l\f$ coodinates
@param[in] *d_ms Pointer to \f$m\f$ coodinates
@param[in] *d_ns Pointer to \f$n\f$ coodinates
@param[in] iBaseline Index for \f$u,v,w\f$
@param[in] iComponent Index for \f$l,m,n\f$
@return `visi`, a `cuFloatComplex` of the visibility
*/
__device__ cuFloatComplex calc_measurement_equation(float *d_us,
           float *d_vs, float *d_ws, float *d_ls, float *d_ms, float *d_ns,
           int iBaseline, int iComponent);

/**
@brief Given primary beam gains and leakage terms for antenna 1
`g1x, D1x, D1y, gy` and antenna 2 `g1x, D1x, D1y, gy`,the complex visibility
phase across those two antennas `visi`, and the Stokes parameters of a source,
simulate the observed XX,XY,YX,YY instrumental cross-correlated visibilities,
where 'x' means aligned north-south, 'y' means east-west.

@details Performs the following calculations:

\f{eqnarray*}{
\mathrm{V}^{XX}_{12} = (g_{1x}g_{2x}^{\ast} + D_{1x}D_{2x}^{\ast})\mathrm{V}^{I}_{12}
  +  (g_{1x}g_{2x}^{\ast} - D_{1x}D_{2x}^{\ast})\mathrm{V}^{Q}_{12} \\
  +  (g_{1x}D_{2x}^{\ast} + D_{1x}g_{2x}^{\ast})\mathrm{V}^{U}_{12}
  +  i(g_{1x}D_{2x}^{\ast} - D_{1x}g_{2x}^{\ast})\mathrm{V}^{V}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{XY}_{12} =
     (g_{1x}D_{2y}^{\ast} + D_{1x}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
  +  (g_{1x}D_{2y}^{\ast} - D_{1x}g_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
  +  (g_{1x}g_{2y}^{\ast} + D_{1x}D_{2y}^{\ast})\mathrm{V}^{U}_{12}
  +  i(g_{1x}g_{2y}^{\ast} - D_{1x}D_{2y}^{\ast})\mathrm{V}^{V}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YX}_{12} =
     (D_{1y}g_{2x}^{\ast} + g_{1y}D_{2x}^{\ast})\mathrm{V}^{I}_{12}
  +  (D_{1y}g_{2x}^{\ast} - g_{1y}D_{2x}^{\ast})\mathrm{V}^{Q}_{12} \\
  +  (D_{1y}D_{2x}^{\ast} + g_{1y}g_{2x}^{\ast})\mathrm{V}^{U}_{12}
  +  i(D_{1y}D_{2x}^{\ast} - g_{1y}g_{2x}^{\ast})\mathrm{V}^{V}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YY}_{12} =
     (D_{1y}D_{2y}^{\ast} + g_{1y}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
  +  (D_{1y}D_{2y}^{\ast} - g_{1y}g_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
  +  (D_{1y}g_{2y}^{\ast} + g_{1y}D_{2y}^{\ast})\mathrm{V}^{U}_{12}
  +  i(D_{1y}g_{2y}^{\ast} - g_{1y}D_{2y}^{\ast})\mathrm{V}^{V}_{12}
\f}

where \f${\ast}\f$ means complex conjugate, and

\f{eqnarray*}{
\mathrm{V}_{12} &=& \mathrm{V}_{\mathrm{env}}\exp \left( 2\pi i\left( u_{12}l + v_{12}m + w_{12}(n-1) \right) \right) \\
\mathrm{V}^{I}_{12} &=& I(l,m) \mathrm{V}_{12} \\
\mathrm{V}^{Q}_{12} &=& Q(l,m) \mathrm{V}_{12} \\
\mathrm{V}^{U}_{12} &=& U(l,m) \mathrm{V}_{12} \\
\mathrm{V}^{V}_{12} &=& V(l,m) \mathrm{V}_{12}
\f}

where \f$ \mathrm{V}_{\mathrm{env}} \f$ is an visibility envelope that turns a
component into a GAUSSIAN or SHAPELET component, and has been applied in
`visi_component`.

@param[in] g1x Beam gain antenna 1 in north-south
@param[in] D1x Beam leakage antenna 1 from north-south
@param[in] D1y Beam gain antenna 1 in east-west
@param[in] g1y Beam leakage antenna 1 from east-west
@param[in] g2x Beam gain antenna 2 in north-south
@param[in] D2x Beam leakage antenna 2 from north-south
@param[in] D2y Beam gain antenna 2 in east-west
@param[in] g2y Beam leakage antenna 2 from east-west
@param[in] flux_I Stokes I flux density (Jy)
@param[in] flux_Q Stokes Q flux density (Jy)
@param[in] flux_U Stokes U flux density (Jy)
@param[in] flux_V Stokes V flux density (Jy)
@param[in] visi_component Complex visibility across antennas 1 and 2
@param[in,out] visi_XX Output XX instrumental visibility
@param[in,out] visi_XY Output XY instrumental visibility
@param[in,out] visi_YX Output YX instrumental visibility
@param[in,out] visi_YY Output YY instrumental visibility
*/
__device__ void apply_beam_gains(cuFloatComplex g1x, cuFloatComplex D1x,
          cuFloatComplex D1y, cuFloatComplex g1y,
          cuFloatComplex g2x, cuFloatComplex D2x,
          cuFloatComplex D2y, cuFloatComplex g2y,
          float flux_I, float flux_Q,
          float flux_U, float flux_V,
          cuFloatComplex visi_component,
          cuFloatComplex * visi_XX, cuFloatComplex * visi_XY,
          cuFloatComplex * visi_YX, cuFloatComplex * visi_YY);

/**
@brief Given the type of primary beam simulated `beamtype`, select the beam
gains and leakages from the arrays `d_primay_beam_J*` that match the indexes
`iBaseline`, `iComponent`. NOTE that currently, primary beams are assumed
to be the same for all recieving elements, and so this function only takes in
one set of primary beam arrays.

@details This function is built to return the correct beam gain for a given
component on the sky, at a given time, at a given frequency, for a given baseline.
The 4 arrays `d_primay_beam_J00`, `d_primay_beam_J01`, `d_primay_beam_J10`,
`d_primay_beam_J11` should contain the primary beam settings for all times,
all frequencies, and COPMONENTs on the sky, and so should be
`num_freqs*num_components*num_times` long. The order elements should increment
through COMPONENT (fastest changing), freqeuncy, and time (slowest changing).
`iBaseline` is a combined index of baseline, frequency, and time.

If `beamtype == NO_BEAM`, set `g1x = g1y = g2x = g2y = 1` and
`D1x = D1y = D2x = D2y = 0`

@todo Prepare this function for two different primary beam responses

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_primay_beam_J00 Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *d_primay_beam_J01 Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_primay_beam_J10 Pointer towards array of primary beam J[1,0]
(east-west gain)
@param[in] *d_primay_beam_J11 Pointer towards array of primary beam J[1,1]
(east-west leakage)
@param[in,out] *g1x Beam gain antenna 1 in north-south
@param[in,out] *D1x Beam leakage antenna 1 from north-south
@param[in,out] *D1y Beam gain antenna 1 in east-west
@param[in,out] *g1y Beam leakage antenna 1 from east-west
@param[in,out] *g2x Beam gain antenna 2 in north-south
@param[in,out] *D2x Beam leakage antenna 2 from north-south
@param[in,out] *D2y Beam gain antenna 2 in east-west
@param[in,out] *g2y Beam leakage antenna 2 from east-west

*/
__device__ void get_beam_gains(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           cuFloatComplex * g1x, cuFloatComplex * D1x,
           cuFloatComplex * D1y, cuFloatComplex * g1y,
           cuFloatComplex * g2x, cuFloatComplex * D2x,
           cuFloatComplex * D2y, cuFloatComplex * g2y);

/**
@brief Given the visibility between two recieving elements for a COMPONENT
`visi_component`, at the baseline/freq/time index given by `iBaseline` and
COMPONENT index `iComponent`, apply the instrumental primary beam and
COMPONENT Stokes parameters to create instrumental XX,XY,YX,YY visibilities,
and sum them into real and imaginary XX,XY,YX,YY visibilities arrays
`d_sim_visi_*_real` and `d_sim_visi_*_imag`.

@details Uses `get_beam_gains` and `apply_beam_gains` as described above to
apply the gains - see descriptions for what should be the arguments to them.

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_primay_beam_J00 Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *d_primay_beam_J01 Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_primay_beam_J10 Pointer towards array of primary beam J[1,0]
(east-west gain)
@param[in] *d_primay_beam_J11 Pointer towards array of primary beam J[1,1]
(east-west leakage)
@param[in] visi_component Complex visibility across antennas 1 and 2
@param[in] flux_I Stokes I flux density (Jy)
@param[in] flux_Q Stokes Q flux density (Jy)
@param[in] flux_U Stokes U flux density (Jy)
@param[in] flux_V Stokes V flux density (Jy)
@param[in,out] *d_sum_visi_XX_real Pointer to array to sum real XX visibility
into
@param[in,out] *d_sum_visi_XX_imag Pointer to array to sum imaginary XX
visibility into
@param[in,out] *d_sum_visi_XY_real Pointer to array to sum real XY visibility
into
@param[in,out] *d_sum_visi_XY_imag Pointer to array to sum imaginary XY
visibility into
@param[in,out] *d_sum_visi_YX_real Pointer to array to sum real YX visibility
into
@param[in,out] *d_sum_visi_YX_imag Pointer to array to sum imaginary YX
visibility into
@param[in,out] *d_sum_visi_YY_real Pointer to array to sum real YY visibility
into
@param[in,out] *d_sum_visi_YY_imag Pointer to array to sum imaginary YY
visibility into
*/
__device__ void update_sum_visis(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           cuFloatComplex visi_component,
           float flux_I, float flux_Q, float flux_U, float flux_V,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag);

/**
@brief Performs necessary calculations that are common to all POINT, GAUSSIAN,
and SHAPELET component types, calculating the \f$l,m,n\f$ coordinates and
primary beam values.

@details Uses `kern_calc_lmn` to calculate \f$l,m,n\f$. Uses either
`primary_beam_cuda.calculate_gaussian_beam`,
`FEE_primary_beam_cuda.calc_CUDA_FEE_beam`, or
`primary_beam_cuda.calculate_analytic_dipole_beam` to calculate the primary
beam.

Device arrays `d_primay_beam_J*` are used to store primary beam outputs, and
should have
`num_components*woden_settings->num_time_steps*woden_settings->num_freqs`
elements memory allocated. \f$l,m,n\f$ coordinate arrays `d_ls`, `d_ms`,`d_ns`
should have `num_components` elements. `sin_para_angs` and `cos_para_angs` are
sine and cosine of the parallactic angle, and only needed if
`beam_settings.beamtype == FEE_BEAM`. They should be
`num_components*woden_settings->num_time_steps` long. `beam_has`,`beam_decs`
are hour angles and declinations of components at all time steps, and only used
when `beam_settings.beamtype == GAUSS_BEAM`. They should be also be
`num_components*woden_settings->num_time_steps`.

@param[in] num_components Number of COMPONENTs
@param[in] *d_primay_beam_J00 Pointer towards array of primary beam J[0,0]
(north-south gain) to be populated
@param[in] *d_primay_beam_J01 Pointer towards array of primary beam J[0,1]
(north-south leakage) to be populated
@param[in] *d_primay_beam_J10 Pointer towards array of primary beam J[1,0]
(east-west gain) to be populated
@param[in] *d_primay_beam_J11 Pointer towards array of primary beam J[1,1]
(east-west leakage) to be populated
@param[in] d_freqs Frequencies used in the simulation
@param[in,out] d_ls Pointer to array to store l coords in
@param[in,out] d_ms Pointer to array to store m coords in
@param[in,out] d_ns Pointer to array to store d coords in
@param[in] d_ras Array of Right Ascensions of COPMONENTs (radians)
@param[in] d_decs Array of Declinations of COPMONENTs (radians)
@param[in] azs Azimuth angles for all time steps for all COMPONENTs (radians)
@param[in] zas Zenith Angles for all time steps for all COMPONENTs
@param[in] sin_para_angs Sine of parallactic angle for all time steps for all
COMPONENTs
@param[in] cos_para_angs Sine of parallactic angle for all time steps for all
COMPONENTs
@param[in] beam_has Hour angles for all time steps for all COMPONENTs (radians)
@param[in] beam_decs Declinations for all time steps for all COMPONENTs (radians)
@param[in] woden_settings A populated `woden_settings_t` struct
@param[in] beam_settings A populated `beam_settings_t` struct
@param[in] FEE_beam An initialised `RTS_MWA_FEE_beam_t` struct

*/
void source_component_common(int num_components,
           cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
           cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
           float *d_freqs, float *d_ls, float *d_ms, float *d_ns,
           float *d_ras, float *d_decs, float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           float *beam_has, float *beam_decs,
           woden_settings_t *woden_settings,
           beam_settings_t beam_settings,
           RTS_MWA_FEE_beam_t *FEE_beam);

/**
@brief Kernel to calculate the visibility response to a number `num_point` of
POINT COMPONENTs, and sum the outputs to `d_sum_visi_*_real`,
`d_sum_visi_*_imag`.

@details Uses the functions `extrap_stokes`, `calc_measurement_equation`,
`update_sum_visis` as detailed above to calculate the visibilities. Sets off
a thread for each visibility to be calculate, with each thread looping over
all COMPONENTs. This seems to keep the memory access low enough to be faster
than having a second dimension over COMPONENT.

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
defined, where:
 - grid.x * threads.x >= `num_visi`

@param[in] *d_point_ras Right Ascensions of all POINTSs (radians)
@param[in] *d_point_decs Declinations of all POINTS (radians)
@param[in] *d_point_freqs Frequencies in the simulation (Hz)
@param[in] *d_point_stokesI Array of reference Stokes I flux densities for all
POINTs (Jy)
@param[in] *d_point_stokesQ Array of reference Stokes Q flux densities for all
POINTs (Jy)
@param[in] *d_point_stokesU Array of reference Stokes U flux densities for all
POINTs (Jy)
@param[in] *d_point_stokesV Array of reference Stokes V flux densities for all
POINTs (Jy)
@param[in] d_point_SIs Array of spectral indicies for all POINTs
@param[in] *d_us Visibility coord \f$u\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_vs Visibility coord \f$v\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_ws Visibility coord \f$w\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in,out] *d_sum_visi_XX_real Pointer to array to sum real XX visibility
into
@param[in,out] *d_sum_visi_XX_imag Pointer to array to sum imaginary XX
visibility into
@param[in,out] *d_sum_visi_XY_real Pointer to array to sum real XY visibility
into
@param[in,out] *d_sum_visi_XY_imag Pointer to array to sum imaginary XY
visibility into
@param[in,out] *d_sum_visi_YX_real Pointer to array to sum real YX visibility
into
@param[in,out] *d_sum_visi_YX_imag Pointer to array to sum imaginary YX
visibility into
@param[in,out] *d_sum_visi_YY_real Pointer to array to sum real YY visibility
into
@param[in,out] *d_sum_visi_YY_imag Pointer to array to sum imaginary YY
visibility into
@param[in] *d_allsteps_wavelengths Wavelength for every baseline, frequency, and
time step in the simulation
@param[in] *d_ls \f$l\f$ sky coord for all time steps for all POINTs
@param[in] *d_ms \f$m\f$ sky coord for all time steps for all POINTs
@param[in] *d_ns \f$n\f$ sky coord for all time steps for all POINTs
@param[in] num_points Number of POINTs
@param[in] num_baselines Number of baselines for one time, one frequency step
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_visis Overall number of visibilities (baselines*freqs*times) in
simulation
@param[in] num_times Number of time steps in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_primay_beam_J00 Array of primary beam J[0,0]
(north-south gain)
@param[in] *d_primay_beam_J01 Array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_primay_beam_J10 Array of primary beam J[1,0]
(east-west leakage)
@param[in] *d_primay_beam_J11 Array of primary beam J[1,1]
(east-west gain)
*/
__global__ void kern_calc_visi_point(float *d_point_ras, float *d_point_decs,
      float *d_point_freqs, float *d_point_stokesI, float *d_point_stokesQ,
      float *d_point_stokesU, float *d_point_stokesV, float *d_point_SIs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_allsteps_wavelengths,
      float *d_ls, float *d_ms, float *d_ns,
      int num_points, int num_baselines, int num_freqs, int num_visis,
      int num_times, int beamtype,
      cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
      cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11);

/**
@brief Kernel to calculate the visibility response to a number `num_gauss` of
GAUSS COMPONENTs, and sum the outputs to `d_sum_visi_*_real`,
`d_sum_visi_*_imag`.

@details Uses the functions `extrap_stokes`, `calc_measurement_equation`,
`update_sum_visis` as detailed above to calculate the visibilities. Furthermore
calculates the visibility envelope \f$\mathrm{V}_{\mathrm{env}}\f$ to convert
the basic visibility into a GAUSSIAN:

\f{eqnarray*}{
\mathrm{V}_{\mathrm{env}}&=& \exp\left( -\frac{\pi^2}{4\ln(2)} \left( k_x^2\theta_{\mathrm{maj}}^2 + k_y^2\theta_{\mathrm{min}}^2\right) \right) \\
k_x &=&  \cos(\phi_{PAj})v + \sin(\phi_{PA})u \\
k_y &=& -\sin(\phi_{PAj})v + \cos(\phi_{PA})u
\f}

where \f$ \theta_{\mathrm{maj}}, \theta_{\mathrm{min}}, \phi_{PA} \f$ are the
major axis, minor axis, and position angle.

Sets off a thread for each visibility to be calculate, with each thread
looping over all COMPONENTs. This seems to keep the memory access low enough
to be faster than having a second dimension over COMPONENT.

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
defined, where:
 - grid.x * threads.x >= `num_visi`

@param[in] *d_gauss_ras Right Ascensions of all GAUSSIANs (radians)
@param[in] *d_gauss_decs Declinations of all GAUSSIANs (radians)
@param[in] *d_gauss_freqs Frequencies in the simulation (Hz)
@param[in] *d_gauss_stokesI Array of reference Stokes I flux densities for all
GAUSSIANs (Jy)
@param[in] *d_gauss_stokesQ Array of reference Stokes Q flux densities for all
GAUSSIANs (Jy)
@param[in] *d_gauss_stokesU Array of reference Stokes U flux densities for all
GAUSSIANs (Jy)
@param[in] *d_gauss_stokesV Array of reference Stokes V flux densities for all
GAUSSIANs (Jy)
@param[in] d_gauss_SIs Array of spectral indicies for all GAUSSIANs
@param[in] *d_us Visibility coord \f$u\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_vs Visibility coord \f$v\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_ws Visibility coord \f$w\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in,out] *d_sum_visi_XX_real Pointer to array to sum real XX visibility
into
@param[in,out] *d_sum_visi_XX_imag Pointer to array to sum imaginary XX
visibility into
@param[in,out] *d_sum_visi_XY_real Pointer to array to sum real XY visibility
into
@param[in,out] *d_sum_visi_XY_imag Pointer to array to sum imaginary XY
visibility into
@param[in,out] *d_sum_visi_YX_real Pointer to array to sum real YX visibility
into
@param[in,out] *d_sum_visi_YX_imag Pointer to array to sum imaginary YX
visibility into
@param[in,out] *d_sum_visi_YY_real Pointer to array to sum real YY visibility
into
@param[in,out] *d_sum_visi_YY_imag Pointer to array to sum imaginary YY
visibility into
@param[in] *d_allsteps_wavelengths Wavelength for every baseline, frequency, and
time step in the simulation
@param[in] *d_ls \f$l\f$ sky coord for all time steps for all GAUSSIANs
@param[in] *d_ms \f$m\f$ sky coord for all time steps for all GAUSSIANs
@param[in] *d_ns \f$n\f$ sky coord for all time steps for all GAUSSIANs
@param[in] *d_gauss_pas Position angles for all GAUSSIANs (radians)
@param[in] *d_gauss_majors Major axis for all GAUSSIANs (FWHM, radians)
@param[in] *d_gauss_minors Minor axis for all GAUSSIANs (FWHM, radians)
@param[in] num_gauss Number of GAUSSIANs components
@param[in] num_baselines Number of baselines for one time, one frequency step
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_visis Overall number of visibilities (baselines*freqs*times) in
simulation
@param[in] num_times Number of time steps in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_primay_beam_J00 Array of primary beam J[0,0]
(north-south gain)
@param[in] *d_primay_beam_J01 Array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_primay_beam_J10 Array of primary beam J[1,0]
(east-west leakage)
@param[in] *d_primay_beam_J11 Array of primary beam J[1,1]
(east-west gain)
*/
__global__ void kern_calc_visi_gaussian(float *d_gauss_ras, float *d_gauss_decs,
      float *d_gauss_freqs, float *d_gauss_stokesI, float *d_gauss_stokesQ,
      float *d_gauss_stokesU, float *d_gauss_stokesV, float *d_gauss_SIs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_allsteps_wavelengths,
      float *d_ls, float *d_ms, float *d_ns,
      float *d_gauss_pas, float *d_gauss_majors, float *d_gauss_minors,
      int num_gauss, int num_baselines, int num_freqs, int num_visis,
      int num_times, int beamtype,
      cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
      cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11);

/**
@brief Kernel to calculate the visibility response to a number `num_shapes` of
SHAPELET COMPONENTs, and sum the outputs to `d_sum_visi_*_real`,
`d_sum_visi_*_imag`.

@details Uses the functions `extrap_stokes`, `calc_measurement_equation`,
`update_sum_visis` as detailed above to calculate the visibilities. Furthermore
calculates the visibility envelope \f$\mathrm{V}_{\mathrm{env}}\f$ to convert
the basic visibility into a SHAPELET:

\f{eqnarray*}{
\mathrm{V}_{\mathrm{env}} &=& \sum^{n_k +n_l < n_\mathrm{max}}_{k,l} C_{n_k,n_l} B_{n_k,n_l}(k_x,k_y)\\
k_x &=& \frac{\pi}{\sqrt{2\ln(2)}} \left[\cos(\phi_{PA})v_{\mathrm{comp}} + \sin(\phi_{PA})u_{\mathrm{comp}} \right] \\
k_y &=& \frac{\pi}{\sqrt{2\ln(2)}} \left[-\sin(\phi_{PA})v_{\mathrm{comp}} + \cos(\phi_{PA})u_{\mathrm{comp}} \right]
\f}

where \f$ B_{n_k,n_l} \f$ is a stored 2D shapelet basis function from a look up
table, of order \f$ n_k,n_l \f$ and scaled by the major and minor axis,
\f$C_{n_k,n_l}\f$ is a scaling coefficient, and \f$ u_{\mathrm{comp}},
v_{\mathrm{comp}} \f$ are visibility coordinates with the SHAPELET ra,dec as the
phase centre. These should have been calculated using
`fundamental_coords.kern_calc_uvw_shapelet`.

This kernel differs from `kern_calc_visi_point` and `kern_calc_visi_gaussian`
in that each SHAPELET component is built up from multiple shapelet basis
functions, and so the kernel (and sometimes the sky model) is split over the
shapelet basis functions, meaning multiple basis function calculations will
use the same COMPONENT \f$l,m,n\f$. The array `d_shape_param_indexes` is
used to match the basis function information `d_shape_n1s`, `d_shape_n2s`,
`d_shape_coeffs` to the COPMONENT information (e.g. `d_shape_ls`).

Sets off a thread for each visibility to be calculate, with each thread
looping over all cofficients. This seems to keep the memory access low enough
to be faster than having a second dimension over coefficient.

`d_sbf` is the shapelet basis function lookup table, and should have been
generated with `shapelet_basis.create_sbf` and copied into device memory.

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
defined, where:
 - grid.x * threads.x >= `num_visi`

@param[in] *d_shape_freqs Frequencies in the simulation (Hz)
@param[in] *d_shape_stokesI Array of reference Stokes I flux densities for all
SHAPELETs (Jy)
@param[in] *d_shape_stokesQ Array of reference Stokes Q flux densities for all
SHAPELETs (Jy)
@param[in] *d_shape_stokesU Array of reference Stokes U flux densities for all
SHAPELETs (Jy)
@param[in] *d_shape_stokesV Array of reference Stokes V flux densities for all
SHAPELETs (Jy)
@param[in] d_shape_SIs Array of spectral indicies for all SHAPELETs
@param[in] *d_us Visibility coord \f$u\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_vs Visibility coord \f$v\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_ws Visibility coord \f$w\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)

@param[in] *d_allsteps_wavelengths Wavelength for every baseline, frequency,
and time step in the simulation
@param[in] *d_u_shapes Array of \f$ u_{\mathrm{comp}} \f$ for every baseline,
frequency, and time step in the simulation (wavelengths)
@param[in] *d_v_shapes Array of \f$ v_{\mathrm{comp}} \f$ for every baseline,
frequency, and time step in the simulation (wavelengths)
@param[in] *d_w_shapes Array of \f$ w_{\mathrm{comp}} \f$ for every baseline,
frequency, and time step in the simulation (wavelengths)

@param[in,out] *d_sum_visi_XX_real Pointer to array to sum real XX visibility
into
@param[in,out] *d_sum_visi_XX_imag Pointer to array to sum imaginary XX
visibility into
@param[in,out] *d_sum_visi_XY_real Pointer to array to sum real XY visibility
into
@param[in,out] *d_sum_visi_XY_imag Pointer to array to sum imaginary XY
visibility into
@param[in,out] *d_sum_visi_YX_real Pointer to array to sum real YX visibility
into
@param[in,out] *d_sum_visi_YX_imag Pointer to array to sum imaginary YX
visibility into
@param[in,out] *d_sum_visi_YY_real Pointer to array to sum real YY visibility
into
@param[in,out] *d_sum_visi_YY_imag Pointer to array to sum imaginary YY
visibility into
time step in the simulation

@param[in] *d_shape_pas Position angles for all SHAPELETs (radians)
@param[in] *d_shape_majors \f$\beta_1 \f$ scaling (major axis) for all SHAPELETs
(radians)
@param[in] *d_shape_minors \f$\beta_2 \f$ scaling (minor axis for all SHAPELETs
(radians)

@param[in] *d_shape_n1s Array of first order of 2D shapelet basis function
@param[in] *d_shape_n2s Array of second order of 2D shapelet basis function
@param[in] *d_shape_coeffs Array of scaling coefficients for basis functions
@param[in] *d_shape_param_indexes Array of index mapping basis functions to
COMPONENT parameters

@param[in] *d_shape_ls \f$l\f$ sky coord for all time steps for all SHAPELETs
@param[in] *d_shape_ms \f$m\f$ sky coord for all time steps for all SHAPELETs
@param[in] *d_shape_ns \f$n\f$ sky coord for all time steps for all SHAPELETs
@param[in] *d_sbf Shapelet basis function lookup table
@param[in] num_shapes Number of SHAPELETs components
@param[in] num_baselines Number of baselines for one time, one frequency step
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_visis Overall number of visibilities (baselines*freqs*times) in
simulation
@param[in] num_coeffs Number of shapelet basis functions and coefficents
@param[in] num_times Number of time steps in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_primay_beam_J00 Array of primary beam J[0,0]
(north-south gain)
@param[in] *d_primay_beam_J01 Array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_primay_beam_J10 Array of primary beam J[1,0]
(east-west leakage)
@param[in] *d_primay_beam_J11 Array of primary beam J[1,1]
(east-west gain)
*/
__global__ void kern_calc_visi_shapelets(float *d_shape_freqs,
      float *d_shape_stokesI, float *d_shape_stokesQ,
      float *d_shape_stokesU, float *d_shape_stokesV, float *d_shape_SIs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_allsteps_wavelengths,
      float *d_u_shapes, float *d_v_shapes, float *d_w_shapes,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_shape_pas, float *d_shape_majors, float *d_shape_minors,
      float *d_shape_n1s, float *d_shape_n2s, float *d_shape_coeffs,
      float *d_shape_param_indexes,
      float *d_shape_ls, float *d_shape_ms, float *d_shape_ns,
      float *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_visis,
      const int num_coeffs, int num_times, int beamtype,
      cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
      cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11);

// __global__ void kern_extrap_stokes(int num_visis, int num_components,
//            float *d_allsteps_wavelengths, float *d_ref_freqs, float *d_SIs,
//            float *d_ref_stokesI, float *d_ref_stokesQ,
//            float *d_ref_stokesU, float *d_ref_stokesV,
//            float *d_flux_I, float *d_flux_Q,
//            float *d_flux_U, float *d_flux_V );

// extern "C" void test_extrap_flux(catsource_t *catsource,
//            const int num_visis, const int num_components,
//            float *allsteps_wavelengths, float *flux_I, float *flux_Q,
//            float *flux_U, float *flux_V );
