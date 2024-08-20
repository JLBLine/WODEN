#pragma once
#include "calculate_visibilities.h"
/*! \file
  Device methods to extrapolate flux densities, assign and apply primary beam
  gains, and visibility responses for different sky model COMPONENT types.
*/

/*!
A struct to contain primary beam values for a give COMPONENT. `d_gxs,d_Dxs,d_Dys,d_gys`
should be used when all antennas have the same primary beam, and `d_gxs,d_Dxs,d_Dys,d_gys` used when all primary beams are different.
*/
typedef struct _d_beam_gains_t {

  gpuUserComplex *d_gxs = NULL; /*!< Device copy of North-South Beam gain values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  gpuUserComplex *d_Dxs = NULL; /*!< Device copy of North-South Beam leakage values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  gpuUserComplex *d_Dys = NULL; /*!< Device copy of East-West Beam leakage values
  for all beams, directions, frequencies, and times for these COMPONENTS*/
  gpuUserComplex *d_gys = NULL; /*!< Device copy of East-West Beam gain values
  for all beams, directions, frequencies, and times for these COMPONENTS*/

  int *d_ant1_to_baseline_map = NULL; /*!< The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1 */
  int *d_ant2_to_baseline_map = NULL; /*!< The index of antenna 2 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 2 */
  int use_twobeams; /*!< The beam gains were made with unique primary beams so
  should use two antenna patterns per visibility */

} d_beam_gains_t;

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
@return `visi`, a `gpuUserComplex` of the visibility
*/
__device__ gpuUserComplex calc_measurement_equation(user_precision_t *d_us,
           user_precision_t *d_vs, user_precision_t *d_ws,
           double *d_ls, double *d_ms, double *d_ns,
           const int iBaseline, const int iComponent);

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
__device__ void apply_beam_gains_stokesIQUV(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY);

/**
@brief Given primary beam gains and leakage terms for antenna 1
`g1x, D1x, D1y, gy` and antenna 2 `g1x, D1x, D1y, gy`,the complex visibility
phase across those two antennas `visi`, and the Stokes I parameter of a source,
simulate the observed XX,XY,YX,YY instrumental cross-correlated visibilities,
where 'x' means aligned north-south, 'y' means east-west.

@details Performs the following calculations:

\f{eqnarray*}{
\mathrm{V}^{XX}_{12} = (g_{1x}g_{2x}^{\ast} + D_{1x}D_{2x}^{\ast})\mathrm{V}^{I}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{XY}_{12} =
     (g_{1x}D_{2y}^{\ast} + D_{1x}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YX}_{12} =
     (D_{1y}g_{2x}^{\ast} + g_{1y}D_{2x}^{\ast})\mathrm{V}^{I}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YY}_{12} =
     (D_{1y}D_{2y}^{\ast} + g_{1y}g_{2y}^{\ast})\mathrm{V}^{I}_{12}
\f}

where \f${\ast}\f$ means complex conjugate, and

\f{eqnarray*}{
\mathrm{V}_{12} &=& \mathrm{V}_{\mathrm{env}}\exp \left( 2\pi i\left( u_{12}l + v_{12}m + w_{12}(n-1) \right) \right) \\
\mathrm{V}^{I}_{12} &=& I(l,m) \mathrm{V}_{12}
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
@param[in] visi_component Complex visibility across antennas 1 and 2
@param[in,out] visi_XX Output XX instrumental visibility
@param[in,out] visi_XY Output XY instrumental visibility
@param[in,out] visi_YX Output YX instrumental visibility
@param[in,out] visi_YY Output YY instrumental visibility
*/
__device__ void apply_beam_gains_stokesI(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY);

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

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *d_Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *d_gys Pointer towards array of primary beam J[1,1]
(east-west gain)
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
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
           gpuUserComplex * g1x, gpuUserComplex * D1x,
           gpuUserComplex * D1y, gpuUserComplex * g1y,
           gpuUserComplex * g2x, gpuUserComplex * D2x,
           gpuUserComplex * D2y, gpuUserComplex * g2y);

/**
@brief Given the type of primary beam simulated `beamtype`, select the beam
gains and leakages from the arrays `d_primay_beam_J*` that match the indexes
`iBaseline`, `iComponent`. This function assumes the primary beams are different
for every antenna, so returns different values for each antenna. 

@todo Currently this is only set up to work with the MWA_FEE and MWA_FEE_INTERP
primary beams

@details This function is built to return the correct beam gain for a given
component on the sky, at a given time, at a given frequency, for a given baseline.
The 4 arrays `d_primay_beam_J00`, `d_primay_beam_J01`, `d_primay_beam_J10`,
`d_primay_beam_J11` should contain the primary beam settings for all antennas,
all times, all frequencies, and COPMONENTs on the sky, and so should be
`num_ants*num_freqs*num_components*num_times` long. The order elements should increment
through COMPONENT (fastest changing), freqeuncy, time, and antenna (slowest changing).
`iBaseline` is a combined index of baseline, frequency, and time.

If `beamtype == NO_BEAM`, set `g1x = g1y = g2x = g2y = 1` and
`D1x = D1y = D2x = D2y = 0`

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *d_Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *d_gys Pointer towards array of primary beam J[1,1]
(east-west gain)
@param[in] *d_ant1_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1
@param[in] *d_ant2_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 2
@param[in,out] *g1x Beam gain antenna 1 in north-south
@param[in,out] *D1x Beam leakage antenna 1 from north-south
@param[in,out] *D1y Beam gain antenna 1 in east-west
@param[in,out] *g1y Beam leakage antenna 1 from east-west
@param[in,out] *g2x Beam gain antenna 2 in north-south
@param[in,out] *D2x Beam leakage antenna 2 from north-south
@param[in,out] *D2y Beam gain antenna 2 in east-west
@param[in,out] *g2y Beam leakage antenna 2 from east-west

*/
__device__ void get_beam_gains_multibeams(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
           int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map,
           gpuUserComplex * g1x, gpuUserComplex * D1x,
           gpuUserComplex * D1y, gpuUserComplex * g1y,
           gpuUserComplex * g2x, gpuUserComplex * D2x,
           gpuUserComplex * D2y, gpuUserComplex * g2y);

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
@param[in] *d_gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *d_Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *d_gys Pointer towards array of primary beam J[1,1]
(east-west gain)
@param[in] *d_ant1_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1 (only used when use_twobeams == 1)
@param[in] *d_ant2_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 2 (only used when use_twobeams == 1)
@param[in] use_twobeams If 1 (True), assume all primary beams are different for each
antenna. If 0 (False), assume all primary beams are the same for all antennas.
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
__device__ void update_sum_visis_stokesIQUV(int iBaseline, int iComponent, int num_freqs,
    int num_baselines, int num_components, int num_times, int beamtype,
    gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
    gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
    int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
    gpuUserComplex visi_component,
    user_precision_t flux_I, user_precision_t flux_Q,
    user_precision_t flux_U, user_precision_t flux_V,
    user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
    user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
    user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
    user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag);

/**
@brief Given the visibility between two recieving elements for a COMPONENT
`visi_component`, at the baseline/freq/time index given by `iBaseline` and
COMPONENT index `iComponent`, apply the instrumental primary beam and
COMPONENT Stokes I parameter to create instrumental XX,XY,YX,YY visibilities,
and sum them into real and imaginary XX,XY,YX,YY visibilities arrays
`d_sim_visi_*_real` and `d_sim_visi_*_imag`.

@details Uses `get_beam_gains` or `get_beam_gains_multibeams`,
 and `apply_beam_gains` as described above to apply the gains -
 see descriptions for what should be the arguments to them.

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] *d_gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *d_Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *d_Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *d_gys Pointer towards array of primary beam J[1,1]
(east-west gain)
@param[in] *d_ant1_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1 (only used when use_twobeams == 1)
@param[in] *d_ant2_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 2 (only used when use_twobeams == 1)
@param[in] use_twobeams If 1 (True), assume all primary beams are different for each
antenna. If 0 (False), assume all primary beams are the same for all antennas.
@param[in] visi_component Complex visibility across antennas 1 and 2
@param[in] flux_I Stokes I flux density (Jy)
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
__device__ void update_sum_visis_stokesI(int iBaseline, int iComponent, int num_freqs,
    int num_baselines, int num_components, int num_times, int beamtype,
    gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
    gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
    int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
    gpuUserComplex visi_component,
    user_precision_t flux_I,
    user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
    user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
    user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
    user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag);



/**
@brief Allocate device memory of extrapolated Stokes arrays int `d_components`

@details It does this if d_components.do_QUV == 1:

    d_components->extrap_stokesI = NULL;
    gpuMalloc( (void**)&d_components->extrap_stokesI, num_comps*num_freqs*sizeof(double) );
    d_components->extrap_stokesQ = NULL;
    gpuMalloc( (void**)&d_components->extrap_stokesQ, num_comps*num_freqs*sizeof(double) );
    d_components->extrap_stokesU = NULL;
    gpuMalloc( (void**)&d_components->extrap_stokesU, num_comps*num_freqs*sizeof(double) );
    d_components->extrap_stokesV = NULL;
    gpuMalloc( (void**)&d_components->extrap_stokesV, num_comps*num_freqs*sizeof(double) );

If do_QUV == 0, only allocate the StokesI array.

@param[in,out] *d_components A populated `components_t` struct
@param[in] num_comps Number of components
@param[in] num_freqs Number of frequencies
*/
void malloc_extrapolated_flux_arrays(components_t *d_components, int num_comps,
                                     int num_freqs);


/**
@brief Assuming a simple power-law SED of \f$ S \propto \nu^\alpha \f$,
extrapolate flux density parameters to the requested frequencies

@details Assumes a referece frequency of 200MHz, and uses the
reference flux density `d_ref_fluxes[iFluxComp]`, and `spectral index d_SIs[iFluxComp]` to calculate the flux density `extrap_flux` at the requested frequency d_extrap_freqs[iFreq].

@param[in] d_ref_fluxes Reference fluxes for all components
@param[in] d_SIs Spectral indexes for all components
@param[in] d_extrap_freqs Pointer to array of frequencies to extrapolate to
@param[in] iFluxComp Index of which component to extrapolate
@param[in] iFreq Index of which frequency to extrapolate to
@param[in,out] *extrap_flux Pointer to extrapolated flux (Jy)
*/
__device__ void extrap_stokes_power_law(user_precision_t *d_ref_fluxes,
           user_precision_t *d_SIs,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux);

/**
@brief Kernel to run `extrap_stokes_power_law` for all Stokes I components in `d_components`.

@details Fills the array `d_components.extrap_stokesI` with the extrapolated Stokes I flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);

/**
@brief Kernel to run `extrap_stokes_power_law` for all Stokes V components in `d_components`.

@details Fills the array `d_components.extrap_stokesV` with the extrapolated Stokes V flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_power_laws_stokesV(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);

/**
@brief Kernel to run `extrap_stokes_power_law` for all linear polarisation components in `d_components`.

@details Temporarily fills the array `d_components.extrap_stokesQ` with the extrapolated linear polarisation flux densities. The intention is to the use `kern_apply_rotation_measure` after the fact, which will use `d_components.extrap_stokesQ` as input flux, do rotations, and then fill `d_components.extrap_stokesU` and `d_components.extrap_stokesV`.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_power_laws_linpol(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);


/**
@brief Assuming a curved power-law SED of \f$ S \propto \nu^\alpha \exp(q \ln (\nu)^2 )  \f$,
extrapolate Stokes I flux density parameters to the requested frequencies

@details Assumes a referece frequency of 200MHz, and uses the
reference flux density `d_ref_fluxes[iFluxComp]`, `spectral index d_SIs[iFluxComp]`,
curvature `q` param `spectral index d_qs[iFluxComp]` to calculate the flux density `extrap_flux` at the requested frequency d_extrap_freqs[iFreq].

@param[in] d_ref_fluxes Reference fluxes for all components
@param[in] d_SIs Spectral indexes for all components
@param[in] d_qs Curvature `q` param for all components
@param[in] d_extrap_freqs Pointer to array of frequencies to extrapolate to
@param[in] iFluxComp Index of which component to extrapolate
@param[in] iFreq Index of which frequency to extrapolate to
@param[in,out] *extrap_flux Pointer to extrapolated flux (Jy)
*/
__device__ void extrap_stokes_curved_power_law(user_precision_t *d_ref_fluxes,
           user_precision_t *d_SIs, user_precision_t *d_qs,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux);

/**
@brief Kernel to run `extrap_stokes_curved_power_law` for all Stokes I components in `d_components`.

@details Fills the array `d_components.extrap_stokesI` with the extrapolated Stokes I flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_curved_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);


/**
@brief Kernel to run `extrap_stokes_curved_power_law` for all Stokes I components in `d_components`.

@details Fills the array `d_components.extrap_stokesV` with the extrapolated Stokes V flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_curved_power_laws_stokesV(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);

/**
@brief Kernel to run `extrap_stokes_curved_power_law` for all linear polarisation components in `d_components`.

@details Temporarily fills the array `d_components.extrap_stokesQ` with the extrapolated linear polarisation flux densities. The intention is to the use `kern_apply_rotation_measure` after the fact, which will use `d_components.extrap_stokesQ` as input flux, do rotations, and then fill `d_components.extrap_stokesU` and `d_components.extrap_stokesV`.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_curved_power_laws_linpol(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);

/**
@brief Kernel to calculate polarisation fractions for all Stokes V components in `d_components`.

@details MUST have the `d_components.extrap_stokesI` filled already. Fills the array `d_components.extrap_stokesV` using `d_components.extrap_stokesI` and `d_components.stokesV_pol_fracs`.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_polarisation_fraction_stokesV(int num_extrap_freqs, 
             double *d_extrap_freqs, int num_comps, components_t d_components);

/**
@brief Kernel to calculate polarisation fractions for all linear polarisation components in `d_components`.

@details MUST have the `d_components.extrap_stokesI` filled already. Temporarily fills the array `d_components.extrap_stokesQ` using `d_components.extrap_stokesI` and `d_components.linpol_pol_fracs`. The intention is to the use `kern_apply_rotation_measure` after the fact, which will use `d_components.extrap_stokesQ` as input flux, do rotations, and then fill `d_components.extrap_stokesU` and `d_components.extrap_stokesV`.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_polarisation_fraction_linpol(int num_extrap_freqs, 
             double *d_extrap_freqs, int num_comps, components_t d_components);



/**
Assuming a list-type spectral model, extrapolates the flux for a given component and
frequency.

@param[in] user_precision_t *list_stokes Array containing all the list-fluxes used for the extrapolation.
@param[in] double *list_freqs Array containing all the list-frequencies used for the extrapolation; must match order of `list_stokes`.
@param[in] int *num_list_values Array containing the number of list values for each component.
@param[in] int *list_start_indexes Array containing the starting index of each component within `list_stokes` and `list_freqs`.
@param[in] d_extrap_freqs The frequencies to use for the extrapolation.
@param[in] iFluxComp The index of the flux component to store the result in.
@param[in] iFreq The index of the frequency component to extrapolate.
@param[in,out]  extrap_flux The extrapolated flux.
 */
__device__ void extrap_stokes_list_fluxes(user_precision_t *list_stokes,
           double *list_freqs, int *arr_num_list_values, int *list_start_indexes,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux);

/**
@brief Extrapolates list-type spectral model fluxes to given frequencies.

@details Fills `extrap_stokes` with the extrapolated stokes flux densities. Runs the function `extrap_stokes_list_fluxes`. 

@param[in] user_precision_t *list_stokes Array containing all the list-fluxes used for the extrapolation.
@param[in] double *list_freqs Array containing all the list-frequencies used for the extrapolation; must match order of `list_stokes`.
@param[in] int *num_list_values Array containing the number of list values for each component.
@param[in] int *list_start_indexes Array containing the starting index of each component within `list_stokes` and `list_freqs`.
@param[in] int *list_comp_inds Array containing the component index for list component; used to index the extrapolated fluxes to `extrap_stokes`.
@param[in] num_extrap_freqs The number of extrapolation frequencies.
@param[in] d_extrap_freqs Pointer to the array of extrapolation frequencies.
@param[in] num_comps The number of components.
@param[in,out] extrap_stokes Output extrapolated fluxes
 */
__global__ void kern_extrap_list_fluxes(user_precision_t *list_stokes, double *list_freqs,
                                        int *num_list_values, int *list_start_indexes,
                                        int *list_comp_inds,
                                        int num_extrap_freqs, double *d_extrap_freqs,
                                        int num_comps, user_precision_t *extrap_stokes);



/**
@brief Extrapolate the flux densities of certain spectral model type in a source to a set of frequencies.

@details This function is a wrapper for the various extrapolation functions.
It will extrapolate all spectral model types for a give component type as specified
by `comptype`. It will extrapolate to all the frequencies specified in `d_extrap_freqs`,
and store the results within `d_chunked_source`


@param[in,out] d_chunked_source Pointer to the source data.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to.
@param[in] num_extrap_freqs Number of frequencies in the `d_extrap_freqs` array.
@param[in] comptype The type of component to extrapolate (e.g. POINT, GAUSSIAN, SHAPELET).
 */
extern "C" void extrapolate_Stokes(source_t *d_chunked_source,
                                   double *d_extrap_freqs, int num_extrap_freqs,
                                   e_component_type comptype);

/**
@brief Performs necessary calculations that are common to all POINT, GAUSSIAN,
and SHAPELET component types: calculating the \f$l,m,n\f$ coordinates;
extrapolating flux densities to given frequencies; calculating
primary beam values.

@details *chunked_source should be a populated `source_t` struct as filled by
`chunk_sky_model::create_chunked_sky_models`, and then remapped by
`chunk_sky_model::remap_source_for_gpu`. The remapping has a memory cost,
hence should be done iteratively over chunks (this is done inside
`calculate_visibilities.cu`).

This function uses `fundamental_coords::kern_calc_lmn` to calculate
\f$l,m,n\f$, and stores the ouputs
in the `components_t` struct `d_components`. Uses the `d_beam_gains_t` struct
`d_component_beam_gains` to store the calculated primary beam responses. Each
gain and leakage will have
`num_components*woden_settings->num_time_steps*woden_settings->num_freqs`
elements of memory allocated in the process.

If the primary beam requested (set via beam_settings->beamtype ) is either
GAUSS_BEAM or MWA_ANALY, both `chunked_source->components->beam_has` and
`chunked_source->components->beam_decs`
need to be set.

If `woden_settings->do_autos` is True, will use `source_components::kern_calc_autos`
to calculate all the auto-correlations, and stores them in `d_visibility_set`.


@param[in] woden_settings Populated `woden_settings_t` struct
@param[in] beam_settings Populated `beam_settings_t` struct
@param[in] d_freqs Frequencies to calculate beam responses to (stored on the device)
@param[in] *chunked_source A populated `source_t` struct as filled by
`chunk_sky_model::create_chunked_sky_models` and then remapped by
`chunk_sky_model::remap_source_for_gpu`
@param[in, out] *d_chunked_source A populated `source_t` struct with device memory
as filled by `source_components::copy_chunked_source_to_GPU`
@param[in,out] d_component_beam_gains Pointer to `d_beam_gains_t` struct in
which to malloc and store results on the device in
@param[in] comptype Either POINT, GAUSSIAN, SHAPELET - selects which set of
COMPONENTS stored in `chunked_source` and `chunked_source` to work with
@param[in] d_visibility_set If `woden_settings->do_autos` is True, use this to store auto
correlations in. Needs to have device memory allocated. Otherwise, `d_visibility_set`
is unused, so can just be an empty pointer
*/
extern "C" void source_component_common(woden_settings_t *woden_settings,
           beam_settings_t *beam_settings, double *d_freqs,
           source_t *chunked_source, source_t *d_chunked_source,
           d_beam_gains_t *d_component_beam_gains,
           e_component_type comptype,
           visibility_set_t *d_visibility_set);

/**
@brief Kernel to calculate the visibility response to a number `num_components`
of either POINT or GAUSSIAN COMPONENTs, and sum the outputs to `d_sum_visi_*_real`,
`d_sum_visi_*_imag`.

@details Uses the functions `extrap_stokes`, `calc_measurement_equation`,
`update_sum_visis` as detailed above to calculate the visibilities. Sets off
a thread for each visibility to be calculated, with each thread looping over
all COMPONENTs. This seems to keep the memory access low enough to be faster
than having a second dimension over COMPONENT. Specify if you want POINT or
GAUSSIAN using `comptype`.

The code is almost exactly the same for POINT or GAUSSIANs, except when
running with a Gaussian, a visiblity envelope is calculated as

\f{eqnarray*}{
\mathrm{V}_{\mathrm{env}}&=& \exp\left( -\frac{\pi^2}{4\ln(2)} \left( k_x^2\theta_{\mathrm{maj}}^2 + k_y^2\theta_{\mathrm{min}}^2\right) \right) \\
k_x &=&  \cos(\phi_{PAj})v + \sin(\phi_{PAj})u \\
k_y &=& -\sin(\phi_{PAj})v + \cos(\phi_{PAj})u
\f}

where \f$ \theta_{\mathrm{maj}}, \theta_{\mathrm{min}}, \phi_{PA} \f$ are the
major axis, minor axis, and position angle. The visiblity equation is multipled
by this envelope to modify from a POINT source into a GAUSSIAN

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
defined, where:
 - grid.x * threads.x >= `num_visi`

@param[in] d_components d_components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] d_component_beam_gains Pointer to a populated `d_beam_gains_t` struct as filled by `source_components::source_component_common`
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
@param[in] num_components Either number of POINT of GAUSSIANS components
@param[in] num_baselines Number of baselines for one time, one frequency step
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_cross Overall number of cross-correlations  (baselines*freqs*times) in
simulation
@param[in] num_times Number of time steps in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] comptype Component type, either POINT or GAUSSIAN
*/
__global__ void kern_calc_visi_point_or_gauss(components_t d_components,
           d_beam_gains_t d_component_beam_gains,
           user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
           user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
           user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
           user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
           user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
           int num_components, int num_baselines, int num_freqs, int num_cross,
           int num_times, e_beamtype beamtype, e_component_type comptype);

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
`fundamental_coords::kern_calc_uvw_shapelet`.

This kernel differs from `kern_calc_visi_point_or_gauss`
in that each SHAPELET component is built up from multiple shapelet basis
functions, and so the kernel (and sometimes the sky model) is split over the
shapelet basis functions, meaning multiple basis function calculations will
use the same COMPONENT \f$l,m,n\f$. The array `d_components.param_indexes` is
used to match the basis function information `d_components.n1s`, `d_components.n2s`,
`d_components.shape_coeffs` to the COPMONENT information (e.g. `d_components.ls`).

A further difference is that the shapelet \f$ u_{\mathrm{comp}},
v_{\mathrm{comp}} \f$ used in the visibility envelope equation should
be stored in metres only - for every baseline (fastest
changing), every time step, and every SHAPELET component (slowest changing).
They will be scaled by wavelength inside the kernel, as memory costs
become prohibitive to store for every wavelength as well. The visibility \f$u,v,w\f$
are stored for every baseline, wavelength, and time step.

Sets off a thread for each visibility to be calculate, with each thread
looping over all cofficients. This seems to keep the memory access low enough
to be faster than having a second dimension over coefficient.

`d_sbf` is the shapelet basis function lookup table, and should have been
generated with `shapelet_basis.create_sbf` and copied into device memory.

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
defined, where:
 - grid.x * threads.x >= `num_visi`

@param[in] d_components d_components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] d_component_beam_gains Pointer to a populated `d_beam_gains_t` struct as filled by `source_components::source_component_common`
@param[in] *d_us Visibility coord \f$u\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_vs Visibility coord \f$v\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_ws Visibility coord \f$w\f$ for every baseline, frequency, and time
step in the simulation (wavelengths)
@param[in] *d_allsteps_wavelengths Wavelength for every baseline, frequency,
and time step in the simulation
@param[in] *d_u_shapes Array of \f$ u_{\mathrm{comp}} \f$ for every baseline,
frequency, and SHAPELET component in the simulation (metres)
@param[in] *d_v_shapes Array of \f$ v_{\mathrm{comp}} \f$ for every baseline,
frequency, and SHAPELET component in the simulation (metres)
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
@param[in] *d_sbf Shapelet basis function lookup table
@param[in] num_shapes Number of SHAPELETs components
@param[in] num_baselines Number of baselines for one time, one frequency step
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_cross Overall number of cross-correlations (baselines*freqs*times) in
simulation
@param[in] num_coeffs Number of shapelet basis functions and coefficents
@param[in] num_times Number of time steps in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
*/
__global__ void kern_calc_visi_shapelets(components_t d_components,
      d_beam_gains_t d_component_beam_gains,
      user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
      user_precision_t *d_allsteps_wavelengths,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
      user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
      user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
      user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
      user_precision_t *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_cross,
      const int num_coeffs, int num_times, e_beamtype beamtype);

/**
@brief Copies the specified type of source components from host memory to device memory.

@details This function copies the specified type of source components from the host memory pointed to by 
 `chunked_source` to the device memory pointed to by `d_chunked_source`, allocating the device
 memory in the process The `comptype` parameter specifies the type of components to copy. 
 
@param[in] chunked_source Pointer to the host memory containing the source components to copy.
@param[in,out] d_chunked_source Pointer to the device memory where the source components will be copied.
@param[in] comptype The type of source components to copy.
 */
void copy_components_to_GPU(source_t *chunked_source, source_t *d_chunked_source,
                            e_component_type comptype);

/**
@brief Frees device memory associated with `d_chunked_source`, depending on
what type of COMPONENT you are freeing

@details Either frees all device memory that will have been allocated from
either `d_chunked_source->point_components`, `d_chunked_source->gauss_components`,
or `d_chunked_source->shape_components`, depending on `comptype`.

@param[in,out] *d_chunked_source A `source_t` with device memory to be
freed
@param[in] comptype Which type of COMPONENT, either POINT, GAUSSIAN, or
 SHAPELET

*/
extern "C" void free_d_components(source_t *d_chunked_source,
                                  e_component_type comptype);


/**
@brief Frees the device memory on `d_beam_gains`, given the type of
primary beam `beamtype`.

@details Only certain primary beam models have leakage terms, so need
to know the `beamtype` to free the correct locations.

@param[in,out] d_beam_gains A `d_beam_gains_t` with device memory to be
freed
@param[in] beamtype An `e_beamtype` describing the primary beam type

*/
extern "C" void free_beam_gains(d_beam_gains_t d_beam_gains, e_beamtype beamtype);


/**
@brief Copies a populated `source_t` from host to device.

@details Depending on what is in chunked_source (e.g. POINT, GAUSSIAN, SHAPELET,
POWER_LAW, CURVED_POWER_LAW, LIST) type information, only what needs to be
copied across to the GPU is (no empty pointers arrays are copied across)

@param[in,out] *chunked_source A populated `source_t` struct
*/
source_t * copy_chunked_source_to_GPU(source_t *chunked_source);

/**
@brief Free device memory of extrapolated Stokes arrays from `d_components`

@details It does this if d_components.do_QUV == 1:

   gpuFree d_components->extrap_stokesI );
   gpuFree d_components->extrap_stokesQ );
   gpuFree d_components->extrap_stokesU );
   gpuFree d_components->extrap_stokesV );

If do_QUV == 0, only free the StokesI array.

@param[in,out] *d_components A populated `components_t` struct
*/
void free_extrapolated_flux_arrays(components_t *d_components);

/**
@brief Calculate the auto-correlations for all antennas given the fluxes
already calculated in `d_components` and beam gains in `d_component_beam_gains`.
Stores the outputs at the end of `d_sum_visi*`, after the cross-correlations.

@details Currently, the primary beam for all antennas (or what the MWA calls
tiles) are identical. So `d_component_beam_gains` should have gains for
one primary beam. This function then just does the dot product of the
component fluxes and beam gains to produce linear stokes polarisation
auto-correlations.

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
and grid.y defined, where:
 - grid.x * threads.x >= `num_freqs*num_times`
 - grid.y * threads.y >= `num_ants`

@param[in] d_components d_components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] d_component_beam_gains Pointer to a populated `d_beam_gains_t` struct as filled by `source_components::source_component_common`
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] num_components Number of components in `d_components`
@param[in] num_baselines Number of baselines for one time, one frequency step
@param[in] num_freqs Number of frequncies in the simulation
@param[in] num_times Number of times in the simulation
@param[in] num_ants Number of antennas in the array
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
@param[in] use_twobeams If True, use a two primary beams per visibility.
Otherwise, assume all primary beams are identical
@param[in] *d_ant1_to_auto_map An index of all primary beams to auto-correlations
Currently this is just an index of all antennas. Gets passed to `get_beam_gains_multibeams`
@param[in] *d_ant2_to_auto_map An index of all primary beams to auto-correlations
Currently this is just an index of all antennas. Gets passed to `get_beam_gains_multibeams`

*/
__global__ void kern_calc_autos(components_t d_components,
                                d_beam_gains_t d_component_beam_gains,
                                int beamtype,
                                int num_components, int num_baselines,
                                int num_freqs, int num_times, int num_ants,
                                user_precision_t *d_sum_visi_XX_real,
                                user_precision_t *d_sum_visi_XX_imag,
                                user_precision_t *d_sum_visi_XY_real,
                                user_precision_t *d_sum_visi_XY_imag,
                                user_precision_t *d_sum_visi_YX_real,
                                user_precision_t *d_sum_visi_YX_imag,
                                user_precision_t *d_sum_visi_YY_real,
                                user_precision_t *d_sum_visi_YY_imag,
                                int use_twobeams,
                                int *d_ant1_to_auto_map,
                                int *d_ant2_to_auto_map);

/**
@brief Fill the `d_ant1_to_baseline_map` and `d_ant2_to_baseline_map` device arrays
with indexes corresponding to ant1 and ant2 for all unique baselines in an
array of `num_ants` antennas.

@details The `d_ant1_to_baseline_map` and `d_ant2_to_baseline_map` should
already have their memory allocated

@param[in] num_ants Number of antennas in the array
@param[in,out] *d_ant1_to_baseline_map Device memory-allocated array of size `((num_ants - 1)*num_ants) / 2`
@param[in,out] *d_ant2_to_baseline_map Device memory-allocated array of size `((num_ants - 1)*num_ants) / 2`

*/
extern "C" void fill_ant_to_baseline_mapping(int num_ants, int *d_ant1_to_baseline_map,
                                               int *d_ant2_to_baseline_map);