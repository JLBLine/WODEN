/*! \file
  GPU methods to extrapolate flux densities, assign and apply primary beam
  gains, and visibility responses for different sky model COMPONENT types.
*/
#pragma once
#include "calculate_visibilities_gpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "gpucomplex.h"
#include "fundamental_coords_gpu.h"
#include "constants.h"
#include "source_components_gpu.h"
#include "woden_struct_defs.h"
#include "primary_beam_gpu.h"
#include "woden_precision_defs.h"
#include "gpu_macros.h"
#include "source_components_common.h"


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
__device__ gpuUserComplex calc_measurement_equation_gpu(user_precision_t *d_us,
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
__device__ void apply_beam_gains_stokesIQUV_on_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
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
__device__ void apply_beam_gains_stokesI_on_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY);

/**
@brief Given primary beam gains and leakage terms for antenna 1
`g1x, D1x, D1y, gy` and antenna 2 `g1x, D1x, D1y, gy`,the complex visibility
phase across those two antennas `visi`, and the Stokes parameters of a source,
simulate the observed XX,XY,YX,YY instrumental cross-correlated visibilities,
where 'x' means aligned north-east south-west (45 deg),
'y' means south-east north-west (135 deg).

@details Performs the following calculations:

\f{eqnarray*}{
\mathrm{V}^{XX}_{12} = (g_{1x} g_{2x}^{\ast} + D_{1x} D_{2x}^{\ast})\mathrm{V}^{I}_{12}
                     - (g_{1x} D_{2x}^{\ast} + D_{1x} g_{2x}^{\ast})\mathrm{V}^{Q}_{12} \\
                     + (g_{1x} g_{2x}^{\ast} - D_{1x} D_{2x}^{\ast})\mathrm{V}^{U}_{12}
                     + i(g_{1x} D_{2x}^{\ast} - D_{1x} g_{2x}^{\ast})\mathrm{V}^{V}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{XY}_{12} =
      (g_{1x} D_{2y}^{\ast} + D_{1x} g_{2y}^{\ast})\mathrm{V}^{I}_{12}
  -  (g_{1x} g_{2y}^{\ast} + D_{1x} D_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
  +   (g_{1x} D_{2y}^{\ast} - D_{1x} g_{2y}^{\ast})\mathrm{V}^{U}_{12}
  +  i(g_{1x} g_{2y}^{\ast} -D_{1x} D_{2y}^{\ast})\mathrm{V}^{V}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YX}_{12} =
     (D_{1y} g_{2x}^{\ast} + G_{1y} D_{2x}^{\ast})\mathrm{V}^{I}_{12}
  -  (D_{1y} D_{2x}^{\ast} + g_{1y} g_{2x}^{\ast})\mathrm{V}^{Q}_{12} \\
  +  (D_{1y} g_{2x}^{\ast} - G_{1y} D_{2x}^{\ast})\mathrm{V}^{U}_{12}
  +  i(D_{1y} D_{2x}^{\ast} -g_{1y} g_{2x}^{\ast})\mathrm{V}^{V}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YY}_{12} =
     (D_{1y} D_{2y}^{\ast} + G_{1y} g_{2y}^{\ast})\mathrm{V}^{I}_{12}
  -  (D_{1y} g_{2y}^{\ast} + g_{1y} D_{2y}^{\ast})\mathrm{V}^{Q}_{12} \\
  +  (D_{1y} D_{2y}^{\ast} - G_{1y} g_{2y}^{\ast})\mathrm{V}^{U}_{12}
  +  i(D_{1y} g_{2y}^{\ast} -g_{1y} D_{2y}^{\ast})\mathrm{V}^{V}_{12}
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

@param[in] g1x Beam gain antenna 1 in north-east south-west (45 deg)
@param[in] D1x Beam leakage antenna 1 from north-east south-west (45 deg)
@param[in] D1y Beam gain antenna 1 in south-east north-west (135 deg)
@param[in] g1y Beam leakage antenna 1 from south-east north-west (135 deg)
@param[in] g2x Beam gain antenna 2 in north-east south-west (45 deg)
@param[in] D2x Beam leakage antenna 2 from north-east south-west (45 deg)
@param[in] D2y Beam gain antenna 2 in south-east north-west (135 deg)
@param[in] g2y Beam leakage antenna 2 from south-east north-west (135 deg)
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
__device__ void apply_beam_gains_stokesIQUV_off_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
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
where 'x' means aligned north-east south-west (45 deg),
'y' means south-east north-west (135 deg).

@details Performs the following calculations:

\f{eqnarray*}{
\mathrm{V}^{XX}_{12} = (g_{1x} g_{2x}^{\ast} + D_{1x} D_{2x}^{\ast})\mathrm{V}^{I}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{XY}_{12} =
      (g_{1x} D_{2y}^{\ast} + D_{1x} g_{2y}^{\ast})\mathrm{V}^{I}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YX}_{12} =
     (D_{1y} g_{2x}^{\ast} + G_{1y} D_{2x}^{\ast})\mathrm{V}^{I}_{12}
\f}

\f{eqnarray*}{
\mathrm{V}^{YY}_{12} =
     (D_{1y} D_{2y}^{\ast} + G_{1y} g_{2y}^{\ast})\mathrm{V}^{I}_{12}
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

@param[in] g1x Beam gain antenna 1 in north-east south-west (45 deg)
@param[in] D1x Beam leakage antenna 1 from north-east south-west (45 deg)
@param[in] D1y Beam gain antenna 1 in south-east north-west (135 deg)
@param[in] g1y Beam leakage antenna 1 from south-east north-west (135 deg)
@param[in] g2x Beam gain antenna 2 in north-east south-west (45 deg)
@param[in] D2x Beam leakage antenna 2 from north-east south-west (45 deg)
@param[in] D2y Beam gain antenna 2 in south-east north-west (135 deg)
@param[in] g2y Beam leakage antenna 2 from south-east north-west (135 deg)
@param[in] flux_I Stokes I flux density (Jy)
@param[in] visi_component Complex visibility across antennas 1 and 2
@param[in,out] visi_XX Output XX instrumental visibility
@param[in,out] visi_XY Output XY instrumental visibility
@param[in,out] visi_YX Output YX instrumental visibility
@param[in,out] visi_YY Output YY instrumental visibility
*/
__device__ void apply_beam_gains_stokesI_off_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY);

/**
@brief Given the type of primary beam simulated `beamtype`, select the beam
gains and leakages from the arrays `d_g*, d_D*` that match the indexes
`iBaseline`, `iComponent`. NOTE that currently, primary beams are assumed
to be the same for all recieving elements, and so this function only takes in
one set of primary beam arrays.

@details This function is built to return the correct beam gain for a given
component on the sky, at a given time, at a given frequency, for a given baseline.
The 4 arrays `d_gxs`, `d_Dxs`, `d_Dys`,`d_gys` should contain the primary beam
values for all times, all frequencies, and COPMONENTs on the sky, and so should be
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
__device__ void get_beam_gains_gpu(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
           gpuUserComplex * g1x, gpuUserComplex * D1x,
           gpuUserComplex * D1y, gpuUserComplex * g1y,
           gpuUserComplex * g2x, gpuUserComplex * D2x,
           gpuUserComplex * D2y, gpuUserComplex * g2y);

/**
@brief Given the type of primary beam simulated `beamtype`, select the beam
gains and leakages from the arrays `d_g*, d_D*` that match the indexes
`iBaseline`, `iComponent`. This function assumes the primary beams are different
for every antenna, so returns different values for each antenna. 

@details This function is built to return the correct beam gain for a given
component on the sky, at a given time, at a given frequency, for a given baseline.
The 4 arrays  `d_gxs`, `d_Dxs`, `d_Dys`,`d_gys`
should contain the primary beam settings for all antennas,
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
__device__ void get_beam_gains_multibeams_gpu(int iBaseline, int iComponent, int num_freqs,
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

@details Uses `get_beam_gains_gpu` and `apply_beam_gains` as described above to
apply the gains - see descriptions for what should be the arguments to them.

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] off_cardinal Boolean to indicate if the dipoles are off-cardinal i.e.
if off_cardinal == 1, then the dipoles are aligned at 45 and 135 degrees,
if off_cardinal == 0, then the dipoles are aligned at 0 and 90 degrees
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
__device__ void update_sum_visis_stokesIQUV_gpu(int iBaseline, int iComponent,
    int num_freqs, int num_baselines, int num_components, int num_times,
    int beamtype, int off_cardinal,
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

@details Uses `get_beam_gains_gpu` or `get_beam_gains_multibeams_gpu`,
 and `apply_beam_gains` as described above to apply the gains -
 see descriptions for what should be the arguments to them.

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] off_cardinal Boolean to indicate if the dipoles are off-cardinal i.e.
if off_cardinal == 1, then the dipoles are aligned at 45 and 135 degrees,
if off_cardinal == 0, then the dipoles are aligned at 0 and 90 degrees
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
__device__ void update_sum_visis_stokesI_gpu(int iBaseline, int iComponent,
    int num_freqs, int num_baselines, int num_components, int num_times,
    int beamtype, int off_cardinal,
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
@brief Kernel to set all elements of `array` to zero.

@details Call using something like
     - dim3 grid, threads;
     - threads.x = 128;
     - grid.x = (int)ceil( (float)num_arr / (float)threads.x );

@param[in,out] *array Array to set to zeros
@param[in] num_arr Number of elements in the array
*/
__global__ void kern_make_zeros_user_precision(user_precision_t *array, int num_arr);

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
extern "C" void malloc_extrapolated_flux_arrays_gpu(components_t *d_components, int num_comps,
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
__device__ void extrap_stokes_power_law_gpu(user_precision_t *d_ref_fluxes,
           user_precision_t *d_SIs,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux);

/**
@brief Kernel to run `extrap_stokes_power_law_gpu` for all Stokes I components in `d_components`.

@details Fills the array `d_components.extrap_stokesI` with the extrapolated Stokes I flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);

/**
@brief Kernel to run `extrap_stokes_power_law_gpu` for all Stokes V components in `d_components`.

@details Fills the array `d_components.extrap_stokesV` with the extrapolated Stokes V flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_power_laws_stokesV(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);

/**
@brief Kernel to run `extrap_stokes_power_law_gpu` for all linear polarisation components in `d_components`.

@details Temporarily fills the array `d_components.extrap_stokesQ` with the extrapolated linear polarisation flux densities. The intention is to the use `kern_apply_rotation_measure` after the fact, which will use `d_components.extrap_stokesQ` as input flux, do rotations, and then fill `d_components.extrap_stokesQ` and `d_components.extrap_stokesU`.

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
__device__ void extrap_stokes_curved_power_law_gpu(user_precision_t *d_ref_fluxes,
           user_precision_t *d_SIs, user_precision_t *d_qs,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux);

/**
@brief Kernel to run `extrap_stokes_curved_power_law_gpu` for all Stokes I components in `d_components`.

@details Fills the array `d_components.extrap_stokesI` with the extrapolated Stokes I flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_curved_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);


/**
@brief Kernel to run `extrap_stokes_curved_power_law_gpu` for all Stokes I components in `d_components`.

@details Fills the array `d_components.extrap_stokesV` with the extrapolated Stokes V flux densities.

@param[in] num_extrap_freqs The number of frequencies to extrapolate to.
@param[in] d_extrap_freqs Pointer to an array of frequencies to extrapolate to on the device.
@param[in] num_comps The number of components to extrapolate.
@param[in,out] d_components The components to extrapolate on the device.
*/
__global__ void kern_extrap_curved_power_laws_stokesV(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components);

/**
@brief Kernel to run `extrap_stokes_curved_power_law_gpu` for all linear polarisation components in `d_components`.

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

@param[in] list_stokes Array containing all the list-fluxes used for the extrapolation.
@param[in] list_freqs Array containing all the list-frequencies used for the extrapolation; must match order of `list_stokes`.
@param[in] arr_num_list_values Array containing the number of list values for each component.
@param[in] list_start_indexes Array containing the starting index of each component within `list_stokes` and `list_freqs`.
@param[in] d_extrap_freqs The frequencies to use for the extrapolation.
@param[in] iFluxComp The index of the flux component to store the result in.
@param[in] iFreq The index of the frequency component to extrapolate.
@param[in,out]  extrap_flux The extrapolated flux.
 */
__device__ void extrap_stokes_list_fluxes_gpu(user_precision_t *list_stokes,
     double *list_freqs, int *arr_num_list_values, int *list_start_indexes,
     double *d_extrap_freqs, int iFluxComp, int iFreq,
     user_precision_t * extrap_flux);

/**
@brief Extrapolates list-type spectral model fluxes to given frequencies. To be clear,
`list_stokes` and `list_freqs` are arrays of flux densities and frequencies, respectively,
for `num_comps` number of components. Each component can have any number of flux and freq
entries, given by `num_list_values`. 

@details Fills `extrap_stokes` with the extrapolated stokes flux densities. Runs the function `extrap_stokes_list_fluxes_gpu`. 

@param[in] list_stokes Array containing all the list-fluxes used for the extrapolation.
@param[in] list_freqs Array containing all the list-frequencies used for the extrapolation; must match order of `list_stokes`.
@param[in] num_list_values Array containing the number of list values for each component.
@param[in] list_start_indexes Array containing the starting index of each component within `list_stokes` and `list_freqs`.
@param[in] list_comp_inds Array containing the component index for list component; used to index the extrapolated fluxes to `extrap_stokes`.
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
* @brief Apply rotation measure to existing fluxes in `components.extrap_stokesQ`.
 * 
 * @details Fills `d_components.extrap_stokesQ` and `d_components.extrap_stokesU` with the rotated linear polarisation flux densities. Assumes the model
 * 
  \f[
  Q = \mathrm{P}(\lambda) \cos(2\chi_0 + 2\phi_{\textrm{RM}} \lambda^2), \\
  U = \mathrm{P}(\lambda) \sin(2\chi_0 + 2\phi_{\textrm{RM}} \lambda^2).
  \f]

  where \f$\mathrm{P}(\lambda)\f$ is the polarised flux, stored in `components.extrap_stokesQ`, \f$\chi_0\f$ is the intrinsic polarisation angle (`components.intr_pol_angle`), and \f$\phi_{\textrm{RM}}\f$ is the rotation measure (`components.rm_values`).

 * @param num_extrap_freqs The number of extrapolation frequencies.
 * @param d_extrap_freqs The frequencies to extrapolate to.
 * @param num_comps The number of components.
 * @param d_components The components to apply rotation measure to
 
 
 */
__global__ void kern_apply_rotation_measure(int num_extrap_freqs, double *d_extrap_freqs,
     int num_comps, components_t d_components);


/**
 * @brief Extrapolates power law SEDSs for Stokes I flux densities to 
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_power_laws_stokesI` to fill `d_components.extrap_stokesI`
 * with the extrapolated flux densities with inputs from `d_components`.
 *
 * @param d_components The device components structure containing necessary data.
 * @param n_powers The number of power laws to extrapolate.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of extrapolated frequencies.
 */
extern "C" void extrap_power_laws_stokesI_gpu(components_t d_components,
                                   int n_powers, double *d_extrap_freqs,
                                   int num_extrap_freqs);

/**
 * @brief Extrapolates curved power law SEds for Stokes I flux densities to 
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_curved_power_laws_stokesI` to fill `d_components.extrap_stokesI`
 * with the extrapolated flux densities with inputs from `d_components`.
 *
 * @param d_components The device components structure containing necessary data for extrapolation.
 * @param n_curves The number of curves to be extrapolated.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */
extern "C" void extrap_curved_power_laws_stokesI_gpu(components_t d_components,
                                   int n_curves, double *d_extrap_freqs,
                                   int num_extrap_freqs);


/**
 * @brief Extrapolate list-type SEDs for Stokes I flux densities to 
 * given frequencies on GPU.
 *
 * Calls `kern_extrap_list_fluxes` to fill `d_components.extrap_stokesI`
 * with the extrapolated flux densities with inputs from `d_components`.
 *
 * @param d_components The device pointer to the components data structure.
 * @param n_lists The number of lists to process.
 * @param d_extrap_freqs The device pointer to the array of frequencies at which to extrapolate.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */
extern "C" void extrap_list_fluxes_stokesI_gpu(components_t d_components,
                                   int n_lists, double *d_extrap_freqs,
                                   int num_extrap_freqs);

/**
 * @brief Extrapolates power law SEDSs for Stokes V flux densities to 
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_power_laws_stokesV` to fill `d_components.extrap_stokesV`
 * with the extrapolated flux densities with inputs from `d_components`.
 *
 * @param d_components The device components structure containing necessary data.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of extrapolated frequencies.
 */
extern "C" void extrap_power_laws_stokesV_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);

/**
 * @brief Extrapolates curved power law SEds for Stokes V flux densities to 
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_curved_power_laws_stokesV` to fill `d_components.extrap_stokesV`
 * with the extrapolated flux densities with inputs from `d_components`.
 *
 * @param d_components The device components structure containing necessary data for extrapolation.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */
extern "C" void extrap_curved_power_laws_stokesV_gpu(components_t d_components,
                                                       double *d_extrap_freqs,
                                                       int num_extrap_freqs);

/**
 * @brief Calculates polarisations fraction for Stokes V flux densities to 
 * given frequencies on the GPU.
 *
 * Calls `kern_polarisation_fraction_stokesV` to fill `d_components.extrap_stokesV`
 * with the extrapolated flux densities with inputs from `d_components`. Requires
 * `d_components.extrap_stokesI` to be filled first.
 *
 * @param d_components The device components structure containing necessary data for extrapolation.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */                                                    
extern "C" void polarisation_fraction_stokesV_gpu(components_t d_components,
                                                       double *d_extrap_freqs,
                                                       int num_extrap_freqs);

/**
 * @brief Extrapolate list-type SEDs for Stokes V flux densities to 
 * given frequencies on GPU.
 *
 * Calls `kern_extrap_list_fluxes` to fill `d_components.extrap_stokesV`
 * with the extrapolated flux densities with inputs from `d_components`.
 *
 * @param d_components The device pointer to the components data structure.
 * @param d_extrap_freqs The device pointer to the array of frequencies at which to extrapolate.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */
extern "C" void extrap_list_fluxes_stokesV_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);

/**
 * @brief Extrapolates power law SEDSs for linear polarisation flux 
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_power_laws_linpol` to fill `d_components.extrap_stokesQ`
 * with the extrapolated flux densities with inputs from `d_components`.
 * Once called, the function `kern_apply_rotation_measure` should be called to rotate
 * the linear polarisation fluxes into both `d_components.extrap_stokesQ` and
 * `d_components.extrap_stokesU`.
 *
 * @param d_components The device components structure containing necessary data.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of extrapolated frequencies.
 */
extern "C" void extrap_power_laws_linpol_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);

/**
 * @brief Extrapolates curved power law SEDs for linear polarisation flux for
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_curved_power_laws_linpol` to fill to fill `d_components.extrap_stokesQ`
 * with the extrapolated flux densities with inputs from `d_components`.
 * Once called, the function `kern_apply_rotation_measure` should be called to rotate
 * the linear polarisation fluxes into both `d_components.extrap_stokesQ` and
 * `d_components.extrap_stokesU`.
 *
 * @param d_components The device components structure containing necessary data for extrapolation.
 * @param n_curves The number of curves to be extrapolated.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */
extern "C" void extrap_curved_power_laws_linpol_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);

/**
 * @brief Calculates polarisations fraction for linear polarisation flux on
 * the GPU
 *
 * Calls `kern_polarisation_fraction_linpol` to fill `d_components.extrap_stokesQ`
 * with the extrapolated flux densities with inputs from `d_components`. Requires
 * `d_components.extrap_stokesI` to be filled first. Once called, the function
 * `kern_apply_rotation_measure` should be called to rotate the linear polarisation
 * fluxes into both `d_components.extrap_stokesQ` and `d_components.extrap_stokesU`.
 * 
 * @param d_components The device components structure containing necessary data for extrapolation.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */  
extern "C" void polarisation_fraction_linpol_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);

/**
 * @brief Extrapolate list-type SEDs for linear polarisation flux for
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_list_fluxes` twice to separately fill
 * `d_components.extrap_stokesQ` and `d_components.extrap_stokesU`. Assumes
 * information for boths Q and U list type fluxes are stored in `d_components`.
 * Does not required rotation measure to be applied as this model type
 * just requires lists of Q and U fluxes.
 *
 * @param d_components The device pointer to the components data structure.
 * @param d_extrap_freqs The device pointer to the array of frequencies at which to extrapolate.
 * @param num_extrap_freqs The number of frequencies in the d_extrap_freqs array.
 */
extern "C" void extrap_list_fluxes_linpol_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);


/**
 * @brief Extrapolate list-type SEDs for linear polarisation flux for
 * given frequencies on the GPU.
 *
 * Calls `kern_extrap_list_fluxes` to fill `d_components.extrap_stokesQ`
 * with the extrapolated flux densities with inputs from `d_components`.
 * Assumes `linpol_p` type list flux information exists in `d_components`.
 * Once called, the function `kern_apply_rotation_measure` should be called to rotate
 * the linear polarisation fluxes into both `d_components.extrap_stokesQ` and
 * `d_components.extrap_stokesU`.
 *
 * @param d_components The device components structure containing necessary data.
 * @param d_extrap_freqs Pointer to the device array of extrapolation frequencies.
 * @param num_extrap_freqs The number of extrapolated frequencies.
 */
extern "C" void extrap_p_list_fluxes_linpol_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);

/**
 * @brief Applies rotation measure to the given components on the GPU.
 *
 * Calls `kern_apply_rotation_measure`. Assumes `d_components.extrap_stokesQ` has
 * been filled with the total linear polarisation flux. The function
 * will fill `d_components.extrap_stokesQ` and `d_components.extrap_stokesU` with
 * the rotated linear polarisation fluxes.
 *
 * @param d_components The device pointer to the components structure.
 * @param d_extrap_freqs The device pointer to the array of extrapolated frequencies.
 * @param num_extrap_freqs The number of extrapolated frequencies.
 */
extern "C" void apply_rotation_measure_gpu(components_t d_components,
                                                  double *d_extrap_freqs,
                                                  int num_extrap_freqs);

/**
@brief Kernel to calculate the visibility response to a number `num_components`
of either POINT or GAUSSIAN COMPONENTs, and sum the outputs to `d_sum_visi_*_real`,
`d_sum_visi_*_imag`. Assumes all fluxes have been extrapolated to the requested
frequencies.

@details Uses the functions `calc_measurement_equation_gpu` and
`update_sum_visis_stokes*gpu` as detailed above to calculate the visibilities. Sets off
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
@param[in] d_component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
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
@param[in] off_cardinal_dipoles Boolean to indicate if the dipoles are off-cardinal i.e.
if off_cardinal == 1, then the dipoles are aligned at 45 and 135 degrees,
if off_cardinal == 0, then the dipoles are aligned at 0 and 90 degrees
*/
__global__ void kern_calc_visi_point_or_gauss(components_t d_components,
           beam_gains_t d_component_beam_gains,
           user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
           user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
           user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
           user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
           user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
           int num_components, int num_baselines, int num_freqs, int num_cross,
           int num_times, e_beamtype beamtype, e_component_type comptype,
           int off_cardinal_dipoles);

/**
@brief Kernel to calculate the visibility response to a number `num_shapes` of
SHAPELET COMPONENTs, and sum the outputs to `d_sum_visi_*_real`,
`d_sum_visi_*_imag`.

@details Uses the functions `extrap_stokes`, `calc_measurement_equation_gpu`,
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
`fundamental_coords_gpu::kern_calc_uvw_shapelet`.

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
@param[in] d_component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
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
@param[in] off_cardinal_dipoles Boolean to specify if the dipoles in the beam are off the cardinal axes
(45 and 135 degrees) or aligned with north-south and east-west (0 and 90 degrees). Effects what visibilities are calculated. 
*/
__global__ void kern_calc_visi_shapelets(components_t d_components,
      beam_gains_t d_component_beam_gains,
      user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
      user_precision_t *d_allsteps_wavelengths,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
      user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
      user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
      user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
      user_precision_t *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_cross,
      const int num_coeffs, int num_times, e_beamtype beamtype, int off_cardinal_dipoles);

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
extern "C" void free_components_gpu(source_t *d_chunked_source,
                                  e_component_type comptype);


/**
@brief Frees the device memory on `d_beam_gains`, given the type of
primary beam `beamtype`.

@details Only certain primary beam models have leakage terms, so need
to know the `beamtype` to free the correct locations.

@param[in,out] d_beam_gains A `beam_gains_t` with device memory to be
freed
@param[in] beamtype An `e_beamtype` describing the primary beam type

*/
extern "C" void free_beam_gains_gpu(beam_gains_t *d_beam_gains, e_beamtype beamtype);


/**
@brief Copies a populated `source_t` from host to device.

@details Depending on what is in chunked_source (e.g. POINT, GAUSSIAN, SHAPELET,
POWER_LAW, CURVED_POWER_LAW, LIST) type information, only what needs to be
copied across to the GPU is (no empty pointers arrays are copied across)

@param[in,out] *chunked_source A populated `source_t` struct
*/
extern "C" source_t * copy_chunked_source_to_GPU(source_t *chunked_source);



/**
 * @brief Copies component gains from CPU to GPU beam gains. This function
 * transfers the component gains stored in the CPU memory to the 
 * GPU memory, specifically to the beam gains structure. It is used to ensure 
 * that the GPU has the necessary data to perform computations involving beam gains.
 * 
 * @details Specifically, these member arrays `gxs, Dxs, Dys, gys` have memory allocated
 * on GPU in `d_beam_gains` and are copied from the `components` CPU memory.
 *
 * @param components Pointer to the components structure containing the gains on the CPU.
 * @param d_beam_gains Pointer to the beam gains structure on the GPU where the gains will be copied.
 * @param num_gains The number of gains to be copied from the CPU to the GPU.
 */
extern "C" void copy_CPU_component_gains_to_GPU_beam_gains(components_t *components,
                beam_gains_t *d_beam_gains, int num_gains);


/**
 * @brief Copies beam gains data from CPU memory to GPU memory.
 * 
 * @details Specifically, these member arrays `gxs, Dxs, Dys, gys` have memory allocated
 * on GPU in `d_beam_gains` and are copied from the `beam_gains` CPU memory.
 *
 * @param beam_gains Pointer to the beam gains data in CPU memory.
 * @param d_beam_gains Pointer to the beam gains data in GPU memory.
 * @param num_gains The number of beam gains to copy.
 */
extern "C" void copy_CPU_beam_gains_to_GPU_beam_gains(beam_gains_t *beam_gains,
  beam_gains_t *d_beam_gains, int num_gains);


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
extern "C" void free_extrapolated_flux_arrays_gpu(components_t *d_components);

/**
@brief Calculate the auto-correlations for all antennas given the fluxes
already calculated in `d_components` and beam gains in `d_component_beam_gains`.
Stores the outputs at the end of `d_sum_visi*`, after the cross-correlations.

@details `d_component_beam_gains` should have gains for
as many antennas (a.k.a stations/tiles) as detailed by `num_ants`.
This function then just does the dot product of the component fluxes and beam gains specific
to each baseline.

When called with `dim3 grid, threads`, kernel should be called with `grid.x`
and grid.y defined, where:
 - grid.x * threads.x >= `num_freqs*num_times`
 - grid.y * threads.y >= `num_ants`

@param[in] d_components d_components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] d_component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
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
Currently this is just an index of all antennas. Gets passed to `get_beam_gains_multibeams_gpu`
@param[in] *d_ant2_to_auto_map An index of all primary beams to auto-correlations
Currently this is just an index of all antennas. Gets passed to `get_beam_gains_multibeams_gpu`
@param[in] off_cardinal_dipoles Boolean to specify if the dipoles in the beam are off the cardinal axes
(45 and 135 degrees) or aligned with north-south and east-west (0 and 90 degrees). Effects what visibilities are calculated.

*/
__global__ void kern_calc_autos(components_t d_components,
                                beam_gains_t d_component_beam_gains,
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
                                int *d_ant2_to_auto_map,
                                int off_cardinal_dipoles);

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
extern "C" void fill_ant_to_baseline_mapping_gpu(int num_ants, int *d_ant1_to_baseline_map,
                                               int *d_ant2_to_baseline_map);


/**
 * @brief Calculate the visibility response to a number `num_components`
of either POINT or GAUSSIAN COMPONENTs, and sum the outputs to `sum_visi_*_real`,
`sum_visi_*_imag` in `d_visibility_set`. Assumes all fluxes have been
extrapolated to the requested frequencies.

@details Runs `kern_calc_visi_point_or_gauss`.

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
 *
@param[in] d_components components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] d_component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
@param d_calc_visi_inouts The input/output parameters for visibility calculation.
@param d_visibility_set The visibility set to update with the results.
@param num_components The number of components.
@param beamtype Beam type, see `woden_struct_defs.e_beamtype`
@param comptype Component type, either POINT or GAUSSIAN
@param woden_settings The filled WODEN settings struct
 */
extern "C" void calc_visi_point_or_gauss_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_components, e_beamtype beamtype,
                                        e_component_type comptype,
                                        woden_settings_t *woden_settings);

/**
@brief the visibility response to a number `num_shapes` of
SHAPELET COMPONENTs, and sum the outputs to `sum_visi_*_real`,
`sum_visi_*_imag` in `d_visibility_set`. Assumes all fluxes have been
extrapolated to the requested frequencies.

@details Uses the functions `calc_measurement_equation_cpu` and `update_sum_visis_stokesI*_cpu`
as detailed in `calc_visi_point_or_gauss_cpu` to calculate the visibilities in
serial on the CPU. Furthermore calculates the visibility envelope
\f$\mathrm{V}_{\mathrm{env}}\f$ to convert the basic visibility into a SHAPELET:

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
`fundamental_coords_gpu::kern_calc_uvw_shapelet`.

This differs from `calc_visi_point_or_gauss_gpu`
in that each SHAPELET component is built up from multiple shapelet basis
functions, and so the kernel (and sometimes the sky model) is split over the
shapelet basis functions, meaning multiple basis function calculations will
use the same COMPONENT \f$l,m,n\f$. The array `components.param_indexes` is
used to match the basis function information `components.n1s`, `components.n2s`,
`components.shape_coeffs` to the COPMONENT information (e.g. `components.ls`).

A further difference is that the shapelet \f$ u_{\mathrm{comp}},
v_{\mathrm{comp}} \f$ used in the visibility envelope equation should
be stored in metres only - for every baseline (fastest
changing), every time step, and every SHAPELET component (slowest changing).
They will be scaled by wavelength inside the kernel, as memory costs
become prohibitive to store for every wavelength as well. The visibility \f$u,v,w\f$
are stored for every baseline, wavelength, and time step.

`sbf` is the shapelet basis function lookup table, and should have been
generated with `create_sbf`
 *
@param[in] d_components components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] d_component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
@param d_calc_visi_inouts The input/output parameters for visibility calculation.
@param d_visibility_set The visibility set to update.
@param num_shapes The number of shapelets.
@param num_shape_coeffs The number of shapelet coefficients.
@param beamtype The type of beam.
@param woden_settings The WODEN settings.
 */
extern "C" void calc_visi_shapelets_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_shapes, int num_shape_coeffs,
                                        e_beamtype beamtype,
                                        woden_settings_t *woden_settings);


/**
 * @brief Allocate memory for beam gains on GPU.
 * 
 * @details Allocates the following memory:
  - `d_component_beam_gains->Dxs`
  - `d_component_beam_gains->Dys`
  - `d_component_beam_gains->gxs`
  - `d_component_beam_gains->gys`

  depending on what beamtype is given, the leakgages may or may not be allocated.
 *
 * @param d_component_beam_gains The beam gains struct to allocate memory in
 * @param beamtype The type of beam.
 * @param num_gains The number of gains.
 */
extern "C" void malloc_beam_gains_gpu(beam_gains_t *d_component_beam_gains,
                                     int beamtype, int num_gains);


/**
 * @brief Calculate the auto-correlations for all antennas given the fluxes
 * already calculated in `d_components` and beam gains in `d_component_beam_gains`.
 * Stores the outputs in `d_visibility_set`.
 * 
 * @details `d_component_beam_gains` should have gains for
 * as many antennas (a.k.a stations/tiles). This function
 * then just does the dot product of the component fluxes and beam gains specific
 * to each baseline.
 *
 * @param d_components Pointer to the device memory holding the components.
 * @param beam_settings Pointer to the beam settings structure.
 * @param d_component_beam_gains Pointer to the device memory holding the component beam gains.
 * @param d_visibility_set Pointer to the device memory holding the visibility set.
 * @param woden_settings Pointer to the WODEN settings structure.
 * @param num_components The number of components to process.
 * @param use_twobeams If True (1), use two primary beams per visibility.
 * Otherwise, assume all primary beams are identical
 */
extern "C" void calc_autos_gpu(components_t *d_components,
                               beam_settings_t *beam_settings,
                               beam_gains_t *d_component_beam_gains,
                               visibility_set_t *d_visibility_set,
                               woden_settings_t *woden_settings,
                               int num_components, int use_twobeams);

//------------------------------------------------------------------------------
//Wrapper kernels used in unit tests to test __device__ functions externally
//These are not used in the main code

/**
 @brief Wrapper kernel to test `calc_measurement_equation_gpu` in unit tests. 
*/
__global__ void kern_calc_measurement_equation(int num_components, int num_baselines,
          user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
          double *d_ls, double *d_ms, double *d_ns, gpuUserComplex *d_visis);

/**
 * @brief Wrapper kernel to test `apply_beam_gains_stokesI*_off_cardinal_gpu`
 * and `apply_beam_gains_stokesI*_on_cardinal_gpu` in unit tests.
 */
__global__ void kern_apply_beam_gains(int num_gains,
          gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
          gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
          gpuUserComplex *d_g2xs, gpuUserComplex *d_D2xs,
          gpuUserComplex *d_D2ys, gpuUserComplex *d_g2ys,
          user_precision_t *d_flux_Is, user_precision_t *d_flux_Qs,
          user_precision_t *d_flux_Us, user_precision_t *d_flux_Vs,
          gpuUserComplex *d_visi_components,
          gpuUserComplex *d_visi_XXs, gpuUserComplex *d_visi_XYs,
          gpuUserComplex *d_visi_YXs, gpuUserComplex *d_visi_YYs, 
          int off_cardinal_dipoles, int do_QUV);

/**
 * @brief Wrapper kernel to test `get_beam_gains_gpu` and
 * `get_beam_gains_multibeams_gpu` in unit tests.
 */
__global__ void kern_get_beam_gains(int num_components, int num_baselines,
           int num_freqs, int num_cross, int num_times, int beamtype,
           gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
           gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
           gpuUserComplex *d_recov_g1x, gpuUserComplex *d_recov_D1x,
           gpuUserComplex *d_recov_D1y, gpuUserComplex *d_recov_g1y,
           gpuUserComplex *d_recov_g2x, gpuUserComplex *d_recov_D2x,
           gpuUserComplex *d_recov_D2y, gpuUserComplex *d_recov_g2y,
           int use_twobeams, int num_ants,
           int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map);

/**
 * @brief Wrapper kernel to test `update_sum_visis_stokesIQUV_gpu` in unit tests.
 */
__global__ void kern_update_sum_visis_stokesIQUV(int num_freqs,
     int num_baselines, int num_components, int num_times,
     int beamtype, int off_cardinal_dipoles,
     gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
     gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
     int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
     gpuUserComplex *d_visi_components,
     user_precision_t *d_flux_I, user_precision_t *d_flux_Q,
     user_precision_t *d_flux_U, user_precision_t *d_flux_V,
     user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
     user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
     user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
     user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag);


