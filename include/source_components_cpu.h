/*! \file
  GPU methods to extrapolate flux densities, assign and apply primary beam
  gains, and visibility responses for different sky model COMPONENT types.
*/
#pragma once
// #include "calculate_visibilities_cpu.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

#include "fundamental_coords_cpu.h"
#include "primary_beam_cpu.h"
#include "source_components_common.h"

/**
@brief Calculate the complex exponential phase part of the visibility function
from the given u,v,w and l,m,n coords

@details Calculates the following complex exponential:

\f[
\exp \left( 2\pi i\left( ul + vm + w(n-1) \right) \right)
\f]

returning a complex number.

Note there is no negative sign here - testing with the current implementation
through imaging software, and from comparisons to the OSKAR and RTS simulators
this yields the correct output visibilities.

@param[in] u \f$u\f$ coordinate (wavelengths)
@param[in] v \f$v\f$ coordinate (wavelengths)
@param[in] w \f$w\f$ coordinate (wavelengths)
@param[in] l \f$l\f$ coordinates
@param[in] m \f$m\f$ coordinates
@param[in] n \f$n\f$ coordinates
@return `visi`, a `user_precision_complex_t` of the visibility
*/
user_precision_complex_t calc_measurement_equation_cpu(user_precision_t u,
           user_precision_t v, user_precision_t w,
           double l, double m, double n);


/**
@brief Calculate the complex exponential phase part of the visibility function
from the given u,v,w and l,m,n coord arrays. Puts the results in the `visis`
array

@details Calls `calc_measurement_equation_cpu` for each baseline and component
pair, storing the results in the `visis` array. Index of output visibilities
is `num_components*iBaseline + iComponent`.

@param[in] num_cross Total number of `u,v,w` coords to be calculated (number of cross correlations)
@param[in] num_components Number of components to calculate visibilities for (number of `l,m,n` coords)
@param[in] *us Pointer to \f$u\f$ coodinates (wavelengths)
@param[in] *vs Pointer to \f$v\f$ coodinates (wavelengths)
@param[in] *ws Pointer to \f$w\f$ coodinates (wavelengths)
@param[in] *ls Pointer to \f$l\f$ coodinates
@param[in] *ms Pointer to \f$m\f$ coodinates
@param[in] *ns Pointer to \f$n\f$ coodinates
@param[in, out] *visis Pointer to the visibility array of size `num_cross*num_components`

*/
void calc_measurement_equation_arrays_cpu(int num_cross, int num_components,
           user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
           double *ls, double *ms, double *ns, user_precision_complex_t *visis);

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
void apply_beam_gains_stokesIQUV_on_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY);


/**
@brief Given primary beam gains and leakage terms for antenna 1
`g1xs, D1xs, D1ys, g1ys` and antenna 2 `g2x, D2x, D2y, g2ys`, the complex visibility
across those two antennas `visi`, and the Stokes parameters `flux*`,
simulate the observed XX,XY,YX,YY instrumental cross-correlated visibilities,
where 'x' means aligned north-south, 'y' means east-west.

@details Calls `apply_beam_gains_stokesIQUV_on_cardinal_cpu` for `num_gains`
visibilities, storing the results in the `visi_XXs`, `visi_XYs`, `visi_YXs`,
and `visi_YYs` arrays.

@param[in] num_gains Number of gains
@param[in] *g1xs Pointer to beam gain antenna 1 in north-south
@param[in] *D1xs Pointer to beam leakage antenna 1 from north-south
@param[in] *D1ys Pointer to beam gain antenna 1 in east-west
@param[in] *g1ys Pointer to beam leakage antenna 1 from east-west
@param[in] *g2xs Pointer to beam gain antenna 2 in north-south
@param[in] *D2xs Pointer to beam leakage antenna 2 from north-south
@param[in] *D2ys Pointer to beam gain antenna 2 in east-west
@param[in] *g2ys Pointer to beam leakage antenna 2 from east-west
@param[in] *flux_Is Pointer to Stokes I flux density (Jy)
@param[in] *flux_Qs Pointer to Stokes Q flux density (Jy)
@param[in] *flux_Us Pointer to Stokes U flux density (Jy)
@param[in] *flux_Vs Pointer to Stokes V flux density (Jy)
@param[in] *visi_components Pointer to complex visibility across antennas 1 and 2
@param[in,out] *visi_XXs Pointer to output XX instrumental visibility
@param[in,out] *visi_XYs Pointer to output XY instrumental visibility
@param[in,out] *visi_YXs Pointer to output YX instrumental visibility
@param[in,out] *visi_YYs Pointer to output YY instrumental visibility
*/
void apply_beam_gains_stokesIQUV_on_cardinal_arrays_cpu(int num_gains,
          user_precision_complex_t *g1xs, user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *g2xs, user_precision_complex_t *D2xs,
          user_precision_complex_t *D2ys, user_precision_complex_t *g2ys,
          user_precision_t *flux_Is, user_precision_t *flux_Qs,
          user_precision_t *flux_Us, user_precision_t *flux_Vs,
          user_precision_complex_t *visi_components,
          user_precision_complex_t *visi_XXs, user_precision_complex_t *visi_XYs,
          user_precision_complex_t *visi_YXs, user_precision_complex_t *visi_YYs);


void apply_beam_gains_stokesIQUV_off_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY);

void apply_beam_gains_stokesIQUV_off_cardinal_arrays_cpu(int num_gains,
          user_precision_complex_t *g1xs, user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *g2xs, user_precision_complex_t *D2xs,
          user_precision_complex_t *D2ys, user_precision_complex_t *g2ys,
          user_precision_t *flux_Is, user_precision_t *flux_Qs,
          user_precision_t *flux_Us, user_precision_t *flux_Vs,
          user_precision_complex_t *visi_components,
          user_precision_complex_t *visi_XXs, user_precision_complex_t *visi_XYs,
          user_precision_complex_t *visi_YXs, user_precision_complex_t *visi_YYs);

void malloc_extrapolated_flux_arrays_cpu(components_t *components, int num_comps,
                                     int num_freqs);

void free_extrapolated_flux_arrays_cpu(components_t *components);

void extrap_stokes_power_law_cpu(user_precision_t ref_flux, user_precision_t SI,
                             double extrap_freq,  user_precision_t * extrap_flux);

void extrap_stokes_curved_power_law_cpu(user_precision_t ref_flux,
           user_precision_t SI, user_precision_t q,
           double extrap_freq, user_precision_t * extrap_flux);

void extrap_stokes_list_flux_arrays_cpu(user_precision_t *list_stokes, double *list_freqs,
                                int *num_list_values, int *list_start_indexes,
                                int *list_comp_inds,
                                int num_extrap_freqs, double *extrap_freqs,
                                int num_comps, user_precision_t *extrap_stokes);

void extrap_power_laws_stokesI_cpu(components_t components, int num_components,
                                   double *extrap_freqs, int num_extrap_freqs);

void extrap_power_laws_stokesV_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs);

void extrap_power_laws_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs);

void extrap_curved_power_laws_stokesI_cpu(components_t components, int num_components,
                                   double *extrap_freqs, int num_extrap_freqs);

void extrap_curved_power_laws_stokesV_cpu(components_t components,
                                    double *extrap_freqs, int num_extrap_freqs);
void extrap_curved_power_laws_linpol_cpu(components_t components,

                                    double *extrap_freqs, int num_extrap_freqs);

void extrap_list_fluxes_stokesI_cpu(components_t components,
                                   int n_lists, double *extrap_freqs,
                                   int num_extrap_freqs);

void extrap_list_fluxes_stokesV_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs);

void extrap_list_fluxes_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs);

void extrap_p_list_fluxes_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs);

void polarisation_fraction_stokesV_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs);

void polarisation_fraction_linpol_cpu(components_t components, double *extrap_freqs,
                                   int num_extrap_freqs);

void apply_rotation_measure_cpu(components_t components,
                                double *extrap_freqs, int num_extrap_freqs);

void get_beam_gains_cpu(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           user_precision_complex_t * g1x, user_precision_complex_t * D1x,
           user_precision_complex_t * D1y, user_precision_complex_t * g1y,
           user_precision_complex_t * g2x, user_precision_complex_t * D2x,
           user_precision_complex_t * D2y, user_precision_complex_t * g2y);

void get_beam_gains_multibeams_cpu(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
           user_precision_complex_t *Dys, user_precision_complex_t *gys,
           int *ant1_to_baseline_map, int *ant2_to_baseline_map,
           user_precision_complex_t * g1x, user_precision_complex_t * D1x,
           user_precision_complex_t * D1y, user_precision_complex_t * g1y,
           user_precision_complex_t * g2x, user_precision_complex_t * D2x,
           user_precision_complex_t * D2y, user_precision_complex_t * g2y);

void update_sum_visis_stokesIQUV_cpu(int iBaseline, int iComponent,
    int num_freqs, int num_baselines, int num_components, int num_times,
    int beamtype, int off_cardinal_dipoles,
    user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
    user_precision_complex_t *Dys, user_precision_complex_t *gys,
    int *ant1_to_baseline_map, int *ant2_to_baseline_map, int use_twobeams,
    user_precision_complex_t visi_component,
    user_precision_t flux_I, user_precision_t flux_Q,
    user_precision_t flux_U, user_precision_t flux_V,
    user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
    user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
    user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
    user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag);