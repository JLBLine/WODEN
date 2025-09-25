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
           double l, double m, double n, double offset);


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
void apply_beam_gains_stokesI_on_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY);

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
void apply_beam_gains_stokesI_off_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
          user_precision_complex_t D1y, user_precision_complex_t g1y,
          user_precision_complex_t g2x, user_precision_complex_t D2x,
          user_precision_complex_t D2y, user_precision_complex_t g2y,
          user_precision_t flux_I, user_precision_complex_t visi_component,
          user_precision_complex_t * visi_XX, user_precision_complex_t * visi_XY,
          user_precision_complex_t * visi_YX, user_precision_complex_t * visi_YY);


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
void apply_beam_gains_stokesIQUV_off_cardinal_cpu(user_precision_complex_t g1x, user_precision_complex_t D1x,
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
calculate the observed XX,XY,YX,YY instrumental cross-correlated visibilities.

If `off_cardinal_dipoles` is set to 1, the off-cardinal terms are calculated, 
where we label north-east south-west (45 deg) aligned dipoles as 'x' and
south-east north-west (135 deg) as 'y'.

If `off_cardinal_dipoles` is set to 0, the on-cardinal terms are calculated, 
meaning 'x' is aligned north-south, 'y' aligned east-west.

If `do_QUV` is set to 1, the Q,U,V arrays are used to calculate XX,XY,YX,YY.
If `do_QUV` is set to 0, only I array is used to calculate XX,XY,YX,YY. Just
pass NULL pointers for Q,U,V arrays if `do_QUV` is set to 0.

@details Calls either `apply_beam_gains_stokesIQUV_off_cardinal_cpu`, 
`apply_beam_gains_stokesIQUV_on_cardinal_cpu`, `apply_beam_gains_stokesI_off_cardinal_cpu`,
or `apply_beam_gains_stokesI_on_cardinal_cpu`, depending on the boolean flags provided.

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
@param[in] off_cardinal_dipoles Flag to calculate off-cardinal dipoles
@param[in] do_QUV Flag to include Stokes QUV in visibility calculations
*/
void apply_beam_gains_arrays_cpu(int num_gains,
          user_precision_complex_t *g1xs, user_precision_complex_t *D1xs,
          user_precision_complex_t *D1ys, user_precision_complex_t *g1ys,
          user_precision_complex_t *g2xs, user_precision_complex_t *D2xs,
          user_precision_complex_t *D2ys, user_precision_complex_t *g2ys,
          user_precision_t *flux_Is, user_precision_t *flux_Qs,
          user_precision_t *flux_Us, user_precision_t *flux_Vs,
          user_precision_complex_t *visi_components,
          user_precision_complex_t *visi_XXs, user_precision_complex_t *visi_XYs,
          user_precision_complex_t *visi_YXs, user_precision_complex_t *visi_YYs,
          int off_cardinal_dipoles, int do_QUV);

/**
@brief Allocate host memory of extrapolated Stokes arrays int `components`

@details It does this if components.do_QUV == 1:

    components->extrap_stokesI = malloc(num_comps*num_freqs*sizeof(user_precision_t));
    components->extrap_stokesQ = malloc(num_comps*num_freqs*sizeof(user_precision_t));
    components->extrap_stokesU = malloc(num_comps*num_freqs*sizeof(user_precision_t));
    components->extrap_stokesV = malloc(num_comps*num_freqs*sizeof(user_precision_t));

If do_QUV == 0, only allocate the StokesI array.

@param[in,out] *components A populated `components_t` struct
@param[in] num_comps Number of components
@param[in] num_freqs Number of frequencies
*/
void malloc_extrapolated_flux_arrays_cpu(components_t *components, int num_comps,
                                     int num_freqs);

/**
@brief Free device memory of extrapolated Stokes arrays from `components`

@details It does this if components.do_QUV == 1:

   free(components->extrap_stokesI);
   free(components->extrap_stokesQ);
   free(components->extrap_stokesU);
   free(components->extrap_stokesV);

If do_QUV == 0, only free the StokesI array.

@param[in,out] *components A populated `components_t` struct
*/
void free_extrapolated_flux_arrays_cpu(components_t *components);

/**
@brief Assuming a simple power-law SED of \f$ S \propto \nu^\alpha \f$,
extrapolate flux density parameters to the requested frequencies

@details Assumes a referece frequency of 200MHz, and uses the
reference flux density `ref_flux` and spectral index `SI` to calculate the flux
density `extrap_flux` at the requested frequency `extrap_freq`.

@param[in] ref_flux Reference fluxes at 200MHz (Jy)
@param[in] SI Spectral index
@param[in] extrap_freq Frequency to extrapolate to
@param[in,out] *extrap_flux Pointer to extrapolated flux (Jy)
*/
void extrap_stokes_power_law_cpu(user_precision_t ref_flux, user_precision_t SI,
                             double extrap_freq,  user_precision_t * extrap_flux);

/**
@brief Assuming a curved power-law SED of \f$ S \propto \nu^\alpha \exp(q \ln (\nu)^2 )  \f$,
extrapolate flux density to the requested frequency

@details Assumes a referece frequency of 200MHz, and uses the
reference flux density `ref_flux`, spectral index `SI`, and curvature `q` param
to calculate the flux density `extrap_flux` at the requested frequency `extrap_freq`.

@param[in] ref_flux Reference fluxes at 200MHz (Jy)
@param[in] SI Spectral index
@param[in] q Curvature parameter
@param[in] extrap_freq Frequency to extrapolate to
@param[in,out] *extrap_flux Pointer to extrapolated flux (Jy)
*/
void extrap_stokes_curved_power_law_cpu(user_precision_t ref_flux,
           user_precision_t SI, user_precision_t q,
           double extrap_freq, user_precision_t * extrap_flux);


/**
@brief Extrapolates list-type spectral model fluxes to given frequencies. To be clear,
`list_stokes` and `list_freqs` are arrays of flux densities and frequencies, respectively,
for `num_comps` number of components. Each component can have any number of flux and freq
entries, given by `num_list_values`. 

@details Fills `extrap_stokes` with `num_comps` extrapolated stokes flux densities. Runs the function `extrap_stokes_list_flux_cpu`. 

@param[in] list_stokes Array containing all the list-fluxes used for the extrapolation.
@param[in] list_freqs Array containing all the list-frequencies used for the extrapolation; must match order of `list_stokes`.
@param[in] num_list_values Array containing the number of list values for each component.
@param[in] list_start_indexes Array containing the starting index of each component within `list_stokes` and `list_freqs`.
@param[in] list_comp_inds Array containing the component index for list component; used to index the extrapolated fluxes to `extrap_stokes`.
@param[in] num_extrap_freqs The number of extrapolation frequencies.
@param[in] extrap_freqs Pointer to the array of extrapolation frequencies.
@param[in] num_comps The number of components.
@param[in,out] extrap_stokes Output extrapolated fluxes
 */
void extrap_stokes_list_flux_arrays_cpu(user_precision_t *list_stokes, double *list_freqs,
                                int *num_list_values, int *list_start_indexes,
                                int *list_comp_inds,
                                int num_extrap_freqs, double *extrap_freqs,
                                int num_comps, user_precision_t *extrap_stokes);


/**
 * @brief Extrapolate power laws for Stokes I in `components`.
 * 
 * @details Fills `components.extrap_stokesI` with `num_comps` extrapolated stokes I flux densities. Runs the function `extrap_stokes_power_law_cpu`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param num_components The number of components.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_power_laws_stokesI_cpu(components_t components, int num_components,
                   double *extrap_freqs, int num_extrap_freqs);

/**
 * @brief Extrapolate power laws for Stokes V component in `components`.
 * 
 * @details Fills `components.extrap_stokesV` with `num_comps` extrapolated stokes V flux densities. Runs the function `extrap_stokes_power_law_cpu`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_power_laws_stokesV_cpu(components_t components, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Extrapolate power laws for linear polarisation in `components`.
 * 
 * @details Temporarily fills the array `components.extrap_stokesQ` with the extrapolated linear polarisation flux densities. The intention is to the use `apply_rotation_measure_cpu` after the fact, which will use `components.extrap_stokesQ` as input flux, do rotations, and then fill `components.extrap_stokesQ` and `components.extrap_stokesU`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_power_laws_linpol_cpu(components_t components, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Extrapolate curved power laws for Stokes I component in `components`.
 * 
 * @details Fills `components.extrap_stokesI` with `num_comps` extrapolated stokes I flux densities. Runs the function `extrap_stokes_curved_power_law_cpu`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param num_components The number of components.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_curved_power_laws_stokesI_cpu(components_t components, int num_components,
                   double *extrap_freqs, int num_extrap_freqs);

/**
 * @brief Extrapolate curved power laws for Stokes V component in `components`.
 * 
 * @details Fills `components.extrap_stokesV` with `num_comps` extrapolated stokes V flux densities. Runs the function `extrap_stokes_curved_power_law_cpu`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_curved_power_laws_stokesV_cpu(components_t components,
                  double *extrap_freqs, int num_extrap_freqs);

/**
 * @brief Extrapolate curved power laws for linear polarisation in `components`.
 * 
 * @details Temporarily fills the array `components.extrap_stokesQ` with the extrapolated linear polarisation flux densities. The intention is to the use `apply_rotation_measure_cpu` after the fact, which will use `components.extrap_stokesQ` as input flux, do rotations, and then fill `components.extrap_stokesQ` and `components.extrap_stokesU`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_curved_power_laws_linpol_cpu(components_t components,
                  double *extrap_freqs, int num_extrap_freqs);

/**
 * @brief Extrapolate list fluxes for Stokes I in `components`.
 * 
 * @details Fills `components.extrap_stokesI` with `num_comps` extrapolated stokes I flux densities. Runs the function `extrap_stokes_list_flux_arrays_cpu`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param n_lists The number of lists.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_list_fluxes_stokesI_cpu(components_t components,
                   int n_lists, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Extrapolate list fluxes for Stokes V component in `components`.
 * 
 * @details Fills `components.extrap_stokesV` with `num_comps` extrapolated stokes V flux densities. Runs the function `extrap_stokes_list_flux_arrays_cpu`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_list_fluxes_stokesV_cpu(components_t components, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Extrapolate list fluxes for linear polarisation in `components`.
 * 
 * @details Assumes list entries for both Stokes Q and U are present in `components`
 * Runs the function `extrap_stokes_list_flux_arrays_cpu` for both, filling `components.extrap_stokesQ` and `components.extrap_stokesU`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_list_fluxes_linpol_cpu(components_t components, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Extrapolate p-list fluxes for linear polarisation on CPU. p-lists are lists of fluxes for linear polarisation, which is then used in combination with the rotation measure to calculate the final linear polarisation fluxes.
 * 
 * @details Temporarily fills the array `components.extrap_stokesQ` with the extrapolated linear polarisation flux densities. The intention is to the use `apply_rotation_measure_cpu` after the fact, which will use `components.extrap_stokesQ` as input flux, do rotations, and then fill `components.extrap_stokesQ` and `components.extrap_stokesU`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void extrap_p_list_fluxes_linpol_cpu(components_t components, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Calculate polarisation fraction for Stokes V in `components`.
 * 
 * @details Fills `components.extrap_stokesV` with `num_comps` polarisation fractions for stokes V. Requires `components.extrap_stokesI` to be filled.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void polarisation_fraction_stokesV_cpu(components_t components, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Calculate polarisation fraction for linear polarisation in `components`.
 * 
 * @details Fills `components.extrap_stokesQ` with `num_comps` polarisation fractions for linear polarisation. Requires `components.extrap_stokesI` to be filled. Intention is to use `apply_rotation_measure_cpu` after the fact, which will use `components.extrap_stokesQ` as input flux, do rotations, and then fill `components.extrap_stokesQ` and `components.extrap_stokesU`.
 *
 * @param components The filled `components` object to run extrapolation on.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void polarisation_fraction_linpol_cpu(components_t components, double *extrap_freqs,
                   int num_extrap_freqs);

/**
 * @brief Apply rotation measure to existing fluxes in `components.extrap_stokesQ`.
 * 
 * @details Fills `components.extrap_stokesQ` and `components.extrap_stokesU` with the rotated linear polarisation flux densities. Assumes the model
 * 
  \f[
  Q = \mathrm{P}(\lambda) \cos(2\chi_0 + 2\phi_{\textrm{RM}} \lambda^2), \\
  U = \mathrm{P}(\lambda) \sin(2\chi_0 + 2\phi_{\textrm{RM}} \lambda^2).
  \f]

  where \f$\mathrm{P}(\lambda)\f$ is the polarised flux, stored in `components.extrap_stokesQ`, \f$\chi_0\f$ is the intrinsic polarisation angle (`components.intr_pol_angle`), and \f$\phi_{\textrm{RM}}\f$ is the rotation measure (`components.rm_values`).

 *
 * @param components The components to apply rotation measure.
 * @param extrap_freqs The frequencies to extrapolate to.
 * @param num_extrap_freqs The number of extrapolation frequencies.
 */
void apply_rotation_measure_cpu(components_t components,
                double *extrap_freqs, int num_extrap_freqs);

/**
@brief Given the type of primary beam simulated `beamtype`, select the beam
gains and leakages from the arrays `g*s, D*s` that match the indexes
`iBaseline`, `iComponent`. NOTE that currently, primary beams are assumed
to be the same for all recieving elements, and so this function only takes in
one set of primary beam arrays.

@details This function is built to return the correct beam gain for a given
component on the sky, at a given time, at a given frequency, for a given baseline.
The 4 arrays `gxs`, `Dxs`, `Dys`,`gys` should contain the primary beam
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
@param[in] *gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *gys Pointer towards array of primary beam J[1,1]
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
void get_beam_gains_cpu(int iBaseline, int iComponent, int num_freqs,
       int num_baselines, int num_components, int num_times, int beamtype,
       user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
       user_precision_complex_t *Dys, user_precision_complex_t *gys,
       user_precision_complex_t * g1x, user_precision_complex_t * D1x,
       user_precision_complex_t * D1y, user_precision_complex_t * g1y,
       user_precision_complex_t * g2x, user_precision_complex_t * D2x,
       user_precision_complex_t * D2y, user_precision_complex_t * g2y);

/**
@brief Given the type of primary beam simulated `beamtype`, select the beam
gains and leakages from the arrays `g*s, D*s` that match the indexes
`iBaseline`, `iComponent`. This function assumes the primary beams are different
for every antenna, so returns different values for each antenna. 

@details This function is built to return the correct beam gain for a given
component on the sky, at a given time, at a given frequency, for a given baseline.
The 4 arrays  `gxs`, `Dxs`, `Dys`,`gys`
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
@param[in] *gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *gys Pointer towards array of primary beam J[1,1]
(east-west gain)
@param[in] *ant1_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1
@param[in] *ant2_to_baseline_map The index of antenna 1 in all unique pairs of
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
void get_beam_gains_multibeams_cpu(int iBaseline, int iComponent, int num_freqs,
       int num_baselines, int num_components, int num_times, int beamtype,
       user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
       user_precision_complex_t *Dys, user_precision_complex_t *gys,
       int *ant1_to_baseline_map, int *ant2_to_baseline_map,
       user_precision_complex_t * g1x, user_precision_complex_t * D1x,
       user_precision_complex_t * D1y, user_precision_complex_t * g1y,
       user_precision_complex_t * g2x, user_precision_complex_t * D2x,
       user_precision_complex_t * D2y, user_precision_complex_t * g2y);

/**
@brief Given the visibility between two recieving elements for a COMPONENT
`visi_component`, at the baseline/freq/time index given by `iBaseline` and
COMPONENT index `iComponent`, apply the instrumental primary beam and
COMPONENT Stokes I parameter to create instrumental XX,XY,YX,YY visibilities,
and sum them into real and imaginary XX,XY,YX,YY visibilities arrays
`sim_visi_*_real` and `sim_visi_*_imag`.

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
@param[in] off_cardinal_dipoles Boolean to indicate if the dipoles are off-cardinal i.e.
if off_cardinal == 1, then the dipoles are aligned at 45 and 135 degrees,
if off_cardinal == 0, then the dipoles are aligned at 0 and 90 degrees
@param[in] *gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *gys Pointer towards array of primary beam J[1,1]
(east-west gain)
@param[in] *ant1_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1 (only used when use_twobeams == 1)
@param[in] *ant2_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 2 (only used when use_twobeams == 1)
@param[in] use_twobeams If 1 (True), assume all primary beams are different for each
antenna. If 0 (False), assume all primary beams are the same for all antennas.
@param[in] visi_component Complex visibility across antennas 1 and 2
@param[in] flux_I Stokes I flux density (Jy)
@param[in,out] *sum_visi_XX_real Pointer to array to sum real XX visibility
into
@param[in,out] *sum_visi_XX_imag Pointer to array to sum imaginary XX
visibility into
@param[in,out] *sum_visi_XY_real Pointer to array to sum real XY visibility
into
@param[in,out] *sum_visi_XY_imag Pointer to array to sum imaginary XY
visibility into
@param[in,out] *sum_visi_YX_real Pointer to array to sum real YX visibility
into
@param[in,out] *sum_visi_YX_imag Pointer to array to sum imaginary YX
visibility into
@param[in,out] *sum_visi_YY_real Pointer to array to sum real YY visibility
into
@param[in,out] *sum_visi_YY_imag Pointer to array to sum imaginary YY
visibility into
*/
void update_sum_visis_stokesI_cpu(int iBaseline, int iComponent,
  int num_freqs, int num_baselines, int num_components, int num_times,
  int beamtype, int off_cardinal_dipoles,
  user_precision_complex_t *gxs, user_precision_complex_t *Dxs,
  user_precision_complex_t *Dys, user_precision_complex_t *gys,
  int *ant1_to_baseline_map, int *ant2_to_baseline_map, int use_twobeams,
  user_precision_complex_t visi_component, user_precision_t flux_I,
  user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
  user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
  user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
  user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag);

/**
@brief Given the visibility between two recieving elements for a COMPONENT
`visi_component`, at the baseline/freq/time index given by `iBaseline` and
COMPONENT index `iComponent`, apply the instrumental primary beam and
COMPONENT Stokes parameters to create instrumental XX,XY,YX,YY visibilities,
and sum them into real and imaginary XX,XY,YX,YY visibilities arrays
`sim_visi_*_real` and `sim_visi_*_imag`.

@details Uses `get_beam_gains_gpu` and `apply_beam_gains` as described above to
apply the gains - see descriptions for what should be the arguments to them.

@param[in] iBaseline Index of which baseline, freq, and time we are on
@param[in] iComponent COMPONENT index
@param[in] num_freqs Number of frequencies in simulation
@param[in] num_baselines Number of baselines for one time and one frequency step
@param[in] num_components Number of COMPONENTs
@param[in] num_times Number of times in simulation
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] off_cardinal_dipoles Boolean to indicate if the dipoles are off-cardinal i.e.
if off_cardinal == 1, then the dipoles are aligned at 45 and 135 degrees,
if off_cardinal == 0, then the dipoles are aligned at 0 and 90 degrees
@param[in] *gxs Pointer towards array of primary beam J[0,0]
(north-south gain)
@param[in] *Dxs Pointer towards array of primary beam J[0,1]
(north-south leakage)
@param[in] *Dys Pointer towards array of primary beam J[1,0]
(east-west leakage)
@param[in] *gys Pointer towards array of primary beam J[1,1]
(east-west gain)
@param[in] *ant1_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 1 (only used when use_twobeams == 1)
@param[in] *ant2_to_baseline_map The index of antenna 1 in all unique pairs of
antennas. Used to map iBaseline to the correct antenna 2 (only used when use_twobeams == 1)
@param[in] use_twobeams If 1 (True), assume all primary beams are different for each
antenna. If 0 (False), assume all primary beams are the same for all antennas.
@param[in] visi_component Complex visibility across antennas 1 and 2
@param[in] flux_I Stokes I flux density (Jy)
@param[in] flux_Q Stokes Q flux density (Jy)
@param[in] flux_U Stokes U flux density (Jy)
@param[in] flux_V Stokes V flux density (Jy)
@param[in,out] *sum_visi_XX_real Pointer to array to sum real XX visibility
into
@param[in,out] *sum_visi_XX_imag Pointer to array to sum imaginary XX
visibility into
@param[in,out] *sum_visi_XY_real Pointer to array to sum real XY visibility
into
@param[in,out] *sum_visi_XY_imag Pointer to array to sum imaginary XY
visibility into
@param[in,out] *sum_visi_YX_real Pointer to array to sum real YX visibility
into
@param[in,out] *sum_visi_YX_imag Pointer to array to sum imaginary YX
visibility into
@param[in,out] *sum_visi_YY_real Pointer to array to sum real YY visibility
into
@param[in,out] *sum_visi_YY_imag Pointer to array to sum imaginary YY
visibility into
*/
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

/**
 * @brief Allocate memory for beam gains on CPU.
 * 
 * @details Allocates the following memory:
  - `component_beam_gains->Dxs`
  - `component_beam_gains->Dys`
  - `component_beam_gains->gxs`
  - `component_beam_gains->gys`

  depending on what beamtype is given, the leakgages may or may not be allocated.
 *
 * @param component_beam_gains The beam gains struct to allocate memory in
 * @param beamtype The type of beam.
 * @param num_gains The number of gains.
 */
void malloc_beam_gains_cpu(beam_gains_t *component_beam_gains,
               int beamtype, int num_gains);

/**
 * @brief Calculate the visibility response to a number `num_components`
of either POINT or GAUSSIAN COMPONENTs, and sum the outputs to `sum_visi_*_real`,
`sum_visi_*_imag` in `visibility_set`. Assumes all fluxes have been
extrapolated to the requested frequencies.

@details Uses the functions `calc_measurement_equation_cpu` and `update_sum_visis_stokesI*_cpu`
as detailed above to calculate the visibilities in serial on the CPU.

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
@param[in] components components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
@param calc_visi_inouts The input/output parameters for visibility calculation.
@param visibility_set The visibility set to update with the results.
@param num_components The number of components.
@param beamtype Beam type, see `woden_struct_defs.e_beamtype`
@param comptype Component type, either POINT or GAUSSIAN
@param woden_settings The filled WODEN settings struct
 */
void calc_visi_point_or_gauss_cpu(components_t components,
                  beam_gains_t component_beam_gains,
                  calc_visi_inouts_t *calc_visi_inouts,
                  visibility_set_t *visibility_set, 
                  int num_components, e_beamtype beamtype,
                  e_component_type comptype,
                  woden_settings_t *woden_settings);

/**
@brief the visibility response to a number `num_shapes` of
SHAPELET COMPONENTs, and sum the outputs to `sum_visi_*_real`,
`sum_visi_*_imag` in `visibility_set`. Assumes all fluxes have been
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

This differs from `calc_visi_point_or_gauss_cpu`
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
@param[in] components components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
@param calc_visi_inouts The input/output parameters for visibility calculation.
@param visibility_set The visibility set to update.
@param num_shapes The number of shapelets.
@param num_shape_coeffs The number of shapelet coefficients.
@param beamtype The type of beam.
@param woden_settings The WODEN settings.
 */
void calc_visi_shapelets_cpu(components_t components,
               beam_gains_t component_beam_gains,
               calc_visi_inouts_t *calc_visi_inouts,
               visibility_set_t *visibility_set,
               int num_shapes, int num_shape_coeffs,
               e_beamtype beamtype,
               woden_settings_t *woden_settings);

/**
@brief Calculate the auto-correlations for all antennas given the fluxes
already calculated in `components` and beam gains in `component_beam_gains`.
Stores the outputs at the end of `sum_visi*`, after the cross-correlations.

@details `component_beam_gains` should have gains for
as many antennas (a.k.a stations/tiles) as detailed by `num_ants`. This function
then just does the dot product of the component fluxes and beam gains specific
to each baseline.

@param[in] components components Pointer to a populated `components_t` struct as filled by `source_components::source_component_common`
@param[in] component_beam_gains Pointer to a populated `beam_gains_t` struct as filled by `source_components::source_component_common`
@param[in] beamtype Beam type see `woden_struct_defs.e_beamtype`
@param[in] num_components Number of components in `components`
@param[in] num_baselines Number of baselines for one time, one frequency step
@param[in] num_freqs Number of frequncies in the simulation
@param[in] num_times Number of times in the simulation
@param[in] num_ants Number of antennas in the array
@param[in,out] *sum_visi_XX_real Pointer to array to sum real XX visibility
into
@param[in,out] *sum_visi_XX_imag Pointer to array to sum imaginary XX
visibility into
@param[in,out] *sum_visi_XY_real Pointer to array to sum real XY visibility
into
@param[in,out] *sum_visi_XY_imag Pointer to array to sum imaginary XY
visibility into
@param[in,out] *sum_visi_YX_real Pointer to array to sum real YX visibility
into
@param[in,out] *sum_visi_YX_imag Pointer to array to sum imaginary YX
visibility into
@param[in,out] *sum_visi_YY_real Pointer to array to sum real YY visibility
into
@param[in,out] *sum_visi_YY_imag Pointer to array to sum imaginary YY
visibility into
time step in the simulation
@param[in] use_twobeams If True, use a two primary beams per visibility.
Otherwise, assume all primary beams are identical
@param[in] *ant1_to_auto_map An index of all primary beams to auto-correlations
Currently this is just an index of all antennas. Gets passed to `get_beam_gains_multibeams_gpu`
@param[in] *ant2_to_auto_map An index of all primary beams to auto-correlations
Currently this is just an index of all antennas. Gets passed to `get_beam_gains_multibeams_gpu`
@param[in] off_cardinal_dipoles Boolean to specify if the dipoles in the beam are off the cardinal axes
(45 and 135 degrees) or aligned with north-south and east-west (0 and 90 degrees). Effects what visibilities are calculated.

*/
void calc_autos_cpu(components_t components, beam_gains_t component_beam_gains,
          int beamtype, int num_components, int num_baselines,
          int num_freqs, int num_times, int num_ants,
          user_precision_t *sum_visi_XX_real, user_precision_t *sum_visi_XX_imag,
          user_precision_t *sum_visi_XY_real, user_precision_t *sum_visi_XY_imag,
          user_precision_t *sum_visi_YX_real, user_precision_t *sum_visi_YX_imag,
          user_precision_t *sum_visi_YY_real, user_precision_t *sum_visi_YY_imag,
          int use_twobeams, int *ant1_to_auto_map, int *ant2_to_auto_map,
          int off_cardinal_dipoles);