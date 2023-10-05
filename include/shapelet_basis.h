/*! \file shapelet_basis.h
  A function to create a look-up shapelet basis function array, with relevant
  constants definitions.
*/
#include "woden_precision_defs.h"
#include "constants.h"

/**
@brief Create the 1D shapelet basis function look-up array `sbf`

@details `sbf` should have `sbf_N * sbf_L * sizeof(user_precision_t)` of memory allocated.

Each 1D basis function is \f$B_n(x, \beta)\f$, with
\f$n\f$ the basis function order, and \f$\beta \f$ the scaling factor. For all
stored arrays, I have set \f$\beta=1\f$. Pre-written arrays of
basis function orders \f$n=0\f$ to \f$n=100\f$ are then copied into `sbf`,
each of length `sbf_L`. Each samples over the range \f$-500\leq x \leq 500 \f$,
with a resolution \f$ \Delta x = 0.01 \f$

This single lookup array can then be used to generate 2D shapelet basis
functions later on varying \f$\beta\f$ values via a scaling factor.

@param[in] sbf An array with `sbf_N * sbf_L * sizeof(user_precision_t)` of memory allocated
@returns The populated `sbf` array
*/
user_precision_t * create_sbf(user_precision_t *sbf);
