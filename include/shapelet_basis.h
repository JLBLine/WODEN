/*! \file shapelet_basis.h
  A function to create a look-up shapelet basis function array, with relevant
  constants definitions.
*/

//Somthing WEIRD is up with the way I'm using the documentation package
//'breathe', so having to include the definition value in the documentation
//string to get values to appear (sigh)

/*! 101 \n
Number of orders of basis functions (from 0 to 100 inclusive) */
#define sbf_N 101
/*! 10001 \n
Number of samples of each order of the basis functions */
#define sbf_L 10001
/*! 5000 \n
If shapelet basis function is B(x), this is the array index where x=0 */
#define sbf_c 5000
/*! 0.01 \n
If shapelet basis function is B(x), this is the x sampling resolution */
#define sbf_dx 0.01

/**
@brief Create the 1D shapelet basis function look-up array `sbf`

@details `sbf` should have `sbf_N * sbf_L * sizeof(float)` of memory allocated.

Each 1D basis function is \f$B_n(x, \beta)\f$, with
\f$n\f$ the basis function order, and \f$\beta \f$ the scaling factor. For all
stored arrays, I have set \f$\beta=1\f$. Pre-written arrays of
basis function orders \f$n=0\f$ to \f$n=100\f$ are then copied into `sbf`,
each of length `sbf_L`. Each samples over the range \f$-500\leq x \leq 500 \f$,
with a resolution \f$ \Delta x = 0.01 \f$

This single lookup array can then be used to generate 2D shapelet basis
functions later on varying \f$\beta\f$ values via a scaling factor.

@param[in] sbf An array with `sbf_N * sbf_L * sizeof(float)` of memory allocated
@returns The populated `sbf` array
*/
float * create_sbf(float *sbf);
