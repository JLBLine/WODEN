/*! \file woden_precision_defs.h
  A switch that sets two custom data types, `user_precision_t` and
  `user_precision_complex_t` to either be float or double. If the flag
  `-DDOUBLE_PRECISION` is added at compilation time, set to double, otherwise
  set to float
*/

#pragma once

#ifdef DOUBLE_PRECISION
  /*! If -DDOUBLE_PRECISION flag is added a compliation,
  then user_precision_t is set to double */
  typedef double user_precision_t;
  /*! If -DDOUBLE_PRECISION flag is added a compliation,
  then user_precision_complex_t is set to double _Complex */
  typedef double _Complex user_precision_complex_t;
#else
  /*! If -DDOUBLE_PRECISION flag is NOT added a compliation,
  then user_precision_t defaults to float */
  typedef float user_precision_t;
  /*! If -DDOUBLE_PRECISION flag is NOTE added a compliation,
  then user_precision_complex_t defaults to float _Complex */
  typedef float _Complex user_precision_complex_t;
#endif
