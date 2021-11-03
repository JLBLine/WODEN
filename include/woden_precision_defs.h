#pragma once

#ifdef DOUBLE_PRECISION
typedef double user_precision_t;
typedef double _Complex user_precision_complex_t;
#else
typedef float user_precision_t;
typedef float _Complex user_precision_complex_t;
#endif
