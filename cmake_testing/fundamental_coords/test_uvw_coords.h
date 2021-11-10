#include "woden_precision_defs.h"

typedef struct _uvw_settings_t {
  double *X_diff;
  double *Y_diff;
  double *Z_diff;
  user_precision_t *wavelengths;
  double *cha0s;
  double *sha0s;
  user_precision_t *us;
  user_precision_t *vs;
  user_precision_t *ws;
  user_precision_t *u_metres;
  user_precision_t *v_metres;
  user_precision_t *w_metres;
  double *lsts;  /*!< lst value for every visibility */

} uvw_settings_t;
