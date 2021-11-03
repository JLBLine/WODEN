#include "woden_precision_defs.h"

typedef struct _uvw_settings_t {
  user_precision_t *X_diff;
  user_precision_t *Y_diff;
  user_precision_t *Z_diff;
  user_precision_t *wavelengths;
  user_precision_t *cha0s;
  user_precision_t *sha0s;
  user_precision_t *us;
  user_precision_t *vs;
  user_precision_t *ws;
  user_precision_t *u_metres;
  user_precision_t *v_metres;
  user_precision_t *w_metres;
  user_precision_t *lsts;  /*!< lst value for every visibility */

} uvw_settings_t;
