typedef struct _uvw_settings_t {
  float *X_diff;
  float *Y_diff;
  float *Z_diff;
  float *wavelengths;
  float *cha0s;
  float *sha0s;
  float *us;
  float *vs;
  float *ws;
  float *u_metres;
  float *v_metres;
  float *w_metres;
  float *lsts;  /*!< lst value for every visibility */

} uvw_settings_t;
