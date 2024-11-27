#include <stdio.h>
#include <stdlib.h>
// #include "gpu_macros.h"
#include <complex.h>
#include <math.h>
#include "constants.h"
#include "woden_precision_defs.h"


void calc_uvw_cpu(double *X_diff, double *Y_diff, double *Z_diff,
                  user_precision_t *u_metres, user_precision_t *v_metres, user_precision_t *w_metres,
                  user_precision_t *us, user_precision_t *vs, user_precision_t *ws,
                  user_precision_t *wavelengths,
                  double sdec0, double cdec0,
                  double *cha0s, double *sha0s,
                  int num_cross, int num_baselines, int num_times, int num_freqs) {

  user_precision_t u, v, w, wavelength;
  double sha0, cha0;

  for (int iBaseline = 0; iBaseline < num_cross; iBaseline++)
  {
    int mod_baseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);
    int time_ind = floorf(((float)iBaseline - (float)mod_baseline) / ((float)num_freqs*(float)num_baselines));
    int xyz_ind = time_ind*num_baselines + mod_baseline;

    sha0 = sha0s[iBaseline];
    cha0 = cha0s[iBaseline];

    u = (sha0*X_diff[xyz_ind]) + (cha0*Y_diff[xyz_ind]);
    v = -(sdec0*cha0*X_diff[xyz_ind]) + (sdec0*sha0*Y_diff[xyz_ind]) + (cdec0*Z_diff[xyz_ind]);
    w = (cdec0*cha0*X_diff[xyz_ind]) - (cdec0*sha0*Y_diff[xyz_ind]) + (sdec0*Z_diff[xyz_ind]);

    u_metres[iBaseline] = u;
    v_metres[iBaseline] = v;
    w_metres[iBaseline] = w;

    wavelength = wavelengths[iBaseline];

    us[iBaseline] = u / wavelength;
    vs[iBaseline] = v / wavelength;
    ws[iBaseline] = w / wavelength;

  }
}


void calc_lmn_cpu(double ra0, double sdec0, double cdec0, 
                  double *ras, double *decs,
                  double *ls, double *ms, double *ns, int num_components){

  for (int iComponent = 0; iComponent < num_components; iComponent++){

    double cdec;
    double sdec;
    double cdra;
    double sdra;

    double ra = ras[iComponent];
    double dec = decs[iComponent];

    cdec = cos(dec);
    sdec = sin(dec);
    cdra = cos((ra - ra0));
    sdra = sin((ra - ra0));

    ls[iComponent] = cdec*sdra;
    ms[iComponent] = sdec*cdec0 - cdec*sdec0*cdra;
    ns[iComponent] = sdec*sdec0 + cdec*cdec0*cdra;

  }
}

void set_auto_uvw_to_zero(int num_cross, int num_autos,
                          user_precision_t *us, user_precision_t *vs,
                          user_precision_t *ws) {

  for (int iAuto = 0; iAuto < num_autos; iAuto++) {
    us[num_cross + iAuto] = 0.0;
    vs[num_cross + iAuto] = 0.0;
    ws[num_cross + iAuto] = 0.0;
  }
}

void calc_uv_shapelet_cpu(double *X_diff, double *Y_diff, double *Z_diff,
      user_precision_t *u_shapes, user_precision_t *v_shapes,
      double *lsts, double *ras, double *decs,
      const int num_baselines, const int num_times, const int num_shapes) {

  for (int iBaseline = 0; iBaseline < num_baselines*num_times; iBaseline++){
    for (int iComponent = 0; iComponent < num_shapes; iComponent++){
      int mobaseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);
      int time_ind = floorf(((float)iBaseline - (float)mobaseline) / (float)num_baselines);
      int xyz_ind = time_ind*num_baselines + mobaseline;

      double sdec0 = sin(decs[iComponent]);
      double cdec0 = cos(decs[iComponent]);
      double sha0 = sin(lsts[time_ind] - ras[iComponent]);
      double cha0 = cos(lsts[time_ind] - ras[iComponent]);

      user_precision_t u_shape = (sha0*X_diff[xyz_ind]) + (cha0*Y_diff[xyz_ind]);
      user_precision_t v_shape = -(sdec0*cha0*X_diff[xyz_ind]) + (sdec0*sha0*Y_diff[xyz_ind]) + (cdec0*Z_diff[xyz_ind]);

      int stripe = num_baselines*num_times*iComponent + time_ind*num_baselines + mobaseline;

      u_shapes[stripe] = u_shape;
      v_shapes[stripe] = v_shape;
    }
  }
}

