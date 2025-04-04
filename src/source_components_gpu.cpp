#include "source_components_gpu.h"

__device__  gpuUserComplex calc_measurement_equation_gpu(user_precision_t *d_us,
           user_precision_t *d_vs, user_precision_t *d_ws,
           double *d_ls, double *d_ms, double *d_ns,
           const int iBaseline, const int iComponent){

  gpuUserComplex visi;

  double u, v, w;
  double l, m, n;

  u = (double)d_us[iBaseline];
  v = (double)d_vs[iBaseline];
  w = (double)d_ws[iBaseline];

  // printf("%d u: %f, v: %f, w: %f\n", iComponent, u, v, w);

  // if (iBaseline == 0) {
  //   printf("l: %f\n", d_ls[iComponent]);
  // }


  l = d_ls[iComponent];
  m = d_ms[iComponent];
  n = d_ns[iComponent];

  //Not sure why, but get match with OSKAR/RTS sims, and correct location
  //on sky through WSClean, without negative infront on 2pi
  double temp = 2*M_PI*( u*l + v*m + w*(n-1) );

  visi.y = (user_precision_t)sin(temp);
  visi.x = (user_precision_t)cos(temp);

  return visi;
}

__device__ void apply_beam_gains_stokesIQUV_on_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY) {

  //Conjugate the second beam gains
  gpuUserComplex g2x_conj = make_gpuUserComplex(g2x.x,-g2x.y);
  gpuUserComplex D2x_conj = make_gpuUserComplex(D2x.x,-D2x.y);
  gpuUserComplex D2y_conj = make_gpuUserComplex(D2y.x,-D2y.y);
  gpuUserComplex g2y_conj = make_gpuUserComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  gpuUserComplex visi_I = make_gpuUserComplex(flux_I, 0.0)*visi_component;
  gpuUserComplex visi_Q = make_gpuUserComplex(flux_Q, 0.0)*visi_component;
  gpuUserComplex visi_U = make_gpuUserComplex(flux_U, 0.0)*visi_component;
  gpuUserComplex visi_V = make_gpuUserComplex(flux_V, 0.0)*visi_component;

  gpuUserComplex this_XX;
  gpuUserComplex this_XY;
  gpuUserComplex this_YX;
  gpuUserComplex this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XX += (g1x*g2x_conj - D1x*D2x_conj)*visi_Q;
  this_XX += (g1x*D2x_conj + D1x*g2x_conj)*visi_U;
  this_XX += (make_gpuUserComplex(0.0,1.0)*visi_V)*(g1x*D2x_conj - D1x*g2x_conj);

  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_XY += (g1x*D2y_conj - D1x*g2y_conj)*visi_Q;
  this_XY += (g1x*g2y_conj + D1x*D2y_conj)*visi_U;
  this_XY += (make_gpuUserComplex(0.0,1.0)*visi_V)*(g1x*g2y_conj - D1x*D2y_conj);

  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YX += (D1y*g2x_conj - g1y*D2x_conj)*visi_Q;
  this_YX += (D1y*D2x_conj + g1y*g2x_conj)*visi_U;
  this_YX += (make_gpuUserComplex(0.0,1.0)*visi_V)*(D1y*D2x_conj - g1y*g2x_conj);

  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;
  this_YY += (D1y*D2y_conj - g1y*g2y_conj)*visi_Q;
  this_YY += (D1y*g2y_conj + g1y*D2y_conj)*visi_U;
  this_YY += (make_gpuUserComplex(0.0,1.0)*visi_V)*(D1y*g2y_conj - g1y*D2y_conj);

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void apply_beam_gains_stokesIQUV_off_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I, user_precision_t flux_Q,
          user_precision_t flux_U, user_precision_t flux_V,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY) {

  //Conjugate the second beam gains
  gpuUserComplex g2x_conj = make_gpuUserComplex(g2x.x,-g2x.y);
  gpuUserComplex D2x_conj = make_gpuUserComplex(D2x.x,-D2x.y);
  gpuUserComplex D2y_conj = make_gpuUserComplex(D2y.x,-D2y.y);
  gpuUserComplex g2y_conj = make_gpuUserComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  gpuUserComplex visi_I = make_gpuUserComplex(flux_I, 0.0)*visi_component;
  gpuUserComplex visi_Q = make_gpuUserComplex(flux_Q, 0.0)*visi_component;
  gpuUserComplex visi_U = make_gpuUserComplex(flux_U, 0.0)*visi_component;
  gpuUserComplex visi_V = make_gpuUserComplex(flux_V, 0.0)*visi_component;

  gpuUserComplex this_XX;
  gpuUserComplex this_XY;
  gpuUserComplex this_YX;
  gpuUserComplex this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XX -= (g1x*D2x_conj + D1x*g2x_conj)*visi_Q;
  this_XX += (g1x*g2x_conj + D1x*D2x_conj)*visi_U;
  this_XX += (make_gpuUserComplex(0.0,1.0)*visi_V)*(g1x*D2x_conj - D1x*g2x_conj);

  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_XY -= (g1x*g2y_conj + D1x*D2y_conj)*visi_Q;
  this_XY += (g1x*D2y_conj - D1x*g2y_conj)*visi_U;
  this_XY += (make_gpuUserComplex(0.0,1.0)*visi_V)*(g1x*g2y_conj - D1x*D2y_conj);

  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YX -= (D1y*D2x_conj + g1y*g2x_conj)*visi_Q;
  this_YX += (D1y*g2x_conj - g1y*D2x_conj)*visi_U;
  this_YX += (make_gpuUserComplex(0.0,1.0)*visi_V)*(D1y*D2x_conj - g1y*g2x_conj);

  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;
  this_YY -= (D1y*g2y_conj + g1y*D2y_conj)*visi_Q;
  this_YY += (D1y*D2y_conj - g1y*g2y_conj)*visi_U;
  this_YY += (make_gpuUserComplex(0.0,1.0)*visi_V)*(D1y*g2y_conj - g1y*D2y_conj);

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void get_beam_gains_gpu(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
           gpuUserComplex * g1x, gpuUserComplex * D1x,
           gpuUserComplex * D1y, gpuUserComplex * g1y,
           gpuUserComplex * g2x, gpuUserComplex * D2x,
           gpuUserComplex * D2y, gpuUserComplex * g2y){

  int beam_ind = 0;
  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
  freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
  beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    //Set gains to one if no beam
  if (beamtype == NO_BEAM) {
    * g1x = make_gpuUserComplex(1.0, 0.0);
    * g2x = make_gpuUserComplex(1.0, 0.0);
    * g1y = make_gpuUserComplex(1.0, 0.0);
    * g2y = make_gpuUserComplex(1.0, 0.0);
  }

  //Get gains if using a beam
  else {
    * g1x = d_gxs[beam_ind];
    * g2x = d_gxs[beam_ind];
    * g1y = d_gys[beam_ind];
    * g2y = d_gys[beam_ind];

  }

  // Set leakage to zero if no leakage
  if (beamtype == NO_BEAM || beamtype == GAUSS_BEAM || beamtype == ANALY_DIPOLE) {
    * D1x = make_gpuUserComplex(0.0, 0.0);
    * D2x = make_gpuUserComplex(0.0, 0.0);
    * D1y = make_gpuUserComplex(0.0, 0.0);
    * D2y = make_gpuUserComplex(0.0, 0.0);
  }
  
  else {
    * D1x = d_Dxs[beam_ind];
    * D2x = d_Dxs[beam_ind];
    * D1y = d_Dys[beam_ind];
    * D2y = d_Dys[beam_ind];
  }

} //end __device__ get_beam_gains_gpu


__device__ void get_beam_gains_multibeams_gpu(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
           gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
           int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map,
           gpuUserComplex * g1x, gpuUserComplex * D1x,
           gpuUserComplex * D1y, gpuUserComplex * g1y,
           gpuUserComplex * g2x, gpuUserComplex * D2x,
           gpuUserComplex * D2y, gpuUserComplex * g2y){

  int time_ind = 0;
  int freq_ind = 0;

  time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
  freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
  // beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

  int baseline_ind = iBaseline % num_baselines;

  int ant1_ind = d_ant1_to_baseline_map[baseline_ind];
  int ant2_ind = d_ant2_to_baseline_map[baseline_ind];

  int beam1 = ant1_ind*num_freqs*num_components*num_times + num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

  int beam2 = ant2_ind*num_freqs*num_components*num_times + num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    //Set gains to one if no beam
  if (beamtype == NO_BEAM) {
    * g1x = make_gpuUserComplex(1.0, 0.0);
    * g2x = make_gpuUserComplex(1.0, 0.0);
    * g1y = make_gpuUserComplex(1.0, 0.0);
    * g2y = make_gpuUserComplex(1.0, 0.0);
  }

  //Get gains if using a beam
  else {
    * g1x = d_gxs[beam1];
    * g2x = d_gxs[beam2];
    * g1y = d_gys[beam1];
    * g2y = d_gys[beam2];

  }

  //Only MWA models have leakge terms at the moment
  if (beamtype == NO_BEAM || beamtype == GAUSS_BEAM || beamtype == ANALY_DIPOLE) {
    * D1x = make_gpuUserComplex(0.0, 0.0);
    * D2x = make_gpuUserComplex(0.0, 0.0);
    * D1y = make_gpuUserComplex(0.0, 0.0);
    * D2y = make_gpuUserComplex(0.0, 0.0);
  }
  // Set leakage to zero if no leakage
  else {
    * D1x = d_Dxs[beam1];
    * D2x = d_Dxs[beam2];
    * D1y = d_Dys[beam1];
    * D2y = d_Dys[beam2];
  }

} //end __device__ get_beam_gains_multibeams_gpu

__device__ void apply_beam_gains_stokesI_on_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY) {

  //Conjugate the second beam gains
  gpuUserComplex g2x_conj = make_gpuUserComplex(g2x.x,-g2x.y);
  gpuUserComplex D2x_conj = make_gpuUserComplex(D2x.x,-D2x.y);
  gpuUserComplex D2y_conj = make_gpuUserComplex(D2y.x,-D2y.y);
  gpuUserComplex g2y_conj = make_gpuUserComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  gpuUserComplex visi_I = make_gpuUserComplex(flux_I, 0.0)*visi_component;

  gpuUserComplex this_XX;
  gpuUserComplex this_XY;
  gpuUserComplex this_YX;
  gpuUserComplex this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void apply_beam_gains_stokesI_off_cardinal_gpu(gpuUserComplex g1x, gpuUserComplex D1x,
          gpuUserComplex D1y, gpuUserComplex g1y,
          gpuUserComplex g2x, gpuUserComplex D2x,
          gpuUserComplex D2y, gpuUserComplex g2y,
          user_precision_t flux_I,
          gpuUserComplex visi_component,
          gpuUserComplex * visi_XX, gpuUserComplex * visi_XY,
          gpuUserComplex * visi_YX, gpuUserComplex * visi_YY) {

  //Conjugate the second beam gains
  gpuUserComplex g2x_conj = make_gpuUserComplex(g2x.x,-g2x.y);
  gpuUserComplex D2x_conj = make_gpuUserComplex(D2x.x,-D2x.y);
  gpuUserComplex D2y_conj = make_gpuUserComplex(D2y.x,-D2y.y);
  gpuUserComplex g2y_conj = make_gpuUserComplex(g2y.x,-g2y.y);

  //Create the Stokes visibilities
  gpuUserComplex visi_I = make_gpuUserComplex(flux_I, 0.0)*visi_component;

  gpuUserComplex this_XX;
  gpuUserComplex this_XY;
  gpuUserComplex this_YX;
  gpuUserComplex this_YY;

  this_XX = (g1x*g2x_conj + D1x*D2x_conj)*visi_I;
  this_XY = (g1x*D2y_conj + D1x*g2y_conj)*visi_I;
  this_YX = (D1y*g2x_conj + g1y*D2x_conj)*visi_I;
  this_YY = (D1y*D2y_conj + g1y*g2y_conj)*visi_I;

  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;

}

__device__ void update_sum_visis_stokesIQUV_gpu(int iBaseline, int iComponent,
    int num_freqs, int num_baselines, int num_components, int num_times,
    int beamtype, int off_cardinal_dipoles,
    gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
    gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
    int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
    gpuUserComplex visi_component,
    user_precision_t flux_I, user_precision_t flux_Q,
    user_precision_t flux_U, user_precision_t flux_V,
    user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
    user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
    user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
    user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag){

    gpuUserComplex g1x;
    gpuUserComplex D1x;
    gpuUserComplex D1y;
    gpuUserComplex g1y;
    gpuUserComplex g2x;
    gpuUserComplex D2x;
    gpuUserComplex D2y;
    gpuUserComplex g2y;

    if (use_twobeams == 1){
      get_beam_gains_multibeams_gpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               d_ant1_to_baseline_map, d_ant2_to_baseline_map,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }
    else {
      get_beam_gains_gpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }

    // if (iBaseline == 0){
    //   printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g1x.x, g1x.y, D1x.x, D1x.y, D1y.x, D1y.y, g1y.x, g1y.y);
    //   printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g2x.x, g2x.y, D2x.x, D2x.y, D2y.x, D2y.y, g2y.x, g2y.y);
    // }

    gpuUserComplex visi_XX;
    gpuUserComplex visi_XY;
    gpuUserComplex visi_YX;
    gpuUserComplex visi_YY;

    // printf("iComponent IQUV: %d %f %f %f %f\n", iComponent, flux_I, flux_Q, flux_U, flux_V);

    if (off_cardinal_dipoles == 1) {
      apply_beam_gains_stokesIQUV_off_cardinal_gpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I, flux_Q, flux_U, flux_V,
                    visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    } else {
      apply_beam_gains_stokesIQUV_on_cardinal_gpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                    flux_I, flux_Q, flux_U, flux_V,
                    visi_component, &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    }

    // if (iBaseline == 0){
    //   printf("Visibilities: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", visi_XX.x, visi_XX.y, visi_XY.x, visi_XY.y,
    //     visi_YX.x, visi_YX.y, visi_YY.x, visi_YY.y);
    // }

    d_sum_visi_XX_real[iBaseline] += visi_XX.x;
    d_sum_visi_XX_imag[iBaseline] += visi_XX.y;

    d_sum_visi_XY_real[iBaseline] += visi_XY.x;
    d_sum_visi_XY_imag[iBaseline] += visi_XY.y;

    d_sum_visi_YX_real[iBaseline] += visi_YX.x;
    d_sum_visi_YX_imag[iBaseline] += visi_YX.y;

    d_sum_visi_YY_real[iBaseline] += visi_YY.x;
    d_sum_visi_YY_imag[iBaseline] += visi_YY.y;

    // if (iBaseline == 0){
    //   printf("Visibilities: %.3e %.3e %.3e %.3e\n", visi_XX.x, visi_XX.y,
    //                         d_sum_visi_XX_real[iBaseline], d_sum_visi_XX_imag[iBaseline]);
    // }
}

__device__ void update_sum_visis_stokesI_gpu(int iBaseline, int iComponent,
    int num_freqs, int num_baselines, int num_components, int num_times,
    int beamtype, int off_cardinal_dipoles,
    gpuUserComplex *d_gxs, gpuUserComplex *d_Dxs,
    gpuUserComplex *d_Dys, gpuUserComplex *d_gys,
    int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
    gpuUserComplex visi_component,
    user_precision_t flux_I,
    user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
    user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
    user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
    user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag){

    gpuUserComplex g1x;
    gpuUserComplex D1x;
    gpuUserComplex D1y;
    gpuUserComplex g1y;
    gpuUserComplex g2x;
    gpuUserComplex D2x;
    gpuUserComplex D2y;
    gpuUserComplex g2y;

    if (use_twobeams == 1){
      get_beam_gains_multibeams_gpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               d_ant1_to_baseline_map, d_ant2_to_baseline_map,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }
    else {
      get_beam_gains_gpu(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gxs, d_Dxs,
               d_Dys, d_gys,
               &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
    }

    gpuUserComplex visi_XX;
    gpuUserComplex visi_XY;
    gpuUserComplex visi_YX;
    gpuUserComplex visi_YY;

    if (off_cardinal_dipoles == 1) {
      apply_beam_gains_stokesI_off_cardinal_gpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                                           flux_I, visi_component,
                                           &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    } else {
      apply_beam_gains_stokesI_on_cardinal_gpu(g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y,
                                           flux_I, visi_component,
                                           &visi_XX, &visi_XY, &visi_YX, &visi_YY);
    }

    // if (iBaseline == 0) {
    //   printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g1x.x, g1x.y, D1x.x, D1x.y, D1y.x, D1y.y, g1y.x, g1y.y);
    //   // printf("Beam gains: %.12f %.12f %.12f %.12f %.12f %.12f %.12f %.12f\n", g2x.x, g2x.y, D2x.x, D2x.y, D2y.x, D2y.y, g2y.x, g2y.y);
    //   printf("Other bits %.4f %.4f %.4f\n", flux_I, visi_component.x, visi_component.y);
    // }

    d_sum_visi_XX_real[iBaseline] += visi_XX.x;
    d_sum_visi_XX_imag[iBaseline] += visi_XX.y;

    d_sum_visi_XY_real[iBaseline] += visi_XY.x;
    d_sum_visi_XY_imag[iBaseline] += visi_XY.y;

    d_sum_visi_YX_real[iBaseline] += visi_YX.x;
    d_sum_visi_YX_imag[iBaseline] += visi_YX.y;

    d_sum_visi_YY_real[iBaseline] += visi_YY.x;
    d_sum_visi_YY_imag[iBaseline] += visi_YY.y;

}

//just make things zero pls
__global__ void kern_make_zeros_user_precision(user_precision_t *array, int num_arr) {

  const int iComp = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iComp < num_arr)
  {
    array[iComp] = 0.0;
  }
}


//Allocate space for the extrapolated Stokes parameters
extern "C" void malloc_extrapolated_flux_arrays_gpu(components_t *d_components, int num_comps,
                                     int num_freqs){
  d_components->extrap_stokesI = NULL;
  gpuMalloc( (void**)&d_components->extrap_stokesI,
                                   num_comps*num_freqs*sizeof(user_precision_t) );

  if (d_components->do_QUV == 1)
  {
      // printf("Doing full polarisation\n");
      d_components->extrap_stokesQ = NULL;
      gpuMalloc( (void**)&d_components->extrap_stokesQ,
                                      num_comps*num_freqs*sizeof(user_precision_t) );
      d_components->extrap_stokesU = NULL;
      gpuMalloc( (void**)&d_components->extrap_stokesU,
                                      num_comps*num_freqs*sizeof(user_precision_t) );
      d_components->extrap_stokesV = NULL;
      gpuMalloc( (void**)&d_components->extrap_stokesV,
                                      num_comps*num_freqs*sizeof(user_precision_t) );

      //set everything to zero at not every component has to have full polarisation

      dim3 grid, threads;

      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)(num_comps*num_freqs) / (float)threads.x );
      grid.y = 1;

      gpuErrorCheckKernel("kern_make_zeros_user_precision",
              kern_make_zeros_user_precision, grid, threads,
              d_components->extrap_stokesQ, num_comps*num_freqs);

      gpuErrorCheckKernel("kern_make_zeros_user_precision",
              kern_make_zeros_user_precision, grid, threads,
              d_components->extrap_stokesU, num_comps*num_freqs);

      gpuErrorCheckKernel("kern_make_zeros_user_precision",
              kern_make_zeros_user_precision, grid, threads,
              d_components->extrap_stokesV, num_comps*num_freqs);

  }
}

__device__ void extrap_stokes_power_law_gpu(user_precision_t *d_ref_fluxes,
           user_precision_t *d_SIs,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux){

  double d_freq = d_extrap_freqs[iFreq];
  user_precision_t flux_ratio = pow(d_freq / REF_FREQ, d_SIs[iFluxComp]);

  * extrap_flux = d_ref_fluxes[iFluxComp] * flux_ratio;
}

__global__ void kern_extrap_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_I;

    extrap_stokes_power_law_gpu(d_components.power_ref_stokesI,
                            d_components.power_SIs,
                            d_extrap_freqs,
                            iFluxComp, iFreq, &flux_I);

    int iComponent = d_components.power_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesI[extrap_ind] = flux_I;
  }
}

__global__ void kern_extrap_power_laws_stokesV(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_V;

    extrap_stokes_power_law_gpu(d_components.stokesV_power_ref_flux,
                            d_components.stokesV_power_SIs,
                            d_extrap_freqs,
                            iFluxComp, iFreq, &flux_V);

    int iComponent = d_components.stokesV_power_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesV[extrap_ind] = flux_V;
  }
}

__global__ void kern_extrap_power_laws_linpol(int num_extrap_freqs, double *d_extrap_freqs,
                                       int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_linpol;

    extrap_stokes_power_law_gpu(d_components.linpol_power_ref_flux,
                            d_components.linpol_power_SIs,
                            d_extrap_freqs,
                            iFluxComp, iFreq, &flux_linpol);

    int iComponent = d_components.linpol_power_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesQ[extrap_ind] = flux_linpol;
  }
}

__device__ void extrap_stokes_curved_power_law_gpu(user_precision_t *d_ref_fluxes,
           user_precision_t *d_SIs, user_precision_t *d_qs,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux){

  double d_freq = d_extrap_freqs[iFreq];

  user_precision_t si_ratio = pow(d_freq / REF_FREQ, d_SIs[iFluxComp]);

  double log_freq_ratio = log(d_freq / REF_FREQ);

  double q = (double)d_qs[iFluxComp];
  double exp_bit = exp(q*log_freq_ratio*log_freq_ratio);

  user_precision_t flux_ratio = si_ratio * exp_bit;

  * extrap_flux = d_ref_fluxes[iFluxComp] * flux_ratio;
}

__global__ void kern_extrap_curved_power_laws_stokesI(int num_extrap_freqs, double *d_extrap_freqs,
                                              int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_I;

    extrap_stokes_curved_power_law_gpu(d_components.curve_ref_stokesI,
                            d_components.curve_SIs, d_components.curve_qs, 
                            d_extrap_freqs,
                            iFluxComp, iFreq,
                            &flux_I);

    int iComponent = d_components.curve_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesI[extrap_ind] = flux_I;

  }
}

__global__ void kern_extrap_curved_power_laws_stokesV(int num_extrap_freqs, double *d_extrap_freqs,
                                              int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_V;

    extrap_stokes_curved_power_law_gpu(d_components.stokesV_curve_ref_flux,
                            d_components.stokesV_curve_SIs,
                            d_components.stokesV_curve_qs, 
                            d_extrap_freqs, iFluxComp, iFreq, &flux_V);

    int iComponent = d_components.stokesV_curve_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesV[extrap_ind] = flux_V;

  }
}

__global__ void kern_extrap_curved_power_laws_linpol(int num_extrap_freqs, double *d_extrap_freqs,
                                              int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t flux_linpol;

    extrap_stokes_curved_power_law_gpu(d_components.linpol_curve_ref_flux,
                            d_components.linpol_curve_SIs,
                            d_components.linpol_curve_qs, 
                            d_extrap_freqs, iFluxComp, iFreq, &flux_linpol);

    int iComponent = d_components.linpol_curve_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesQ[extrap_ind] = flux_linpol;

  }
}

__global__ void kern_polarisation_fraction_stokesV(int num_extrap_freqs, 
             double *d_extrap_freqs, int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t pol_frac;
    pol_frac = d_components.stokesV_pol_fracs[iFluxComp];

    int iComponent = d_components.stokesV_pol_frac_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesV[extrap_ind] = pol_frac*d_components.extrap_stokesI[extrap_ind];
  }
}

__global__ void kern_polarisation_fraction_linpol(int num_extrap_freqs, 
             double *d_extrap_freqs, int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t pol_frac;
    pol_frac = d_components.linpol_pol_fracs[iFluxComp];

    int iComponent = d_components.linpol_pol_frac_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    d_components.extrap_stokesQ[extrap_ind] = pol_frac*d_components.extrap_stokesI[extrap_ind];
  }
}

__device__ user_precision_t calc_gradient_extrap_list_gpu(user_precision_t *list_fluxes,
          double *list_freqs, double desired_freq, int low_ind_1, int low_ind_2) {

  user_precision_t gradient;
  user_precision_t extrap_flux;

  //If both zero, just stick to zero
  if (list_fluxes[low_ind_1] == 0 && list_fluxes[low_ind_2] == 0) {
   extrap_flux = 0.0;
  }

  //If one is negative, do interpolation in linear space
  else if (list_fluxes[low_ind_1] <= 0 || list_fluxes[low_ind_2] <= 0) {
    gradient = (list_fluxes[low_ind_2] - list_fluxes[low_ind_1]) / (list_freqs[low_ind_2] - list_freqs[low_ind_1]);
    extrap_flux = list_fluxes[low_ind_1] + gradient*(desired_freq - list_freqs[low_ind_1]);
  }

  else {

    user_precision_t logflux1, logflux2, logfreq1, logfreq2, log_des_freq;

    logflux1 = log10(list_fluxes[low_ind_1]);
    logflux2 = log10(list_fluxes[low_ind_2]);
    logfreq1 = log10(list_freqs[low_ind_1]);
    logfreq2 = log10(list_freqs[low_ind_2]);
    log_des_freq = log10(desired_freq);

    gradient = (logflux2 - logflux1) / (logfreq2 - logfreq1);
    extrap_flux = logflux1 + gradient*(log_des_freq - logfreq1);

    extrap_flux = pow(10, extrap_flux);

  }
  return extrap_flux;
}


__device__ void extrap_stokes_list_fluxes_gpu(user_precision_t *list_stokes,
           double *list_freqs, int *arr_num_list_values, int *list_start_indexes,
           double *d_extrap_freqs, int iFluxComp, int iFreq,
           user_precision_t * extrap_flux){

  int num_list_values = arr_num_list_values[iFluxComp];
  int list_start_ind = list_start_indexes[iFluxComp];

  double d_extrap_freq = d_extrap_freqs[iFreq];

  int low_ind_1 = -1;
  int low_ind_2 = -1;

  double low_val_1 = 1e16;
  // double low_val_2 = 1e16;

  double ref_freq;
  double abs_diff_freq;

  if (num_list_values == 1) {
    * extrap_flux = list_stokes[list_start_ind];
    return;
  }

  //First loop finds the absolute closest frequency
  for (int i = 0; i < num_list_values; i++) {
    ref_freq = list_freqs[list_start_ind + i];
    abs_diff_freq = abs(ref_freq - d_extrap_freq);

    if (abs_diff_freq < low_val_1) {
      low_val_1 = abs_diff_freq;
      low_ind_1 = i;
    }
  }

  //Depending on the closest frequency, we either want to search above or
  //below the target frequency to find points either side of the target freq

  //We happen to need the reference frequency; just return the refs
  if (list_freqs[list_start_ind + low_ind_1] == d_extrap_freq) {
    * extrap_flux = list_stokes[list_start_ind + low_ind_1];
    return;
  }
  else {
    //The closest freq is the first index, so set the second index to the second
    if (low_ind_1 == 0) {
      low_ind_2 = 1;
    }
    //closest freq the highest list entry - set second index to one below
    //(order of indexes doesn't matter, as the calculated gradient is pos/neg
    //as needed)
    else if (low_ind_1 == num_list_values - 1){
      low_ind_2 = low_ind_1 - 1;
    }
    else {
      //closest freq is higher than desired - set second index to one below
      //(order of indexes doesn't matter, as the calculated gradient is pos/neg
      //as needed)
      if (list_freqs[list_start_ind + low_ind_1] > d_extrap_freq){
        low_ind_2 = low_ind_1 - 1;
      }
      else {
        low_ind_2 = low_ind_1 + 1;
      }
        //We are extrapolating to a frequency that is lower than all list entries
        //so just stick low_ind_2 to one above low_ind_1
    }
  }

  * extrap_flux = calc_gradient_extrap_list_gpu(list_stokes,
            list_freqs, d_extrap_freq,
            list_start_ind + low_ind_1, list_start_ind + low_ind_2);

  if (low_ind_2 == -1){

    printf("wrong range %.3e %.3e iFreq %d %.3e low %d %.3e\n", list_freqs[list_start_ind],
    list_freqs[list_start_ind + num_list_values-1],
    iFreq, d_extrap_freq,
    low_ind_1, list_freqs[list_start_ind + low_ind_1]);
    printf("The flooxes %.3e \n",* extrap_flux);
  }
}

__global__ void kern_extrap_list_fluxes(user_precision_t *list_stokes, double *list_freqs,
                                        int *num_list_values, int *list_start_indexes,
                                        int *list_comp_inds,
                                        int num_extrap_freqs, double *d_extrap_freqs,
                                        int num_comps, user_precision_t *extrap_stokes) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    user_precision_t extrap_flux;

    int iComponent = list_comp_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    extrap_stokes_list_fluxes_gpu(list_stokes, list_freqs, num_list_values,
                              list_start_indexes, d_extrap_freqs,
                              iFluxComp, iFreq,
                              &extrap_flux);

    extrap_stokes[extrap_ind] = extrap_flux;

  }
}

__global__ void kern_apply_rotation_measure(int num_extrap_freqs, double *d_extrap_freqs,
                                            int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iFluxComp = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iFreq = threadIdx.y + (blockDim.y*blockIdx.y);
  if(iFluxComp < num_comps && iFreq < num_extrap_freqs) {

    int iComponent = d_components.linpol_angle_inds[iFluxComp];
    int extrap_ind = num_extrap_freqs*iComponent + iFreq;

    user_precision_t rm = d_components.rm_values[iFluxComp];
    user_precision_t intr_pol_angle = d_components.intr_pol_angle[iFluxComp];

    //We should have calculated the linear polarisation flux before, and
    //shoved it in the Stokes Q extrap array. No point in wasting precious
    //GPU memory as we'll immediate overwrite it
    user_precision_t linpol_flux = d_components.extrap_stokesQ[extrap_ind];

    double wavelength = VELC / d_extrap_freqs[iFreq];
    double angle = 2*(intr_pol_angle + rm*wavelength*wavelength);

    d_components.extrap_stokesQ[extrap_ind] = linpol_flux*cos(angle);
    d_components.extrap_stokesU[extrap_ind] = linpol_flux*sin(angle);

  }
}

__global__ void kern_print_extrap_fluxes(int freq_ind, int num_extrap_freqs,
                                         int num_comps, components_t d_components) {

  // Start by computing which baseline we're going to do
  const int iComp = threadIdx.x + (blockDim.x*blockIdx.x);
  if(iComp < num_comps) {

    int extrap_ind = num_extrap_freqs*iComp + freq_ind;

    printf("iComp: %d %.5f %.5f %.5f %.5f\n", iComp, d_components.extrap_stokesI[extrap_ind],
                                                     d_components.extrap_stokesQ[extrap_ind],
                                                     d_components.extrap_stokesU[extrap_ind],
                                                     d_components.extrap_stokesV[extrap_ind]);

    

  }
}

//wrap all flux extrap functions so we can call them from `source_components_common.c`
//I've done this so if someone is adding funnctionality to `source_components_common.c`,
//it should be clear where you need to create both a GPU and CPU version of a 
//particular function
extern "C" void extrap_power_laws_stokesI_gpu(components_t d_components,
                                   int n_powers, double *d_extrap_freqs,
                                   int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.x = (int)ceilf( (float)n_powers / (float)threads.x );
  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );

  gpuErrorCheckKernel("kern_extrap_power_laws_stokesI",
                        kern_extrap_power_laws_stokesI, grid, threads,
                        num_extrap_freqs, d_extrap_freqs,
                        n_powers, d_components);
}

extern "C" void extrap_curved_power_laws_stokesI_gpu(components_t d_components,
                                   int n_curves, double *d_extrap_freqs,
                                   int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)n_curves / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_curved_power_laws_stokesI",
                    kern_extrap_curved_power_laws_stokesI, grid, threads,
                    num_extrap_freqs, d_extrap_freqs,
                    n_curves, d_components);
}

extern "C" void extrap_list_fluxes_stokesI_gpu(components_t d_components,
                                   int n_lists, double *d_extrap_freqs,
                                   int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)n_lists / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_list_fluxes",
                       kern_extrap_list_fluxes, grid, threads,
                       d_components.list_stokesI, d_components.list_freqs,
                       d_components.num_list_values, d_components.list_start_indexes,
                       d_components.list_comp_inds,
                       num_extrap_freqs, d_extrap_freqs,
                       n_lists, d_components.extrap_stokesI);
}

extern "C" void extrap_power_laws_stokesV_gpu(components_t d_components,
                                              double *d_extrap_freqs,
                                              int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_stokesV_power / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_power_laws_stokesV",
                       kern_extrap_power_laws_stokesV, grid, threads,
                       num_extrap_freqs, d_extrap_freqs,
                       d_components.n_stokesV_power, d_components);
}

extern "C" void extrap_curved_power_laws_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_stokesV_curve / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_curved_power_laws_stokesV",
                          kern_extrap_curved_power_laws_stokesV, grid, threads,
                          num_extrap_freqs, d_extrap_freqs,
                          d_components.n_stokesV_curve, d_components);
}

extern "C" void polarisation_fraction_stokesV_gpu(components_t d_components,
                                                     double *d_extrap_freqs,
                                                     int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_stokesV_pol_frac / (float)threads.x );
  gpuErrorCheckKernel("kern_polarisation_fraction_stokesV",
                          kern_polarisation_fraction_stokesV, grid, threads,
                          num_extrap_freqs, d_extrap_freqs,
                          d_components.n_stokesV_pol_frac, d_components);
}

extern "C" void extrap_list_fluxes_stokesV_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_stokesV_list / (float)threads.x );
  // printf("Extrapolating Stokes V list fluxes\n");
  gpuErrorCheckKernel("kern_extrap_list_fluxes",
                      kern_extrap_list_fluxes, grid, threads,
                      d_components.stokesV_list_ref_flux, d_components.stokesV_list_ref_freqs,
                      d_components.stokesV_num_list_values, d_components.stokesV_list_start_indexes,
                      d_components.stokesV_list_comp_inds,
                      num_extrap_freqs, d_extrap_freqs,
                      d_components.n_stokesV_list, d_components.extrap_stokesV);
}

extern "C" void extrap_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_linpol_power / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_power_laws_linpol",
                        kern_extrap_power_laws_linpol, grid, threads,
                        num_extrap_freqs, d_extrap_freqs,
                        d_components.n_linpol_power, d_components);
}

extern "C" void extrap_curved_power_laws_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_linpol_curve / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_curved_power_laws_linpol",
                        kern_extrap_curved_power_laws_linpol, grid, threads,
                        num_extrap_freqs, d_extrap_freqs,
                        d_components.n_linpol_curve, d_components);
}

extern "C" void polarisation_fraction_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_linpol_pol_frac / (float)threads.x );
  gpuErrorCheckKernel("kern_polarisation_fraction_linpol",
                        kern_polarisation_fraction_linpol, grid, threads,
                        num_extrap_freqs, d_extrap_freqs,
                        d_components.n_linpol_pol_frac, d_components);
}

extern "C" void extrap_list_fluxes_linpol_gpu(components_t d_components,
                                               double *d_extrap_freqs,
                                               int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_linpol_list / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_list_fluxes",
                      kern_extrap_list_fluxes, grid, threads,
                      d_components.stokesQ_list_ref_flux, d_components.stokesQ_list_ref_freqs,
                      d_components.stokesQ_num_list_values, d_components.stokesQ_list_start_indexes,
                      d_components.stokesQ_list_comp_inds,
                      num_extrap_freqs, d_extrap_freqs,
                      d_components.n_linpol_list, d_components.extrap_stokesQ);

  grid.x = (int)ceilf( (float)d_components.n_linpol_list / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_list_fluxes",
                      kern_extrap_list_fluxes, grid, threads,
                      d_components.stokesU_list_ref_flux, d_components.stokesU_list_ref_freqs,
                      d_components.stokesU_num_list_values, d_components.stokesU_list_start_indexes,
                      d_components.stokesU_list_comp_inds,
                      num_extrap_freqs, d_extrap_freqs,
                      d_components.n_linpol_list, d_components.extrap_stokesU);
}

extern "C" void extrap_p_list_fluxes_linpol_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_linpol_p_list / (float)threads.x );
  gpuErrorCheckKernel("kern_extrap_list_fluxes",
                      kern_extrap_list_fluxes, grid, threads,
                      d_components.linpol_p_list_ref_flux, d_components.linpol_p_list_ref_freqs,
                      d_components.linpol_p_num_list_values, d_components.linpol_p_list_start_indexes,
                      d_components.linpol_p_list_comp_inds,
                      num_extrap_freqs, d_extrap_freqs,
                      d_components.n_linpol_p_list, d_components.extrap_stokesQ);
}

extern "C" void apply_rotation_measure_gpu(components_t d_components,
                                                double *d_extrap_freqs,
                                                int num_extrap_freqs){
  dim3 grid, threads;
  threads.x = 16;
  threads.y = 16;

  grid.y = (int)ceilf( (float)num_extrap_freqs / (float)threads.y );
  grid.x = (int)ceilf( (float)d_components.n_linpol_angles / (float)threads.x );
  gpuErrorCheckKernel("kern_apply_rotation_measure",
                        kern_apply_rotation_measure, grid, threads,
                        num_extrap_freqs, d_extrap_freqs,
                        d_components.n_linpol_angles, d_components);
}

extern "C" void malloc_beam_gains_gpu(beam_gains_t *d_component_beam_gains,
                                     int beamtype, int num_gains){

  //If we're using an everybeam model, all memory and values have already
  //been copied to GPU, so no need to allocate here
  if (beamtype == FEE_BEAM || beamtype == MWA_ANALY || beamtype == FEE_BEAM_INTERP
      || beamtype == GAUSS_BEAM || beamtype == ANALY_DIPOLE || beamtype == NO_BEAM) {

    //Only some models would have had leakage terms malloced
    if (beamtype == FEE_BEAM || beamtype == MWA_ANALY || beamtype == FEE_BEAM_INTERP) {
      gpuMalloc( (void**)&d_component_beam_gains->Dxs,
                      num_gains*sizeof(user_precision_complex_t) );
      gpuMalloc( (void**)&d_component_beam_gains->Dys,
                      num_gains*sizeof(user_precision_complex_t) );
    }
    gpuMalloc( (void**)&d_component_beam_gains->gxs,
                      num_gains*sizeof(user_precision_complex_t) );
    gpuMalloc( (void**)&d_component_beam_gains->gys,
                      num_gains*sizeof(user_precision_complex_t) );

  }
}





__global__ void kern_calc_visi_point_or_gauss(components_t d_components,
           beam_gains_t d_component_beam_gains,
           user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
           user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
           user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
           user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
           user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
           int num_components, int num_baselines, int num_freqs, int num_cross,
           int num_times, e_beamtype beamtype, e_component_type comptype,
           int off_cardinal_dipoles) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  if(iBaseline < num_cross) {

    int use_twobeams = d_component_beam_gains.use_twobeams;

    user_precision_t flux_I;
    user_precision_t flux_Q;
    user_precision_t flux_U;
    user_precision_t flux_V;

    gpuUserComplex visi_comp;
    gpuUserComplex V_envelop;

    user_precision_t pa, sinpa, cospa, u, v, x, y, invsig_x, invsig_y;

    //Find out what time and freq index this baseline corresponds to
    int time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    int freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);

    for (int iComponent = 0; iComponent < num_components; iComponent++) {
      int extrap_ind = num_freqs*iComponent + freq_ind;
      // printf("INSIDE KERN d_components.extrap_stokesI %p\n", d_components.extrap_stokesI);

      flux_I = d_components.extrap_stokesI[extrap_ind];
      if (d_components.do_QUV == 1) {
        flux_Q = d_components.extrap_stokesQ[extrap_ind];
        flux_U = d_components.extrap_stokesU[extrap_ind];
        flux_V = d_components.extrap_stokesV[extrap_ind];
      }

      // if (iBaseline == 0 && d_components.do_QUV == 1) {
      //   printf("Fluxes %.3e %.3e %.3e %.3e\n", flux_I, flux_Q, flux_U, flux_V);
      // }
      
      visi_comp = calc_measurement_equation_gpu(d_us, d_vs, d_ws,
                             d_components.ls, d_components.ms, d_components.ns,
                             iBaseline, iComponent);

      // if (iBaseline == 0 && d_components.do_QUV == 1) {
      //   printf("visi_comp %.3e %.3e \n", visi_comp.x, visi_comp.y);
      // }
      // printf("iComponent %d d_components.ls[iComponent] %f \n", iComponent, d_components.ls[iComponent]);

      if (comptype == GAUSSIAN) {

        V_envelop = make_gpuUserComplex( 1.0, 0.0 );

        pa = d_components.pas[iComponent];
        sinpa = sin(pa);
        cospa = cos(pa);
        u = d_us[iBaseline];
        v = d_vs[iBaseline];

        x =  cospa*v + sinpa*u; // major axis
        y = -sinpa*v + cospa*u; // minor axis
        invsig_x = d_components.majors[iComponent];
        invsig_y = d_components.minors[iComponent];

        V_envelop = make_gpuUserComplex( exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ), 0.0 );

        visi_comp = visi_comp*V_envelop;
      }

      if (d_components.do_QUV == 1)
      {
        update_sum_visis_stokesIQUV_gpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype, off_cardinal_dipoles,
             (gpuUserComplex *)d_component_beam_gains.gxs,
             (gpuUserComplex *)d_component_beam_gains.Dxs,
             (gpuUserComplex *)d_component_beam_gains.Dys,
             (gpuUserComplex *)d_component_beam_gains.gys,
             d_component_beam_gains.ant1_to_baseline_map,
             d_component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_comp, flux_I, flux_Q, flux_U, flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      } else {
        update_sum_visis_stokesI_gpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype, off_cardinal_dipoles,
             (gpuUserComplex *)d_component_beam_gains.gxs,
             (gpuUserComplex *)d_component_beam_gains.Dxs,
             (gpuUserComplex *)d_component_beam_gains.Dys,
             (gpuUserComplex *)d_component_beam_gains.gys,
             d_component_beam_gains.ant1_to_baseline_map,
             d_component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_comp, flux_I,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      }
    }
  }
}


extern "C" void calc_visi_point_or_gauss_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_components, e_beamtype beamtype,
                                        e_component_type comptype,
                                        woden_settings_t *woden_settings){

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)woden_settings->num_cross / (float)threads.x );

  gpuErrorCheckKernel("kern_calc_visi_point_or_gauss",
                  kern_calc_visi_point_or_gauss, grid, threads,
                  d_components, d_component_beam_gains,
                  d_calc_visi_inouts->us, d_calc_visi_inouts->vs,
                  d_calc_visi_inouts->ws,
                  d_visibility_set->sum_visi_XX_real,
                  d_visibility_set->sum_visi_XX_imag,
                  d_visibility_set->sum_visi_XY_real,
                  d_visibility_set->sum_visi_XY_imag,
                  d_visibility_set->sum_visi_YX_real,
                  d_visibility_set->sum_visi_YX_imag,
                  d_visibility_set->sum_visi_YY_real,
                  d_visibility_set->sum_visi_YY_imag,
                  num_components, woden_settings->num_baselines,
                  woden_settings->num_freqs, woden_settings->num_cross,
                  woden_settings->num_time_steps, beamtype, comptype,
                  woden_settings->off_cardinal_dipoles);

}


__global__ void kern_calc_visi_shapelets(components_t d_components,
      beam_gains_t d_component_beam_gains,
      user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
      user_precision_t *d_allsteps_wavelengths,
      user_precision_t *d_u_shapes, user_precision_t *d_v_shapes,
      user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
      user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
      user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
      user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag,
      user_precision_t *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_cross,
      const int num_coeffs, int num_times, e_beamtype beamtype,
      int off_cardinal_dipoles) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iBaseline < num_cross) {
    int use_twobeams = d_component_beam_gains.use_twobeams;

    user_precision_t shape_flux_I;
    user_precision_t shape_flux_Q;
    user_precision_t shape_flux_U;
    user_precision_t shape_flux_V;
    gpuUserComplex visi_shape;

    int mod_baseline = iBaseline - num_baselines*floorf((float)iBaseline / (float)num_baselines);

    //Find out what time and freq index this baseline corresponds to
    int time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    int freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);

    for (int iCoeff = 0; iCoeff < num_coeffs; iCoeff++) {

      //We have multiple coefficients per SHAPELET component - reference
      //them via this array. We chunk over coeffs so might have any
      //number of components here
      int iComponent = d_components.param_indexes[iCoeff];
      int extrap_ind = num_freqs*iComponent + freq_ind;

      // if (iBaseline == 0) {
      //   printf("iComponent %d iCoeff %d extrap_ind %d\n", iComponent, iCoeff, extrap_ind);
      // }
      // printf("d_components.extrap_stokesI %p\n", d_components.extrap_stokesI);

      shape_flux_I = d_components.extrap_stokesI[extrap_ind];

      if (d_components.do_QUV == 1) {
        shape_flux_Q = d_components.extrap_stokesQ[extrap_ind];
        shape_flux_U = d_components.extrap_stokesU[extrap_ind];
        shape_flux_V = d_components.extrap_stokesV[extrap_ind];
      }

      visi_shape = calc_measurement_equation_gpu(d_us, d_vs, d_ws,
                            d_components.ls, d_components.ms, d_components.ns,
                            iBaseline, iComponent);

      user_precision_t pa = d_components.pas[iComponent];
      user_precision_t sinpa = sin(pa);
      user_precision_t cospa = cos(pa);

      int uv_stripe = num_baselines*num_times*iComponent + time_ind*num_baselines + mod_baseline;

      user_precision_t u_shape = d_u_shapes[uv_stripe] / d_allsteps_wavelengths[iBaseline];
      user_precision_t v_shape = d_v_shapes[uv_stripe] / d_allsteps_wavelengths[iBaseline];

      user_precision_t x = (cospa*v_shape + sinpa*u_shape); // major axis
      user_precision_t y = (-sinpa*v_shape + cospa*u_shape); // minor axis

      //Scales the FWHM to std to match basis functions, and account for the
      //basis functions being stored with beta = 1.0
      //Basis functions have been stored in such a way that x is in the same
      //direction as on sky, but y is opposite, so include negative here
      user_precision_t const_x = (d_components.majors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
      user_precision_t const_y = -(d_components.minors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

      // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
      gpuUserComplex Ipow_lookup[] = { make_gpuUserComplex(  1.0,  0.0 ),
                                       make_gpuUserComplex(  0.0,  1.0 ),
                                       make_gpuUserComplex( -1.0,  0.0 ),
                                       make_gpuUserComplex(  0.0, -1.0 ) };

      user_precision_t xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

      // find the indices in the basis functions for u*beta_u and v*beta_v

      user_precision_t xpos = x*const_x + sbf_c;
      user_precision_t ypos = y*const_y + sbf_c;

      int xindex = (int)floor(xpos);
      int yindex = (int)floor(ypos);
      //
      int n1 = (int)d_components.n1s[iCoeff];
      int n2 = (int)d_components.n2s[iCoeff];

      f_hat = d_components.shape_coeffs[iCoeff];

      sbf_n = &d_sbf[n1*sbf_L];
      xlow  = sbf_n[xindex];
      xhigh = sbf_n[xindex+1];
      u_value = xlow + (xhigh-xlow)*(xpos-xindex);

      sbf_n = &d_sbf[n2*sbf_L];
      ylow  = sbf_n[yindex];
      yhigh = sbf_n[yindex+1];
      v_value = ylow + (yhigh-ylow)*(ypos-yindex);

      // accumulate the intensity model for baseline pair (u,v)
      gpuUserComplex V_envelop = make_gpuUserComplex( 0.0, 0.0 );
      V_envelop = V_envelop + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

      visi_shape = visi_shape*V_envelop;

      if (d_components.do_QUV == 1) {
        update_sum_visis_stokesIQUV_gpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype, off_cardinal_dipoles,
             (gpuUserComplex *)d_component_beam_gains.gxs,
             (gpuUserComplex *)d_component_beam_gains.Dxs,
             (gpuUserComplex *)d_component_beam_gains.Dys,
             (gpuUserComplex *)d_component_beam_gains.gys,
             d_component_beam_gains.ant1_to_baseline_map,
             d_component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_shape,
             shape_flux_I, shape_flux_Q, shape_flux_U, shape_flux_V,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      } else {
        update_sum_visis_stokesI_gpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_shapes, num_times, beamtype, off_cardinal_dipoles,
             (gpuUserComplex *)d_component_beam_gains.gxs,
             (gpuUserComplex *)d_component_beam_gains.Dxs,
             (gpuUserComplex *)d_component_beam_gains.Dys,
             (gpuUserComplex *)d_component_beam_gains.gys,
             d_component_beam_gains.ant1_to_baseline_map,
             d_component_beam_gains.ant2_to_baseline_map, use_twobeams,
             visi_shape, shape_flux_I,
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);
      }
    }
  }
}

extern "C" void calc_visi_shapelets_gpu(components_t d_components,
                                        beam_gains_t d_component_beam_gains,
                                        calc_visi_inouts_t *d_calc_visi_inouts,
                                        visibility_set_t *d_visibility_set,
                                        int num_shapes, int num_shape_coeffs,
                                        e_beamtype beamtype,
                                        woden_settings_t *woden_settings){

  dim3 grid, threads;

  threads.x = 128;
  grid.x = (int)ceil( (float)woden_settings->num_cross / (float)threads.x );

  gpuErrorCheckKernel("kern_calc_visi_shapelets",
                  kern_calc_visi_shapelets, grid, threads,
                  d_components, d_component_beam_gains,
                  d_calc_visi_inouts->us, d_calc_visi_inouts->vs,
                  d_calc_visi_inouts->ws,
                  d_calc_visi_inouts->allsteps_wavelengths,
                  d_calc_visi_inouts->u_shapes, d_calc_visi_inouts->v_shapes,
                  d_visibility_set->sum_visi_XX_real,
                  d_visibility_set->sum_visi_XX_imag,
                  d_visibility_set->sum_visi_XY_real,
                  d_visibility_set->sum_visi_XY_imag,
                  d_visibility_set->sum_visi_YX_real,
                  d_visibility_set->sum_visi_YX_imag,
                  d_visibility_set->sum_visi_YY_real,
                  d_visibility_set->sum_visi_YY_imag,
                  d_calc_visi_inouts->sbf,  num_shapes,
                  woden_settings->num_baselines, woden_settings->num_freqs,
                  woden_settings->num_cross,
                  num_shape_coeffs, woden_settings->num_time_steps,
                  beamtype, woden_settings->off_cardinal_dipoles);
}



//Copy the sky model info from a set of components from the CPU to the GPU
void copy_components_to_GPU(source_t *chunked_source, source_t *d_chunked_source,
                            e_component_type comptype) {

  components_t *components=NULL;
  components_t *d_components=NULL;
  int num_comps = 0, num_shape_coeffs = 0;
  int num_powers = 0, num_curves = 0, num_lists = 0;

  if (comptype == POINT) {
    components = &chunked_source->point_components;
    d_components = &d_chunked_source->point_components;

    num_comps = chunked_source->n_points;
    num_shape_coeffs = 0;
    num_powers = chunked_source->n_point_powers;
    num_curves = chunked_source->n_point_curves;
    num_lists = chunked_source->n_point_lists;

  }
  else if (comptype == GAUSSIAN) {
    components = &chunked_source->gauss_components;
    d_components = &d_chunked_source->gauss_components;

    num_comps = chunked_source->n_gauss;
    num_shape_coeffs = 0;
    num_powers = chunked_source->n_gauss_powers;
    num_curves = chunked_source->n_gauss_curves;
    num_lists = chunked_source->n_gauss_lists;

  }
  // else if (comptype == SHAPELET) {
  else {
    components = &chunked_source->shape_components;
    d_components = &d_chunked_source->shape_components;

    num_comps = chunked_source->n_shapes;
    num_shape_coeffs = chunked_source->n_shape_coeffs;
    num_powers = chunked_source->n_shape_powers;
    num_curves = chunked_source->n_shape_curves;
    num_lists = chunked_source->n_shape_lists;

  }

  //Common attributes between all flux types and components types
  gpuMalloc( (void**)&d_components->ras, num_comps*sizeof(double) );
  gpuMemcpy( d_components->ras, components->ras, num_comps*sizeof(double),
                                                        gpuMemcpyHostToDevice );

  gpuMalloc( (void**)&d_components->decs, num_comps*sizeof(double) );
  gpuMemcpy( d_components->decs, components->decs, num_comps*sizeof(double),
                                                        gpuMemcpyHostToDevice );

  d_components->num_primarybeam_values = components->num_primarybeam_values;

  //GAUSSIAN and SHAPELET only attributes
  if (comptype == GAUSSIAN || comptype == SHAPELET ) {
    gpuMalloc( (void**)&d_components->pas, num_comps*sizeof(user_precision_t) );
    gpuMemcpy( d_components->pas, components->pas, num_comps*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->majors, num_comps*sizeof(user_precision_t) );
    gpuMemcpy( d_components->majors, components->majors, num_comps*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->minors, num_comps*sizeof(user_precision_t) );
    gpuMemcpy( d_components->minors, components->minors, num_comps*sizeof(user_precision_t),
                                                        gpuMemcpyHostToDevice );
  }

  //SHAPELET only attributes
  if (comptype == SHAPELET) {
    gpuMalloc( (void**)&d_components->shape_coeffs,
                        num_shape_coeffs*sizeof(user_precision_t) );
    gpuMemcpy( d_components->shape_coeffs, components->shape_coeffs,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->n1s,
                        num_shape_coeffs*sizeof(user_precision_t) );
    gpuMemcpy( d_components->n1s, components->n1s,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->n2s,
                        num_shape_coeffs*sizeof(user_precision_t) );
    gpuMemcpy( d_components->n2s, components->n2s,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->param_indexes,
                        num_shape_coeffs*sizeof(user_precision_t) );
    gpuMemcpy( d_components->param_indexes, components->param_indexes,
                        num_shape_coeffs*sizeof(user_precision_t),
                        gpuMemcpyHostToDevice );
  }

  //POWER_LAW flux things
  if (num_powers > 0) {
    gpuMalloc( (void**)&d_components->power_comp_inds,
                        num_powers*sizeof(int) );
    gpuMemcpy( d_components->power_comp_inds, components->power_comp_inds,
                        num_powers*sizeof(int), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->power_ref_stokesI,
                        num_powers*sizeof(user_precision_t) );
    gpuMemcpy( d_components->power_ref_stokesI, components->power_ref_stokesI,
                        num_powers*sizeof(user_precision_t), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->power_SIs,
                        num_powers*sizeof(user_precision_t) );
    gpuMemcpy( d_components->power_SIs, components->power_SIs,
                        num_powers*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  }

  //CURVED_POWER_LAW things
  if (num_curves > 0) {
    gpuMalloc( (void**)&d_components->curve_comp_inds,
                        num_curves*sizeof(int) );
    gpuMemcpy( d_components->curve_comp_inds, components->curve_comp_inds,
                        num_curves*sizeof(int), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->curve_ref_stokesI,
                        num_curves*sizeof(user_precision_t) );
    gpuMemcpy( d_components->curve_ref_stokesI, components->curve_ref_stokesI,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->curve_SIs,
                        num_curves*sizeof(user_precision_t) );
    gpuMemcpy( d_components->curve_SIs, components->curve_SIs,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->curve_qs,
                        num_curves*sizeof(user_precision_t) );
    gpuMemcpy( d_components->curve_qs, components->curve_qs,
                        num_curves*sizeof(user_precision_t), gpuMemcpyHostToDevice );
  }

  //LIST things
  if (num_lists > 0) {
    int num_list_values = components->total_num_flux_entires;

    gpuMalloc( (void**)&d_components->list_comp_inds,
                        num_lists*sizeof(int) );
    gpuMemcpy( d_components->list_comp_inds,
                        components->list_comp_inds,
                        num_lists*sizeof(int), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->num_list_values,
                        num_lists*sizeof(int) );
    gpuMemcpy( d_components->num_list_values,
                        components->num_list_values,
                        num_lists*sizeof(int), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->list_start_indexes,
                        num_lists*sizeof(int) );
    gpuMemcpy( d_components->list_start_indexes,
                        components->list_start_indexes,
                        num_lists*sizeof(int), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->list_freqs,
                        num_list_values*sizeof(double) );
    gpuMemcpy( d_components->list_freqs, components->list_freqs,
                        num_list_values*sizeof(double), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->list_stokesI,
                        num_list_values*sizeof(user_precision_t) );
    gpuMemcpy( d_components->list_stokesI, components->list_stokesI,
                        num_list_values*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  }

  int n_stokesV_pol_frac = components->n_stokesV_pol_frac;
  int n_stokesV_power = components->n_stokesV_power;
  int n_stokesV_curve = components->n_stokesV_curve;
  int n_stokesV_list = components->n_stokesV_list;
  int n_stokesV_list_flux_entries = components->n_stokesV_list_flux_entries;
  int n_linpol_pol_frac = components->n_linpol_pol_frac;
  int n_linpol_power = components->n_linpol_power;
  int n_linpol_curve = components->n_linpol_curve;
  int n_linpol_angles = components->n_linpol_angles;
  int n_linpol_list = components->n_linpol_list;
  int n_stokesQ_list_flux_entries = components->n_stokesQ_list_flux_entries;
  int n_stokesU_list_flux_entries = components->n_stokesU_list_flux_entries;
  int n_linpol_p_list = components->n_linpol_p_list;
  int n_linpol_p_list_flux_entries = components->n_linpol_p_list_flux_entries;

  d_components->n_stokesV_pol_frac = n_stokesV_pol_frac;
  d_components->n_stokesV_power = n_stokesV_power;
  d_components->n_stokesV_curve = n_stokesV_curve;
  d_components->n_stokesV_list = n_stokesV_list;
  d_components->n_stokesV_list_flux_entries = n_stokesV_list_flux_entries;
  d_components->n_linpol_pol_frac = n_linpol_pol_frac;
  d_components->n_linpol_power = n_linpol_power;
  d_components->n_linpol_curve = n_linpol_curve;
  d_components->n_linpol_list = n_linpol_list;
  d_components->n_stokesQ_list_flux_entries = n_stokesQ_list_flux_entries;
  d_components->n_stokesU_list_flux_entries = n_stokesU_list_flux_entries;
  d_components->n_linpol_p_list = n_linpol_p_list;
  d_components->n_linpol_p_list_flux_entries = n_linpol_p_list_flux_entries;

  d_components->n_linpol_angles = n_linpol_angles;
  d_components->do_QUV = components->do_QUV;

  if (n_stokesV_pol_frac > 0) {
    gpuMalloc( (void**)&d_components->stokesV_pol_fracs,
                        n_stokesV_pol_frac*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesV_pol_fracs, components->stokesV_pol_fracs,
                n_stokesV_pol_frac*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_pol_frac_comp_inds,
                        n_stokesV_pol_frac*sizeof(int) );
    gpuMemcpy( d_components->stokesV_pol_frac_comp_inds, components->stokesV_pol_frac_comp_inds,
                n_stokesV_pol_frac*sizeof(int), gpuMemcpyHostToDevice );
  }
  if (n_stokesV_power > 0){
    gpuMalloc( (void**)&d_components->stokesV_power_ref_flux,
                        n_stokesV_power*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesV_power_ref_flux, components->stokesV_power_ref_flux,
                n_stokesV_power*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_power_SIs,
                        n_stokesV_power*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesV_power_SIs, components->stokesV_power_SIs,
                n_stokesV_power*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_power_comp_inds,
                        n_stokesV_power*sizeof(int) );
    gpuMemcpy( d_components->stokesV_power_comp_inds, components->stokesV_power_comp_inds,
                n_stokesV_power*sizeof(int), gpuMemcpyHostToDevice );
  }

  if (n_stokesV_curve > 0){
    gpuMalloc( (void**)&d_components->stokesV_curve_ref_flux,
                        n_stokesV_curve*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesV_curve_ref_flux, components->stokesV_curve_ref_flux,
                n_stokesV_curve*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_curve_SIs,
                        n_stokesV_curve*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesV_curve_SIs, components->stokesV_curve_SIs,
                n_stokesV_curve*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_curve_qs,
                        n_stokesV_curve*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesV_curve_qs, components->stokesV_curve_qs,
                n_stokesV_curve*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_curve_comp_inds,
                        n_stokesV_curve*sizeof(int) );
    gpuMemcpy( d_components->stokesV_curve_comp_inds, components->stokesV_curve_comp_inds,
                n_stokesV_curve*sizeof(int), gpuMemcpyHostToDevice );
  }

  if (n_stokesV_list > 0) {
    gpuMalloc( (void**)&d_components->stokesV_num_list_values,
                        n_stokesV_list*sizeof(int) );
    gpuMemcpy( d_components->stokesV_num_list_values, components->stokesV_num_list_values,
                n_stokesV_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_list_start_indexes,
                        n_stokesV_list*sizeof(int) );
    gpuMemcpy( d_components->stokesV_list_start_indexes, components->stokesV_list_start_indexes,
                n_stokesV_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_list_comp_inds,
                        n_stokesV_list*sizeof(int) );
    gpuMemcpy( d_components->stokesV_list_comp_inds, components->stokesV_list_comp_inds,
                n_stokesV_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_list_ref_freqs,
                        n_stokesV_list_flux_entries*sizeof(double) );
    gpuMemcpy( d_components->stokesV_list_ref_freqs, components->stokesV_list_ref_freqs,
                n_stokesV_list_flux_entries*sizeof(double), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesV_list_ref_flux,
                        n_stokesV_list_flux_entries*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesV_list_ref_flux, components->stokesV_list_ref_flux,
                n_stokesV_list_flux_entries*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  }

  if (n_linpol_pol_frac > 0) {
    gpuMalloc( (void**)&d_components->linpol_pol_fracs,
                        n_linpol_pol_frac*sizeof(user_precision_t) );
    gpuMemcpy( d_components->linpol_pol_fracs, components->linpol_pol_fracs,
                n_linpol_pol_frac*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_pol_frac_comp_inds,
                        n_linpol_pol_frac*sizeof(int) );
    gpuMemcpy( d_components->linpol_pol_frac_comp_inds, components->linpol_pol_frac_comp_inds,
                n_linpol_pol_frac*sizeof(int), gpuMemcpyHostToDevice );
  }
  if (n_linpol_power > 0){
    gpuMalloc( (void**)&d_components->linpol_power_ref_flux,
                        n_linpol_power*sizeof(user_precision_t) );
    gpuMemcpy( d_components->linpol_power_ref_flux, components->linpol_power_ref_flux,
                n_linpol_power*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_power_SIs,
                        n_linpol_power*sizeof(user_precision_t) );
    gpuMemcpy( d_components->linpol_power_SIs, components->linpol_power_SIs,
                n_linpol_power*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_power_comp_inds,
                        n_linpol_power*sizeof(int) );
    gpuMemcpy( d_components->linpol_power_comp_inds, components->linpol_power_comp_inds,
                n_linpol_power*sizeof(int), gpuMemcpyHostToDevice );
  }

  if (n_linpol_curve > 0){
    gpuMalloc( (void**)&d_components->linpol_curve_ref_flux,
                        n_linpol_curve*sizeof(user_precision_t) );
    gpuMemcpy( d_components->linpol_curve_ref_flux, components->linpol_curve_ref_flux,
                n_linpol_curve*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_curve_SIs,
                        n_linpol_curve*sizeof(user_precision_t) );
    gpuMemcpy( d_components->linpol_curve_SIs, components->linpol_curve_SIs,
                n_linpol_curve*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_curve_qs,
                        n_linpol_curve*sizeof(user_precision_t) );
    gpuMemcpy( d_components->linpol_curve_qs, components->linpol_curve_qs,
                n_linpol_curve*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_curve_comp_inds,
                        n_linpol_curve*sizeof(int) );
    gpuMemcpy( d_components->linpol_curve_comp_inds, components->linpol_curve_comp_inds,
                n_linpol_curve*sizeof(int), gpuMemcpyHostToDevice );
  }

  if (n_linpol_list > 0) {
    gpuMalloc( (void**)&d_components->stokesQ_num_list_values,
                        n_linpol_list*sizeof(int) );
    gpuMemcpy( d_components->stokesQ_num_list_values, components->stokesQ_num_list_values,
                n_linpol_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesQ_list_start_indexes,
                        n_linpol_list*sizeof(int) );
    gpuMemcpy( d_components->stokesQ_list_start_indexes, components->stokesQ_list_start_indexes,
                n_linpol_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesQ_list_comp_inds,
                        n_linpol_list*sizeof(int) );
    gpuMemcpy( d_components->stokesQ_list_comp_inds, components->stokesQ_list_comp_inds,
                n_linpol_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesQ_list_ref_freqs,
                        n_stokesQ_list_flux_entries*sizeof(double) );
    gpuMemcpy( d_components->stokesQ_list_ref_freqs, components->stokesQ_list_ref_freqs,
                n_stokesQ_list_flux_entries*sizeof(double), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesQ_list_ref_flux,
                        n_stokesQ_list_flux_entries*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesQ_list_ref_flux, components->stokesQ_list_ref_flux,
                n_stokesQ_list_flux_entries*sizeof(user_precision_t), gpuMemcpyHostToDevice );

    gpuMalloc( (void**)&d_components->stokesU_num_list_values,
                        n_linpol_list*sizeof(int) );
    gpuMemcpy( d_components->stokesU_num_list_values, components->stokesU_num_list_values,
                n_linpol_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesU_list_start_indexes,
                        n_linpol_list*sizeof(int) );
    gpuMemcpy( d_components->stokesU_list_start_indexes, components->stokesU_list_start_indexes,
                n_linpol_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesU_list_comp_inds,
                        n_linpol_list*sizeof(int) );
    gpuMemcpy( d_components->stokesU_list_comp_inds, components->stokesU_list_comp_inds,
                n_linpol_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesU_list_ref_freqs,
                        n_stokesU_list_flux_entries*sizeof(double) );
    gpuMemcpy( d_components->stokesU_list_ref_freqs, components->stokesU_list_ref_freqs,
                n_stokesU_list_flux_entries*sizeof(double), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->stokesU_list_ref_flux,
                        n_stokesU_list_flux_entries*sizeof(user_precision_t) );
    gpuMemcpy( d_components->stokesU_list_ref_flux, components->stokesU_list_ref_flux,
                n_stokesU_list_flux_entries*sizeof(user_precision_t), gpuMemcpyHostToDevice );

  }

  if (n_linpol_p_list > 0) {
    gpuMalloc( (void**)&d_components->linpol_p_num_list_values,
                        n_linpol_p_list*sizeof(int) );
    gpuMemcpy( d_components->linpol_p_num_list_values, components->linpol_p_num_list_values,
                n_linpol_p_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_p_list_start_indexes,
                        n_linpol_p_list*sizeof(int) );
    gpuMemcpy( d_components->linpol_p_list_start_indexes, components->linpol_p_list_start_indexes,
                n_linpol_p_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_p_list_comp_inds,
                        n_linpol_p_list*sizeof(int) );
    gpuMemcpy( d_components->linpol_p_list_comp_inds, components->linpol_p_list_comp_inds,
                n_linpol_p_list*sizeof(int), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_p_list_ref_freqs,
                        n_linpol_p_list_flux_entries*sizeof(double) );
    gpuMemcpy( d_components->linpol_p_list_ref_freqs, components->linpol_p_list_ref_freqs,
                n_linpol_p_list_flux_entries*sizeof(double), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_p_list_ref_flux,
                        n_linpol_p_list_flux_entries*sizeof(user_precision_t) );
    gpuMemcpy( d_components->linpol_p_list_ref_flux, components->linpol_p_list_ref_flux,
                n_linpol_p_list_flux_entries*sizeof(user_precision_t), gpuMemcpyHostToDevice );
                
  }


  if (n_linpol_angles > 0){
    gpuMalloc( (void**)&d_components->rm_values,
                        n_linpol_angles*sizeof(user_precision_t) );
    gpuMemcpy( d_components->rm_values, components->rm_values,
                n_linpol_angles*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->intr_pol_angle,
                        n_linpol_angles*sizeof(user_precision_t) );
    gpuMemcpy( d_components->intr_pol_angle, components->intr_pol_angle,
                n_linpol_angles*sizeof(user_precision_t), gpuMemcpyHostToDevice );
    gpuMalloc( (void**)&d_components->linpol_angle_inds,
                        n_linpol_angles*sizeof(int) );
    gpuMemcpy( d_components->linpol_angle_inds, components->linpol_angle_inds,
                n_linpol_angles*sizeof(int), gpuMemcpyHostToDevice );
  }
}

extern "C" source_t * copy_chunked_source_to_GPU(source_t *chunked_source){

  source_t *d_chunked_source = (source_t*)malloc(sizeof(source_t));

  if (chunked_source->n_points > 0) {
    copy_components_to_GPU(chunked_source, d_chunked_source, POINT);
  }
  if (chunked_source->n_gauss > 0) {
    copy_components_to_GPU(chunked_source, d_chunked_source, GAUSSIAN);
  }
  if (chunked_source->n_shapes > 0) {
    copy_components_to_GPU(chunked_source, d_chunked_source, SHAPELET);
  }

  //copy across the component counters

  d_chunked_source->n_points = chunked_source->n_points;
  d_chunked_source->n_point_lists = chunked_source->n_point_lists;
  d_chunked_source->n_point_powers = chunked_source->n_point_powers;
  d_chunked_source->n_point_curves = chunked_source->n_point_curves;

  d_chunked_source->n_gauss = chunked_source->n_gauss;
  d_chunked_source->n_gauss_lists = chunked_source->n_gauss_lists;
  d_chunked_source->n_gauss_powers = chunked_source->n_gauss_powers;
  d_chunked_source->n_gauss_curves = chunked_source->n_gauss_curves;

  d_chunked_source->n_shapes = chunked_source->n_shapes;
  d_chunked_source->n_shape_lists = chunked_source->n_shape_lists;
  d_chunked_source->n_shape_powers = chunked_source->n_shape_powers;
  d_chunked_source->n_shape_curves = chunked_source->n_shape_curves;
  d_chunked_source->n_shape_coeffs = chunked_source->n_shape_coeffs;

  return d_chunked_source;
}

//TODO might be circumstances in the future where there are no leakages
//allocated in `components`, so need a switch to not copy them
extern "C" void copy_CPU_component_gains_to_GPU_beam_gains(components_t *components,
  beam_gains_t *d_beam_gains, int num_gains) {


    gpuMalloc( (void**)&d_beam_gains->gxs,
                      num_gains*sizeof(user_precision_complex_t) );
    gpuMalloc( (void**)&d_beam_gains->Dxs,
                      num_gains*sizeof(user_precision_complex_t) );
    gpuMalloc( (void**)&d_beam_gains->Dys,
                      num_gains*sizeof(user_precision_complex_t) );
    gpuMalloc( (void**)&d_beam_gains->gys,
                      num_gains*sizeof(user_precision_complex_t) );

    // printf("Copying beam gains to GPU %p %p\n", components->gxs, d_beam_gains->gxs);

    gpuMemcpy( d_beam_gains->gxs, components->gxs,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

    gpuMemcpy( d_beam_gains->Dxs, components->Dxs,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

    gpuMemcpy( d_beam_gains->Dys, components->Dys,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

    gpuMemcpy( d_beam_gains->gys, components->gys,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
}


extern "C" void copy_CPU_beam_gains_to_GPU_beam_gains(beam_gains_t *beam_gains,
  beam_gains_t *d_beam_gains, int num_gains) {


    gpuMalloc( (void**)&d_beam_gains->gxs,
                      num_gains*sizeof(user_precision_complex_t) );
    gpuMalloc( (void**)&d_beam_gains->Dxs,
                      num_gains*sizeof(user_precision_complex_t) );
    gpuMalloc( (void**)&d_beam_gains->Dys,
                      num_gains*sizeof(user_precision_complex_t) );
    gpuMalloc( (void**)&d_beam_gains->gys,
                      num_gains*sizeof(user_precision_complex_t) );

    // printf("Copying beam gains to GPU %p %p\n", components->gxs, d_beam_gains->gxs);

    gpuMemcpy( d_beam_gains->gxs, beam_gains->gxs,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

    gpuMemcpy( d_beam_gains->Dxs, beam_gains->Dxs,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

    gpuMemcpy( d_beam_gains->Dys, beam_gains->Dys,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );

    gpuMemcpy( d_beam_gains->gys, beam_gains->gys,
               num_gains*sizeof(user_precision_complex_t), gpuMemcpyHostToDevice );
}

extern "C" void free_extrapolated_flux_arrays_gpu(components_t *d_components){
  gpuFree( d_components->extrap_stokesI );

  if (d_components->do_QUV) {
    gpuFree( d_components->extrap_stokesQ );
    gpuFree( d_components->extrap_stokesU );
    gpuFree( d_components->extrap_stokesV );
  }
}



extern "C" void free_components_gpu(source_t *d_chunked_source,
                                  e_component_type comptype){
  components_t d_components;
  int n_powers = 0;
  int n_curves = 0;
  int n_lists = 0;

  if (comptype == POINT) {
    d_components = d_chunked_source->point_components;
    n_powers = d_chunked_source->n_point_powers;
    n_curves = d_chunked_source->n_point_curves;
    n_lists = d_chunked_source->n_point_lists;
  }
  else if (comptype == GAUSSIAN) {
    d_components = d_chunked_source->gauss_components;
    n_powers = d_chunked_source->n_gauss_powers;
    n_curves = d_chunked_source->n_gauss_curves;
    n_lists = d_chunked_source->n_gauss_lists;
  }
  else {
    d_components = d_chunked_source->shape_components;
    n_powers = d_chunked_source->n_shape_powers;
    n_curves = d_chunked_source->n_shape_curves;
    n_lists = d_chunked_source->n_shape_lists;
  }

  gpuFree( d_components.decs);
  gpuFree( d_components.ras);
  gpuFree( d_components.ls);
  gpuFree( d_components.ms);
  gpuFree( d_components.ns);

  //The az,za,beam_has,beam_decs are handled by other functions

  if (n_powers > 0) {
    gpuFree( d_components.power_ref_stokesI );
    gpuFree( d_components.power_SIs );
    gpuFree( d_components.power_comp_inds );
  }

  if (n_curves > 0) {
    gpuFree( d_components.curve_ref_stokesI );
    gpuFree( d_components.curve_SIs );
    gpuFree( d_components.curve_qs );
    gpuFree( d_components.curve_comp_inds );
  }
  if (n_lists > 0) {
    gpuFree( d_components.list_comp_inds );
    gpuFree( d_components.list_freqs );
    gpuFree( d_components.list_stokesI );
    gpuFree( d_components.num_list_values );
    gpuFree( d_components.list_start_indexes );
  }

  if (comptype == GAUSSIAN || comptype == SHAPELET) {
    gpuFree( d_components.pas );
    gpuFree( d_components.majors );
    gpuFree( d_components.minors );
  }

  if (comptype == SHAPELET) {
    gpuFree( d_components.shape_coeffs );
    gpuFree( d_components.n1s );
    gpuFree( d_components.n2s );
    gpuFree( d_components.param_indexes );
  }
  //Free whatever polarisation information we have
  if (d_components.n_stokesV_pol_frac > 0) {
    gpuFree(d_components.stokesV_pol_fracs);
    gpuFree(d_components.stokesV_pol_frac_comp_inds);
  }
  if (d_components.n_stokesV_power > 0){
    gpuFree(d_components.stokesV_power_ref_flux);
    gpuFree(d_components.stokesV_power_SIs);
    gpuFree(d_components.stokesV_power_comp_inds);
  }

  if (d_components.n_stokesV_curve > 0){
    gpuFree(d_components.stokesV_curve_ref_flux);
    gpuFree(d_components.stokesV_curve_SIs);
    gpuFree(d_components.stokesV_curve_qs);
    gpuFree(d_components.stokesV_curve_comp_inds);
  }

  if (d_components.n_stokesV_list > 0) {
    gpuFree(d_components.stokesV_num_list_values);
    gpuFree(d_components.stokesV_list_start_indexes);
    gpuFree(d_components.stokesV_list_comp_inds);
    gpuFree(d_components.stokesV_list_ref_freqs);
    gpuFree(d_components.stokesV_list_ref_flux);
  }

  if (d_components.n_linpol_pol_frac > 0) {
    gpuFree(d_components.linpol_pol_fracs);
    gpuFree(d_components.linpol_pol_frac_comp_inds);
  }
  if (d_components.n_linpol_power > 0){
    gpuFree(d_components.linpol_power_ref_flux);
    gpuFree(d_components.linpol_power_SIs);
    gpuFree(d_components.linpol_power_comp_inds);
  }

  if (d_components.n_linpol_curve > 0){
    gpuFree(d_components.linpol_curve_ref_flux);
    gpuFree(d_components.linpol_curve_SIs);
    gpuFree(d_components.linpol_curve_qs);
    gpuFree(d_components.linpol_curve_comp_inds);
  }

  if (d_components.n_linpol_list > 0) {
    gpuFree(d_components.stokesQ_num_list_values);
    gpuFree(d_components.stokesQ_list_start_indexes);
    gpuFree(d_components.stokesQ_list_comp_inds);
    gpuFree(d_components.stokesQ_list_ref_freqs);
    gpuFree(d_components.stokesQ_list_ref_flux);
    gpuFree(d_components.stokesU_num_list_values);
    gpuFree(d_components.stokesU_list_start_indexes);
    gpuFree(d_components.stokesU_list_comp_inds);
    gpuFree(d_components.stokesU_list_ref_freqs);
    gpuFree(d_components.stokesU_list_ref_flux);
  }

  if (d_components.n_linpol_p_list > 0) {
    gpuFree(d_components.linpol_p_num_list_values);
    gpuFree(d_components.linpol_p_list_start_indexes);
    gpuFree(d_components.linpol_p_list_comp_inds);
    gpuFree(d_components.linpol_p_list_ref_freqs);
    gpuFree(d_components.linpol_p_list_ref_flux);
  }

  if (d_components.n_linpol_angles > 0){
    gpuFree(d_components.rm_values);
    gpuFree(d_components.intr_pol_angle);
    gpuFree(d_components.linpol_angle_inds);
  }
}

extern "C" void free_beam_gains_gpu(beam_gains_t *d_beam_gains, e_beamtype beamtype){

  gpuFree( d_beam_gains->gxs );
  gpuFree( d_beam_gains->gys );
  // gpuFree( d_beam_gains->Dxs );
  // gpuFree( d_beam_gains->Dys );

  if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY || beamtype == EB_OSKAR || beamtype == EB_LOFAR || beamtype == EB_MWA){
    gpuFree( d_beam_gains->Dxs );
    gpuFree( d_beam_gains->Dys );
  }
}

//Calculate auto-correlations
__global__ void kern_calc_autos(components_t d_components,
                                beam_gains_t d_component_beam_gains,
                                int beamtype,
                                int num_components, int num_baselines,
                                int num_freqs, int num_times, int num_ants,
                                user_precision_t *d_sum_visi_XX_real,
                                user_precision_t *d_sum_visi_XX_imag,
                                user_precision_t *d_sum_visi_XY_real,
                                user_precision_t *d_sum_visi_XY_imag,
                                user_precision_t *d_sum_visi_YX_real,
                                user_precision_t *d_sum_visi_YX_imag,
                                user_precision_t *d_sum_visi_YY_real,
                                user_precision_t *d_sum_visi_YY_imag,
                                int use_twobeams,
                                int *d_ant1_to_auto_map,
                                int *d_ant2_to_auto_map,
                                int off_cardinal_dipoles) {

  // Start by computing which baseline we're going to do
  const int iTimeFreq = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iAnt = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iAnt < num_ants && iTimeFreq < num_times*num_freqs) {

    int time_ind = (int)floorf( (float)iTimeFreq / (float)num_freqs);
    int freq_ind = iTimeFreq - time_ind*num_freqs;

    //Set up iBaseline to be a cross-pol of the correct time
    //and frequency step, that also correpsonds to the correct antenna
    //get_beam_gains_gpu and get_beam_gains_multibeams_gpu will use this to access the
    //correct beam gains.
    int iBaseline = num_baselines*num_freqs*time_ind + num_baselines*freq_ind + iAnt;

    int num_visis = num_baselines*num_freqs*num_times;
    int iAuto = num_visis + num_ants*num_freqs*time_ind + num_ants*freq_ind + iAnt;

    gpuUserComplex auto_XX, auto_XY, auto_YX, auto_YY;
    gpuUserComplex g1x, D1x, D1y, g1y, g2x, D2x, D2y, g2y;

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      if (use_twobeams == 1){
        get_beam_gains_multibeams_gpu(iBaseline, iComponent, num_freqs,
                num_baselines, num_components, num_times, beamtype,
                (gpuUserComplex *)d_component_beam_gains.gxs,
                (gpuUserComplex *)d_component_beam_gains.Dxs,
                (gpuUserComplex *)d_component_beam_gains.Dys,
                (gpuUserComplex *)d_component_beam_gains.gys,
                d_ant1_to_auto_map, d_ant2_to_auto_map,
                &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      }
      else {
        // printf("WE IS ONLY DOING THIS TING\n");
        get_beam_gains_gpu(iBaseline, iComponent, num_freqs,
                num_baselines, num_components, num_times, beamtype,
                (gpuUserComplex *)d_component_beam_gains.gxs,
                (gpuUserComplex *)d_component_beam_gains.Dxs,
                (gpuUserComplex *)d_component_beam_gains.Dys,
                (gpuUserComplex *)d_component_beam_gains.gys,
                &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      }

      gpuUserComplex visi_component;
      visi_component = make_gpuUserComplex(1.0, 0.0);

      int extrap_ind = num_freqs*iComponent + freq_ind;

      user_precision_t flux_I = d_components.extrap_stokesI[extrap_ind];

      if (d_components.do_QUV == 1) {
        user_precision_t flux_Q = d_components.extrap_stokesQ[extrap_ind];
        user_precision_t flux_U = d_components.extrap_stokesU[extrap_ind];
        user_precision_t flux_V = d_components.extrap_stokesV[extrap_ind];

        if (off_cardinal_dipoles == 1) {
          apply_beam_gains_stokesIQUV_off_cardinal_gpu(g1x, D1x, D1y, g1y,
                                    g2x, D2x, D2y, g2y,
                                    flux_I, flux_Q, flux_U, flux_V,
                                    visi_component,
                                    &auto_XX, &auto_XY, &auto_YX, &auto_YY);
        } else {
          apply_beam_gains_stokesIQUV_on_cardinal_gpu(g1x, D1x, D1y, g1y,
                                    g2x, D2x, D2y, g2y,
                                    flux_I, flux_Q, flux_U, flux_V,
                                    visi_component,
                                    &auto_XX, &auto_XY, &auto_YX, &auto_YY);
        }
      } else {
        if (off_cardinal_dipoles == 1) {
          apply_beam_gains_stokesI_off_cardinal_gpu(g1x, D1x, D1y, g1y,
                                 g2x, D2x, D2y, g2y,
                                 flux_I, visi_component,
                                 &auto_XX, &auto_XY, &auto_YX, &auto_YY);
        } else {
          apply_beam_gains_stokesI_on_cardinal_gpu(g1x, D1x, D1y, g1y,
                                 g2x, D2x, D2y, g2y,
                                 flux_I, visi_component,
                                 &auto_XX, &auto_XY, &auto_YX, &auto_YY);
        }
      }
      
      d_sum_visi_XX_real[iAuto] += auto_XX.x;
      d_sum_visi_XX_imag[iAuto] += auto_XX.y;

      d_sum_visi_XY_real[iAuto] += auto_XY.x;
      d_sum_visi_XY_imag[iAuto] += auto_XY.y;

      d_sum_visi_YX_real[iAuto] += auto_YX.x;
      d_sum_visi_YX_imag[iAuto] += auto_YX.y;

      d_sum_visi_YY_real[iAuto] += auto_YY.x;
      d_sum_visi_YY_imag[iAuto] += auto_YY.y;

    }
  }
}

extern "C" void fill_ant_to_baseline_mapping_gpu(int num_ants, int *d_ant1_to_baseline_map,
                                                 int *d_ant2_to_baseline_map){

  int num_baselines = ((num_ants - 1)*num_ants) / 2;

  // gpuMalloc( (void**)&d_ant1_to_baseline_map, num_baselines*sizeof(int) );
  // gpuMalloc( (void**)&d_ant2_to_baseline_map, num_baselines*sizeof(int) );

  int *ant1_to_baseline_map = NULL;
  int *ant2_to_baseline_map = NULL;

  ant1_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));
  ant2_to_baseline_map = (int *)malloc(num_baselines*sizeof(int));

  //These functions only do cross correlations, so create all combos of antennas
  //that make up all the crosses
  int cross_index = 0;
  for (int ant1 = 0; ant1 < num_ants-1; ant1++)
  {
    for (int ant2 = ant1 + 1; ant2 < num_ants; ant2++)
    {
      ant1_to_baseline_map[cross_index] = ant1;
      ant2_to_baseline_map[cross_index] = ant2;

      cross_index += 1;
    }
  }

  gpuMemcpy(d_ant1_to_baseline_map, ant1_to_baseline_map,
                                  num_baselines*sizeof(int), gpuMemcpyHostToDevice );
  gpuMemcpy(d_ant2_to_baseline_map, ant2_to_baseline_map,
                                  num_baselines*sizeof(int), gpuMemcpyHostToDevice );

  free(ant1_to_baseline_map);
  free(ant2_to_baseline_map);

}

extern "C" void calc_autos_gpu(components_t *d_components,
                               beam_settings_t *beam_settings,
                               beam_gains_t *d_component_beam_gains,
                               visibility_set_t *d_visibility_set,
                               woden_settings_t *woden_settings,
                               int num_components, int use_twobeams){
  int num_freqs = woden_settings->num_freqs;
  int num_times = woden_settings->num_time_steps;
  int num_ants = woden_settings->num_ants;
  int num_baselines = woden_settings->num_baselines;

  dim3 threads, grid;

  threads.x = 64;
  threads.y = 2;
  threads.z = 1;
  grid.x = (int)ceil( (float)(num_freqs*num_times) / (float)threads.x );
  grid.y = (int)ceil( (float)(num_ants) / (float)threads.y );
  grid.z = 1;

  int *d_ant_to_auto_map = NULL;

  if (use_twobeams == 1) {
    int *ant_to_auto_map = NULL;
    ant_to_auto_map = (int *)malloc(num_ants*sizeof(int));
    for (int ant = 0; ant < num_ants; ant++){
        ant_to_auto_map[ant] = ant;
    }
    gpuMalloc( (void**)&d_ant_to_auto_map, num_ants*sizeof(int) );
    gpuMemcpy(d_ant_to_auto_map, ant_to_auto_map, num_ants*sizeof(int),
                                                        gpuMemcpyHostToDevice );
    free(ant_to_auto_map);
  }

  gpuErrorCheckKernel("kern_calc_autos",
                kern_calc_autos, grid, threads,
                *d_components, *d_component_beam_gains,
                beam_settings->beamtype,
                num_components, num_baselines,
                num_freqs, num_times, num_ants,
                d_visibility_set->sum_visi_XX_real,
                d_visibility_set->sum_visi_XX_imag,
                d_visibility_set->sum_visi_XY_real,
                d_visibility_set->sum_visi_XY_imag,
                d_visibility_set->sum_visi_YX_real,
                d_visibility_set->sum_visi_YX_imag,
                d_visibility_set->sum_visi_YY_real,
                d_visibility_set->sum_visi_YY_imag,
                use_twobeams, d_ant_to_auto_map,
                d_ant_to_auto_map,
                woden_settings->off_cardinal_dipoles);

  if (use_twobeams == 1) {
     gpuFree( d_ant_to_auto_map );
  }
}

//------------------------------------------------------------------------------
//Functions below to be used in unit tests. Unless you want to do separate linking
//and much faffing, it's hard to import any __device__ functions to separate
//test files. So below, we write kernels (__global__) around anything that we want
//to test so we can call them from the test files.
//------------------------------------------------------------------------------

__global__ void kern_calc_measurement_equation(int num_components, int num_baselines,
          user_precision_t *d_us, user_precision_t *d_vs, user_precision_t *d_ws,
          double *d_ls, double *d_ms, double *d_ns, gpuUserComplex *d_visis) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iComponent < num_components && iBaseline < num_baselines) {

    gpuUserComplex visi;
    visi = calc_measurement_equation_gpu(d_us, d_vs, d_ws, d_ls, d_ms, d_ns,
                                     iBaseline, iComponent);

    int visi_ind = num_components*iBaseline + iComponent;
    d_visis[visi_ind] = visi;

  }
}

__global__ void kern_apply_beam_gains(int num_gains,
          gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
          gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
          gpuUserComplex *d_g2xs, gpuUserComplex *d_D2xs,
          gpuUserComplex *d_D2ys, gpuUserComplex *d_g2ys,
          user_precision_t *d_flux_Is, user_precision_t *d_flux_Qs,
          user_precision_t *d_flux_Us, user_precision_t *d_flux_Vs,
          gpuUserComplex *d_visi_components,
          gpuUserComplex *d_visi_XXs, gpuUserComplex *d_visi_XYs,
          gpuUserComplex *d_visi_YXs, gpuUserComplex *d_visi_YYs, 
          int off_cardinal_dipoles, int do_QUV) {

  const int iGain = threadIdx.x + (blockDim.x*blockIdx.x);
  if (iGain < num_gains) {

    gpuUserComplex visi_XX;
    gpuUserComplex visi_XY;
    gpuUserComplex visi_YX;
    gpuUserComplex visi_YY;

    if (do_QUV == 1) {
      if (off_cardinal_dipoles == 1){
        apply_beam_gains_stokesIQUV_off_cardinal_gpu(d_g1xs[iGain], d_D1xs[iGain],
              d_D1ys[iGain], d_g1ys[iGain],
              d_g2xs[iGain], d_D2xs[iGain],
              d_D2ys[iGain], d_g2ys[iGain],
              d_flux_Is[iGain], d_flux_Qs[iGain],
              d_flux_Us[iGain], d_flux_Vs[iGain],
              d_visi_components[iGain],
              &visi_XX, &visi_XY,
              &visi_YX, &visi_YY);
      } else {
        apply_beam_gains_stokesIQUV_on_cardinal_gpu(d_g1xs[iGain], d_D1xs[iGain],
              d_D1ys[iGain], d_g1ys[iGain],
              d_g2xs[iGain], d_D2xs[iGain],
              d_D2ys[iGain], d_g2ys[iGain],
              d_flux_Is[iGain], d_flux_Qs[iGain],
              d_flux_Us[iGain], d_flux_Vs[iGain],
              d_visi_components[iGain],
              &visi_XX, &visi_XY,
              &visi_YX, &visi_YY);
      }
    } else {
      if (off_cardinal_dipoles == 1){
        apply_beam_gains_stokesI_off_cardinal_gpu(d_g1xs[iGain], d_D1xs[iGain],
              d_D1ys[iGain], d_g1ys[iGain],
              d_g2xs[iGain], d_D2xs[iGain],
              d_D2ys[iGain], d_g2ys[iGain],
              d_flux_Is[iGain],
              d_visi_components[iGain],
              &visi_XX, &visi_XY,
              &visi_YX, &visi_YY);
      } else {
        apply_beam_gains_stokesI_on_cardinal_gpu(d_g1xs[iGain], d_D1xs[iGain],
              d_D1ys[iGain], d_g1ys[iGain],
              d_g2xs[iGain], d_D2xs[iGain],
              d_D2ys[iGain], d_g2ys[iGain],
              d_flux_Is[iGain],
              d_visi_components[iGain],
              &visi_XX, &visi_XY,
              &visi_YX, &visi_YY);
      }
    }

    d_visi_XXs[iGain] = visi_XX;
    d_visi_XYs[iGain] = visi_XY;
    d_visi_YXs[iGain] = visi_YX;
    d_visi_YYs[iGain] = visi_YY;

  }
}

__global__ void kern_get_beam_gains(int num_components, int num_baselines,
           int num_freqs, int num_cross, int num_times, int beamtype,
           gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
           gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
           gpuUserComplex *d_recov_g1x, gpuUserComplex *d_recov_D1x,
           gpuUserComplex *d_recov_D1y, gpuUserComplex *d_recov_g1y,
           gpuUserComplex *d_recov_g2x, gpuUserComplex *d_recov_D2x,
           gpuUserComplex *d_recov_D2y, gpuUserComplex *d_recov_g2y,
           int use_twobeams, int num_ants,
           int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map) {

  // Start by computing which baseline we're going to do
  // This iBaseline means all baselines for all times and freqs
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  if(iBaseline < num_cross) {

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      gpuUserComplex g1x;
      gpuUserComplex D1x;
      gpuUserComplex D1y;
      gpuUserComplex g1y;
      gpuUserComplex g2x;
      gpuUserComplex D2x;
      gpuUserComplex D2y;
      gpuUserComplex g2y;

      if (use_twobeams == 1) {
        get_beam_gains_multibeams_gpu(iBaseline, iComponent, num_freqs,
                 num_baselines, num_components, num_times, beamtype,
                 d_g1xs, d_D1xs,
                 d_D1ys, d_g1ys,
                 d_ant1_to_baseline_map, d_ant2_to_baseline_map,
                 &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      } else {
        get_beam_gains_gpu(iBaseline, iComponent, num_freqs,
                 num_baselines, num_components, num_times, beamtype,
                 d_g1xs, d_D1xs,
                 d_D1ys, d_g1ys,
                 &g1x, &D1x, &D1y, &g1y, &g2x, &D2x, &D2y, &g2y);
      }

      int out_ind = num_cross*iComponent + iBaseline;

      d_recov_g1x[out_ind] = g1x;
      d_recov_D1x[out_ind] = D1x;
      d_recov_D1y[out_ind] = D1y;
      d_recov_g1y[out_ind] = g1y;
      d_recov_g2x[out_ind] = g2x;
      d_recov_D2x[out_ind] = D2x;
      d_recov_D2y[out_ind] = D2y;
      d_recov_g2y[out_ind] = g2y;

    }
  }
}

__global__ void kern_update_sum_visis_stokesIQUV(int num_freqs,
     int num_baselines, int num_components, int num_times,
     int beamtype, int off_cardinal_dipoles,
     gpuUserComplex *d_g1xs, gpuUserComplex *d_D1xs,
     gpuUserComplex *d_D1ys, gpuUserComplex *d_g1ys,
     int *d_ant1_to_baseline_map, int *d_ant2_to_baseline_map, int use_twobeams,
     gpuUserComplex *d_visi_components,
     user_precision_t *d_flux_I, user_precision_t *d_flux_Q,
     user_precision_t *d_flux_U, user_precision_t *d_flux_V,
     user_precision_t *d_sum_visi_XX_real, user_precision_t *d_sum_visi_XX_imag,
     user_precision_t *d_sum_visi_XY_real, user_precision_t *d_sum_visi_XY_imag,
     user_precision_t *d_sum_visi_YX_real, user_precision_t *d_sum_visi_YX_imag,
     user_precision_t *d_sum_visi_YY_real, user_precision_t *d_sum_visi_YY_imag) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);

  if(iBaseline < num_freqs*num_baselines*num_times) {

    int time_ind = (int)floorf( (user_precision_t)iBaseline / ((user_precision_t)num_baselines * (user_precision_t)num_freqs));
    int freq_ind = (int)floorf( ((user_precision_t)iBaseline - ((user_precision_t)time_ind*(user_precision_t)num_baselines * (user_precision_t)num_freqs)) / (user_precision_t)num_baselines);

    for (int iComponent = 0; iComponent < num_components; iComponent++) {

      //There is a flux for every frequnecy and component
      int flux_ind = num_components*freq_ind + iComponent;

      update_sum_visis_stokesIQUV_gpu(iBaseline, iComponent, num_freqs,
             num_baselines, num_components, num_times, beamtype, off_cardinal_dipoles,
             d_g1xs, d_D1xs,
             d_D1ys, d_g1ys,
             d_ant1_to_baseline_map, d_ant2_to_baseline_map, use_twobeams,
             d_visi_components[iBaseline],
             d_flux_I[flux_ind], d_flux_Q[flux_ind],
             d_flux_U[flux_ind], d_flux_V[flux_ind],
             d_sum_visi_XX_real, d_sum_visi_XX_imag,
             d_sum_visi_XY_real, d_sum_visi_XY_imag,
             d_sum_visi_YX_real, d_sum_visi_YX_imag,
             d_sum_visi_YY_real, d_sum_visi_YY_imag);

    }
  }
}