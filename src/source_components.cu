#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "cudacomplex.h"
#include "fundamental_coords.h"
#include "constants.h"
#include "shapelet_basis.h"
#include "read_and_write.h"
#include "source_components.h"

__device__ void extrap_flux(float *d_wavelengths, float *d_freqs,
           float *d_fluxes, int iComponent, int iBaseline,
           float * extrap_flux){

  float d_wavelength = d_wavelengths[iBaseline];
  float cat_wavelength = VELC / d_freqs[iComponent];
  * extrap_flux = d_fluxes[iComponent] * powf(cat_wavelength / d_wavelength,DEFAULT_SI);
}

__device__ void extrap_stokes(float *d_wavelengths, float *d_ref_freqs,
           float *d_ref_stokesI, float *d_ref_stokesQ,
           float *d_ref_stokesU, float *d_ref_stokesV,
           float *d_SIs, int iComponent, int iBaseline,
           float * flux_I, float * flux_Q, float * flux_U, float * flux_V){

  float d_freq = VELC / d_wavelengths[iBaseline];
  float d_ref_freq = d_ref_freqs[iComponent];

  float flux_ratio = powf(d_freq / d_ref_freq, d_SIs[iComponent]);

  * flux_I = d_ref_stokesI[iComponent] * flux_ratio;
  * flux_Q = d_ref_stokesQ[iComponent] * flux_ratio;
  * flux_U = d_ref_stokesU[iComponent] * flux_ratio;
  * flux_V = d_ref_stokesV[iComponent] * flux_ratio;

  if (iBaseline == 0) {
    // printf("Extrap fluxes %.5f %.5f %.5f %.5f\n", d_freq, d_ref_freq, d_SIs[iComponent], flux_ratio);
    // printf("Extrap fluxes %d %.5f %.5f %.5f %.5f\n",iComponent,* flux_I,* flux_Q,* flux_U,* flux_V );
  }
}

__global__ void kern_extrap_stokes(int num_visis, int num_components,
           float *d_wavelengths, float *d_ref_freqs, float *d_SIs,
           float *d_ref_stokesI, float *d_ref_stokesQ,
           float *d_ref_stokesU, float *d_ref_stokesV,
           float *d_flux_I, float *d_flux_Q,
           float *d_flux_U, float *d_flux_V ) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_components) {
    // printf("Here 2? %d %d\n",iBaseline,iComponent );

    float flux_I, flux_Q, flux_U, flux_V;

    extrap_stokes(d_wavelengths, d_ref_freqs,
                 d_ref_stokesI, d_ref_stokesQ,
                 d_ref_stokesU, d_ref_stokesV,
                 d_SIs, iComponent, iBaseline,
                 &flux_I, &flux_Q, &flux_U, &flux_V);

      d_flux_I[iComponent] = flux_I;
      d_flux_Q[iComponent] = flux_Q;
      d_flux_U[iComponent] = flux_U;
      d_flux_V[iComponent] = flux_V;

  }
}

extern "C" void test_extrap_flux(catsource_t *catsource,
           int num_visis, int num_components,
           float *wavelengths, float *flux_I, float *flux_Q,
           float *flux_U, float *flux_V ){

  float *d_wavelengths = NULL;
  float *d_ref_freqs = NULL;
  float *d_SIs = NULL;
  float *d_ref_stokesI = NULL;
  float *d_ref_stokesQ = NULL;
  float *d_ref_stokesU = NULL;
  float *d_ref_stokesV = NULL;
  float *d_flux_I = NULL;
  float *d_flux_Q = NULL;
  float *d_flux_U = NULL;
  float *d_flux_V = NULL;

  cudaMalloc( (void**)&d_wavelengths, num_visis*sizeof(float));
  cudaMalloc( (void**)&d_ref_freqs, num_components*sizeof(float));
  cudaMalloc( (void**)&d_SIs, num_components*sizeof(float));
  cudaMalloc( (void**)&d_ref_stokesI, num_components*sizeof(float));
  cudaMalloc( (void**)&d_ref_stokesQ, num_components*sizeof(float));
  cudaMalloc( (void**)&d_ref_stokesU, num_components*sizeof(float));
  cudaMalloc( (void**)&d_ref_stokesV, num_components*sizeof(float));
  cudaMalloc( (void**)&d_flux_I, num_components*sizeof(float));
  cudaMalloc( (void**)&d_flux_Q, num_components*sizeof(float));
  cudaMalloc( (void**)&d_flux_U, num_components*sizeof(float));
  cudaMalloc( (void**)&d_flux_V, num_components*sizeof(float));

  // printf("CUDA error 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
  cudaMemcpy(d_wavelengths, wavelengths, num_visis*sizeof(float), cudaMemcpyHostToDevice );

  // printf("CUDA error 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );
  cudaMemcpy(d_ref_stokesI, catsource->point_ref_stokesI, num_components*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_ref_stokesQ, catsource->point_ref_stokesQ, num_components*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_ref_stokesU, catsource->point_ref_stokesU, num_components*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_ref_stokesV, catsource->point_ref_stokesV, num_components*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_ref_freqs, catsource->point_ref_freqs, num_components*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_SIs, catsource->point_SIs, num_components*sizeof(float), cudaMemcpyHostToDevice );

  // float *zero_array = NULL;
  // // cudaMallocHost( (void**)&zero_array, num_time_steps*num_components*MAX_POLS*sizeof(float _Complex) );
  // zero_array = (float *)malloc( num_components*sizeof(float) );
  //
  // for (int i = 0; i < num_components; i++) {
  //   zero_array[i] = 0.0;
  // }
  //
  // cudaMemcpy(d_ref_stokesQ, zero_array, num_components*sizeof(float), cudaMemcpyHostToDevice );
  // cudaMemcpy(d_ref_stokesU, zero_array, num_components*sizeof(float), cudaMemcpyHostToDevice );
  // cudaMemcpy(d_ref_stokesV, zero_array, num_components*sizeof(float), cudaMemcpyHostToDevice );
  //
  // free(zero_array);

  // printf("Here 4?\n");

  dim3 grid, threads;

  threads.x = 64;
  threads.y = 2;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );
  grid.y = (int)ceil( (float)num_components / (float)threads.y );

  kern_extrap_stokes<<< grid, threads >>>(num_visis, num_components,
                     d_wavelengths, d_ref_freqs, d_SIs,
                     d_ref_stokesI, d_ref_stokesQ,
                     d_ref_stokesU, d_ref_stokesV,
                     d_flux_I, d_flux_Q,
                     d_flux_U, d_flux_V);

  cudaMemcpy(flux_I, d_flux_I, num_components*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(flux_Q, d_flux_Q, num_components*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(flux_U, d_flux_U, num_components*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(flux_V, d_flux_V, num_components*sizeof(float),cudaMemcpyDeviceToHost);

  cudaFree( d_wavelengths );
  cudaFree( d_ref_freqs );
  cudaFree( d_SIs );
  cudaFree( d_ref_stokesI );
  cudaFree( d_ref_stokesQ );
  cudaFree( d_ref_stokesU );
  cudaFree( d_ref_stokesV );
  cudaFree( d_flux_I );
  cudaFree( d_flux_Q );
  cudaFree( d_flux_U );
  cudaFree( d_flux_V );


}

__device__ cuFloatComplex calc_measurement_equation(float *d_us,
           float *d_vs, float *d_ws, float *d_ls, float *d_ms, float *d_ns,
           const int iBaseline, const int iComponent){

  float u, v, w;
  float l, m, n;

  u = d_us[iBaseline];
  v = d_vs[iBaseline];
  w = d_ws[iBaseline];

  l = d_ls[iComponent];
  m = d_ms[iComponent];
  n = d_ns[iComponent];

  cuFloatComplex visi;

  //Not sure why, but get exact match with oskar sims and correct location
  //on sky through wsclean without negative infront on 2pi
  float temp = 2*M_PI*( u*l + v*m + w*(n-1) );
  sincosf(temp, &(visi.y), &(visi.x));

  return visi;
}

__device__ void apply_beam_gains(cuFloatComplex g1xx, cuFloatComplex g1xy,
          cuFloatComplex g1yx, cuFloatComplex g1yy,
          cuFloatComplex g2xx, cuFloatComplex g2xy,
          cuFloatComplex g2yx, cuFloatComplex g2yy,
          float flux_I, float flux_Q,
          float flux_U, float flux_V,
          cuFloatComplex visi,
          cuFloatComplex * visi_XX, cuFloatComplex * visi_XY,
          cuFloatComplex * visi_YX, cuFloatComplex * visi_YY,
          int beamtype ) {

  //Conjugate the second beam gains
  cuFloatComplex g2xx_conj = make_cuFloatComplex(g2xx.x,-g2xx.y);
  cuFloatComplex g2xy_conj = make_cuFloatComplex(g2xy.x,-g2xy.y);
  cuFloatComplex g2yx_conj = make_cuFloatComplex(g2yx.x,-g2yx.y);
  cuFloatComplex g2yy_conj = make_cuFloatComplex(g2yy.x,-g2yy.y);

  //Create the Stokes visibilities
  cuFloatComplex visi_I = cuCmulf(make_cuComplex(flux_I,0.0),visi );
  cuFloatComplex visi_Q = cuCmulf(make_cuComplex(flux_Q,0.0),visi );
  cuFloatComplex visi_U = cuCmulf(make_cuComplex(flux_U,0.0),visi );
  cuFloatComplex visi_V = cuCmulf(make_cuComplex(flux_V,0.0),visi );

  cuFloatComplex this_XX;
  cuFloatComplex this_XY;
  cuFloatComplex this_YX;
  cuFloatComplex this_YY;

  //Convert the Stokes into instrumental visibilities
  this_XX = cuCmulf(cuCmulf(g1xx,g2xx_conj) + cuCmulf(g1xy,g2xy_conj),visi_I);
  this_XX += cuCmulf(cuCmulf(g1xx,g2xx_conj) - cuCmulf(g1xy,g2xy_conj),visi_Q);
  this_XX += cuCmulf(cuCmulf(g1xx,g2xy_conj) + cuCmulf(g1xy,g2xx_conj),visi_U);
  this_XX += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(g1xx,g2xy_conj) - cuCmulf(g1xy,g2xx_conj) );

  this_XY = cuCmulf(cuCmulf(g1xx,g2yx_conj) + cuCmulf(g1xy,g2yy_conj),visi_I);
  this_XY += cuCmulf(cuCmulf(g1xx,g2yx_conj) - cuCmulf(g1xy,g2yy_conj),visi_Q);
  this_XY += cuCmulf(cuCmulf(g1xx,g2yy_conj) + cuCmulf(g1xy,g2yx_conj),visi_U);
  this_XY += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(g1xx,g2yy_conj) - cuCmulf(g1xy,g2yx_conj) );

  this_YX = cuCmulf(cuCmulf(g1yx,g2xx_conj) + cuCmulf(g1yy,g2xy_conj),visi_I);
  this_YX += cuCmulf(cuCmulf(g1yx,g2xx_conj) - cuCmulf(g1yy,g2xy_conj),visi_Q);
  this_YX += cuCmulf(cuCmulf(g1yx,g2xy_conj) + cuCmulf(g1yy,g2xx_conj),visi_U);
  this_YX += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(g1yx,g2xy_conj) - cuCmulf(g1yy,g2xx_conj) );

  this_YY = cuCmulf(cuCmulf(g1yx,g2yx_conj) + cuCmulf(g1yy,g2yy_conj),visi_I);
  this_YY += cuCmulf(cuCmulf(g1yx,g2yx_conj) - cuCmulf(g1yy,g2yy_conj),visi_Q);
  this_YY += cuCmulf(cuCmulf(g1yx,g2yy_conj) + cuCmulf(g1yy,g2yx_conj),visi_U);
  this_YY += cuCmulf(cuCmulf(make_cuFloatComplex(0.0,1.0),visi_V), cuCmulf(g1yx,g2yy_conj) - cuCmulf(g1yy,g2yx_conj) );

  // if (beamtype == FEE_BEAM) {
  //   * visi_XX = this_YY;
  //   * visi_XY = this_YX;
  //   * visi_YX = this_XY;
  //   * visi_YY = this_XX;
  // }
  //
  // else {
  * visi_XX = this_XX;
  * visi_XY = this_XY;
  * visi_YX = this_YX;
  * visi_YY = this_YY;
  // }


}

__device__ void get_beam_gains(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           float *d_gauss_beam_reals, float *d_gauss_beam_imags,
           cuFloatComplex *d_analy_beam_X, cuFloatComplex *d_analy_beam_Y,
           cuFloatComplex *d_FEE_beam_gain_matrices,
           cuFloatComplex * g1xx, cuFloatComplex * g1xy,
           cuFloatComplex * g1yx, cuFloatComplex * g1yy,
           cuFloatComplex * g2xx, cuFloatComplex * g2xy,
           cuFloatComplex * g2yx, cuFloatComplex * g2yy){

  int beam_ind = 0;
  int time_ind = 0;
  int freq_ind = 0;

  if (beamtype == GAUSS_BEAM) {
    //Do some epic indexing to work out which beam value
    time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    cuFloatComplex gauss_beam_complex = make_cuFloatComplex(d_gauss_beam_reals[beam_ind],d_gauss_beam_imags[beam_ind]);

    * g1xx = gauss_beam_complex;
    * g2xx = gauss_beam_complex;
    * g1yy = gauss_beam_complex;
    * g2yy = gauss_beam_complex;
    * g1xy = make_cuComplex(0.0, 0.0);
    * g2xy = make_cuComplex(0.0, 0.0);
    * g1yx = make_cuComplex(0.0, 0.0);
    * g2yx = make_cuComplex(0.0, 0.0);

    // printf("%d %d %d %d %d %d %f %f\n",iBaseline,num_baselines,num_freqs,time_ind,freq_ind,beam_ind,beam_real,beam_imag);
  }

  else if (beamtype == FEE_BEAM) {
    // printf("Like, here? %d\n", iBaseline);
    //

    time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    beam_ind = iComponent*num_times + time_ind;

    * g1xx = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 0];
    * g1xy = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 1];
    * g1yx = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 2];
    * g1yy = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 3];

    * g2xx = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 0];
    * g2xy = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 1];
    * g2yx = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 2];
    * g2yy = d_FEE_beam_gain_matrices[beam_ind*MAX_POLS + 3];

    // cuFloatComplex thing1 = * g1xx;
    // cuFloatComplex thing2 = * g1xy;
    // cuFloatComplex thing3 = * g1yx;
    // cuFloatComplex thing4 = * g1yy;
    //
    // if (iBaseline == 0) {
    //   if (thing1.x > 2.0 || thing2.x > 2.0 || thing3.x > 2.0 || thing4.x > 2.0) {
    //     printf("%d %.5f %.5f %.5f %.5f\n",iComponent,thing1.x,thing2.x,thing3.x,thing4.x );
    //   }
    // }

  }

  else if (beamtype == ANALY_DIPOLE) {

    time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    * g1xx = d_analy_beam_X[beam_ind];
    * g2xx = d_analy_beam_X[beam_ind];
    * g1yy = d_analy_beam_Y[beam_ind];
    * g2yy = d_analy_beam_Y[beam_ind];

    * g1xy = make_cuComplex(0.0, 0.0);
    * g2xy = make_cuComplex(0.0, 0.0);
    * g1yx = make_cuComplex(0.0, 0.0);
    * g2yx = make_cuComplex(0.0, 0.0);

    cuFloatComplex thing1 = * g1xx;
    cuFloatComplex thing2 = * g1xy;
    cuFloatComplex thing3 = * g1yx;
    cuFloatComplex thing4 = * g1yy;

    // if (iBaseline == 0) {
    //   if (thing1.x > 2.0 || thing2.x > 2.0 || thing3.x > 2.0 || thing4.x > 2.0 || thing1.y > 2.0 || thing2.y > 2.0 || thing3.y > 2.0 || thing4.y > 2.0) {
    //         printf("%d %.5f %.5f %.5f %.5f\n",iComponent,thing1.x,thing2.x,thing3.x,thing4.x );
    //   }
    // }

  }

  else {
    * g1xx = make_cuComplex(1.0, 0.0);
    * g2xx = make_cuComplex(1.0, 0.0);
    * g1yy = make_cuComplex(1.0, 0.0);
    * g2yy = make_cuComplex(1.0, 0.0);
    * g1xy = make_cuComplex(0.0, 0.0);
    * g2xy = make_cuComplex(0.0, 0.0);
    * g1yx = make_cuComplex(0.0, 0.0);
    * g2yx = make_cuComplex(0.0, 0.0);
  }

} //end __device__ get_beam_gains

__device__ void update_sum_visis(int iBaseline, int iComponent, int num_freqs,
           int num_baselines, int num_components, int num_times, int beamtype,
           float *d_gauss_beam_reals, float *d_gauss_beam_imags,
           cuFloatComplex *d_FEE_beam_gain_matrices,
           cuFloatComplex *d_analy_beam_X, cuFloatComplex *d_analy_beam_Y,
           cuFloatComplex visi,
           float flux_I, float flux_Q, float flux_U, float flux_V,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag){

    cuFloatComplex g1xx;
    cuFloatComplex g1xy;
    cuFloatComplex g1yx;
    cuFloatComplex g1yy;
    cuFloatComplex g2xx;
    cuFloatComplex g2xy;
    cuFloatComplex g2yx;
    cuFloatComplex g2yy;

    get_beam_gains(iBaseline, iComponent, num_freqs,
               num_baselines, num_components, num_times, beamtype,
               d_gauss_beam_reals, d_gauss_beam_imags,
               d_analy_beam_X, d_analy_beam_Y,
               d_FEE_beam_gain_matrices,
               &g1xx, &g1xy, &g1yx, &g1yy, &g2xx, &g2xy, &g2yx, &g2yy);

    cuFloatComplex visi_XX;
    cuFloatComplex visi_XY;
    cuFloatComplex visi_YX;
    cuFloatComplex visi_YY;

    apply_beam_gains(g1xx, g1xy, g1yx, g1yy, g2xx, g2xy, g2yx, g2yy,
                    flux_I, flux_Q, flux_U, flux_V,
                    visi, &visi_XX, &visi_XY, &visi_YX, &visi_YY, beamtype );

    atomicAdd(&d_sum_visi_XX_real[iBaseline],visi_XX.x);
    atomicAdd(&d_sum_visi_XX_imag[iBaseline],visi_XX.y);

    atomicAdd(&d_sum_visi_XY_real[iBaseline],visi_XY.x);
    atomicAdd(&d_sum_visi_XY_imag[iBaseline],visi_XY.y);

    atomicAdd(&d_sum_visi_YX_real[iBaseline],visi_YX.x);
    atomicAdd(&d_sum_visi_YX_imag[iBaseline],visi_YX.y);

    atomicAdd(&d_sum_visi_YY_real[iBaseline],visi_YY.x);
    atomicAdd(&d_sum_visi_YY_imag[iBaseline],visi_YY.y);

}

__global__ void kern_calc_visi_point(float *d_point_ras, float *d_point_decs,
           float *d_point_freqs, float *d_point_stokesI, float *d_point_stokesQ,
           float *d_point_stokesU, float *d_point_stokesV, float *d_point_SIs,
           float *d_us, float *d_vs, float *d_ws,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
           float *d_angles_array, float *d_wavelengths,
           float *d_ls, float *d_ms, float *d_ns,
           int num_points, int num_baselines, int num_freqs, int num_visis,
           int num_times,
           float *d_gauss_beam_reals, float *d_gauss_beam_imags, int beamtype,
           cuFloatComplex *d_FEE_beam_gain_matrices,
           cuFloatComplex *d_analy_beam_X, cuFloatComplex *d_analy_beam_Y) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_points) {
    // float point_flux_I;
    // extrap_flux(d_wavelengths, d_point_freqs, d_point_fluxes,
    //               iComponent, iBaseline, &point_flux_I);
    //
    // float point_flux_Q = 0.0;
    // float point_flux_U = 0.0;
    // float point_flux_V = 0.0;

    float point_flux_I;
    float point_flux_Q;
    float point_flux_U;
    float point_flux_V;

    extrap_stokes(d_wavelengths, d_point_freqs,
                 d_point_stokesI, d_point_stokesQ,
                 d_point_stokesU, d_point_stokesV,
                 d_point_SIs, iComponent, iBaseline,
                 &point_flux_I, &point_flux_Q, &point_flux_U, &point_flux_V);

    // if (iBaseline == 0) {
    //   printf("FLUXES %.5f %.5f %.5f %.5f\n", point_flux_I, point_flux_Q, point_flux_U, point_flux_V);
    // }



    cuFloatComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws,
                           d_ls, d_ms, d_ns,
                           iBaseline, iComponent);

    update_sum_visis(iBaseline, iComponent, num_freqs,
           num_baselines, num_points, num_times, beamtype,
           d_gauss_beam_reals, d_gauss_beam_imags,
           d_FEE_beam_gain_matrices,
           d_analy_beam_X, d_analy_beam_Y,
           visi,
           point_flux_I, point_flux_Q, point_flux_U, point_flux_V,
           d_sum_visi_XX_real, d_sum_visi_XX_imag,
           d_sum_visi_XY_real, d_sum_visi_XY_imag,
           d_sum_visi_YX_real, d_sum_visi_YX_imag,
           d_sum_visi_YY_real, d_sum_visi_YY_imag);

  }
}

__global__ void kern_calc_visi_gaussian(float *d_gauss_ras, float *d_gauss_decs,
           float *d_gauss_freqs, float *d_gauss_stokesI, float *d_gauss_stokesQ,
           float *d_gauss_stokesU, float *d_gauss_stokesV, float *d_gauss_SIs,
           float *d_us, float *d_vs, float *d_ws,
           float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
           float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
           float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
           float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
           float *d_angles_array, float *d_wavelengths,
           float *d_ls, float *d_ms, float *d_ns,
           float *d_gauss_pas, float *d_gauss_majors, float *d_gauss_minors,
           int num_gauss, int num_baselines, int num_freqs, int num_visis,
           int num_times,
           float *d_gauss_beam_reals, float *d_gauss_beam_imags, int beamtype,
           cuFloatComplex *d_FEE_beam_gain_matrices) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_gauss) {

    float gauss_flux_I;
    float gauss_flux_Q;
    float gauss_flux_U;
    float gauss_flux_V;

    extrap_stokes(d_wavelengths, d_gauss_freqs,
                 d_gauss_stokesI, d_gauss_stokesQ,
                 d_gauss_stokesU, d_gauss_stokesV,
                 d_gauss_SIs, iComponent, iBaseline,
                 &gauss_flux_I, &gauss_flux_Q, &gauss_flux_U, &gauss_flux_V);

    cuFloatComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws,
                           d_ls, d_ms, d_ns,
                           iBaseline, iComponent);

    cuFloatComplex V_envelop = make_cuFloatComplex( 1.0, 0.0 );

    float pa = d_gauss_pas[iComponent];
    float u = d_us[iBaseline];
    float v = d_vs[iBaseline];
    float sinpa = sin(pa);
    float cospa = cos(pa);

    float x =  cospa*v + sinpa*u; // major axis
    float y = -sinpa*v + cospa*u; // minor axis
    float invsig_x = d_gauss_majors[iComponent];
    float invsig_y = d_gauss_minors[iComponent];

    V_envelop = make_cuFloatComplex( exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ), 0.0 );

    visi = cuCmulf(visi, V_envelop);

    cuFloatComplex *d_analy_beam_X = NULL;
    cuFloatComplex *d_analy_beam_Y = NULL;

    update_sum_visis(iBaseline, iComponent, num_freqs,
           num_baselines, num_gauss, num_times, beamtype,
           d_gauss_beam_reals, d_gauss_beam_imags,
           d_FEE_beam_gain_matrices,
           d_analy_beam_X, d_analy_beam_Y,
           visi,
           gauss_flux_I, gauss_flux_Q, gauss_flux_U, gauss_flux_V,
           d_sum_visi_XX_real, d_sum_visi_XX_imag,
           d_sum_visi_XY_real, d_sum_visi_XY_imag,
           d_sum_visi_YX_real, d_sum_visi_YX_imag,
           d_sum_visi_YY_real, d_sum_visi_YY_imag);

  }
}

__global__ void kern_calc_visi_shapelets(float *d_shape_ras,
      float *d_shape_decs, float *d_shape_fluxes, float *d_shape_freqs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_wavelengths,
      float *d_u_s_metres, float *d_v_s_metres, float *d_w_s_metres,
      float *d_sum_visi_XX_real, float *d_sum_visi_XX_imag,
      float *d_sum_visi_XY_real, float *d_sum_visi_XY_imag,
      float *d_sum_visi_YX_real, float *d_sum_visi_YX_imag,
      float *d_sum_visi_YY_real, float *d_sum_visi_YY_imag,
      float *d_angles_array, float *d_shape_pas, float *d_shape_majors,
      float *d_shape_minors,
      float *d_shape_n1s, float *d_shape_n2s, float *d_shape_coeffs,
      float *d_shape_param_indexes,
      float *d_shape_ls, float *d_shape_ms, float *d_shape_ns,
      float *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_visis,
      const int num_coeffs, int num_times,
      float *d_gauss_beam_reals, float *d_gauss_beam_imags, int beamtype,
      cuFloatComplex *d_FEE_beam_gain_matrices) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iCoeff = threadIdx.y + (blockDim.y*blockIdx.y);


  if (iBaseline < num_visis && iCoeff < num_coeffs) {
    // printf("Made it here all g %d %d %d %d\n", iBaseline,num_visis,iCoeff,num_coeffs);
    int iComponent = d_shape_param_indexes[iCoeff];

    float shape_flux_I;
    //Use param index below and not iCoeff as there
    extrap_flux(d_wavelengths, d_shape_freqs,
               d_shape_fluxes, iComponent, iBaseline,
               &shape_flux_I);

    float shape_flux_Q = 0.0;
    float shape_flux_U = 0.0;
    float shape_flux_V = 0.0;

    cuFloatComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws,
                          d_shape_ls, d_shape_ms, d_shape_ns,
                          iBaseline, iComponent);

    float pa = d_shape_pas[iComponent];
    float sinpa = sin(pa);
    float cospa = cos(pa);

    float d_wavelength = d_wavelengths[iBaseline];

    float u_s = d_u_s_metres[iComponent*num_visis + iBaseline] / d_wavelength;
    float v_s = d_v_s_metres[iComponent*num_visis + iBaseline] / d_wavelength;
    //
    float x = (cospa*v_s + sinpa*u_s); // major axis
    float y = (-sinpa*v_s + cospa*u_s); // minor axis

    //Scales the FWHM to std to match basis functions, and account for the
    //basis functions being stored with beta = 1.0
    //Basis functions have been stored in such a way that x is in the same
    //direction as on sky, but y is opposite, so include negative here
    float const_x = (d_shape_majors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
    float const_y = -(d_shape_minors[iComponent]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

    // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
    cuFloatComplex Ipow_lookup[] = { make_cuFloatComplex(  1.0,  0.0 ),
                                     make_cuFloatComplex(  0.0,  1.0 ),
                                     make_cuFloatComplex( -1.0,  0.0 ),
                                     make_cuFloatComplex(  0.0, -1.0 ) };

    float xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

    // find the indices in the basis functions for u*beta_u and v*beta_v

    float xpos = x*const_x + sbf_c;
    float ypos = y*const_y + sbf_c;

    int xindex = (int)floor(xpos);
    int yindex = (int)floor(ypos);
    //
    int n1 = (int)d_shape_n1s[iCoeff];
    int n2 = (int)d_shape_n2s[iCoeff];

    // if ( n1 < 0 || n2 < 0 || n1 >= sbf_N || n2 >= sbf_N ) continue;

    f_hat = d_shape_coeffs[iCoeff];
    //
    sbf_n = &d_sbf[n1*sbf_L];
    xlow  = sbf_n[xindex];
    xhigh = sbf_n[xindex+1];
    u_value = xlow + (xhigh-xlow)*(xpos-xindex);

    sbf_n = &d_sbf[n2*sbf_L];
    ylow  = sbf_n[yindex];
    yhigh = sbf_n[yindex+1];
    v_value = ylow + (yhigh-ylow)*(ypos-yindex);

    // accumulate the intensity model for baseline pair (u,v)
    cuFloatComplex V_envelop = make_cuFloatComplex( 0.0, 0.0 );
    V_envelop = V_envelop + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

    visi = cuCmulf(visi, V_envelop);

    cuFloatComplex *d_analy_beam_X = NULL;
    cuFloatComplex *d_analy_beam_Y = NULL;

    update_sum_visis(iBaseline, iComponent, num_freqs,
           num_baselines, num_shapes, num_times, beamtype,
           d_gauss_beam_reals, d_gauss_beam_imags,
           d_FEE_beam_gain_matrices,
           d_analy_beam_X, d_analy_beam_Y,
           visi,
           shape_flux_I, shape_flux_Q, shape_flux_U, shape_flux_V,
           d_sum_visi_XX_real, d_sum_visi_XX_imag,
           d_sum_visi_XY_real, d_sum_visi_XY_imag,
           d_sum_visi_YX_real, d_sum_visi_YX_imag,
           d_sum_visi_YY_real, d_sum_visi_YY_imag);

  }

}
