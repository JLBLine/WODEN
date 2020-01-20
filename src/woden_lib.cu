#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "woden_lib.h"
#include "shapelet_basis.h"
#include "cudacomplex.h"

__global__ void calc_uvw(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
      float *d_u_metres, float *d_v_metres, float *d_w_metres,
      float *d_angles_array, float *d_cha0s, float *d_sha0s) {
  //

  float u, v, w;

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  float d_sdec0 = d_angles_array[0];
  float d_cdec0 = d_angles_array[1];
  float d_sha0 = d_sha0s[iBaseline];
  float d_cha0 = d_cha0s[iBaseline];

  int mod_baseline = iBaseline - 8128*floorf(iBaseline / 8128.0);

  u = (d_sha0*d_X_diff[mod_baseline]) + (d_cha0*d_Y_diff[mod_baseline]);
  v = (d_sdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_cdec0*d_Z_diff[mod_baseline]) - (d_sdec0*d_cha0*d_X_diff[mod_baseline]);
  w = (d_cdec0*d_cha0*d_X_diff[mod_baseline]) - (d_cdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_sdec0*d_Z_diff[mod_baseline]);

  d_u_metres[iBaseline] = u;
  d_v_metres[iBaseline] = v;
  d_w_metres[iBaseline] = w;

}

__global__ void calc_uvw_shapelet(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
      float *d_u_s_metres, float *d_v_s_metres, float *d_w_s_metres,
      float *d_lsts, float *d_ras, float *d_decs,
      const int num_baselines, const int num_visis) {

  float u_s, v_s, w_s;

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  //TODO do the sin/cos outside of the GPU kernel?
  float d_sdec0 = sinf(d_decs[iComponent]);
  float d_cdec0 = cosf(d_decs[iComponent]);
  float d_sha0 = sinf(d_lsts[iBaseline] - d_ras[iComponent]);
  float d_cha0 = cosf(d_lsts[iBaseline] - d_ras[iComponent]);

  int mod_baseline = iBaseline - 8128*floorf(iBaseline / 8128.0);

  u_s = (d_sha0*d_X_diff[mod_baseline]) + (d_cha0*d_Y_diff[mod_baseline]);
  v_s = (d_sdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_cdec0*d_Z_diff[mod_baseline]) - (d_sdec0*d_cha0*d_X_diff[mod_baseline]);
  w_s = (d_cdec0*d_cha0*d_X_diff[mod_baseline]) - (d_cdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_sdec0*d_Z_diff[mod_baseline]);

  d_u_s_metres[num_visis*iComponent + iBaseline] = u_s;
  d_v_s_metres[num_visis*iComponent + iBaseline] = v_s;
  d_w_s_metres[num_visis*iComponent + iBaseline] = w_s;

}

__global__ void calc_lmn(float *d_angles_array, float *d_ras, float *d_decs,
                         float *d_l, float *d_m, float *d_n){
  const int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);

  float d_sdec0 = d_angles_array[0];
  float d_cdec0 = d_angles_array[1];
  float d_ra0 = d_angles_array[2];

  float cdec;
  float sdec;
  float cdra;
  float sdra;

  cdec = cosf(d_decs[iComponent]);
  sdec = sinf(d_decs[iComponent]);
  cdra = cosf((d_ras[iComponent] - d_ra0));
  sdra = sinf((d_ras[iComponent] - d_ra0));

  d_l[iComponent] = cdec*sdra;
  d_m[iComponent] = sdec*d_cdec0 - cdec*d_sdec0*cdra;
  d_n[iComponent] = sdec*d_sdec0 + cdec*d_cdec0*cdra;

}

__device__ void extrap_uvw_flux_calc_lmn(float *d_angles_array,
                float *d_ras, float *d_decs, float *d_fluxes, float *d_freqs,
                float *d_u_metres, float *d_v_metres, float *d_w_metres, float *d_wavelengths,
                int iComponent, int iBaseline,
                float * l, float * m, float * n, float * u, float * v, float * w, float * extrap_flux){
  float cdec;
  float sdec;
  float cdra;
  float sdra;

  float d_sdec0 = d_angles_array[0];
  float d_cdec0 = d_angles_array[1];
  float d_ra0 = d_angles_array[2];
  float d_wavelength = d_wavelengths[iBaseline];

  cdec = cosf(d_decs[iComponent]);
  sdec = sinf(d_decs[iComponent]);
  cdra = cosf((d_ras[iComponent]-d_ra0));
  sdra = sinf((d_ras[iComponent]-d_ra0));

  * l = cdec*sdra;
  * m = sdec*d_cdec0 - cdec*d_sdec0*cdra;
  * n = sdec*d_sdec0 + cdec*d_cdec0*cdra;

  * u = d_u_metres[iBaseline] / d_wavelength;
  * v = d_v_metres[iBaseline] / d_wavelength;
  * w = d_w_metres[iBaseline] / d_wavelength;

  float cat_wavelength = VELC / d_freqs[iComponent];

  * extrap_flux = d_fluxes[iComponent] * powf(cat_wavelength / d_wavelength,DEFAULT_SI);

}

__device__ void extrap_uvw_flux(float *d_angles_array,
                float *d_u_metres, float *d_v_metres, float *d_w_metres, float *d_wavelengths,
                float *d_freqs, float *d_fluxes,
                int iComponent, int iBaseline, int param_index,
                float * u, float * v, float * w, float * extrap_flux){

  float d_wavelength = d_wavelengths[iBaseline];

  * u = d_u_metres[iBaseline] / d_wavelength;
  * v = d_v_metres[iBaseline] / d_wavelength;
  * w = d_w_metres[iBaseline] / d_wavelength;

  float cat_wavelength = VELC / d_freqs[param_index];
  * extrap_flux = d_fluxes[param_index] * powf(cat_wavelength / d_wavelength,DEFAULT_SI);

}


__global__ void calc_visi_point(float *d_point_ras, float *d_point_decs, float *d_point_fluxes, float *d_point_freqs,
      float *d_u_metres, float *d_v_metres, float *d_w_metres,
      float *d_sum_visi_real, float *d_sum_visi_imag,
      float *d_angles_array, float *d_wavelengths) {
  //

  float u, v, w;
  float l, m, n;
  float extrap_flux;

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  //
  extrap_uvw_flux_calc_lmn(d_angles_array,
                  d_point_ras, d_point_decs, d_point_fluxes, d_point_freqs,
                  d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
                  iComponent, iBaseline,
                  &l, &m, &n, &u, &v, &w, &extrap_flux);

  cuFloatComplex visi;

  //Not sure why, but get exact match with oskar sims and correct location
  //on sky through wsclean without negative infront on 2pi
  float temp = 2*M_PI*( u*l + v*m + w*(n-1) );
  sincosf(temp, &(visi.y), &(visi.x));

  register float i = atomicAdd(&d_sum_visi_real[iBaseline],visi.x*extrap_flux);
  register float j = atomicAdd(&d_sum_visi_imag[iBaseline],visi.y*extrap_flux);

  //HERE
  // register float i = atomicAdd(&d_sum_visi_real[iBaseline],0.0);
  // register float j = atomicAdd(&d_sum_visi_imag[iBaseline],0.0);

}

__global__ void calc_visi_gaussian(float *d_gauss_ras, float *d_gauss_decs, float *d_gauss_fluxes, float *d_gauss_freqs,
      float *d_u_metres, float *d_v_metres, float *d_w_metres, float *d_wavelengths,
      float *d_sum_visi_real, float *d_sum_visi_imag,
      float *d_angles_array,
      float *d_gauss_pas, float *d_gauss_majors, float *d_gauss_minors ) {
  //

  float u, v, w;
  float l, m, n;
  float extrap_flux;

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  //
  extrap_uvw_flux_calc_lmn(d_angles_array,
                  d_gauss_ras, d_gauss_decs, d_gauss_fluxes, d_gauss_freqs,
                  d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
                  iComponent, iBaseline,
                  &l, &m, &n, &u, &v, &w, &extrap_flux);

  cuFloatComplex visi;

  //Not sure why, but get exact match with oskar sims and correct location
  //on sky through wsclean without negative infront on 2pi

  float temp = 2*M_PI*( u*l + v*m + w*(n-1) );
  sincosf(temp, &(visi.y), &(visi.x));

  cuFloatComplex V_envelop = make_cuFloatComplex( 1.0, 0.0 );

  float pa = d_gauss_pas[iComponent];
  float sinpa = sin(pa);
  float cospa = cos(pa);

  float x =  cospa*v + sinpa*u; // major axis
  float y = -sinpa*v + cospa*u; // minor axis
  float invsig_x = d_gauss_majors[iComponent];
  float invsig_y = d_gauss_minors[iComponent];

  V_envelop = make_cuFloatComplex( exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ), 0.0 );

  visi = cuCmulf(visi, V_envelop);

  register float i = atomicAdd(&d_sum_visi_real[iBaseline],visi.x * extrap_flux);
  register float j = atomicAdd(&d_sum_visi_imag[iBaseline],visi.y * extrap_flux);

  //HERE
  // register float i = atomicAdd(&d_sum_visi_real[iBaseline],0.0);
  // register float j = atomicAdd(&d_sum_visi_imag[iBaseline],0.0);

}

__global__ void calc_visi_shapelets(float *d_shape_ras, float *d_shape_decs, float *d_shape_fluxes, float *d_shape_freqs,
      float *d_u_metres, float *d_v_metres, float *d_w_metres, float *d_wavelengths,
      float *d_u_s_metres, float *d_v_s_metres, float *d_w_s_metres,
      float *d_sum_visi_real, float *d_sum_visi_imag,
      float *d_angles_array,float *d_shape_pas, float *d_shape_majors, float *d_shape_minors,
      float *d_shape_n1s, float *d_shape_n2s, float *d_shape_coeffs,float *d_shape_param_indexes,
      float *d_shape_ls, float *d_shape_ms, float *d_shape_ns,
      float *d_sbf,
      const int num_baselines, const int num_visis){


  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);


  float u, v, w;
  float l, m, n;
  float extrap_flux;

  int param_index = d_shape_param_indexes[iComponent];

  extrap_uvw_flux(d_angles_array,
                  d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
                  d_shape_freqs, d_shape_fluxes,
                  iComponent, iBaseline, param_index,
                  &u, &v, &w, &extrap_flux);
  float pa = d_shape_pas[param_index];
  float sinpa = sin(pa);
  float cospa = cos(pa);

  float d_wavelength = d_wavelengths[iBaseline];

  float u_s = d_u_s_metres[param_index*num_visis + iBaseline] / d_wavelength;
  float v_s = d_v_s_metres[param_index*num_visis + iBaseline] / d_wavelength;
  //
  float x = (cospa*v_s + sinpa*u_s); // major axis
  float y = (-sinpa*v_s + cospa*u_s); // minor axis

  //Scales the FWHM to std to match basis functions, and account for the
  //basis functions being stored with beta = 1.0
  //Basis functions have been stored in such a way that x is in the same
  //direction as on sky, but y is opposite, so include negative here
  float const_x = (d_shape_majors[param_index]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
  float const_y = -(d_shape_minors[param_index]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

  // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
  cuFloatComplex Ipow_lookup[] = { make_cuFloatComplex(  1.0,  0.0 ),
                                   make_cuFloatComplex(  0.0,  1.0 ),
                                   make_cuFloatComplex( -1.0,  0.0 ),
                                   make_cuFloatComplex(  0.0, -1.0 ) };
  //
  float xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

  // find the indices in the basis functions for u*beta_u and v*beta_v

  float xpos = x*const_x + sbf_c;
  float ypos = y*const_y + sbf_c;

  int xindex = (int)floor(xpos);
  int yindex = (int)floor(ypos);
  //
  // loop over shapelet coefficients, building up the intensity model, if this baseline is in range

  // if ( xindex >= 0 && yindex >= 0 && xindex+1 < sbf_L && yindex+1 < sbf_L ) continue;

  int n1 = (int)d_shape_n1s[iComponent];
  int n2 = (int)d_shape_n2s[iComponent];

  // if ( n1 < 0 || n2 < 0 || n1 >= sbf_N || n2 >= sbf_N ) continue;

  f_hat = d_shape_coeffs[iComponent];
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
  //

  l = d_shape_ls[param_index];
  m = d_shape_ms[param_index];
  n = d_shape_ns[param_index];

  cuFloatComplex visi;
  //Not sure why, but get exact match with oskar sims and correct location
  //on sky through wsclean without negative infront on 2pi

  float temp = 2*M_PI*( u*l + v*m + w*(n-1) );
  sincosf(temp, &(visi.y), &(visi.x));

  visi = cuCmulf(visi, V_envelop);

  register float i = atomicAdd(&d_sum_visi_real[iBaseline],visi.x * extrap_flux);
  register float j = atomicAdd(&d_sum_visi_imag[iBaseline],visi.y * extrap_flux);
  //
  //HERE
  // register float i = atomicAdd(&d_sum_visi_real[iBaseline],0.0);
  // register float j = atomicAdd(&d_sum_visi_imag[iBaseline],0.0);

}

extern "C" void Atomic_time_step(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array,
                    const int num_baselines, const int num_visis,
                    visibility_set_t *visibility_set,
                    float *sbf2) {

  /*START We should be able to do all this outside of this function and transfer in--------------*/
  //==========================================================================================
  float *d_X_diff = NULL;
  float *d_Y_diff = NULL;
  float *d_Z_diff = NULL;

  cudaMalloc( (void**)&d_X_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_X_diff, X_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Y_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Y_diff, Y_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Z_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Z_diff, Z_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );

  float *d_angles_array;
  cudaMalloc( (void**)&d_angles_array, 3*sizeof(float) );
  cudaMemcpy( d_angles_array, angles_array, 3*sizeof(float), cudaMemcpyHostToDevice );

  float *d_sha0s = NULL;
  float *d_cha0s = NULL;
  float *d_wavelengths = NULL;
  cudaMalloc( (void**)&d_sha0s, num_visis*sizeof(float) );
  cudaMemcpy( d_sha0s, visibility_set->sha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_cha0s, num_visis*sizeof(float) );
  cudaMemcpy( d_cha0s, visibility_set->cha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_wavelengths, num_visis*sizeof(float) );
  cudaMemcpy( d_wavelengths, visibility_set->wavelengths, num_visis*sizeof(float), cudaMemcpyHostToDevice );

  //
  /* END We should be able to do all this outside of this function and transfer in--------------*/

  float *d_sum_visi_real = NULL;
  float *d_sum_visi_imag = NULL;
  float *d_u_metres = NULL;
  float *d_v_metres = NULL;
  float *d_w_metres = NULL;

  cudaMalloc( (void**)&d_sum_visi_real, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_sum_visi_imag, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_u_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_v_metres, num_visis*sizeof(float) );
  cudaMalloc( (void**)&d_w_metres, num_visis*sizeof(float) );


  dim3 grid, threads;

  threads.x = 128;
  threads.y = 1;
  threads.z = 1;
  grid.x = (int)ceil( (float)num_visis / (float)threads.x );
  // grid.x = 127;
  grid.y = 1;
  grid.z = 1;

  calc_uvw<<< grid , threads, 0 >>>( d_X_diff,
          d_Y_diff, d_Z_diff,
          d_u_metres, d_v_metres, d_w_metres,
          d_angles_array, d_cha0s, d_sha0s);


  //
  int num_points = catsource.n_points;
  int num_gauss = catsource.n_gauss;
  int num_shapes = catsource.n_shapes;

  if (num_points > 0) {

    float *d_point_ras=NULL;
    float *d_point_decs=NULL;
    float *d_point_fluxes=NULL;
    float *d_point_freqs=NULL;

    cudaMalloc( (void**)&(d_point_ras), num_points*sizeof(float) );
    cudaMemcpy( d_point_ras, catsource.point_ras, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_decs), num_points*sizeof(float) );
    cudaMemcpy( d_point_decs, catsource.point_decs, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_fluxes), num_points*sizeof(float) );
    cudaMemcpy( d_point_fluxes, catsource.point_fluxes, num_points*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_point_freqs), num_points*sizeof(float) );
    cudaMemcpy( d_point_freqs, catsource.point_freqs, num_points*sizeof(float), cudaMemcpyHostToDevice );

    //TODO need to put some kind of check inside calc_visi_point to skip a grid
    //thread point if it's outside the realms of sanity

    if (num_points == 1) {
      threads.x = 64;
      threads.y = 1;
      grid.x = grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;

    }

    else {
      threads.x = 64;
      threads.y = 2;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)num_points) / ((float)threads.y) );

    }

    calc_visi_point<<< grid , threads, 0 >>>(d_point_ras,
            d_point_decs, d_point_fluxes,d_point_freqs,
            d_u_metres, d_v_metres, d_w_metres,
            d_sum_visi_real, d_sum_visi_imag,
            d_angles_array, d_wavelengths);

    cudaFree( d_point_freqs );
    cudaFree( d_point_fluxes );
    cudaFree( d_point_decs);
    cudaFree( d_point_ras);

  }//if point sources

  if (num_gauss > 0) {

    float *d_gauss_ras=NULL;
    float *d_gauss_decs=NULL;
    float *d_gauss_fluxes=NULL;
    float *d_gauss_freqs=NULL;
    float *d_gauss_pas=NULL;
    float *d_gauss_majors=NULL;
    float *d_gauss_minors=NULL;

    cudaMalloc( (void**)&(d_gauss_ras), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_ras, catsource.gauss_ras, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_decs), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_decs, catsource.gauss_decs, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_fluxes), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_fluxes, catsource.gauss_fluxes, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_pas), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_pas, catsource.gauss_pas, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_majors), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_majors, catsource.gauss_majors, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_minors), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_minors, catsource.gauss_minors, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_gauss_freqs), num_gauss*sizeof(float) );
    cudaMemcpy( d_gauss_freqs, catsource.gauss_freqs, num_gauss*sizeof(float), cudaMemcpyHostToDevice );

    //TODO need to put some kind of check inside calc_visi_point to skip a grid
    //thread point if it's outside the realms of sanity
    threads.x = 64;
    threads.y = 2;
    grid.x = (int)ceil( (float)num_visis / (float)threads.x );
    grid.y = (int)ceil( ((float)num_gauss) / ((float)threads.y) );

    calc_visi_gaussian<<< grid , threads, 0 >>>(d_gauss_ras,
            d_gauss_decs, d_gauss_fluxes, d_gauss_freqs,
            d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
            d_sum_visi_real, d_sum_visi_imag,
            d_angles_array,
            d_gauss_pas, d_gauss_majors, d_gauss_minors);

    cudaFree( d_gauss_minors );
    cudaFree( d_gauss_majors);
    cudaFree( d_gauss_pas);
    cudaFree( d_gauss_freqs );
    cudaFree( d_gauss_fluxes );
    cudaFree( d_gauss_decs);
    cudaFree( d_gauss_ras);

  }//if gauss sources

  if (num_shapes > 0) {
    float *d_shape_ras=NULL;
    float *d_shape_decs=NULL;
    float *d_shape_fluxes=NULL;
    float *d_shape_freqs=NULL;
    float *d_shape_pas=NULL;
    float *d_shape_majors=NULL;
    float *d_shape_minors=NULL;

    float *d_shape_coeffs=NULL;
    float *d_shape_n1s=NULL;
    float *d_shape_n2s=NULL;
    float *d_shape_param_indexes=NULL;

    float *d_sbf=NULL;
    float *d_lsts=NULL;

    float *d_u_s_metres = NULL;
    float *d_v_s_metres = NULL;
    float *d_w_s_metres = NULL;

    float *d_shape_ls=NULL;
    float *d_shape_ms=NULL;
    float *d_shape_ns=NULL;

    //Who likes cudaMalloc cudaMallocs? We like cudaMalloc cudaMallocs
    cudaMalloc( (void**)&(d_shape_ras), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_ras, catsource.shape_ras, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_decs), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_decs, catsource.shape_decs, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_fluxes), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_fluxes, catsource.shape_fluxes, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_pas), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_pas, catsource.shape_pas, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_majors), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_majors, catsource.shape_majors, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_minors), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_minors, catsource.shape_minors, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_freqs), num_shapes*sizeof(float) );
    cudaMemcpy( d_shape_freqs, catsource.shape_freqs, num_shapes*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_coeffs), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_coeffs, catsource.shape_coeffs, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_n1s), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_n1s, catsource.shape_n1s, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_n2s), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_n2s, catsource.shape_n2s, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_param_indexes), catsource.n_shape_coeffs*sizeof(float) );
    cudaMemcpy( d_shape_param_indexes, catsource.shape_param_indexes, catsource.n_shape_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_sbf), sbf_N*sbf_L*sizeof(float) );
    cudaMemcpy( d_sbf, sbf2, sbf_N*sbf_L*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_shape_ls), num_shapes*sizeof(float) );
    cudaMalloc( (void**)&(d_shape_ms), num_shapes*sizeof(float) );
    cudaMalloc( (void**)&(d_shape_ns), num_shapes*sizeof(float) );

    cudaMalloc( (void**)&(d_lsts), num_visis*sizeof(float) );
    cudaMemcpy( d_lsts, visibility_set->lsts, num_visis*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&d_u_s_metres, num_shapes*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_v_s_metres, num_shapes*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_w_s_metres, num_shapes*num_visis*sizeof(float) );


    if (num_shapes == 1) {
      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;

      calc_uvw_shapelet<<< grid , threads, 0 >>>(d_X_diff, d_Y_diff, d_Z_diff,
            d_u_s_metres, d_v_s_metres, d_w_s_metres,
            d_lsts, d_shape_ras, d_shape_decs,
            num_baselines, num_visis);

    }
    else {
      threads.x = 64;
      threads.y = 2;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)num_shapes) / ((float)threads.y) );

      calc_uvw_shapelet<<< grid , threads, 0 >>>(d_X_diff, d_Y_diff, d_Z_diff,
            d_u_s_metres, d_v_s_metres, d_w_s_metres,
            d_lsts, d_shape_ras, d_shape_decs,
            num_baselines, num_visis);

    }


    threads.x = 64;
    threads.y = 1;
    grid.x = (int)ceil( ((float)num_shapes / (float)threads.x) );
    grid.y = 1;

    calc_lmn<<< grid , threads, 0 >>>(d_angles_array, d_shape_ras, d_shape_decs,
                             d_shape_ls, d_shape_ms, d_shape_ns);

    //TODO need to put some kind of check inside calc_visi_point to skip a grid
    //thread point if it's outside the realms of sanity
    if (catsource.n_shape_coeffs == 1) {
      threads.x = 64;
      threads.y = 1;
      threads.z = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
      grid.z = 1;

      calc_visi_shapelets<<< grid , threads, 0 >>>(d_shape_ras,
              d_shape_decs, d_shape_fluxes, d_shape_freqs,
              d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
              d_u_s_metres, d_v_s_metres, d_w_s_metres,
              d_sum_visi_real, d_sum_visi_imag,
              d_angles_array, d_shape_pas, d_shape_majors, d_shape_minors,
              d_shape_n1s, d_shape_n2s, d_shape_coeffs, d_shape_param_indexes,
              d_shape_ls, d_shape_ms, d_shape_ns,
              d_sbf,
              num_baselines, num_visis);
    }

    else {

      threads.x = 64;
      threads.y = 2;
      threads.z = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)catsource.n_shape_coeffs) / ((float)threads.y) );
      grid.z = 1;

      calc_visi_shapelets<<< grid , threads, 0 >>>(d_shape_ras,
              d_shape_decs, d_shape_fluxes, d_shape_freqs,
              d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
              d_u_s_metres, d_v_s_metres, d_w_s_metres,
              d_sum_visi_real, d_sum_visi_imag,
              d_angles_array, d_shape_pas, d_shape_majors, d_shape_minors,
              d_shape_n1s, d_shape_n2s, d_shape_coeffs, d_shape_param_indexes,
              d_shape_ls, d_shape_ms, d_shape_ns,
              d_sbf,
              num_baselines, num_visis);
    }

    cudaFree( d_shape_ns );
    cudaFree( d_shape_ms );
    cudaFree( d_shape_ls );
    cudaFree(d_w_s_metres);
    cudaFree(d_v_s_metres);
    cudaFree(d_u_s_metres);
    cudaFree(d_lsts);
    cudaFree(d_sbf);
    cudaFree(d_shape_param_indexes);
    cudaFree(d_shape_n2s);
    cudaFree(d_shape_n1s);
    cudaFree(d_shape_coeffs);
    cudaFree( d_shape_minors );
    cudaFree( d_shape_majors);
    cudaFree( d_shape_pas);
    cudaFree( d_shape_freqs );
    cudaFree( d_shape_fluxes );
    cudaFree( d_shape_decs);
    cudaFree( d_shape_ras);

  }//if shapelet

  //Get the results into host memory
  cudaMemcpy(visibility_set->sum_visi_real,d_sum_visi_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->sum_visi_imag,d_sum_visi_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set->us_metres,d_u_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->vs_metres,d_v_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set->ws_metres,d_w_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  //Free up the GPU memory
  cudaFree( d_w_metres );
  cudaFree( d_v_metres );
  cudaFree( d_u_metres );
  cudaFree( d_sum_visi_imag );
  cudaFree( d_sum_visi_real );
  cudaFree( d_wavelengths );
  cudaFree( d_cha0s );
  cudaFree( d_sha0s );
  cudaFree( d_angles_array );
  cudaFree( d_Z_diff );
  cudaFree( d_Y_diff );
  cudaFree( d_X_diff );

}

extern "C" void assign_pointsource_on_GPU(catsource_t src){

  src.d_point_ras=NULL;
  src.d_point_decs=NULL;
  src.d_point_fluxes=NULL;
  src.d_point_freqs=NULL;

  int num_comps = src.n_points;

  cudaMalloc( (void**)&(src.d_point_ras), num_comps*sizeof(float) );
  cudaMemcpy( src.d_point_ras, src.point_ras, num_comps*sizeof(float), cudaMemcpyHostToDevice );

  cudaMalloc( (void**)&(src.d_point_decs), num_comps*sizeof(float) );
  cudaMemcpy( src.d_point_decs, src.point_decs, num_comps*sizeof(float), cudaMemcpyHostToDevice );

  cudaMalloc( (void**)&(src.d_point_fluxes), num_comps*sizeof(float) );
  cudaMemcpy( src.d_point_fluxes, src.point_fluxes, num_comps*sizeof(float), cudaMemcpyHostToDevice );

  cudaMalloc( (void**)&(src.d_point_freqs), num_comps*sizeof(float) );
  cudaMemcpy( src.d_point_freqs, src.point_freqs, num_comps*sizeof(float), cudaMemcpyHostToDevice );

}

//

extern "C" void copy_XYZ_to_GPU(float *d_X_diff, float *d_Y_diff, float *d_Z_diff,
                                float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                                const int num_baselines) {
  // printf("What the fuck chico man\n");

  cudaError_t err;
  err = cudaMalloc( (void**)&d_X_diff, num_baselines*sizeof(float) );
  if( cudaSuccess != err ) { \
    fprintf( stderr, "CUDA Error in Malloc\n");	\
  }

  if( cudaSuccess != err ) { \
    fprintf( stderr, "CUDA Error in Memcpy\n");	\
  }
  err = cudaMemcpy( d_X_diff, X_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );

  cudaMalloc( (void**)&d_Y_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Y_diff, Y_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_Z_diff, num_baselines*sizeof(float) );
  cudaMemcpy( d_Z_diff, Z_diff_metres, num_baselines*sizeof(float), cudaMemcpyHostToDevice );

}
