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
  // float mod_baseline = (float)floorf(iBaseline / 8128.0);
  // int mod_baseline = 0;

  u = (d_sha0*d_X_diff[mod_baseline]) + (d_cha0*d_Y_diff[mod_baseline]);
  v = (d_sdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_cdec0*d_Z_diff[mod_baseline]) - (d_sdec0*d_cha0*d_X_diff[mod_baseline]);
  w = (d_cdec0*d_cha0*d_X_diff[mod_baseline]) - (d_cdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_sdec0*d_Z_diff[mod_baseline]);

  // printf("%d\n", iBaseline);

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
  // float d_sdec0 = d_angles_array[0];
  // float d_cdec0 = d_angles_array[1];
  // float d_sha0 = d_sha0s[iBaseline];
  // float d_cha0 = d_cha0s[iBaseline];

  //TODO do the sin/cos outside of the GPU kernel?
  float d_sdec0 = sinf(d_decs[iComponent]);
  float d_cdec0 = cosf(d_decs[iComponent]);
  float d_sha0 = sinf(d_lsts[iBaseline] - d_ras[iComponent]);
  float d_cha0 = cosf(d_lsts[iBaseline] - d_ras[iComponent]);

  int mod_baseline = iBaseline - 8128*floorf(iBaseline / 8128.0);
  // float mod_baseline = (float)floorf(iBaseline / 8128.0);
  // int mod_baseline = 0;

  u_s = (d_sha0*d_X_diff[mod_baseline]) + (d_cha0*d_Y_diff[mod_baseline]);
  v_s = (d_sdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_cdec0*d_Z_diff[mod_baseline]) - (d_sdec0*d_cha0*d_X_diff[mod_baseline]);
  w_s = (d_cdec0*d_cha0*d_X_diff[mod_baseline]) - (d_cdec0*d_sha0*d_Y_diff[mod_baseline]) + (d_sdec0*d_Z_diff[mod_baseline]);

  // printf("%d\n", iBaseline);

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

  // calc_lmn(d_angles_array, d_ras[iComponent], d_decs[iComponent],
  //                          l, m, n);

  * u = d_u_metres[iBaseline] / d_wavelength;
  * v = d_v_metres[iBaseline] / d_wavelength;
  * w = d_w_metres[iBaseline] / d_wavelength;

  float cat_wavelength = VELC / d_freqs[iComponent];
  // float cat_flux = d_fluxes[iComponent];
  // float wave_ratio = cat_wavelength / d_wavelength;
  // float SI_extrap = powf(wave_ratio,DEFAULT_SI);
  // printf("%f %f %f %f %f\n",flux,d_wavelength, extrap_wavelength, wave_ratio,SI_extrap);
  //
  * extrap_flux = d_fluxes[iComponent] * powf(cat_wavelength / d_wavelength,DEFAULT_SI);
  // * extrap_flux = flux * SI_extrap;

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

  // register float i = atomicAdd(&d_sum_visi_real[iBaseline],v);
  // register float j = atomicAdd(&d_sum_visi_imag[iBaseline],w);

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

}

__global__ void calc_visi_shapelets2(float *d_S2_ras, float *d_S2_decs, float *d_S2_fluxes, float *d_S2_freqs,
      float *d_u_metres, float *d_v_metres, float *d_w_metres, float *d_wavelengths,
      float *d_u_s_metres, float *d_v_s_metres, float *d_w_s_metres,
      float *d_sum_visi_real, float *d_sum_visi_imag,
      float *d_angles_array,float *d_S2_pas, float *d_S2_majors, float *d_S2_minors,
      float *d_S2_n1s, float *d_S2_n2s, float *d_S2_coeffs,float *d_S2_param_indexes,
      float *d_S2_ls, float *d_S2_ms, float *d_S2_ns,
      float *d_sbf2,
      const int num_baselines, const int num_visis){


  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);


  float u, v, w;
  float l, m, n;
  float extrap_flux;



  int param_index = d_S2_param_indexes[iComponent];

  extrap_uvw_flux(d_angles_array,
                  d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
                  d_S2_freqs, d_S2_fluxes,
                  iComponent, iBaseline, param_index,
                  &u, &v, &w, &extrap_flux);
  //
  //
  float pa = d_S2_pas[param_index];
  float sinpa = sin(pa);
  float cospa = cos(pa);

  float d_wavelength = d_wavelengths[iBaseline];
  //
  // //Again, need negs of what's in the RTS - dunno why BAD
  float u_s = -d_u_s_metres[param_index*num_visis + iBaseline] / d_wavelength;
  float v_s = -d_v_s_metres[param_index*num_visis + iBaseline] / d_wavelength;
  //
  float x = -( cospa*v_s + sinpa*u_s); // major axis
  float y =  (-sinpa*v_s + cospa*u_s); // minor axis

  //Scales the FWHM to std to match basis functions, and account for the
  //basis functions being stored with beta = 1.0
  float const_x = (d_S2_majors[param_index]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
  float const_y = (d_S2_minors[param_index]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

  // I^(n1+n2) = Ipow_lookup[(n1+n2) % 4]
  cuFloatComplex Ipow_lookup[] = { make_cuFloatComplex(  1.0,  0.0 ),
                                   make_cuFloatComplex(  0.0,  1.0 ),
                                   make_cuFloatComplex( -1.0,  0.0 ),
                                   make_cuFloatComplex(  0.0, -1.0 ) };
  //
  float xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

  // find the indices in the basis functions for u*beta_u and v*beta_v

  float xpos = x*const_x + sbf_c2;
  float ypos = y*const_y + sbf_c2;

  int xindex = (int)floor(xpos);
  int yindex = (int)floor(ypos);
  //
  // loop over shapelet coefficients, building up the intensity model, if this baseline is in range

  // if ( xindex >= 0 && yindex >= 0 && xindex+1 < sbf_L2 && yindex+1 < sbf_L2 ) continue;

  int n1 = (int)d_S2_n1s[iComponent];
  int n2 = (int)d_S2_n2s[iComponent];

  // if ( n1 < 0 || n2 < 0 || n1 >= sbf_N2 || n2 >= sbf_N2 ) continue;

  f_hat = d_S2_coeffs[iComponent];
  //
  sbf_n = &d_sbf2[n1*sbf_L2];
  xlow  = sbf_n[xindex];
  xhigh = sbf_n[xindex+1];
  u_value = xlow + (xhigh-xlow)*(xpos-xindex);

  sbf_n = &d_sbf2[n2*sbf_L2];
  ylow  = sbf_n[yindex];
  yhigh = sbf_n[yindex+1];
  v_value = ylow + (yhigh-ylow)*(ypos-yindex);

  // accumulate the intensity model for baseline pair (u,v)
  cuFloatComplex V_envelop = make_cuFloatComplex( 0.0, 0.0 );
  V_envelop = V_envelop + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;
  //

  l = d_S2_ls[param_index];
  m = d_S2_ms[param_index];
  n = d_S2_ns[param_index];

  cuFloatComplex visi;
  //Not sure why, but get exact match with oskar sims and correct location
  //on sky through wsclean without negative infront on 2pi

  float temp = 2*M_PI*( u*l + v*m + w*(n-1) );
  sincosf(temp, &(visi.y), &(visi.x));

  visi = cuCmulf(visi, V_envelop);

  register float i = atomicAdd(&d_sum_visi_real[iBaseline],visi.x * extrap_flux);
  register float j = atomicAdd(&d_sum_visi_imag[iBaseline],visi.y * extrap_flux);
  //
  // register float i = atomicAdd(&d_sum_visi_real[iBaseline],extrap_flux);
  // register float j = atomicAdd(&d_sum_visi_imag[iBaseline],param_index);

}

extern "C" void Atomic_time_step(float *X_diff_metres, float *Y_diff_metres, float *Z_diff_metres,
                    catsource_t catsource, float *angles_array,
                    const int num_baselines, const int num_visis,
                    visibility_set_t visibility_set,
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
  cudaMemcpy( d_sha0s, visibility_set.sha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_cha0s, num_visis*sizeof(float) );
  cudaMemcpy( d_cha0s, visibility_set.cha0s, num_visis*sizeof(float), cudaMemcpyHostToDevice );
  cudaMalloc( (void**)&d_wavelengths, num_visis*sizeof(float) );
  cudaMemcpy( d_wavelengths, visibility_set.wavelengths, num_visis*sizeof(float), cudaMemcpyHostToDevice );

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
  int num_S2s = catsource.n_S2s;

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

    cudaFree( d_point_ras);
    cudaFree( d_point_decs);
    cudaFree( d_point_fluxes );
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

    cudaFree( d_gauss_ras);
    cudaFree( d_gauss_decs);
    cudaFree( d_gauss_fluxes );
    cudaFree( d_gauss_pas);
    cudaFree( d_gauss_majors);
    cudaFree( d_gauss_minors );

  }//if gauss sources

  if (num_S2s > 0) {
    float *d_S2_ras=NULL;
    float *d_S2_decs=NULL;
    float *d_S2_fluxes=NULL;
    float *d_S2_freqs=NULL;
    float *d_S2_pas=NULL;
    float *d_S2_majors=NULL;
    float *d_S2_minors=NULL;

    float *d_S2_coeffs=NULL;
    float *d_S2_n1s=NULL;
    float *d_S2_n2s=NULL;
    float *d_S2_param_indexes=NULL;

    float *d_sbf2=NULL;
    float *d_lsts=NULL;

    float *d_u_s_metres = NULL;
    float *d_v_s_metres = NULL;
    float *d_w_s_metres = NULL;

    float *d_S2_ls=NULL;
    float *d_S2_ms=NULL;
    float *d_S2_ns=NULL;

    //Who likes cudaMalloc cudaMallocs? We like cudaMalloc cudaMallocs
    cudaMalloc( (void**)&(d_S2_ras), num_S2s*sizeof(float) );
    cudaMemcpy( d_S2_ras, catsource.S2_ras, num_S2s*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_decs), num_S2s*sizeof(float) );
    cudaMemcpy( d_S2_decs, catsource.S2_decs, num_S2s*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_fluxes), num_S2s*sizeof(float) );
    cudaMemcpy( d_S2_fluxes, catsource.S2_fluxes, num_S2s*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_pas), num_S2s*sizeof(float) );
    cudaMemcpy( d_S2_pas, catsource.S2_pas, num_S2s*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_majors), num_S2s*sizeof(float) );
    cudaMemcpy( d_S2_majors, catsource.S2_majors, num_S2s*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_minors), num_S2s*sizeof(float) );
    cudaMemcpy( d_S2_minors, catsource.S2_minors, num_S2s*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_freqs), num_S2s*sizeof(float) );
    cudaMemcpy( d_S2_freqs, catsource.S2_freqs, num_S2s*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_coeffs), catsource.n_S2_coeffs*sizeof(float) );
    cudaMemcpy( d_S2_coeffs, catsource.S2_coeffs, catsource.n_S2_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_n1s), catsource.n_S2_coeffs*sizeof(float) );
    cudaMemcpy( d_S2_n1s, catsource.S2_n1s, catsource.n_S2_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_n2s), catsource.n_S2_coeffs*sizeof(float) );
    cudaMemcpy( d_S2_n2s, catsource.S2_n2s, catsource.n_S2_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_param_indexes), catsource.n_S2_coeffs*sizeof(float) );
    cudaMemcpy( d_S2_param_indexes, catsource.S2_param_indexes, catsource.n_S2_coeffs*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_sbf2), sbf_N2*sbf_L2*sizeof(float) );
    cudaMemcpy( d_sbf2, sbf2, sbf_N2*sbf_L2*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&(d_S2_ls), num_S2s*sizeof(float) );
    cudaMalloc( (void**)&(d_S2_ms), num_S2s*sizeof(float) );
    cudaMalloc( (void**)&(d_S2_ns), num_S2s*sizeof(float) );

    cudaMalloc( (void**)&(d_lsts), num_visis*sizeof(float) );
    cudaMemcpy( d_lsts, visibility_set.lsts, num_visis*sizeof(float), cudaMemcpyHostToDevice );

    cudaMalloc( (void**)&d_u_s_metres, num_S2s*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_v_s_metres, num_S2s*num_visis*sizeof(float) );
    cudaMalloc( (void**)&d_w_s_metres, num_S2s*num_visis*sizeof(float) );


    if (num_S2s == 1) {
      threads.x = 128;
      threads.y = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;

      calc_uvw_shapelet<<< grid , threads, 0 >>>(d_X_diff, d_Y_diff, d_Z_diff,
            d_u_s_metres, d_v_s_metres, d_w_s_metres,
            d_lsts, d_S2_ras, d_S2_decs,
            num_baselines, num_visis);

    }
    else {
      threads.x = 64;
      threads.y = 2;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)num_S2s) / ((float)threads.y) );

      calc_uvw_shapelet<<< grid , threads, 0 >>>(d_X_diff, d_Y_diff, d_Z_diff,
            d_u_s_metres, d_v_s_metres, d_w_s_metres,
            d_lsts, d_S2_ras, d_S2_decs,
            num_baselines, num_visis);

    }


    threads.x = 64;
    threads.y = 1;
    grid.x = (int)ceil( ((float)num_S2s / (float)threads.x) );
    grid.y = 1;

    calc_lmn<<< grid , threads, 0 >>>(d_angles_array, d_S2_ras, d_S2_decs,
                             d_S2_ls, d_S2_ms, d_S2_ns);

    //TODO need to put some kind of check inside calc_visi_point to skip a grid
    //thread point if it's outside the realms of sanity
    if (catsource.n_S2_coeffs == 1) {
      threads.x = 64;
      threads.y = 1;
      threads.z = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = 1;
      grid.z = 1;

      calc_visi_shapelets2<<< grid , threads, 0 >>>(d_S2_ras,
              d_S2_decs, d_S2_fluxes, d_S2_freqs,
              d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
              d_u_s_metres, d_v_s_metres, d_w_s_metres,
              d_sum_visi_real, d_sum_visi_imag,
              d_angles_array, d_S2_pas, d_S2_majors, d_S2_minors,
              d_S2_n1s, d_S2_n2s, d_S2_coeffs, d_S2_param_indexes,
              d_S2_ls, d_S2_ms, d_S2_ns,
              d_sbf2,
              num_baselines, num_visis);
    }

    else {

      threads.x = 64;
      threads.y = 2;
      threads.z = 1;
      grid.x = (int)ceil( (float)num_visis / (float)threads.x );
      grid.y = (int)ceil( ((float)catsource.n_S2_coeffs) / ((float)threads.y) );
      grid.z = 1;

      calc_visi_shapelets2<<< grid , threads, 0 >>>(d_S2_ras,
              d_S2_decs, d_S2_fluxes, d_S2_freqs,
              d_u_metres, d_v_metres, d_w_metres, d_wavelengths,
              d_u_s_metres, d_v_s_metres, d_w_s_metres,
              d_sum_visi_real, d_sum_visi_imag,
              d_angles_array, d_S2_pas, d_S2_majors, d_S2_minors,
              d_S2_n1s, d_S2_n2s, d_S2_coeffs, d_S2_param_indexes,
              d_S2_ls, d_S2_ms, d_S2_ns,
              d_sbf2,
              num_baselines, num_visis);
    }

    cudaFree( d_S2_ras);
    cudaFree( d_S2_decs);
    cudaFree( d_S2_fluxes );
    cudaFree( d_S2_pas);
    cudaFree( d_S2_majors);
    cudaFree( d_S2_minors );

    cudaFree(d_S2_coeffs);
    cudaFree(d_S2_n1s);
    cudaFree(d_S2_n2s);
    cudaFree(d_S2_param_indexes);

    cudaFree(d_sbf2);
    cudaFree(d_lsts);

    cudaFree(d_u_s_metres);
    cudaFree(d_v_s_metres);
    cudaFree(d_w_s_metres);

  }//if S2s



  //Get the results into host memory
  cudaMemcpy(visibility_set.sum_visi_real,d_sum_visi_real,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set.sum_visi_imag,d_sum_visi_imag,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  cudaMemcpy(visibility_set.us_metres,d_u_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set.vs_metres,d_v_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(visibility_set.ws_metres,d_w_metres,num_visis*sizeof(float),cudaMemcpyDeviceToHost);

  //Free up the GPU memory
  cudaFree( d_sum_visi_real);
  cudaFree( d_sum_visi_imag);
  cudaFree( d_u_metres );
  cudaFree( d_v_metres );
  cudaFree( d_w_metres );

  cudaFree( d_X_diff);
  cudaFree( d_Y_diff);
  cudaFree( d_Z_diff );

  cudaFree( d_cha0s);
  cudaFree( d_sha0s);

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
