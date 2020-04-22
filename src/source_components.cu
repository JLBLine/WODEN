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

__device__ void extrap_flux(float *d_wavelengths, float *d_freqs,
           float *d_fluxes, int iComponent, int iBaseline,
           float * extrap_flux){

  float d_wavelength = d_wavelengths[iBaseline];
  float cat_wavelength = VELC / d_freqs[iComponent];
  * extrap_flux = d_fluxes[iComponent] * powf(cat_wavelength / d_wavelength,DEFAULT_SI);
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

__global__ void kern_calc_visi_point(float *d_point_ras,
           float *d_point_decs, float *d_point_fluxes, float *d_point_freqs,
           float *d_us, float *d_vs, float *d_ws,
           float *d_sum_visi_real, float *d_sum_visi_imag,
           float *d_angles_array, float *d_wavelengths,
           float *d_ls, float *d_ms, float *d_ns,
           int num_points, int num_baselines, int num_freqs, int num_visis,
           float *d_beam_reals, float *d_beam_imags, int beamtype) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_points) {

    float point_flux;
    extrap_flux(d_wavelengths, d_point_freqs, d_point_fluxes,
                  iComponent, iBaseline, &point_flux);

    cuFloatComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws,
                           d_ls, d_ms, d_ns,
                           iBaseline, iComponent);

    float beam_real;
    float beam_imag;

    int beam_ind = 0.0;
    int time_ind = 0.0;
    int freq_ind = 0.0;

    if (beamtype == GAUSS_BEAM) {
      //Do some epic indexing to work out which beam value
      time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
      freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
      beam_ind = num_freqs*time_ind*num_points + (num_points*freq_ind) + iComponent;

      beam_real = d_beam_reals[beam_ind];
      beam_imag = d_beam_imags[beam_ind];

      // printf("%d %d %d %d %d %d\n",iBaseline,num_baselines,num_freqs,time_ind,freq_ind,beam_ind);

    }
    else {
      beam_real = 1.0;
      beam_imag = 0.0;
    }

    cuFloatComplex beam_complex = make_cuFloatComplex(beam_real,beam_imag);

    visi = cuCmulf(visi, beam_complex);

    atomicAdd(&d_sum_visi_real[iBaseline],visi.x*point_flux);
    atomicAdd(&d_sum_visi_imag[iBaseline],visi.y*point_flux);

    // atomicAdd(&d_sum_visi_real[iBaseline],beam_real);
    // atomicAdd(&d_sum_visi_imag[iBaseline],beam_imag);
  }
}

__global__ void kern_calc_visi_gaussian(float *d_gauss_ras,
           float *d_gauss_decs, float *d_gauss_fluxes, float *d_gauss_freqs,
           float *d_us, float *d_vs, float *d_ws,
           float *d_sum_visi_real, float *d_sum_visi_imag,
           float *d_angles_array, float *d_wavelengths,
           float *d_ls, float *d_ms, float *d_ns,
           float *d_gauss_pas, float *d_gauss_majors, float *d_gauss_minors,
           int num_gauss, int num_baselines, int num_freqs, int num_visis,
           float *d_beam_reals, float *d_beam_imags, int beamtype) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_gauss) {

    float gauss_flux;
    extrap_flux(d_wavelengths, d_gauss_freqs, d_gauss_fluxes,
                  iComponent, iBaseline, &gauss_flux);


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

    float beam_real;
    float beam_imag;

    int beam_ind = 0.0;
    int time_ind = 0.0;
    int freq_ind = 0.0;

    if (beamtype == GAUSS_BEAM) {
      //Do some epic indexing to work out which beam value
      time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
      freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
      beam_ind = num_freqs*time_ind*num_gauss + (num_gauss*freq_ind) + iComponent;

      beam_real = d_beam_reals[beam_ind];
      beam_imag = d_beam_imags[beam_ind];

    }
    else {
      beam_real = 1.0;
      beam_imag = 0.0;
    }

    cuFloatComplex beam_complex = make_cuFloatComplex(beam_real,beam_imag);

    visi = cuCmulf(visi, beam_complex);

    atomicAdd(&d_sum_visi_real[iBaseline],visi.x*gauss_flux);
    atomicAdd(&d_sum_visi_imag[iBaseline],visi.y*gauss_flux);

  }
}

__global__ void kern_calc_visi_shapelets(float *d_shape_ras,
      float *d_shape_decs, float *d_shape_fluxes, float *d_shape_freqs,
      float *d_us, float *d_vs, float *d_ws,
      float *d_wavelengths,
      float *d_u_s_metres, float *d_v_s_metres, float *d_w_s_metres,
      float *d_sum_visi_real, float *d_sum_visi_imag,
      float *d_angles_array, float *d_shape_pas, float *d_shape_majors,
      float *d_shape_minors,
      float *d_shape_n1s, float *d_shape_n2s, float *d_shape_coeffs,
      float *d_shape_param_indexes,
      float *d_shape_ls, float *d_shape_ms, float *d_shape_ns,
      float *d_sbf,
      int num_shapes, int num_baselines, int num_freqs, int num_visis,
      const int num_coeffs,
      float *d_beam_reals, float *d_beam_imags, int beamtype) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);


  if (iBaseline < num_visis && iComponent < num_coeffs) {
    // printf("Made it here all g %d %d %d %d\n", iBaseline,num_visis,iComponent,num_coeffs);
    int param_index = d_shape_param_indexes[iComponent];

    float shape_flux;
    //Use param index below and not iComponent as there
    extrap_flux(d_wavelengths, d_shape_freqs,
               d_shape_fluxes, param_index, iBaseline,
               &shape_flux);

    cuFloatComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws,
                          d_shape_ls, d_shape_ms, d_shape_ns,
                          iBaseline, param_index);

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

    float xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat, *sbf_n;

    // find the indices in the basis functions for u*beta_u and v*beta_v

    float xpos = x*const_x + sbf_c;
    float ypos = y*const_y + sbf_c;

    int xindex = (int)floor(xpos);
    int yindex = (int)floor(ypos);
    //
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

    visi = cuCmulf(visi, V_envelop);

    float beam_real;
    float beam_imag;

    int beam_ind = 0.0;
    int time_ind = 0.0;
    int freq_ind = 0.0;

    if (beamtype == GAUSS_BEAM) {
      //Do some epic indexing to work out which beam value
      time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
      freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
      beam_ind = num_freqs*time_ind*num_shapes + (num_shapes*freq_ind) + param_index;

      beam_real = d_beam_reals[beam_ind];
      beam_imag = d_beam_imags[beam_ind];

    }
    else {
      beam_real = 1.0;
      beam_imag = 0.0;
    }

    cuFloatComplex beam_complex = make_cuFloatComplex(beam_real,beam_imag);

    visi = cuCmulf(visi, beam_complex);

    atomicAdd(&d_sum_visi_real[iBaseline],visi.x * shape_flux);
    atomicAdd(&d_sum_visi_imag[iBaseline],visi.y * shape_flux);

  }

}
