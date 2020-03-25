#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "cudacomplex.h"
#include "fundamental_coords.h"
#include "constants.h"


__device__ void extrap_flux(float *d_wavelengths,
                float *d_freqs, float *d_fluxes,
                int iComponent, int iBaseline,
                float * extrap_flux ){

  float d_wavelength = d_wavelengths[iBaseline];
  float cat_wavelength = VELC / d_freqs[iComponent];
  * extrap_flux = d_fluxes[iComponent] * powf(cat_wavelength / d_wavelength,DEFAULT_SI);

}

__device__ cuFloatComplex calc_point_source_visi(float *d_u, float *d_v, float *d_w,
                                       float *d_ls, float *d_ms, float *d_ns,
                                       const int iBaseline, const int iComponent){

  float u, v, w;
  float l, m, n;

  u = d_u[iBaseline];
  v = d_v[iBaseline];
  w = d_w[iBaseline];

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



__global__ void calc_visi_point(float *d_point_ras, float *d_point_decs, float *d_point_fluxes, float *d_point_freqs,
      float *d_u_metres, float *d_v_metres, float *d_w_metres,
      float *d_u, float *d_v, float *d_w,
      float *d_sum_visi_real, float *d_sum_visi_imag,
      float *d_angles_array, float *d_wavelengths,
      float *d_ls, float *d_ms, float *d_ns,
      int num_points, int num_visis) {
  //
  float point_flux;

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis || iComponent < num_visis) {

    extrap_flux(d_wavelengths, d_point_freqs, d_point_fluxes,
                  iComponent, iBaseline, &point_flux);


    cuFloatComplex visi;
    visi = calc_point_source_visi(d_u, d_v, d_w,
                           d_ls, d_ms, d_ns,
                           iBaseline, iComponent);


    register float i = atomicAdd(&d_sum_visi_real[iBaseline],visi.x*point_flux);
    register float j = atomicAdd(&d_sum_visi_imag[iBaseline],visi.y*point_flux);

    //HERE
    // register float i = atomicAdd(&d_sum_visi_real[iBaseline],u);
    // register float j = atomicAdd(&d_sum_visi_imag[iBaseline],v);
  }
}
