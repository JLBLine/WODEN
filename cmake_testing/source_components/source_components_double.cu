#include <stdio.h>
#include <stdlib.h>
#include <cuComplex.h>
#include <complex.h>
#include <math.h>
#include "cudacomplex.h"
// #include "fundamental_coords.h"
// #include "constants.h"
// #include "shapelet_basis.h"
// #include "source_components.h"
#include "cudacheck.h"
// #include "woden_struct_defs.h"
// #include "primary_beam_cuda.h"
// #include "FEE_primary_beam_cuda.h"

__device__  cuDoubleComplex calc_measurement_equation(double *d_us,
           double *d_vs, double *d_ws, double *d_ls, double *d_ms, double *d_ns,
           const int iBaseline, const int iComponent){

  double u, v, w;
  double l, m, n;

  u = d_us[iBaseline];
  v = d_vs[iBaseline];
  w = d_ws[iBaseline];

  l = d_ls[iComponent];
  m = d_ms[iComponent];
  n = d_ns[iComponent];

  cuDoubleComplex visi;

  //Not sure why, but get match with OSKAR/RTS sims, and correct location
  //on sky through WSClean, without negative infront on 2pi
  double temp = 2*M_PI*( u*l + v*m + w*(n-1) );
  // sincosf(temp, &(visi.y), &(visi.x));

  visi.y = sin(temp);
  visi.x = cos(temp);

  return visi;
}



__global__ void kern_calc_measurement_equation(int num_components, int num_baselines,
                  double *d_us, double *d_vs, double *d_ws,
                  double *d_ls, double *d_ms, double *d_ns,
                  cuDoubleComplex *d_visis) {

  // Start by computing which baseline we're going to do
  const int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);
  // if(iBaseline < num_visis && iComponent < num_points) {
  if(iComponent < num_components && iBaseline < num_baselines) {

    cuDoubleComplex visi;
    visi = calc_measurement_equation(d_us, d_vs, d_ws, d_ls, d_ms, d_ns,
                                     iBaseline, iComponent);

    int visi_ind = num_components*iBaseline + iComponent;
    d_visis[visi_ind] = visi;

  }
}

extern "C" void test_kern_calc_measurement_equation(int num_components,
                  int num_baselines,
                  double *us, double *vs, double *ws,
                  double *ls, double *ms, double *ns,
                  double _Complex *visis){

  double *d_us = NULL;
  double *d_vs = NULL;
  double *d_ws = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_us, num_baselines*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_vs, num_baselines*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ws, num_baselines*sizeof(double) ));
  cudaErrorCheckCall( cudaMemcpy(d_us, us, num_baselines*sizeof(double), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_vs, vs, num_baselines*sizeof(double), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ws, ws, num_baselines*sizeof(double), cudaMemcpyHostToDevice ));

  double *d_ls = NULL;
  double *d_ms = NULL;
  double *d_ns = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ls, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ms, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMalloc( (void**)&d_ns, num_components*sizeof(double) ));
  cudaErrorCheckCall( cudaMemcpy(d_ls, ls, num_components*sizeof(double), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ms, ms, num_components*sizeof(double), cudaMemcpyHostToDevice ));
  cudaErrorCheckCall( cudaMemcpy(d_ns, ns, num_components*sizeof(double), cudaMemcpyHostToDevice ));

  double _Complex *d_visis = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_visis, num_baselines*num_components*sizeof(double _Complex) ));

  dim3 grid, threads;

  threads.x = 16;
  threads.y = 16;
  grid.x = (int)ceil( (double)num_baselines / (double)threads.x );
  grid.y = (int)ceil( (double)num_components / (double)threads.y );

  cudaErrorCheckKernel("kern_calc_measurement_equation",
                      kern_calc_measurement_equation, grid, threads,
                      num_components, num_baselines,
                      d_us, d_vs, d_ws,
                      d_ls, d_ms, d_ns,
                      (cuDoubleComplex*)d_visis );

  cudaErrorCheckCall( cudaMemcpy(visis, (double _Complex*)d_visis, num_components*num_baselines*sizeof(double _Complex),cudaMemcpyDeviceToHost ));

  cudaErrorCheckCall( cudaFree( d_us ) );
  cudaErrorCheckCall( cudaFree( d_vs ) );
  cudaErrorCheckCall( cudaFree( d_ws ) );
  cudaErrorCheckCall( cudaFree( d_ls ) );
  cudaErrorCheckCall( cudaFree( d_ms ) );
  cudaErrorCheckCall( cudaFree( d_ns ) );
  cudaErrorCheckCall( cudaFree(d_visis ) );

}
