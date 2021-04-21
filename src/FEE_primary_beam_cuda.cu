#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuComplex.h>
#include <complex.h>
#include <assert.h>
#include <stdio.h>

#include "cudacomplex.h"
#include "cudacheck.h"
#include "FEE_primary_beam_cuda.h"

/**
Factorials used in calculating legendre polynomials
*/
__constant__ float ffactorials[] = {
  1.0,
  1.0,
  2.0,
  6.0,
  24.0,
  120.0,
  720.0,
  5040.0,
  40320.0,
  362880.0,
  3628800.0,
  39916800.0,
  479001600.0,
  6227020800.0,
  87178291200.0,
  1.307674368e+12,
  2.0922789888e+13,
  3.55687428096e+14,
  6.40237370573e+15,
  1.21645100409e+17,
  2.43290200818e+18,
  5.10909421717e+19,
  1.12400072778e+21,
  2.58520167389e+22,
  6.20448401733e+23,
  1.55112100433e+25,
  4.03291461127e+26,
  1.08888694504e+28,
  3.04888344612e+29,
  8.84176199374e+30,
  2.65252859812e+32,
  8.22283865418e+33,
  2.63130836934e+35,
  8.68331761881e+36,
  2.9523279904e+38,
  1.03331479664e+40,
  3.7199332679e+41,
  1.37637530912e+43,
  5.23022617467e+44,
  2.03978820812e+46,
  8.15915283248e+47,
  3.34525266132e+49,
  1.40500611775e+51,
  6.04152630634e+52,
  2.65827157479e+54,
  1.19622220865e+56,
  5.50262215981e+57,
  2.58623241511e+59,
  1.24139155925e+61,
  6.08281864034e+62,
  3.04140932017e+64,
  1.55111875329e+66,
  8.06581751709e+67,
  4.27488328406e+69,
  2.30843697339e+71,
  1.26964033537e+73,
  7.10998587805e+74,
  4.05269195049e+76,
  2.35056133128e+78,
  1.38683118546e+80,
  8.32098711274e+81,
  5.07580213877e+83,
  3.14699732604e+85,
  1.9826083154e+87,
  1.26886932186e+89,
  8.24765059208e+90,
  5.44344939077e+92,
  3.64711109182e+94,
  2.48003554244e+96,
  1.71122452428e+98,
  1.197857167e+100,
  8.50478588568e+101,
  6.12344583769e+103,
  4.47011546151e+105,
  3.30788544152e+107,
  2.48091408114e+109,
  1.88549470167e+111,
  1.45183092028e+113,
  1.13242811782e+115,
  8.94618213078e+116
};

extern "C" void copy_FEE_primary_beam_to_GPU(RTS_MWA_FEE_beam_t *FEE_beam){

  const int n_pols=2; // instrumental pols
  const int nMN = FEE_beam->nMN; // max_length HDFBeamInit
  int arrSize = n_pols * nMN;

  //Currently only calculate for one station - expand to more in the future?
  int nStations = 1;

  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_M, arrSize*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_N, arrSize*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_Q1, arrSize*sizeof(float _Complex) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_Q2, arrSize*sizeof(float _Complex) ) );

  float *h_params;

  h_params = (float *)malloc(arrSize*sizeof(float) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (float)(FEE_beam->M[i][j]);
    }
  }

  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_M, h_params,
                      arrSize*sizeof(float), cudaMemcpyHostToDevice ) );

  for(unsigned int i=0; i<n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (float)(FEE_beam->N[i][j]);
    }
  }

  cudaErrorCheckCall( cudaMemcpy( FEE_beam->d_N, h_params,
                      arrSize*nStations*sizeof(float), cudaMemcpyHostToDevice ) );
  // cudaFreeHost( h_params );
  free(h_params);

  float _Complex *h_Qdata;
  h_Qdata = (float _Complex *)malloc(arrSize*sizeof(float _Complex) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = FEE_beam->Q1[i][j];
    }
  }

  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_Q1, h_Qdata,
                      arrSize*nStations*sizeof(float _Complex), cudaMemcpyHostToDevice ) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = FEE_beam->Q2[i][j];
    }
  }

  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_Q2, h_Qdata,
                      arrSize*nStations*sizeof(float _Complex) , cudaMemcpyHostToDevice ) );
  free( h_Qdata );

}

// // GPU version of legendre polynomial function. This version only returns the value
// // at a single point as required in HDF beam calculation
// // See legendre_polynomial.c for source and info
__device__ void RTS_CUDA_pm_polynomial_value_singlef(int n, int m, float x, float *values  ){

  float fact;
  int i;
  int j;
  int k;
  int mm=1;

  //v = ( float * ) malloc ( mm*(n+1) * sizeof ( float ) );;

  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      values[i+j*mm] = 0.0;
    }
  }
  /*
    J = M is the first nonzero function.
  */
  if ( m <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      values[i+m*mm] = 1.0;
    }

    fact = 1.0;
    for ( k = 0; k < m; k++ )
    {
      for ( i = 0; i < mm; i++ )
      {
        values[i+m*mm] = - values[i+m*mm] * fact * sqrt ( 1.0 - x * x );
      }
      fact = fact + 2.0;
    }
  }
  /*
    J = M + 1 is the second nonzero function.
  */
  if ( m + 1 <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      values[i+(m+1)*mm] = x * ( float ) ( 2 * m + 1 ) * values[i+m*mm];
    }
  }
  /*
    Now we use a three term recurrence.
  */
  for ( j = m + 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      values[i+j*mm] = ( ( float ) ( 2 * j     - 1 ) * x * values[i+(j-1)*mm]
                  + ( float ) (   - j - m + 1 ) *        values[i+(j-2)*mm] )
                  / ( float ) (     j - m     );
    }
  }
}

__global__ void RTS_P1SINfKernel(float *d_theta, cuFloatComplex *rts_P_sin,
           cuFloatComplex *rts_p1, int nmax, int num_coords){

  int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  float theta = d_theta[iCoord];

  //RTS style variables
  int i = threadIdx.y + (blockDim.y*blockIdx.y);
  int n = floor(sqrt((float)(i+1)));
  int index = i - ((n-1)*(n-1) + 2*(n-1));
  int pm_index;
  float u = cos(theta);
  float sin_th = sin(theta);
  float delu = 1.0e-6;

  int nMN = nmax*nmax + 2*nmax;

  if(iCoord < num_coords && i < nMN) {

    if(index < n){
      pm_index = n-index;
    } else {
      pm_index = index-n;
    }

    float p;
    float Pm1;
    float Pm_sin;

    float pm_vals[100];

    RTS_CUDA_pm_polynomial_value_singlef(n, pm_index, u, pm_vals);
    p = pm_vals[n];
    // Diverging??
    RTS_CUDA_pm_polynomial_value_singlef(n, pm_index+1, u, pm_vals);
    if(pm_index==n){
      Pm1 = 0.0;
    } else {
      Pm1 = pm_vals[n];
    }
    Pm_sin = 0.0;

    if(u==1.0){
      float pm_val;
      RTS_CUDA_pm_polynomial_value_singlef(n, 0, u-delu, pm_vals);
      pm_val = pm_vals[n];
      float Pm_sin1,p0;
      RTS_CUDA_pm_polynomial_value_singlef(0, 0, u, pm_vals);
      p0 = pm_vals[0];
      Pm_sin1 = -(p0 - pm_val) / delu;
      if(pm_index==1) Pm_sin = Pm_sin1;
    } else if (u == -1.0){
      float pm_val;
      RTS_CUDA_pm_polynomial_value_singlef(n, 0, u-delu, pm_vals);
      pm_val = pm_vals[n];
      float Pm_sin1,p0;
      RTS_CUDA_pm_polynomial_value_singlef(0, 0, u, pm_vals);
      p0 = pm_vals[0];
      Pm_sin1 = -(pm_val - p0) / delu;
      if(pm_index==1) Pm_sin = Pm_sin1;
    } else {
      Pm_sin = p / sin_th;
    }

    rts_P_sin[iCoord*nMN + i] = make_cuFloatComplex(Pm_sin,0);
    rts_p1[iCoord*nMN + i] = make_cuFloatComplex(Pm1,0);

  } //end IF thread lies within number of calculations we want to do
}

__global__ void RTS_getTileGainsKernel( float *d_phi, float *d_theta, int nMN, int num_coords,
           float *pb_M, float *pb_N,
           cuFloatComplex *pb_Q1, cuFloatComplex *pb_Q2,
           cuFloatComplex *rts_P_sin, cuFloatComplex *rts_P1,
           cuFloatComplex *emn_T, cuFloatComplex *emn_P){


  int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  int i = threadIdx.y + (blockDim.y*blockIdx.y);

  int pol = blockIdx.z;
  int n_pols = gridDim.z;

  int index;
  float factor1,factor2;
  float C_MN, M_absM;

  int beam_index = i+(pol*nMN);

  if(iCoord < num_coords && i < nMN) {
    float phi = d_phi[iCoord];
    float theta = d_theta[iCoord];

    index = (int) (pb_N[beam_index] - fabs(pb_M[beam_index]));
    if(index >=80){
      printf("Maximum factorial exceeded in RTS_getTileGainsKernel (attempted %d)\n", index );
    }
    factor1 = ffactorials[index];
    index = (int) (pb_N[beam_index] + fabs(pb_M[beam_index]));
    if(index >=80){
      printf("Maximum factorial exceeded  in RTS_getTileGainsKernel (attempted %d)\n", index );
    }
    factor2 = ffactorials[index];

    C_MN = sqrt(0.5 * (2 * pb_N[beam_index] + 1) * factor1 / factor2);

    if (pb_M[beam_index] > 0) {
      if ((int)pb_M[beam_index] % 2 == 0) {
        M_absM = 1;
      }
      else {
        M_absM = -1;
      }
    }
    else {
      M_absM = 1;
    }

    cuFloatComplex phi_comp;

    phi_comp = U1polar(pb_M[beam_index] * phi) * C_MN * M_absM / sqrt(pb_N[beam_index] * (pb_N[beam_index] + 1.0));

    cuFloatComplex this_emn_T, this_emn_P;
    float T_exp, P_exp;
    float u = cos(theta);

    T_exp = (pb_N[beam_index]);
    P_exp = (pb_N[beam_index])+1.0;

    this_emn_T = make_cuComplex(cos(T_exp * M_PI/2.0),sin(T_exp * M_PI/2.0));
    this_emn_P = make_cuComplex(cos(P_exp * M_PI/2.0),sin(P_exp * M_PI/2.0));

    this_emn_T *= (rts_P_sin[iCoord*nMN + i]*(fabs(pb_M[beam_index])*pb_Q2[beam_index]*u - pb_M[beam_index] * pb_Q1[beam_index]) + pb_Q2[beam_index] * rts_P1[iCoord*nMN + i]);
    this_emn_P *= (rts_P_sin[iCoord*nMN + i]*(pb_M[beam_index]*pb_Q2[beam_index] - fabs(pb_M[beam_index])*pb_Q1[beam_index]*u) - pb_Q1[beam_index] * rts_P1[iCoord*nMN + i]);

    emn_T[iCoord*n_pols*nMN + beam_index] = cuCmulf(this_emn_T,phi_comp);
    emn_P[iCoord*n_pols*nMN + beam_index] = cuCmulf(this_emn_P,phi_comp);

  }
}

__global__ void kern_sum_emn_PT_by_M(cuFloatComplex *emn_T, cuFloatComplex *emn_P,
           float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
           float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
           float *d_m_range, float *d_M, int nMN, int nmax, int num_coords){

  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iM_index = threadIdx.y + (blockDim.y*blockIdx.y);
  const int iM_value = threadIdx.z + (blockDim.z*blockIdx.z);


  int iPol = (int)floorf((float)iM_index / (float)nMN);
  int n_pols = 2;

  int num_sum_M = 2*nmax + 1;

  if(iM_index < n_pols*nMN && iM_value < num_sum_M && iCoord < num_coords) {
    if (d_M[iM_index] == d_m_range[iM_value]) {
      cuFloatComplex emn_P_value = emn_P[n_pols*nMN*iCoord + iM_index];
      cuFloatComplex emn_T_value = emn_T[n_pols*nMN*iCoord + iM_index];

      atomicAdd(&d_emn_P_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_P_value.x);
      atomicAdd(&d_emn_P_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_P_value.y);

      atomicAdd(&d_emn_T_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_T_value.x);
      atomicAdd(&d_emn_T_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_T_value.y);

    }
  }
}

__global__ void kern_sum_emn_PT_by_M_loop(cuFloatComplex *emn_T, cuFloatComplex *emn_P,
           float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
           float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
           float *d_m_range, float *d_M, int nMN, int nmax, int num_coords){

  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  // const int iM_index = threadIdx.y + (blockDim.y*blockIdx.y);
  // const int iM_value = threadIdx.z + (blockDim.z*blockIdx.z);

  int n_pols = 2;
  int num_sum_M = 2*nmax + 1;

  // if(iM_index < n_pols*nMN && iM_value < num_sum_M && iCoord < num_coords) {
  if(iCoord < num_coords) {

    for (int iM_index = 0; iM_index < n_pols*nMN; iM_index++) {
      int iPol = (int)floorf((float)iM_index / (float)nMN);

      for (int iM_value = 0; iM_value < num_sum_M; iM_value++) {
        if (d_M[iM_index] == d_m_range[iM_value]) {
          cuFloatComplex emn_P_value = emn_P[n_pols*nMN*iCoord + iM_index];
          cuFloatComplex emn_T_value = emn_T[n_pols*nMN*iCoord + iM_index];

          // atomicAdd(&d_emn_P_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_P_value.x);
          // atomicAdd(&d_emn_P_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_P_value.y);
          //
          // atomicAdd(&d_emn_T_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_T_value.x);
          // atomicAdd(&d_emn_T_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax], emn_T_value.y);

          d_emn_P_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_P_value.x;
          d_emn_P_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_P_value.y;

          d_emn_T_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_T_value.x;
          d_emn_T_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_T_value.y;
        }
      }
    }
  }
}

__global__ void kern_map_emn(cuFloatComplex *TileGainMatrices,
                float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
                float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
                int nmax, int num_coords) {

  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iPol = threadIdx.y + (blockDim.y*blockIdx.y);
  int n_pols=2;

  int numM = 2*nmax + 1;

  if (iCoord < num_coords && iPol < n_pols){

    for(int i=0; i < numM; i++){

      int array_loc = n_pols*numM*iCoord + iPol*numM + i;

      cuFloatComplex d_emn_T = make_cuComplex(d_emn_T_sum_real[array_loc],d_emn_T_sum_imag[array_loc]);
      cuFloatComplex d_emn_P = make_cuComplex(d_emn_P_sum_real[array_loc],d_emn_P_sum_imag[array_loc]);

      TileGainMatrices[iCoord*MAX_POLS + iPol*n_pols] += d_emn_T;
      TileGainMatrices[iCoord*MAX_POLS + iPol*n_pols + 1] += d_emn_P;

    }

  }
}

__global__ void kern_apply_FEE_norm(cuFloatComplex *TileGainMatrices,
           cuFloatComplex *d_norm_fac, int num_coords ){

  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iCoord < num_coords) {

    cuFloatComplex before0 = TileGainMatrices[iCoord*MAX_POLS + 0];
    cuFloatComplex before1 = TileGainMatrices[iCoord*MAX_POLS + 1];
    cuFloatComplex before2 = TileGainMatrices[iCoord*MAX_POLS + 2];
    cuFloatComplex before3 = TileGainMatrices[iCoord*MAX_POLS + 3];

    //Again, reorder polarisations here to follow the RTS
    TileGainMatrices[iCoord*MAX_POLS + 0] = cuCdivf(before3, d_norm_fac[0]);
    TileGainMatrices[iCoord*MAX_POLS + 1] = cuCdivf(before2, d_norm_fac[1]);
    TileGainMatrices[iCoord*MAX_POLS + 2] = cuCdivf(before1, d_norm_fac[2]);
    TileGainMatrices[iCoord*MAX_POLS + 3] = cuCdivf(before0, d_norm_fac[3]);

  }
}

__global__ void kern_rotate_FEE_beam(cuFloatComplex *d_FEE_beam_gain_matrices,
                                float *d_sin_para_angs, float *d_cos_para_angs,
                                int num_coords) {
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iCoord < num_coords) {

    cuFloatComplex prerot0 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 0];
    cuFloatComplex prerot1 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 1];
    cuFloatComplex prerot2 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 2];
    cuFloatComplex prerot3 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 3];

    float sinrot = d_sin_para_angs[iCoord];
    float cosrot = d_cos_para_angs[iCoord];

    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 0] = -prerot1*sinrot + prerot0*cosrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 1] = prerot1*cosrot + prerot0*sinrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 2] = -prerot3*sinrot + prerot2*cosrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 3] = prerot3*cosrot + prerot2*sinrot;

  }
}

extern "C" void RTS_CUDA_get_TileGains(float *phi, float *theta,
           float *sin_para_angs, float *cos_para_angs,
           int num_time_steps, int num_components,
           float rotation, RTS_MWA_FEE_beam_t *primary_beam,
           float _Complex *TileGainMatrices, int scaling){

  int n_pols = 2;
  int num_coords = num_time_steps*num_components;

  for (size_t phi_ind = 0; phi_ind < num_coords; phi_ind++) {
    phi[phi_ind] = (M_PI/2.0) - phi[phi_ind];
    if(phi[phi_ind] < 0){
      phi[phi_ind] += 2.0*M_PI;
    }
  }

  float *d_phi=NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_phi, num_coords*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy(d_phi, phi, num_coords*sizeof(float),
                      cudaMemcpyHostToDevice ) );

  float *d_theta=NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_theta, num_coords*sizeof(float) ) );
  cudaErrorCheckCall( cudaMemcpy(d_theta, theta, num_coords*sizeof(float),
                      cudaMemcpyHostToDevice ) );

  int nmax = primary_beam->nmax; // Assumes all beams have same dimensions
  int nMN = primary_beam->nMN;

  int P_size = sizeof(float _Complex)*(nmax*nmax + 2*nmax)*num_coords;

  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->rts_P_sin, P_size) );
  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->rts_P1, P_size) );

  dim3 grid, threads;

  threads.y = 2;
  grid.y = (int)ceil((float)(nmax*nmax + 2*nmax) / (float)threads.y);

  threads.x = 64;
  grid.x = (int)ceil((float)num_coords / (float)threads.x);

  cudaErrorCheckKernel("RTS_P1SINfKernel",
                       RTS_P1SINfKernel, grid, threads,
                       d_theta, (cuFloatComplex*)primary_beam->rts_P_sin,
                       (cuFloatComplex*)primary_beam->rts_P1, nmax, num_coords);

  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->emn_P,
                    num_coords * nMN * n_pols * sizeof(float _Complex)) );
  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->emn_T,
                    num_coords * nMN * n_pols * sizeof(float _Complex)) );
  //
  //
  threads.x = 128;
  threads.y = 2;
  threads.z = 1;
  //
  grid.x = (int)ceil((float)num_coords / (float)threads.x);
  grid.y = (int)ceil((float)(nMN) / (float)threads.y);
  grid.z = n_pols;

  cudaErrorCheckKernel("RTS_getTileGainsKernel",
                        RTS_getTileGainsKernel, grid, threads,
                        d_phi, d_theta, nMN, num_coords,
                        primary_beam->d_M,primary_beam->d_N,
                        (cuFloatComplex*)primary_beam->d_Q1,
                        (cuFloatComplex*)primary_beam->d_Q2,
                        (cuFloatComplex*)primary_beam->rts_P_sin,
                        (cuFloatComplex*)primary_beam->rts_P1,
                        (cuFloatComplex*)primary_beam->emn_T,
                        (cuFloatComplex*)primary_beam->emn_P);

  float *d_emn_P_sum_real = NULL;
  float *d_emn_T_sum_real = NULL;
  float *d_emn_P_sum_imag = NULL;
  float *d_emn_T_sum_imag = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_P_sum_real,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_T_sum_real,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_P_sum_imag,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(float) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_T_sum_imag,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(float) ) );

  //Ensure gains are zero before summing results to them
  float *zero_array = NULL;
  zero_array = (float *)malloc(num_coords*(2*nmax + 1)*n_pols*sizeof(float) );

  for (int i = 0; i < num_coords*(2*nmax + 1)*n_pols; i++) {
    zero_array[i] = 0.0;
  }

  cudaErrorCheckCall( cudaMemcpy(d_emn_P_sum_real, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice) );
  cudaErrorCheckCall( cudaMemcpy(d_emn_P_sum_imag, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice) );
  cudaErrorCheckCall( cudaMemcpy(d_emn_T_sum_real, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice) );
  cudaErrorCheckCall( cudaMemcpy(d_emn_T_sum_imag, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice) );

  float *d_m_range = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_m_range, (2*nmax + 1)*sizeof(float)) );
  cudaErrorCheckCall( cudaMemcpy(d_m_range, primary_beam->m_range,
                      (2*nmax + 1)*sizeof(float), cudaMemcpyHostToDevice ) );

  threads.x = 16;
  threads.y = 8;
  threads.z = 8;

  grid.x = (int)ceil( (float)num_coords / threads.x);
  grid.y = (int)ceil( (2.0*(float)primary_beam->nMN) / threads.y);
  grid.z = (int)ceil( ((float)nmax*2.0 + 1.0) / threads.z);

  // threads.x = 128;
  // threads.y = 1;
  // threads.z = 1;
  //
  // grid.x = (int)ceil( (float)num_coords / threads.x);
  // grid.y = 1;
  // grid.z = 1;

  cudaErrorCheckKernel("kern_sum_emn_PT_by_M",
                        kern_sum_emn_PT_by_M, grid, threads,
                        (cuFloatComplex*)primary_beam->emn_T,
                        (cuFloatComplex*)primary_beam->emn_P,
                        d_emn_T_sum_real, d_emn_T_sum_imag,
                        d_emn_P_sum_real, d_emn_P_sum_imag,
                        d_m_range, primary_beam->d_M, primary_beam->nMN,
                        nmax, num_coords);

  threads.x = 64;
  threads.y = n_pols;
  threads.z = 1;

  grid.x = (int)ceil( (float)num_coords / threads.x);
  grid.y = grid.z = 1;

  free(zero_array);

  cudaErrorCheckKernel("kern_map_emn",
                        kern_map_emn, grid, threads,
                        (cuFloatComplex*)TileGainMatrices,
                        d_emn_T_sum_real, d_emn_T_sum_imag,
                        d_emn_P_sum_real, d_emn_P_sum_imag,
                        nmax, num_coords );

  //If we want to normalise the beam, do it here
  if (scaling == 1.0) {

    threads.x = 128;
    grid.x = (int)ceil( (float)num_coords / threads.x);

    threads.y = threads.z = 1;
    grid.y = grid.z = 1;

    float _Complex *d_norm_fac=NULL;
    cudaErrorCheckCall( cudaMalloc( (void**)&d_norm_fac,
                        MAX_POLS*sizeof(float _Complex) ) );
    cudaErrorCheckCall( cudaMemcpy(d_norm_fac, primary_beam->norm_fac,
               MAX_POLS*sizeof(float _Complex), cudaMemcpyHostToDevice ) );

    cudaErrorCheckKernel("kern_apply_FEE_norm",
                         kern_apply_FEE_norm, grid, threads,
                         (cuFloatComplex*)TileGainMatrices,
                         (cuFloatComplex*)d_norm_fac, num_coords );

    cudaErrorCheckCall( cudaFree(d_norm_fac) );

    //If rotating by parallactic angle - only ever want to do this
    //after normalistation
    if (rotation == 1.0) {

      float *d_cos_para_angs=NULL;
      float *d_sin_para_angs=NULL;

      cudaErrorCheckCall( cudaMalloc( (void**)&d_cos_para_angs,
                          num_components*num_time_steps*sizeof(float)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_sin_para_angs,
                          num_components*num_time_steps*sizeof(float)) );


      cudaErrorCheckCall( cudaMemcpy( d_cos_para_angs, cos_para_angs,
         num_components*num_time_steps*sizeof(float), cudaMemcpyHostToDevice) );
      cudaErrorCheckCall( cudaMemcpy( d_sin_para_angs, sin_para_angs,
         num_components*num_time_steps*sizeof(float), cudaMemcpyHostToDevice) );

      cudaErrorCheckKernel("kern_rotate_FEE_beam",
                            kern_rotate_FEE_beam, grid, threads,
                            (cuFloatComplex*)TileGainMatrices,
                            d_sin_para_angs, d_cos_para_angs,
                            num_coords);

      cudaErrorCheckCall( cudaFree(d_cos_para_angs) );
      cudaErrorCheckCall( cudaFree(d_sin_para_angs) );
    }
  }

  cudaErrorCheckCall( cudaFree(d_emn_T_sum_real) );
  cudaErrorCheckCall( cudaFree(d_emn_T_sum_imag) );
  cudaErrorCheckCall( cudaFree(d_emn_P_sum_real) );
  cudaErrorCheckCall( cudaFree(d_emn_P_sum_imag) );

  cudaErrorCheckCall( cudaFree(d_m_range) );

  cudaErrorCheckCall( cudaFree( d_phi ) );
  cudaErrorCheckCall( cudaFree( d_theta) );

  cudaErrorCheckCall( cudaFree( primary_beam->emn_T ) );
  cudaErrorCheckCall( cudaFree( primary_beam->emn_P ) );

  cudaErrorCheckCall( cudaFree( primary_beam->rts_P_sin ) );
  cudaErrorCheckCall( cudaFree( primary_beam->rts_P1 ) );

}

extern "C" void calc_CUDA_FEE_beam(float *azs, float *zas,
                                   float *sin_para_angs, float *cos_para_angs,
                                   int num_components, int num_time_steps,
                                   RTS_MWA_FEE_beam_t *FEE_beam,
                                   int rotation, int scaling) {

  // printf("\tDoing FEE beam tings\n");

  cudaErrorCheckCall( cudaMalloc( (void **)&FEE_beam->d_FEE_beam_gain_matrices, num_time_steps*num_components*MAX_POLS*sizeof(float _Complex)) );

  //Ensure gains are zero before summing results to them
  float _Complex *zero_array = NULL;
  zero_array = (float _Complex *)malloc( num_time_steps*num_components*MAX_POLS*sizeof(float _Complex) );

  for (int i = 0; i < num_time_steps*num_components*MAX_POLS; i++) {
    zero_array[i] = {0.0, 0.0};
  }
  //
  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_FEE_beam_gain_matrices,
                     zero_array,
                     num_time_steps*num_components*MAX_POLS*sizeof(float _Complex), cudaMemcpyHostToDevice ) );
  free(zero_array);

  RTS_CUDA_get_TileGains(azs, zas,
    sin_para_angs, cos_para_angs,
    num_time_steps, num_components, rotation, FEE_beam,
    FEE_beam->d_FEE_beam_gain_matrices, scaling);

}

__global__ void kern_make_norm_abs(cuFloatComplex *d_norm_fac) {

  const int pol = threadIdx.x + (blockDim.x*blockIdx.x);

  //If pol is less than the total number of polarisations * 4 normalisation
  //directions
  if (pol < MAX_POLS*MAX_POLS){
    cuFloatComplex norm = d_norm_fac[pol];
    cuFloatComplex abs_norm;

    abs_norm.x = sqrt(norm.x*norm.x + norm.y*norm.y);
    abs_norm.y = 0.0;

    d_norm_fac[pol] = abs_norm;

  }
}

extern "C" void get_HDFBeam_normalisation(RTS_MWA_FEE_beam_t *FEE_beam_zenith,
                RTS_MWA_FEE_beam_t *FEE_beam) {

  //The FEE beam is in theta-phi polarisations, which means we need a different
  //azimuth direction in north, east, south, west
  float norm_azs[4] = {M_PI/2.0, 0.0, M_PI, M_PI/2.0};
  float norm_zas[4] = {0.0, 0.0, 0.0, 0.0};
  int num_azza = 4;

  //When running the FEE beam, we later rotate by parrallactic angle. We
  //will normalise before this rotation, so only need to calculate the
  //beam values at zenith once, and the beam is stationary, so only need
  //to calculate for one time step here
  int num_time_steps = 1;

  //Not rotating normalisation by parallactic angle, so leave as NULL
  float *para_sinrot = NULL;
  float *para_cosrot = NULL;

  copy_FEE_primary_beam_to_GPU(FEE_beam_zenith);

  //Scaling 0 means don't apply normalistaion, as we are calculating it here
  int scaling = 0;
  //Also don't rotate as we will scale before rotating later
  int rotation = 0;

  calc_CUDA_FEE_beam(norm_azs, norm_zas,
         para_sinrot, para_cosrot,
         num_azza, num_time_steps, FEE_beam_zenith, rotation, scaling);

  //The FEE code spits out 4 polarisation values per direction on the sky,
  //so we have to calculate 16 complex beam gains, and later select the
  //correct outputs to get our 4 normlisation values
  float _Complex *all_norm_gains=NULL;
  all_norm_gains = (float _Complex*)malloc(num_time_steps*num_azza*MAX_POLS*sizeof(float _Complex));

  //Convert the normalisation results to their absolute values
  dim3 grid, threads;
  threads.x = num_time_steps*num_azza*MAX_POLS;
  grid.x = grid.y = grid.z = threads.y = threads.z = 1;

  cudaErrorCheckKernel("kern_make_norm_abs",
                       kern_make_norm_abs, grid, threads,
                       (cuFloatComplex*)FEE_beam_zenith->d_FEE_beam_gain_matrices);

  cudaErrorCheckCall( cudaMemcpy(all_norm_gains, FEE_beam_zenith->d_FEE_beam_gain_matrices,  num_time_steps*num_azza*MAX_POLS*sizeof(float _Complex), cudaMemcpyDeviceToHost )) ;

  for (size_t i = 0; i < MAX_POLS; i++) {
    //Do a polarisation reordering here, to follow the RTS code
    FEE_beam->norm_fac[i] = all_norm_gains[i*MAX_POLS + MAX_POLS - 1 - i];

  }
  free(all_norm_gains);
}

__global__ void kern_map_FEE_beam_gains(cuFloatComplex *d_FEE_beam_gain_matrices,
    cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
    cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
    int num_freqs, int num_components, int num_visis, int num_baselines,
    int num_times){

  //All baselines at all freqs and all times
  int iBaseline = threadIdx.x + (blockDim.x*blockIdx.x);
  //Direction on sky
  int iComponent = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iBaseline < num_visis && iComponent < num_components ) {

    //I hate indexing, good grief
    int time_ind = (int)floorf( (float)iBaseline / ((float)num_baselines * (float)num_freqs));
    int freq_ind = (int)floorf( ((float)iBaseline - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);

    int current_ind = iComponent*num_times + time_ind;
    int new_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + iComponent;

    d_primay_beam_J00[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 0];
    d_primay_beam_J01[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 1];
    d_primay_beam_J10[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 2];
    d_primay_beam_J11[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 3];

  }
}

extern "C" void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam){
  cudaErrorCheckCall( cudaFree(primary_beam->d_M) );
  cudaErrorCheckCall( cudaFree(primary_beam->d_N) );
  cudaErrorCheckCall( cudaFree(primary_beam->d_Q1) );
  cudaErrorCheckCall( cudaFree(primary_beam->d_Q2) );
}

extern "C" void test_RTS_CUDA_FEE_beam(int num_components,
           float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           RTS_MWA_FEE_beam_t *FEE_beam_zenith,
           RTS_MWA_FEE_beam_t *FEE_beam,
           int rotation, int scaling,
           float _Complex *FEE_beam_gains){

  int num_time_steps = 1;

  printf("Getting FEE beam normalisation...\n");
  get_HDFBeam_normalisation(FEE_beam_zenith, FEE_beam);
  printf(" done.\n");

  free_FEE_primary_beam_from_GPU(FEE_beam_zenith);

  printf("Copying the FEE beam across to the GPU...");
  copy_FEE_primary_beam_to_GPU(FEE_beam);
  printf(" done.\n");

  printf("Calculating FEE beam for %d directions...\n", num_components);
  calc_CUDA_FEE_beam(azs, zas, sin_para_angs, cos_para_angs,
         num_components, num_time_steps, FEE_beam,
         rotation, scaling);

  cudaErrorCheckCall( cudaMemcpy(FEE_beam_gains,
            FEE_beam->d_FEE_beam_gain_matrices,
            num_components*num_time_steps*MAX_POLS*sizeof(float _Complex),
            cudaMemcpyDeviceToHost) );

  // free_FEE_primary_beam_from_GPU(FEE_beam_zenith);
  // printf("GPU beam realeased, calculation complete\n");

}
