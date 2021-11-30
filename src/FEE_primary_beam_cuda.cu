#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuComplex.h>
#include <complex.h>
#include <assert.h>
#include <stdio.h>
#include <erfa.h>

#include "cudacomplex.h"
#include "cudacheck.h"
#include "FEE_primary_beam_cuda.h"

/**
Factorials used in calculating legendre polynomials
*/
__constant__ user_precision_t ffactorials[] = {
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

  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_M, arrSize*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_N, arrSize*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_Q1, arrSize*sizeof(user_precision_complex_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&FEE_beam->d_Q2, arrSize*sizeof(user_precision_complex_t) ) );

  user_precision_t *h_params;

  h_params = (user_precision_t *)malloc(arrSize*sizeof(user_precision_t) );

  for(int i=0; i< n_pols; i++){
    for(int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (user_precision_t)(FEE_beam->M[i][j]);
    }
  }

  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_M, h_params,
                      arrSize*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  for(int i=0; i<n_pols; i++){
    for(int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (user_precision_t)(FEE_beam->N[i][j]);
    }
  }

  cudaErrorCheckCall( cudaMemcpy( FEE_beam->d_N, h_params,
                      arrSize*nStations*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );
  // cudaFreeHost( h_params );
  free(h_params);

  user_precision_complex_t *h_Qdata;
  h_Qdata = (user_precision_complex_t *)malloc(arrSize*sizeof(user_precision_complex_t) );

  for(int i=0; i< n_pols; i++){
    for(int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = FEE_beam->Q1[i][j];
    }
  }

  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_Q1, h_Qdata,
                      arrSize*nStations*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ) );

  for(int i=0; i< n_pols; i++){
    for(int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = FEE_beam->Q2[i][j];
    }
  }

  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_Q2, h_Qdata,
                      arrSize*nStations*sizeof(user_precision_complex_t) , cudaMemcpyHostToDevice ) );
  free( h_Qdata );

}

// // GPU version of legendre polynomial function. This version only returns the value
// // at a single point as required in HDF beam calculation
// // See legendre_polynomial.c for source and info
__device__ void RTS_CUDA_pm_polynomial_value_single(int n, int m,
                               user_precision_t x, user_precision_t *values  ){

  user_precision_t fact;
  int i;
  int j;
  int k;
  int mm=1;

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
      values[i+(m+1)*mm] = x * ( user_precision_t ) ( 2 * m + 1 ) * values[i+m*mm];
    }
  }
  /*
    Now we use a three term recurrence.
  */
  for ( j = m + 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      values[i+j*mm] = ( ( user_precision_t ) ( 2 * j     - 1 ) * x * values[i+(j-1)*mm]
                  + ( user_precision_t ) (   - j - m + 1 ) *        values[i+(j-2)*mm] )
                  / ( user_precision_t ) (     j - m     );
    }
  }
}

__global__ void RTS_P1SINfKernel(user_precision_t *d_theta, cuUserComplex *rts_P_sin,
           cuUserComplex *rts_p1, int nmax, int num_coords){

  int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  user_precision_t theta = d_theta[iCoord];

  //RTS style variables
  int i = threadIdx.y + (blockDim.y*blockIdx.y);
  int n = floor(sqrt((user_precision_t)(i+1)));
  int index = i - ((n-1)*(n-1) + 2*(n-1));
  int pm_index;
  user_precision_t u = cos(theta);
  user_precision_t sin_th = sin(theta);
  user_precision_t delu = 1.0e-6;

  int nMN = nmax*nmax + 2*nmax;

  if(iCoord < num_coords && i < nMN) {

    if(index < n){
      pm_index = n-index;
    } else {
      pm_index = index-n;
    }

    user_precision_t p;
    user_precision_t Pm1;
    user_precision_t Pm_sin;

    user_precision_t pm_vals[100];

    RTS_CUDA_pm_polynomial_value_single(n, pm_index, u, pm_vals);
    p = pm_vals[n];
    // Diverging??
    RTS_CUDA_pm_polynomial_value_single(n, pm_index+1, u, pm_vals);
    if(pm_index==n){
      Pm1 = 0.0;
    } else {
      Pm1 = pm_vals[n];
    }
    Pm_sin = 0.0;

    if(u==1.0){
      user_precision_t pm_val;
      RTS_CUDA_pm_polynomial_value_single(n, 0, u-delu, pm_vals);
      pm_val = pm_vals[n];
      user_precision_t Pm_sin1,p0;
      RTS_CUDA_pm_polynomial_value_single(0, 0, u, pm_vals);
      p0 = pm_vals[0];
      Pm_sin1 = -(p0 - pm_val) / delu;
      if(pm_index==1) Pm_sin = Pm_sin1;
    } else if (u == -1.0){
      user_precision_t pm_val;
      RTS_CUDA_pm_polynomial_value_single(n, 0, u-delu, pm_vals);
      pm_val = pm_vals[n];
      user_precision_t Pm_sin1,p0;
      RTS_CUDA_pm_polynomial_value_single(0, 0, u, pm_vals);
      p0 = pm_vals[0];
      Pm_sin1 = -(pm_val - p0) / delu;
      if(pm_index==1) Pm_sin = Pm_sin1;
    } else {
      Pm_sin = p / sin_th;
    }

    rts_P_sin[iCoord*nMN + i] = make_cuUserComplex(Pm_sin,0);
    rts_p1[iCoord*nMN + i] = make_cuUserComplex(Pm1,0);

  } //end IF thread lies within number of calculations we want to do
}

__global__ void RTS_getTileGainsKernel( user_precision_t *d_phi, user_precision_t *d_theta, int nMN, int num_coords,
           user_precision_t *pb_M, user_precision_t *pb_N,
           cuUserComplex *pb_Q1, cuUserComplex *pb_Q2,
           cuUserComplex *rts_P_sin, cuUserComplex *rts_P1,
           cuUserComplex *emn_T, cuUserComplex *emn_P){


  int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  int i = threadIdx.y + (blockDim.y*blockIdx.y);

  int pol = blockIdx.z;
  int n_pols = gridDim.z;

  int index;
  user_precision_t factor1,factor2;
  user_precision_t C_MN, M_absM;

  int beam_index = i+(pol*nMN);

  if(iCoord < num_coords && i < nMN) {
    user_precision_t phi = d_phi[iCoord];
    user_precision_t theta = d_theta[iCoord];

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

    cuUserComplex phi_comp;

    phi_comp = U1polar(pb_M[beam_index] * phi) * C_MN * M_absM / sqrt(pb_N[beam_index] * (pb_N[beam_index] + 1.0));

    cuUserComplex this_emn_T, this_emn_P;
    user_precision_t T_exp, P_exp;
    user_precision_t u = cos(theta);

    T_exp = (pb_N[beam_index]);
    P_exp = (pb_N[beam_index])+1.0;

    this_emn_T = make_cuUserComplex(cos(T_exp * M_PI/2.0),sin(T_exp * M_PI/2.0));
    this_emn_P = make_cuUserComplex(cos(P_exp * M_PI/2.0),sin(P_exp * M_PI/2.0));

    this_emn_T *= (rts_P_sin[iCoord*nMN + i]*(abs(pb_M[beam_index])*pb_Q2[beam_index]*u - pb_M[beam_index] * pb_Q1[beam_index]) + pb_Q2[beam_index] * rts_P1[iCoord*nMN + i]);
    this_emn_P *= (rts_P_sin[iCoord*nMN + i]*(pb_M[beam_index]*pb_Q2[beam_index] - abs(pb_M[beam_index])*pb_Q1[beam_index]*u) - pb_Q1[beam_index] * rts_P1[iCoord*nMN + i]);

    // emn_T[iCoord*n_pols*nMN + beam_index] = cuCmulf(this_emn_T,phi_comp);
    // emn_P[iCoord*n_pols*nMN + beam_index] = cuCmulf(this_emn_P,phi_comp);

    emn_T[iCoord*n_pols*nMN + beam_index] = this_emn_T*phi_comp;
    emn_P[iCoord*n_pols*nMN + beam_index] = this_emn_P*phi_comp;

  }
}

__global__ void kern_sum_emn_PT_by_M(cuUserComplex *emn_T, cuUserComplex *emn_P,
           user_precision_t *d_emn_T_sum_real, user_precision_t *d_emn_T_sum_imag,
           user_precision_t *d_emn_P_sum_real, user_precision_t *d_emn_P_sum_imag,
           user_precision_t *d_m_range, user_precision_t *d_M, int nMN, int nmax, int num_coords){

  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iM_value = threadIdx.y + (blockDim.y*blockIdx.y);

  int n_pols = 2;
  int num_sum_M = 2*nmax + 1;

  if(iCoord < num_coords && iM_value < num_sum_M) {

    for (int iM_index = 0; iM_index < n_pols*nMN; iM_index++) {
      int iPol = (int)floorf((user_precision_t)iM_index / (user_precision_t)nMN);

      if (d_M[iM_index] == d_m_range[iM_value]) {
        cuUserComplex emn_P_value = emn_P[n_pols*nMN*iCoord + iM_index];
        cuUserComplex emn_T_value = emn_T[n_pols*nMN*iCoord + iM_index];

        d_emn_P_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_P_value.x;
        d_emn_P_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_P_value.y;

        d_emn_T_sum_real[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_T_value.x;
        d_emn_T_sum_imag[n_pols*num_sum_M*iCoord + iPol*num_sum_M + (int)d_m_range[iM_value] + nmax] += emn_T_value.y;
      }
    }
  }
}

__global__ void kern_map_emn(cuUserComplex *TileGainMatrices,
                user_precision_t *d_emn_T_sum_real, user_precision_t *d_emn_T_sum_imag,
                user_precision_t *d_emn_P_sum_real, user_precision_t *d_emn_P_sum_imag,
                int nmax, int num_coords) {

  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iPol = threadIdx.y + (blockDim.y*blockIdx.y);
  int n_pols=2;

  int numM = 2*nmax + 1;

  if (iCoord < num_coords && iPol < n_pols){

    for(int i=0; i < numM; i++){

      int array_loc = n_pols*numM*iCoord + iPol*numM + i;

      cuUserComplex d_emn_T = make_cuUserComplex(d_emn_T_sum_real[array_loc],d_emn_T_sum_imag[array_loc]);
      cuUserComplex d_emn_P = make_cuUserComplex(d_emn_P_sum_real[array_loc],d_emn_P_sum_imag[array_loc]);

      TileGainMatrices[iCoord*MAX_POLS + iPol*n_pols] += d_emn_T;
      TileGainMatrices[iCoord*MAX_POLS + iPol*n_pols + 1] += d_emn_P;

    }

  }
}

__global__ void kern_apply_FEE_norm(cuUserComplex *TileGainMatrices,
           cuUserComplex *d_norm_fac, int num_coords ){

  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iCoord < num_coords) {

    cuUserComplex before0 = TileGainMatrices[iCoord*MAX_POLS + 0];
    cuUserComplex before1 = TileGainMatrices[iCoord*MAX_POLS + 1];
    cuUserComplex before2 = TileGainMatrices[iCoord*MAX_POLS + 2];
    cuUserComplex before3 = TileGainMatrices[iCoord*MAX_POLS + 3];

    //Apply normalistation
    TileGainMatrices[iCoord*MAX_POLS + 0] = before3 / d_norm_fac[0];
    TileGainMatrices[iCoord*MAX_POLS + 1] = before2 / d_norm_fac[1];
    TileGainMatrices[iCoord*MAX_POLS + 2] = before1 / d_norm_fac[2];
    TileGainMatrices[iCoord*MAX_POLS + 3] = before0 / d_norm_fac[3];

  }
}

__global__ void kern_rotate_FEE_beam(cuUserComplex *d_FEE_beam_gain_matrices,
                                user_precision_t *d_sin_para_angs, user_precision_t *d_cos_para_angs,
                                int num_coords) {
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iCoord < num_coords) {

    cuUserComplex prerot0 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 0];
    cuUserComplex prerot1 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 1];
    cuUserComplex prerot2 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 2];
    cuUserComplex prerot3 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 3];

    user_precision_t sinrot = d_sin_para_angs[iCoord];
    user_precision_t cosrot = d_cos_para_angs[iCoord];

    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 0] = -prerot1*sinrot + prerot0*cosrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 1] = prerot1*cosrot + prerot0*sinrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 2] = -prerot3*sinrot + prerot2*cosrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 3] = prerot3*cosrot + prerot2*sinrot;

  }
}

extern "C" void RTS_CUDA_get_TileGains(user_precision_t *phi, user_precision_t *theta,
           user_precision_t *sin_para_angs, user_precision_t *cos_para_angs,
           int num_time_steps, int num_components,
           user_precision_t rotation, RTS_MWA_FEE_beam_t *primary_beam,
           user_precision_complex_t *TileGainMatrices, int scaling){

  int n_pols = 2;
  int num_coords = num_time_steps*num_components;

  for (int phi_ind = 0; phi_ind < num_coords; phi_ind++) {
    phi[phi_ind] = (M_PI/2.0) - phi[phi_ind];
    if(phi[phi_ind] < 0){
      phi[phi_ind] += 2.0*M_PI;
    }
  }

  user_precision_t *d_phi=NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_phi, num_coords*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy(d_phi, phi, num_coords*sizeof(user_precision_t),
                      cudaMemcpyHostToDevice ) );

  user_precision_t *d_theta=NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_theta, num_coords*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMemcpy(d_theta, theta, num_coords*sizeof(user_precision_t),
                      cudaMemcpyHostToDevice ) );

  int nmax = primary_beam->nmax; // Assumes all beams have same dimensions
  int nMN = primary_beam->nMN;

  int P_size = sizeof(user_precision_complex_t)*(nmax*nmax + 2*nmax)*num_coords;

  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->rts_P_sin, P_size) );
  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->rts_P1, P_size) );

  dim3 grid, threads;

  threads.y = 2;
  grid.y = (int)ceil((user_precision_t)(nmax*nmax + 2*nmax) / (user_precision_t)threads.y);

  threads.x = 64;
  grid.x = (int)ceil((user_precision_t)num_coords / (user_precision_t)threads.x);

  cudaErrorCheckKernel("RTS_P1SINfKernel",
                       RTS_P1SINfKernel, grid, threads,
                       d_theta, (cuUserComplex*)primary_beam->rts_P_sin,
                       (cuUserComplex*)primary_beam->rts_P1, nmax, num_coords);

  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->emn_P,
                    num_coords * nMN * n_pols * sizeof(user_precision_complex_t)) );
  cudaErrorCheckCall( cudaMalloc( (void **)&primary_beam->emn_T,
                    num_coords * nMN * n_pols * sizeof(user_precision_complex_t)) );
  //
  //
  // threads.x = 128;
  threads.x = 128;
  threads.y = 2;
  threads.z = 1;
  //
  grid.x = (int)ceil((user_precision_t)num_coords / (user_precision_t)threads.x);
  grid.y = (int)ceil((user_precision_t)(nMN) / (user_precision_t)threads.y);
  grid.z = n_pols;

  cudaErrorCheckKernel("RTS_getTileGainsKernel",
                        RTS_getTileGainsKernel, grid, threads,
                        d_phi, d_theta, nMN, num_coords,
                        primary_beam->d_M,primary_beam->d_N,
                        (cuUserComplex*)primary_beam->d_Q1,
                        (cuUserComplex*)primary_beam->d_Q2,
                        (cuUserComplex*)primary_beam->rts_P_sin,
                        (cuUserComplex*)primary_beam->rts_P1,
                        (cuUserComplex*)primary_beam->emn_T,
                        (cuUserComplex*)primary_beam->emn_P);

  user_precision_t *d_emn_P_sum_real = NULL;
  user_precision_t *d_emn_T_sum_real = NULL;
  user_precision_t *d_emn_P_sum_imag = NULL;
  user_precision_t *d_emn_T_sum_imag = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_P_sum_real,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_T_sum_real,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_P_sum_imag,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_emn_T_sum_imag,
                      num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t) ) );

  //Ensure gains are zero before summing results to them
  user_precision_t *zero_array = NULL;
  zero_array = (user_precision_t *)malloc(num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t) );

  for (int i = 0; i < num_coords*(2*nmax + 1)*n_pols; i++) {
    zero_array[i] = 0.0;
  }

  cudaErrorCheckCall( cudaMemcpy(d_emn_P_sum_real, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t), cudaMemcpyHostToDevice) );
  cudaErrorCheckCall( cudaMemcpy(d_emn_P_sum_imag, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t), cudaMemcpyHostToDevice) );
  cudaErrorCheckCall( cudaMemcpy(d_emn_T_sum_real, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t), cudaMemcpyHostToDevice) );
  cudaErrorCheckCall( cudaMemcpy(d_emn_T_sum_imag, zero_array,
         num_coords*(2*nmax + 1)*n_pols*sizeof(user_precision_t), cudaMemcpyHostToDevice) );

  user_precision_t *d_m_range = NULL;
  cudaErrorCheckCall( cudaMalloc( (void**)&d_m_range, (2*nmax + 1)*sizeof(user_precision_t)) );
  cudaErrorCheckCall( cudaMemcpy(d_m_range, primary_beam->m_range,
                      (2*nmax + 1)*sizeof(user_precision_t), cudaMemcpyHostToDevice ) );

  threads.x = 8;
  threads.y = 32;
  threads.z = 1;
  //
  grid.x = (int)ceil( (user_precision_t)num_coords / threads.x);
  grid.z = (int)ceil( ((user_precision_t)nmax*2.0 + 1.0) / threads.y);
  grid.z = 1;

  cudaErrorCheckKernel("kern_sum_emn_PT_by_M",
                        kern_sum_emn_PT_by_M, grid, threads,
                        (cuUserComplex*)primary_beam->emn_T,
                        (cuUserComplex*)primary_beam->emn_P,
                        d_emn_T_sum_real, d_emn_T_sum_imag,
                        d_emn_P_sum_real, d_emn_P_sum_imag,
                        d_m_range, primary_beam->d_M, primary_beam->nMN,
                        nmax, num_coords);

  threads.x = 64;
  threads.y = n_pols;
  threads.z = 1;

  grid.x = (int)ceil( (user_precision_t)num_coords / threads.x);
  grid.y = grid.z = 1;

  free(zero_array);

  cudaErrorCheckKernel("kern_map_emn",
                        kern_map_emn, grid, threads,
                        (cuUserComplex*)TileGainMatrices,
                        d_emn_T_sum_real, d_emn_T_sum_imag,
                        d_emn_P_sum_real, d_emn_P_sum_imag,
                        nmax, num_coords );

  //If we want to normalise the beam, do it here
  if (scaling == 1.0) {

    threads.x = 128;
    grid.x = (int)ceil( (user_precision_t)num_coords / threads.x);

    threads.y = threads.z = 1;
    grid.y = grid.z = 1;

    user_precision_complex_t *d_norm_fac=NULL;
    cudaErrorCheckCall( cudaMalloc( (void**)&d_norm_fac,
                        MAX_POLS*sizeof(user_precision_complex_t) ) );
    cudaErrorCheckCall( cudaMemcpy(d_norm_fac, primary_beam->norm_fac,
               MAX_POLS*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ) );

    cudaErrorCheckKernel("kern_apply_FEE_norm",
                         kern_apply_FEE_norm, grid, threads,
                         (cuUserComplex*)TileGainMatrices,
                         (cuUserComplex*)d_norm_fac, num_coords );

    cudaErrorCheckCall( cudaFree(d_norm_fac) );

    //If rotating by parallactic angle - only ever want to do this
    //after normalistation
    if (rotation == 1.0) {

      user_precision_t *d_cos_para_angs=NULL;
      user_precision_t *d_sin_para_angs=NULL;

      cudaErrorCheckCall( cudaMalloc( (void**)&d_cos_para_angs,
                          num_components*num_time_steps*sizeof(user_precision_t)) );
      cudaErrorCheckCall( cudaMalloc( (void**)&d_sin_para_angs,
                          num_components*num_time_steps*sizeof(user_precision_t)) );


      cudaErrorCheckCall( cudaMemcpy( d_cos_para_angs, cos_para_angs,
         num_components*num_time_steps*sizeof(user_precision_t), cudaMemcpyHostToDevice) );
      cudaErrorCheckCall( cudaMemcpy( d_sin_para_angs, sin_para_angs,
         num_components*num_time_steps*sizeof(user_precision_t), cudaMemcpyHostToDevice) );

      cudaErrorCheckKernel("kern_rotate_FEE_beam",
                            kern_rotate_FEE_beam, grid, threads,
                            (cuUserComplex*)TileGainMatrices,
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

  //undo the azimuth rotation in case user keeps using the same az/za
  //otherwise we'll just keep spinning around and around and around
  for (int phi_ind = 0; phi_ind < num_coords; phi_ind++) {
    phi[phi_ind] = (M_PI/2.0) - phi[phi_ind];
    if(phi[phi_ind] < 0){
      phi[phi_ind] += 2.0*M_PI;
    }
  }

}

extern "C" void calc_CUDA_FEE_beam(user_precision_t *azs, user_precision_t *zas,
                                   user_precision_t *sin_para_angs, user_precision_t *cos_para_angs,
                                   int num_components, int num_time_steps,
                                   RTS_MWA_FEE_beam_t *FEE_beam,
                                   int rotation, int scaling) {

  // printf("\tDoing FEE beam tings\n");
  cudaErrorCheckCall( cudaMalloc( (void **)&FEE_beam->d_FEE_beam_gain_matrices, num_time_steps*num_components*MAX_POLS*sizeof(user_precision_complex_t)) );

  //Ensure gains are zero before summing results to them
  user_precision_complex_t *zero_array = NULL;
  zero_array = (user_precision_complex_t *)malloc( num_time_steps*num_components*MAX_POLS*sizeof(user_precision_complex_t) );

  for (int i = 0; i < num_time_steps*num_components*MAX_POLS; i++) {
    zero_array[i] = {0.0, 0.0};
  }
  //
  cudaErrorCheckCall( cudaMemcpy(FEE_beam->d_FEE_beam_gain_matrices,
                     zero_array,
                     num_time_steps*num_components*MAX_POLS*sizeof(user_precision_complex_t), cudaMemcpyHostToDevice ) );
  free(zero_array);

  RTS_CUDA_get_TileGains(azs, zas,
    sin_para_angs, cos_para_angs,
    num_time_steps, num_components, rotation, FEE_beam,
    FEE_beam->d_FEE_beam_gain_matrices, scaling);

}

__global__ void kern_make_norm_abs(cuUserComplex *d_norm_fac) {

  const int pol = threadIdx.x + (blockDim.x*blockIdx.x);

  //If pol is less than the total number of polarisations * 4 normalisation
  //directions
  if (pol < MAX_POLS*MAX_POLS){
    cuUserComplex norm = d_norm_fac[pol];
    cuUserComplex abs_norm;

    abs_norm.x = sqrt(norm.x*norm.x + norm.y*norm.y);
    abs_norm.y = 0.0;

    d_norm_fac[pol] = abs_norm;

  }
}

extern "C" void get_HDFBeam_normalisation(RTS_MWA_FEE_beam_t *FEE_beam_zenith,
                RTS_MWA_FEE_beam_t *FEE_beam) {

  //The FEE beam is in theta-phi polarisations, which means we need a different
  //azimuth direction in north, east, south, west
  user_precision_t norm_azs[4] = {M_PI/2.0, 0.0, M_PI, M_PI/2.0};
  user_precision_t norm_zas[4] = {0.0, 0.0, 0.0, 0.0};
  int num_azza = 4;

  //When running the FEE beam, we later rotate by parrallactic angle. We
  //will normalise before this rotation, so only need to calculate the
  //beam values at zenith once, and the beam is stationary, so only need
  //to calculate for one time step here
  int num_time_steps = 1;

  //Not rotating normalisation by parallactic angle, so leave as NULL
  user_precision_t *para_sinrot = NULL;
  user_precision_t *para_cosrot = NULL;

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
  user_precision_complex_t *all_norm_gains=NULL;
  all_norm_gains = (user_precision_complex_t*)malloc(num_time_steps*num_azza*MAX_POLS*sizeof(user_precision_complex_t));

  //Convert the normalisation results to their absolute values
  dim3 grid, threads;
  threads.x = num_time_steps*num_azza*MAX_POLS;
  grid.x = grid.y = grid.z = threads.y = threads.z = 1;

  cudaErrorCheckKernel("kern_make_norm_abs",
                       kern_make_norm_abs, grid, threads,
                       (cuUserComplex*)FEE_beam_zenith->d_FEE_beam_gain_matrices);

  cudaErrorCheckCall( cudaMemcpy(all_norm_gains, FEE_beam_zenith->d_FEE_beam_gain_matrices,  num_time_steps*num_azza*MAX_POLS*sizeof(user_precision_complex_t), cudaMemcpyDeviceToHost )) ;

  for (size_t i = 0; i < MAX_POLS; i++) {
    //Do a polarisation reordering here, to follow the RTS code
    FEE_beam->norm_fac[i] = all_norm_gains[i*MAX_POLS + MAX_POLS - 1 - i];

  }
  free(all_norm_gains);
}

__global__ void kern_map_FEE_beam_gains(cuUserComplex *d_FEE_beam_gain_matrices,
    cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
    cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11,
    int num_freqs, int num_components, int num_visis, int num_baselines,
    int num_times){

  //All baselines at all freqs and all times
  int iComponent = threadIdx.x + (blockDim.x*blockIdx.x);
  //time step
  int iTime = threadIdx.y + (blockDim.y*blockIdx.y);

  if(iComponent < num_components && iTime < num_times ) {

    //Where this component and time step lives in the current array
    int current_ind = iComponent*num_times + iTime;

    //For the FEE beam currently, only have resolution of 1.28MHz, so assign
    //the same beam value to all frequencies in this band (hence the loop)
    for (int freq_ind = 0; freq_ind < num_freqs; freq_ind++) {
      //Get an index to shove into the generic beam jones containers
      int new_ind = num_freqs*iTime*num_components + (num_components*freq_ind) + iComponent;
      //Split up the polarisations
      d_primay_beam_J00[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 0];
      d_primay_beam_J01[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 1];
      d_primay_beam_J10[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 2];
      d_primay_beam_J11[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 3];
    }
  }
}

extern "C" void free_FEE_primary_beam_from_GPU(RTS_MWA_FEE_beam_t *primary_beam){
  cudaErrorCheckCall( cudaFree(primary_beam->d_M) );
  cudaErrorCheckCall( cudaFree(primary_beam->d_N) );
  cudaErrorCheckCall( cudaFree(primary_beam->d_Q1) );
  cudaErrorCheckCall( cudaFree(primary_beam->d_Q2) );
}

extern "C" void multifreq_get_MWAFEE_normalisation(beam_settings_t *beam_settings){

  printf("Getting normalisations for MWA FEE beams...");

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];
    RTS_MWA_FEE_beam_t *FEE_beam_zenith = &beam_settings->FEE_beam_zeniths[freq_ind];

    //Only needs to be done once per frequency as zenith never moves, so do
    //this outside the sky model chunk loop
    get_HDFBeam_normalisation(FEE_beam_zenith, FEE_beam);
    //Free the zenith pointing as done with it now
    free_FEE_primary_beam_from_GPU(FEE_beam_zenith);

    copy_FEE_primary_beam_to_GPU(FEE_beam);

  }

  printf(" done\n");
}

extern "C" void multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
           user_precision_t *azs, user_precision_t *zas, int num_time_steps,
           user_precision_t *sin_para_angs, user_precision_t *cos_para_angs,
           int num_components, int rotation, int scaling){

  // printf("Calculating MWA FEE beams...");

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];

    calc_CUDA_FEE_beam(azs, zas, sin_para_angs, cos_para_angs,
             num_components, num_time_steps, FEE_beam,
             rotation, scaling);

  }
  // printf(" done\n");
}

__global__ void kern_map_FEE_beam_gains_multi_freq(cuUserComplex *d_FEE_beam_gain_matrices,
           cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
           cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11,
           int num_freqs, int num_components, int num_times, int iFreq){

  int iComp = threadIdx.x + (blockDim.x*blockIdx.x);
  int iTime = threadIdx.y + (blockDim.y*blockIdx.y);

  if ((iComp < num_components) && (iTime < num_times)) {
    int current_ind = iComp*num_times + iTime;
    int new_ind = num_freqs*iTime*num_components + (num_components*iFreq) + iComp;

    d_primay_beam_J00[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 0];
    d_primay_beam_J01[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 1];
    d_primay_beam_J10[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 2];
    d_primay_beam_J11[new_ind] = d_FEE_beam_gain_matrices[current_ind*MAX_POLS + 3];
  }
}


extern "C" void map_FEE_beam_gains_multi_freq(beam_settings_t *beam_settings,
    cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
    cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11,
    int num_freqs, int num_components, int num_times){

  dim3 grid, threads;
  threads.x = 64;
  threads.y = 2;

  grid.x = (int)ceil( (float)num_components / (float)threads.x );
  grid.y = (int)ceil( (float)num_times / (float)threads.y );
  grid.z = threads.z = 1;

  //Loop over frequencies. This way only doing the conversion of user_precision_complex_t
  //to cuUserComplex on d_FEE_beam_gain_matrices outside the kernel
  for (int iFreq = 0; iFreq < num_freqs; iFreq++) {

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[iFreq];

    cudaErrorCheckKernel("kern_map_FEE_beam_gains_multi_freq",
                        kern_map_FEE_beam_gains_multi_freq, grid, threads,
                        (cuUserComplex *)FEE_beam->d_FEE_beam_gain_matrices,
                        d_primay_beam_J00, d_primay_beam_J01,
                        d_primay_beam_J10, d_primay_beam_J11,
                        num_freqs, num_components, num_times, iFreq);
  }
}

extern "C" void run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
        user_precision_t *azs, user_precision_t *zas,
        user_precision_t *sin_para_angs, user_precision_t *cos_para_angs,
        int num_components, int num_freqs, int num_times,
        int rotation, int scaling,
        cuUserComplex *d_primay_beam_J00, cuUserComplex *d_primay_beam_J01,
        cuUserComplex *d_primay_beam_J10, cuUserComplex *d_primay_beam_J11){

  // printf("Calculating MWA FEE beams...");

  multifreq_calc_CUDA_FEE_beam(beam_settings, azs, zas, num_times,
                              sin_para_angs, cos_para_angs,
                              num_components, rotation, scaling);
  // printf(" done\n");

  // printf("Mapping MWA FEE beams...\n");

  map_FEE_beam_gains_multi_freq(beam_settings,
                                d_primay_beam_J00, d_primay_beam_J01,
                                d_primay_beam_J10, d_primay_beam_J11,
                                num_freqs, num_components, num_times);

  // printf(" done\n");
}








/*******************************************************************************
Testing code below
*******************************************************************************/

extern "C" void test_RTS_CUDA_FEE_beam(int num_components,
           user_precision_t *azs, user_precision_t *zas, double latitude,
           RTS_MWA_FEE_beam_t *FEE_beam_zenith,
           RTS_MWA_FEE_beam_t *FEE_beam,
           int rotation, int scaling,
           user_precision_complex_t *FEE_beam_gains){

  //Just taking in a list of az/za, so set to single time step
  //Don;t need any fancy ordering of outputs
  int num_time_steps = 1;

  //Get the parallatic angles
  user_precision_t *sin_para_angs = NULL;
  user_precision_t *cos_para_angs = NULL;

  sin_para_angs = (user_precision_t *)malloc(num_components*sizeof(user_precision_t));
  cos_para_angs = (user_precision_t *)malloc(num_components*sizeof(user_precision_t));

  double ha, dec, el, para_angle;

  for (int comp = 0; comp < num_components; comp++) {
    el = M_PI/2.0 - (double)zas[comp];

    eraAe2hd((double)azs[comp], el, latitude, &ha, &dec);
    para_angle = eraHd2pa(ha, dec, latitude);

    sin_para_angs[comp] = (user_precision_t)sin(para_angle + M_PI/2.0);
    cos_para_angs[comp] = (user_precision_t)cos(para_angle + M_PI/2.0);

  }

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
            num_components*num_time_steps*MAX_POLS*sizeof(user_precision_complex_t),
            cudaMemcpyDeviceToHost) );

  free_FEE_primary_beam_from_GPU(FEE_beam);
  printf("GPU beam realeased, calculation complete\n");

}


extern "C" void test_multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
                            user_precision_t *azs, user_precision_t *zas,
                            int num_time_steps, double latitude,
                            int num_components, int rotation, int scaling,
                            user_precision_complex_t *all_FEE_beam_gains){


  //Get the parallatic angles
  user_precision_t *sin_para_angs = NULL;
  user_precision_t *cos_para_angs = NULL;

  sin_para_angs = (user_precision_t *)malloc(num_components*sizeof(user_precision_t));
  cos_para_angs = (user_precision_t *)malloc(num_components*sizeof(user_precision_t));

  double ha, dec, el, para_angle;

  for (int comp = 0; comp < num_components; comp++) {
    el = M_PI/2.0 - zas[comp];

    eraAe2hd((double)azs[comp], el, latitude, &ha, &dec);
    para_angle = eraHd2pa(ha, dec, latitude);

    sin_para_angs[comp] = (user_precision_t)sin(para_angle + M_PI/2.0);
    cos_para_angs[comp] = (user_precision_t)cos(para_angle + M_PI/2.0);

  }

  multifreq_calc_CUDA_FEE_beam(beam_settings, azs, zas, num_time_steps,
                              sin_para_angs, cos_para_angs,
                              num_components, rotation, scaling);

  for (int freq_ind = 0; freq_ind < beam_settings->num_MWAFEE; freq_ind++) {

    RTS_MWA_FEE_beam_t *FEE_beam = &beam_settings->FEE_beams[freq_ind];

    calc_CUDA_FEE_beam(azs, zas, sin_para_angs, cos_para_angs,
             num_components, num_time_steps, FEE_beam,
             rotation, scaling);

    //Iterate the pointer to all_FEE_beam_gains and copy across the gains from
    //the GPU into the chunk of memory
    int pointer_iter = num_time_steps*num_components*MAX_POLS*freq_ind;

    cudaErrorCheckCall( cudaMemcpy(all_FEE_beam_gains + pointer_iter,
              FEE_beam->d_FEE_beam_gain_matrices,
              num_time_steps*num_components*MAX_POLS*sizeof(user_precision_complex_t),
              cudaMemcpyDeviceToHost) );

    cudaErrorCheckCall( cudaFree(FEE_beam->d_FEE_beam_gain_matrices) );

  }
  free(sin_para_angs);
  free(cos_para_angs);
  printf(" done\n");
}



extern "C" void test_run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings_t *beam_settings,
    user_precision_t *azs, user_precision_t *zas, double latitude,
    user_precision_complex_t *primay_beam_J00,
    user_precision_complex_t *primay_beam_J01,
    user_precision_complex_t *primay_beam_J10,
    user_precision_complex_t *primay_beam_J11,
    int num_freqs, int num_components, int num_times,
    int rotation, int scaling){

  //Get the parallatic angles
  user_precision_t *sin_para_angs = NULL;
  user_precision_t *cos_para_angs = NULL;

  sin_para_angs = (user_precision_t *)malloc(num_times*num_components*sizeof(user_precision_t));
  cos_para_angs = (user_precision_t *)malloc(num_times*num_components*sizeof(user_precision_t));

  double ha, dec, el, para_angle;

  for (int comp = 0; comp < num_components*num_times; comp++) {
    el = M_PI/2.0 - zas[comp];

    eraAe2hd((double)azs[comp], el, latitude, &ha, &dec);
    para_angle = eraHd2pa(ha, dec, latitude);

    sin_para_angs[comp] = (user_precision_t)sin(para_angle + M_PI/2.0);
    cos_para_angs[comp] = (user_precision_t)cos(para_angle + M_PI/2.0);

  }

  int num_beam_values = num_freqs*num_times*num_components;

  user_precision_complex_t *d_primay_beam_J00 = NULL;
  user_precision_complex_t *d_primay_beam_J01 = NULL;
  user_precision_complex_t *d_primay_beam_J10 = NULL;
  user_precision_complex_t *d_primay_beam_J11 = NULL;

  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J00, num_beam_values*sizeof(user_precision_complex_t ) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J01, num_beam_values*sizeof(user_precision_complex_t ) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J10, num_beam_values*sizeof(user_precision_complex_t ) ) );
  cudaErrorCheckCall( cudaMalloc( (void**)&d_primay_beam_J11, num_beam_values*sizeof(user_precision_complex_t ) ) );

  run_and_map_multifreq_calc_CUDA_FEE_beam(beam_settings,
                                          azs, zas,
                                          sin_para_angs, cos_para_angs,
                                          num_components, num_freqs, num_times,
                                          rotation, scaling,
                                          (cuUserComplex *)d_primay_beam_J00,
                                          (cuUserComplex *)d_primay_beam_J01,
                                          (cuUserComplex *)d_primay_beam_J10,
                                          (cuUserComplex *)d_primay_beam_J11);

  cudaErrorCheckCall( cudaMemcpy(primay_beam_J00, d_primay_beam_J00,
                num_beam_values*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J01, d_primay_beam_J01,
                num_beam_values*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J10, d_primay_beam_J10,
                num_beam_values*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost) );
  cudaErrorCheckCall( cudaMemcpy(primay_beam_J11, d_primay_beam_J11,
                num_beam_values*sizeof(user_precision_complex_t),cudaMemcpyDeviceToHost) );

  cudaErrorCheckCall( cudaFree(d_primay_beam_J00) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J01) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J10) );
  cudaErrorCheckCall( cudaFree(d_primay_beam_J11) );

}
