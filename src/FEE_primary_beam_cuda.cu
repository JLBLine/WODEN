#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuComplex.h>
#include <complex.h>
#include <assert.h>
#include <stdio.h>

#include "cudacomplex.h"
#include "legendre_polynomial.h"
#include "cudacheck.h"
#include "FEE_primary_beam_cuda.h"
#include "rts_cube.h"


// #include <iostream>
// using namespace std;

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

__global__ void kern_rotate_FEE_beam(cuFloatComplex *d_FEE_beam_gain_matrices,
                                float *d_sin_para_angs, float *d_cos_para_angs,
                                int num_components, int num_time_steps) {
  const int iCoord = threadIdx.x + (blockDim.x*blockIdx.x);

  if (iCoord < num_components*num_time_steps) {

    int time_ind = (int)floorf((float)iCoord / (float)num_components);

    cuFloatComplex prerot0 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 0];
    cuFloatComplex prerot1 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 1];
    cuFloatComplex prerot2 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 2];
    cuFloatComplex prerot3 = d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 3];

    float sinrot = d_sin_para_angs[iCoord];
    float cosrot = d_cos_para_angs[iCoord];

    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 3] = prerot3*cosrot + prerot2*sinrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 2] = -prerot3*sinrot + prerot2*cosrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 1] = prerot1*cosrot + prerot0*sinrot;
    d_FEE_beam_gain_matrices[iCoord*MAX_POLS + 0] = -prerot1*sinrot + prerot0*cosrot;

  }
}

extern "C" void calc_CUDA_FEE_beam(float *d_beam_reals, float *d_beam_imags,
                                   float *azs, float *zas,
                                   float *sin_para_angs, float *cos_para_angs,
                                   int num_components, int num_time_steps,
                                   copy_primary_beam_t *FEE_beam) {

  // float _Complex *d_FEE_beam_gain_matrices = NULL;

  printf("\tDoing FEE beam tings\n");

  cudaMalloc( (void **)&FEE_beam->d_FEE_beam_gain_matrices, num_time_steps*num_components*MAX_POLS*sizeof(float _Complex));

  // printf("calc_CUDA_FEE_beam 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  //Ensure gains are zero before summing results to them
  float _Complex *zero_array = NULL;
  // cudaMallocHost( (void**)&zero_array, num_time_steps*num_components*MAX_POLS*sizeof(float _Complex) );
  zero_array = (float _Complex *)malloc( num_time_steps*num_components*MAX_POLS*sizeof(float _Complex) );

  for (int i = 0; i < num_time_steps*num_components*MAX_POLS; i++) {
    zero_array[i] = {0.0, 0.0};
  }
  // printf("calc_CUDA_FEE_beam 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );
  //
  cudaMemcpy(FEE_beam->d_FEE_beam_gain_matrices, zero_array, num_time_steps*num_components*MAX_POLS*sizeof(float _Complex), cudaMemcpyHostToDevice );
  // cudaFreeHost(zero_array);
  free(zero_array);

  // float az = 0.0;
  // float za = 0.0;
  float rotation = 0.0;
  float scaling = 1.0;

  // printf("Location 3\n");


  RTS_CUDA_get_TileGains(azs, zas,
    sin_para_angs, cos_para_angs,
    num_time_steps, num_components, rotation, FEE_beam,
    FEE_beam->d_FEE_beam_gain_matrices, scaling);

  // printf("calc_CUDA_FEE_beam 3: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  cudaFree( FEE_beam->emn_T );
  cudaFree( FEE_beam->emn_P );
  cudaFree( FEE_beam->d_emn_T_sum );
  cudaFree( FEE_beam->d_emn_P_sum );

  cudaFree(FEE_beam->rts_P_sin);
  cudaFree(FEE_beam->rts_P1);

}

extern "C" void free_FEE_primary_beam_from_GPU(copy_primary_beam_t *primary_beam){
  cudaFree(primary_beam->d_M);
  cudaFree(primary_beam->d_N);
  cudaFree(primary_beam->d_Q1);
  cudaFree(primary_beam->d_Q2);

  // cudaFree(primary_beam->cos_para_angs);
  // cudaFree(primary_beam->sin_para_angs);


}



extern "C" void copy_FEE_primary_beam_to_GPU(beam_settings_t beam_settings,
                                             int num_time_steps){

  const int n_pols=2; // instrumental pols
  const int nMN = beam_settings.FEE_beam->nMN; // max_length HDFBeamInit
  int arrSize = n_pols * nMN;

  int nStations = 1;

  cudaMalloc( (void**)&beam_settings.FEE_beam->d_M, arrSize*sizeof(float) );
  cudaMalloc( (void**)&beam_settings.FEE_beam->d_N, arrSize*sizeof(float) );
  cudaMalloc( (void**)&beam_settings.FEE_beam->d_Q1, arrSize*sizeof(float _Complex) );
  cudaMalloc( (void**)&beam_settings.FEE_beam->d_Q2, arrSize*sizeof(float _Complex) );

  float *h_params;

  // cudaMallocHost( (void**)&h_params, arrSize*sizeof(float) );
  h_params = (float *)malloc(arrSize*sizeof(float) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (float)(beam_settings.FEE_beam->M[i][j]);
    }
  }

  cudaMemcpy(beam_settings.FEE_beam->d_M, h_params, arrSize*sizeof(float), cudaMemcpyHostToDevice );

  for(unsigned int i=0; i<n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (float)(beam_settings.FEE_beam->N[i][j]);
    }
  }

  cudaMemcpy( beam_settings.FEE_beam->d_N, h_params, arrSize*nStations*sizeof(float), cudaMemcpyHostToDevice );
  // cudaFreeHost( h_params );
  free(h_params);

  float _Complex *h_Qdata;
  // cudaMallocHost( (void**)&h_Qdata, arrSize*nStations*sizeof(float _Complex) );
  h_Qdata = (float _Complex *)malloc(arrSize*sizeof(float _Complex) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = beam_settings.FEE_beam->Q1[i][j];
    }
  }

  cudaMemcpy(beam_settings.FEE_beam->d_Q1, h_Qdata, arrSize*nStations*sizeof(float _Complex), cudaMemcpyHostToDevice );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = beam_settings.FEE_beam->Q2[i][j];
    }
  }

  cudaMemcpy(beam_settings.FEE_beam->d_Q2, h_Qdata ,arrSize*nStations*sizeof(float _Complex) , cudaMemcpyHostToDevice );
  // cudaFreeHost( h_Qdata );
  free( h_Qdata );


  // cudaMalloc( (void**)&beam_settings.FEE_beam->cos_para_angs, num_time_steps*sizeof(float) );
  // cudaMalloc( (void**)&beam_settings.FEE_beam->sin_para_angs, num_time_steps*sizeof(float) );
  //
  //
  // cudaMemcpy( beam_settings.FEE_beam->cos_para_angs, beam_settings.para_cosrot, num_time_steps*sizeof(float), cudaMemcpyHostToDevice );
  // cudaMemcpy( beam_settings.FEE_beam->sin_para_angs, beam_settings.para_sinrot, num_time_steps*sizeof(float), cudaMemcpyHostToDevice );

}

extern "C" void calc_FEE_beam(float *az, float *za, int num_azza,
           float *sin_para_angs, float *cos_para_angs,
           int num_time_steps, copy_primary_beam_t *primary_beam,
           float _Complex *h_FEE_beam_gains,
           int scaling) {
  // printf("We are beginning the calc_FEE_beam\n");

  const int n_pols=2; // instrumental pols
  const int nMN = primary_beam->nMN; // max_length HDFBeamInit
  int arrSize = n_pols * nMN;
  int nStations = 1;

  cudaMalloc( (void**)&primary_beam->d_M, arrSize*sizeof(float) );
  cudaMalloc( (void**)&primary_beam->d_N, arrSize*sizeof(float) );
  cudaMalloc( (void**)&primary_beam->d_Q1, arrSize*sizeof(float _Complex) );
  cudaMalloc( (void**)&primary_beam->d_Q2, arrSize*sizeof(float _Complex) );

  // printf("CUDA error 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  float *h_params=NULL;

  // cudaMallocHost( (void**)&h_params, arrSize*sizeof(float) );
  h_params = (float *)malloc(arrSize*sizeof(float));

  // printf("CUDA error 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (float)(primary_beam->M[i][j]);
    }
  }

  // printf("CUDA error 3: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  cudaMemcpy(primary_beam->d_M, h_params, arrSize*sizeof(float), cudaMemcpyHostToDevice );

  for(unsigned int i=0; i<n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      h_params[j+(i*nMN)] = (float)(primary_beam->N[i][j]);
    }
  }

  cudaMemcpy( primary_beam->d_N, h_params, arrSize*sizeof(float), cudaMemcpyHostToDevice );
  // cudaFreeHost( h_params );
  free( h_params );

  float _Complex *h_Qdata=NULL;
  // cudaMallocHost( (void**)&h_Qdata, arrSize*nStations*sizeof(float _Complex) );
  h_Qdata = (float _Complex *)malloc(arrSize*sizeof(float _Complex));

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = primary_beam->Q1[i][j];
    }
  }

  cudaMemcpy(primary_beam->d_Q1, h_Qdata, arrSize*nStations*sizeof(float _Complex), cudaMemcpyHostToDevice );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
        h_Qdata[j+(i*nMN)] = primary_beam->Q2[i][j];
    }
  }

  cudaMemcpy(primary_beam->d_Q2, h_Qdata ,arrSize*nStations*sizeof(float _Complex) , cudaMemcpyHostToDevice );
  // cudaFreeHost( h_Qdata );
  free( h_Qdata );

  float _Complex *TileGainMatrices = NULL;
  cudaMalloc( (void **)&TileGainMatrices, num_azza*MAX_POLS*sizeof(float _Complex));

  //Ensure gains are zero before summing results to them
  for (size_t i = 0; i < num_azza*MAX_POLS; i++) {
    h_FEE_beam_gains[i] = {0.0,0.0};
  }
  cudaMemcpy(TileGainMatrices, h_FEE_beam_gains, num_azza*MAX_POLS*sizeof(float _Complex), cudaMemcpyHostToDevice );

  float rotation = 0.0;
  int nmax = primary_beam->nmax;

  RTS_CUDA_get_TileGains(az, za, sin_para_angs, cos_para_angs,
    num_time_steps, num_azza, rotation, primary_beam, TileGainMatrices, scaling);

  // printf("CUDA error 4: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  cudaMemcpy(h_FEE_beam_gains, TileGainMatrices, num_azza*MAX_POLS*sizeof(float _Complex), cudaMemcpyDeviceToHost );

  cudaFree( primary_beam->emn_T );
  cudaFree( primary_beam->emn_P );
  cudaFree( primary_beam->d_emn_T_sum );
  cudaFree( primary_beam->d_emn_P_sum );

  cudaFree(primary_beam->rts_P_sin);
  cudaFree(primary_beam->rts_P1);

  cudaFree( TileGainMatrices );

}



extern "C" void get_HDFBeam_normalisation(beam_settings_t beam_settings) {

  //The FEE beam is in theta-phi polarisations, which means we need a different
  //norm factor in north, east, south, west
  float norm_azs[4] = {M_PI/2.0, 0.0, M_PI, M_PI/2.0};
  float norm_zas[4] = {0.0, 0.0, 0.0, 0.0};

  //When running the FEE beam, we later rotate by parrallactic angle. We
  //will normalise before this rotation, so only need to calculate the
  //beam values at zenith once, and the beam is stationary, so only need
  //to calculate for one time step here
  int num_time_steps = 1;



  float _Complex *all_norm_gains=NULL;
  all_norm_gains = (float _Complex *)malloc(num_time_steps*MAX_POLS*MAX_POLS*sizeof(float _Complex));
  printf("Before calc_FEE_beam: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  //Scaling means don't apply normalistaion, as we are calculating it here
  int scaling = 0;

  calc_FEE_beam(norm_azs, norm_zas,  MAX_POLS*num_time_steps,
    beam_settings.para_sinrot, beam_settings.para_cosrot,
    num_time_steps, beam_settings.FEE_beam_zenith, all_norm_gains, scaling);

  printf("After calc_FEE_beam: %s\n", cudaGetErrorString( cudaGetLastError() ) );

  for (size_t i = 0; i < MAX_POLS; i++) {

    // float real_norm = creal(all_norm_gains[i*MAX_POLS + MAX_POLS - 1 - i]);
    // float imag_norm = cimag(all_norm_gains[i*MAX_POLS + MAX_POLS - 1 - i]);
    // float abs_norm_value = sqrtf(real_norm*real_norm + imag_norm*imag_norm);
    // printf("NORM BEAM STOOF %.6f %.6f %.6f\n",real_norm, imag_norm, abs_norm_value );

    //Do a polarisation reordering here, to follow the RTS code
    beam_settings.FEE_beam->norm_fac[i] = all_norm_gains[i*MAX_POLS + MAX_POLS - 1 - i];

    // beam_settings.FEE_beam->norm_fac[i] = {abs(real_norm),abs(imag_norm)};
  }
  free(all_norm_gains);
}

extern "C" void RTS_CUDA_get_TileGains(float *phi, float *theta,
           float *sin_para_angs, float *cos_para_angs,
           int num_time_steps, int num_components,
           float rotation, copy_primary_beam_t *primary_beam,
           float _Complex *TileGainMatrices, int scaling){

  int n_pols = 2;
  int num_coords = num_time_steps*num_components;

  for (size_t phi_ind = 0; phi_ind < num_coords; phi_ind++) {
    phi[phi_ind] = (M_PI/2.0) - phi[phi_ind];
    if(phi[phi_ind] < 0){
      phi[phi_ind] += 2.0*M_PI;
    }
    // printf("Input phi theta %d %.5f %.5f\n",(int)phi_ind, phi[phi_ind], theta[phi_ind] );
  }

  // printf("Loation 6\n");

  float *d_phi=NULL;
  cudaMalloc( (void**)&d_phi, num_coords*sizeof(float) );
  cudaMemcpy(d_phi, phi, num_coords*sizeof(float), cudaMemcpyHostToDevice );

  float *d_theta=NULL;
  cudaMalloc( (void**)&d_theta, num_coords*sizeof(float) );
  cudaMemcpy(d_theta, theta, num_coords*sizeof(float), cudaMemcpyHostToDevice );

  int nmax = primary_beam->nmax; // Assumes all beams have same dimensions
  int nMN = primary_beam->nMN;

  int P_size = sizeof(float _Complex)*(nmax*nmax + 2*nmax)*num_coords;

  cudaMalloc( (void **)&primary_beam->rts_P_sin, P_size);
  cudaMalloc( (void **)&primary_beam->rts_P1, P_size);

  dim3 grid, threads;

  threads.x = kP1BlockSize;
  grid.x = (int)ceil((float)(nmax*nmax + 2*nmax) / (float)threads.x);

  threads.y = 2;
  grid.y = (int)ceil((float)num_coords / (float)threads.y);

  // printf("WE BE HERE INNIT %d %d\n",grid.x,grid.y);

  RTS_P1SINfKernel<<<grid, threads >>>( d_theta, (cuFloatComplex*)primary_beam->rts_P_sin,
                                      (cuFloatComplex*)primary_beam->rts_P1, nmax, num_coords);

  cudaMalloc( (void **)&primary_beam->emn_P, num_coords * nMN * n_pols * sizeof(float _Complex));
  cudaMalloc( (void **)&primary_beam->emn_T, num_coords * nMN * n_pols * sizeof(float _Complex));
  //
  //
  threads.x = kP1BlockSize;
  threads.y = 2;
  threads.z = 1;
  //
  grid.x = (int)ceil((float)(nMN) / (float)threads.x);
  grid.y = (int)ceil((float)num_coords / (float)threads.y);
  grid.z = n_pols;

  RTS_getTileGainsKernel<<<grid, threads >>>(d_phi, d_theta, nMN, num_coords,
                                primary_beam->d_M,primary_beam->d_N,
                                (cuFloatComplex*)primary_beam->d_Q1, (cuFloatComplex*)primary_beam->d_Q2,
                                (cuFloatComplex*)primary_beam->rts_P_sin, (cuFloatComplex*)primary_beam->rts_P1,
                                (cuFloatComplex*)primary_beam->emn_T, (cuFloatComplex*)primary_beam->emn_P);

  float *d_emn_P_sum_real = NULL;
  float *d_emn_T_sum_real = NULL;
  float *d_emn_P_sum_imag = NULL;
  float *d_emn_T_sum_imag = NULL;

  cudaMalloc( (void**)&d_emn_P_sum_real, num_coords*(2*nmax + 1)*n_pols*sizeof(float) );
  cudaMalloc( (void**)&d_emn_T_sum_real, num_coords*(2*nmax + 1)*n_pols*sizeof(float) );
  cudaMalloc( (void**)&d_emn_P_sum_imag, num_coords*(2*nmax + 1)*n_pols*sizeof(float) );
  cudaMalloc( (void**)&d_emn_T_sum_imag, num_coords*(2*nmax + 1)*n_pols*sizeof(float) );

  //Ensure gains are zero before summing results to them
  float *zero_array = NULL;
  zero_array = (float *)malloc(num_coords*(2*nmax + 1)*n_pols*sizeof(float) );

  for (int i = 0; i < num_coords*(2*nmax + 1)*n_pols; i++) {
    zero_array[i] = 0.0;
  }

  cudaMemcpy(d_emn_P_sum_real, zero_array, num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_emn_P_sum_imag, zero_array, num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_emn_T_sum_real, zero_array, num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_emn_T_sum_imag, zero_array, num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice );


  float *d_m_range = NULL;
  cudaMalloc( (void**)&d_m_range, (2*nmax + 1)*sizeof(float) );
  cudaMemcpy(d_m_range, primary_beam->m_range, (2*nmax + 1)*sizeof(float), cudaMemcpyHostToDevice );

  threads.x = 16;
  threads.y = 8;
  threads.z = 8;

  grid.x = (int)ceil( (float)num_coords / threads.x);
  grid.y = (int)ceil( (2.0*(float)primary_beam->nMN) / threads.y);
  grid.z = (int)ceil( ((float)nmax*2.0 + 1.0) / threads.z);

  kern_sum_emn_PT_by_M<<< grid, threads >>>((cuFloatComplex*)primary_beam->emn_T,
                      (cuFloatComplex*)primary_beam->emn_P,
                      d_emn_T_sum_real, d_emn_T_sum_imag,
                      d_emn_P_sum_real, d_emn_P_sum_imag,
                      d_m_range, primary_beam->d_M, primary_beam->nMN, nmax, num_coords);

  threads.x = 64;
  threads.y = n_pols;
  threads.z = 1;

  grid.x = (int)ceil( (float)num_coords / threads.x);
  grid.y = grid.z = 1;

  cudaMalloc( (void**)&primary_beam->d_emn_T_sum, num_coords*n_pols*(2*nmax + 1)*sizeof(float _Complex) );
  cudaMalloc( (void**)&primary_beam->d_emn_P_sum, num_coords*n_pols*(2*nmax + 1)*sizeof(float _Complex) );

  cudaMemcpy(primary_beam->d_emn_T_sum, zero_array, num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(primary_beam->d_emn_P_sum, zero_array, num_coords*(2*nmax + 1)*n_pols*sizeof(float), cudaMemcpyHostToDevice );

  free(zero_array);

  kern_calc_sigmaTP<<< grid, threads >>>((cuFloatComplex*)TileGainMatrices,
                  d_emn_T_sum_real, d_emn_T_sum_imag,
                  d_emn_P_sum_real, d_emn_P_sum_imag,
                  (cuFloatComplex*)primary_beam->d_emn_T_sum, (cuFloatComplex*)primary_beam->d_emn_P_sum,
                  nmax, num_coords );

  if (scaling == 1.0) {

    threads.x = 128;
    grid.x = (int)ceil( (float)num_coords / threads.x);

    threads.y = threads.z = 1;
    grid.y = grid.z = 1;

    float _Complex *d_norm_fac=NULL;
    cudaMalloc( (void**)&d_norm_fac, MAX_POLS*sizeof(float _Complex) );
    cudaMemcpy(d_norm_fac, primary_beam->norm_fac, MAX_POLS*sizeof(float _Complex), cudaMemcpyHostToDevice );

    kern_apply_FEE_norm<<< grid, threads >>>((cuFloatComplex*)TileGainMatrices,
               (cuFloatComplex*)d_norm_fac, num_coords );

    cudaFree(d_norm_fac);

    float *d_cos_para_angs=NULL;
    float *d_sin_para_angs=NULL;

    cudaMalloc( (void**)&d_cos_para_angs, num_components*num_time_steps*sizeof(float) );
    cudaMalloc( (void**)&d_sin_para_angs, num_components*num_time_steps*sizeof(float) );


    cudaMemcpy( d_cos_para_angs, cos_para_angs, num_components*num_time_steps*sizeof(float), cudaMemcpyHostToDevice );
    cudaMemcpy( d_sin_para_angs, sin_para_angs, num_components*num_time_steps*sizeof(float), cudaMemcpyHostToDevice );

    kern_rotate_FEE_beam<<< grid, threads >>>( (cuFloatComplex*)TileGainMatrices,
                                  d_sin_para_angs, d_cos_para_angs,
                                  num_components, num_time_steps);

    cudaFree(d_cos_para_angs);
    cudaFree(d_sin_para_angs);

  }

  cudaFree(d_emn_T_sum_real);
  cudaFree(d_emn_T_sum_imag);
  cudaFree(d_emn_P_sum_real);
  cudaFree(d_emn_P_sum_imag);

  cudaFree(d_m_range);

  cudaFree( d_phi );
  cudaFree( d_theta);

  cudaFree( primary_beam->emn_T );
  cudaFree( primary_beam->emn_P );
  cudaFree( primary_beam->d_emn_T_sum );
  cudaFree( primary_beam->d_emn_P_sum );

  cudaFree( primary_beam->rts_P_sin );
  cudaFree( primary_beam->rts_P1 );

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

__global__ void kern_calc_sigmaTP(cuFloatComplex *TileGainMatrices,
                float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
                float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
                cuFloatComplex *d_emn_T_sum, cuFloatComplex *d_emn_P_sum,
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

      d_emn_T_sum[array_loc] = d_emn_T;
      d_emn_P_sum[array_loc] = d_emn_P;

      TileGainMatrices[iCoord*MAX_POLS + iPol*n_pols] += d_emn_T;
      TileGainMatrices[iCoord*MAX_POLS + iPol*n_pols + 1] += d_emn_P;

    }
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


//
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

  int iCoord = threadIdx.y + (blockDim.y*blockIdx.y);
  float theta = d_theta[iCoord];

  int i = threadIdx.x + (blockDim.x*blockIdx.x);
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


  int i = threadIdx.x + (blockDim.x*blockIdx.x);
  int iCoord = threadIdx.y + (blockDim.y*blockIdx.y);
  int pol = blockIdx.z;
  int n_pols = gridDim.z; // 2

  // int iStation = blockIdx.z;
  int index;
  float factor1,factor2;
  float C_MN, M_absM;

  int beam_index = i+(pol*nMN);

  if(iCoord < num_coords && i < nMN) {
    float phi = d_phi[iCoord];
    float theta = d_theta[iCoord];

    index = (int) (pb_N[beam_index] - fabs(pb_M[beam_index]));
    if(index >=80){
      printf("Maximum factorial exceeded %d\n", __FUNCTION__,index );
    }
    factor1 = ffactorials[index];
    index = (int) (pb_N[beam_index] + fabs(pb_M[beam_index]));
    if(index >=80){
      printf("Maximum factorial exceeded %d\n", __FUNCTION__,index );
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

    phi_comp = CUBE_DEVICE_CALL(U1polar,(pb_M[beam_index] * phi)) * C_MN * M_absM / sqrt(pb_N[beam_index] * (pb_N[beam_index] + 1.0));

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
