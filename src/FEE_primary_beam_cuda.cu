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

__global__ void kern_print_summin(cuFloatComplex *d_here_is_an_array,
           int array_length, cuFloatComplex *d_output_array) {

  // Start by computing which baseline we're going to do
  const int ind = threadIdx.x + (blockDim.x*blockIdx.x);

  if (ind < array_length) {
    d_output_array[ind] = cuCmulf(d_here_is_an_array[ind],d_here_is_an_array[ind]);
    // printf("%d\n",ind );

  }
  // printf("%d %f \n", ind, (d_Q1[0].x) );

}

extern "C" void test_FEE_beam(float *az, float *za, copy_primary_beam_t *primary_beam,
           float _Complex *copy_hostGains,
           float _Complex *h_rts_P_sin, float _Complex *h_rts_P1,
           float _Complex *h_emn_T, float _Complex *h_emn_P,
           float _Complex *h_emn_T_sum, float _Complex *h_emn_P_sum) {

  const int n_pols=2; // instrumental pols
  const int nMN = primary_beam->nMN; // max_length HDFBeamInit
  int arrSize = n_pols * nMN;

  // az = 0.17453;
  // za = 0.17453;

  float single_az = az[0];
  float single_za = za[0];

  // printf("INPUT GPU3 %.5f %.5f\n", single_az, single_za );
  // printf("INPUT GPU4 %.5f %.5f\n", az, za );

  int nStations = 1;

  cudaMalloc( (void**)&primary_beam->d_M, nStations*arrSize*sizeof(float) );
  cudaMalloc( (void**)&primary_beam->d_N, nStations*arrSize*sizeof(float) );
  cudaMalloc( (void**)&primary_beam->d_Q1, nStations*arrSize*sizeof(float _Complex) );
  cudaMalloc( (void**)&primary_beam->d_Q2, nStations*arrSize*sizeof(float _Complex) );

  float *h_params;

  cudaMallocHost( (void**)&h_params, arrSize*nStations*sizeof(float) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      for(unsigned int n=0; n< nStations; n++){
      h_params[j+(i*nMN)+(n*nMN*n_pols)] = (float)(primary_beam->M[i][j]);
      }
    }
  }

  cudaMemcpy(primary_beam->d_M, h_params, arrSize*nStations*sizeof(float), cudaMemcpyHostToDevice );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      for(unsigned int n=0; n< nStations; n++){
      h_params[j+(i*nMN)+(n*nMN*n_pols)] = (float)(primary_beam->N[i][j]);
      }
    }
  }

  cudaMemcpy( primary_beam->d_N, h_params, arrSize*nStations*sizeof(float), cudaMemcpyHostToDevice );
  cudaFreeHost( h_params );

  float _Complex *h_Qdata;

  cudaMallocHost( (void**)&h_Qdata, arrSize*nStations*sizeof(float _Complex) );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      for(unsigned int n=0; n< nStations; n++){
        h_Qdata[j+(i*nMN)+(n*nMN*n_pols)] = primary_beam->Q1[i][j];
      }
    }
  }

  cudaMemcpy(primary_beam->d_Q1, h_Qdata, arrSize*nStations*sizeof(float _Complex), cudaMemcpyHostToDevice );

  for(unsigned int i=0; i< n_pols; i++){
    for(unsigned int j=0; j<nMN; j++){
      for(unsigned int n=0; n< nStations; n++){
        h_Qdata[j+(i*nMN)+(n*nMN*n_pols)] = primary_beam->Q2[i][j];
      }
    }
  }

  cudaMemcpy(primary_beam->d_Q2, h_Qdata ,arrSize*nStations*sizeof(float _Complex) , cudaMemcpyHostToDevice );
  cudaFreeHost( h_Qdata );

  float _Complex *TileGainMatrices=NULL;
  // cudaMalloc( (void **)TileGainMatrices, nStations*MAX_POLS*sizeof(float _Complex));
  cudaMalloc( (void **)&TileGainMatrices, 4*sizeof(float _Complex));

  // printf("FINISHED %d\n",grid.x );



  // float az = 0.0;
  // float za = 0.0;
  float rotation = 0.0;
  float scaling = 1.0;

  int nmax = primary_beam->nmax;

  // printf("INPUT GPU2 %.5f %.5f\n",single_az,single_za );
  RTS_CUDA_get_TileGains(single_az, single_za, rotation, primary_beam, nStations, TileGainMatrices, scaling);
  //

  cudaMemcpy(h_rts_P1, primary_beam->rts_P1, (nmax*nmax + 2*nmax)*sizeof(float _Complex), cudaMemcpyDeviceToHost );
  cudaMemcpy(h_rts_P_sin, primary_beam->rts_P_sin, (nmax*nmax + 2*nmax)*sizeof(float _Complex), cudaMemcpyDeviceToHost );

  cudaMemcpy(h_emn_T, primary_beam->emn_T, 2*primary_beam->nMN*sizeof(float _Complex), cudaMemcpyDeviceToHost );
  cudaMemcpy(h_emn_P, primary_beam->emn_P, 2*primary_beam->nMN*sizeof(float _Complex), cudaMemcpyDeviceToHost );

  cudaMemcpy(h_emn_T_sum, primary_beam->d_emn_T_sum, n_pols*(2*nmax + 1)*sizeof(float _Complex), cudaMemcpyDeviceToHost );
  cudaMemcpy(h_emn_P_sum, primary_beam->d_emn_P_sum, n_pols*(2*nmax + 1)*sizeof(float _Complex), cudaMemcpyDeviceToHost );

  cudaMemcpy(copy_hostGains, TileGainMatrices, nStations*MAX_POLS*sizeof(float _Complex), cudaMemcpyDeviceToHost );


  cudaFree( primary_beam->d_Q1 );
  cudaFree( primary_beam->d_Q2 );
  // cudaFree( d_norm_fac );
  cudaFree( primary_beam->d_N );
  cudaFree( primary_beam->d_M );
  cudaFree( primary_beam->emn_T );
  cudaFree( primary_beam->emn_P );
  cudaFree( primary_beam->d_emn_T_sum );
  cudaFree( primary_beam->d_emn_P_sum );

  cudaFree(primary_beam->rts_P_sin);
  cudaFree(primary_beam->rts_P1);

  cudaFree( TileGainMatrices );


}

extern "C" void RTS_CUDA_get_TileGains(float phi, float theta, float rotation, copy_primary_beam_t *primary_beam, int nStations, float _Complex *TileGainMatrices, int scaling){

  int n_pols = 2;

  phi = (M_PI/2.0 - phi);
  if(phi < 0) phi += 2.0*M_PI;

  copy_primary_beam_t *pb;

  pb = primary_beam;

  int nmax = pb->nmax; // Assumes all beams have same dimensions
  int nMN = pb->nMN;

  //CUDA_SAFE_CALL( cudaMalloc( (void **)&TileGainMatrices, nStations*MAX_POLS*sizeof(float _Complex)));

  // Call P1SIN. This only needs to happen once per last_theta
  // IF nMN == nmax*nmax + 2*nmax then this could be part of the
  // same kernel, but this is somewhat obscure

  // float _Complex *rts_P_sin=NULL;
  // float _Complex *rts_P1=NULL;

  int P_size = sizeof(float _Complex)*(nmax*nmax + 2*nmax);
  // printf("WE GOOD WE DID THIS %d\n",nmax*nmax + 2*nmax);
  CUDA_SAFE_CALL( cudaMalloc( (void **)&pb->rts_P_sin, P_size));
  CUDA_SAFE_CALL( cudaMalloc( (void **)&pb->rts_P1, P_size));

  // dim3 blockDim;
  // blockDim.x = kP1BlockSize;
  // blockDim.y = blockDim.z = 1;
  //
  // dim3 gridDim;
  // gridDim.x = (int)ceil((float)(nmax*nmax + 2*nmax) / (float)(kP1BlockSize));
  // gridDim.y = gridDim.z = 1;

  dim3 grid, threads;

  // CUBE_KERNEL_CALL(RTS_P1SINfKernel, gridDim, blockDim, 0, theta, (cuFloatComplex*)rts_P_sin, (cuFloatComplex*)rts_P1, nmax );
  threads.x = kP1BlockSize;
  grid.x = (int)ceil((float)(nmax*nmax + 2*nmax) / (float)(kP1BlockSize));

  // printf("%d %d\n",threads.x,grid.x );

  printf("GPU INPUT %.5f\n",theta );

  RTS_P1SINfKernel<<<grid, threads >>>( theta, (cuFloatComplex*)pb->rts_P_sin, (cuFloatComplex*)pb->rts_P1, nmax );


  // CUDA_CHECK_KERNEL( "P1SINKernel execution failed\n" );

  // Now we
  // Now call main get_FF kernal. In CPU this is a Loop
  // over 2 instrumental pols and nMN (modes)
  // Here dimensions are nMN by npols by nStations

  // blockDim.x = kP1BlockSize;
  // blockDim.y = blockDim.z = 1;
  //
  //
  // gridDim.x = (int)ceil((float)(nMN) / (float)(kP1BlockSize));
  // gridDim.y = n_pols;
  // gridDim.z = nStations;


  //


  float _Complex *Sigma_P;
  float _Complex *Sigma_T;
  //
  CUDA_SAFE_CALL( cudaMalloc( (void **)&Sigma_P, nMN * n_pols * nStations * sizeof(float _Complex)));
  CUDA_SAFE_CALL( cudaMalloc( (void **)&Sigma_T, nMN * n_pols * nStations * sizeof(float _Complex)));

  CUDA_SAFE_CALL( cudaMalloc( (void **)&primary_beam->emn_P, nMN * n_pols * nStations * sizeof(float _Complex)));
  CUDA_SAFE_CALL( cudaMalloc( (void **)&primary_beam->emn_T, nMN * n_pols * nStations * sizeof(float _Complex)));


  threads.x = kP1BlockSize;
  threads.y = threads.z = 1;

  grid.x = (int)ceil((float)(nMN) / (float)(kP1BlockSize));
  grid.y = n_pols;
  grid.z = nStations;

  RTS_getTileGainsKernel<<<grid, threads >>>(phi, theta, nMN,
                                primary_beam->d_M,primary_beam->d_N,
                                (cuFloatComplex*)primary_beam->d_Q1, (cuFloatComplex*)primary_beam->d_Q2,
                                (cuFloatComplex*)pb->rts_P_sin, (cuFloatComplex*)pb->rts_P1,
                                (cuFloatComplex*)Sigma_T, (cuFloatComplex*)Sigma_P,
                                (cuFloatComplex*)primary_beam->emn_T, (cuFloatComplex*)primary_beam->emn_P);
  //
  // // CUDA_CHECK_KERNEL( "RTS_getTileGainsKernel execution failed\n" );
  //
  // threads.x = nStations;
  // threads.y = threads.z = 1;
  //
  // grid.x = n_pols;
  // grid.y = grid.z = 1;


  float _Complex *d_norm_fac=NULL;
  cudaMalloc( (void**)&d_norm_fac, MAX_POLS*sizeof(float _Complex) );
  cudaMemcpy(d_norm_fac, primary_beam->norm_fac, MAX_POLS*sizeof(float _Complex), cudaMemcpyHostToDevice );

  // printf("NORM FAC AFTER mem copy %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n",
  //     creal(primary_beam->norm_fac[0]), cimag(primary_beam->norm_fac[0]),
  //     creal(primary_beam->norm_fac[1]), cimag(primary_beam->norm_fac[1]),
  //     creal(primary_beam->norm_fac[2]), cimag(primary_beam->norm_fac[2]),
  //     creal(primary_beam->norm_fac[3]), cimag(primary_beam->norm_fac[3]) );


  //
  // RTS_reduce_SigmaPTKernel<<<grid, threads >>>(nMN,
  //                   (cuFloatComplex*)Sigma_T, (cuFloatComplex*)Sigma_P,
  //                   (cuFloatComplex*)TileGainMatrices,
  //                   (cuFloatComplex*)d_norm_fac, scaling,
  //                   rotation);

  float *d_emn_P_sum_real = NULL;
  float *d_emn_T_sum_real = NULL;
  float *d_emn_P_sum_imag = NULL;
  float *d_emn_T_sum_imag = NULL;

  cudaMalloc( (void**)&d_emn_P_sum_real, (2*nmax + 1)*n_pols*sizeof(float) );
  cudaMalloc( (void**)&d_emn_T_sum_real, (2*nmax + 1)*n_pols*sizeof(float) );
  cudaMalloc( (void**)&d_emn_P_sum_imag, (2*nmax + 1)*n_pols*sizeof(float) );
  cudaMalloc( (void**)&d_emn_T_sum_imag, (2*nmax + 1)*n_pols*sizeof(float) );
  // int num_theta = 1;



  float *d_m_range = NULL;
  cudaMalloc( (void**)&d_m_range, (2*nmax + 1)*sizeof(float) );

  float *m_range = NULL;
  cudaMallocHost( (void**)&m_range, (2*nmax + 1)*sizeof(float) );


  float *d_phi_comp = NULL;
  cudaMalloc( (void**)&d_phi_comp, (2*nmax + 1)*sizeof(float _Complex) );

  float _Complex *phi_comp = NULL;
  cudaMallocHost( (void**)&phi_comp, (2*nmax + 1)*sizeof(float _Complex) );

  printf("Managed to do the cudaMalloc\n");

  for (int m = 0; m < 2*nmax + 1; m++) {
    m_range[m] = -(float)nmax + (float)m;

    // float _Complex for_exp = 0 + (float)m*_Complex_I;
    float _Complex for_exp = {0.0, (float)m*phi};
    phi_comp[m] = cexpf(for_exp);
    // printf("The looops %d %.5f %.5f %.5f\n", m, m_range[m], creal(phi_comp[m]), cimag(phi_comp[m]));
  }

  cudaMemcpy(d_m_range, m_range, (2*nmax + 1)*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(d_phi_comp, phi_comp, (2*nmax + 1)*sizeof(float _Complex), cudaMemcpyHostToDevice );

  printf("Managed to do the loop\n");

  threads.x = 32;
  threads.y = 4;
  threads.z = n_pols;

  printf("pb->nM is %d, nmax is %d\n",pb->nMN,nmax );

  grid.x = (int)ceil( (float)pb->nMN / threads.x);
  grid.y = (int)ceil( ((float)nmax*2.0 + 1.0) / threads.y);
  grid.z = 1;

  kern_sum_emn_PT_by_M<<< grid, threads >>>((cuFloatComplex*)primary_beam->emn_T,
                      (cuFloatComplex*)primary_beam->emn_P,
                      d_emn_T_sum_real, d_emn_T_sum_imag,
                      d_emn_P_sum_real, d_emn_P_sum_imag,
                      d_m_range, primary_beam->d_M, pb->nMN, nmax);

  threads.x = 2;
  threads.y = threads.z = 1;

  grid.x = grid.y = grid.z = 1;

  cudaMalloc( (void**)&primary_beam->d_emn_T_sum, n_pols*(2*nmax + 1)*sizeof(float _Complex) );
  cudaMalloc( (void**)&primary_beam->d_emn_P_sum, n_pols*(2*nmax + 1)*sizeof(float _Complex) );

  kern_calc_sigmaTP<<< grid, threads >>>((cuFloatComplex*)TileGainMatrices,
                  (cuFloatComplex*)d_phi_comp,
                  d_emn_T_sum_real, d_emn_T_sum_imag,
                  d_emn_P_sum_real, d_emn_P_sum_imag,
                  (cuFloatComplex*)primary_beam->d_emn_T_sum, (cuFloatComplex*)primary_beam->d_emn_P_sum,
                  nmax );

  // //
  // // // CUDA_SAFE_CALL( cudaFree(Sigma_P));
  // // // CUDA_SAFE_CALL( cudaFree(Sigma_T));


  threads.x = 1;
  threads.y = threads.z = 1;

  grid.x = grid.y = grid.z = 1;

  kern_apply_norm<<< grid, threads >>>((cuFloatComplex*)TileGainMatrices,
             (cuFloatComplex*)d_norm_fac );

  cudaFree(d_emn_T_sum_real);
  cudaFree(d_emn_T_sum_imag);
  cudaFree(d_emn_P_sum_real);
  cudaFree(d_emn_P_sum_imag);

  cudaFree(d_m_range);
  cudaFreeHost( m_range);

  cudaFree(d_phi_comp);
  cudaFreeHost( phi_comp);

  cudaFree(d_norm_fac);
  // //
  cudaFree(Sigma_P);
  cudaFree(Sigma_T);

}

__global__ void kern_apply_norm(cuFloatComplex *TileGainMatrices,
           cuFloatComplex *d_norm_fac ){

  // const int iM_index = threadIdx.x + (blockDim.x*blockIdx.x);

  TileGainMatrices[0] = TileGainMatrices[0] / d_norm_fac[0];
  TileGainMatrices[1] = TileGainMatrices[1] / d_norm_fac[1];
  TileGainMatrices[2] = TileGainMatrices[2] / d_norm_fac[2];
  TileGainMatrices[3] = TileGainMatrices[3] / d_norm_fac[3];

}

__global__ void kern_calc_sigmaTP(cuFloatComplex *TileGainMatrices,
                cuFloatComplex *d_phi_comp,
                float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
                float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
                cuFloatComplex *d_emn_T_sum, cuFloatComplex *d_emn_P_sum,
                int nmax) {

  const int iPol = threadIdx.x + (blockDim.x*blockIdx.x);
  int n_pols=2;

  int numM = 2*nmax + 1;

  for(int i=0; i < numM; i++){
    cuFloatComplex d_emn_T = make_cuComplex(d_emn_T_sum_real[iPol*numM + i],d_emn_T_sum_imag[iPol*numM + i]);
    cuFloatComplex d_emn_P = make_cuComplex(d_emn_P_sum_real[iPol*numM + i],d_emn_P_sum_imag[iPol*numM + i]);
    // cuFloatComplex d_emn_P = d_phi_comp[i];

    // printf("%d %.5f %.5f %.5f %.5f\n",i,d_emn_T_sum_real[iPol*numM + i],d_emn_T_sum_imag[iPol*numM + i],d_emn_P.x,d_emn_P.y );

    d_emn_T_sum[iPol*numM + i] = d_emn_T;
    d_emn_P_sum[iPol*numM + i] = d_emn_P;

    // TileGainMatrices[iPol*n_pols] += cuCmulf(d_emn_T , d_phi_comp[i]);
    // TileGainMatrices[iPol*n_pols + 1] += cuCmulf(d_emn_P , d_phi_comp[i]);

    TileGainMatrices[iPol*n_pols] += d_emn_T;
    TileGainMatrices[iPol*n_pols + 1] += d_emn_P;

  }

}
//
//
//
//   // // order of instrumental pols is reversed between FEKO and RTS
//   const int rtspol = abs(pol-1);
//   // const int rtspol = pol;
//   const int iGainIndex = iStation*MAX_POLS;
//   TileGainMatrices[iGainIndex+2*rtspol] = sum_P;
//   TileGainMatrices[iGainIndex+2*rtspol+1] = sum_T;
//
//   // printf("CHANGES? %d %d %.5f %.5f %.5f %.5f\n",iGainIndex+2*rtspol, iGainIndex+2*rtspol+1, sum_P.x,sum_P.y,sum_T.x,sum_T.y );
//
//   // printf("INDEXES %d %d\n",iGainIndex+2*rtspol, iGainIndex+2*rtspol+1);
//
//   cuFloatComplex to_print;
//   to_print = norm_fac[0];
//
//   //
//
//   if(scaling==1){
//     // printf("Norm fac inside GPU %.5f %.5f\n",to_print.x, to_print.y );
//
//   }
//   //
//   // if(rotation!=0.0){
//   //   float crot, srot;
//   //   crot = cosf(rotation);
//   //   srot = sinf(rotation);
//   //   cuFloatComplex prerot[2];
//   //
//   //   prerot[0] = TileGainMatrices[iGainIndex+2*rtspol];
//   //   prerot[1] = TileGainMatrices[iGainIndex+2*rtspol+1];
//   //
//   //   TileGainMatrices[iGainIndex+2*rtspol+1] = prerot[1]*crot + prerot[0]*srot;
//   //   TileGainMatrices[iGainIndex+2*rtspol] = -prerot[1]*srot + prerot[0]*crot;
//
//   // }
//
// }


__global__ void kern_sum_emn_PT_by_M(cuFloatComplex *emn_T, cuFloatComplex *emn_P,
           float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
           float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
           float *d_m_range, float *d_M, int nMN, int nmax){

  const int iM_index = threadIdx.x + (blockDim.x*blockIdx.x);
  const int iM_value = threadIdx.y + (blockDim.y*blockIdx.y);
  const int iPol = threadIdx.z + (blockDim.z*blockIdx.z);

  if(iM_index < nMN && iM_value < (2*nmax + 1)) {
    if (d_M[iPol*nMN + iM_index] == d_m_range[iM_value]) {
      cuFloatComplex emn_P_value = emn_P[iPol*nMN + iM_index];
      cuFloatComplex emn_T_value = emn_T[iPol*nMN + iM_index];

      // printf("%.8f %.8f %.8f %.8f\n",emn_P_value.x,emn_P_value.y,emn_T_value.x,emn_T_value.y );

      //

      // if (emn_P_value.x > 0) {
      //   printf("%d %d %d %d %.5f %d\n", iM_index, iM_value, iPol*nMN + iM_index,
      //                                   iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax,
      //                                   d_M[iPol*nMN + iM_index], (int)d_m_range[iM_value] );
      // }

      atomicAdd(&d_emn_P_sum_real[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], emn_P_value.x);
      atomicAdd(&d_emn_P_sum_imag[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], emn_P_value.y);

      atomicAdd(&d_emn_T_sum_real[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], emn_T_value.x);
      atomicAdd(&d_emn_T_sum_imag[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], emn_T_value.y);

      // printf("%d %d %d\n",iM_index,iM_value,iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax );
      //
      // atomicAdd(&d_emn_P_sum_real[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], 1.0);
      // atomicAdd(&d_emn_T_sum_real[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], 1.0);
      // atomicAdd(&d_emn_P_sum_imag[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], 1.0);
      // atomicAdd(&d_emn_T_sum_imag[iPol*(2*nmax + 1) + (int)d_m_range[iM_value] + nmax], 1.0);

    }

  }

}

// extern "C" void RTS_CUDA_SendPrimaryBeam_local( copy_cal_context_t *calContext,
//                            const int nStations ) {
//
// // FEE Beam models
//
//   float *h_params;
//
//   const int n_pols=2; // instrumental pols
//   const int nMN = calContext->primary_beams[0].nMN; // max_length HDFBeamInit
//   int arrSize = n_pols * nMN;
//
//   CUDA_SAFE_CALL( cudaMallocHost( (void**)&h_params, arrSize*nStations*sizeof(float) ) );
//
//   for(unsigned int i=0; i< n_pols; i++){
//     for(unsigned int j=0; j<nMN; j++){
//       for(unsigned int n=0; n< nStations; n++){
//       h_params[j+(i*nMN)+(n*nMN*n_pols)] = (float)(calContext->primary_beams[n].M[i][j]);
//       }
//     }
//   }
//
//   CUDA_SAFE_CALL( cudaMemcpy( calContext->d_M, h_params,arrSize*nStations*sizeof(float) ,
//                               cudaMemcpyHostToDevice ) );
//
//   for(unsigned int i=0; i< n_pols; i++){
//     for(unsigned int j=0; j<nMN; j++){
//       for(unsigned int n=0; n< nStations; n++){
//       h_params[j+(i*nMN)+(n*nMN*n_pols)] = (float)(calContext->primary_beams[n].N[i][j]);
//       }
//     }
//   }
//
//   CUDA_SAFE_CALL( cudaMemcpy( calContext->d_N, h_params,arrSize*nStations*sizeof(float) ,
//                               cudaMemcpyHostToDevice ) );
//
//   CUDA_SAFE_CALL( cudaFreeHost( h_params ) );
//
//   float _Complex *h_Qdata;
//
//   CUDA_SAFE_CALL( cudaMallocHost( (void**)&h_Qdata, arrSize*nStations*sizeof(float _Complex) ) );
//
//   for(unsigned int i=0; i< n_pols; i++){
//     for(unsigned int j=0; j<nMN; j++){
//       for(unsigned int n=0; n< nStations; n++){
//       //h_Qdata[j+(i*nMN)+(n*nMN*n_pols)] = {crealf(calContext->primary_beams[n].Q1[i][j]),cimagf(calContext->primary_beams[n].Q1[i][j])};
//         h_Qdata[j+(i*nMN)+(n*nMN*n_pols)] = calContext->primary_beams[n].Q1[i][j];
//       }
//     }
//   }
//
//   CUDA_SAFE_CALL( cudaMemcpy( calContext->d_Q1, h_Qdata,arrSize*nStations*sizeof(float _Complex) ,
//                               cudaMemcpyHostToDevice ) );
//
//   for(unsigned int i=0; i< n_pols; i++){
//     for(unsigned int j=0; j<nMN; j++){
//       for(unsigned int n=0; n< nStations; n++){
//       //h_Qdata[j+(i*nMN)+(n*nMN*n_pols)] = {crealf(calContext->primary_beams[n].Q2[i][j]),cimagf(calContext->primary_beams[n].Q2[i][j])};
//         h_Qdata[j+(i*nMN)+(n*nMN*n_pols)] = calContext->primary_beams[n].Q2[i][j];
//       }
//     }
//   }
//
//   CUDA_SAFE_CALL( cudaMemcpy( calContext->d_Q2, h_Qdata,arrSize*nStations*sizeof(float _Complex) ,
//                               cudaMemcpyHostToDevice ) );
//
//
//   CUDA_SAFE_CALL( cudaFreeHost( h_Qdata ) );
//
//   // Copy zenith normalization factor. The same for all beams
//
//   //float *norm_fac;
//
//
//   //CUDA_SAFE_CALL( cudaMallocHost( (void**)&norm_fac, MAX_POLS*sizeof(float _Complex) ));
//
//   CUDA_SAFE_CALL( cudaMemcpy( calContext->d_norm_fac, calContext->primary_beams[0].norm_fac,MAX_POLS*sizeof(float _Complex) ,
//                               cudaMemcpyHostToDevice ) );
//
//   //CUDA_SAFE_CALL( cudaFreeHost( norm_fac ) );
//
// }
//
// extern "C" void RTS_CUDA_allocate_TileGains(float _Complex **TileGainMatrices, int nStations){
//   CUDA_SAFE_CALL( cudaMalloc( (void **)TileGainMatrices, nStations*MAX_POLS*sizeof(float _Complex)));
// }
//
//
//
// float _Complex *rts_P_sin=NULL;
// float _Complex *rts_P1=NULL;
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

  //return v;

}
//
//
//
__global__ void RTS_P1SINfKernel(float theta, cuFloatComplex *rts_P_sin, cuFloatComplex *rts_p1, int nmax){

  // CUBE_START;
  //int i = threadIdx.x;
  int i = threadIdx.x + (blockDim.x*blockIdx.x);
  int n = floor(sqrt((float)(i+1)));
  int index = i - ((n-1)*(n-1) + 2*(n-1));
  int pm_index;
  float u = cos(theta);
  float sin_th = sin(theta);
  float delu = 1.0e-6;

  if(i >= nmax*nmax + 2*nmax){
    return;
  }

  // printf("THIS THIS THIS %.5f\n",theta );

  if(index < n){
    pm_index = n-index;
  } else {
    pm_index = index-n;
  }

  float p;
  float Pm1;
  float Pm_sin;

  //float *pm_vals;

  //pm_vals = ( float * ) malloc ((n+1) * sizeof ( float ) );

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

  // rts_P_sin[i] = make_float _Complex(Pm_sin,0);
  // rts_p1[i] = make_float _Complex(Pm1,0);

  // printf("dis %d %.1f\n",i, Pm_sin );
  rts_P_sin[i] = make_cuFloatComplex(Pm_sin,0);
  rts_p1[i] = make_cuFloatComplex(Pm1,0);

  // printf("RTS_P1SINfKernel %.3f %.3f\n",Pm_sin,Pm1 );

  // CUBE_END;
}
//
//
__global__ void RTS_getTileGainsKernel( float phi, float theta, int nMN,
           float *pb_M, float *pb_N,
           cuFloatComplex *pb_Q1, cuFloatComplex *pb_Q2,
           cuFloatComplex *rts_P_sin, cuFloatComplex *rts_P1,
           cuFloatComplex *Sigma_T, cuFloatComplex *Sigma_P,
           cuFloatComplex *emn_T, cuFloatComplex *emn_P){


  int i = threadIdx.x + (blockDim.x*blockIdx.x);
  int pol = blockIdx.y;
  // int iStation = blockIdx.z;
  int index;
  float factor1,factor2;
  float C_MN, M_absM;

  int iStation = 0;

  int n_pols = gridDim.y; // 2
  int beam_index = i+(pol*nMN)+(iStation*n_pols*nMN);

  if(i >= nMN){ // Thread out of bounds
    return;
  }

 // CHANGE INDEXING FOR STATIONS!!!!!!!!!!!!!!!!



  index = (int) (pb_N[beam_index] - fabs(pb_M[beam_index]));
  if(index >=80){
    printf("Maximum factorial exceeded\n", __FUNCTION__ );
    //LOG( LOG_ERR,"Maximum factorial exceeded\n", __FUNCTION__ );
    //exit(-1);
  }
  factor1 = ffactorials[index];
  index = (int) (pb_N[beam_index] + fabs(pb_M[beam_index]));
  if(index >=80){
    printf("Maximum factorial exceeded\n", __FUNCTION__ );
    //LOG( LOG_ERR,"Maximum factorial exceeded\n", __FUNCTION__ );
    //exit(-1);
  }
  factor2 = ffactorials[index];

  C_MN = sqrt(0.5 * (2 * pb_N[beam_index] + 1) * factor1 / factor2);

  // printf("%d %d %d %.5f \n",i,pol,iStation,C_MN);

  // if(pb_M[beam_index] == 0){ //diverging??
  //   M_absM = 1;
  // } else {
  //   M_absM = -(pb_M[beam_index] / fabs(pb_M[beam_index]));
  // }
  // M_absM = pow(M_absM,pb_M[beam_index]);

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

  // double real_phi_const;
  // real_phi_const = C_MN * M_absM / sqrt(pb_N[beam_index] * (pb_N[beam_index] + 1.0));
  //
  // cuFloatComplex phi_const = make_cuComplex(real_phi_const, 0.0);

  cuFloatComplex this_emn_T, this_emn_P;
  float T_exp, P_exp;
  float u = cos(theta);

  T_exp = (pb_N[beam_index]);
  P_exp = (pb_N[beam_index])+1.0;

  this_emn_T = make_cuComplex(cos(T_exp * M_PI/2.0),sin(T_exp * M_PI/2.0));
  this_emn_P = make_cuComplex(cos(P_exp * M_PI/2.0),sin(P_exp * M_PI/2.0));

  this_emn_T *= (rts_P_sin[i]*(fabs(pb_M[beam_index])*pb_Q2[beam_index]*u - pb_M[beam_index] * pb_Q1[beam_index]) + pb_Q2[beam_index] * rts_P1[i]);
  this_emn_P *= (rts_P_sin[i]*(pb_M[beam_index]*pb_Q2[beam_index] - fabs(pb_M[beam_index])*pb_Q1[beam_index]*u) - pb_Q1[beam_index] * rts_P1[i]);

  emn_T[beam_index] = cuCmulf(this_emn_T,phi_comp);
  emn_P[beam_index] = cuCmulf(this_emn_P,phi_comp);

  Sigma_P[i+(nMN*pol)+(iStation*n_pols*nMN)] = cuCmulf(this_emn_P,phi_comp);
  Sigma_T[i+(nMN*pol)+(iStation*n_pols*nMN)] = cuCmulf(this_emn_T,phi_comp);

  // emn_T[beam_index] = cuCmulf(this_emn_T,phi_const);
  // emn_P[beam_index] = cuCmulf(this_emn_P,phi_const);

  // Sigma_P[i+(nMN*pol)+(iStation*n_pols*nMN)] = cuCmulf(this_emn_P,phi_const);
  // Sigma_T[i+(nMN*pol)+(iStation*n_pols*nMN)] = cuCmulf(this_emn_T,phi_const);

  // Now have one long Sigma_P for all modes, pols and tiles. Want to reduce (sum) per tile
  // and assign to TileGainMatrices

}
//
//



// __global__ void RTS_reduce_SigmaPTKernel(int nMN, cuFloatComplex *Sigma_T,
//            cuFloatComplex *Sigma_P, cuFloatComplex *TileGainMatrices,
//            cuFloatComplex *norm_fac, int scaling, float rotation) {
//
//   int iStation = threadIdx.x;
//   int pol = blockIdx.x;
//   int n_pols=2;
//
//   cuFloatComplex sum_P;
//   cuFloatComplex sum_T;
//
//   sum_P = sum_T = make_cuComplex(0.0,0.0);
//
//   for(int i=0;i<nMN;i++){
//     sum_P += Sigma_P[i+(pol*nMN)+(iStation*n_pols*nMN)];
//     sum_T += Sigma_T[i+(pol*nMN)+(iStation*n_pols*nMN)];
//   }
//
//
//
//   // // order of instrumental pols is reversed between FEKO and RTS
//   const int rtspol = abs(pol-1);
//   // const int rtspol = pol;
//   const int iGainIndex = iStation*MAX_POLS;
//   TileGainMatrices[iGainIndex+2*rtspol] = sum_P;
//   TileGainMatrices[iGainIndex+2*rtspol+1] = sum_T;
//
//   // printf("CHANGES? %d %d %.5f %.5f %.5f %.5f\n",iGainIndex+2*rtspol, iGainIndex+2*rtspol+1, sum_P.x,sum_P.y,sum_T.x,sum_T.y );
//
//   // printf("INDEXES %d %d\n",iGainIndex+2*rtspol, iGainIndex+2*rtspol+1);
//
//   // cuFloatComplex to_print;
//   // to_print = norm_fac[0];
//
//   //
//
//   if(scaling==1){
//     // printf("Norm fac inside GPU %.5f %.5f\n",to_print.x, to_print.y );
//     TileGainMatrices[iGainIndex+2*rtspol] = TileGainMatrices[iGainIndex+2*rtspol] / norm_fac[2*rtspol];
//     TileGainMatrices[iGainIndex+2*rtspol+1] = TileGainMatrices[iGainIndex+2*rtspol+1] / norm_fac[2*rtspol+1];
//   }
//   //
//   // if(rotation!=0.0){
//   //   float crot, srot;
//   //   crot = cosf(rotation);
//   //   srot = sinf(rotation);
//   //   cuFloatComplex prerot[2];
//   //
//   //   prerot[0] = TileGainMatrices[iGainIndex+2*rtspol];
//   //   prerot[1] = TileGainMatrices[iGainIndex+2*rtspol+1];
//   //
//   //   TileGainMatrices[iGainIndex+2*rtspol+1] = prerot[1]*crot + prerot[0]*srot;
//   //   TileGainMatrices[iGainIndex+2*rtspol] = -prerot[1]*srot + prerot[0]*crot;
//
//   // }
//
// }
//
//

//
// extern "C" void RTS_copy_GPU_TileGains(float _Complex *TileGainMatrices, float _Complex *hostGains, int n_gains){
//
//   CUDA_SAFE_CALL( cudaMemcpy( hostGains, TileGainMatrices , n_gains*sizeof(float _Complex) ,
//                               cudaMemcpyDeviceToHost ) );
// }
//
// extern "C" void RTS_CUDA_ReleasePrimaryBeam_local( copy_cal_context_t *calContext ) {
//
//   if( calContext->d_M != NULL ) {
//     CUDA_SAFE_CALL( cudaFree( calContext->d_M ) );
//     calContext->d_M = NULL;
//   }
//
//   if( calContext->d_N != NULL ) {
//     CUDA_SAFE_CALL( cudaFree( calContext->d_N ) );
//     calContext->d_N = NULL;
//   }
//
//   if( calContext->d_Q1 != NULL ) {
//     CUDA_SAFE_CALL( cudaFree( calContext->d_Q1 ) );
//     calContext->d_Q1 = NULL;
//   }
//
//   if( calContext->d_Q2 != NULL ) {
//     CUDA_SAFE_CALL( cudaFree( calContext->d_Q2 ) );
//     calContext->d_Q2 = NULL;
//   }
//
//   if( calContext->d_norm_fac != NULL ) {
//     CUDA_SAFE_CALL( cudaFree( calContext->d_norm_fac ) );
//     calContext->d_norm_fac = NULL;
//   }
//
// }

// extern "C" void RTS_CUDA_AllocatePrimaryBeam_local( copy_cal_context_t *calContext,
//                                const int nStations ) {
//
//   const int n_pols=2; // instrumental pols
//   const int nMN = calContext->primary_beams[0].nMN; // max_length HDFBeamInit
//   int arrSize = n_pols * nMN;
//
//   if( calContext->d_M == NULL ) {
//     CUDA_SAFE_CALL( cudaMalloc( (void**)&calContext->d_M, nStations*arrSize*sizeof(float) ) );
//   } else {
//     cerr << __FUNCTION__ << ": Non-null d_M" << endl;
//     exit( EXIT_FAILURE );
//   }
//
//   if( calContext->d_N == NULL ) {
//     CUDA_SAFE_CALL( cudaMalloc( (void**)&calContext->d_N, nStations*arrSize*sizeof(float) ) );
//   } else {
//     cerr << __FUNCTION__ << ": Non-null d_N" << endl;
//     exit( EXIT_FAILURE );
//   }
//
//   if( calContext->d_Q1 == NULL ) {
//     CUDA_SAFE_CALL( cudaMalloc( (void**)&calContext->d_Q1, nStations*arrSize*sizeof(float _Complex) ) );
//   } else {
//     cerr << __FUNCTION__ << ": Non-null d_Q1" << endl;
//     exit( EXIT_FAILURE );
//   }
//
//   if( calContext->d_Q2 == NULL ) {
//     CUDA_SAFE_CALL( cudaMalloc( (void**)&calContext->d_Q2, nStations*arrSize*sizeof(float _Complex) ) );
//   } else {
//     cerr << __FUNCTION__ << ": Non-null d_Q2" << endl;
//     exit( EXIT_FAILURE );
//   }
//
//   if( calContext->d_norm_fac == NULL ) {
//     CUDA_SAFE_CALL( cudaMalloc( (void**)&calContext->d_norm_fac, MAX_POLS*sizeof(float _Complex) ) );
//   } else {
//     cerr << __FUNCTION__ << ": Non-null d_norm_fac" << endl;
//     exit( EXIT_FAILURE );
//   }
//
// }
