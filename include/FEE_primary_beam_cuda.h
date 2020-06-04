#include "FEE_primary_beam.h"
#define kP1BlockSize 64
#include "cudacheck.h"

// extern "C" void RTS_CUDA_AllocatePrimaryBeam_local( copy_cal_context_t *calContext,
//                                const int nStations );
//
// extern "C" void RTS_CUDA_SendPrimaryBeam_local( copy_cal_context_t *calContext,
                           // const int nStations );
//
// extern "C" void RTS_CUDA_allocate_TileGains(float _Complex **TileGainMatrices, int nStations);

//
// extern "C" void RTS_copy_GPU_TileGains(float _Complex *TileGainMatrices, float _Complex *hostGains, int n_gains);
//
// extern "C" void RTS_CUDA_ReleasePrimaryBeam_local( copy_cal_context_t *calContext );

//
extern "C" void RTS_CUDA_get_TileGains(float phi, float theta, float rotation, copy_primary_beam_t *primary_beam, int nStations, float _Complex *TileGainMatrices, int scaling);

extern "C" void test_FEE_beam(float *az, float *za, copy_primary_beam_t *primary_beam,
           float _Complex *copy_hostGains,
           float _Complex *h_rts_P_sin, float _Complex *h_rts_P1,
           float _Complex *h_emn_T, float _Complex *h_emn_P,
           float _Complex *h_emn_T_sum, float _Complex *h_emn_P_sum);

__global__ void RTS_P1SINfKernel( float theta, cuFloatComplex *rts_P_sin, cuFloatComplex *rts_p1, int nmax);

__global__ void RTS_getTileGainsKernel( float phi, float theta, int nMN,
           float *pb_M, float *pb_N,
           cuFloatComplex *pb_Q1, cuFloatComplex *pb_Q2,
           cuFloatComplex *rts_P_sin, cuFloatComplex *rts_P1,
           cuFloatComplex *Sigma_T, cuFloatComplex *Sigma_P,
           cuFloatComplex *emn_T, cuFloatComplex *emn_P);

// __global__ void RTS_reduce_SigmaPTKernel(int nMN, cuFloatComplex *Sigma_T,
//            cuFloatComplex *Sigma_P, cuFloatComplex *TileGainMatrices,
//            cuFloatComplex *norm_fac, int scaling, float rotation);

__global__ void kern_sum_emn_PT_by_M(cuFloatComplex *emn_T, cuFloatComplex *emn_P,
           float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
           float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
           float *d_m_range, float *d_M, int nMN, int nmax);

__global__ void kern_calc_sigmaTP(cuFloatComplex *TileGainMatrices,
                cuFloatComplex *d_phi_comp,
                float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
                float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
                cuFloatComplex *d_emn_T_sum, cuFloatComplex *d_emn_P_sum,
                int nmax);

__global__ void kern_apply_norm(cuFloatComplex *TileGainMatrices,
           cuFloatComplex *d_norm_fac );
