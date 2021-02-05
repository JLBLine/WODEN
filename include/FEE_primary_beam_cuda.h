#define kP1BlockSize 64
#include "cudacheck.h"
#include "woden_struct_defs.h"

extern "C" void calc_CUDA_FEE_beam(float *azs, float *zas,
                                   float *sin_para_angs, float *cos_para_angs,
                                   int num_components, int num_time_steps,
                                   copy_primary_beam_t *FEE_beam,
                                   int rotation, int scaling);

__global__ void kern_rotate_FEE_beam(cuFloatComplex *d_FEE_beam_gain_matrices,
                                float *d_sin_para_angs, float *d_cos_para_angs,
                                int num_components, int num_time_steps);

extern "C" void get_HDFBeam_normalisation(copy_primary_beam_t *FEE_beam_zenith,
                copy_primary_beam_t *FEE_beam);

extern "C" void free_FEE_primary_beam_from_GPU(copy_primary_beam_t *primary_beam);

extern "C" void copy_FEE_primary_beam_to_GPU(copy_primary_beam_t *FEE_beam,
                                             int num_time_steps);

extern "C" void calc_FEE_beam(float *az, float *za, int num_azza,
           float *sin_para_angs, float *cos_para_angs,
           int num_time_steps, copy_primary_beam_t *primary_beam,
           float _Complex *h_FEE_beam_gains,
           int scaling);

extern "C" void RTS_CUDA_get_TileGains(float *phi, float *theta,
           float *sin_para_angs, float *cos_para_angs,
           int num_time_steps, int num_components,
           float rotation, copy_primary_beam_t *primary_beam,
           float _Complex *TileGainMatrices, int scaling);

__global__ void RTS_P1SINfKernel( float *d_theta, cuFloatComplex *rts_P_sin,
           cuFloatComplex *rts_p1, int nmax, int num_coords);

__global__ void RTS_getTileGainsKernel( float *d_phi, float *d_theta,
           int nMN, int num_coords,
           float *pb_M, float *pb_N,
           cuFloatComplex *pb_Q1, cuFloatComplex *pb_Q2,
           cuFloatComplex *rts_P_sin, cuFloatComplex *rts_P1,
           cuFloatComplex *emn_T, cuFloatComplex *emn_P);

__global__ void kern_sum_emn_PT_by_M(cuFloatComplex *emn_T, cuFloatComplex *emn_P,
           float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
           float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
           float *d_m_range, float *d_M, int nMN, int nmax, int num_coords);

__global__ void kern_calc_sigmaTP(cuFloatComplex *TileGainMatrices,
                float *d_emn_T_sum_real, float *d_emn_T_sum_imag,
                float *d_emn_P_sum_real, float *d_emn_P_sum_imag,
                cuFloatComplex *d_emn_T_sum, cuFloatComplex *d_emn_P_sum,
                int nmax, int num_coords);

__global__ void kern_apply_FEE_norm(cuFloatComplex *TileGainMatrices,
           cuFloatComplex *d_norm_fac, int num_coords );

__global__ void kern_map_FEE_beam_gains(cuFloatComplex *d_FEE_beam_gain_matrices,
    cuFloatComplex *d_primay_beam_J00, cuFloatComplex *d_primay_beam_J01,
    cuFloatComplex *d_primay_beam_J10, cuFloatComplex *d_primay_beam_J11,
    int num_freqs, int num_components, int num_visis, int num_baselines,
    int num_times);

extern "C" void test_RTS_CUDA_FEE_beam(int num_components,
           float *azs, float *zas,
           float *sin_para_angs, float *cos_para_angs,
           copy_primary_beam_t *FEE_beam_zenith,
           copy_primary_beam_t *FEE_beam,
           int rotation, int scaling,
           float _Complex *FEE_beam_gains);
