#pragma once
#define N_LNA_FREQS 451
#define MAX_ZMATRIX_FREQS 256
#define NUM_DIPOLES 16
#define SIZ_Z_MATRIX (NUM_DIPOLES*2*NUM_DIPOLES*2)
#define VEL_LIGHT 299792458.0
#define DQ (435e-12*VEL_LIGHT)  // delay quantum of the MWA beamformer in meters.
#define MAX_POLS 4
#define N_COPOL 2
#define NUM_PARAMS_PER_DIPOLE 2

#include "hdf5.h"
#include "hdf5_hl.h"

typedef struct _copy_primary_beam {
  float *parameters;
  float ref_freq;               //!< if used to generate reference J matrices, some applications need a frequency.
  double _Complex **Q1, **Q2;      //!< Beam modes used for Spherical Harmonic model
  double _Complex **p_T, **p_P;    //!< Some pre-computed values used in tile response
  double **M,**N;
  int nmax;
  int nMN;
  float _Complex norm_fac[MAX_POLS];

  // BP 2019: All the Spherical Harmonic Beam data are double
  // so we will use them on the GPUs as well or there will be all kinds
  // of issues with copying

  float _Complex *d_Q1, *d_Q2;
  float _Complex *d_p_T, *d_p_P; // Precomputed values used in CPU optimization
  float *d_M, *d_N;
  int d_nmax;
  int d_nMN;
  float _Complex d_norm_fac[MAX_POLS];

  float _Complex *emn_P;
  float _Complex *emn_T;

  float _Complex *d_emn_T_sum;
  float _Complex *d_emn_P_sum;

  float _Complex *rts_P1;
  float _Complex *rts_P_sin;

  float *m_range;

  float _Complex *d_FEE_beam_gain_matrices;

  float *d_para_cosrot;
  float *d_para_sinrot;

} copy_primary_beam_t;

/***************************************************************
Data Structure used by H5Literate.
***************************************************************/
struct opdata {
  float freq_in;
  float freq_out;
  float least_diff;
};

/*
 * Operator function to be called by H5Literate.
 */
herr_t RTS_op_func (hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data);

// void RTS_reorderDelays(float in[NUM_DIPOLES], float out[NUM_DIPOLES]);
//
// void RTS_reorderDelays2RTS(float in[NUM_DIPOLES], float out[NUM_DIPOLES]);

int RTS_HDFBeamInit(char *h5filename, float freq_Hz, copy_primary_beam_t *pb,
                float *FEE_delays, int stn);

void RTS_P1SIN(int nmax, float theta, double _Complex **P_sin, double _Complex **P1);
//
int RTS_get_FF2(float phi, float theta, copy_primary_beam_t *pb, float _Complex result[4], int scaling, int clean_up);
//
int RTS_getJonesSphHarm(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);
//
int RTS_getJonesDipole(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);
//
int RTS_getTileResponse(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4], int scaling, float rotation);
//
void RTS_freeHDFBeam(copy_primary_beam_t *pb);
//
void RTS_HDFBeamCleanUp();
