#pragma once
#include <stdlib.h>
#include "constants.h"

enum component_type {POINT=0, GAUSSIAN, SHAPELET, SHAPELET2};
typedef enum {NO_BEAM, GAUSS_BEAM, FEE_BEAM, ANALY_DIPOLE}e_beamtype;

typedef struct _catsource_t {
  //General source info
  char name[32];
  //Components info
  int n_comps;
  int n_points;
  int n_gauss;
  int n_shapes;
  int n_shape_coeffs;

  //Pointsource params
  float *point_ras;
  float *point_decs;
  float *point_ref_freqs;
  float *point_ref_stokesI;
  float *point_ref_stokesQ;
  float *point_ref_stokesU;
  float *point_ref_stokesV;
  float *point_SIs;
  float *point_freqs;
  float *point_azs;
  float *point_zas;
  float *sin_point_para_angs;
  float *cos_point_para_angs;

  //Gaussian params
  float *gauss_ras;
  float *gauss_decs;
  float *gauss_ref_freqs;
  float *gauss_ref_stokesI;
  float *gauss_ref_stokesQ;
  float *gauss_ref_stokesU;
  float *gauss_ref_stokesV;
  float *gauss_SIs;
  float *gauss_majors;
  float *gauss_minors;
  float *gauss_pas;
  float *gauss_azs;
  float *gauss_zas;
  float *sin_gauss_para_angs;
  float *cos_gauss_para_angs;

  //Shapelet params
  float *shape_ras;
  float *shape_decs;
  float *shape_ref_freqs;
  float *shape_ref_stokesI;
  float *shape_ref_stokesQ;
  float *shape_ref_stokesU;
  float *shape_ref_stokesV;
  float *shape_SIs;
  float *shape_coeffs;
  float *shape_n1s;
  float *shape_n2s;
  float *shape_majors;
  float *shape_minors;
  float *shape_pas;
  float *shape_param_indexes;
  float *shape_azs;
  float *shape_zas;
  float *sin_shape_para_angs;
  float *cos_shape_para_angs;

} catsource_t;

typedef struct _RTS_MWA_FEE_beam {
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

} RTS_MWA_FEE_beam_t;

typedef struct _beam_settings_t {
    // float *beam_angles_array;

    float gauss_sdec;
    float gauss_cdec;
    float gauss_ha;

    float *beam_point_has;
    float *beam_point_decs;
    int num_point_beam_values;

    float *beam_gausscomp_has;
    float *beam_gausscomp_decs;
    int num_gausscomp_beam_values;

    float *beam_shape_has;
    float *beam_shape_decs;
    int num_shape_beam_values;

    float beam_FWHM_rad;
    // float *beam_ref_freq_array;
    float beam_ref_freq;
    int beamtype;

    float *para_cosrot;
    float *para_sinrot;

    RTS_MWA_FEE_beam_t *FEE_beam;
    RTS_MWA_FEE_beam_t *FEE_beam_zenith;

} beam_settings_t;

typedef struct _source_catalogue_t {
    int num_sources;
    int num_shapelets;
    catsource_t *catsources;
    beam_settings_t *beam_settings;
} source_catalogue_t;

typedef struct _visibility_set_t {
  float *sum_visi_real;
  float *sum_visi_imag;
  float *us_metres;
  float *vs_metres;
  float *ws_metres;
  float *allsteps_sha0s;
  float *allsteps_cha0s;
  float *allsteps_lsts;
  float *allsteps_wavelengths;
  float *channel_frequencies;

  float *beam_has;
  float *beam_decs;
  float *beam_ls;
  float *beam_ms;
  float *beam_reals;
  float *beam_imags;

  float *sum_visi_XX_real;
  float *sum_visi_XX_imag;
  float *sum_visi_XY_real;
  float *sum_visi_XY_imag;
  float *sum_visi_YX_real;
  float *sum_visi_YX_imag;
  float *sum_visi_YY_real;
  float *sum_visi_YY_imag;

} visibility_set_t;

typedef struct _woden_settngs_t {
  float lst_base;
  float ra0;
  float dec0;
  float sdec0;
  float cdec0;
  int num_baselines;
  int num_freqs;
  float frequency_resolution;
  float base_low_freq;
  int num_time_steps;
  float time_res;
  const char* cat_filename;
  int num_bands;
  int *band_nums;
  int sky_crop_type;
  e_beamtype beamtype;
  float gauss_beam_FWHM;
  float gauss_beam_ref_freq;
  int chunking_size;
  const char* hdf5_beam_path;
  double jd_date;
  int EDA2_sim;
  int array_layout_file;
  const char* array_layout_file_path;
  float latitude;
  float longitude;
  float FEE_ideal_delays[16];
  float coarse_band_width;
  float gauss_ra_point;
  float gauss_dec_point;
  int num_visis;
  float base_band_freq;

} woden_settings_t;

typedef struct _array_layout_t {
    float *ant_X;
    float *ant_Y;
    float *ant_Z;
    float *X_diff_metres;
    float *Y_diff_metres;
    float *Z_diff_metres;
    float *ant_east;
    float *ant_north;
    float *ant_height;
    float latitude;
    int num_baselines;
    int num_tiles;
    float lst_base;

} array_layout_t;

// typedef struct _d_XYZ_t {
//     float *d_X_diff;
//     float *d_Y_diff;
//     float *d_Z_diff;
// } d_XYZ_t;
