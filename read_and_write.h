#include <math.h>
#include <stdint.h>
#include <fitsio.h>

#define CS_LEN_SRC_NAME 16
#define SRC_KEY         "SOURCE"
#define SRC_END         "ENDSOURCE"
#define COMP_KEY        "COMPONENT"
#define COMP_END        "ENDCOMPONENT"
#define FREQ_KEY        "FREQ"
#define POINT_KEY       "POINT"
#define GAUSSIAN_KEY    "GAUSSIAN"
#define GPARAMS_KEY     "GPARAMS"
#define SHAPELET_KEY    "SHAPELET"
#define S1PARAMS_KEY    "S1PARAMS"
#define S1COEFF_KEY     "S1COEFF"
#define SHAPELET2_KEY   "SHAPELET2"
#define S2PARAMS_KEY    "S2PARAMS"
#define S2COEFF_KEY     "S2COEFF"
#define DH2R 0.26179938779914943653855361527329190701643078328126
#define DD2R 0.017453292519943295769236907684886127134428718885417
#define VELC 299792458.0
#define SOLAR2SIDEREAL 1.00274
#define D2R M_PI/180.0
#define DEFAULT_SI -0.8
#define MWA_LAT -26.703319

enum component_type {POINT=0, GAUSSIAN, SHAPELET, SHAPELET2};

typedef struct _catsource_t {
  //General source info
  char name[32];
  //Components info
  int n_comps;
  int n_points;
  int n_gauss;
  int n_S1s;
  int n_S1_coeffs;
  int n_S2s;
  int n_S2_coeffs;

  //Pointsource params
  float *point_ras;
  float *point_decs;
  float *point_fluxes;
  float *point_freqs;

  float *d_point_ras;
  float *d_point_decs;
  float *d_point_fluxes;
  float *d_point_freqs;

  //Gaussian params
  float *gauss_ras;
  float *gauss_decs;
  float *gauss_fluxes;
  float *gauss_freqs;
  float *gauss_majors;
  float *gauss_minors;
  float *gauss_pas;

  float *S1_ras;
  float *S1_decs;
  float *S1_fluxes;
  float *S1_freqs;
  float *S1_coeffs;
  float *S1_n1s;
  float *S1_n2s;
  float *S1_majors;
  float *S1_minors;
  float *S1_pas;
  float *S1_param_indexes;

  float *S2_ras;
  float *S2_decs;
  float *S2_fluxes;
  float *S2_freqs;
  float *S2_coeffs;
  float *S2_n1s;
  float *S2_n2s;
  float *S2_majors;
  float *S2_minors;
  float *S2_pas;
  float *S2_param_indexes;

  // float *n_coeffs

} catsource_t;

typedef struct _source_catalogue_t {
    int num_sources;
    catsource_t *catsources;
} source_catalogue_t;

typedef struct _visibility_set_t {
  float *sum_visi_real;
  float *sum_visi_imag;
  float *us_metres;
  float *vs_metres;
  float *ws_metres;
  float *sha0s;
  float *cha0s;
  float *lsts;


  // float *d_sum_visi_real;
  // float *d_sum_visi_imag;
  // float *d_us_metres;
  // float *d_vs_metres;
  // float *d_ws_metres;
} visibility_set_t;

typedef struct _woden_settngs_t {
  float lst_base;
  float ra0;
  float dec0;
  int num_baselines;
  int num_freqs;
  float frequency_resolution;
  float base_low_freq;
  int num_time_steps;
  float time_res;
  const char* cat_filename;
  const char* metafits_filename;
  int num_bands;
  int *band_nums;

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
    int num_antennas;
    float lst_base;
} array_layout_t;

//Stolen from the RTS
typedef struct _Meta_Ffile {
    char recvrs[38];
    char calib[1];
    char calsrc[16];
    char tileflg[16];
    char version[4];
    char mwaVersion;
    long naxes[2];
    int inps[256], ants[256], tile[256], flags[256];
    float E[256], N[256], H[256];
    char leng[256][14];
    int ddlys[256][16];
    int dig_gains[256][24];
    int centchan;
    float lst_base;
    float frequency_resolution;
    float frequency_cent;
    float bandwidth;
    float base_low_freq;
    float time_res;

} MetaFfile_t;

source_catalogue_t * read_source_catalogue(const char *filename);

woden_settings_t * read_json_settings(const char *filename);

int init_meta_file(fitsfile *mfptr, MetaFfile_t *metafits, const char *nome);

array_layout_t * calc_XYZ_diffs(MetaFfile_t *metafits);
