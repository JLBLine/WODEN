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


source_catalogue_t * read_source_catalogue(char *filename);
