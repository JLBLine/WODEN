#pragma once
#include <math.h>
#include <stdint.h>
#include <fitsio.h>
// #include "FEE_primary_beam.h"
#include "woden_struct_defs.h"

#define CS_LEN_SRC_NAME 16
#define SRC_KEY         "SOURCE"
#define SRC_END         "ENDSOURCE"
#define COMP_KEY        "COMPONENT"
#define COMP_END        "ENDCOMPONENT"
#define FREQ_KEY        "FREQ"
#define LINEAR_KEY      "LINEAR"
#define POINT_KEY       "POINT"
#define GAUSSIAN_KEY    "GAUSSIAN"
#define GPARAMS_KEY     "GPARAMS"
#define SHAPELET_KEY    "SHAPELET"
#define SPARAMS_KEY     "SPARAMS"
#define SCOEFF_KEY      "SCOEFF"

source_catalogue_t * read_source_catalogue(const char *filename);

woden_settings_t * read_json_settings(const char *filename);

array_layout_t * calc_XYZ_diffs(woden_settings_t *woden_settings);

void RTS_precXYZ(double rmat[3][3], double x, double y, double z, double lmst,
         double *xp, double *yp, double *zp, double lmst2000);

void RTS_PrecessXYZtoJ2000( array_layout_t *array_layout,
                       woden_settings_t *woden_settings);
