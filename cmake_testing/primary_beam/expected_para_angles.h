#pragma once
#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "primary_beam.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

/*
Bunch of settings and expected output arrays. Defined in expected_para_angles.c
*/
extern double lsts[];
extern double point_ras[];
extern double point_decs[];
extern double gauss_ras[];
extern double gauss_decs[];
extern double shape_ras[];
extern double shape_decs[];
extern double expec_point_sin_para[];
extern double expec_point_cos_para[];
extern double expec_gauss_sin_para[];
extern double expec_gauss_cos_para[];
extern double expec_shape_sin_para[];
extern double expec_shape_cos_para[];

/*
Make the polpulated catsource_t struct. Stick in some necessary values
*/
catsource_t * make_sky_model(void);
