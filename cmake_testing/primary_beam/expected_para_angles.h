#pragma once
#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "primary_beam.h"
#include "woden_struct_defs.h"

/*
Bunch of settings and expected output arrays. Defined in expected_para_angles.c
*/
extern float lsts[];
extern float point_ras[];
extern float point_decs[];
extern float gauss_ras[];
extern float gauss_decs[];
extern float shape_ras[];
extern float shape_decs[];
extern float expec_point_sin_para[];
extern float expec_point_cos_para[];
extern float expec_gauss_sin_para[];
extern float expec_gauss_cos_para[];
extern float expec_shape_sin_para[];
extern float expec_shape_cos_para[];

/*
Make the polpulated catsource_t struct. Stick in some necessary values
*/
catsource_t * make_sky_model(void) ;
