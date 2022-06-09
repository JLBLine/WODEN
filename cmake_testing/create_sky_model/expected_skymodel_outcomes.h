#pragma once

#include <math.h>
#include "constants.h"
#include "woden_struct_defs.h"
#include "read_yaml_skymodel.h"
#include "read_text_skymodel.h"
#include "create_sky_model.h"

//First point component values--------------------------------------------------
extern double POINT0_RA;
extern double POINT0_DEC;
//Power law / curved power law expected outcomes
extern double POINT0_POW_FREQ;
extern double POINT0_POW_I;
extern double POINT0_POW_Q;
extern double POINT0_POW_U;
extern double POINT0_POW_V;
extern double POINT0_POW_SI;
extern double POINT0_CURVE_Q;

//List type flux model
extern double POINT0_LIST_FREQ[];
extern double POINT0_LIST_I[];
extern double POINT0_LIST_Q[];
extern double POINT0_LIST_U[];
extern double POINT0_LIST_V[];
extern int POINT0_NUM_LIST_ENTRIES;

//First gaussian component values-----------------------------------------------
extern double GAUSS0_RA;
extern double GAUSS0_DEC;
extern double GAUSS0_MAJ;
extern double GAUSS0_MIN;
extern double GAUSS0_PA;

//Power law / curved power law expected outcomes
extern double GAUSS0_POW_FREQ;
extern double GAUSS0_POW_I;
extern double GAUSS0_POW_Q;
extern double GAUSS0_POW_U;
extern double GAUSS0_POW_V;
extern double GAUSS0_POW_SI;
extern double GAUSS0_CURVE_Q;

//List type flux model
extern double GAUSS0_LIST_FREQ[];
extern double GAUSS0_LIST_I[];
extern double GAUSS0_LIST_Q[];
extern double GAUSS0_LIST_U[];
extern double GAUSS0_LIST_V[];
extern int GAUSS0_NUM_LIST_ENTRIES;

//First shapelet component values-----------------------------------------------
extern double SHAPE0_RA;
extern double SHAPE0_DEC;
extern double SHAPE0_MAJ;
extern double SHAPE0_MIN;
extern double SHAPE0_PA;
extern double SHAPE0_N1S[];
extern double SHAPE0_N2S[];
extern double SHAPE0_COEFFS[];

//Power law / curved power law expected outcomes
extern double SHAPE0_POW_FREQ;
extern double SHAPE0_POW_I;
extern double SHAPE0_POW_Q;
extern double SHAPE0_POW_U;
extern double SHAPE0_POW_V;
extern double SHAPE0_POW_SI;
extern double SHAPE0_CURVE_Q;


//List type flux model
extern double SHAPE0_LIST_FREQ[];
extern double SHAPE0_LIST_I[];
extern double SHAPE0_LIST_Q[];
extern double SHAPE0_LIST_U[];
extern double SHAPE0_LIST_V[];
extern int SHAPE0_NUM_LIST_ENTRIES;

/*
Check that the counts of number of SOURCEs and COMPONENT types in the
source_catalogue_t matches those provide through argument
*/
void check_single_source_numbers(source_catalogue_t *raw_srccat,
       int num_sources,  int num_shapelets,  int n_comps,
       int n_points, int n_point_lists, int n_point_powers, int n_point_curves,
       int n_gauss, int n_gauss_lists, int n_gauss_powers, int n_gauss_curves,
       int n_shapes, int n_shape_lists, int n_shape_powers, int n_shape_curves,
       int n_shape_coeffs, int source_index);

/*
Checks that the POINT values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_point(source_t source, int comp_index,
                                  e_flux_type flux_type, int flux_ind);

/*
Checks that the GAUSSIAN values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_gauss(source_t source, int comp_index,
                                  e_flux_type flux_type, int flux_ind);

/*
Checks that the SHAPELET values for the COMPONENT at the given index are correct
Assumes the catalogue being read in is one of the specific sky models made
for testing in the WODEN/cmake_testing/read_and_write dir
*/
void check_single_component_shapelet(source_t source, int comp_index,
                                     e_flux_type flux_type, int flux_ind);

/*
Test whether a single point source is read in correctly
Use this to test in some allmost broken srclists
*/
void test_read_skymodel_SinglePoint(char *srclist, e_flux_type flux_type);

/*
Test whether a single Gaussian source is read in correctly
*/
void test_read_skymodel_SingleGaussian(char *srclist, e_flux_type flux_type);

/*
Test whether a single shapelet source is read in correctly
*/
void test_read_skymodel_SingleShapelet(char *srclist, e_flux_type flux_type);

/*
Test whether three separate SOURCEs, each with a single COMPONENT,
ordered as POINT, GAUSSIAN, SHAPELET, read in correctly
*/
void test_read_skymodel_ThreeSources(char *srclist, e_flux_type flux_type);

/*
Test whether a single SOURCE with three COMPONENTS,
ordered as POINT, GAUSSIAN, SHAPELET, read in correctly
*/
void test_read_skymodel_ThreeComponents(char *srclist, e_flux_type flux_type);
