/*******************************************************************************
`crop_sky_model` takes in a sky catalogue, and crops everything below the
horizon for the first time step. Has to crop three different COMPONENT
types, with SHAPELETs being the trickiest as mulitple basis function parameters
from a single COMPONENT have to be matched correctly in the cropped sky model
Can also crop in two different ways - if one COMPONENT in a SOURCE is below
horizon, through out the whole source, or just crop out the COMPONENTs that are
below the horizon. Former is good if you have a discrete set of foreground like
sources. The latter is best if you have a massive diffuse sky SOURCE.
*******************************************************************************/
#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "create_sky_model.h"
#include "woden_struct_defs.h"
#include "test_crop_sky_model.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL = 1e-12;
#else
  double TOL = 1e-6;
#endif

/*
This makes a dummy sky model, including as many of the three COMPONENT types
as requested. This produces an example of what `read_and_write::read_source_catalogue`
would output.

We'll always make two SOURCEs here, and add COMPONENTs of different types
based on the point, gauss, shape args. Dummy values are stored in
test_crop_sky_model.h
*/
source_catalogue_t * create_input_srccat(int point, int gauss, int shape){

  source_catalogue_t *raw_srccat = malloc(sizeof(source_catalogue_t));

  //Initially set number of shapelets to zero

  raw_srccat->num_shapelets = 0;

  //Allocate two catsources (two SOURCEs)
  raw_srccat->num_sources = 2;
  raw_srccat->sources = malloc(2*sizeof(source_t));

  //set the COMPONENT counters to zero for both SOURCEs

  for (int src = 0; src < 2; src++) {
    raw_srccat->sources[src].n_comps = 0;
    raw_srccat->sources[src].n_points = 0;
    raw_srccat->sources[src].n_gauss = 0;
    raw_srccat->sources[src].n_shapes = 0;
    raw_srccat->sources[src].n_shape_coeffs = 0;

    raw_srccat->sources[src].n_point_lists = 0;
    raw_srccat->sources[src].n_point_powers = 0;
    raw_srccat->sources[src].n_point_curves = 0;
    raw_srccat->sources[src].n_gauss_lists = 0;
    raw_srccat->sources[src].n_gauss_powers = 0;
    raw_srccat->sources[src].n_gauss_curves = 0;
    raw_srccat->sources[src].n_shape_lists = 0;
    raw_srccat->sources[src].n_shape_powers = 0;
    raw_srccat->sources[src].n_shape_curves = 0;


  }

  if (point) {
    //Assign the first SOURCE some POINT COMPONENTs
    raw_srccat->sources[0].point_components.ras = comp0_ras;
    raw_srccat->sources[0].point_components.decs = comp0_decs;

    //Power law stuff
    raw_srccat->sources[0].point_components.power_ref_freqs = point0_ref_freqs;
    raw_srccat->sources[0].point_components.power_ref_stokesI = point0_ref_stokesI;
    raw_srccat->sources[0].point_components.power_ref_stokesQ = point0_ref_stokesQ;
    raw_srccat->sources[0].point_components.power_ref_stokesU = point0_ref_stokesU;
    raw_srccat->sources[0].point_components.power_ref_stokesV = point0_ref_stokesV;
    raw_srccat->sources[0].point_components.power_SIs = point0_SIs;
    raw_srccat->sources[0].point_components.power_comp_inds = point0_power_inds;

    //curved power law stuff
    raw_srccat->sources[0].point_components.curve_ref_freqs = point0_ref_freqs;
    raw_srccat->sources[0].point_components.curve_ref_stokesI = point0_ref_stokesI;
    raw_srccat->sources[0].point_components.curve_ref_stokesQ = point0_ref_stokesQ;
    raw_srccat->sources[0].point_components.curve_ref_stokesU = point0_ref_stokesU;
    raw_srccat->sources[0].point_components.curve_ref_stokesV = point0_ref_stokesV;
    raw_srccat->sources[0].point_components.curve_SIs = point0_SIs;
    raw_srccat->sources[0].point_components.curve_comp_inds = point0_power_inds;
    raw_srccat->sources[0].point_components.curve_qs = point0_curve_qs;

    //List flux type things
    raw_srccat->sources[0].point_components.list_comp_inds = point0_list_comp_inds;
    raw_srccat->sources[0].point_components.num_list_values = point0_num_list_values;
    raw_srccat->sources[0].point_components.list_start_indexes = point0_list_start_indexes;
    raw_srccat->sources[0].point_components.list_freqs = point0_list_freqs;
    raw_srccat->sources[0].point_components.list_stokesI = point0_list_stokesI;
    raw_srccat->sources[0].point_components.list_stokesQ = point0_list_stokesQ;
    raw_srccat->sources[0].point_components.list_stokesU = point0_list_stokesU;
    raw_srccat->sources[0].point_components.list_stokesV = point0_list_stokesV;

    raw_srccat->sources[0].n_comps += 9;
    raw_srccat->sources[0].n_points += 9;
    raw_srccat->sources[0].n_point_powers += 3;
    raw_srccat->sources[0].n_point_curves += 3;
    raw_srccat->sources[0].n_point_lists += 3;

    //Assign the second SOURCE some POINT COMPONENTs
    raw_srccat->sources[1].point_components.ras = comp1_ras;
    raw_srccat->sources[1].point_components.decs = comp1_decs;
    raw_srccat->sources[1].point_components.power_ref_freqs = point1_ref_freqs;
    raw_srccat->sources[1].point_components.power_ref_stokesI = point1_ref_stokesI;
    raw_srccat->sources[1].point_components.power_ref_stokesQ = point1_ref_stokesQ;
    raw_srccat->sources[1].point_components.power_ref_stokesU = point1_ref_stokesU;
    raw_srccat->sources[1].point_components.power_ref_stokesV = point1_ref_stokesV;
    raw_srccat->sources[1].point_components.power_SIs = point1_SIs;
    raw_srccat->sources[1].point_components.power_comp_inds = point1_power_inds;

    //curved power law stuff
    raw_srccat->sources[1].point_components.curve_ref_freqs = point1_ref_freqs;
    raw_srccat->sources[1].point_components.curve_ref_stokesI = point1_ref_stokesI;
    raw_srccat->sources[1].point_components.curve_ref_stokesQ = point1_ref_stokesQ;
    raw_srccat->sources[1].point_components.curve_ref_stokesU = point1_ref_stokesU;
    raw_srccat->sources[1].point_components.curve_ref_stokesV = point1_ref_stokesV;
    raw_srccat->sources[1].point_components.curve_SIs = point1_SIs;
    raw_srccat->sources[1].point_components.curve_comp_inds = point1_power_inds;
    raw_srccat->sources[1].point_components.curve_qs = point1_curve_qs;

    raw_srccat->sources[1].n_comps += 9;
    raw_srccat->sources[1].n_points += 9;
    raw_srccat->sources[1].n_point_powers += 3;
    raw_srccat->sources[1].n_point_curves += 3;
    raw_srccat->sources[1].n_point_lists += 3;

    //List flux type things
    raw_srccat->sources[1].point_components.list_comp_inds = point1_list_comp_inds;
    raw_srccat->sources[1].point_components.num_list_values = point1_num_list_values;
    raw_srccat->sources[1].point_components.list_start_indexes = point1_list_start_indexes;
    raw_srccat->sources[1].point_components.list_freqs = point1_list_freqs;
    raw_srccat->sources[1].point_components.list_stokesI = point1_list_stokesI;
    raw_srccat->sources[1].point_components.list_stokesQ = point1_list_stokesQ;
    raw_srccat->sources[1].point_components.list_stokesU = point1_list_stokesU;
    raw_srccat->sources[1].point_components.list_stokesV = point1_list_stokesV;

  }

  if (gauss) {
    //Assign the first SOURCE some GAUSSIAN COMPONENTs
    raw_srccat->sources[0].gauss_components.ras = comp0_ras;
    raw_srccat->sources[0].gauss_components.decs = comp0_decs;
    raw_srccat->sources[0].gauss_components.power_SIs = gauss0_SIs;
    raw_srccat->sources[0].gauss_components.majors = gauss0_majors;
    raw_srccat->sources[0].gauss_components.minors = gauss0_minors;
    raw_srccat->sources[0].gauss_components.pas = gauss0_pas;

    //power law stuff
    raw_srccat->sources[0].gauss_components.power_ref_freqs = gauss0_ref_freqs;
    raw_srccat->sources[0].gauss_components.power_ref_stokesI = gauss0_ref_stokesI;
    raw_srccat->sources[0].gauss_components.power_ref_stokesQ = gauss0_ref_stokesQ;
    raw_srccat->sources[0].gauss_components.power_ref_stokesU = gauss0_ref_stokesU;
    raw_srccat->sources[0].gauss_components.power_ref_stokesV = gauss0_ref_stokesV;
    raw_srccat->sources[0].gauss_components.power_comp_inds = gauss0_power_inds;

    //curved power law stuff
    raw_srccat->sources[0].gauss_components.curve_ref_freqs = gauss0_ref_freqs;
    raw_srccat->sources[0].gauss_components.curve_ref_stokesI = gauss0_ref_stokesI;
    raw_srccat->sources[0].gauss_components.curve_ref_stokesQ = gauss0_ref_stokesQ;
    raw_srccat->sources[0].gauss_components.curve_ref_stokesU = gauss0_ref_stokesU;
    raw_srccat->sources[0].gauss_components.curve_ref_stokesV = gauss0_ref_stokesV;
    raw_srccat->sources[0].gauss_components.curve_SIs = gauss0_SIs;
    raw_srccat->sources[0].gauss_components.curve_comp_inds = gauss0_power_inds;
    raw_srccat->sources[0].gauss_components.curve_qs = gauss0_curve_qs;

    //List flux type things
    raw_srccat->sources[0].gauss_components.list_comp_inds = gauss0_list_comp_inds;
    raw_srccat->sources[0].gauss_components.num_list_values = gauss0_num_list_values;
    raw_srccat->sources[0].gauss_components.list_start_indexes = gauss0_list_start_indexes;
    raw_srccat->sources[0].gauss_components.list_freqs = gauss0_list_freqs;
    raw_srccat->sources[0].gauss_components.list_stokesI = gauss0_list_stokesI;
    raw_srccat->sources[0].gauss_components.list_stokesQ = gauss0_list_stokesQ;
    raw_srccat->sources[0].gauss_components.list_stokesU = gauss0_list_stokesU;
    raw_srccat->sources[0].gauss_components.list_stokesV = gauss0_list_stokesV;

    raw_srccat->sources[0].n_comps += 9;
    raw_srccat->sources[0].n_gauss += 9;
    raw_srccat->sources[0].n_gauss_powers += 3;
    raw_srccat->sources[0].n_gauss_curves += 3;
    raw_srccat->sources[0].n_gauss_lists += 3;

    ///SECOND SOURCE TIME-------------------------------------------------------

    //Assign the first SOURCE some GAUSSIAN COMPONENTs
    raw_srccat->sources[1].gauss_components.ras = comp1_ras;
    raw_srccat->sources[1].gauss_components.decs = comp1_decs;
    raw_srccat->sources[1].gauss_components.power_SIs = gauss1_SIs;
    raw_srccat->sources[1].gauss_components.majors = gauss1_majors;
    raw_srccat->sources[1].gauss_components.minors = gauss1_minors;
    raw_srccat->sources[1].gauss_components.pas = gauss1_pas;

    //power law stuff
    raw_srccat->sources[1].gauss_components.power_ref_freqs = gauss1_ref_freqs;
    raw_srccat->sources[1].gauss_components.power_ref_stokesI = gauss1_ref_stokesI;
    raw_srccat->sources[1].gauss_components.power_ref_stokesQ = gauss1_ref_stokesQ;
    raw_srccat->sources[1].gauss_components.power_ref_stokesU = gauss1_ref_stokesU;
    raw_srccat->sources[1].gauss_components.power_ref_stokesV = gauss1_ref_stokesV;
    raw_srccat->sources[1].gauss_components.power_comp_inds = gauss1_power_inds;

    //curved power law stuff
    raw_srccat->sources[1].gauss_components.curve_ref_freqs = gauss1_ref_freqs;
    raw_srccat->sources[1].gauss_components.curve_ref_stokesI = gauss1_ref_stokesI;
    raw_srccat->sources[1].gauss_components.curve_ref_stokesQ = gauss1_ref_stokesQ;
    raw_srccat->sources[1].gauss_components.curve_ref_stokesU = gauss1_ref_stokesU;
    raw_srccat->sources[1].gauss_components.curve_ref_stokesV = gauss1_ref_stokesV;
    raw_srccat->sources[1].gauss_components.curve_SIs = gauss1_SIs;
    raw_srccat->sources[1].gauss_components.curve_comp_inds = gauss1_power_inds;
    raw_srccat->sources[1].gauss_components.curve_qs = gauss1_curve_qs;

    //List flux type things
    raw_srccat->sources[1].gauss_components.list_comp_inds = gauss1_list_comp_inds;
    raw_srccat->sources[1].gauss_components.num_list_values = gauss1_num_list_values;
    raw_srccat->sources[1].gauss_components.list_start_indexes = gauss1_list_start_indexes;
    raw_srccat->sources[1].gauss_components.list_freqs = gauss1_list_freqs;
    raw_srccat->sources[1].gauss_components.list_stokesI = gauss1_list_stokesI;
    raw_srccat->sources[1].gauss_components.list_stokesQ = gauss1_list_stokesQ;
    raw_srccat->sources[1].gauss_components.list_stokesU = gauss1_list_stokesU;
    raw_srccat->sources[1].gauss_components.list_stokesV = gauss1_list_stokesV;

    raw_srccat->sources[1].n_comps += 9;
    raw_srccat->sources[1].n_gauss += 9;
    raw_srccat->sources[1].n_gauss_powers += 3;
    raw_srccat->sources[1].n_gauss_curves += 3;
    raw_srccat->sources[1].n_gauss_lists += 3;

  }

  if (shape) {
    //Assign the first SOURCE some shapeIAN COMPONENTs
    raw_srccat->sources[0].shape_components.ras = comp0_ras;
    raw_srccat->sources[0].shape_components.decs = comp0_decs;
    raw_srccat->sources[0].shape_components.power_SIs = shape0_SIs;
    raw_srccat->sources[0].shape_components.majors = shape0_majors;
    raw_srccat->sources[0].shape_components.minors = shape0_minors;
    raw_srccat->sources[0].shape_components.pas = shape0_pas;

    //power law stuff
    raw_srccat->sources[0].shape_components.power_ref_freqs = shape0_ref_freqs;
    raw_srccat->sources[0].shape_components.power_ref_stokesI = shape0_ref_stokesI;
    raw_srccat->sources[0].shape_components.power_ref_stokesQ = shape0_ref_stokesQ;
    raw_srccat->sources[0].shape_components.power_ref_stokesU = shape0_ref_stokesU;
    raw_srccat->sources[0].shape_components.power_ref_stokesV = shape0_ref_stokesV;
    raw_srccat->sources[0].shape_components.power_comp_inds = shape0_power_inds;

    //curved power law stuff
    raw_srccat->sources[0].shape_components.curve_ref_freqs = shape0_ref_freqs;
    raw_srccat->sources[0].shape_components.curve_ref_stokesI = shape0_ref_stokesI;
    raw_srccat->sources[0].shape_components.curve_ref_stokesQ = shape0_ref_stokesQ;
    raw_srccat->sources[0].shape_components.curve_ref_stokesU = shape0_ref_stokesU;
    raw_srccat->sources[0].shape_components.curve_ref_stokesV = shape0_ref_stokesV;
    raw_srccat->sources[0].shape_components.curve_SIs = shape0_SIs;
    raw_srccat->sources[0].shape_components.curve_comp_inds = shape0_curve_inds;
    raw_srccat->sources[0].shape_components.curve_qs = shape0_curve_qs;

    //List flux type things
    raw_srccat->sources[0].shape_components.list_comp_inds = shape0_list_comp_inds;
    raw_srccat->sources[0].shape_components.num_list_values = shape0_num_list_values;
    raw_srccat->sources[0].shape_components.list_start_indexes = shape0_list_start_indexes;
    raw_srccat->sources[0].shape_components.list_freqs = shape0_list_freqs;
    raw_srccat->sources[0].shape_components.list_stokesI = shape0_list_stokesI;
    raw_srccat->sources[0].shape_components.list_stokesQ = shape0_list_stokesQ;
    raw_srccat->sources[0].shape_components.list_stokesU = shape0_list_stokesU;
    raw_srccat->sources[0].shape_components.list_stokesV = shape0_list_stokesV;

    //shapelet basis function stuff
    raw_srccat->sources[0].shape_components.shape_coeffs = shape0_coeffs;
    raw_srccat->sources[0].shape_components.n1s = shape0_n1s;
    raw_srccat->sources[0].shape_components.n2s = shape0_n2s;
    raw_srccat->sources[0].shape_components.param_indexes = shape0_param_indexes;

    raw_srccat->sources[0].n_comps += 9;
    raw_srccat->sources[0].n_shapes += 9;
    raw_srccat->sources[0].n_shape_powers += 3;
    raw_srccat->sources[0].n_shape_curves += 3;
    raw_srccat->sources[0].n_shape_lists += 3;
    raw_srccat->sources[0].n_shape_coeffs += NUM_SHAPE0_COEFFS;

    ///SECOND SOURCE TIME-------------------------------------------------------

    //Assign the first SOURCE some shapeIAN COMPONENTs
    raw_srccat->sources[1].shape_components.ras = comp1_ras;
    raw_srccat->sources[1].shape_components.decs = comp1_decs;
    raw_srccat->sources[1].shape_components.power_SIs = shape1_SIs;
    raw_srccat->sources[1].shape_components.majors = shape1_majors;
    raw_srccat->sources[1].shape_components.minors = shape1_minors;
    raw_srccat->sources[1].shape_components.pas = shape1_pas;

    //power law stuff
    raw_srccat->sources[1].shape_components.power_ref_freqs = shape1_ref_freqs;
    raw_srccat->sources[1].shape_components.power_ref_stokesI = shape1_ref_stokesI;
    raw_srccat->sources[1].shape_components.power_ref_stokesQ = shape1_ref_stokesQ;
    raw_srccat->sources[1].shape_components.power_ref_stokesU = shape1_ref_stokesU;
    raw_srccat->sources[1].shape_components.power_ref_stokesV = shape1_ref_stokesV;
    raw_srccat->sources[1].shape_components.power_comp_inds = shape1_power_inds;

    //curved power law stuff
    raw_srccat->sources[1].shape_components.curve_ref_freqs = shape1_ref_freqs;
    raw_srccat->sources[1].shape_components.curve_ref_stokesI = shape1_ref_stokesI;
    raw_srccat->sources[1].shape_components.curve_ref_stokesQ = shape1_ref_stokesQ;
    raw_srccat->sources[1].shape_components.curve_ref_stokesU = shape1_ref_stokesU;
    raw_srccat->sources[1].shape_components.curve_ref_stokesV = shape1_ref_stokesV;
    raw_srccat->sources[1].shape_components.curve_SIs = shape1_SIs;
    raw_srccat->sources[1].shape_components.curve_comp_inds = shape1_curve_inds;
    raw_srccat->sources[1].shape_components.curve_qs = shape1_curve_qs;

    //List flux type things
    raw_srccat->sources[1].shape_components.list_comp_inds = shape1_list_comp_inds;
    raw_srccat->sources[1].shape_components.num_list_values = shape1_num_list_values;
    raw_srccat->sources[1].shape_components.list_start_indexes = shape1_list_start_indexes;
    raw_srccat->sources[1].shape_components.list_freqs = shape1_list_freqs;
    raw_srccat->sources[1].shape_components.list_stokesI = shape1_list_stokesI;
    raw_srccat->sources[1].shape_components.list_stokesQ = shape1_list_stokesQ;
    raw_srccat->sources[1].shape_components.list_stokesU = shape1_list_stokesU;
    raw_srccat->sources[1].shape_components.list_stokesV = shape1_list_stokesV;

    //shapelet basis function stuff
    raw_srccat->sources[1].shape_components.shape_coeffs = shape1_coeffs;
    raw_srccat->sources[1].shape_components.n1s = shape1_n1s;
    raw_srccat->sources[1].shape_components.n2s = shape1_n2s;
    raw_srccat->sources[1].shape_components.param_indexes = shape1_param_indexes;

    raw_srccat->sources[1].n_comps += 9;
    raw_srccat->sources[1].n_shapes += 9;
    raw_srccat->sources[1].n_shape_powers += 3;
    raw_srccat->sources[1].n_shape_curves += 3;
    raw_srccat->sources[1].n_shape_lists += 3;
    raw_srccat->sources[1].n_shape_coeffs += NUM_SHAPE1_COEFFS;

  }

  return raw_srccat;
}

/*
This checks that the overall number of retained COMPONENTs is correct. This
depends on whether I've injected POINT, GAUSSIAN, or SHAPELETS into the sky,
whether I'm cropping by SOURCE or COMPONENT, and what the LST is
*/

void check_num_components_after_crop(double lst_base, source_t *cropped_src,
                                    e_sky_crop sky_crop_type,
                                    int point, int gauss, int shape) {

  //I've given the POINT, GAUSSIAN and SHAPELET components the same RA/Decs,
  //so the number of overall retained components is a multiple of how those
  //RA/DECs convert to az/za as to how many overall COMPONENTs are retained

  //For LST0, all COMPONENTs of the first SOURCE
  //(as defined in test_crop_sky_model.h) should have been retained.
  //Only the first COMPONENT of the second SOURCE is above the horizon, so
  //if cropping by SOURCE there should only be three COMPONENTs retained.
  //If cropping by COMPONENT, there should be a 4th component

  int num_component_types = point + gauss + shape;

  if (lst_base == LST0) {

    if (sky_crop_type == CROP_SOURCES) {
      TEST_ASSERT_EQUAL_INT(9*num_component_types, cropped_src->n_comps);

    } else if (sky_crop_type == CROP_COMPONENTS) {
      TEST_ASSERT_EQUAL_INT(12*num_component_types, cropped_src->n_comps);
    }
  }
  else if (lst_base == LST1) { //For this LST, both SOURCEs should be retained
    TEST_ASSERT_EQUAL_INT(18*num_component_types, cropped_src->n_comps);
  }
  else if (lst_base == LST2) { //first source should have been cropped, second retained
    TEST_ASSERT_EQUAL_INT(9*num_component_types, cropped_src->n_comps);
  }
  else if (lst_base == LST3) { //both should be flagged
    TEST_ASSERT_EQUAL_INT(0, cropped_src->n_comps);
  }
}


/*
This tests whether the az/za calculations are correct. As all example POINT,
GAUSSIAN, and SHAPELET sources have been given the same RA/Dec, they should
give the same results so feed in the correct COMPONENT array via functions
below
*/
void test_azza_coords(double lst_base, int num_comps, e_sky_crop sky_crop_type,
                      user_precision_t *calc_azs, user_precision_t *calc_zas){

  //Compare calculated az/za for all time steps to expectation
  //Some of these are small numbers so test within an accuracy
  for (size_t azza_ind = 0; azza_ind < num_comps*3; azza_ind++) {

    if (lst_base == LST0) {
      double *expec_az_LST0;
      double *expec_za_LST0;
      if (sky_crop_type == CROP_SOURCES){
        expec_az_LST0 = expec_az_sourcecrop_LST0;
        expec_za_LST0 = expec_za_sourcecrop_LST0;
      } else {
        expec_az_LST0 = expec_az_compcrop_LST0;
        expec_za_LST0 = expec_za_compcrop_LST0;
      }
      // printf("calc %.5f expec %.5f\n",calc_azs[azza_ind], expec_az_LST0[azza_ind] );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_az_LST0[azza_ind],
                               calc_azs[azza_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_za_LST0[azza_ind],
                               calc_zas[azza_ind]);
    } else if (lst_base == LST1) {
      // printf("calc %.5f expec %.5f\n",calc_azs[azza_ind], expec_az_LST1[azza_ind] );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_az_LST1[azza_ind],
                               calc_azs[azza_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_za_LST1[azza_ind],
                               calc_zas[azza_ind]);
    } else if (lst_base == LST2) {
      // printf("calc %.5f expec %.5f\n",calc_azs[azza_ind], expec_az_LST2[azza_ind] );
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_az_LST2[azza_ind],
                               calc_azs[azza_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_za_LST2[azza_ind],
                               calc_zas[azza_ind]);
    }
  }
}


/*
Checks the values match expected for given LST and COMPONENT type, when cropping
by SOURCE or by COMPONENT
*/
void check_component_values_after_crop(double lst_base, source_t *cropped_src,
                                       e_component_type comptype,
                                       e_sky_crop sky_crop_type) {

  double *comp0_ref_freqs = NULL;
  user_precision_t *comp0_ref_stokesI = NULL;
  user_precision_t *comp0_ref_stokesQ = NULL;
  user_precision_t *comp0_ref_stokesU = NULL;
  user_precision_t *comp0_ref_stokesV = NULL;
  user_precision_t *comp0_SIs = NULL;
  user_precision_t *comp0_curve_qs = NULL;
  int *comp0_power_inds = NULL;
  int *comp0_curve_inds = NULL;
  int *comp0_list_comp_inds = NULL;
  int *comp0_num_list_values = NULL;
  int *comp0_list_start_indexes = NULL;
  double *comp0_list_freqs = NULL;
  user_precision_t *comp0_list_stokesI = NULL;
  user_precision_t *comp0_list_stokesQ = NULL;
  user_precision_t *comp0_list_stokesU = NULL;
  user_precision_t *comp0_list_stokesV = NULL;
  user_precision_t *comp0_majors = NULL;
  user_precision_t *comp0_minors = NULL;
  user_precision_t *comp0_pas = NULL;


  double *comp1_ref_freqs = NULL;
  user_precision_t *comp1_ref_stokesI = NULL;
  user_precision_t *comp1_ref_stokesQ = NULL;
  user_precision_t *comp1_ref_stokesU = NULL;
  user_precision_t *comp1_ref_stokesV = NULL;
  user_precision_t *comp1_SIs = NULL;
  user_precision_t *comp1_curve_qs = NULL;
  int *comp1_power_inds = NULL;
  int *comp1_curve_inds = NULL;
  int *comp1_list_comp_inds = NULL;
  int *comp1_num_list_values = NULL;
  int *comp1_list_start_indexes = NULL;
  double *comp1_list_freqs = NULL;
  user_precision_t *comp1_list_stokesI = NULL;
  user_precision_t *comp1_list_stokesQ = NULL;
  user_precision_t *comp1_list_stokesU = NULL;
  user_precision_t *comp1_list_stokesV = NULL;
  user_precision_t *comp1_majors = NULL;
  user_precision_t *comp1_minors = NULL;
  user_precision_t *comp1_pas = NULL;

  components_t *components = NULL;


  //The total length of all list entires for 3 components in point0
  int len_point0_lists = 9;

  int len_gauss0_lists = 10;
  int len_gauss1_lists = 12;

  int len_shape0_lists = 8;
  int len_shape1_lists = 9;

  int len_comp0_lists = 0;
  int len_comp1_lists = 0;

  // int num_shape0_coeffs = 29;
  // int num_shape1_coeffs = 24;

  int n_comps = 0, n_powers = 0, n_curves = 0, n_lists = 0;

  if (comptype == POINT) {

    components = &cropped_src->point_components;
    n_comps = cropped_src->n_points;
    n_powers = cropped_src->n_point_powers;
    n_curves = cropped_src->n_point_curves;
    n_lists = cropped_src->n_point_lists;

    comp0_ref_freqs = point0_ref_freqs;
    comp0_ref_stokesI = point0_ref_stokesI;
    comp0_ref_stokesQ = point0_ref_stokesQ;
    comp0_ref_stokesU = point0_ref_stokesU;
    comp0_ref_stokesV = point0_ref_stokesV;
    comp0_SIs = point0_SIs;
    comp0_curve_qs = point0_curve_qs;
    comp0_power_inds = point0_power_inds;
    comp0_curve_inds = point0_curve_inds;
    comp0_list_comp_inds = point0_list_comp_inds;
    comp0_num_list_values = point0_num_list_values;
    comp0_list_start_indexes = point0_list_start_indexes;
    comp0_list_freqs = point0_list_freqs;
    comp0_list_stokesI = point0_list_stokesI;
    comp0_list_stokesQ = point0_list_stokesQ;
    comp0_list_stokesU = point0_list_stokesU;
    comp0_list_stokesV = point0_list_stokesV;
    len_comp0_lists = len_point0_lists;

    comp1_ref_freqs = point1_ref_freqs;
    comp1_ref_stokesI = point1_ref_stokesI;
    comp1_ref_stokesQ = point1_ref_stokesQ;
    comp1_ref_stokesU = point1_ref_stokesU;
    comp1_ref_stokesV = point1_ref_stokesV;
    comp1_SIs = point1_SIs;
    comp1_curve_qs = point1_curve_qs;
    comp1_power_inds = point1_power_inds;
    comp1_curve_inds = point1_curve_inds;
    comp1_list_comp_inds = point1_list_comp_inds;
    comp1_num_list_values = point1_num_list_values;
    comp1_list_start_indexes = point1_list_start_indexes;
    comp1_list_freqs = point1_list_freqs;
    comp1_list_stokesI = point1_list_stokesI;
    comp1_list_stokesQ = point1_list_stokesQ;
    comp1_list_stokesU = point1_list_stokesU;
    comp1_list_stokesV = point1_list_stokesV;
    len_comp1_lists = len_point0_lists;

  } else if (comptype == GAUSSIAN) {

    components = &cropped_src->gauss_components;
    n_comps = cropped_src->n_gauss;
    n_powers = cropped_src->n_gauss_powers;
    n_curves = cropped_src->n_gauss_curves;
    n_lists = cropped_src->n_gauss_lists;

    comp0_majors = gauss0_majors;
    comp0_minors = gauss0_minors;
    comp0_pas = gauss0_pas;
    comp0_ref_freqs = gauss0_ref_freqs;
    comp0_ref_stokesI = gauss0_ref_stokesI;
    comp0_ref_stokesQ = gauss0_ref_stokesQ;
    comp0_ref_stokesU = gauss0_ref_stokesU;
    comp0_ref_stokesV = gauss0_ref_stokesV;
    comp0_SIs = gauss0_SIs;
    comp0_curve_qs = gauss0_curve_qs;
    comp0_power_inds = gauss0_power_inds;
    comp0_curve_inds = gauss0_curve_inds;
    comp0_list_comp_inds = gauss0_list_comp_inds;
    comp0_num_list_values = gauss0_num_list_values;
    comp0_list_start_indexes = gauss0_list_start_indexes;
    comp0_list_freqs = gauss0_list_freqs;
    comp0_list_stokesI = gauss0_list_stokesI;
    comp0_list_stokesQ = gauss0_list_stokesQ;
    comp0_list_stokesU = gauss0_list_stokesU;
    comp0_list_stokesV = gauss0_list_stokesV;
    len_comp0_lists = len_gauss0_lists;

    comp1_majors = gauss1_majors;
    comp1_minors = gauss1_minors;
    comp1_pas = gauss1_pas;
    comp1_ref_freqs = gauss1_ref_freqs;
    comp1_ref_stokesI = gauss1_ref_stokesI;
    comp1_ref_stokesQ = gauss1_ref_stokesQ;
    comp1_ref_stokesU = gauss1_ref_stokesU;
    comp1_ref_stokesV = gauss1_ref_stokesV;
    comp1_SIs = gauss1_SIs;
    comp1_curve_qs = gauss1_curve_qs;
    comp1_power_inds = gauss1_power_inds;
    comp1_curve_inds = gauss1_curve_inds;
    comp1_list_comp_inds = gauss1_list_comp_inds;
    comp1_num_list_values = gauss1_num_list_values;
    comp1_list_start_indexes = gauss1_list_start_indexes;
    comp1_list_freqs = gauss1_list_freqs;
    comp1_list_stokesI = gauss1_list_stokesI;
    comp1_list_stokesQ = gauss1_list_stokesQ;
    comp1_list_stokesU = gauss1_list_stokesU;
    comp1_list_stokesV = gauss1_list_stokesV;
    len_comp1_lists = len_gauss1_lists;

  } else if (comptype == SHAPELET) {

    components = &cropped_src->shape_components;
    n_comps = cropped_src->n_shapes;
    n_powers = cropped_src->n_shape_powers;
    n_curves = cropped_src->n_shape_curves;
    n_lists = cropped_src->n_shape_lists;

    comp0_majors = shape0_majors;
    comp0_minors = shape0_minors;
    comp0_pas = shape0_pas;
    comp0_ref_freqs = shape0_ref_freqs;
    comp0_ref_stokesI = shape0_ref_stokesI;
    comp0_ref_stokesQ = shape0_ref_stokesQ;
    comp0_ref_stokesU = shape0_ref_stokesU;
    comp0_ref_stokesV = shape0_ref_stokesV;
    comp0_SIs = shape0_SIs;
    comp0_curve_qs = shape0_curve_qs;
    comp0_power_inds = shape0_power_inds;
    comp0_curve_inds = shape0_curve_inds;
    comp0_list_comp_inds = shape0_list_comp_inds;
    comp0_num_list_values = shape0_num_list_values;
    comp0_list_start_indexes = shape0_list_start_indexes;
    comp0_list_freqs = shape0_list_freqs;
    comp0_list_stokesI = shape0_list_stokesI;
    comp0_list_stokesQ = shape0_list_stokesQ;
    comp0_list_stokesU = shape0_list_stokesU;
    comp0_list_stokesV = shape0_list_stokesV;
    len_comp0_lists = len_shape0_lists;

    comp1_majors = shape1_majors;
    comp1_minors = shape1_minors;
    comp1_pas = shape1_pas;
    comp1_ref_freqs = shape1_ref_freqs;
    comp1_ref_stokesI = shape1_ref_stokesI;
    comp1_ref_stokesQ = shape1_ref_stokesQ;
    comp1_ref_stokesU = shape1_ref_stokesU;
    comp1_ref_stokesV = shape1_ref_stokesV;
    comp1_SIs = shape1_SIs;
    comp1_curve_qs = shape1_curve_qs;
    comp1_power_inds = shape1_power_inds;
    comp1_curve_inds = shape1_curve_inds;
    comp1_list_comp_inds = shape1_list_comp_inds;
    comp1_num_list_values = shape1_num_list_values;
    comp1_list_start_indexes = shape1_list_start_indexes;
    comp1_list_freqs = shape1_list_freqs;
    comp1_list_stokesI = shape1_list_stokesI;
    comp1_list_stokesQ = shape1_list_stokesQ;
    comp1_list_stokesU = shape1_list_stokesU;
    comp1_list_stokesV = shape1_list_stokesV;
    len_comp1_lists = len_shape1_lists;

  }

  //Test if the az/za coords have been calculated correctly
  test_azza_coords(lst_base, n_comps, sky_crop_type,
                   components->azs, components->zas);

  int num_comps = 0;
  int num_powers = 0;
  int num_lists = 0;
  int num_curves = 0;

  //For LST0, all COMPONENTs of the first POINT SOURCE
  //(as defined in test_crop_sky_model.h) should have been retained.
  //Only the first COMPONENT of the second SOURCE is above the horizon, so
  //if cropping by SOURCE there should only be three COMPONENTs retained.
  //If cropping by COMPONENT, there should be a 4th component
  if (lst_base == LST0) {

    num_comps = 9;
    num_powers = 3;
    num_lists = 3;
    num_curves = 3;

    //Here we are testing up the array up until num_comps, so the test covers
    //the first set of components for both the SOURCE and COMPONENT sky cropping cases
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_ras, components->ras, num_comps);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_decs, components->decs, num_comps);

    //power law stuff-----------------------------------------------------------
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_ref_freqs, components->power_ref_freqs, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesI,
                                  components->power_ref_stokesI, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesQ,
                                  components->power_ref_stokesQ, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesU,
                                  components->power_ref_stokesU, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesV,
                                  components->power_ref_stokesV, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_SIs, components->power_SIs, num_powers);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_power_inds,
            components->power_comp_inds, num_powers);

    //curved power law stuff----------------------------------------------------
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_ref_freqs, components->curve_ref_freqs, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesI,
                                  components->curve_ref_stokesI, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesQ,
                                  components->curve_ref_stokesQ, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesU,
                                  components->curve_ref_stokesU, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesV,
                                  components->curve_ref_stokesV, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_SIs, components->curve_SIs, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_curve_qs, components->curve_qs, num_curves);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_curve_inds,
            components->curve_comp_inds, num_curves);

    //list flux type stuff------------------------------------------------------

    //indexing arrays, test up to number of point list components
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_list_comp_inds,
            components->list_comp_inds, num_lists);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_num_list_values,
            components->num_list_values, num_lists);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_list_start_indexes,
            components->list_start_indexes, num_lists);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_list_freqs,
            components->list_freqs, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesI,
            components->list_stokesI, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesQ,
            components->list_stokesQ, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesU,
            components->list_stokesU, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesV,
            components->list_stokesV, len_comp0_lists);

    if (comptype == GAUSSIAN || comptype == SHAPELET) {
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_majors,
              components->majors, num_comps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_minors,
              components->minors, num_comps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_pas,
              components->pas, num_comps);
    }

    if (comptype == SHAPELET) {
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_coeffs, components->shape_coeffs,
                                                             NUM_SHAPE0_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n1s, components->n1s,
                                                             NUM_SHAPE0_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n2s, components->n2s,
                                                             NUM_SHAPE0_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_param_indexes,
                                  components->param_indexes, NUM_SHAPE0_COEFFS);
    }

    //

    if (sky_crop_type == CROP_SOURCES) {
      TEST_ASSERT_EQUAL_INT(num_comps, n_comps);
      TEST_ASSERT_EQUAL_INT(num_powers, n_powers);
      TEST_ASSERT_EQUAL_INT(num_curves, n_curves);
      TEST_ASSERT_EQUAL_INT(num_lists, n_lists);

      if (comptype == SHAPELET){
        TEST_ASSERT_EQUAL_INT(NUM_SHAPE0_COEFFS, cropped_src->n_shape_coeffs);
      }

    } else if (sky_crop_type == CROP_COMPONENTS) {
      num_comps = 12;
      num_powers = 4;
      num_lists = 4;
      num_curves = 4;
      TEST_ASSERT_EQUAL_INT(num_comps, n_comps);
      TEST_ASSERT_EQUAL_INT(num_powers, n_powers);
      TEST_ASSERT_EQUAL_INT(num_curves, n_curves);
      TEST_ASSERT_EQUAL_INT(num_lists, n_lists);

      //test that 10, 11, 12th point COMPONENT ra, dec is equal to the first of comp1 ra,dec
      //Works because I've repeated three ra,dec coords pairs three times
      TEST_ASSERT_EQUAL_DOUBLE(comp1_ras[0], components->ras[9]);
      TEST_ASSERT_EQUAL_DOUBLE(comp1_decs[0], components->decs[9]);
      TEST_ASSERT_EQUAL_DOUBLE(comp1_ras[0], components->ras[10]);
      TEST_ASSERT_EQUAL_DOUBLE(comp1_decs[0], components->decs[10]);
      TEST_ASSERT_EQUAL_DOUBLE(comp1_ras[0], components->ras[11]);
      TEST_ASSERT_EQUAL_DOUBLE(comp1_decs[0], components->decs[11]);


      if (comptype == GAUSSIAN || comptype == SHAPELET) {
        for (int comp_ind = 9; comp_ind < 12; comp_ind++) {
          TEST_ASSERT_EQUAL_FLOAT(comp1_majors[0], components->majors[comp_ind]);
          TEST_ASSERT_EQUAL_FLOAT(comp1_minors[0], components->minors[comp_ind]);
          TEST_ASSERT_EQUAL_FLOAT(comp1_pas[0], components->pas[comp_ind]);
        }
      }


      //Then the 4th power-law component should be equal to the first of comp1 power law values
      int power_ind = 3;
      TEST_ASSERT_EQUAL_DOUBLE(comp1_ref_freqs[0], components->power_ref_freqs[power_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesI[0],
                             components->power_ref_stokesI[power_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesQ[0],
                             components->power_ref_stokesQ[power_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesU[0],
                             components->power_ref_stokesU[power_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesV[0],
                             components->power_ref_stokesV[power_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_SIs[0], components->power_SIs[power_ind]);
      //Should be the 10th component
      TEST_ASSERT_EQUAL_INT(9, components->power_comp_inds[power_ind]);

      //Similarly for the curved power law values
      int curve_ind = 3;
      TEST_ASSERT_EQUAL_DOUBLE(comp1_ref_freqs[0], components->curve_ref_freqs[curve_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesI[0],
                             components->curve_ref_stokesI[curve_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesQ[0],
                             components->curve_ref_stokesQ[curve_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesU[0],
                             components->curve_ref_stokesU[curve_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_ref_stokesV[0],
                             components->curve_ref_stokesV[curve_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_SIs[0], components->curve_SIs[curve_ind]);
      TEST_ASSERT_EQUAL_FLOAT(comp1_curve_qs[0], components->curve_qs[curve_ind]);
      TEST_ASSERT_EQUAL_INT(10, components->curve_comp_inds[curve_ind]);

      //OK, the 4th list type flux component should be equal to the first
      //component of the second source. The indexing arrays change however
      //as they are based on the first three components
      //This should be the 12th component in the cropped source, so
      //list_comp_ind should be the 11th index
      int list_ind = 3;
      //indexing arrays, test up to number of point list components
      TEST_ASSERT_EQUAL_INT(11,
              components->list_comp_inds[list_ind]);
      TEST_ASSERT_EQUAL_INT(comp1_num_list_values[0],
              components->num_list_values[list_ind]);
      TEST_ASSERT_EQUAL_INT(comp1_list_start_indexes[0] + len_comp0_lists,
              components->list_start_indexes[list_ind]);

      //Ok, the amount of flux list entries we should have in this cropped
      //source should be given by the number in the first component of the second
      //source, which is listed in comp1_num_list_values
      TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_list_freqs,
              components->list_freqs + len_comp0_lists, comp1_num_list_values[0]);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesI,
              components->list_stokesI + len_comp0_lists, comp1_num_list_values[0]);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesQ,
              components->list_stokesQ + len_comp0_lists, comp1_num_list_values[0]);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesU,
              components->list_stokesU + len_comp0_lists, comp1_num_list_values[0]);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesV,
              components->list_stokesV + len_comp0_lists, comp1_num_list_values[0]);


      //Ok, the first of each flux type SHAPELET component should have
      //made it through. This results in a complicated mapping from
      //the second SOURCE, so there will be hard-coding here - basically
      //pulling out everything from the arrays in test_crop_sky_model.h
      //that match the indexes 0, 3, 6 in shape1_param_indexes
      if (comptype == SHAPELET){
        TEST_ASSERT_EQUAL_INT(cropped_src->n_shape_coeffs,
                              NUM_SHAPE0_COEFFS + 9);

        user_precision_t expected_param_indexes[] = {9, 9, 9,
                                                      10, 10, 10, 10, 10, 11};
        TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_param_indexes,
                components->param_indexes + NUM_SHAPE0_COEFFS, 9);

        user_precision_t expected_coeffs[] = {0.609, 0.695, 0.561,
        0.628, 0.06, 0.53, 0.772, 0.253, 0.805};
        TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_coeffs,
                components->shape_coeffs + NUM_SHAPE0_COEFFS, 9);

        user_precision_t expected_n1s[] = {18, 24, 81, 50, 88, 95, 63, 33, 89};
        TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_n1s,
                components->n1s + NUM_SHAPE0_COEFFS, 9);

        user_precision_t expected_n2s[] = {51, 47, 10, 45,  5,  8, 90, 44, 55};
        TEST_ASSERT_EQUAL_FLOAT_ARRAY(expected_n2s,
                components->n2s + NUM_SHAPE0_COEFFS, 9);
      }
    }
  }
  else if (lst_base == LST1) { //For this LST, both SOURCEs should be retained
    num_comps = 18;
    num_powers = 6;
    num_lists = 6;
    num_curves = 6;

    int radec_len = 9;
    int freq_len = 3;


    TEST_ASSERT_EQUAL_INT(num_comps, n_comps);
    TEST_ASSERT_EQUAL_INT(num_powers, n_powers);
    TEST_ASSERT_EQUAL_INT(num_curves, n_curves);
    TEST_ASSERT_EQUAL_INT(num_lists, n_lists);

    //First three COMPONENTs should match comp0
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_ras, components->ras, radec_len);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_decs, components->decs, radec_len);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_ref_freqs, components->power_ref_freqs, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesI,
                                  components->power_ref_stokesI, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesQ,
                                  components->power_ref_stokesQ, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesU,
                                  components->power_ref_stokesU, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesV,
                                  components->power_ref_stokesV, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_SIs, components->power_SIs, freq_len);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_ref_freqs, components->curve_ref_freqs, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesI,
                                  components->curve_ref_stokesI, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesQ,
                                  components->curve_ref_stokesQ, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesU,
                                  components->curve_ref_stokesU, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_ref_stokesV,
                                  components->curve_ref_stokesV, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_SIs, components->curve_SIs, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_curve_qs, components->curve_qs, freq_len);

    //indexing arrays, test up to number of point list components
    //Only testing that first half of the indexing arrays match comp0
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_list_comp_inds,
            components->list_comp_inds, num_lists/2);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_num_list_values,
            components->num_list_values, num_lists/2);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp0_list_start_indexes,
            components->list_start_indexes, num_lists/2);

    //The total length of all list entires for 3 components in comp0
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp0_list_freqs,
            components->list_freqs, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesI,
            components->list_stokesI, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesQ,
            components->list_stokesQ, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesU,
            components->list_stokesU, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_list_stokesV,
            components->list_stokesV, len_comp0_lists);

    //Last three COMPONENTs of each flux type should match comp1
    //Gooooooooooo pointer arithmatic
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_ras, components->ras + radec_len, radec_len);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_decs, components->decs + radec_len, radec_len);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_ref_freqs, components->power_ref_freqs + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesI,
                                  components->power_ref_stokesI + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesQ,
                                  components->power_ref_stokesQ + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesU,
                                  components->power_ref_stokesU + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesV,
                                  components->power_ref_stokesV + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_SIs, components->power_SIs + freq_len, freq_len);

    for (int power_ind = 0; power_ind < 3; power_ind++) {
      TEST_ASSERT_EQUAL_INT(comp1_power_inds[power_ind] + radec_len,
                                        components->power_comp_inds[power_ind + 3]);
    }





    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_ref_freqs, components->curve_ref_freqs + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesI,
                                  components->curve_ref_stokesI + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesQ,
                                  components->curve_ref_stokesQ + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesU,
                                  components->curve_ref_stokesU + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesV,
                                  components->curve_ref_stokesV + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_SIs, components->curve_SIs + freq_len, freq_len);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_curve_qs, components->curve_qs + freq_len, freq_len);
    for (int curve_ind = 0; curve_ind < 3; curve_ind++) {
      TEST_ASSERT_EQUAL_INT(comp1_curve_inds[curve_ind] + radec_len,
                                        components->curve_comp_inds[curve_ind + 3]);
    }

    //indexing arrays, test that second half matches values of the second
    //source, along with the fact there are now 9 components from the first
    //source in this cropped source already

    for (int list_ind = 0; list_ind < 3; list_ind++) {
      TEST_ASSERT_EQUAL_INT(comp1_list_comp_inds[list_ind] + radec_len,
              components->list_comp_inds[list_ind + 3]);
      TEST_ASSERT_EQUAL_INT(comp1_num_list_values[list_ind],
              components->num_list_values[list_ind + 3]);
      TEST_ASSERT_EQUAL_INT(comp1_list_start_indexes[list_ind] + len_comp0_lists,
              components->list_start_indexes[list_ind + 3]);
    }

    //Check second half of the cropped model matches comp1 list entries
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_list_freqs,
            components->list_freqs + len_comp0_lists, len_comp1_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesI,
            components->list_stokesI + len_comp0_lists, len_comp1_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesQ,
            components->list_stokesQ + len_comp0_lists, len_comp1_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesU,
            components->list_stokesU + len_comp0_lists, len_comp1_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesV,
            components->list_stokesV + len_comp0_lists, len_comp1_lists);

    if (comptype == GAUSSIAN || comptype == SHAPELET) {
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_majors,
              components->majors, radec_len);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_minors,
              components->minors, radec_len);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp0_pas,
              components->pas, radec_len);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_majors,
              components->majors + radec_len, radec_len);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_minors,
              components->minors + radec_len, radec_len);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_pas,
              components->pas + radec_len, radec_len);
    }

    if (comptype == SHAPELET) {
      TEST_ASSERT_EQUAL_INT(cropped_src->n_shape_coeffs,
                            NUM_SHAPE0_COEFFS + NUM_SHAPE1_COEFFS);
      //First half should match source0 stuff
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_coeffs, components->shape_coeffs,
                                                             NUM_SHAPE0_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n1s, components->n1s,
                                                             NUM_SHAPE0_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n2s, components->n2s,
                                                             NUM_SHAPE0_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_param_indexes,
                                  components->param_indexes, NUM_SHAPE0_COEFFS);

      //Second half should match apart from the indexes that should ne higher
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_coeffs,
               components->shape_coeffs + NUM_SHAPE0_COEFFS, NUM_SHAPE1_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n1s,
               components->n1s + NUM_SHAPE0_COEFFS, NUM_SHAPE1_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n2s,
               components->n2s + NUM_SHAPE0_COEFFS, NUM_SHAPE1_COEFFS);

      for (int coeff = 0; coeff < NUM_SHAPE1_COEFFS; coeff++) {
        TEST_ASSERT_EQUAL_FLOAT(shape1_param_indexes[coeff] + radec_len,
                          components->param_indexes[NUM_SHAPE0_COEFFS + coeff] );
      }
    }


  } else if (lst_base == LST2) { //first source should have been cropped, second retained

    num_comps = 9;
    num_powers = 3;
    num_lists = 3;
    num_curves = 3;

    TEST_ASSERT_EQUAL_INT(num_comps, n_comps);
    TEST_ASSERT_EQUAL_INT(num_powers, n_powers);
    TEST_ASSERT_EQUAL_INT(num_curves, n_curves);
    TEST_ASSERT_EQUAL_INT(num_lists, n_lists);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_ras, components->ras, num_powers);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_decs, components->decs, num_powers);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_ref_freqs, components->power_ref_freqs, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesI,
                                  components->power_ref_stokesI, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesQ,
                                  components->power_ref_stokesQ, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesU,
                                  components->power_ref_stokesU, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesV,
                                  components->power_ref_stokesV, num_powers);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_SIs, components->power_SIs, num_powers);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp1_power_inds, components->power_comp_inds, num_powers);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_ref_freqs, components->curve_ref_freqs, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesI,
                                  components->curve_ref_stokesI, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesQ,
                                  components->curve_ref_stokesQ, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesU,
                                  components->curve_ref_stokesU, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_ref_stokesV,
                                  components->curve_ref_stokesV, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_SIs, components->curve_SIs, num_curves);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_curve_qs, components->curve_qs, num_curves);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp1_curve_inds, components->curve_comp_inds, num_curves);

    //indexing arrays, test up to number of point list components
    TEST_ASSERT_EQUAL_INT_ARRAY(comp1_list_comp_inds,
            components->list_comp_inds, num_lists);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp1_num_list_values,
            components->num_list_values, num_lists);
    TEST_ASSERT_EQUAL_INT_ARRAY(comp1_list_start_indexes,
            components->list_start_indexes, num_lists);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(comp1_list_freqs,
            components->list_freqs, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesI,
            components->list_stokesI, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesQ,
            components->list_stokesQ, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesU,
            components->list_stokesU, len_comp0_lists);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_list_stokesV,
            components->list_stokesV, len_comp0_lists);

    if (comptype == GAUSSIAN || comptype == SHAPELET) {
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_majors, components->majors, num_comps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_minors, components->minors, num_comps);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(comp1_pas, components->pas, num_comps);
    }

    if (comptype == SHAPELET) {
      TEST_ASSERT_EQUAL_INT(cropped_src->n_shape_coeffs, NUM_SHAPE1_COEFFS);

      //First half should match source0 stuff
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_coeffs, components->shape_coeffs,
                                                             NUM_SHAPE1_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n1s, components->n1s,
                                                             NUM_SHAPE1_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n2s, components->n2s,
                                                             NUM_SHAPE1_COEFFS);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_param_indexes,
                                  components->param_indexes, NUM_SHAPE1_COEFFS);
    }

  } else if (lst_base == LST3) { //both should be flagged
    TEST_ASSERT_EQUAL_INT(num_comps, n_comps);
    TEST_ASSERT_EQUAL_INT(num_powers, n_powers);
    TEST_ASSERT_EQUAL_INT(num_curves, n_curves);
    TEST_ASSERT_EQUAL_INT(num_lists, n_lists);
  }
}


void test_crop_sky_model_Point(double lst_base, e_sky_crop sky_crop_type) {

  //Just use POINT type COMPONENTs
  int point = 1;
  int gauss = 0;
  int shape = 0;

  //Setup the lsts - the first of this list is used in calculating az/za
  //for cropping. Once cropped, az/za is calculated for all surviving COMPONENTs
  int num_time_steps = 3;
  double lsts[] = {lst_base, lst_base+60*SOLAR2SIDEREAL*DS2R,
                  lst_base + 120.0*SOLAR2SIDEREAL*DS2R};

  //Create the input sky catalogue
  source_catalogue_t *raw_srccat = create_input_srccat(point, gauss, shape);

  //Call the function being tested
  source_t *cropped_src;
  cropped_src =  crop_sky_model(raw_srccat, lsts, (double)MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  check_component_values_after_crop(lst_base, cropped_src, POINT, sky_crop_type);

  check_num_components_after_crop(lst_base, cropped_src,
                                  sky_crop_type,
                                  point, gauss, shape);


}

void test_crop_sky_model_PointCropSourceLST0(void) {
  test_crop_sky_model_Point(LST0, CROP_SOURCES);
}

void test_crop_sky_model_PointCropSourceLST1(void) {
  test_crop_sky_model_Point(LST1, CROP_SOURCES);
}

void test_crop_sky_model_PointCropSourceLST2(void) {
  test_crop_sky_model_Point(LST2, CROP_SOURCES);
}

void test_crop_sky_model_PointCropSourceLST3(void) {
  test_crop_sky_model_Point(LST3, CROP_SOURCES);
}

void test_crop_sky_model_PointCropComponentLST0(void) {
  test_crop_sky_model_Point(LST0, CROP_COMPONENTS);
}

void test_crop_sky_model_PointCropComponentLST1(void) {
  test_crop_sky_model_Point(LST1, CROP_COMPONENTS);
}

void test_crop_sky_model_PointCropComponentLST2(void) {
  test_crop_sky_model_Point(LST2, CROP_COMPONENTS);
}

void test_crop_sky_model_PointCropComponentLST3(void) {
  test_crop_sky_model_Point(LST3, CROP_COMPONENTS);
}


void test_crop_sky_model_Gauss(double lst_base, e_sky_crop sky_crop_type) {

  //Just use GAUSSIAN type COMPONENTs
  int point = 0;
  int gauss = 1;
  int shape = 0;

  //Setup the lsts - the first of this list is used in calculating az/za
  //for cropping. Once cropped, az/za is calculated for all surviving COMPONENTs
  int num_time_steps = 3;
  double lsts[] = {lst_base, lst_base+60*SOLAR2SIDEREAL*DS2R,
                  lst_base + 120.0*SOLAR2SIDEREAL*DS2R};

  //Create the input sky catalogue
  source_catalogue_t *raw_srccat = create_input_srccat(point, gauss, shape);

  //Call the function being tested
  source_t *cropped_src;
  cropped_src =  crop_sky_model(raw_srccat, lsts, MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  check_component_values_after_crop(lst_base, cropped_src, GAUSSIAN, sky_crop_type);

  check_num_components_after_crop(lst_base, cropped_src,
                                  sky_crop_type,
                                  point, gauss, shape);


}

void test_crop_sky_model_GaussCropSourceLST0(void) {
  test_crop_sky_model_Gauss(LST0, CROP_SOURCES);
}

void test_crop_sky_model_GaussCropSourceLST1(void) {
  test_crop_sky_model_Gauss(LST1, CROP_SOURCES);
}

void test_crop_sky_model_GaussCropSourceLST2(void) {
  test_crop_sky_model_Gauss(LST2, CROP_SOURCES);
}

void test_crop_sky_model_GaussCropSourceLST3(void) {
  test_crop_sky_model_Gauss(LST3, CROP_SOURCES);
}

void test_crop_sky_model_GaussCropComponentLST0(void) {
  test_crop_sky_model_Gauss(LST0, CROP_COMPONENTS);
}

void test_crop_sky_model_GaussCropComponentLST1(void) {
  test_crop_sky_model_Gauss(LST1, CROP_COMPONENTS);
}

void test_crop_sky_model_GaussCropComponentLST2(void) {
  test_crop_sky_model_Gauss(LST2, CROP_COMPONENTS);
}

void test_crop_sky_model_GaussCropComponentLST3(void) {
  test_crop_sky_model_Gauss(LST3, CROP_COMPONENTS);
}

/*******************************************************************************
---------------SHAPELET SOURCE functions------------------------------------------
*******************************************************************************/


void test_crop_sky_model_Shape(double lst_base, e_sky_crop sky_crop_type) {

  //Just use SHAPELET type COMPONENTs
  int point = 0;
  int gauss = 0;
  int shape = 1;

  //Setup the lsts - the first of this list is used in calculating az/za
  //for cropping. Once cropped, az/za is calculated for all surviving COMPONENTs
  int num_time_steps = 3;
  double lsts[] = {lst_base, lst_base+60*SOLAR2SIDEREAL*DS2R,
                  lst_base + 120.0*SOLAR2SIDEREAL*DS2R};

  //Create the input sky catalogue
  source_catalogue_t *raw_srccat = create_input_srccat(point, gauss, shape);

  //Call the function being tested
  source_t *cropped_src;
  cropped_src = crop_sky_model(raw_srccat, lsts, (double)MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  // check_shape_values_after_crop(lst_base, cropped_src, sky_crop_type);
  check_component_values_after_crop(lst_base, cropped_src, SHAPELET, sky_crop_type);

  check_num_components_after_crop(lst_base, cropped_src,
                                  sky_crop_type,
                                  point, gauss, shape);


}

void test_crop_sky_model_ShapeCropSourceLST0(void) {
  test_crop_sky_model_Shape(LST0, CROP_SOURCES);
}

void test_crop_sky_model_ShapeCropSourceLST1(void) {
  test_crop_sky_model_Shape(LST1, CROP_SOURCES);
}

void test_crop_sky_model_ShapeCropSourceLST2(void) {
  test_crop_sky_model_Shape(LST2, CROP_SOURCES);
}

void test_crop_sky_model_ShapeCropSourceLST3(void) {
  test_crop_sky_model_Shape(LST3, CROP_SOURCES);
}

void test_crop_sky_model_ShapeCropComponentLST0(void) {
  test_crop_sky_model_Shape(LST0, CROP_COMPONENTS);
}

void test_crop_sky_model_ShapeCropComponentLST1(void) {
  test_crop_sky_model_Shape(LST1, CROP_COMPONENTS);
}

void test_crop_sky_model_ShapeCropComponentLST2(void) {
  test_crop_sky_model_Shape(LST2, CROP_COMPONENTS);
}

void test_crop_sky_model_ShapeCropComponentLST3(void) {
  test_crop_sky_model_Shape(LST3, CROP_COMPONENTS);
}


/*******************************************************************************
------TEST ALL THREE COMPONENT TYPES AT ONCE
*******************************************************************************/
void test_crop_sky_model_AllTypes(double lst_base, e_sky_crop sky_crop_type) {

  int point = 1;
  int gauss = 1;
  int shape = 1;

  //Setup the lsts - the first of this list is used in calculating az/za
  //for cropping. Once cropped, az/za is calculated for all surviving COMPONENTs
  int num_time_steps = 3;
  double lsts[] = {lst_base, lst_base+60*SOLAR2SIDEREAL*DS2R,
                  lst_base + 120.0*SOLAR2SIDEREAL*DS2R};

  //Create the input sky catalogue
  source_catalogue_t *raw_srccat = create_input_srccat(point, gauss, shape);

  //Call the function being tested
  source_t *cropped_src;
  cropped_src = crop_sky_model(raw_srccat, lsts, (double)MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  // check_point_values_after_crop(lst_base, cropped_src, sky_crop_type);
  check_component_values_after_crop(lst_base, cropped_src, POINT, sky_crop_type);
  check_component_values_after_crop(lst_base, cropped_src, GAUSSIAN, sky_crop_type);
  check_component_values_after_crop(lst_base, cropped_src, SHAPELET, sky_crop_type);

  check_num_components_after_crop(lst_base, cropped_src,
                                  sky_crop_type,
                                  point, gauss, shape);


}

void test_crop_sky_model_AllTypesCropSourceLST0(void) {
  test_crop_sky_model_AllTypes(LST0, CROP_SOURCES);
}

void test_crop_sky_model_AllTypesCropSourceLST1(void) {
  test_crop_sky_model_AllTypes(LST1, CROP_SOURCES);
}

void test_crop_sky_model_AllTypesCropSourceLST2(void) {
  test_crop_sky_model_AllTypes(LST2, CROP_SOURCES);
}

void test_crop_sky_model_AllTypesCropSourceLST3(void) {
  test_crop_sky_model_AllTypes(LST3, CROP_SOURCES);
}

void test_crop_sky_model_AllTypesCropComponentLST0(void) {
  test_crop_sky_model_AllTypes(LST0, CROP_COMPONENTS);
}

void test_crop_sky_model_AllTypesCropComponentLST1(void) {
  test_crop_sky_model_AllTypes(LST1, CROP_COMPONENTS);
}

void test_crop_sky_model_AllTypesCropComponentLST2(void) {
  test_crop_sky_model_AllTypes(LST2, CROP_COMPONENTS);
}

void test_crop_sky_model_AllTypesCropComponentLST3(void) {
  test_crop_sky_model_AllTypes(LST3, CROP_COMPONENTS);
}


//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_crop_sky_model_PointCropSourceLST0);
    RUN_TEST(test_crop_sky_model_PointCropSourceLST1);
    RUN_TEST(test_crop_sky_model_PointCropSourceLST2);
    RUN_TEST(test_crop_sky_model_PointCropSourceLST3);
    RUN_TEST(test_crop_sky_model_PointCropComponentLST0);
    RUN_TEST(test_crop_sky_model_PointCropComponentLST1);
    RUN_TEST(test_crop_sky_model_PointCropComponentLST2);
    RUN_TEST(test_crop_sky_model_PointCropComponentLST3);

    RUN_TEST(test_crop_sky_model_GaussCropSourceLST0);
    RUN_TEST(test_crop_sky_model_GaussCropSourceLST1);
    RUN_TEST(test_crop_sky_model_GaussCropSourceLST2);
    RUN_TEST(test_crop_sky_model_GaussCropSourceLST3);
    RUN_TEST(test_crop_sky_model_GaussCropComponentLST0);
    RUN_TEST(test_crop_sky_model_GaussCropComponentLST1);
    RUN_TEST(test_crop_sky_model_GaussCropComponentLST2);
    RUN_TEST(test_crop_sky_model_GaussCropComponentLST3);

    RUN_TEST(test_crop_sky_model_ShapeCropSourceLST0);
    RUN_TEST(test_crop_sky_model_ShapeCropSourceLST1);
    RUN_TEST(test_crop_sky_model_ShapeCropSourceLST2);
    RUN_TEST(test_crop_sky_model_ShapeCropSourceLST3);
    RUN_TEST(test_crop_sky_model_ShapeCropComponentLST0);
    RUN_TEST(test_crop_sky_model_ShapeCropComponentLST1);
    RUN_TEST(test_crop_sky_model_ShapeCropComponentLST2);
    RUN_TEST(test_crop_sky_model_ShapeCropComponentLST3);

    RUN_TEST(test_crop_sky_model_AllTypesCropSourceLST0);
    RUN_TEST(test_crop_sky_model_AllTypesCropSourceLST1);
    RUN_TEST(test_crop_sky_model_AllTypesCropSourceLST2);
    RUN_TEST(test_crop_sky_model_AllTypesCropSourceLST3);
    RUN_TEST(test_crop_sky_model_AllTypesCropComponentLST0);
    RUN_TEST(test_crop_sky_model_AllTypesCropComponentLST1);
    RUN_TEST(test_crop_sky_model_AllTypesCropComponentLST2);
    RUN_TEST(test_crop_sky_model_AllTypesCropComponentLST3);

    return UNITY_END();
}
