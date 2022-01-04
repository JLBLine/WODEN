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
  raw_srccat->catsources = malloc(2*sizeof(catsource_t));

  //set the COMPONENT counters to zero for both SOURCEs

  for (int src = 0; src < 2; src++) {
    raw_srccat->catsources[src].n_comps = 0;
    raw_srccat->catsources[src].n_points= 0;
    raw_srccat->catsources[src].n_gauss = 0;
    raw_srccat->catsources[src].n_shapes = 0;
    raw_srccat->catsources[src].n_shape_coeffs = 0;
  }

  if (point) {
    //Assign the first SOURCE some POINT COMPONENTs
    raw_srccat->catsources[0].point_ras = point0_ras;
    raw_srccat->catsources[0].point_decs = point0_decs;
    raw_srccat->catsources[0].point_ref_freqs = point0_ref_freqs;
    raw_srccat->catsources[0].point_ref_stokesI = point0_ref_stokesI;
    raw_srccat->catsources[0].point_ref_stokesQ = point0_ref_stokesQ;
    raw_srccat->catsources[0].point_ref_stokesU = point0_ref_stokesU;
    raw_srccat->catsources[0].point_ref_stokesV = point0_ref_stokesV;
    raw_srccat->catsources[0].point_SIs = point0_SIs;

    raw_srccat->catsources[0].n_comps += 3;
    raw_srccat->catsources[0].n_points += 3;

    //Assign the second SOURCE some POINT COMPONENTs
    raw_srccat->catsources[1].point_ras = point1_ras;
    raw_srccat->catsources[1].point_decs = point1_decs;
    raw_srccat->catsources[1].point_ref_freqs = point1_ref_freqs;
    raw_srccat->catsources[1].point_ref_stokesI = point1_ref_stokesI;
    raw_srccat->catsources[1].point_ref_stokesQ = point1_ref_stokesQ;
    raw_srccat->catsources[1].point_ref_stokesU = point1_ref_stokesU;
    raw_srccat->catsources[1].point_ref_stokesV = point1_ref_stokesV;
    raw_srccat->catsources[1].point_SIs = point1_SIs;

    raw_srccat->catsources[1].n_comps += 3;
    raw_srccat->catsources[1].n_points += 3;

  }

  if (gauss) {
    //Assign the first SOURCE some GAUSSIAN COMPONENTs
    raw_srccat->catsources[0].gauss_ras = gauss0_ras;
    raw_srccat->catsources[0].gauss_decs = gauss0_decs;
    raw_srccat->catsources[0].gauss_ref_freqs = gauss0_ref_freqs;
    raw_srccat->catsources[0].gauss_ref_stokesI = gauss0_ref_stokesI;
    raw_srccat->catsources[0].gauss_ref_stokesQ = gauss0_ref_stokesQ;
    raw_srccat->catsources[0].gauss_ref_stokesU = gauss0_ref_stokesU;
    raw_srccat->catsources[0].gauss_ref_stokesV = gauss0_ref_stokesV;
    raw_srccat->catsources[0].gauss_SIs = gauss0_SIs;
    raw_srccat->catsources[0].gauss_majors = gauss0_majors;
    raw_srccat->catsources[0].gauss_minors = gauss0_minors;
    raw_srccat->catsources[0].gauss_pas = gauss0_pas;

    raw_srccat->catsources[0].n_comps += 3;
    raw_srccat->catsources[0].n_gauss += 3;

    //Assign the second SOURCE some GAUSSIAN COMPONENTs
    raw_srccat->catsources[1].gauss_ras = gauss1_ras;
    raw_srccat->catsources[1].gauss_decs = gauss1_decs;
    raw_srccat->catsources[1].gauss_ref_freqs = gauss1_ref_freqs;
    raw_srccat->catsources[1].gauss_ref_stokesI = gauss1_ref_stokesI;
    raw_srccat->catsources[1].gauss_ref_stokesQ = gauss1_ref_stokesQ;
    raw_srccat->catsources[1].gauss_ref_stokesU = gauss1_ref_stokesU;
    raw_srccat->catsources[1].gauss_ref_stokesV = gauss1_ref_stokesV;
    raw_srccat->catsources[1].gauss_SIs = gauss1_SIs;
    raw_srccat->catsources[1].gauss_majors = gauss1_majors;
    raw_srccat->catsources[1].gauss_minors = gauss1_minors;
    raw_srccat->catsources[1].gauss_pas = gauss1_pas;

    raw_srccat->catsources[1].n_comps += 3;
    raw_srccat->catsources[1].n_gauss += 3;

  }

  if (shape) {
    //Assign the first SOURCE some shapeIAN COMPONENTs
    raw_srccat->catsources[0].shape_ras = shape0_ras;
    raw_srccat->catsources[0].shape_decs = shape0_decs;
    raw_srccat->catsources[0].shape_ref_freqs = shape0_ref_freqs;
    raw_srccat->catsources[0].shape_ref_stokesI = shape0_ref_stokesI;
    raw_srccat->catsources[0].shape_ref_stokesQ = shape0_ref_stokesQ;
    raw_srccat->catsources[0].shape_ref_stokesU = shape0_ref_stokesU;
    raw_srccat->catsources[0].shape_ref_stokesV = shape0_ref_stokesV;
    raw_srccat->catsources[0].shape_SIs = shape0_SIs;
    raw_srccat->catsources[0].shape_majors = shape0_majors;
    raw_srccat->catsources[0].shape_minors = shape0_minors;
    raw_srccat->catsources[0].shape_pas = shape0_pas;
    raw_srccat->catsources[0].shape_coeffs = shape0_coeffs;
    raw_srccat->catsources[0].shape_n1s = shape0_n1s;
    raw_srccat->catsources[0].shape_n2s = shape0_n2s;
    raw_srccat->catsources[0].shape_param_indexes = shape0_param_indexes;

    raw_srccat->catsources[0].n_comps += 3;
    raw_srccat->num_shapelets += 3;
    raw_srccat->catsources[0].n_shapes += 3;
    raw_srccat->catsources[0].n_shape_coeffs += 10;

    //Assign the second SOURCE some shapeIAN COMPONENTs
    raw_srccat->catsources[1].shape_ras = shape1_ras;
    raw_srccat->catsources[1].shape_decs = shape1_decs;
    raw_srccat->catsources[1].shape_ref_freqs = shape1_ref_freqs;
    raw_srccat->catsources[1].shape_ref_stokesI = shape1_ref_stokesI;
    raw_srccat->catsources[1].shape_ref_stokesQ = shape1_ref_stokesQ;
    raw_srccat->catsources[1].shape_ref_stokesU = shape1_ref_stokesU;
    raw_srccat->catsources[1].shape_ref_stokesV = shape1_ref_stokesV;
    raw_srccat->catsources[1].shape_SIs = shape1_SIs;
    raw_srccat->catsources[1].shape_majors = shape1_majors;
    raw_srccat->catsources[1].shape_minors = shape1_minors;
    raw_srccat->catsources[1].shape_pas = shape1_pas;
    raw_srccat->catsources[1].shape_coeffs = shape1_coeffs;
    raw_srccat->catsources[1].shape_n1s = shape1_n1s;
    raw_srccat->catsources[1].shape_n2s = shape1_n2s;
    raw_srccat->catsources[1].shape_param_indexes = shape1_param_indexes;

    raw_srccat->catsources[1].n_comps += 3;
    raw_srccat->num_shapelets += 3;
    raw_srccat->catsources[1].n_shapes += 3;
    raw_srccat->catsources[1].n_shape_coeffs += 7;

  }

  return raw_srccat;
}

/*
This checks that the overall number of retained COMPONENTs is correct. This
depends on whether I've injected POINT, GAUSSIAN, or SHAPELETS into the sky,
whether I'm cropping by SOURCE or COMPONENT, and what the LST is
*/

void check_num_components_after_crop(double lst_base, catsource_t *cropped_src,
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
      TEST_ASSERT_EQUAL_INT(3*num_component_types, cropped_src->n_comps);

    } else if (sky_crop_type == CROP_COMPONENTS) {
      TEST_ASSERT_EQUAL_INT(4*num_component_types, cropped_src->n_comps);
    }
  }
  else if (lst_base == LST1) { //For this LST, both SOURCEs should be retained
    TEST_ASSERT_EQUAL_INT(6*num_component_types, cropped_src->n_comps);
  }
  else if (lst_base == LST2) { //first source should have been cropped, second retained
    TEST_ASSERT_EQUAL_INT(3*num_component_types, cropped_src->n_comps);
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
void test_azza_coords(double lst_base, int num_comps,
                      user_precision_t *calc_azs, user_precision_t *calc_zas){

  //Compare calculated az/za for all time steps to expectation
  //Some of these are small numbers so test within an accuracy
  for (size_t azza_ind = 0; azza_ind < num_comps*3; azza_ind++) {

    if (lst_base == LST0) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_az_LST0[azza_ind],
                               calc_azs[azza_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_za_LST0[azza_ind],
                               calc_zas[azza_ind]);
    } else if (lst_base == LST1) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_az_LST1[azza_ind],
                               calc_azs[azza_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_za_LST1[azza_ind],
                               calc_zas[azza_ind]);
    } else if (lst_base == LST2) {
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_az_LST2[azza_ind],
                               calc_azs[azza_ind]);
      TEST_ASSERT_DOUBLE_WITHIN(TOL, expec_za_LST2[azza_ind],
                               calc_zas[azza_ind]);
    }
  }
}

/*******************************************************************************
---------------POINT SOURCE functions------------------------------------------
*******************************************************************************/

/*
Checks the POINT source values match expected for given LST, when cropping
by SOURCE or by COMPONENT
*/
void check_point_values_after_crop(double lst_base, catsource_t *cropped_src,
                                    e_sky_crop sky_crop_type) {

  //Test if the az/za coords have been calculated correctly
  test_azza_coords(lst_base, cropped_src->n_points,
                          cropped_src->point_azs, cropped_src->point_zas);

  //For LST0, all COMPONENTs of the first POINT SOURCE
  //(as defined in test_crop_sky_model.h) should have been retained.
  //Only the first COMPONENT of the second SOURCE is above the horizon, so
  //if cropping by SOURCE there should only be three COMPONENTs retained.
  //If cropping by COMPONENT, there should be a 4th component
  if (lst_base == LST0) {

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point0_ras, cropped_src->point_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point0_decs, cropped_src->point_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point0_ref_freqs, cropped_src->point_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesI,
                                  cropped_src->point_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesQ,
                                  cropped_src->point_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesU,
                                  cropped_src->point_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesV,
                                  cropped_src->point_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_SIs, cropped_src->point_SIs, 3);

    if (sky_crop_type == CROP_SOURCES) {
      TEST_ASSERT_EQUAL_INT(3, cropped_src->n_points);

    } else if (sky_crop_type == CROP_COMPONENTS) {
      TEST_ASSERT_EQUAL_INT(4, cropped_src->n_points);
      //test that 4th COMPONENT is equal to the first of point1
      TEST_ASSERT_EQUAL_DOUBLE(point1_ras[0], cropped_src->point_ras[3]);
      TEST_ASSERT_EQUAL_DOUBLE(point1_decs[0], cropped_src->point_decs[3]);
      TEST_ASSERT_EQUAL_DOUBLE(point1_ref_freqs[0], cropped_src->point_ref_freqs[3]);
      TEST_ASSERT_EQUAL_FLOAT(point1_ref_stokesI[0],
                             cropped_src->point_ref_stokesI[3]);
      TEST_ASSERT_EQUAL_FLOAT(point1_ref_stokesQ[0],
                             cropped_src->point_ref_stokesQ[3]);
      TEST_ASSERT_EQUAL_FLOAT(point1_ref_stokesU[0],
                             cropped_src->point_ref_stokesU[3]);
      TEST_ASSERT_EQUAL_FLOAT(point1_ref_stokesV[0],
                             cropped_src->point_ref_stokesV[3]);
      TEST_ASSERT_EQUAL_FLOAT(point1_SIs[0], cropped_src->point_SIs[3]);

    }
  }
  else if (lst_base == LST1) { //For this LST, both SOURCEs should be retained
    TEST_ASSERT_EQUAL_INT(6, cropped_src->n_points);

    //First three COMPONENTs should match point0
    // TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ras, cropped_src->point_ras, 3);
    // TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_decs, cropped_src->point_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point0_ras, cropped_src->point_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point0_decs, cropped_src->point_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point0_ref_freqs, cropped_src->point_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesI,
                                  cropped_src->point_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesQ,
                                  cropped_src->point_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesU,
                                  cropped_src->point_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_ref_stokesV,
                                  cropped_src->point_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point0_SIs, cropped_src->point_SIs, 3);

    //Last three COMPONENTs should match point1
    //Gooooooooooo pointer arithmatic
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point1_ras, cropped_src->point_ras + 3, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point1_decs, cropped_src->point_decs + 3, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point1_ref_freqs, cropped_src->point_ref_freqs + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesI,
                                  cropped_src->point_ref_stokesI + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesQ,
                                  cropped_src->point_ref_stokesQ + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesU,
                                  cropped_src->point_ref_stokesU + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesV,
                                  cropped_src->point_ref_stokesV + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_SIs, cropped_src->point_SIs + 3, 3);

  } else if (lst_base == LST2) { //first source should have been cropped, second retained
    TEST_ASSERT_EQUAL_INT(3, cropped_src->n_points);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point1_ras, cropped_src->point_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point1_decs, cropped_src->point_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(point1_ref_freqs, cropped_src->point_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesI,
                                  cropped_src->point_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesQ,
                                  cropped_src->point_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesU,
                                  cropped_src->point_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_ref_stokesV,
                                  cropped_src->point_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(point1_SIs, cropped_src->point_SIs, 3);

  } else if (lst_base == LST3) { //both should be flagged
    TEST_ASSERT_EQUAL_INT(0, cropped_src->n_points);
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
  catsource_t *cropped_src;
  cropped_src =  crop_sky_model(raw_srccat, lsts, (double)MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  check_point_values_after_crop(lst_base, cropped_src, sky_crop_type);

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


/*******************************************************************************
---------------GAUSSIAN SOURCE functions----------------------------------------
*******************************************************************************/

/*
Checks the GAUSSIAN source values match expected for given LST, when cropping
by SOURCE or by COMPONENT
*/
void check_gauss_values_after_crop(double lst_base, catsource_t *cropped_src,
                                    e_sky_crop sky_crop_type) {
  //Test if the az/za coords have been calculated correctly
  test_azza_coords(lst_base, cropped_src->n_gauss,
                          cropped_src->gauss_azs, cropped_src->gauss_zas);


  //For LST0, all COMPONENTs of the first GAUSSIAN SOURCE
  //(as defined in test_crop_sky_model.h) should have been retained.
  //Only the first COMPONENT of the second SOURCE is above the horizon, so
  //if cropping by SOURCE there should only be three COMPONENTs retained.
  //If cropping by COMPONENT, there should be a 4th component
  if (lst_base == LST0) {
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss0_ras, cropped_src->gauss_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss0_decs, cropped_src->gauss_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss0_ref_freqs, cropped_src->gauss_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesI,
                                  cropped_src->gauss_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesQ,
                                  cropped_src->gauss_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesU,
                                  cropped_src->gauss_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesV,
                                  cropped_src->gauss_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_SIs, cropped_src->gauss_SIs, 3);

    if (sky_crop_type == CROP_SOURCES) {
      TEST_ASSERT_EQUAL_INT(3, cropped_src->n_gauss);

    } else if (sky_crop_type == CROP_COMPONENTS) {
      TEST_ASSERT_EQUAL_INT(4, cropped_src->n_gauss);
      //test that 4th COMPONENT is equal to the first of gauss1
      TEST_ASSERT_EQUAL_DOUBLE(gauss1_ras[0], cropped_src->gauss_ras[3]);
      TEST_ASSERT_EQUAL_DOUBLE(gauss1_decs[0], cropped_src->gauss_decs[3]);
      TEST_ASSERT_EQUAL_DOUBLE(gauss1_ref_freqs[0], cropped_src->gauss_ref_freqs[3]);
      TEST_ASSERT_EQUAL_FLOAT(gauss1_ref_stokesI[0],
                              cropped_src->gauss_ref_stokesI[3]);
      TEST_ASSERT_EQUAL_FLOAT(gauss1_ref_stokesQ[0],
                              cropped_src->gauss_ref_stokesQ[3]);
      TEST_ASSERT_EQUAL_FLOAT(gauss1_ref_stokesU[0],
                              cropped_src->gauss_ref_stokesU[3]);
      TEST_ASSERT_EQUAL_FLOAT(gauss1_ref_stokesV[0],
                              cropped_src->gauss_ref_stokesV[3]);
      TEST_ASSERT_EQUAL_FLOAT(gauss1_SIs[0], cropped_src->gauss_SIs[3]);
    }
  }
  else if (lst_base == LST1) { //For this LST, both SOURCEs should be retained
    TEST_ASSERT_EQUAL_INT(6, cropped_src->n_gauss);

    //First three COMPONENTs should match gauss0
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss0_ras, cropped_src->gauss_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss0_decs, cropped_src->gauss_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss0_ref_freqs, cropped_src->gauss_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesI,
                                  cropped_src->gauss_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesQ,
                                  cropped_src->gauss_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesU,
                                  cropped_src->gauss_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_ref_stokesV,
                                  cropped_src->gauss_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss0_SIs, cropped_src->gauss_SIs, 3);

    //Last three COMPONENTs should match gauss1
    //Gooooooooooo gausser arithmatic
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss1_ras, cropped_src->gauss_ras + 3, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss1_decs, cropped_src->gauss_decs + 3, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss1_ref_freqs,
                                  cropped_src->gauss_ref_freqs + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesI,
                                  cropped_src->gauss_ref_stokesI + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesQ,
                                  cropped_src->gauss_ref_stokesQ + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesU,
                                  cropped_src->gauss_ref_stokesU + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesV,
                                  cropped_src->gauss_ref_stokesV + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_SIs, cropped_src->gauss_SIs + 3, 3);

  } else if (lst_base == LST2) { //first source should have been cropped, second retained
    TEST_ASSERT_EQUAL_INT(3, cropped_src->n_gauss);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss1_ras, cropped_src->gauss_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss1_decs, cropped_src->gauss_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(gauss1_ref_freqs, cropped_src->gauss_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesI,
                                  cropped_src->gauss_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesQ,
                                  cropped_src->gauss_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesU,
                                  cropped_src->gauss_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_ref_stokesV,
                                  cropped_src->gauss_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(gauss1_SIs, cropped_src->gauss_SIs, 3);

  } else if (lst_base == LST3) { //both should be flagged
    TEST_ASSERT_EQUAL_INT(0, cropped_src->n_gauss);
  }
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
  catsource_t *cropped_src;
  cropped_src =  crop_sky_model(raw_srccat, lsts, MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  check_gauss_values_after_crop(lst_base, cropped_src, sky_crop_type);

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

/*
Checks the SHAPELET source values match expected for given LST, when cropping
by SOURCE or by COMPONENT
*/
void check_shape_values_after_crop(double lst_base, catsource_t *cropped_src,
                                    e_sky_crop sky_crop_type) {

  //Test if the az/za coords have been calculated correctly
  test_azza_coords(lst_base, cropped_src->n_shapes,
                          cropped_src->shape_azs, cropped_src->shape_zas);

  //For LST0, all COMPONENTs of the first SHAPELET SOURCE
  //(as defined in test_crop_sky_model.h) should have been retained.
  //Only the first COMPONENT of the second SOURCE is above the horizon, so
  //if cropping by SOURCE there should only be three COMPONENTs retained.
  //If cropping by COMPONENT, there should be a 4th component
  if (lst_base == LST0) {
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape0_ras, cropped_src->shape_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape0_decs, cropped_src->shape_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape0_ref_freqs, cropped_src->shape_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesI,
                                  cropped_src->shape_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesQ,
                                  cropped_src->shape_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesU,
                                  cropped_src->shape_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesV,
                                  cropped_src->shape_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_SIs, cropped_src->shape_SIs, 3);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_coeffs, cropped_src->shape_coeffs, 10);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n1s, cropped_src->shape_n1s, 10);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n2s, cropped_src->shape_n2s, 10);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_param_indexes, cropped_src->shape_param_indexes, 10);

    if (sky_crop_type == CROP_SOURCES) {
      TEST_ASSERT_EQUAL_INT(3, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_INT(10, cropped_src->n_shape_coeffs);

    } else if (sky_crop_type == CROP_COMPONENTS) {
      TEST_ASSERT_EQUAL_INT(4, cropped_src->n_shapes);
      TEST_ASSERT_EQUAL_INT(14, cropped_src->n_shape_coeffs);
      //test that 4th COMPONENT is equal to the first of shape1
      TEST_ASSERT_EQUAL_DOUBLE(shape1_ras[0], cropped_src->shape_ras[3]);
      TEST_ASSERT_EQUAL_DOUBLE(shape1_decs[0], cropped_src->shape_decs[3]);
      TEST_ASSERT_EQUAL_DOUBLE(shape1_ref_freqs[0], cropped_src->shape_ref_freqs[3]);
      TEST_ASSERT_EQUAL_FLOAT(shape1_ref_stokesI[0], cropped_src->shape_ref_stokesI[3]);
      TEST_ASSERT_EQUAL_FLOAT(shape1_ref_stokesQ[0], cropped_src->shape_ref_stokesQ[3]);
      TEST_ASSERT_EQUAL_FLOAT(shape1_ref_stokesU[0], cropped_src->shape_ref_stokesU[3]);
      TEST_ASSERT_EQUAL_FLOAT(shape1_ref_stokesV[0], cropped_src->shape_ref_stokesV[3]);
      TEST_ASSERT_EQUAL_FLOAT(shape1_SIs[0], cropped_src->shape_SIs[3]);

      //The last 4 coeffs, n1s, n2s should match the first 4 of shape1
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_coeffs, cropped_src->shape_coeffs + 10, 4);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n1s, cropped_src->shape_n1s + 10, 4);
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n2s, cropped_src->shape_n2s + 10, 4);
      //There is now a 4th component, so the param indexes should be 3
      user_precision_t expec_param_indexes1[] = {3, 3, 3, 3};
      TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_param_indexes1,
                                    cropped_src->shape_param_indexes + 10, 4);


    }
  }
  else if (lst_base == LST1) { //For this LST, both SOURCEs should be retained
    TEST_ASSERT_EQUAL_INT(6, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_INT(17, cropped_src->n_shape_coeffs);

    //First three COMPONENTs should match shape0
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape0_ras, cropped_src->shape_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape0_decs, cropped_src->shape_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape0_ref_freqs, cropped_src->shape_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesI,
                                  cropped_src->shape_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesQ,
                                  cropped_src->shape_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesU,
                                  cropped_src->shape_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_ref_stokesV,
                                  cropped_src->shape_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_SIs, cropped_src->shape_SIs, 3);

    //Last three COMPONENTs should match shape1
    //Gooooooooooo shapeer arithmatic
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape1_ras, cropped_src->shape_ras + 3, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape1_decs, cropped_src->shape_decs + 3, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape1_ref_freqs,
                                  cropped_src->shape_ref_freqs + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesI,
                                  cropped_src->shape_ref_stokesI + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesQ,
                                  cropped_src->shape_ref_stokesQ + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesU,
                                  cropped_src->shape_ref_stokesU + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesV,
                                  cropped_src->shape_ref_stokesV + 3, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_SIs, cropped_src->shape_SIs + 3, 3);

    //There were 10 n1, n2, coeffs from the first SOURCE
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_coeffs, cropped_src->shape_coeffs, 10);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n1s, cropped_src->shape_n1s, 10);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_n2s, cropped_src->shape_n2s, 10);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape0_param_indexes, cropped_src->shape_param_indexes, 10);

    //There are 7 n1, n2, coeffs from the second SOURCE
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_coeffs, cropped_src->shape_coeffs + 10, 7);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n1s, cropped_src->shape_n1s + 10, 7);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n2s, cropped_src->shape_n2s + 10, 7);

    //The indexes are additive coeffs from second SOURCE are higher than shape1
    user_precision_t expec_param_indexes2[] = {3, 3, 3, 3, 4, 5, 5};
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(expec_param_indexes2,
                                  cropped_src->shape_param_indexes + 10, 7);

    //There are

  } else if (lst_base == LST2) { //first source should have been cropped, second retained
    TEST_ASSERT_EQUAL_INT(3, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_INT(7, cropped_src->n_shape_coeffs);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape1_ras, cropped_src->shape_ras, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape1_decs, cropped_src->shape_decs, 3);
    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(shape1_ref_freqs, cropped_src->shape_ref_freqs, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesI,
                                  cropped_src->shape_ref_stokesI, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesQ,
                                  cropped_src->shape_ref_stokesQ, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesU,
                                  cropped_src->shape_ref_stokesU, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_ref_stokesV,
                                  cropped_src->shape_ref_stokesV, 3);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_SIs, cropped_src->shape_SIs, 3);

    //There were 7 n1, n2, coeffs from the second SOURCE
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_coeffs, cropped_src->shape_coeffs, 7);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n1s, cropped_src->shape_n1s, 7);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_n2s, cropped_src->shape_n2s, 7);
    TEST_ASSERT_EQUAL_FLOAT_ARRAY(shape1_param_indexes,
                                  cropped_src->shape_param_indexes, 7);

  } else if (lst_base == LST3) { //both should be flagged
    TEST_ASSERT_EQUAL_INT(0, cropped_src->n_shapes);
    TEST_ASSERT_EQUAL_INT(0, cropped_src->n_shape_coeffs);
  }
}

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
  catsource_t *cropped_src;
  cropped_src = crop_sky_model(raw_srccat, lsts, (double)MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  check_shape_values_after_crop(lst_base, cropped_src, sky_crop_type);

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

  //Just use SHAPELET type COMPONENTs
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
  catsource_t *cropped_src;
  cropped_src = crop_sky_model(raw_srccat, lsts, (double)MWA_LAT_RAD,
                                num_time_steps, sky_crop_type);

  check_point_values_after_crop(lst_base, cropped_src, sky_crop_type);
  check_gauss_values_after_crop(lst_base, cropped_src, sky_crop_type);
  check_shape_values_after_crop(lst_base, cropped_src, sky_crop_type);

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
