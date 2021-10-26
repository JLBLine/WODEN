#include "expected_para_angles.h"

float lsts[] = {1*DD2R, 120*DD2R, 240*DD2R};

double point_ras[] = {0 , 5*DD2R, 10*DD2R};
double point_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};

double gauss_ras[] = {115*DD2R , 125*DD2R, 130*DD2R};
double gauss_decs[] = {-15*DD2R , -20*DD2R, -25*DD2R};

double shape_ras[] = {235*DD2R , 250*DD2R, 265*DD2R};
double shape_decs[] = {-35*DD2R , -40*DD2R, -45*DD2R};

float expec_point_sin_para[] = {-0.0039216, -0.6142135, -0.6142136, -0.0156904,
                                -0.5764040, -0.6534451, -0.0353442, -0.5401087,
                                -0.6939089 };

float expec_point_cos_para[] = {-0.9999923, -0.7891399, 0.7891398, 0.9998769,
                                -0.8171648, 0.7569739, 0.9993752, -0.8415953,
                                 0.7200627 };

float expec_gauss_sin_para[] = {-0.5432732, -0.9341043, -0.6122640, -0.6251023,
                                -0.8344397, -0.5628954, -0.6805615, -0.2228374,
                                -0.5384376 };

float expec_gauss_cos_para[] = {0.8395559, -0.3570003, -0.7906534, 0.7805428,
                                0.5510992, -0.8265281, 0.7326910, 0.9748556,
                                -0.8426654 };

float expec_shape_sin_para[] = {-0.6794566, -0.5854309, 0.8773373, -0.5505476,
                                -0.7216252, 0.8188159, -0.3965520, -0.8523725,
                                 0.5593383 };

float expec_shape_cos_para[] = {-0.7337157, 0.8107223, -0.4798742, -0.8348038,
                                 0.6922840, 0.5740561, -0.9180123, 0.5229351,
                                 0.8289395 };

/*
Make the polpulated catsource_t struct. Stick in some necessary values
*/
catsource_t * make_sky_model(void) {

  catsource_t *src = malloc(sizeof(catsource_t));

  src->n_comps = 9;
  src->n_points = 3;
  src->n_gauss = 3;
  src->n_shapes = 3;

  src->point_ras = point_ras;
  src->point_decs = point_decs;

  src->gauss_ras = gauss_ras;
  src->gauss_decs = gauss_decs;

  src->shape_ras = shape_ras;
  src->shape_decs = shape_decs;

  return src;
}
