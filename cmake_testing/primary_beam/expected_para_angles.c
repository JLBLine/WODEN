#include "expected_para_angles.h"
#include "woden_precision_defs.h"

user_precision_t lsts[] = {1*DD2R, 120*DD2R, 240*DD2R};

double point_ras[] = {0 , 5*DD2R, 10*DD2R};
double point_decs[] = {MWA_LAT_RAD, MWA_LAT_RAD, MWA_LAT_RAD};

double gauss_ras[] = {115*DD2R , 125*DD2R, 130*DD2R};
double gauss_decs[] = {-15*DD2R , -20*DD2R, -25*DD2R};

double shape_ras[] = {235*DD2R , 250*DD2R, 265*DD2R};
double shape_decs[] = {-35*DD2R , -40*DD2R, -45*DD2R};

user_precision_t expec_point_sin_para[] = {-0.003921568957513001, -0.6142135609790145,
 -0.6142135609790146, -0.015690440557775406, -0.5764041507979725, -0.6534450388756359,
 -0.03534414809761719, -0.5401088129215882, -0.6939090402675685 };

user_precision_t expec_point_cos_para[] = {-0.9999923106188925, -0.7891398491455608,
    0.7891398491455607, 0.9998768974604338, -0.8171647660924131, 0.7569739633361366,
    0.9993752004103632, -0.8415951937864384, 0.7200626666026659};

user_precision_t expec_gauss_sin_para[] = {-0.5432731274525338, -0.93410419485169,
  -0.6122639942221019, -0.6251022987747572, -0.8344396610646464, -0.5628954663710054,
  -0.6805613903651295, -0.22283738651783844, -0.5384376075227046};

user_precision_t expec_gauss_cos_para[] = {0.8395560189695164, -0.3570004946221726,
  -0.7906534015478578, 0.7805428342291756, 0.5510993123224871, -0.8265280962792301,
  0.7326910630984124, 0.9748556299113729, -0.842665380091781};

user_precision_t expec_shape_sin_para[] = {-0.6794568466308039, -0.5854307441903317,
   0.877337280555733, -0.5505476889106884, -0.7216251259281579, 0.8188159509192331,
   -0.3965520384424639, -0.8523725523916869, 0.5593381941309807};

user_precision_t expec_shape_cos_para[] = {-0.7337154718053343, 0.8107224209041923,
  -0.47987425034801673, -0.8348037147947414, 0.6922840295927462, 0.5740561283709391,
  -0.9180122443666678, 0.5229350169278025, 0.8289395542416186};

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
