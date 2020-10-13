// #include "woden.h"
#include "read_and_write.h"

typedef enum {CROP_SOURCES, CROP_COMPONENTS}e_sky_crop;
typedef enum {BELOW, ABOVE}e_horizon;

void convert_radec2azza(double ra, double dec, double lst,
     double * az, double * za);

void horizon_test(double za, e_sky_crop sky_crop_type,
     e_horizon * all_comps_above_horizon, int * num_comp_retained,
     int * num_shape_coeff_retained, int num_shape_coeff_component,
     float *shape_param_indexes, int shape);

catsource_t * crop_sky_model(source_catalogue_t *raw_srccat, float *lsts,
             int num_time_steps, e_sky_crop sky_crop_type);
