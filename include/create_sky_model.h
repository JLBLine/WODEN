#include "read_and_write.h"

#define MAX_NUM_COMPONENTS 10000

typedef enum {CROP_SOURCES, CROP_COMPONENTS}e_sky_crop;

void convert_radec2azza(double ra, double dec, double lst,
     double * az, double * za);


catsource_t * crop_sky_model(source_catalogue_t *raw_srccat, float *lsts,
             int num_time_steps, e_sky_crop sky_crop_type);
