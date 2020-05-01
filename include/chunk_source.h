void null_point_comps(catsource_t *temp_cropped_src);
void null_gauss_comps(catsource_t *temp_cropped_src);
void null_shapelet_comps(catsource_t *temp_cropped_src);

void increment_point(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int chunk, int chunking_size, int num_time_steps);

void increment_gauss(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int chunk, int chunking_size, int gauss_iter, int num_time_steps);

void increment_shapelet(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int chunk, int chunking_size, int shape_iter, int num_time_steps)

void fill_chunk_src(catsource_t *temp_cropped_src, catsource_t *cropped_src,
     int num_chunks, int chunk, int chunking_size, int num_time_steps );
