#include "woden_struct_defs.h"

/**
@brief Sets all values associated with a SHAPELET component in `temp_cropped_src`
to NULL

@details When splitting up a sky model into chunks, and we know there are no
SHAPELET components, set everything in `temp_cropped_src` that refers to a SHAPELET
source component to a NULL or zero to ensure no SHAPELET sources are simulated

@param[in] *temp_cropped_src Pointer to a `catsource_t` struct
*/
void null_point_comps(catsource_t *temp_cropped_src);

/**
@brief Sets all values associated with a SHAPELET component in `temp_cropped_src`
to NULL

@details When splitting up a sky model into chunks, and we know there are no
SHAPELET components, set everything in `temp_cropped_src` that refers to a SHAPELET
source component to a NULL or zero to ensure no SHAPELET sources are simulated

@param[in] *temp_cropped_src Pointer to a `catsource_t` struct
*/
void null_gauss_comps(catsource_t *temp_cropped_src);

/**
@brief Sets all values associated with a SHAPELET component in `temp_cropped_src`
to NULL

@details When splitting up a sky model into chunks, and we know there are no
SHAPELET components, set everything in `temp_cropped_src` that refers to a
SHAPELET source component to a NULL or zero to ensure no SHAPELET sources are
simulated

@param[in] *temp_cropped_src pointer to a `catsource_t` struct
*/
void null_shapelet_comps(catsource_t *temp_cropped_src);

/**
@brief When chunking the sky model in `cropped_src`, use `point_iter` to
stick the correct parts of POINT components `cropped_src` into `temp_cropped_src`

@details Given the larger sky model `cropped_src` and the index `point_iter`,
use pointer arithmatic to update the POINT component related attributes in
`temp_cropped_src` to point to the correct chunk of `cropped_src`. Also need
`num_time_steps` to correctly index the az/za coords as az/za change per time

@param[in,out] *temp_cropped_src Pointer to a `catsource_t` struct to contain
the 'chunked' sky model
@param[in] *cropped_src Pointer to a `catsource_t` struct that contains the full
sky model
@param[in] *point_iter Pointer to integer to correctly index the information
in `cropped_src` and point to in `temp_cropped_src`
@param[in] num_time_steps Number of time steps in the simulation
*/
void increment_point(catsource_t *temp_cropped_src, catsource_t *cropped_src,
                     int * point_iter, int num_time_steps);

/**
@brief When chunking the sky model in `cropped_src`, use `gauss_iter` to
stick the correct parts of GAUSS components `cropped_src` into `temp_cropped_src`

@details Given the larger sky model `cropped_src` and the index `gauss_iter`,
use pointer arithmatic to update the GAUSS component related attributes in
`temp_cropped_src` to point to the correct chunk of `cropped_src`. Also need
`num_time_steps` to correctly index the az/za coords as az/za change per time

@param[in,out] *temp_cropped_src Pointer to a `catsource_t` struct to contain
the 'chunked' sky model
@param[in] *cropped_src Pointer to a `catsource_t` struct that contains the full
sky model
@param[in] *gauss_iter Pointer to integer to correctly index the information
in `cropped_src` and point to in `temp_cropped_src`
@param[in] num_time_steps Number of time steps in the simulation
*/
void increment_gauss(catsource_t *temp_cropped_src, catsource_t *cropped_src,
                     int * gauss_iter, int num_time_steps);

/**
@brief When chunking the sky model in `cropped_src`, use `shape_iter` to
stick the correct parts of SHAPELET components `cropped_src` into
`temp_cropped_src`

@details Given the larger sky model `cropped_src` and the index `shape_iter`,
use pointer arithmatic to update the SHAPELET component related attributes in
`temp_cropped_src` to point to the correct chunk of `cropped_src`. Also need
`num_time_steps` to correctly index the az/za coords as az/za change per time

@param[in,out] *temp_cropped_src Pointer to a `catsource_t` struct to contain
the 'chunked' sky model
@param[in] *cropped_src Pointer to a `catsource_t` struct that contains the full
sky model
@param[in] *shape_iter Pointer to integer to correctly index the information
in `cropped_src` and point to in `temp_cropped_src`
@param[in] num_time_steps Number of time steps in the simulation
*/
void increment_shapelet(catsource_t *temp_cropped_src, catsource_t *cropped_src,
                     int * shape_iter, int num_time_steps);

// /**
// @brief When splitting the sky model `cropped_src` into `num_chunks` smaller
// sky models each of size `chunking_size`, fill the chunked sky model
// `temp_cropped_src`at index `chunk_ind` with the correct number of POINT,
// GAUSSIAN, and SHAPELET parameters.
//
// @details Here we are splitting the sky model `cropped_src` into `num_chunks`
// smaller sky models, each containing a number (`chunking_size`) of components.
// This function returns a single 'chunked' sky model in `temp_cropped_src`, which
// when iterating through `num_chunks`, corresponds to the chunked model at index
// `chunk_ind`.
//
// The function performs a large amount
// of logic to generate the correct combination of POINT, GAUSSIAN, and SHAPELET
// type components, that add up to the `chunking_size`, and not include any
// previously used components. Furthermore, as each SHAPELET component is expected
// to have multiple SHAPELET_coefficients, it is more efficient to split the sky
// model by the number of SHAPELET_coefficients rather than SHAPELET components.
// `point_iter`, `gauss_iter`, `shape_iter` are all used by the function
// `make_beam_settings_chunk` so returns their updated values.
//
// @param[in,out] *temp_cropped_src Pointer to a `catsource_t` struct to contain
// the 'chunked' sky model
// @param[in] *cropped_src Pointer to a `catsource_t` struct that contains the full
// sky model
// @param[in] num_chunks The number of chunks `cropped_src` is being split into
// @param[in] chunk_ind The index of the chunked sky model to be returned
// @param[in] chunking_size Number of components to put in each chunk
// @param[in] num_time_steps Number of time steps in the simulation
// @param[in,out] *point_iter Pointer to integer to correctly index the POINT
// component information in `cropped_src`
// @param[in,out] *gauss_iter Pointer to integer to correctly index the GAUSSIAN
// component information in `cropped_src`
// @param[in,out] *shape_iter Pointer to integer to correctly index the SHAPELET
// component information in `cropped_src`
// */
// void fill_chunk_src(catsource_t *temp_cropped_src, catsource_t *cropped_src,
//      int num_chunks, int chunk_ind, int chunking_size, int num_time_steps,
//      int * point_iter, int * gauss_iter, int * shape_iter);
//
// /**
// @brief Setup primary beam attributes required to match the 'chunked' sky model
//
// @details This is mostly needed when simulating using a Gaussian primary beam
// as currently the Gaussian primary beam is defined in hour angle/declination
// coordinates, which are stored in `beam_settings`. Need chunk `bmea_settings` to
// match the components in the given chunked sky model `temp_cropped_src`. SHAPELET
// information in the full sky model `cropped_src` is also needed.
//
// @param[in] beam_settings A `beam_settings_t` containing primary beam
// settings for the full sky model `cropped_src`
// @param[in] *temp_cropped_src Pointer to a `catsource_t` struct that contains the
// chunked sky model
// @param[in] *cropped_src Pointer to a `catsource_t` struct that contains the full
// sky model
// @param[in] *woden_settings Pointer to woden_settings_t containing simulation
// settings
// @param[in] point_iter Integer to correctly index the POINT
// component information in `beam_settings` and `cropped_src`
// @param[in] gauss_iter Integer to correctly index the GAUSSIAN
// component information in `beam_settings` and `cropped_src`
// @param[in] shape_iter Integer to correctly index the SHAPELET
// component information in `beam_settings` and `cropped_src`
//
// @return `beam_settings_chunk` A `beam_settings_t` struct containing the
// 'chunked' beam parameters that correspond to the components in `temp_cropped_src`
// */
// beam_settings_t make_beam_settings_chunk(beam_settings_t beam_settings,
//                 catsource_t *temp_cropped_src, catsource_t *cropped_src,
//                 woden_settings_t *woden_settings,
//                 int point_iter, int gauss_iter, int shape_iter);

void fill_chunk_src_with_pointgauss(catsource_t *temp_cropped_src,
     catsource_t *cropped_src, int chunk_ind, int comps_per_chunk,
     int num_time_steps);

source_catalogue_t * create_chunked_sky_models(catsource_t *cropped_src,
                                              woden_settings_t *woden_settings);
