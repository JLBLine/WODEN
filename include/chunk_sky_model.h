#include "woden_struct_defs.h"
#include "create_sky_model.h"

//Going to be calling this code in both C and CUDA, so stick the
//conditional extern C around so linkage by compilers survives all
//the mangling fun

#ifdef __cplusplus
extern "C" {
#endif

/**
@brief Sets all values associated with either the POINT, GAUSSIAN, or SHAPELET
components in `source_t *src` to NULL, depending on which type is given in
`component_type`.

@details When splitting up a sky model into chunks, often pass around a `source_t`
struct to do the splitting up, so good to have a function that ensures certain
fields are set to NULL and avoid copying incorrect parts of a model

@param[in] *src Pointer to a `source_t` struct
@param[in] component_type Can be either POINT, GAUSSIAN, or SHAPELET
*/
void null_components(source_t *src, e_component_type component_type);

/**
@brief When chunking the sky model in `cropped_src`, use `iter` to
stick the correct parts of either POINT or GAUSSIAN components from `cropped_src`
into `temp_cropped_src`

@details Given the larger sky model `cropped_src` and the index `iter`,
use pointer arithmatic to update the POINT/GAUSSIAN component related attributes
in `temp_cropped_src` to point to the correct chunk of `cropped_src`. Choose
which component type using `component_type`. Also need
`num_time_steps` to correctly index the az/za coords as az/za change per time.

@param[in,out] *temp_cropped_src Pointer to a `source_t` struct to contain
the 'chunked' sky model
@param[in] *cropped_src Pointer to a `source_t` struct that contains the full
sky model
@param[in] component_type Either POINT or GAUSSIAN
@param[in] power_iter Used to correctly index the information
in `cropped_src` for POWER_LAW type components
@param[in] curve_iter Used to correctly index the information
in `cropped_src` for CURVED_POWER_LAW type components
@param[in] list_iter Used to correctly index the information
in `cropped_src` for LIST type components
@param[in] num_chunk_power How many POWER_LAW type components are in this chunk
@param[in] num_chunk_curve How many CURVE_POWER_LAW type components are in this chunk
@param[in] num_chunk_list How many LIST type components are in this chunk
*/
void increment_pointgauss(source_t *temp_cropped_src, source_t *cropped_src,
                          e_component_type component_type,
                          int power_iter, int curve_iter, int list_iter,
                          int num_chunk_power, int num_chunk_curve, int num_chunk_list);

/**
@brief When splitting the sky model `cropped_src` into `num_chunks` smaller
sky models each of size `chunking_size`, fill the chunked sky model
`temp_cropped_src` at index `chunk_ind` with the correct number of POINT and
GAUSSIAN component settings.

@details Here we are splitting the sky model `cropped_src` into `num_chunks`
smaller sky models, each containing a number (`comps_per_chunk`) of components.
This function returns a single 'chunked' sky model in `temp_cropped_src`, which
when iterating through the total number of POINT and GAUSSIAN components
in chunks of size `comps_per_chunk`, corresponds to the chunked model at index
`chunk_ind`.


@param[in,out] *temp_cropped_src Pointer to a `source_t` struct to contain
the 'chunked' sky model
@param[in] *cropped_src Pointer to a `source_t` struct that contains the full
sky model
@param[in] chunk_ind The index of the chunked sky model to be returned
@param[in] comps_per_chunk Number of components to put in each chunk
@param[in] *woden_settings A populated `woden_settings_t` containing the
simulation setting
@param[in] comptype Either a POINT or GAUSSIAN component type
*/
void fill_chunk_src_with_pointgauss(source_t *temp_cropped_src,
     source_t *cropped_src, int chunk_ind, int comps_per_chunk,
     woden_settings_t *woden_settings, e_component_type comptype);



/**
@brief When chunking the sky model in `cropped_src`, use `shape_iter` to
stick the correct parts of SHAPELET components `cropped_src` into
`temp_cropped_src`

@details Given the larger sky model `cropped_src` and the index `shape_iter`,
use pointer arithmatic to update the SHAPELET component related attributes in
`temp_cropped_src` to point to the correct chunk of `cropped_src`. Also need
`num_time_steps` to correctly index the az/za coords as az/za change per time

@param[in,out] *temp_cropped_src Pointer to a `source_t` struct to contain
the 'chunked' sky model
@param[in] *cropped_src Pointer to a `source_t` struct that contains the full
sky model
@param[in] *shape_iter Pointer to integer to correctly index the information
in `cropped_src` and point to in `temp_cropped_src`
*/
void increment_shapelet(source_t *temp_cropped_src, source_t *cropped_src,
                     int * shape_iter);



/**
@brief When splitting the sky model `cropped_src` into `num_chunks` smaller
sky models each of size `chunking_size`, fill the chunked sky model
`temp_cropped_src` at index `chunk_ind` with the correct number of SHAPELET
component settings.

@details When splitting the sky model up, we split SHAPELET components
separately to POINT and GAUSSIANs as the besis functions take up more GPU
memory - this can butt heads with primary beam calculations that also take up
memory for POINT and GAUSSIANs. This function returns a single 'chunked' sky
model in `temp_cropped_src`, which when iterating through the total number of
SHAPELET basis functions in chunks of size `comps_per_chunk`, corresponds to
the chunked model at index `chunk_ind`.

@param[in,out] *temp_cropped_src Pointer to a `source_t` struct to contain
the 'chunked' sky model
@param[in] *cropped_src Pointer to a `source_t` struct that contains the full
sky model
@param[in] chunk_ind The index of the chunked sky model to be returned
@param[in] coeffs_per_chunk Number of basis function params to put in each chunk
@param[in] *woden_settings A populated `woden_settings_t` containing the
simulation setting
*/
void fill_chunk_src_with_shapelets(source_t *temp_cropped_src,
     source_t *cropped_src, int chunk_ind, int coeffs_per_chunk,
     woden_settings_t *woden_settings);

/**
@brief Takes the sky model in `cropped_src` and splits it into bitesize pieces
to fit on the GPU. Returns the models inside a `source_catalogue_t` struct,
with the models stored in `source_catalogue_t->sources`.

@details This is basically a wrapper around `fill_chunk_src_with_pointgauss`
and `fill_chunk_src_with_shapelets`

@param[in] *cropped_src Pointer to a `source_t` struct that contains the full
@param[in] *woden_settings A populated `woden_settings_t` containing the
simulation settings, including the chunking size
@returns `source_catalogue_t` A struct containing the chunked sky models
that can be fed into `calculate_visibilities::calculate_visibilities`.
*/
source_catalogue_t * create_chunked_sky_models(source_t *cropped_src,
                                               woden_settings_t *woden_settings);



/**
@brief Remap a number of arrays from pointers to the full sky model to only what
needs to be copied over to the GPU, so we can fit things into GPU memory

@details Should have already run `chunk_sky_model::create_chunked_sky_models`,
which splits a full sky model into smaller chunks that can fit onto the GPU.
Some of this mapping is just copying a pointer, but we want to get things
ready to be copied across to the GPU, and depending on whether we have
POWER_LAW, CURVED_POWER_LAW, or LIST flux types, we need to copy across
any location in arrays such as those that hold RA,Dec, so we can't do a simple
pointer arithamatic. This function should be used while iterating over
the chunks in the sky model to assign enough memory to copy only what we need,
which can be passed over and copied to the GPU, and then freed straight after.

@param[in,out] *remapped_source Pointer to a new `source_t` struct that will hold the
remapped source ready to be given to the GPU
@param[in] *chunked_source Pointer to a populated `source_t` that was filled using
`chunk_sky_model::create_chunked_sky_models`
@param[in] num_time_steps Number of time steps in the simulation, effects the
number of az/za coords that need to be copied
@param[in] beamtype The type of primary beam used in the simulation. Only
MWA_ANALY and GAUSS_BEAM need to have the `beam_has` and `beam_decs` arrays
copied across
*/
void remap_source_for_gpu(source_t *remapped_source, source_t *chunked_source,
                          int num_time_steps, e_beamtype beamtype);




/**
@brief Free the memory in `source` which would have been assigned by `chunk_sky_model::remap_source_for_gpu`

@details

@param[in] *source Pointer to `source_t` struct that was populated using
`chunk_sky_model::remap_source_for_gpu`
@param[in] beamtype What type of primary beam was used in the simulation, effects what memory was allocated in the first place
*/
void free_remapped_source_for_gpu(source_t *source, e_beamtype beamtype);


#ifdef __cplusplus
}
#endif
