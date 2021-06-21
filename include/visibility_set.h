/*! \file
  Helper functions to do malloc-ing, filling, and free-ing of the
  `visibility_set_t` struct, which holds the outputs visiblities.
  @author J.L.B. Line
*/
#include "constants.h"
#include "woden_struct_defs.h"

//Going to be calling this code in both C and CUDA, so stick the
//conditional extern C around so linkage by compilers survives all
//the mangling fun

#ifdef __cplusplus
extern "C" {
#endif

/**
@brief Create a `visibility_set_t` and return it with a number of arrays
malloc-ed to the sizeof `num_visis`.

@details Mallocs `num_visis*sizeof(float)` for the following arrays:

    visibility_set->us_metres
    visibility_set->vs_metres
    visibility_set->ws_metres
    visibility_set->sum_visi_XX_real
    visibility_set->sum_visi_XX_imag
    visibility_set->sum_visi_XY_real
    visibility_set->sum_visi_XY_imag
    visibility_set->sum_visi_YX_real
    visibility_set->sum_visi_YX_imag
    visibility_set->sum_visi_YY_real
    visibility_set->sum_visi_YY_imag

@returns `visibility_set` with a number of malloc-ed arrays
*/
visibility_set_t * setup_visibility_set(int num_visis);

/**
@brief Fills in hour angle, LSTS, and wavelength/frequency information in
`visibility_set`.

@details Mallocs and fills in the following arrays:

    visibility_set->allsteps_sha0s
    visibility_set->allsteps_cha0s
    visibility_set->allsteps_lsts
    visibility_set->allsteps_wavelengths
    visibility_set->channel_frequencies

Anything with `allsteps` will have a value for every baseline (fastest changing),
every frequency, every time step (slowest changing). This makes it easy to
injest into CUDA kernels as a single axis. `channel_frequencies` will list
the frequencies in this 'band'.

@param[in,out] *visibility_set A `visibility_set_t` to populate with arrays
@param[in] *woden_settings Pointer to populated `woden_settings_t` struct
containing simulation settings
@param[in] base_band_freq The lowest frequency to simulate in this band
@param[in] *lsts Array of LST for all time steps in the simulation
*/
void fill_timefreq_visibility_set(visibility_set_t *visibility_set,
                                  woden_settings_t *woden_settings,
                                  float base_band_freq,
                                  float *lsts);

void write_visi_set_binary(visibility_set_t *visibility_set,
                           int band_num, int num_visis);

void write_visi_set_text(visibility_set_t *visibility_set, int band_num,
                         woden_settings_t *woden_settings);

void free_visi_set_inputs(visibility_set_t *visibility_set);

void free_visi_set_outputs(visibility_set_t *visibility_set);



#ifdef __cplusplus
}
#endif
