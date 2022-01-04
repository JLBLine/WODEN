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

@details Mallocs `num_visis*sizeof(user_precision_t)` for the following arrays:

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
                                  double base_band_freq,
                                  double *lsts);

/**
@brief Write out the simulated visibilities to a binary file

@details Dumps u,v,w (metres), XX Real, XX Imag, XY Real, XY Imag,
YX Real, YX Imag, YY Real, YY Imag (Jy) to a binary file
Output order is by baseline (fastest changing), frequency, time (slowest
changing). This means the us_metres, vs_metres, ws_metres are repeated over
frequency, but keeps the dimensions of the output sane. Output is written as
"output_visi_band%02d.dat", `band_num`.

@param[in] *visibility_set `visibility_set_t` populated with visibility
data to be written out
@param[in] band_num Number of the band being processed (used for output name
@param[in] num_visis Number of visibilities to output

*/
void write_visi_set_binary(visibility_set_t *visibility_set,
                           int band_num, int num_visis);
/**
@brief Write out XX visibilities and u,v,w to a text file

@details Used for desperation debugging. Just writes u,v,w (metres),
XX Real, XX Imag (Jy)

@param[in] *visibility_set `visibility_set_t` populated with visibility
data to be written out
@param[in] band_num Number of the band being processed (used for output name
@param[in] *woden_settings Populated `woden_settings_t` used to run the
simulation
*/
void write_visi_set_text(visibility_set_t *visibility_set, int band_num,
                         woden_settings_t *woden_settings);

/**
@brief Free the settings that were used for this visibility set

@details We don't free everything inside visibility_set as inside
`calculate_visibilities`we set these arrays as pointers as don't want to
free the original allocation. This function frees:

    visibility_set->allsteps_sha0s
    visibility_set->allsteps_cha0s
    visibility_set->allsteps_lsts
    visibility_set->allsteps_wavelengths
    visibility_set->channel_frequencies


@param[in,out] *visibility_set Contains arrays that need freeing
*/
void free_visi_set_inputs(visibility_set_t *visibility_set);

/**
@brief Free the output arrays in `visibility_set`

@details This function frees:

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

@param[in,out] *visibility_set Contains arrays that need freeing
*/
void free_visi_set_outputs(visibility_set_t *visibility_set);



#ifdef __cplusplus
}
#endif
