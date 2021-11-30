/*! \file
  Host methods to read in hdf5 stored MWA FEE model, initialise the model
  for the given telescope delay settings, and accumulate the necessary arrays
  to make calculations on the device.
  @author J.L.B. Line
*/
#pragma once
#include "hdf5.h"
#include "hdf5_hl.h"
#include "constants.h"
#include "woden_precision_defs.h"
#include "woden_struct_defs.h"

/**
@brief Operator function to be called by H5Literate and `RTS_op_func`.

@details Used when finding the closest available frequency in the MWA FEE
primary beam model
*/
struct opdata {
  float freq_in; /**< Desired frequency (Hz) */
  float freq_out; /**< Closest available frequency (Hz) */
  float least_diff; /**< Difference between frequencies (Hz) */
};

/**
@brief Operator function to be called by `H5Literate`, finds closest frequency
in the element pattern hdf5 file

@details The MWA FEE beam is stored at a 1.28MHz frequency resolution. Using
this function in conjunction with `H5Literate`, find the frequencies stored
in the MWA FEE hdf5, and return the closest frequency to the frequency in
`operator_data.freq_in` (Hz)

*/
herr_t RTS_op_func(hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data);

/**
@brief Read in the spherical harmonic basis functions and stored coefficients
from `h5filename`, and calculate the tile phases for the given `FEE_delays` on
the host device to intialise the MWA FEE beam at the requested frequency.

@details The MWA primary beam is steered usings quantised delays on each of the
16 dipoles in the tile. Specify these delays using `FEE_delays` - they can be
found in the metafits file associated with an MWA observation

The MWA FEE beam is stored at a 1.28MHz frequency resolution, so the beam is
intialised at the closest stored frequency to `freq_Hz`.

@todo How to do some kind of interpolation over frequency in the future

@param[in] h5filename Path to `mwa_full_embedded_element_pattern.h5`
@param[in] freq_Hz Requested frequency (Hz)
@param[in,out] pb A `RTS_MWA_FEE_beam_t` struct to store outputs in
@param[in] FEE_delays Dipole delay factors to specify beam pointings

*/
int RTS_MWAFEEInit(const char *h5filename, double freq_Hz,
                   RTS_MWA_FEE_beam_t *pb, user_precision_t *FEE_delays);

/**
@brief Free the stored spherical harmonic basis functions and coefficients in
`pb`

@details Frees `pb->Q1, pb->Q2, pb->M, pb->N`.
*/
void RTS_freeHDFBeam(RTS_MWA_FEE_beam_t *pb);



/**
@brief Read in the spherical harmonic basis functions and stored coefficients
from `woden_settings->hdf5_beam_path`, and calculate the tile phases for the
given `FEE_delays` on the host device to intialise the MWA FEE beam at the
requested frequencies in `beam_freqs`.

@details The MWA primary beam is steered usings quantised delays on each of the
16 dipoles in the tile. These should be specified by
`woden_settings->FEE_ideal_delays` - they can be found in the metafits file
associated with an MWA observation. The intialised `RTS_MWA_FEE_beam_t`
strucuts will be stored in `beam_settings->FEE_beams`, and the normalisation
beams stored in `beam_settings->FEE_beam_zeniths`.

If you are using stored spherical harmonic coeffs that are at a 1.28MHz
freq resolution, this function will create multiple instances of the beam at
the closest available frequency. In that case, it's computationally FAR
cheaper to just use `RTS_MWAFEEInit` and call the beam once per coarse channel.
If you have interpolated fine-frequency spherical harmonics, then crack on.

@param[in,out] *beam_settings Initialised `beam_settings_t` struct. Output
beams are stored here
@param[in] *woden_settings Initialised `woden_settings_t` struct
@param[in] *beam_freqs Array of frequencies to set the beam up at (Hz)
@returns status Error code. 0 if all good, 1 if something is wrong
*/
int multifreq_RTS_MWAFEEInit(beam_settings_t *beam_settings,
                             woden_settings_t *woden_settings,
                             double *beam_freqs);
