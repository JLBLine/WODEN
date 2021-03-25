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
#include "woden_struct_defs.h"

/**
@brief Operator function to be called by H5Literate.

@details Used when finding the closest available frequency in the MWA FEE
primary beam model
*/
struct opdata {
  float freq_in;
  float freq_out;
  float least_diff;
};

/**
@brief Operator function to be called by `H5Literate`, finds closest frequency
in the element pattern hdf5 file

@details The MWA FEE beam is stored at a 1.28MHz frequency resolution. Using
this function is conjunction with `H5Literate`, find the frequencies stored
in the MWA FEE hdf5, and return the closest frequency to the frequency in
`operator_data.freq_in` (Hz)

*/
herr_t RTS_op_func(hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data);

int RTS_HDFBeamInit(const char *h5filename, float freq_Hz, copy_primary_beam_t *pb,
                float *FEE_delays, int stn);

void RTS_P1SIN(int nmax, float theta, double _Complex **P_sin, double _Complex **P1);

int RTS_get_FF2(float phi, float theta, copy_primary_beam_t *pb, float _Complex result[4], int scaling, int clean_up);

int RTS_getJonesSphHarm(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);

int RTS_getJonesDipole(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);

int RTS_getTileResponse(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4], int scaling, float rotation);

void RTS_freeHDFBeam(copy_primary_beam_t *pb);

void RTS_HDFBeamCleanUp();
