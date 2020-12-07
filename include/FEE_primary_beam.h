#pragma once
#include "hdf5.h"
#include "hdf5_hl.h"
#include "constants.h"
#include "woden_struct_defs.h"


/***************************************************************
Data Structure used by H5Literate.
***************************************************************/
struct opdata {
  float freq_in;
  float freq_out;
  float least_diff;
};

/*
 * Operator function to be called by H5Literate.
 */
herr_t RTS_op_func (hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data);

int RTS_HDFBeamInit(const char *h5filename, float freq_Hz, copy_primary_beam_t *pb,
                float *FEE_delays, int stn);

void RTS_P1SIN(int nmax, float theta, double _Complex **P_sin, double _Complex **P1);
//
int RTS_get_FF2(float phi, float theta, copy_primary_beam_t *pb, float _Complex result[4], int scaling, int clean_up);
//
int RTS_getJonesSphHarm(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);
//
int RTS_getJonesDipole(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);
//
int RTS_getTileResponse(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4], int scaling, float rotation);
//
void RTS_freeHDFBeam(copy_primary_beam_t *pb);
//
void RTS_HDFBeamCleanUp();
