#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "woden_precision_defs.h"
#include "constants.h"
#include "woden_struct_defs.h"

void read_woden_settings(woden_settings_t *woden_settings)
{

  FILE *output_text;
  output_text = fopen("woden_settings.txt","w");

  fprintf(output_text,"ra0 %.16f\n", woden_settings->ra0);
  fprintf(output_text,"dec0 %.16f\n", woden_settings->dec0);

  fprintf(output_text,"latitude %.12f\n",
                       woden_settings->latitude);
  fprintf(output_text,"latitude_obs_epoch_base %.12f\n",
                       woden_settings->latitude_obs_epoch_base);
  fprintf(output_text,"lst_base %.12f\n",
                       woden_settings->lst_base);
  fprintf(output_text,"lst_obs_epoch_base %.12f\n",
                       woden_settings->lst_obs_epoch_base);
  fprintf(output_text,"num_freqs %d\n",
                       woden_settings->num_freqs);
  fprintf(output_text,"frequency_resolution %.12e\n",
                       woden_settings->frequency_resolution);
  fprintf(output_text,"base_low_freq %.12e\n",
                       woden_settings->base_low_freq);
  fprintf(output_text,"coarse_band_width %.12e\n",
                       woden_settings->coarse_band_width);
  fprintf(output_text,"num_time_steps %d\n",
                       woden_settings->num_time_steps);
  fprintf(output_text,"time_res %.12f\n",
                       woden_settings->time_res);
  // fprintf(output_text,"cat_filename %s\n",
  //                      woden_settings->cat_filename);

  fprintf(output_text,"jd_date %.12f\n",
                       woden_settings->jd_date);
  fprintf(output_text,"sky_crop_type %d\n",
                       woden_settings->sky_crop_type);
  fprintf(output_text,"beamtype %d\n",
                       woden_settings->beamtype);

  if (woden_settings->beamtype == FEE_BEAM || woden_settings->beamtype == FEE_BEAM_INTERP || woden_settings->beamtype == MWA_ANALY) {

    fprintf(output_text, "FEE_ideal_delays %d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", woden_settings->FEE_ideal_delays[0], woden_settings->FEE_ideal_delays[1],
    woden_settings->FEE_ideal_delays[2], woden_settings->FEE_ideal_delays[3],
    woden_settings->FEE_ideal_delays[4], woden_settings->FEE_ideal_delays[5],
    woden_settings->FEE_ideal_delays[6], woden_settings->FEE_ideal_delays[7],
    woden_settings->FEE_ideal_delays[8], woden_settings->FEE_ideal_delays[9],
    woden_settings->FEE_ideal_delays[10], woden_settings->FEE_ideal_delays[11],
    woden_settings->FEE_ideal_delays[12], woden_settings->FEE_ideal_delays[13],
    woden_settings->FEE_ideal_delays[14], woden_settings->FEE_ideal_delays[15]);

  }

  fprintf(output_text,"do_precession %d\n",
                       woden_settings->do_precession);
  fprintf(output_text,"do_autos %d\n",
                       woden_settings->do_autos);
  fprintf(output_text,"chunking_size %ld\n",
                       woden_settings->chunking_size);

  fprintf(output_text,"array_layout_file_path %s\n",
                       woden_settings->array_layout_file_path);

  fprintf(output_text,"num_bands %d\n",
                       woden_settings->num_bands);
  int num_bands = woden_settings->num_bands;
  
  for (int i = 0; i < num_bands; i++)
  {
    fprintf(output_text,"bandnum%02d %d\n", i, 
                       woden_settings->band_nums[i]);
  }
  
  fprintf(output_text,"gauss_beam_FWHM %.12f\n",
           woden_settings->gauss_beam_FWHM);
  fprintf(output_text,"gauss_beam_ref_freq %.12e\n",
           woden_settings->gauss_beam_ref_freq);
  fprintf(output_text,"gauss_ra_point %.12f\n",
           woden_settings->gauss_ra_point);
  fprintf(output_text,"gauss_dec_point %.12f\n",
           woden_settings->gauss_dec_point);
  fprintf(output_text,"hdf5_beam_path %s\n",
                       woden_settings->hdf5_beam_path);

  fprintf(output_text,"use_dipamps %d\n",
                       woden_settings->use_dipamps);

  fprintf(output_text,"off_cardinal_dipoles %d\n",
                       woden_settings->off_cardinal_dipoles);

  fprintf(output_text,"do_gpu %d\n",
                       woden_settings->do_gpu);


  if (woden_settings->use_dipamps == 1) {

    fprintf(output_text, "mwa_dipole_amps %.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f\n", woden_settings->mwa_dipole_amps[0], woden_settings->mwa_dipole_amps[1],
    woden_settings->mwa_dipole_amps[2], woden_settings->mwa_dipole_amps[3],
    woden_settings->mwa_dipole_amps[4], woden_settings->mwa_dipole_amps[5],
    woden_settings->mwa_dipole_amps[6], woden_settings->mwa_dipole_amps[7],
    woden_settings->mwa_dipole_amps[8], woden_settings->mwa_dipole_amps[9],
    woden_settings->mwa_dipole_amps[10], woden_settings->mwa_dipole_amps[11],
    woden_settings->mwa_dipole_amps[12], woden_settings->mwa_dipole_amps[13],
    woden_settings->mwa_dipole_amps[14], woden_settings->mwa_dipole_amps[15]);

  }

  
  fflush(output_text);
  fclose(output_text);
    
}
