#include <stdio.h>

void print_cmdline_help() {
  printf("Must input a .json settings file to run woden, e.g.\n");
  printf("\t $ woden woden_settings.json\n");
  printf("\n");

  printf("This woden_settings.json file must include:\n");
  printf("\t+ ra0: RA phase centre (float in degrees)\n");
  printf("\t+ dec0: Dec phase centre (float in degrees)\n");
  printf("\t+ num_freqs: number of fine frequency channels to simulate (int)\n");
  printf("\t+ num_time_steps: number of time steps to simulate (int)\n");
  printf("\t+ cat_filename: path to and name of WODEN-style srclist (string)\n");
  printf("\t+ time_res: Time resolution (s) of the simulation\n");
  printf("\t+ frequency_resolution: Fine channel frequnecy resolution (Hz)\n");
  printf("\t+ jd_date: Julian date of the start of the observation\n");
  printf("\t+ LST: Local Sidereal Time at the start of the observation (radians)\n");
  printf("\t+ lowest_channel_freq: Frequency of the lowest fine channel of\n");
      printf("\t\tband 1 in the simulation \n");
  printf("\t+ coarse_band_width: the frequency bandwidth of each band (Hz)\n");
  printf("\t+ band_nums: an array of which band numbers to simulate (ints), e.g.\n");
      printf("\t\t band_nums: [1,4,7]; \n");
  printf("\t+ coarse_band_width: the frequency bandwidth of each band (Hz)\n");
  printf("\t+ array_layout: Path to a text file containing local E,N,H antenna\n");
      printf("\t\tlocations (meters)\n");

  printf("\n");

  printf("Optionally, the .json can include:\n");
  printf("\t+ latitude: Latitude of the array (deg) - defaults to MWA\n");
  printf("\t+ use_FEE_beam: use the MWA FEE primary beam model\n");
  printf("\t+ hdf5_beam_path: Location of the hdf5 file holding the FEE beam\n");
      printf("\t\t coefficients\n");
  printf("\t+ FEE_delays: A array of 16 delays to point the MWA FEE primary beam\n");
      printf("\t\t model\n");

  printf("\t+ use_gaussian_beam=True: Use a gaussian primary beam in the\n");
      printf("\t\tsimulation pointed at the RA,DEC specified by --gauss_ra_point\n");
      printf("\t\tand --gauss_dec_point\n");
  printf("\t+ gauss_beam_FWHM: Sets the FWHM of the gaussian primary beam (degrees)\n");
  printf("\t+ gauss_beam_ref_freq: The frequency (Hz) at which gauss_beam_FWHM is\n");
      printf("\t\tsetThe beam FWHM will scale with frequency about this reference\n");
  printf("\t+ gauss_ra_point: The initial RA (deg) to point the Gaussian beam at\n");
  printf("\t+ gauss_dec_point: he initial Dec (deg) to point the Gaussian beam at\n");

  printf("\t+ use_EDA2_beam=True: Use the EDA2 beam (Analytic dipole with a ground mesh)\n");

  printf("\t+ sky_crop_components=True: WODEN crops sources with any component\n");
      printf("\t\tbelow the horizon. Add this arg to include all components\n");
      printf("\t\tabove horizon, regardless of which source they belong to\n");
  printf("\t+ chunking_size: The chunk size to break up the point sources into \n");
      printf("\t\tfor processing\n");
}
