#include <stdio.h>

void print_cmdline_help() {
  printf("Must input a .json settings file to run woden. This file must include:\n");
  printf("\tra0: ra phase centre (float in degrees)\n");
  printf("\tdec0: dec phase centre (float in degrees)\n");
  printf("\tnum_freqs: number of fine frequency channels to simulate (int)\n");
  printf("\tnum_time_steps: number of time steps to simulate (int)\n");
  printf("\tcat_filename: path to and name of WODEN-style srclist (string)\n");
  printf("\tmetafits_filename: path to MWA and name of metafits file to base\n\t\tsimulation on (string)\n");
  printf("\n");
  printf("Optionally, the .json can include:\n");
  printf("\tsky_crop_components=True: WODEN crops sources with any component\n\t\tbelow the horizon. Add this arg to include all components\n\t\tabove horizon, regardless of which source they belong to\n");
  printf("\tuse_gaussian_beam=True: Use a gaussian primary beam in the simulation\n\t\tpointed at the RA,DEC pointing specified in the metafits file\n");
  printf("\tgauss_beam_FWHM: Sets the FWHM of the gaussian primary beam (degrees)\n");
  printf("\tgauss_beam_ref_freq: The frequency (Hz) at which gauss_beam_FWHM is set\n\t\tThe beam FWHM will scale with frequency about this reference\n");
}
