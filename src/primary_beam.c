#include "primary_beam.h"
#include <complex.h>
#include <math.h>

beam_settings_t * fill_primary_beam_settings(woden_settings_t *woden_settings,
                                             double *lsts) {


  //Setup primary beam settings for observation
  beam_settings_t *beam_settings = malloc(sizeof(beam_settings_t));

  //Both GAUSS_BEAM needs a bunch of extra settings
  if (woden_settings->beamtype == GAUSS_BEAM) {
    
    //Extra settings that just GAUSS_BEAM needs
    beam_settings->beamtype = GAUSS_BEAM;

    //Angles used in calculating beam centred l,m,ns
    beam_settings->gauss_sdec = sin(woden_settings->gauss_dec_point);
    beam_settings->gauss_cdec = cos(woden_settings->gauss_dec_point);
    beam_settings->gauss_ha = woden_settings->lst_base - woden_settings->gauss_ra_point;

    printf("Setting up Gaussian primary beam settings\n");
    printf("   pointing at HA, Dec = %.5fdeg, %.5fdeg\n",
                beam_settings->gauss_ha/DD2R, woden_settings->gauss_dec_point/DD2R );
    printf("   setting beam FWHM to %.5fdeg and ref freq to %.3fMHz\n",
            woden_settings->gauss_beam_FWHM,woden_settings->gauss_beam_ref_freq / 1e+6  );

    //Set constants used in beam calculation
    beam_settings->beam_FWHM_rad = woden_settings->gauss_beam_FWHM * DD2R;
    beam_settings->beam_ref_freq = woden_settings->gauss_beam_ref_freq;

  } // End if (woden_settings->gaussian_beam)

  else if (woden_settings->beamtype == MWA_ANALY ){
      beam_settings->beamtype = MWA_ANALY;
  }

  else if (woden_settings->beamtype == FEE_BEAM) {
    beam_settings->beamtype = FEE_BEAM;
  }

  else if (woden_settings->beamtype == FEE_BEAM_INTERP) {
    beam_settings->beamtype = FEE_BEAM_INTERP;
  }

  else if (woden_settings->beamtype == ANALY_DIPOLE) {
    beam_settings->beamtype = ANALY_DIPOLE;
  }

  else if (woden_settings->beamtype == EB_OSKAR) {
    beam_settings->beamtype = EB_OSKAR;
  }

  else if (woden_settings->beamtype == EB_LOFAR) {
    beam_settings->beamtype = EB_LOFAR;
  }

  else if (woden_settings->beamtype == EB_MWA) {
    beam_settings->beamtype = EB_MWA;
  }

  else {
    printf("NO PRIMARY BEAM HAS BEEN SELECTED\n\tWill run without a primary beam\n");
    beam_settings->beamtype = NO_BEAM;
  }

  return beam_settings;

} // end of fill_primary_beam_settings()
