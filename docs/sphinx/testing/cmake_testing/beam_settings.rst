``beam_settings``
=========================
Tests for the functions in ``WODEN/src/beam_settings.c``. These functions
setup primary beam settings, ready to calculate beam responses CPU or GPU.

test_fill_primary_beam_settings.c
***********************************
``primary_beam::fill_primary_beam_settings`` prepares a ``beam_settings_t``
struct to be used by ``calculate_visibilities::calculate_visibilities``. Most
beam types just require that ``beam_settings->beamtype`` is set, but ``GAUSS_BEAM`` requires
some further settings. So in the case of ``GAUSS_BEAM``, check that::

   beam_settings->beam_FWHM_rad
   beam_settings->beam_ref_freq
   beam_settings->gauss_sdec
   beam_settings->gauss_cdec
   beam_settings->gauss_ha

are set correctly.