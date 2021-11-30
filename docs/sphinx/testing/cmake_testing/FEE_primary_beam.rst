``FEE_primary_beam``
=========================
Tests for the functions in ``WODEN/src/FEE_primary_beam.c``.

test_RTS_MWAFEEInit.c
*********************************
``create_sky_model::RTS_MWAFEEInit`` reads in stored spherical harmonic
coefficients from ``mwa_full_embedded_element_pattern.h5`` and stores them in an ``RTS_MWA_FEE_beam_t`` struct, to be used later in MWA FEE beam calculations.
There are four generated arrays that matter, which are::

  RTS_MWA_FEE_beam_t->Q1 (double _Complex)
  RTS_MWA_FEE_beam_t->Q2 (double _Complex)
  RTS_MWA_FEE_beam_t->M (double)
  RTS_MWA_FEE_beam_t->N (double)

The MWA beam pointing direction on the sky is controlled by a set of 16 delays.
A different delay setting is stored in a different table in the ``hdf5`` file.
In these tests, three different delays settings are tested at 50MHz, 150MHz, and
250MHz (a total of nine tests). For each combination of settings, the values
of the four arrays are tested at 6 locations against stored known values (
which can be found in ``test_RTS_MWAFEEInit.h``), and must be within an absolute
tolerance of 1e-7 to pass. Only six array entries each are tested purely
to keep the stored values to compare to down a reasonable number (but suffice
as a test that things are working correctly).

test_multifreq_RTS_MWAFEEInit.c
*********************************
``create_sky_model::multifreq_RTS_MWAFEEInit`` calls ``RTS_MWAFEEInit`` for
multiple frequencies, storing the output ``RTS_MWA_FEE_beam_t`` structs in
``beam_settings->FEE_beams``. This function is designed to be used in conjunction
with the interpolated coefficient file
``MWA_embedded_element_pattern_rev2_interp_167_197MHz.h5`` which is at a 40kHz
resolution (between 167 - 197MHz). The ``mwa_full_embedded_element_pattern.h5``
file is at 1.28MHz resolution, however covers 50 - 300MHz.

This test calls the function for three different delays and frequencies,
covering a range of frequencies between 167 and 197MHz. The first test calls
32 frequencies so takes 30s per call. The first two
elements in the ``Q1`` abd ``Q2`` arrays are tested for both polarisations for
all requested frequencies, and must match to within an absolute tolerance of 1e-7
to stored values.
