``woden_settings``
=========================
Tests for the functions in ``WODEN/src/woden_settings.c``. These functions handle
reading input simulation settings from a ``.json`` file, and filling a
``woden_settings_t`` struct based on them. ``woden_settings_t`` is passed
around by most functions internally to ``WODEN`` to tell them what to do.

``test_read_json_settings.c``
******************************
Tests the function ``woden_settings::read_json_settings``. This function
reads in the simulation settings held in a ``.json`` file (usually written by
``run_woden.py``). This test runs reads multiple ``.json`` files in, each
with different settings. All test ``.json`` files share these common core
settings:

.. code-block:: json

   {
   "ra0": 0.0000000000,
   "dec0": -27.0000000000,
   "num_freqs": 16,
   "num_time_steps": 4,
   "cat_filename": "srclist_singlepoint.txt",
   "time_res": 2.0,
   "frequency_resolution": 40000.0,
   "chunking_size": 5000,
   "jd_date": 2457278.2010995,
   "LST": 0.44312771,
   "array_layout": "example_array_layout.txt",
   "lowest_channel_freq": 1.6703500000e+08,
   "latitude": -26.70331944,
   "coarse_band_width": 1.2800000000e+06,
   "sky_crop_components": "True",
   "band_nums": [1,4,9]
   }

For all tests, these core attributes are tested as being read in correctly into
``woden_settings_t`` after calling ``read_json_settings``. The test
files listed in this table are all run to test different simulation setups and
eventualities.

.. list-table::
   :widths: 25 50
   :header-rows: 1

   * - run_woden_nobeam.json
     - The NO_BEAM primary beam is selected
   * - run_woden_EDA2.json
     - The ANALY_DIPOLE primary beam is selected
   * - run_woden_gaussian_bespoke.json
     - The GAUSSIAN primary beam is selected, using input values to set the FWHM and reference frequency
   * - run_woden_gaussian_default.json
     - The GAUSSIAN primary beam is selected, using default values for the FWHM and reference frequency
   * - run_woden_MWAFEE.json
     - The MWA_FEE primary beam is selected, and associated delays and path to hdf5 file are read in correctly
   * - run_woden_MWAFEE_interp.json
     - The frequency interpolated MWA_FEE primary beam is selected, and associated delays and path to hdf5 file are read in correctly
   * - run_woden_MWA_analy.json
     - The analytic MWA_primary beam is selected, and associated delays are read in correctly
   * - run_woden_MWAFEE_baddelay.json
     - Contains a bad set of MWA FEE delays and should throw an error
   * - run_woden_MWAFEE_nopath.json
     - Has no path to the MWA FEE hdf5 file so should throw an error
   * - run_woden_multiple_beams.json
     - Contains multiple primary beam selections so should throw an error
   * - run_woden_noprecession.json
     - Checks that precession is switched off if requested

``test_setup_lsts_and_phase_centre.c``
*****************************************
Tests ``woden_settings::setup_lsts_and_phase_centre`` which uses the input
simulation parameters to calculate the sine and cosine of the declination of
the phase centre, and the LST at the centre of every time integration. Runs
three tests, two with :math:`\mathrm{Dec}_{\mathrm{phase}} = \pi/2`, where one has a
zero initial LST, another with non-zero initial LST, and third test with
:math:`\mathrm{Dec}_{\mathrm{phase}}` set to the latitude of the MWA and a
non-zero initial LST. All tests have multiple time steps, and are tested against
equivalent calculations made in ``C`` with 64 bit precision. The FLOAT code
outputs are tested with an absolute tolerance of 1e-7 and the DOUBLE a tolerance
of 1e-15.
