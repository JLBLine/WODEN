``calculate_visibilities``
===========================
Tests for the functions in ``WODEN/src/calculate_visibilities_common.c``.
``calculate_visibilities_common::calculate_visibilities`` is the gateway function
to all visibility calculations in ``WODEN``. It takes in simulations settings and
a sky model, and performs the necessary coordindate and measurement equation
calculations, as well as summations over sky model components to generate visibilities.

Depending on the settings, these calculations are be done on the CPU or GPU.
This allows us to use the same test code for the CPU and GPU, as indicated by the
``gpu.c`` and ``cpu.c`` suffixes, ensuring the CPU and GPU versions are consistent.
Each ``*gpu.c/*cpu.c`` uses a different primary beam model, calling the matching
``calculate_visibilities_*common.c`` function. For example,
``test_calculate_visibilities_gaussbeam_cpu.c`` calls
``calculate_visibilities_gaussbeam_common.c``, which sets up running with a
Gaussian. All tests use some common functionality stored in  ``calculate_visibilities_common_common.c``.

The tests below are fairly limited in the number of variables tested, as they
are integration tests to make sure all the functionality expected is called by
``calculate_visibilities_common::calculate_visibilities``, and is able to talk to the
CPU/GPU correctly. More rigorous testing of this functionality is included in other
test suites in ``cmake_testing``.

The tests below all use ``calculate_visibilities_common_common.c`` to setup a
``source_catalogue_t`` sky model struct with the following sky model setups.
These differing sky models should cover all COMPONENT combinations possible, to
make sure ``calculate_visibilities`` is calling the appropriate functions:

.. list-table::
   :widths: 30 30 30 30 30
   :header-rows: 1

   * - Num SOURCEs
     - Num POINT in each SOURCE
     - Num GAUSS in each SOURCE
     - Num SHAPELET in each SOURCE
     - Total COMPONENTS
   * - 1
     - 1
     - 0
     - 0
     - 1
   * - 1
     - 0
     - 1
     - 0
     - 1
   * - 1
     - 0
     - 0
     - 1
     - 1
   * - 1
     - 1
     - 1
     - 1
     - 3
   * - 3
     - 1
     - 0
     - 0
     - 3
   * - 3
     - 0
     - 1
     - 0
     - 3
   * - 3
     - 0
     - 0
     - 1
     - 3
   * - 3
     - 1
     - 1
     - 1
     - 9
   * - 3
     - 5
     - 0
     - 0
     - 15
   * - 3
     - 0
     - 5
     - 0
     - 15
   * - 3
     - 0
     - 0
     - 5
     - 15
   * - 3
     - 5
     - 5
     - 5
     - 45
  
All sources are given an *RA,Dec* equal to the phase centre, and a Stokes I
flux density of 0.3333333333333333 Jy. Each type of flux component (either
POWER_LAW, CURVED_POWER_LAW, or LIST) are set up so they have a perfectly flat
spectrum. 

In addition to the Stokes I parameters above, each component is given a linear
and circular polarisation component. Again, the polarised components are given
flat spectras and a reference flux density of 0.3333333333333333 Jy. The linear
polarisation is given an RM of zero, meaning Stokes Q is equal to the linear
polarisation flux, and Stokes U is zero. All versions of the polarisation
models are tested (power-law, curved power-law, polarisation fraction, list-type).

This way, the expected sum of Stokes IQUV are simply multiples of the number of
components in the model, modulu the beam response (see below). This can be easily
combined with stored beam gains to predict the visibilities. All tests are run
with and without auto-correlations.

When the tests have a single COMPONENT per SOURCE, I use a POWER_LAW type flux
behaviour. When there are three COMPONENTs, I use one of each of POWER_LAW,
CURVED_POWER_LAW, and LIST types.

GAUSSIAN and SHAPELET components with any size cause a reduction of the real part
of the visibility due to the extended size on the sky. For these tests, I've set
their major and minor axes to 1e-10 to make them behave similarly to point sources.

For all tests below, I setup a simulation with three baselines, three frequencies,
and two timesteps. For all sky models, frequencies, time steps, and baselines, I check:

 - The *u,v,w* are as expected
 - The XX,XY,YX,YY visibilities are as expected. The expected visibilities are calculated by predicting the Stokes IQUV values as explained above, and then applying stored beam gain values to predict the visibilities. These beam gain values are taken from the DOUBLE compiled versions of the respective beam models.

While I run with three frequency steps, the frequency resolution is set to close
to zero, so the beam response is the same for all frequencies. This allows me
to test the correct dimensions of the output visibilities, but only store one
set of beam gains.

To keep testing the *u,v,w* straight forward, I've set the baseline lengths in :math:`X` and :math:`Y` equal, (i.e. :math:`X_{\mathrm{diff}} = X_{\mathrm{ant1}} - X_{\mathrm{ant2}} = Y_{\mathrm{diff}}`), and the length in :math:`Z` to zero. With this configuration, the
following is true:

 - :math:`u = X_{\mathrm{diff}}(\cos(ha_0) + \sin(ha_0))`
 - :math:`v = X_{\mathrm{diff}}\sin(\phi_{\mathrm{lat}})(-\cos(ha_0) + \sin(ha_0))`
 - :math:`w = X_{\mathrm{diff}}\cos(\phi_{\mathrm{lat}})(\cos(ha_0) - \sin(ha_0))`

where :math:`ha_0` is the hour angle of the phase centre, and :math:`\phi_{\mathrm{lat}}`
the latitude of the array. The allows us to check the *u,v,w* are changing with time.

For FLOAT compiled code, the absolute tolerance threshold on the *u,v,w*
values is set to 1e-5, and 1e-12 for DOUBLE compiled code.

``test_calculate_visibilities_nobeam*.c``
*********************************************
Runs tests above, but with no primary beam model applied. Predicted gains are 
ones with zero leakage.

``test_calculate_visibilities_gaussbeam*.c``
*********************************************
Runs tests using the Gaussian primary beam model. No leakage in this model.

``test_calculate_visibilities_edabeam*.c``
*********************************************
Runs tests using the EDA2 primary beam model. No leakage in this model.

``test_calculate_visibilities_mwafeebeam*.c``
*********************************************
Runs tests using the MWA FEE primary beam model. Has gains and leakage terms.

``test_calculate_visibilities_mwafeebeaminterp*.c``
****************************************************
Runs tests using the interpolated MWA FEE primary beam model. Has gains and leakage terms.

``test_calculate_visibilities_mwaanalybeam*c``
****************************************************
Runs tests using the interpolated MWA FEE primary beam model. Has gains and leakage terms,
but is a real-only model.

``test_calculate_visibilities_multibeams*.c``
*********************************************
Same as ``test_calculate_visibilities_mwafeebeaminterp.c``, but where every antenna (tile) has a different primary beam. Uses dipole amplitude gains that cause predictable scalar multiplications of expected visibilities from the previous tests. Doing indexing and keeping track of which primary beams should match which visibilities allows predictions of the expected visibilities to be made. 


``test_calculate_visibilities_everybeam_lofar_*.c``
*******************************************************
Runs with the EveryBeam ``LOFAR`` primary beam model. Has gains and leakage terms. Runs
with a different model for each station. Here we intentionally point the 
primary beam slightly off zenith, so that the beam response is different for
each station. If we left it locked to the phase centre, the beam response would
gains of 1 and leakage of zero. Use stored beam gain values to calculate 
expected visibilities to compare results against.


``test_calculate_visibilities_everybeam_mwa_*.c``
*******************************************************
Runs with the EveryBeam ``MWA`` primary beam model. Has gains and leakage terms. Runs
with a single primary beam, and I haven't worked out how to change dipole amplitudes
with EveryBeam. Use stored beam gain values to calculate expected visibilities
to compare results against.


``test_calculate_visibilities_everybeam_oskar_*.c``
*******************************************************
Runs with the EveryBeam ``OSKAR`` primary beam model. Has gains and leakage terms. Runs
with a different model for each station. Here we intentionally point the 
primary beam slightly off zenith, so that the beam response is different for
each station. If we left it locked to the phase centre, the beam response would
gains of 1 and leakage of zero. Use stored beam gain values to calculate 
expected visibilities to compare results against.

.. ``make_exe_to_profile_lofar_everybeam.c``
.. **********************************************
.. This creates an executable that can be used to profile the EveryBeam LOFAR beam model.
.. It calculates visibilities for 200 