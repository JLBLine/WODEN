``calculate_visibilities``
===========================
Tests for the functions in ``WODEN/src/calculate_visibilities.cu`` (there is only one function). ``calculate_visibilities::calculate_visibilities`` is the gateway function
to all ``CUDA`` functionality in ``WODEN``. It takes in simulations settings and
a sky model, and performs the necessary coorindate and measurement equation calculations, as well as summations over sky model components to generate visibilities.

The tests below are fairly limited in the number of variables tested, as they
are integration tests to make sure all the functionality expected is called by ``calculate_visibilities::calculate_visibilities``, and is able to talk to the
GPU correctly. More rigorous testing of this functionality is included in other
test suites in ``cmake_testing``.

The tests below all use ``test_calculate_visibilities_common.c`` to setup a
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
     - 3
     - 0
     - 0
     - 9
   * - 3
     - 0
     - 3
     - 0
     - 9
   * - 3
     - 0
     - 0
     - 3
     - 9
   * - 3
     - 3
     - 3
     - 3
     - 27
  
All sources are given an *RA,Dec* equal to the phase centre, and a Stokes I
flux density of 0.3333333333333333 Jy. Each type of flux component (either
POWER_LAW, CURVED_POWER_LAW, or LIST) are set up so they have a perfectly flat
spectrum. 

In addition to the Stokes I parameters above, each component is given a linear
and circular polarisation component. Again, the polarised components are given
flat spectras and a reference flux density of 0.3333333333333333 Jy. The linear
polarisation is given an RM of zero, meaning Stokes Q is equal to the linear
polarisation flux, and Stokes U is zero. All versions of the polarisation
models are tested (power-law, curved power-law, and polarisation fraction).

This way, the expected sum of Stokes IQUV are simply multiples of the number of
components in the model, modulu the beam response (see below).

.. note::

	Given all sources are at phase centre, and the beam model is the same for all antennas, the cross-correlations and the auto-correlations should have exactly the same values. So in all tests described below, ``calculate_visibilities`` is run with both ``woden_settings->do_autos`` set to 0 and to 1. If the autos are calculated, the number of output visibilities should double (as I only simulate three baselines from three antennas here). The latter half are where the autos are set, and are tested to be equal to the expected values for the crosses.

When the tests have a single COMPONENT per SOURCE, I use a POWER_LAW type flux
behaviour. When there are three COMPONENTs, I use one of each of POWER_LAW,
CURVED_POWER_LAW, and LIST types.

GAUSSIAN and SHAPELET components with any size cause a reduction of the real part
of the visibility due to the extended size on the sky. For these tests, I've set
their major and minor axes to 1e-10 to make them behave similarly to point sources.

For all tests below, I setup a simulation with three baselines, three frequencies,
and two timesteps. For all sky models, frequencies, time steps, and baselines, I check:

 - The *u,v,w* are as expected
 - The XX,XY,YX,YY visibilities are as expected. The expected visibilities are calculated by predicting the Stokes IQUV values as explained above, and then applying stored beam gain values to predict the visibilities. These beam gain values are taken from the DOUBLE compiled versions of the repsective beam models

To keep testing the *u,v,w* straight forward, I've set the baseline lengths in :math:`X` and :math:`Y` equal, (i.e. :math:`X_{\mathrm{diff}} = X_{\mathrm{ant1}} - X_{\mathrm{ant2}} = Y_{\mathrm{diff}}`), and the length in :math:`Z` to zero. With this configuration, the
following is true:

 - :math:`u = X_{\mathrm{diff}}(\cos(ha_0) + \sin(ha_0))`
 - :math:`v = X_{\mathrm{diff}}\sin(\phi_{\mathrm{lat}})(-\cos(ha_0) + \sin(ha_0))`
 - :math:`w = X_{\mathrm{diff}}\cos(\phi_{\mathrm{lat}})(\cos(ha_0) - \sin(ha_0))`

where :math:`ha_0` is the hour angle of the phase centre, and :math:`\phi_{\mathrm{lat}}`
the latitude of the array. The allows us to check the *u,v,w* are changing with time.

For FLOAT compiled code, the absolute tolerance threshold on the *u,v,w*
values is set to 1e-5, and 1e-12 for DOUBLE compiled code.

``test_calculate_visibilities_nobeam.c``
*********************************************
Runs tests above, but with no primary beam model applied. Predicted gains are 
ones with zero leakage.

``test_calculate_visibilities_gaussbeam.c``
*********************************************
Runs tests using the Gaussian primary beam model. No leakage in this model.

``test_calculate_visibilities_edabeam.c``
*********************************************
Runs tests using the EDA2 primary beam model. No leakage in this model.

``test_calculate_visibilities_mwafeebeam.c``
*********************************************
Runs tests using the MWA FEE primary beam model. Has gains and leakage terms.

``test_calculate_visibilities_mwafeebeaminterp.c``
****************************************************
Runs tests using the interpolated MWA FEE primary beam model. Has gains and leakage terms.

``test_calculate_visibilities_mwaanalybeam.c``
****************************************************
Runs tests using the interpolated MWA FEE primary beam model. Has gains and leakage terms,
but is a real-only model.

``test_calculate_visibilities_multibeams.c``
*********************************************
Same as ``test_calculate_visibilities_mwafeebeaminterp.c``, but where every antenna has a different primary beam. Uses dipole amplitude gains that cause predictable scalar multiplications of expected visibilities from the previous tests. Doing indexing and keeping track of which primary beams should match which visibilities allows predictions of the expected visibilities to be made.