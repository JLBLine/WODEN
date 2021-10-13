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

All sources are given an *RA,Dec* equal to the phase centre, a spectral index
of zero, and a flux density of 1.0 Jy. This way, the output visibilities (in the absence
of a primary beam) should be fully real, and equal to the sum of the number of
COMPONENTs in the sky model. GAUSSIAN and SHAPELET components with any size would
cause a loss in flux, so I've set their major and minor axes to a minuscule number so
they behave like point sources.

For all tests below, I setup a simulation with three baselines, three frequencies,
and two timesteps. For all sky models, frequencies, timesteps, and baselines, I check:

 - The *u,v,w* are as expected
 - The visibilities are equal to the sum of the COMPONENTs (modulu the beam repsonse - see below)

To keep testing the *u,v,w* straight forward, I've set the baseline lengths in :math:`X` and :math:`Y` equal, (i.e. :math:`X_{\mathrm{diff}} = X_{\mathrm{ant1}} - X_{\mathrm{ant2}} = Y_{\mathrm{diff}}`), and the length in :math:`Z` to zero. With this configuration, the
following is true:

 - :math:`u = X_{\mathrm{diff}}(\cos(ha_0) + \sin(ha_0))`
 - :math:`v = X_{\mathrm{diff}}\sin(\phi_{\mathrm{lat}})(-\cos(ha_0) + \sin(ha_0))`
 - :math:`w = X_{\mathrm{diff}}\cos(\phi_{\mathrm{lat}})(\cos(ha_0) - \sin(ha_0))`

where :math:`ha_0` is the hour angle of the phase centre, and :math:`\phi_{\mathrm{lat}}`
the latitude of the array. The allows us to check the *u,v,w* are changing with time.

``test_calculate_visibilities_nobeam.c``
*********************************************
This runs the tests explained above, whilst switching the primary beam off. This
really does check that real visibilities are equal to the number of COMPONENTs for all 12
sky model configurations, and the imaginary equal to zero. A tolerance of within 0.1% of the
expected value is applied when checking the visibilities (the same tolerance is applied for
the following test files as well).

``test_calculate_visibilities_gaussbeam.c``
*********************************************
This runs the same tests as ``test_calculate_visibilities_nobeam.c``, but applies
a Gaussian primary beam model. As all sky model COMPONENTs are set at the same location,
only one beam gain per time step should be applied to visibilities, so for each time
step, we can check whether the visibilities equal this gain times the number of
COMPONENTs. The Gaussian beam is fully real and has no cross-pol values, so only
check that the real XX and YY visibilities have value, and ensure all other
visibility information is zero.

``test_calculate_visibilities_edabeam.c``
*********************************************
Runs exactly the same tests as ``test_calculate_visibilities_gaussbeam.c``, but
using the analytic single dipole beam (the EDA2 beam). Obviously tests the
visiblities match the gain values that are expected for the EDA2 beam and not
the Gaussian Beam test.

``test_calculate_visibilities_mwafeebeam.c``
*********************************************
Again, runs the same tests as ``test_calculate_visibilities_gaussbeam.c``, but
this time for the coarse resolution MWA FEE primary beam. As this model is
complex and includes mutual coupling, both the real and imaginary values
for all XX, XY, YX, and YY polarisations are tested.
