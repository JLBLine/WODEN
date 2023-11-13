``add_instrumental_effects_woden.py``
=======================================
``add_instrumental_effects_woden.py`` can be used after the  main ``WODEN`` executable to add instrumental effects like noise and antenna gains.


test_add_instrumental_effects_woden.py
****************************************
This test suite writes out some example ``uvfits`` files with all visibilities equal to $1 + 0j$. This makes it straight forward to predict the effects of adding instrumental effects like noise and antenna gains. All tests are detailed below and are run by simulated command line arguments to feed into ``add_instrumental_effects_woden.main``, and set to run on the example ``uvfits`` file. The instrumentally corrupted ``uvfits`` files are then read in and checked to make sure the correct effects have been added.

 - **Noise** - runs ``add_instrumental_effects_woden.main`` to add noise to the example ``uvfits`` file. Then reads in effected ``uvfits`` file, and measures the standard deviation of the real and imaginary components for all 4 pols. It then checks this matches what is predicted by the radiometric noise equation, to within 5%, which is good enough given the small number of visibilities created in this test. Auto-correlations are also tested but to within 10% as there are less samples than the cross-correlations.
 - **Antenna gains and leakages**: This test generates the gains and leakages it expects ``add_instrumental_effects_woden.py`` to create, and uses those gains and leakages to create expected visibilities. Instead of using matrix multiplication, it uses explicit equations for each instrumental pol, to check the two methods agree. As the gain amplitudes, phases, and leakages can be generated separately, and number of different combination of input arguments are tested through ``add_instrumental_effects_woden.main``. It also checks that antenna gains are added to the auto-correlations as well as the cross-correlations.
 - **Antenna gains and leakages**: