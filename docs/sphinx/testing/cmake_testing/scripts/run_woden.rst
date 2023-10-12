``run_woden.py``
===========================
``run_woden.py`` is the main ``WODEN`` executable. It should parse input command arguments, gather metafits and sky model information, calculate observational parameters, launch GPU simulation code, and output a ``uvfits`` file.


test_run_woden.py
***************************
This test writes out some example ``uvfits`` files with known contents, creates a mock of a command line argument object, and runs ``add_woden_uvfits.main`` to create combined outputs. These outputs are checked to be the sum of the inputs. A number of input arguments that should cause failures, including attempting to add to uvfits with different dimensions, are also tested.

To make things easy to predict, the following inputs are used:
 - The array is placed at latitude of 0.0$^\circ$, longitude 79.53789457709424$^\circ$. This gives an LST of 0.0$^\circ$ on the date J2000 date 2000-01-01 12:00:00
 - Only three antennas, with *e,n,h* positions (0,0,0), (10,0,0), (20,0,0) are used. Not only does this reduce the number of baselines, but making them east-west only means we can predict the *u,v,w*, as only *u* should have values based on the differences in *e* (being -10,-20,-10)
 - Set the sky model to be 5 point sources, each with 1 Jy and a flat SED, located at the phase centre. This way all visibilities should just be real valued and equal to 5. Turn off the primary beam as well to ensure this is the case; this means XX,YY should have value, XY,YZ should all be zero.

As this is more of an integration test, just test the u,v,w coordinates and visibilities in the output uvfits are correct (tested to be within an absolute value of 1e-12 metres or Jy). Other unit tests cover more specific functionality. We use this test to also check coverage for ``wodenpy.use_libwoden.visibility_set``, so we also test with and without the ``--use_auto`` option, and test both the ``float`` and ``double`` precision options.