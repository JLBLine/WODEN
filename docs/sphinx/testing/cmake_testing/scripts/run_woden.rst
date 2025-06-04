``run_woden.py``
===========================
``run_woden.py`` is the main ``WODEN`` executable. It should parse input command arguments, gather metafits and sky model information, calculate observational parameters, launch GPU simulation code, and output a ``uvfits`` file.


test_get_future_result.py
***************************
Within ``run_woden.py``, unless in serial mode, we run two ``concurrent.futures.ProcessPoolExecutor``s to first run sky model reading in parallel, and then either launch a GPU process, or CPU ctypes processes in parallel. This tests check the ``get_future_result``, which handles recovering outputs from the ``future`` objects that are output by the ``ProcessPoolExecutor``. Test by launching simple Python functions in one ``ProcessPoolExecutor``, and some simple ``C`` in a second. Big motivation here is to be able to actually recover error messages, as they are often obscured by the ``concurrent.futures`` module. Test we can get outputs / errors from the following combinations:
 - Python error in first ``ProcessPoolExecutor``
 - Python error in second ``ProcessPoolExecutor``
 - C error in first ``ProcessPoolExecutor``
 - C error in second ``ProcessPoolExecutor``
 - No errors and expected outputs from both ``ProcessPoolExecutor``s for Python functions
 - No errors and expected outputs from both ``ProcessPoolExecutor``s for C functions


test_run_woden.py
***************************
This test runs with a variety of settings to ensure the ``run_woden.py`` script is working correctly.

To make things easy to predict, the following inputs are used:
 - The array is placed at latitude of 0.0$^\circ$, longitude 79.53789457709424$^\circ$. This gives an LST of 0.0$^\circ$ on the date J2000 date 2000-01-01 12:00:00
 - Only three antennas, with *e,n,h* positions (0,0,0), (10,0,0), (20,0,0) are used. Not only does this reduce the number of baselines, but making them east-west only means we can predict the *u,v,w*, as only *u* should have values based on the differences in *e* (being -10,-20,-10)
 - Set the sky model to be 5 point sources, each with 1 Jy and a flat SED, located at the phase centre. This way all visibilities should just be real valued and equal to 5. Turn off the primary beam as well to ensure this is the case; this means XX,YY should have value, XY,YZ should all be zero.

As this is more of an integration test, just test the u,v,w coordinates and (sometimes) visibilities in the output uvfits are correct. Other unit tests cover more specific functionality. Tests are run with the following settings (if visibilities aren't checked, noted below):
 - No primary beam, float/double precision, number threads 1, cpu/gpu mode (so 4 tests)
 - No primary beam with auto-correlations, number threads 1, cpu/gpu mode
 - No primary beam with auto-correlations, number threads 3, multiple coarse bands, cpu/gpu mode
 - No primary beam, profiler turned on, num threads 1
 - No primary beam, num threads 1, verbose logging on, saved log (checks file exists)
 - No primary beam, num threads 6, checks things work if a sky model thread has nothing to process
 - UVBeam MWA primary beam, num threads 1, cpu mode. Does not check outputs, just that it runs
 - No primary beam, num threads 2, gpu mode. Checks that multiple rounds of GPU processing work by setting chunking size to 1.
 - EveryBeam LOFAR primary beam. Checks things run when moving the LOFAR array centre to MWA. Does not check outputs, just that it runs
 - UVBeam HERA beam when using a FITS model. Does not check outputs, just that it runs