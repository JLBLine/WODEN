``primary_beam``
=========================
Tests for the functions in ``wodenpy/primary_beam``, which are either wrappers around ``C++`` to run the EveryBeam primary beam code, or wrappers around the ``pyuvdata UVBeam`` code.


test_calc_uvbeam_for_components.py
***********************************************************
Given an initialised ``Components_Python`` objects, and an array of ``UVBeam`` objects, this function calculates the ``UVBeam`` primary beam response for each components, storing the results in the ``Components_Python`` object. Setup 5 RA/Dec locations, and test output gains against expected values for the ``HERA`` (from both a set of CST files and a FITS file), and ``MWA`` primary beams. One RA/Dec coord is set to be beam centre, so also check the gains here are close to one, and leakages close to zero. 

test_check_ms_telescope_type_matches_element_response.py
***********************************************************
This function checks that the telescope type in a measurement set allows for the requested element response. For example, the only type of primary beam model that can be used with the MWA telescope is the ``MWA`` beam model. Checks by reading in three measurement sets with ``MWA``, ``LOFAR``, and ``OSKAR`` telescope types. Asserts function fails if an incorrect element response is requested, and passes if an appropriate element response is requested. Also writes out an unknown telescope type measurement set to check that fails as well.

test_create_filtered_ms.py
***********************************************************
``create_filtered_ms.py`` is used to create a minimalistic copy of a measurement set, changing the desired phase centre, and optionally moving the array centre to a new latitude/longitude. As measurement sets are used exclusively for running ``EveryBeam``, test by running the ``EveryBeam`` code on the filtered measurement set. First, test by moving a LOFAR measurement set to the same array centre. Run the primary beam with the old and new measurement sets, ensuring the calculations used to re-normalise the faux-zenith of each station don't change the beam response. Then, move the LOFAR measurement set to the MWA array centre, and run the primary beam again with identical azmimuth/zenith angles. Again, ensure the beam response doesn't change. Second test ensures the beam pointing is changed as well as the correct array centre move is performed.


test_run_everybeam_OSKAR.py
***********************************************************
We test the MWA ``EveryBeam`` code against the ``hyperdrive`` code below, which gives us confidence that the ``EveryBeam`` MWA code is working correctly. This test runs the ``EveryBeam`` OSKAR beam code, which uses a logatithmic dipole akin to the SKA dipoles, but on a measurement set with MWA-like stations. This test call both an MWA measurement set and this MWA-like OSKAR measurement set, runs ``EveryBeam`` on both, and checks that the output gains are similar. They should never match exactly as two different dipole models are used, so just check they match within 15%.

test_run_everybeam_over_threads.py
***********************************************************
Checks that the multi-threaded version of calling ``EveryBeam`` gives the same answer as the serial version when running the ``LOFAR LBA`` beam. Don't check actual values here as that's done in ``cmake_testing/GPU_or_C_code/call_everybeam/``


test_run_everybeam_over_threads_MWA.py
***********************************************************
Same as ``test_run_everybeam_over_threads.py``, but using the MWA primary beam, which requires passing in ``parallactic`` angles as ``EveryBeam`` doesn't rotate the MWA primary beam for some reason.


test_run_uvbeam_HERA_CST.py
***********************************************************
Runs the HERA primary beam using ``pyuvdata UVBeam``, when using a list of CST files as input. Runs the primary beam for a grid of coordinates. Checks that the expected zenith gain is close to one. Currently doesn't compare to a set of expected values; might be worth doing that in the future. Checks that the code fails/warns in some expected circumstances: fails with bad file paths; fails if length of list of CST files is not equal to associated frequency array; warns if less that three CST files are provided (as frequency interpolation  must be switched off). There are commented lines within that can be uncommented to make beam plots.

test_run_uvbeam_HERA_FITS.py
***********************************************************
Runs the HERA primary beam using ``pyuvdata UVBeam``, when using a beam FITS file as input. Runs the primary beam for a grid of coordinates. Checks that the expected zenith gain is close to one. Currently doesn't compare to a set of expected values; might be worth doing that in the future. Checks that the code fails if a bad path to the FITS file is give.

test_run_uvbeam_MWA.py
***********************************************************
Runs the MWA primary beam using ``pyuvdata UVBeam``. Runs a grid of directions through both ``hyperbeam`` and ``pyuvdata UVBeam`` and checks the gains are close. Tests for a zenith and off-zenith MWA pointing, and when setting a single dipole as dead. Checks the code fails is the user passes bad amplitude or bad delays. There are commented lines within that can be uncommented to make beam plots.