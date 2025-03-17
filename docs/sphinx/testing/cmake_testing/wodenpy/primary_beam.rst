``primary_beam``
=========================
Tests for the functions in ``wodenpy/primary_beam``, which are wrappers around ``C++`` to run the EveryBeam primary beam code.

test_check_ms_telescope_type_matches_element_response.py
***********************************************************
This function checks that the telescope type in a measurement set allows for the requested element response. For example, the only type of primary beam model that can be used with the MWA telescope is the ``MWA`` beam model. Checks by reading in three measurement sets with ``MWA``, ``LOFAR``, and ``OSKAR`` telescope types. Asserts function fails if an incorrect element response is requested, and passes if an appropriate element response is requested. Also writes out an unknown telescope type measurement set to check that fails as well.


test_run_everybeam_over_threads.py
***********************************************************
Checks that the multi-threaded version of calling ``EveryBeam`` gives the same answer as the serial version when running the ``LOFAR LBA`` beam. Don't check actual values here as that's done in ``cmake_testing/GPU_or_C_code/call_everybeam/``


test_run_everybeam_over_threads_MWA.py
***********************************************************
Same as ``test_run_everybeam_over_threads.py``, but using the MWA primary beam, which requires passing in ``parallactic`` angles as ``EveryBeam`` doesn't rotate the MWA primary beam for some reason.