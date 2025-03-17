``call_everybeam``
=========================
Tests for the functions in ``WODEN/src/call_everybeam.c``. Obviously, these
tests are only pertinent if you have compiled against the EveryBeam library, 
and I've setup CMake to only compile these tests if you have. 

To create values to test against, I've written a number of Python notebooks
that call the native EveryBeam python wrapper. These notebooks dump out
values to ``.h`` header files, to directly include in the tests. If you ever
need to regenerate these values, you'll need to install the EveryBeam Python
wrapper, which means installing with ``cmake .. -DBUILD_WITH_PYTHON=ON``.



``test_check_ms_telescope_type.c``
********************************************
``check_ms_telescope_type`` reads the telescope type from a Measurement Set.
Test by reading in three different MS, and checking the correct telescope type
is returned. The MS are from the MWA, LOFAR, and OSKAR.

``test_load_everybeam_telescope.c``
********************************************
``load_everybeam_telescope`` reads in the EveryBeam telescope object from a 
measurement set. Check we can read in an MWA, LOFAR, and OSKAR MS with no error.


``test_run_*_beam.c``
********************************************
All tests ``test_run_*_beam.c`` work the same, barring ``test_run_mwa_beam.c``,
which is explained below. These tests run the relevant EveryBeam primary beam
model for a given telescope (LBA and HBA of LOFAR and an SKA model of OSKAR).
A 31 by 31 grid of ``ra,dec`` directions are given to the LOFAR tests, and a 21
by 21 grid to the OSKAR set (OSKAR beam model is much slower). All tests run the
beam at two frequencies, two time steps, for two stations. The tests are run;
 - one with parallactic rotation and without normalisation
 - one with parallactic rotation and normalisation
 - one reordering into the IAU polarisation convention

All values are testing against stored values in the header files.

``test_run_mwa_beam.c``
*********************************
Similar to other ``test_run_*_beam.c`` tests, however the MWA model takes in
``az,za`` as directions. Only one station is run as I have no way of controlling
dipole amplitudes so all ``MWA`` stations are identical. The MWA model
also seems to always come out normalised, so only test the parallactic rotation
and reordering into the IAU polarisation convention.