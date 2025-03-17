``wodenpy_setup``
=====================
Contains code to parse arguments from the command line and read in settings from a metafits file.

test_argument_inputs.py
************************
Tests the argparser by creating an ``inputs`` object to mimic command line inputs. This test suite is run on both the ``wodenpy_setup.run_setup.get_parser`` and ``wodenpy_setup.run_setup.check_args`` functions. The latter should check certain inputs and fail out with an error message, so a number of the tests are designed to assert that an error message is returned. An example metafits file is bundled with the tests, and values read in by the code are compared to stored values to ensure the reader. Also checks values from a measurement set are read in, needed by EveryBeam.

test_get_code_version.py
***************************
This is a token test to just check this function runs; the output depends on what version of ``WODEN`` you are using and how you installed it, so very difficult to test with a fixed value.


test_log_chosen_beamtype.py
***************************
This function logs information about what beam type has been selected. Test by running the ``WODEN`` parser with various beam models in the input arguments, and checking the looger writes out the correct information to a log file.

test_make_logger.py
***************************
Test the functionality of the ``WODEN`` logger when running multiple threads and ``C`` code. Launches a couple of threads, calls some ``C`` code that takes input numbers, and checks the logger writes out the correct information to a log file.