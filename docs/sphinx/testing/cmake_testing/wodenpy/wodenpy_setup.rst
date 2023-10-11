``wodenpy_setup``
=====================
Contains code to parse arguments from the command line and read in settings from a metafits file.

test_argument_inputs.py
************************
Tests the argparser by creating an ``inputs`` object to mimic command line inputs. This test suite is run on both the ``wodenpy_setup.run_setup.get_parser`` and ``wodenpy_setup.run_setup.check_args`` functions. The latter should check certain inputs and fail out with an error message, so a number of the tests are designed to assert that an error message is returned. An example metafits file is bundled with the tests, and values read in by the code are compared to stored values to ensure the reader.