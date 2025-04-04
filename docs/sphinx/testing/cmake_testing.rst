Testing via ``ctest``
======================

This is totally optional, but you can run the unit / integration tests I use for
developing ``WODEN`` to check your system / the code runs as expected. The tests
are compiled using the same ``cmake`` script as the main code.

The test code lives in the ``WODEN/cmake_testing`` directory. See below on how to compile. Everything is run through ``ctest``, a testing tool that comes with ``cmake``. The Python code is techincally run via ``pytest``, but the tests are still called through ``ctest``.

ALL tests are compiled in both FLOAT and DOUBLE precision, meaning the same
test code is used to test the two different precision versions, with the compiler
dropping in the necessary precision. See below in :ref:`What do the tests actually do?` for an overview of what the tests actually do.

Dependencies
-------------

The tests use the following C library:

* **Unity** - https://github.com/ThrowTheSwitch/Unity

You don't need to install Unity, just clone it somewhere::

  $ git clone https://github.com/ThrowTheSwitch/Unity.git

You don't need to install Unity anywhere, as ``WODEN`` compiles the ``C`` code directly. You just need to tell ``WODEN`` where Unity lives in your system (for example, you could download a release version e.g. version 2.5.2 - see example below for how to link without installation).

You'll also need to initiate the ``git submodule`` that runs code coverage. Simply navigate to the ``WODEN`` directory and run::

  $ git submodule init
  $ git submodule update

which will pull in the relevant ``CMake-codecov`` dependencies. This allows us to track code coverage for the ``python`` and ``C`` code (no free tools exist for ``CUDA`` at the time of writing, boooo).

To tell ``cmake`` to build tests, you add ``TARGET_GROUP=test`` to your command to tell CMake to build tests instead of the main code::

  $ cd $WODEN_DIR
  $ cmake .. -DTARGET_GROUP=test -DUNITY_ROOT=/usr/local/Unity-2.5.2
  $ make -j 4

(where ``-DUNITY_ROOT=`` is needed if you didn't install ``unity``).

.. warning:: once you have done this, to go back to compiling the main ``woden`` executable, you need to run::

    $ cmake .. -DTARGET_GROUP=production

    otherwise you'll just keep building the tests.

Running tests
--------------

Once that compiles, you can run the tests by running::

  $ ctest

You should see something like the following if successful::

  $ ctest
  Test project /home/jline/software/WODEN/build
        Start  1: C_test_fill_primary_beam_settings_float
   1/85 Test  #1: C_test_fill_primary_beam_settings_float ....................   Passed    0.00 sec
        Start  2: C_test_fill_primary_beam_settings_double
   2/85 Test  #2: C_test_fill_primary_beam_settings_double ...................   Passed    0.00 sec
        Start  3: C_test_fill_timefreq_visibility_set_float
   3/85 Test  #3: C_test_fill_timefreq_visibility_set_float ..................   Passed    0.00 sec
        Start  4: C_test_malloc_and_free_float
   4/85 Test  #4: C_test_malloc_and_free_float ...............................   Passed    0.00 sec
        Start  5: C_test_write_visi_set_binary_float
   5/85 Test  #5: C_test_write_visi_set_binary_float .........................   Passed    0.00 sec
        Start  6: C_test_write_visi_set_text_float
   6/85 Test  #6: C_test_write_visi_set_text_float ...........................   Passed    0.00 sec
        Start  7: C_test_fill_timefreq_visibility_set_double
   7/85 Test  #7: C_test_fill_timefreq_visibility_set_double .................   Passed    0.00 sec
        Start  8: C_test_malloc_and_free_double
   8/85 Test  #8: C_test_malloc_and_free_double ..............................   Passed    0.00 sec
   ...etc etc

.. note:: To test MWA Fully Embedded Element beam code, you must have the environment variable::

    MWA_FEE_HDF5=/path/to/mwa_full_embedded_element_pattern.h5

  declared. If you don't have that, the MWA FEE beam code tests will just be skipped.

If you want more detail of what the tests are doing, run::

  $ ctest --verbose

What do the tests actually do?
---------------------------------

The tests are all located in ``WODEN/cmake_testing``, and each directory therein contains tests
for a different file from ``WODEN/wodenpy``, ``WODEN/src``, and ``WODEN/scripts``. Within each test directory, there are separate directories for different files, and in each directory, a different test file for testing a specific function, which include the function name. As an example, the directory ``WODEN/cmake_testing/GPU_or_C_code/beam_settings`` contains tests for the file ``WODEN/src/beam_settings.c``, and contains test files that test the following functions::

  test_fill_primary_beam_settings.c -> src/beam_settings.c::fill_primary_beam_settings


The ``C`` and ``GPU`` functions are tested using the `Unity`_ library, which has useful functions like::

  TEST_ASSERT_FLOAT_EQUAL();
  TEST_ASSERT_DOUBLE_WITHIN();
  TEST_ASSERT_NULL();

allowing a simple testing of values. If a test says outputs are tested to be
equal, it refers to the ``TEST_ASSERT_FLOAT_EQUAL`` or ``TEST_ASSERT_DOUBLE_EQUAL``
function.

.. _`Unity`: https://github.com/ThrowTheSwitch/Unity

.. note:: For those unfamiliar with ``CMake`` testing, even though the tests are located in ``WODEN/cmake_testing/``, when you run ``ctest``, the test files are copied and run in ``WODEN/build/cmake_testing``, so any output from the tests will be located there.

The sections below give an outline of the tests performed in each directory.

``wodenpy`` tests
~~~~~~~~~~~~~~~~~~~~
These tests are for the ``python`` module in ``WODEN/wodenpy``. They test the ``python`` code using the ``pytest`` library, and are run using the ``ctest`` command. The tests are located in the ``WODEN/cmake_testing/wodenpy`` directory.

.. toctree::
   :maxdepth: 1

   cmake_testing/wodenpy/array_layout
   cmake_testing/wodenpy/observational
   cmake_testing/wodenpy/phase_rotate
   cmake_testing/wodenpy/primary_beam
   cmake_testing/wodenpy/skymodel
   cmake_testing/wodenpy/use_libwoden
   cmake_testing/wodenpy/uvfits
   cmake_testing/wodenpy/wodenpy_setup


``GPU_or_C_code`` tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Tests in some of these sections test code that has both ``C`` and ``GPU`` versions; by using the same test code for both, we can ensure the ``C`` and ``GPU`` versions are consistent. The tests are located in the ``WODEN/cmake_testing/GPU_or_C_code`` directory.

.. toctree::
   :maxdepth: 1

   cmake_testing/beam_settings
   cmake_testing/calculate_visibilities
   cmake_testing/call_everybeam
   cmake_testing/fundamental_coords
   cmake_testing/logger
   cmake_testing/primary_beam
   cmake_testing/run_woden
   cmake_testing/source_components
   cmake_testing/visibility_set

``script`` code:
~~~~~~~~~~~~~~~~~~~~~~~~~
These integration tests check the major ``python`` scripts as working as expected. The tests are located in the ``WODEN/cmake_testing/scripts`` directory.

.. toctree::
   :maxdepth: 1

   cmake_testing/scripts/add_instrumental_effects_woden
   cmake_testing/scripts/add_woden_uvfits
   cmake_testing/scripts/concat_woden_uvfits
   cmake_testing/scripts/run_woden
   cmake_testing/scripts/woden_uv2ms
