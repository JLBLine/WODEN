Testing via ``ctest``
======================

This is totally optional, but you can run the unit / integration tests I use for
developing ``WODEN`` to check your system / the code runs as expected. The tests
are compiled using the same ``cmake`` script as the main code.

ALL tests are compiled in both FLOAT and DOUBLE precision, meaning the same
test code is used to test the two different precision versions, with the compiler
dropping in the necessary precision. Unless explicitly noted in the details
in :ref:`What do the tests actually do?`, the tests require the same level
of accuracy from the FLOAT and DOUBLE precision versions.

Dependencies
-------------

The tests use the following C library:

* **Unity** - https://github.com/ThrowTheSwitch/Unity

I installed ``unity`` via::

  $ git clone https://github.com/ThrowTheSwitch/Unity.git
  $ cd Unity
  $ mkdir build && cd build
  $ cmake ..
  $ make -j 4
  $ sudo make install

That way, ``cmake`` can find ``unity``. However, you don't need to install Unity anywhere, as ``WODEN`` uses the ``C`` code directly. You just need to tell ``WODEN`` where Unity lives in your system (for example, you could download a release version e.g. version 2.5.2 - see example below for how to link without installation).

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
         Start  1: C_test_RTS_ENH2XYZ_local_float
    1/87 Test  #1: C_test_RTS_ENH2XYZ_local_float .......................   Passed    0.00 sec
         Start  2: C_test_calc_XYZ_diffs_float
    2/87 Test  #2: C_test_calc_XYZ_diffs_float ..........................   Passed    0.00 sec
         Start  3: C_test_RTS_PrecessXYZtoJ2000_float
    3/87 Test  #3: C_test_RTS_PrecessXYZtoJ2000_float ...................   Passed    0.00 sec
         Start  4: C_test_RTS_ENH2XYZ_local_double
    4/87 Test  #4: C_test_RTS_ENH2XYZ_local_double ......................   Passed    0.00 sec
         Start  5: C_test_calc_XYZ_diffs_double
    5/87 Test  #5: C_test_calc_XYZ_diffs_double .........................   Passed    0.00 sec
         Start  6: C_test_RTS_PrecessXYZtoJ2000_double
    6/87 Test  #6: C_test_RTS_PrecessXYZtoJ2000_double ..................   Passed    0.00 sec
         Start  7: C_test_null_comps_float
    7/87 Test  #7: C_test_null_comps_float ..............................   Passed    0.00 sec
         Start  8: C_test_fill_chunk_src_with_pointgauss_float
    8/87 Test  #8: C_test_fill_chunk_src_with_pointgauss_float ..........   Passed    0.03 sec
         Start  9: C_test_fill_chunk_src_with_shapelets_float
   ...etc etc

.. note:: To test MWA Fully Embedded Element beam code, you must have the environment variable::

    MWA_FEE_HDF5=/path/to/mwa_full_embedded_element_pattern.h5

  declared. If you don't have that, the MWA FEE beam code tests will just be skipped.

If you want more detail of what the tests are doing, run::

  $ ctest --verbose

What do the tests actually do?
---------------------------------

The tests are all located in ``WODEN/cmake_testing``, and each directory within contains tests
for a different file from ``WODEN/src``. Within each test directory, there are separate files for testing different functions, which include the function name. As an example, the directory ``WODEN/cmake_testing/array_layout`` contains tests for the file ``WODEN/src/array_layout.c``, and contains test files that test the following functions::

  cmake_testing/array_layout/test_calc_XYZ_diffs.c -> src/array_layout.c::calc_XYZ_diffs
  cmake_testing/array_layout/test_RTS_ENH2XYZ_local.c -> src/array_layout.c::RTS_ENH2XYZ_local
  cmake_testing/array_layout/test_RTS_PrecessXYZtoJ2000.c -> src/array_layout.c::RTS_PrecessXYZtoJ2000

.. note:: For those unfamiliar with ``CMake`` testing, even though the tests are located in ``WODEN/cmake_testing/``, when you run ``ctest``, the test files are copied and run in ``WODEN/build/cmake_testing``, so any output from the tests will be located there.

The sections below give an outline of the tests performed in each directory.

.. toctree::
   :maxdepth: 1

   cmake_testing/array_layout
   cmake_testing/chunk_sky_model
   cmake_testing/create_sky_model
   cmake_testing/FEE_primary_beam
   cmake_testing/fundamental_coords
   cmake_testing/primary_beam
   cmake_testing/shapelet_basis

``CUDA`` code tests:

.. toctree::
   :maxdepth: 1

   cmake_testing/calculate_visibilities
   cmake_testing/FEE_primary_beam_cuda
   cmake_testing/primary_beam_cuda
   cmake_testing/source_components

.. note:: To be able to test ``CUDA`` functions that are designed to work solely in GPU memory, it's necessary to write wrapper functions that allocate GPU memory, pass the data into the ``CUDA`` code to be tested, and then copy the results back into host memory. I've kept these 'intermediate' test functions inside the ``*.cu`` files that contain the code being tested, as it's not straight forward / performance degrading to have them in separate files. On casual inspection it looks like there are many functions in the ``*.cu`` files I haven't written tests for, but the extra functions are there *because* of testing. Sigh.
