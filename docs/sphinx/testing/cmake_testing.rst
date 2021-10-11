Testing via ``ctest``
======================

This is totally optional, but you can run the unit / integration tests I use for developing ``WODEN`` to check your system / the code runs as expected. The tests are compiled using the same ``cmake`` script as the main code.

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
        Start  1: C_test_RTS_ENH2XYZ_local
   1/49 Test  #1: C_test_RTS_ENH2XYZ_local ......................   Passed    0.00 sec
        Start  2: C_test_calc_XYZ_diffs
   2/49 Test  #2: C_test_calc_XYZ_diffs .........................   Passed    0.00 sec
        Start  3: C_test_RTS_PrecessXYZtoJ2000
   3/49 Test  #3: C_test_RTS_PrecessXYZtoJ2000 ..................   Passed    0.00 sec
        Start  4: C_test_null_comps
   4/49 Test  #4: C_test_null_comps .............................   Passed    0.00 sec
        Start  5: C_test_fill_chunk_src_with_pointgauss
   5/49 Test  #5: C_test_fill_chunk_src_with_pointgauss .........   Passed    0.03 sec
        Start  6: C_test_fill_chunk_src_with_shapelets
   6/49 Test  #6: C_test_fill_chunk_src_with_shapelets ..........   Passed    0.01 sec
        Start  7: C_test_create_chunked_sky_models
   7/49 Test  #7: C_test_create_chunked_sky_models ..............   Passed    0.27 sec
        Start  8: C_test_read_source_catalogue
   8/49 Test  #8: C_test_read_source_catalogue ..................   Passed    0.00 sec
   ...etc etc

.. note:: To test MWA Fully Embedded Element beam code, you must have the environment variable::

    MWA_FEE_HDF5=/path/to/mwa_full_embedded_element_pattern.h5

  declared. If you don't have that, the MWA FEE beam code tests will just be skipped.

If you want more detail of what the tests are doing, run::

  $ ctest --verbose
