Testing
----------

There are two ways to test WODEN. The first is to use the CMake testing suite which is more for developing ``WODEN``, and contains unit/integration tests for ``Python``, ``C``, and ``GPU`` code. The second is via the ``test_installation`` directory which contains a number of commands that can be used to test the functionality of the local installation of ``WODEN``. See below for details.

.. toctree::
   :maxdepth: 2

   cmake_testing
   script_testing

EveryBeam Testing
-----------------------

For testing ``everybeam`` functionality , I've made a number of notebooks that live in ``WODEN/test_installation/everybeam/``. This is on top of normal unit testing, to do proper sanity checks and for full disclosure of what is happening.  The integration tests in the notebooks are:

.. toctree::
   :maxdepth: 2

   everybeam_testing


``pyuvdata UVBeam`` Testing
-----------------------

For testing ``UVBeam`` functionality , I've made a number of notebooks that live in ``WODEN/test_installation/pyuvbeam/``. This is on top of normal unit testing, to do proper sanity checks and for full disclosure of what is happening.  The integration tests in the notebooks are:

.. toctree::
   :maxdepth: 2

   uvbeam_testing.rst