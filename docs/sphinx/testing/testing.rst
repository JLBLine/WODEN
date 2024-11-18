Testing
----------

There are two ways to test WODEN. The first is to use the CMake testing suite which is more for developing ``WODEN``, and contains unit/integration tests for ``Python``, ``C``, and ``GPU`` code. The second is via the ``test_installation`` directory which contains a number of commands that can be used to test the functionality of the local installation of ``WODEN``. See below for details.

.. toctree::
   :maxdepth: 2

   cmake_testing
   script_testing

EveryBeam Testing
-----------------------

For testing ``everybeam`` functionality , I've made a number of notebooks that live in ``WODEN/test_installation/everybeam/``. This is in lieu of full unit testing to speed up development. At some point the correct thing would be to go back and break down the testing into units. However, the notebooks are fairly rigorous, so suffice for now. The integration tests in the notebooks are:

.. toctree::
   :maxdepth: 2

   everybeam_testing