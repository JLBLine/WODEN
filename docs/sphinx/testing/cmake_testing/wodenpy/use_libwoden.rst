``use_libwoden``
=========================
This is a collection of modules that interact with the C/CUDA code via ``ctypes``. Much of the code here defined ``ctypes`` ``structures`` classes that are equivalent to a ``C`` ``struct``. This allows us to pass data back and forth between ``Python`` and ``C``. Much of the functionality in ``wodenpy.use_libwoden`` is tested by other unit tests that use the aforementioned classes. The tests here attempt to fill in any gaps.

test_create_sbf.py 
**************************************
This tests that the shapelet basis function array ``sbf`` is created correctly. Test by simply by calling the function and sampling at 20 random array indexes, and comparing to stored expected values.


test_make_woden_settings.py
**************************************
Tests ``woden_settings.create_woden_settings``, which populates a ``Woden_Settings_Double`` or ``Woden_Settings_Float``. The function takes in arguments directly from ``argparse``, so test by inputting a number of fake argument tests. Go the whole hog to verify we understand how ``ctypes`` works by writing ``C`` code that we can pass the populated ``Woden_Setting`` class to. This ``C`` code then writes the contents to a text file, which is read back in by ``Python`` and compared to the original ``Woden_Settings`` class. This is a bit of a pain, but it's the only way to be sure that the ``ctypes`` ``structure`` is being populated correctly.

test_setup_lsts_and_phase_centre.py
**************************************
This populates the LST values in a ``woden_settings`` instance. Test by generating a ``woden_settings`` instance, and calling with and without precession set. Compare to a set of expected values.