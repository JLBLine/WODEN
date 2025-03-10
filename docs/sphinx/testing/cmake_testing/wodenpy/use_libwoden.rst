``use_libwoden``
=========================
This is a collection of modules that interact with the C/GPU code via ``ctypes``. Much of the code here defined ``ctypes`` ``structures`` classes that are equivalent to a ``C`` ``struct``. This allows us to pass data back and forth between ``Python`` and ``C``. Much of the functionality in ``wodenpy.use_libwoden`` is tested by other unit tests that use the aforementioned classes. The tests here attempt to fill in any gaps.

test_check_for_everybeam.py
*************************************
``check_for_everybeam`` interrogates ``libwoden_double.so`` to see if it was compiled with the ``HAVE_EVERBEAM`` flag. This lets ``wodenpy`` know if it can call EveryBeam functions or not. Test by inputting two different test dynamic libraries, one with the flag and one without, and checking the output.

test_convert_woden_settings_to_ctypes.py
****************************************
``woden_settings.convert_woden_settings_to_ctypes`` takes a filled ``Woden_Settings_Python`` class, and copies it across into an equivalent ``Woden_Settings_Double`` or ``Woden_Settings_Float`` ``ctypes`` class. Test by creating a ``Woden_Settings_Python`` class, filling it with some values, and then calling the function. We then pass the ``ctypes`` class to a ``C`` function that writes the contents to a text file. We then read the text file back in and compare to the original ``Woden_Settings_Python`` class. Test with various settings to check different combinations of settings are copied across correctly.  

test_create_sbf.py 
**************************************
This tests that the shapelet basis function array ``sbf`` is created correctly. Test by simply by calling the function and sampling at 20 random array indexes, and comparing to stored expected values.

test_fill_woden_settings_python.py
**************************************
Tests ``woden_settings.fill_woden_settings_python``, which populates a ``Woden_Settings_Python`` class. The function takes in arguments directly from ``argparse``, so test by inputting a number of fake argument tests. Whack in many many combinations of arguments.

test_make_skymodel_structs.py
************************************
This is a very basic test. ``make_skymodel_structs`` creates a ``ctypes`` struct that is used to pass a sky model to the C code. We test by creating an empty ``ctypes`` ``Source_Catalogue`` sky model, filling in the ``RA`` values, passing that sky model across to some C code which writes the ``RA`` values to a text file. We then read the text file back in and compare to the original ``RA`` values. This function is used more in other integration tests, but this test should probably be expanded. 


test_setup_lsts_and_phase_centre.py
**************************************
This populates the LST values in a ``woden_settings`` instance. Test by generating a ``woden_settings`` instance, and calling with and without precession set. Compare to a set of expected values.
