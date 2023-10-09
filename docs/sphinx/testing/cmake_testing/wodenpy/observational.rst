``observational``
=========================
Tests for the functions in ``wodenpy.observational``. These functions handle calculating observational parameters given input dates.

test_calc_jdcal.py
*******************************************************
Tests the ``observational.calc_jdcal`` function, which should calculate the Julian Date and
split into a integer day and fractional day values. Just test by calling
``observational.calc_jdcal`` with two known date strings, and checking the output values
match expectations.

test_get_uvfits_date_and_position_constants.py
*******************************************************
Tests the ``observational.get_uvfits_date_and_position_constants`` function,
which should calculate the LST, GST0 (greenwich sidereal time at 0 hours
of the given date), DEGPDY (rotational speed of the Earth) and UT1UTC (
difference between UT1 and UTC) for a given Long/Lat/Height array location and
UTC date. Test with two combinations of different UTC/Long/Lat/Height and
check the returned values are as expected.

