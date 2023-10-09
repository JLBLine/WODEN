``array_layout``
=========================
Tests for the functions in ``wodenpy.array_layout``. These functions handle
reading in array coords of *e,n,h* and converting them into *X,Y,Z*, used later
to calculate *u,v,w*. Also handles precession of the array layout and LST from the
simulation date back to J2000, so the array is in the same coord frame as the
sky model. Should be the equivalent to rotating the sky model forward to current
day and ensures the output visibilities are in the J2000 frame.

.. _test_RTS_ENH2XYZ_local.c:

``test_enh2xyz.py``
*****************************
``array_layout.create_array_layout.enh2xyz`` transforms a local east, north, height coord
into *X,Y,Z* coords, which can be used to calculate *u,v,w*. Tests here
generate some test *e,n,h* coords, and then tests for two latitude test cases:

 - latitude = :math:`0^\circ`, should result in:
    - :math:`X = h`
    - :math:`Y = e`
    - :math:`Z = n`
 - latitude = :math:`-30^\circ`, should result in:
    - :math:`X = \frac{n}{2} + \frac{\sqrt(3)}{2}h`
    - :math:`Y = e`
    - :math:`Z = \frac{\sqrt(3)}{2}n - \frac{h}{2}`

Tested to an absolute tolerance of 1e-15.

``test_RTS_PrecessXYZtoJ2000.py``
*********************************
``array_layout.create_layout.RTS_PrecessXYZtoJ2000`` takes observation date ``X,Y,Z`` coords (via a ``woden_settings`` and ``args`` object)
and precesses them back to J2000 for all time steps in
``woden_settings.lsts``. Test by giving a known set of *X,Y,Z* coordinates
lsts and mjds, checks output precessed *X,Y,Z* s match the expected values
stored in ``test_RTS_PrecessXYZtoJ2000.py``. Test covers two time steps.

``test_calc_XYZ_diffs.py``
****************************
Tests the function ``array_layout.create_array_layout.calc_XYX_diffs``, which takes in
*e,n,h* coords from ``woden_settings`` and ``args`` objects, and transforms to *X,Y,Z*, precesses the locations
back to J2000 if requested, and then calculates the baseline length in *X,Y,Z*.

These tests simply provide a set of known input coords for both
the precession and no precession cases, and check the output values match the
expected values stored in ``test_calc_XYZ_diffs.py``.

``test_RTS_precess.py``
*************************
Tests the function ``array_layout.precession.RTS_Precess_LST_Lat_to_J2000``, which takes in a current LST, latitude, and modified julian date, and precesses the LST and latitude back to J2000. Test by inputting known LSTs, latitudes, and mjds, and checking outputs match.