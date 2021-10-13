``array_layout``
=========================
Tests for the functions in ``WODEN/src/array_layout.c``. These functions handle
reading in array coords of *e,n,h* and converting them into *X,Y,Z*, used later
to calculate *u,v,w*. Also handles precession of the array layout and LST from the
simulation date back to J2000, so the array is in the same coord frame as the
sky model. Should be the equivalent to rotating the sky model forward to current
day.

``test_RTS_ENH2XYZ_local.c``
****************************
``array_layout::RTS_ENH2XYZ_local`` transforms a local east, north, height coord
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

``test_RTS_PrecessXYZtoJ2000.c``
*********************************
``array_layout::RTS_PrecessXYZtoJ2000`` takes current date ``X,Y,Z`` coords
and precesses them back to J2000, and updates the LST accordingly. This test
reads. Test by giving a known set of *X,Y,Z* and julian date, and
checks outputs match the expected values stored in ``test_RTS_XYZ_common.h``.

``test_calc_XYX_diffs.c``
****************************
``array_layout::calc_XYZ_diffs`` reads in an array text file in *e,n,h* coords
and transforms to *X,Y,Z*, precesses the locations back to J2000 if requested,
and then calculates the baseline length in *X,Y,Z*. These tests simply provide
a set of known input coords (in ``example_array_layout.txt``), run
``array_layout::calc_XYZ_diffs`` for both the precession and no precssion cases,
and check the output values match the expected values stored in ``test_RTS_XYZ_common.h``.
