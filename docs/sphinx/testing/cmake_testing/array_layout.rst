``array_layout``
=========================
Tests for the functions in ``WODEN/src/array_layout.c``. These functions handle
reading in array coords of *e,n,h* and converting them into *X,Y,Z*, used later
to calculate *u,v,w*. Also handles precession of the array layout and LST from the
simulation date back to J2000, so the array is in the same coord frame as the
sky model. Should be the equivalent to rotating the sky model forward to current
day and ensures the output visibilities are in the J2000 frame.

.. _test_RTS_ENH2XYZ_local.c:

``test_RTS_ENH2XYZ_local.c``
*****************************
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

Both the FLOAT and DOUBLE compiled versions are 64 bit precision. Both are
tested to match a 64 bit ``C`` calculation of the above equations to within an
absolute tolerance of 1e-13.

``test_RTS_PrecessXYZtoJ2000.c``
*********************************
``array_layout::RTS_PrecessXYZtoJ2000`` takes observation date ``X,Y,Z`` coords
and precesses them back to J2000 for all time steps in
``woden_settings->lsts``. Test by giving a known set of *X,Y,Z* coordinates
lsts and mjds, checks output precessed *X,Y,Z* s match the expected values
stored in ``test_RTS_XYZ_common.h``. Test covers two time steps.

``test_calc_XYX_diffs.c``
****************************
Tests the function ``array_layout::calc_XYZ_diffs``, which reads in an array
text file in *e,n,h* coords and transforms to *X,Y,Z*, precesses the locations
back to J2000 if requested, and then calculates the baseline length in *X,Y,Z*.

These tests simply provide a set of known input coords (read in from
``example_array_layout.txt``), runs ``array_layout::calc_XYZ_diffs`` for both
the precession and no precession cases, and check the output values match the
expected values stored in ``test_RTS_XYZ_common.h``.

``calc_XYZ_diffs`` also sets the basic values of::

   woden_settings->num_visis
   woden_settings->num_ants
   woden_settings->num_baselines
   woden_settings->num_cross
   woden_settings->num_autos

which are used at various times during the rest of the simulation. The
value of ``num_autos`` depends on whether the user has elicited to calculate
auto-correlations, so should be 0 if ``woden_settings->do_autos`` if False,
and correctly set if True. The numbers above are tested in both the False
and True cases.
