``shapelet_basis``
=========================
Tests for the functions in ``WODEN/src/shapelet_basis.c``.

``test_create_sbf.c``
****************************
``shapelet_basis::create_sbf`` just creates a massive look-up table of shapelet
basis function values. Here we call the function and test 20 different
array indexes and assert that they are equal to expected values. Exciting stuff.

For FLOAT compiled code, the absolute tolerance threshold on values is set to
1e-7, and 1e-15 for DOUBLE compiled code.
