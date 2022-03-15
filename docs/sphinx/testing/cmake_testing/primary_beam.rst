``primary_beam``
=========================
Tests for the functions in ``WODEN/src/primary_beam.c``. These functions
setup primary beam settings, ready to calculate beam responses on the GPU.

test_calc_para_angle.c
*********************************
``primary_beam::calc_para_angle`` calculates the parallactic angle for
all COMPONENTs, for all time steps. This test calls ``calc_para_angle`` for
three POINT, three GAUSSIAN, and three SHAPELET COMPONENTs, for three different
LSTs. The different COMPONENT types are given different RA/Decs. The
results are compared to expected values, which are stored in
``expected_para_angles.h``. For FLOAT compiled code, the absolute tolerance
threshold is set to 1e-7, and 1e-15 for DOUBLE compiled code.


test_fill_primary_beam_settings.c
***********************************
``primary_beam::fill_primary_beam_settings`` prepares a ``beam_settings_t``
struct to be used by ``calculate_visibilities::calculate_visibilities``. The
az,za coords have already been calculated by
``chunk_sky_model::create_chunked_sky_models``, which is sufficient for some
beam models. There are a number of beam models that need further inputs:

   - ``GAUSS_BEAM`` and ``MWA_ANALY``: both need the hour angle of all COMPONENTs for all time steps, so ``fill_primary_beam_settings`` makes this calculation
   - ``FEE_BEAM`` and ``FEE_BEAM_INTERP``: need the parallactic angle to rotate the telescope-based coords into the Stokes frame

The tests here call ``fill_primary_beam_settings`` for all primary
beam types, and perform the following checks:

 - ``GAUSS_BEAM``:
    - Assert that ``beam_settings->beamtype == GAUSS_BEAM``
    - Assert that a number of constants are copied from ``woden_settings`` into ``beam_settings``

 - ``GAUSS_BEAM`` or ``MWA_ANALY``
    - Assert that::

        src->point_gaussbeam_decs
        src->point_gaussbeam_has
        src->gauss_gaussbeam_decs
        src->gauss_gaussbeam_has
        src->shape_gaussbeam_decs
        src->shape_gaussbeam_has

      have been set to the correct values for all COMPONENTs and time steps.
   .. note::
     Even though the arrays say 'gaussbeam', really the HA/Dec are what's important, as the analytic MWA beam also needs them. In the next release, the ``catsource_t`` and ``source_catalogue_t`` structs will be overhauled with something far more sensible.
 - ``FEE_BEAM`` or ``FEE_BEAM_INTERP``:
    - Assert that::

        src->sin_point_para_angs
        src->cos_point_para_angs
        src->sin_gauss_para_angs
        src->cos_gauss_para_angs
        src->sin_shape_para_angs
        src->cos_shape_para_angs

      have been set to the correct values for all COMPONENTs and time steps
 - ``ANALY_DIPOLE``:
    - Assert that ``beam_settings->beamtype == ANALY_DIPOLE``
 - ``NO_BEAM``:
    - Assert that ``beam_settings->beamtype == NO_BEAM``

For FLOAT compiled code, the absolute tolerance threshold on values is set to
1e-7, and 1e-15 for DOUBLE compiled code. For values that are 64 bit in both the
FLOAT and DOUBLE versions the values are tested using ``TEST_ASSERT_EQUAL_DOUBLE``.
