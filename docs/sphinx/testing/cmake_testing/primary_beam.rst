``primary_beam``
=========================
Tests for the functions in ``WODEN/src/primary_beam.c``. These functions
setup primary beam settings, ready to calculate beam responses on the GPU.

test_fill_primary_beam_settings.c
***********************************
``primary_beam::fill_primary_beam_settings`` prepares a ``beam_settings_t``
struct to be used by ``calculate_visibilities::calculate_visibilities``. The
az,za coords have already been calculated by
``chunk_sky_model::create_chunked_sky_models``, which is sufficient for some
beam models. The following beam models need further inputs:

   - ``GAUSS_BEAM`` and ``MWA_ANALY``: both need the hour angle of all COMPONENTs for all time steps, so ``fill_primary_beam_settings`` makes this calculation

The tests here call ``fill_primary_beam_settings`` for all primary
beam types, and perform the following checks:

 - ``GAUSS_BEAM``:
    - Assert that ``beam_settings->beamtype == GAUSS_BEAM``
    - Assert that a number of constants are copied from ``woden_settings`` into ``beam_settings``

 - ``GAUSS_BEAM`` or ``MWA_ANALY``
    - Assert that::

        src->point_components.beam_decs
        src->point_components.beam_has
        src->gauss_components.beam_decs
        src->gauss_components.beam_has
        src->shape_components.beam_decs
        src->shape_components.beam_has

      have been set to the correct values for all COMPONENTs and time steps.
 - All other beams:
    - Assert that ``beam_settings->beamtype`` is set correctly

For FLOAT compiled code, the absolute tolerance threshold on values is set to
1e-7, and 1e-15 for DOUBLE compiled code. For values that are 64 bit in both the
FLOAT and DOUBLE versions the values are tested using ``TEST_ASSERT_EQUAL_DOUBLE``.
