``primary_beam``
=========================
Tests for the functions in ``WODEN/src/primary_beam.c``. These functions
setup primary beam settings, ready to calculate beam responses on the GPU.

test_calc_para_angle.c
*********************************
``primary_beam::calc_para_angle`` calculates the parallactic angle for
all COMPONENTs, for all time steps. This test calls ``calc_para_angle`` for
three POINT, three GAUSSIAN, and three SHAPELET COMPONENTs, for three different
LSTs.. The different COMPONENT types are given different RA/Decs. The
results are compared to expected values, which are stored in
``expected_para_angles.h``


test_fill_primary_beam_settings.c
***********************************
``primary_beam::fill_primary_beam_settings`` prepares a ``beam_settings_t``
struct to be used by ``calculate_visibilities::calculate_visibilities``. The
az,za coords have already been calculated by
``chunk_sky_model::create_chunked_sky_models``, which is sufficient for some
beam models. There are two beam models that need further inputs:

   - ``GAUSS_BEAM``: the Gaussian beam function uses *l,m,n* coords to      incorporate projection effects, so ``fill_primary_beam_settings`` calculates the hour angle of all COMPONENTs for all time steps, to feed into ``fundamental_coords::kern_calc_lmn`` later down the line
   - ``FEE_BEAM``: needs the parallactic angle to rotate the telescope-based coords into the Stokes frame

The tests here call ``fill_primary_beam_settings`` for the four primary
beam types, and perform the following checks:

 - ``GAUSS_BEAM``:
    - Assert that ``beam_settings->beamtype == GAUSS_BEAM``
    - Assert that a number of constants are copied from ``woden_settings`` into ``beam_settings``
    - Assert that::

        src->point_gaussbeam_decs
        src->point_gaussbeam_has
        src->gauss_gaussbeam_decs
        src->gauss_gaussbeam_has
        src->shape_gaussbeam_decs
        src->shape_gaussbeam_has

      have been set to the correct values for all COMPONENTs and time steps
 - ``FEE_BEAM``:
    - Assert that ``beam_settings->beamtype == FEE_BEAM``
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
