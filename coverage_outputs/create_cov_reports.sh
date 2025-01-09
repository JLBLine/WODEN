##Delete any old reports hanging around
rm coverage.xml *.gcov *.gcno *.gcda coverage.info #.coverage

##Use the ctest C code outputs and run gcov over them. Converts them into
##something that codecov can read

cov_src_dir=../build/cmake_testing/GPU_or_C_code/CMakeFiles/calculate_visibilities_CPU_float.dir/__/__/src

for file in ${cov_src_dir}/*.gcda; do
    fileroot=$(basename $file .gcda)
    gcov ${cov_src_dir}/${fileroot}.gcda ${cov_src_dir}/${fileroot}.gcno
done

# ##Run the python tests using python 'coverage'
# ##Can be grabbed with 'pip install coverage'

##array_layout
coverage run --source=wodenpy.array_layout.create_array_layout ../cmake_testing/wodenpy/array_layout/test_calc_XYZ_diffs.py
coverage run --source=wodenpy.array_layout.create_array_layout ../cmake_testing/wodenpy/array_layout/test_enh2xyz.py
coverage run --source=wodenpy.array_layout.create_array_layout ../cmake_testing/wodenpy/array_layout/test_RTS_PrecessXYZtoJ2000.py
coverage run --source=wodenpy.use_libwoden.array_layout_struct ../cmake_testing/wodenpy/array_layout/test_RTS_PrecessXYZtoJ2000.py
coverage run --source=wodenpy.array_layout.precession ../cmake_testing/wodenpy/array_layout/test_RTS_precess.py

##observational
coverage run --source=wodenpy.observational.calc_obs ../cmake_testing/wodenpy/observational/test_calc_jdcal.py
coverage run --source=wodenpy.observational.calc_obs ../cmake_testing/wodenpy/observational/test_get_uvfits_date_and_position_constants.py

##phase_rotate
coverage run --source=wodenpy.phase_rotate.remove_phase_track ../cmake_testing/wodenpy/phase_rotate/test_remove_phase_tracking.py

##primary_beam
coverage run --source=wodenpy.primary_beam.use_everybeam \
         ../cmake_testing/wodenpy/skymodel/test_calc_everybeam_for_components.py
coverage run --source=wodenpy.primary_beam.use_everybeam \
         ../cmake_testing/wodenpy/primary_beam/test_run_everybeam_over_threads.py
coverage run --source=wodenpy.primary_beam.use_everybeam \
         ../cmake_testing/wodenpy/primary_beam/test_run_woden_lofar_para_rotate.py
coverage run --source=wodenpy.primary_beam.use_everybeam \
         ../cmake_testing/wodenpy/primary_beam/test_run_fast_everybeam.py

##skymodel
coverage run --source=wodenpy.skymodel.woden_skymodel ../cmake_testing/wodenpy/skymodel/test_crop_below_horizon.py
coverage run --source=wodenpy.skymodel.chunk_sky_model ../cmake_testing/wodenpy/skymodel/test_create_skymodel_chunk_map.py
coverage run --source=wodenpy.skymodel.chunk_sky_model ../cmake_testing/wodenpy/skymodel/test_map_chunk_pointgauss.py
coverage run --source=wodenpy.skymodel.chunk_sky_model ../cmake_testing/wodenpy/skymodel/test_map_chunk_shapelets.py
coverage run --source=wodenpy.skymodel.read_fits_skymodel ../cmake_testing/wodenpy/skymodel/test_read_FITS_skymodel_chunk.py
coverage run --source=wodenpy.skymodel.read_fits_skymodel ../cmake_testing/wodenpy/skymodel/test_check_columns_fits.py
coverage run --source=wodenpy.use_libwoden.skymodel_structs ../cmake_testing/wodenpy/skymodel/test_read_FITS_skymodel_chunk.py
coverage run --source=wodenpy.skymodel.read_text_skymodel ../cmake_testing/wodenpy/skymodel/test_read_text_skymodel_chunk.py
coverage run --source=wodenpy.skymodel.woden_skymodel ../cmake_testing/wodenpy/skymodel/test_read_text_skymodel_chunk.py
coverage run --source=wodenpy.skymodel.read_yaml_skymodel ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
coverage run --source=wodenpy.skymodel.read_skymodel ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
coverage run --source=wodenpy.skymodel.woden_skymodel ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
coverage run --source=wodenpy.use_libwoden.beam_settings ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
coverage run --source=wodenpy.skymodel.read_fits_skymodel ../cmake_testing/wodenpy/skymodel/test_calc_everybeam_for_components.py

##use_libwoden
export CMAKE_CURRENT_SOURCE_DIR=../cmake_testing/wodenpy/use_libwoden/
coverage run --source=wodenpy.use_libwoden.shapelets ../cmake_testing/wodenpy/use_libwoden/test_create_sbf.py
coverage run --source=wodenpy.wodenpy_setup.run_setup ../cmake_testing/wodenpy/use_libwoden/test_make_woden_settings.py
coverage run --source=wodenpy.use_libwoden.woden_settings ../cmake_testing/wodenpy/use_libwoden/test_make_woden_settings.py
coverage run --source=wodenpy.wodenpy_setup.run_setup ../cmake_testing/wodenpy/use_libwoden/test_setup_lsts_and_phase_centre.py
coverage run --source=wodenpy.use_libwoden.woden_settings ../cmake_testing/wodenpy/use_libwoden/test_setup_lsts_and_phase_centre.py

##uvfits
coverage run --source=wodenpy.uvfits.wodenpy_uvfits ../cmake_testing/wodenpy/uvfits/test_RTS_encoding.py
coverage run --source=wodenpy.uvfits.wodenpy_uvfits ../cmake_testing/wodenpy/uvfits/test_create_uvfits.py
coverage run --source=wodenpy.uvfits.wodenpy_uvfits ../cmake_testing/wodenpy/uvfits/test_make_antenna_table.py
coverage run --source=wodenpy.uvfits.wodenpy_uvfits ../cmake_testing/wodenpy/uvfits/test_make_baseline_date_arrays.py
coverage run --source=wodenpy.uvfits.wodenpy_uvfits ../cmake_testing/wodenpy/uvfits/test_read_uvfits_into_pyuvdata.py

##run_setup
coverage run --source=wodenpy.wodenpy_setup.run_setup ../cmake_testing/wodenpy/wodenpy_setup/test_argument_inputs.py
coverage run --source=wodenpy.wodenpy_setup.run_setup ../cmake_testing/wodenpy/wodenpy_setup/test_get_code_version.py

##scripts
coverage run --source=add_woden_uvfits ../cmake_testing/scripts/add_woden_uvfits/test_add_woden_uvfits.py
coverage run --source=concat_woden_uvfits ../cmake_testing/scripts/concat_woden_uvfits/test_concat_woden_uvfits.py
coverage run --source=run_woden ../cmake_testing/scripts/run_woden/test_run_woden.py
coverage run --source=wodenpy.use_libwoden.visibility_set ../cmake_testing/scripts/run_woden/test_run_woden.py
coverage run --source=wodenpy.use_libwoden.use_libwoden ../cmake_testing/scripts/run_woden/test_run_woden.py
coverage run --source=woden_uv2ms ../cmake_testing/scripts/woden_uv2ms/test_woden_uv2ms.py
coverage run --source=add_instrumental_effects_woden ../cmake_testing/scripts/add_instrumental_effects_woden/test_add_instrumental_effects_woden.py

#convert output to something that codecov accepts
coverage combine
coverage xml

##delete things that were written out by running tests
rm WODEN_array_layout.txt *.uvfits *.json \
    example.txt test_load_data.dat test_full_skymodel* \
    woden_settings.txt *.so *.png *.npz 


##Use this to create two local reports, one for C and one for python
##Only uncomment if you want to see the reports locally without pushing to codecov

lcov --capture --directory ${cov_src_dir} --output-file coverage.info
genhtml coverage.info --output-directory C_coverage
coverage html -d python_coverage

