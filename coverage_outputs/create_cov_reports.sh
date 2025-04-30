##Delete any old reports hanging around
rm coverage.xml *.gcov *.gcno *.gcda coverage.info .coverage.*

export COVERAGE_PROCESS_START=.coveragerc

##Use the ctest C code outputs and run gcov over them. Converts them into
##something that codecov can read

cov_src_dir=../build/cmake_testing/GPU_or_C_code/CMakeFiles/calculate_visibilities_CPU_double.dir/__/__/src

for file in ${cov_src_dir}/*.gcda; do
    fileroot=$(basename $file .gcda)
    gcov ${cov_src_dir}/${fileroot}.gcda ${cov_src_dir}/${fileroot}.gcno
done

eb_src_dir=../build/cmake_testing/GPU_or_C_code/CMakeFiles/use_everybeam.dir/__/__/src

for file in ${eb_src_dir}/call_everybeam*.gcda; do
    fileroot=$(basename $file .gcda)
    gcov ${eb_src_dir}/${fileroot}.gcda ${cov_src_dir}/${fileroot}.gcno
done

do_python=True

if [ "$do_python" = "True" ]; then
    

    # ##Run the python tests using python 'coverage'
    # ##Can be grabbed with 'pip install coverage'
    ##array_layout
    coverage run --source=wodenpy ../cmake_testing/wodenpy/array_layout/test_calc_XYZ_diffs.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/array_layout/test_enh2xyz.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/array_layout/test_RTS_PrecessXYZtoJ2000.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/array_layout/test_RTS_precess.py

    ##observational
    coverage run --source=wodenpy ../cmake_testing/wodenpy/observational/test_calc_jdcal.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/observational/test_get_uvfits_date_and_position_constants.py

    ##phase_rotate
    coverage run --source=wodenpy ../cmake_testing/wodenpy/phase_rotate/test_remove_phase_tracking.py

    ##primary_beam
    coverage run --source=wodenpy ../cmake_testing/wodenpy/primary_beam/test_run_everybeam_over_threads.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/primary_beam/test_run_everybeam_over_threads_MWA.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/primary_beam/test_check_ms_telescope_type_matches_element_response.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/primary_beam/test_calc_uvbeam_for_components.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/primary_beam/test_run_uvbeam_MWA.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/primary_beam/test_create_filtered_ms.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/primary_beam/test_run_everybeam_OSKAR.py

    ##skymodel
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_crop_below_horizon.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_create_skymodel_chunk_map.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_map_chunk_pointgauss.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_map_chunk_shapelets.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_read_FITS_skymodel_chunk.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_check_columns_fits.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_read_text_skymodel_chunk.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_read_skymodel_chunk.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/skymodel/test_calc_everybeam_for_components.py

    ##use_libwoden
    export CMAKE_CURRENT_SOURCE_DIR=../cmake_testing/wodenpy/use_libwoden/
    coverage run --source=wodenpy ../cmake_testing/wodenpy/use_libwoden/test_create_sbf.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/use_libwoden/test_fill_woden_settings_python.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/use_libwoden/test_convert_woden_settings_to_ctypes.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/use_libwoden/test_setup_lsts_and_phase_centre.py

    ##uvfits
    coverage run --source=wodenpy ../cmake_testing/wodenpy/uvfits/test_RTS_encoding.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/uvfits/test_create_uvfits.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/uvfits/test_make_antenna_table.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/uvfits/test_make_baseline_date_arrays.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/uvfits/test_read_uvfits_into_pyuvdata.py

    ##run_setup
    coverage run --source=wodenpy ../cmake_testing/wodenpy/wodenpy_setup/test_argument_inputs.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/wodenpy_setup/test_get_code_version.py
    # coverage run --source=wodenpy ../cmake_testing/wodenpy/wodenpy_setup/test_make_logger.py
    coverage run --source=wodenpy ../cmake_testing/wodenpy/wodenpy_setup/test_log_chosen_beamtype.py

    #scripts
    coverage run --source=add_woden_uvfits ../cmake_testing/scripts/add_woden_uvfits/test_add_woden_uvfits.py
    coverage run --source=concat_woden_uvfits ../cmake_testing/scripts/concat_woden_uvfits/test_concat_woden_uvfits.py
    coverage run --source=wodenpy ../cmake_testing/scripts/run_woden/test_run_woden.py
    coverage run --source=run_woden ../cmake_testing/scripts/run_woden/test_run_woden.py
    coverage run --source=run_woden ../cmake_testing/scripts/run_woden/test_run_woden.py Test.test_runs_with_profiler_on
    coverage run --source=run_woden ../cmake_testing/scripts/run_woden/test_run_woden.py Test.test_runs_with_uvbeam


    coverage run --source=wodenpy ../cmake_testing/scripts/run_woden/test_run_woden.py
    coverage run --source=woden_uv2ms ../cmake_testing/scripts/woden_uv2ms/test_woden_uv2ms.py
    coverage run --source=add_instrumental_effects_woden ../cmake_testing/scripts/add_instrumental_effects_woden/test_add_instrumental_effects_woden.py

    #convert output to something that codecov accepts
    coverage combine #--keep
    coverage xml
fi



##delete things that were written out by running tests
rm -r WODEN_array_layout.txt *.uvfits *.json \
    example.txt test_load_data.dat test_full_skymodel* \
    woden_settings.txt *.so *.png *.npz  *.lprof *.log *.ms

# ##Use this to create two local reports, one for C and one for python
# ##Only uncomment if you want to see the reports locally without pushing to codecov

lcov --gcov-tool gcov-12 --capture --directory ${cov_src_dir} --output-file coverage.info
lcov --gcov-tool gcov-12  --capture --directory ${eb_src_dir} --output-file eb_coverage.info
lcov --extract eb_coverage.info 'src/*'  --output-file filtered_coverage.info
lcov --add-tracefile coverage.info --add-tracefile filtered_coverage.info --output-file merged_coverage.info
mv merged_coverage.info coverage.info
rm eb_coverage.info filtered_coverage.info

genhtml coverage.info --output-directory C_coverage

if [ "$do_python" = "True" ]; then
    coverage html -d python_coverage
fi

