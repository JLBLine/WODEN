# ## Create tests for both float and double precision
# # foreach(PRECISION IN LISTS FLOAT DOUBLE)
foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)
  ##----------------------------------------------------------------------------
  ##Begin adding tests
  ##----------------------------------------------------------------------------

  add_executable(test_load_everybeam_telescope_${PRECISION}_app
      test_load_everybeam_telescope.c
  )
  target_link_libraries(test_load_everybeam_telescope_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  target_compile_options(test_load_everybeam_telescope_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(CXX_test_load_everybeam_telescope_${PRECISION}
           test_load_everybeam_telescope_${PRECISION}_app)

  add_executable(test_load_everybeam_MWABeam_${PRECISION}_app
      test_load_everybeam_MWABeam.c
  )
  target_link_libraries(test_load_everybeam_MWABeam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  target_compile_options(test_load_everybeam_MWABeam_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(CXX_test_load_everybeam_MWABeam_${PRECISION}
           test_load_everybeam_MWABeam_${PRECISION}_app)


  add_executable(test_check_ms_telescope_type_${PRECISION}_app
      test_check_ms_telescope_type.c
  )
  target_link_libraries(test_check_ms_telescope_type_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  target_compile_options(test_check_ms_telescope_type_${PRECISION}_app
                                                              PRIVATE ${C_FLAGS})
  add_test(CXX_test_check_ms_telescope_type_${PRECISION}
          test_check_ms_telescope_type_${PRECISION}_app)

  
  add_executable(test_run_hba_beam_${PRECISION}_app
      test_run_hba_beam.c
  )
  target_link_libraries(test_run_hba_beam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  target_compile_options(test_run_hba_beam_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(CXX_test_run_hba_beam_${PRECISION}
           test_run_hba_beam_${PRECISION}_app)


  add_executable(test_run_lba_beam_${PRECISION}_app
      test_run_lba_beam.c
  )
  target_link_libraries(test_run_lba_beam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  target_compile_options(test_run_lba_beam_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(CXX_test_run_lba_beam_${PRECISION}
           test_run_lba_beam_${PRECISION}_app)


  add_executable(test_run_mwa_beam_${PRECISION}_app
      test_run_mwa_beam.c
  )
  target_link_libraries(test_run_mwa_beam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  target_compile_options(test_run_mwa_beam_${PRECISION}_app
                                                          PRIVATE ${C_FLAGS})
  add_test(CXX_test_run_mwa_beam_${PRECISION}
          test_run_mwa_beam_${PRECISION}_app)

  add_executable(test_run_oskar_beam_${PRECISION}_app
                 test_run_oskar_beam.c
  )
  target_link_libraries(test_run_oskar_beam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  target_compile_options(test_run_oskar_beam_${PRECISION}_app
                                                          PRIVATE ${C_FLAGS})
  add_test(CXX_test_run_oskar_beam_${PRECISION}
          test_run_oskar_beam_${PRECISION}_app)

endforeach()
