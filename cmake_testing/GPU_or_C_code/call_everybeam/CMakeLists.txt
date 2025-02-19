# ## Create tests for both float and double precision
# # foreach(PRECISION IN LISTS FLOAT DOUBLE)
foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)
  ##----------------------------------------------------------------------------
  ##Begin adding tests
  ##----------------------------------------------------------------------------

  ##Test primary beam filling function works
  add_executable(test_load_everybeam_telescope_${PRECISION}_app
      test_load_everybeam_telescope.c
  )
  target_link_libraries(test_load_everybeam_telescope_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_load_everybeam_telescope_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(CXX_test_load_everybeam_telescope_${PRECISION}
           test_load_everybeam_telescope_${PRECISION}_app)

  ##Test primary beam filling function works
  add_executable(test_run_hba_beam_${PRECISION}_app
      test_run_hba_beam.c
  )
  target_link_libraries(test_run_hba_beam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_run_hba_beam_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(CXX_test_run_hba_beam_${PRECISION}
           test_run_hba_beam_${PRECISION}_app)


  ##Test primary beam filling function works
  add_executable(test_run_mwa_beam_${PRECISION}_app
      test_run_mwa_beam.c
  )
  target_link_libraries(test_run_mwa_beam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_run_mwa_beam_${PRECISION}_app
                                                          PRIVATE ${C_FLAGS})
  add_test(CXX_test_run_mwa_beam_${PRECISION}
          test_run_mwa_beam_${PRECISION}_app)

endforeach()
