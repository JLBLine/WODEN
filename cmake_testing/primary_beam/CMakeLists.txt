## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(primary_beam_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/primary_beam.c")
  target_link_libraries(primary_beam_${PRECISION} gcov)
  target_compile_options(primary_beam_${PRECISION} PRIVATE ${C_COVER_FLAGS})

  ##----------------------------------------------------------------------------
  ##Begin adding tests
  ##----------------------------------------------------------------------------

  ##Test the parallactic angle calculation works
  add_executable(test_calc_para_angle_${PRECISION}_app
      test_calc_para_angle.c
      expected_para_angles.c
  )
  target_link_libraries(test_calc_para_angle_${PRECISION}_app
      Unity gcov
      primary_beam_${PRECISION}
      ${CC_LINKLIBS}
  )

  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_calc_para_angle_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_calc_para_angle_${PRECISION} test_calc_para_angle_${PRECISION}_app)

  ##Test primary beam filling function works
  add_executable(test_fill_primary_beam_settings_${PRECISION}_app
      test_fill_primary_beam_settings.c
      expected_para_angles.c
  )
  target_link_libraries(test_fill_primary_beam_settings_${PRECISION}_app
      Unity gcov
      primary_beam_${PRECISION}
      ${CC_LINKLIBS}
  )
  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_fill_primary_beam_settings_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})

  add_test(C_test_fill_primary_beam_settings_${PRECISION}
           test_fill_primary_beam_settings_${PRECISION}_app)

endforeach()
