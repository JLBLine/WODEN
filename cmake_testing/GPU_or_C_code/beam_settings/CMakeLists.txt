## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)
  ##----------------------------------------------------------------------------
  ##Begin adding tests
  ##----------------------------------------------------------------------------

  ##Test primary beam filling function works
  add_executable(test_fill_primary_beam_settings_${PRECISION}_app
      test_fill_primary_beam_settings.c
  )
  ##In the link below, we link against the testing GPU library, even though we don't
  ##call it. This is because the CPU library has references to the GPU library,
  ##and it will poop the bed if it can't find them
  target_link_libraries(test_fill_primary_beam_settings_${PRECISION}_app
      Unity gcov
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_fill_primary_beam_settings_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(C_test_fill_primary_beam_settings_${PRECISION}
           test_fill_primary_beam_settings_${PRECISION}_app)

endforeach()
