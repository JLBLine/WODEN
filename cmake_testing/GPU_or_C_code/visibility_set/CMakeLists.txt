## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  ##Test settings for visibility_set are added correctly
  add_executable(test_fill_timefreq_visibility_set_${PRECISION}_app
      test_fill_timefreq_visibility_set.c
  )
  target_link_libraries(test_fill_timefreq_visibility_set_${PRECISION}_app
      Unity gcov
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_fill_timefreq_visibility_set_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_fill_timefreq_visibility_set_${PRECISION}
           test_fill_timefreq_visibility_set_${PRECISION}_app)

  ##Test some malloc and freeing functions
  add_executable(test_malloc_and_free_${PRECISION}_app
      test_malloc_and_free.c
  )
  target_link_libraries(test_malloc_and_free_${PRECISION}_app
      Unity gcov
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_malloc_and_free_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_malloc_and_free_${PRECISION}
           test_malloc_and_free_${PRECISION}_app)

  ##Test writing the binary file
  add_executable(test_write_visi_set_binary_${PRECISION}_app
      test_write_visi_set_binary.c
  )
  target_link_libraries(test_write_visi_set_binary_${PRECISION}_app
      Unity gcov
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_write_visi_set_binary_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_write_visi_set_binary_${PRECISION}
           test_write_visi_set_binary_${PRECISION}_app)

  ##Test writing to a text file
  add_executable(test_write_visi_set_text_${PRECISION}_app
      test_write_visi_set_text.c
  )
  target_link_libraries(test_write_visi_set_text_${PRECISION}_app
      Unity gcov
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_write_visi_set_text_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_write_visi_set_text_${PRECISION}
           test_write_visi_set_text_${PRECISION}_app)

endforeach()
