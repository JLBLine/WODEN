foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)
  ##----------------------------------------------------------------------------
  ##Begin adding tests
  ##----------------------------------------------------------------------------

  ##Test primary beam filling function works
  add_executable(test_woden_logger_${PRECISION}_app
      test_woden_logger.c
  )
  ##In the link below, we link against the testing GPU library, even though we don't
  ##call it. This is because the CPU library has references to the GPU library,
  ##and it will poop the bed if it can't find them
  target_link_libraries(test_woden_logger_${PRECISION}_app
      Unity gcov
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
  )
  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_woden_logger_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_woden_logger_${PRECISION}
           test_woden_logger_${PRECISION}_app)


  ##Test primary beam filling function works
  add_executable(test_hyperbeam_error_${PRECISION}_app
    test_hyperbeam_error.c
    # ${CMAKE_SOURCE_DIR}/src/logger.c
    # ${CMAKE_SOURCE_DIR}/src/hyperbeam_error.c
  )
  ##In the link below, we link against the testing GPU library, even though we don't
  ##call it. This is because the CPU library has references to the GPU library,
  ##and it will poop the bed if it can't find them
  target_link_libraries(test_hyperbeam_error_${PRECISION}_app
      Unity gcov
      calculate_visibilities_CPU_${PRECISION}
      calculate_visibilities_GPU_${PRECISION}
      ${CC_LINKLIBS}
  )

  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_hyperbeam_error_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_hyperbeam_error_${PRECISION}
          test_hyperbeam_error_${PRECISION}_app)



  if(USE_CUDA)
      set_source_files_properties(bad_gpu_code.cpp
                                  PROPERTIES LANGUAGE CUDA)
            
  endif()

  add_library(bad_cpu_code_${PRECISION} SHARED
              bad_gpu_code.cpp)

  target_compile_options(bad_cpu_code_${PRECISION} PRIVATE ${GPU_FLAGS})

  ##Test primary beam filling function works
  add_executable(test_log_gpu_errors_${PRECISION}_app
    test_log_gpu_errors.c
    ${CMAKE_SOURCE_DIR}/src/logger.c
  )
  ##In the link below, we link against the testing GPU library, even though we don't
  ##call it. This is because the CPU library has references to the GPU library,
  ##and it will poop the bed if it can't find them
  target_link_libraries(test_log_gpu_errors_${PRECISION}_app
      Unity gcov bad_cpu_code_${PRECISION}
      ${CC_LINKLIBS}
  )

  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_log_gpu_errors_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_log_gpu_errors_${PRECISION}
          test_log_gpu_errors_${PRECISION}_app)

endforeach()
