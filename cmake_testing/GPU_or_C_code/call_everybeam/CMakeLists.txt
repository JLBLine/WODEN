# ## Create tests for both float and double precision

# set_source_files_properties(${CMAKE_SOURCE_DIR}/src/call_everybeam.cc
#                             PROPERTIES LANGUAGE CXX)

# ##Test primary beam filling function works
# add_executable(fucking_eh
#     ${CMAKE_SOURCE_DIR}/src/call_everybeam.cc
#     ${CMAKE_SOURCE_DIR}/src/logger.c
# )
# #   ##In the link below, we link against the testing GPU library, even though we don't
# #   ##call it. This is because the CPU library has references to the GPU library,
# #   ##and it will poop the bed if it can't find them
# #   target_link_libraries(test_load_everybeam_telescope_${PRECISION}_app
# #       gcov
# #       use_everybeam
# #   )
# ##Compile using the -DDOUBLE_PRECISION flag if needed
# target_link_libraries(fucking_eh
#                     ${EVERYBEAM_LIB}
#                     ${CASACORE_LIBRARIES})
# #   target_link_options(fucking_eh PRIVATE -fstack-protector-all -fsanitize=address)
# # target_link_options(fucking_eh PRIVATE -fsanitize=thread)

# target_include_directories(fucking_eh PRIVATE
#                 ${EVERYBEAM_INCLUDE_DIRS}
#                 ${AOCOM_BASE}
#                 ${CASACORE_INCLUDE_DIRS})
# target_compile_options(fucking_eh PRIVATE -g -Wall) #-fstack-protector-all -fsanitize=address)
# # target_compile_options(fucking_eh PRIVATE -g -Wall -fsanitize=thread)

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
  ##In the link below, we link against the testing GPU library, even though we don't
  ##call it. This is because the CPU library has references to the GPU library,
  ##and it will poop the bed if it can't find them
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
  add_executable(test_run_lofar_beam_${PRECISION}_app
      test_run_lofar_beam.c
  )
  ##In the link below, we link against the testing GPU library, even though we don't
  ##call it. This is because the CPU library has references to the GPU library,
  ##and it will poop the bed if it can't find them
  target_link_libraries(test_run_lofar_beam_${PRECISION}_app
      Unity gcov
      use_everybeam
  )
  ##Compile using the -DDOUBLE_PRECISION flag if needed
  target_compile_options(test_run_lofar_beam_${PRECISION}_app
                                                             PRIVATE ${C_FLAGS})
  add_test(CXX_test_run_lofar_beam_${PRECISION}
           test_run_lofar_beam_${PRECISION}_app)

  ##Test primary beam filling function works
  add_executable(test_run_mwa_beam_${PRECISION}_app
      test_run_mwa_beam.c
  )

#   target_include_directories(test_run_mwa_beam_${PRECISION}_app PRIVATE 
#       ${CMAKE_SOURCE_DIR}/cmake_testing/GPU_or_C_code/primary_beam/
#   )
  ##In the link below, we link against the testing GPU library, even though we don't
  ##call it. This is because the CPU library has references to the GPU library,
  ##and it will poop the bed if it can't find them
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
