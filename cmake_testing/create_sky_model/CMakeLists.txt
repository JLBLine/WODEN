# ## Will be doing a bunch of sym linking to files for testing. Make sure we
# ##get the absolute path OK here
find_path(TXT_PATH "srclist_singlepoint.txt" HINTS "${CMAKE_SOURCE_DIR}/cmake_testing/read_text_skymodel/")
find_path(YAML_PATH "srclist_singlepoint_power.yaml" HINTS "${CMAKE_SOURCE_DIR}/cmake_testing/read_text_skymodel/")

ADD_CUSTOM_TARGET(link_srclist_singlepoint.txt ALL
   COMMAND ${CMAKE_COMMAND} -E create_symlink "${TXT_PATH}/srclist_singlepoint.txt"
 "srclist_singlepoint_power.txt")

ADD_CUSTOM_TARGET(link_srclist_singlepoint.yaml ALL
   COMMAND ${CMAKE_COMMAND} -E create_symlink "${YAML_PATH}/srclist_singlepoint_power.yaml"
 "srclist_singlepoint.yaml")

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(create_sky_model_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/create_sky_model.c"
  "${CMAKE_SOURCE_DIR}/src/read_text_skymodel.c"
  "${CMAKE_SOURCE_DIR}/src/read_yaml_skymodel.c")
  target_link_libraries(create_sky_model_${PRECISION} gcov)
  target_include_directories(create_sky_model_${PRECISION} PUBLIC "${CMAKE_SOURCE_DIR}/cmake_testing/create_sky_model/"
  "${CMAKE_SOURCE_DIR}/include")
  target_compile_options(create_sky_model_${PRECISION} PRIVATE ${C_COVER_FLAGS})

  ##Test reading of WODEN catalogue works
  add_executable(test_read_source_catalogue_${PRECISION}_app
      test_read_source_catalogue.c
      "${CMAKE_SOURCE_DIR}/cmake_testing/create_sky_model/expected_skymodel_outcomes.c"
  )
  target_link_libraries(test_read_source_catalogue_${PRECISION}_app
      Unity gcov
      create_sky_model_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_read_source_catalogue_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_read_source_catalogue_${PRECISION} test_read_source_catalogue_${PRECISION}_app)

  # ##Test the horizon test function works
  # add_executable(test_horizon_test_${PRECISION}_app
  #     test_horizon_test.c
  # )
  # target_link_libraries(test_horizon_test_${PRECISION}_app
  #     Unity gcov
  #     create_sky_model_${PRECISION}
  #     ${CC_LINKLIBS}
  # )
  # target_compile_options(test_horizon_test_${PRECISION}_app PRIVATE ${C_FLAGS})
  # add_test(C_test_horizon_test_${PRECISION} test_horizon_test_${PRECISION}_app)

  ##Test the horizon cropping function creates sky models correctly
  add_executable(test_crop_sky_model_${PRECISION}_app
      test_crop_sky_model.c
  )
  target_link_libraries(test_crop_sky_model_${PRECISION}_app
      Unity gcov
      create_sky_model_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_crop_sky_model_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_crop_sky_model_${PRECISION} test_crop_sky_model_${PRECISION}_app)

endforeach()
