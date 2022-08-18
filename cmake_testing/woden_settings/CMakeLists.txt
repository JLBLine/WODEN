## Will be doing a bunch of sym linking to files for testing. Make sure we
##get the absolute path OK here
find_path(JSON_PATH "run_woden_MWAFEE.json" HINTS "${CMAKE_CURRENT_LIST_DIR}")

##Make a bunch of symlinks to sky models for testing
if(JSON_PATH)
  message(STATUS "JSON_PATH path: ${JSON_PATH}")
endif()

## Name of json inputs used for testing
set(J1 "run_woden_MWAFEE.json")
set(J2 "run_woden_gaussian_default.json")
set(J3 "run_woden_gaussian_bespoke.json")
set(J4 "run_woden_EDA2.json")
set(J5 "run_woden_nobeam.json")
set(J6 "run_woden_multiple_beams.json")
set(J7 "run_woden_MWAFEE_baddelay.json")
set(J7 "run_woden_MWAFEE_baddelay.json")
set(J8 "run_woden_MWAFEE_nopath.json")
set(J9 "example_array_layout.txt")
set(J10 "run_woden_noprecession.json")
set(J11 "run_woden_MWA_analy.json")
set(J12 "run_woden_MWAFEE_interp.json")
set(J12 "run_woden_doautos.json")

## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS J1 J2 J3 J4 J5 J6 J7 J8 J9 J10 J11 J12 J13)
    ADD_CUSTOM_TARGET(targetthis_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${JSON_PATH}/${X}" "${X}")
endforeach()

## Create tests for both float and double precision
foreach(PRECISION IN LISTS FLOAT DOUBLE)
# foreach(PRECISION IN LISTS DOUBLE)

  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(woden_settings_${PRECISION} SHARED "${CMAKE_SOURCE_DIR}/src/woden_settings.c" "${CMAKE_SOURCE_DIR}/src/array_layout.c")
  target_link_libraries(woden_settings_${PRECISION} gcov)
  target_compile_options(woden_settings_${PRECISION} PRIVATE ${C_COVER_FLAGS})
  ## Accumulate the coverage from the tests
  # add_coverage(woden_settings_${PRECISION})

  # ##Test reading of WODEN catalogue works
  add_executable(test_read_json_settings_${PRECISION}_app
      test_read_json_settings.c
  )
  target_link_libraries(test_read_json_settings_${PRECISION}_app
      Unity gcov
      woden_settings_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_read_json_settings_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_read_json_settings_${PRECISION}
           test_read_json_settings_${PRECISION}_app)


  # ##Test reading of WODEN catalogue works
  add_executable(test_setup_lsts_and_phase_centre_${PRECISION}_app
      test_setup_lsts_and_phase_centre.c
  )
  target_link_libraries(test_setup_lsts_and_phase_centre_${PRECISION}_app
      Unity gcov
      woden_settings_${PRECISION}
      ${CC_LINKLIBS}
  )
  target_compile_options(test_setup_lsts_and_phase_centre_${PRECISION}_app PRIVATE ${C_FLAGS})
  add_test(C_test_setup_lsts_and_phase_centre_${PRECISION}
           test_setup_lsts_and_phase_centre_${PRECISION}_app)


endforeach()
