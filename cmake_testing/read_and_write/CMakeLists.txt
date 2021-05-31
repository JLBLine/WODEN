## Will be doing a bunch of sym linking to files for testing. Make sure we
##get the absolute path OK here
find_path(SRCLIST_PATH "srclist_singlepoint.txt" HINTS "${CMAKE_CURRENT_LIST_DIR}")

# ##Test the east,north,height to X,Y,Z works
add_executable(test_RTS_ENH2XYZ_local_app
    test_RTS_ENH2XYZ_local.c
    ../../src/read_and_write.c
)
target_link_libraries(test_RTS_ENH2XYZ_local_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_RTS_ENH2XYZ_local test_RTS_ENH2XYZ_local_app)



##Make a bunch of symlinks to sky models for testing
if(SRCLIST_PATH)
  message(STATUS "SRCLIST_PATH path: ${SRCLIST_PATH}")
endif()
#
# ## Name a srclists used for testing
set(S1 "srclist_singlepoint.txt")
set(S2 "srclist_singlegauss.txt")
set(S3 "srclist_singleshape.txt")
set(S4 "srclist_threesources.txt")
set(S5 "srclist_threecomponents.txt")
set(S6 "srclist_mulitple_source-components.txt")
set(S7 "srclist_empty_line.txt")
set(S8 "srclist_comment.txt")
set(S9 "srclist_badspell.txt")
set(S10 "srclist_badcoeff.txt")
set(S11 "srclist_no-comp_numbers.txt")


## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRCLIST_PATH}/${X}" "${X}")
endforeach()



# ##Test reading of WODEN catalogue works
add_executable(test_read_source_catalogue_app
    test_read_source_catalogue.c
    ../../src/read_and_write.c
)
target_link_libraries(test_read_source_catalogue_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_read_source_catalogue test_read_source_catalogue_app)



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


## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS J1 J2 J3 J4 J5 J6 J7 J8)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRCLIST_PATH}/${X}" "${X}")
endforeach()



# ##Test reading of WODEN catalogue works
add_executable(test_read_json_settings_app
    test_read_json_settings.c
    ../../src/read_and_write.c
)
target_link_libraries(test_read_json_settings_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_read_json_settings test_read_json_settings_app)

set(A1 "example_array_layout.txt")
## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS A1)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${SRCLIST_PATH}/${X}" "${X}")
endforeach()
#
# ##Test calculating baseline length in X,Y,Z works
add_executable(test_calc_XYZ_diffs_app
    test_calc_XYZ_diffs.c
    ../../src/read_and_write.c
)
target_link_libraries(test_calc_XYZ_diffs_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_calc_XYZ_diffs test_calc_XYZ_diffs_app)

add_executable(test_RTS_PrecessXYZtoJ2000_app
    test_RTS_PrecessXYZtoJ2000.c
    ../../src/read_and_write.c
)
target_link_libraries(test_RTS_PrecessXYZtoJ2000_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_RTS_PrecessXYZtoJ2000 test_RTS_PrecessXYZtoJ2000_app)
