## Will be doing a bunch of sym linking to files for testing. Make sure we
##get the absolute path OK here
find_path(LAYOUT_PATH "example_array_layout.txt" HINTS "${CMAKE_CURRENT_LIST_DIR}")

# ##Test the east,north,height to X,Y,Z works
add_executable(test_RTS_ENH2XYZ_local_app
    test_RTS_ENH2XYZ_local.c
    ../../src/array_layout.c
)
target_link_libraries(test_RTS_ENH2XYZ_local_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_RTS_ENH2XYZ_local test_RTS_ENH2XYZ_local_app)

set(A1 "example_array_layout.txt")
## Loop through them and write a symlink so the testing executable can find them
foreach(X IN LISTS A1)
    ADD_CUSTOM_TARGET(target_${X} ALL
        COMMAND ${CMAKE_COMMAND} -E create_symlink "${LAYOUT_PATH}/${X}" "${X}")
endforeach()
#
# ##Test calculating baseline length in X,Y,Z works
add_executable(test_calc_XYZ_diffs_app
    test_calc_XYZ_diffs.c
    ../../src/array_layout.c
)
target_link_libraries(test_calc_XYZ_diffs_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_calc_XYZ_diffs test_calc_XYZ_diffs_app)

add_executable(test_RTS_PrecessXYZtoJ2000_app
    test_RTS_PrecessXYZtoJ2000.c
    ../../src/array_layout.c
)
target_link_libraries(test_RTS_PrecessXYZtoJ2000_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_RTS_PrecessXYZtoJ2000 test_RTS_PrecessXYZtoJ2000_app)
