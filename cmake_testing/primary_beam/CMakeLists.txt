add_library(primary_beam SHARED "${CMAKE_SOURCE_DIR}/src/primary_beam.c")
## Accumulate the coverage from the tests
add_coverage(primary_beam)

##Test the parallactic angle calculation works
add_executable(test_calc_para_angle_app
    test_calc_para_angle.c
    expected_para_angles.c
)
target_link_libraries(test_calc_para_angle_app
    Unity
    primary_beam
    ${CC_LINKLIBS}
)
add_test(C_test_calc_para_angle test_calc_para_angle_app)

##Test primary beam filling function works
add_executable(test_fill_primary_beam_settings_app
    test_fill_primary_beam_settings.c
    expected_para_angles.c
)
target_link_libraries(test_fill_primary_beam_settings_app
    Unity
    primary_beam
    ${CC_LINKLIBS}
)
add_test(C_test_fill_primary_beam_settings test_fill_primary_beam_settings_app)
