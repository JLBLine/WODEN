# ##Test the horizon test function works
add_executable(test_fill_primary_beam_settings_app
    test_fill_primary_beam_settings.c
    ../../src/primary_beam.c
)
target_link_libraries(test_fill_primary_beam_settings_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_fill_primary_beam_settings test_fill_primary_beam_settings_app)
