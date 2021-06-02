# ##Test the horizon test function works
add_executable(test_null_comps_app
    test_null_comps.c
    ../../src/chunk_sky_model.c
)
target_link_libraries(test_null_comps_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_null_comps test_null_comps_app)
