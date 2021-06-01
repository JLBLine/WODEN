# ##Test the horizon test function works
add_executable(test_horizon_test_app
    test_horizon_test.c
    ../../src/create_sky_model.c
)
target_link_libraries(test_horizon_test_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_horizon_test test_horizon_test_app)

# ##Test the horizon cropping function creates sky models correctly
add_executable(test_crop_sky_model_app
    test_crop_sky_model.c
    ../../src/create_sky_model.c
)
target_link_libraries(test_crop_sky_model_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_crop_sky_model test_crop_sky_model_app)
