# ##Test Initialisation of the MWA FEE bem works
add_executable(test_RTS_MWAFEEInit_app
    test_RTS_MWAFEEInit.c
    ../../src/FEE_primary_beam.c
)
target_link_libraries(test_RTS_MWAFEEInit_app
    Unity
    ${CC_LINKLIBS}
)
add_test(C_test_RTS_MWAFEEInit test_RTS_MWAFEEInit_app)
