add_library(FEE_primary_beam SHARED "${CMAKE_SOURCE_DIR}/src/FEE_primary_beam.c")
## Accumulate the coverage from the tests
add_coverage(FEE_primary_beam)

# ##Test Initialisation of the MWA FEE bem works
add_executable(test_RTS_MWAFEEInit_app
    test_RTS_MWAFEEInit.c
)
target_link_libraries(test_RTS_MWAFEEInit_app
    Unity
    FEE_primary_beam
    ${CC_LINKLIBS}
)
add_test(C_test_RTS_MWAFEEInit test_RTS_MWAFEEInit_app)
