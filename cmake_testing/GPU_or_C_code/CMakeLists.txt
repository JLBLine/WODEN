
##To make coverage reporting more straight forward, compile all the code
##into one big ol' library here. The link any tests against this library. As
##we call all our tests in serial, each test should just add to the coverage

##Currently, the C++ code we use is always double, so only need to 
##compile one library. DO that first, but only if we have the EveryBeam
##library to link against. The variable USE_EVERYBEAM is set in the
##top level CMakeLists.txt file

if (USE_EVERYBEAM)
    add_library(use_everybeam SHARED
                ${CXXSOURCES}
                ${CMAKE_SOURCE_DIR}/src/logger.c
                )
    target_link_libraries(use_everybeam gcov ${EVERYBEAM_LIB} ${CASACORE_LIBRARIES})

    target_include_directories(use_everybeam PRIVATE
                            ${EVERYBEAM_INCLUDE_DIRS}
                            ${AOCOM_BASE}
                            ${CASACORE_INCLUDE_DIRS})
    target_compile_options(use_everybeam PRIVATE -DHAVE_EVERYBEAM -coverage -g -Wall)
endif()


foreach(PRECISION IN LISTS FLOAT DOUBLE)
    DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS GPU_FLAGS C_COVER_FLAGS)

    ##If we are using CUDA, we have to tell CMake that the files are CUDA
    ##and to use the CUDA compiler. We've already to CMake to use hipcc
    ##in the case of HIP.
    if(USE_CUDA)
        set_source_files_properties(${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp
              ${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp
              ${CMAKE_SOURCE_DIR}/src/source_components_gpu.cpp
              ${CMAKE_SOURCE_DIR}/src/calculate_visibilities_gpu.cpp
              PROPERTIES LANGUAGE CUDA)
              
    endif()

    add_library(calculate_visibilities_GPU_${PRECISION} SHARED
                ${CMAKE_SOURCE_DIR}/src/fundamental_coords_gpu.cpp
                ${CMAKE_SOURCE_DIR}/src/primary_beam_gpu.cpp
                ${CMAKE_SOURCE_DIR}/src/source_components_gpu.cpp
                ${CMAKE_SOURCE_DIR}/src/calculate_visibilities_gpu.cpp)
    target_compile_options(calculate_visibilities_GPU_${PRECISION} PRIVATE ${GPU_FLAGS})
    if (USE_EVERYBEAM)
        target_compile_options(calculate_visibilities_GPU_${PRECISION} PRIVATE
                              -DHAVE_EVERYBEAM)
    endif()

    add_library(calculate_visibilities_CPU_${PRECISION} SHARED
                ${CSOURCES})
    if (USE_EVERYBEAM)
        target_link_libraries(calculate_visibilities_CPU_${PRECISION}
                              gcov ${CC_LINKLIBS} use_everybeam)
        target_compile_options(calculate_visibilities_CPU_${PRECISION} PRIVATE
                               ${C_COVER_FLAGS} -DHAVE_EVERYBEAM)
    else()
        target_link_libraries(calculate_visibilities_CPU_${PRECISION}
                              gcov ${CC_LINKLIBS})
    endif()
    target_compile_options(calculate_visibilities_CPU_${PRECISION} PRIVATE ${C_COVER_FLAGS})

endforeach()

add_subdirectory(beam_settings)
add_subdirectory(calculate_visibilities)
add_subdirectory(fundamental_coords)
add_subdirectory(logger)
add_subdirectory(primary_beam)
add_subdirectory(run_woden)
add_subdirectory(source_components)
add_subdirectory(visibility_set)

if (USE_EVERYBEAM)
    add_subdirectory(call_everybeam)
endif()