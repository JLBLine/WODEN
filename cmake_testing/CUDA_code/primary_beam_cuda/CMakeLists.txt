
## Create tests for both float and double precision
# foreach(PRECISION IN LISTS FLOAT DOUBLE)
foreach(PRECISION IN LISTS DOUBLE)
  ##Set the flags depending on the precision
  DEFINE_COMP_FLAGS(${PRECISION} C_FLAGS CUDA_FLAGS C_COVER_FLAGS)

  add_library(primary_beam_cudaCUDA_${PRECISION} SHARED
              "${CMAKE_SOURCE_DIR}/src/primary_beam_cuda.cu"
              "${CMAKE_SOURCE_DIR}/src/fundamental_coords.cu")

  target_compile_options(primary_beam_cudaCUDA_${PRECISION} PRIVATE ${CUDA_FLAGS})

#   ## For CMake versions > 3.18, the arch flag settings are handled by
#   ## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
#   ## replace the CUDA_FLAG lines above with something like
#   set_property(TARGET primary_beam_cudaCUDA_${PRECISION} PROPERTY CUDA_ARCHITECTURES OFF)

#   ##Test gaussian beam functions
#   add_executable(test_gaussian_beam_${PRECISION}_app
#       test_gaussian_beam.c
#   )
#   target_link_libraries(test_gaussian_beam_${PRECISION}_app
#       Unity
#       primary_beam_cudaCUDA_${PRECISION}
#       ${CC_LINKLIBS}
#   )
#   target_compile_options(test_gaussian_beam_${PRECISION}_app PRIVATE ${C_FLAGS})
#   add_test(CUDA_test_gaussian_beam_${PRECISION} test_gaussian_beam_${PRECISION}_app)
  
  
  
#   ##Test analytic dipole beam functions
#   add_executable(test_analytic_dipole_beam_${PRECISION}_app
#       test_analytic_dipole_beam.c
#   )
#   target_link_libraries(test_analytic_dipole_beam_${PRECISION}_app
#       Unity
#       primary_beam_cudaCUDA_${PRECISION}
#       ${CC_LINKLIBS}
#   )
#   target_compile_options(test_analytic_dipole_beam_${PRECISION}_app PRIVATE ${C_FLAGS})
#   add_test(CUDA_test_analytic_dipole_beam_${PRECISION} test_analytic_dipole_beam_${PRECISION}_app)
  
  
  
#   ##Test analytic MWA beam functions
#   add_executable(test_MWA_analytic_${PRECISION}_app
#       test_MWA_analytic.c
#   )
#   target_link_libraries(test_MWA_analytic_${PRECISION}_app
#       Unity
#       primary_beam_cudaCUDA_${PRECISION}
#       ${CC_LINKLIBS}
#   )
#   target_compile_options(test_MWA_analytic_${PRECISION}_app PRIVATE ${C_FLAGS})
#   add_test(CUDA_test_MWA_analytic_${PRECISION} test_MWA_analytic_${PRECISION}_app)


#   #Test hyperbeam with the coarse hdf5 file
#   add_executable(test_run_hyperbeam_${PRECISION}_app
#       test_run_hyperbeam.c
#   )
#   target_link_libraries(test_run_hyperbeam_${PRECISION}_app
#       Unity
#       primary_beam_cudaCUDA_${PRECISION}
#       ${CC_LINKLIBS}
#       # "-L/usr/local/mwa_hyperbeam -lmwa_hyperbeam"
#   )
#   target_compile_options(test_run_hyperbeam_${PRECISION}_app PRIVATE ${C_FLAGS})
#   add_test(CUDA_test_run_hyperbeam_${PRECISION} test_run_hyperbeam_${PRECISION}_app)


#   #Test hyperbeam with the interpolated beam
#   add_executable(test_run_hyperbeam_interp_${PRECISION}_app
#       test_run_hyperbeam_interp.c
#   )
#   target_link_libraries(test_run_hyperbeam_interp_${PRECISION}_app
#       Unity
#       primary_beam_cudaCUDA_${PRECISION}
#       ${CC_LINKLIBS}
#       # "-L/usr/local/mwa_hyperbeam -lmwa_hyperbeam"
#   )
#   target_compile_options(test_run_hyperbeam_interp_${PRECISION}_app PRIVATE ${C_FLAGS})
#   add_test(CUDA_test_run_hyperbeam_interp_${PRECISION} test_run_hyperbeam_interp_${PRECISION}_app)


    #Test hyperbeam with the interpolated beam
    add_executable(test_run_hyperbeam_two_antennas_${PRECISION}_app
    test_run_hyperbeam_two_antennas.c
    )
    target_link_libraries(test_run_hyperbeam_two_antennas_${PRECISION}_app
        Unity
        primary_beam_cudaCUDA_${PRECISION}
        ${CC_LINKLIBS}
        # "-L/usr/local/mwa_hyperbeam -lmwa_hyperbeam"
    )
    target_compile_options(test_run_hyperbeam_two_antennas_${PRECISION}_app PRIVATE ${C_FLAGS})
    add_test(CUDA_test_run_hyperbeam_two_antennas_${PRECISION} test_run_hyperbeam_two_antennas_${PRECISION}_app)

endforeach()
