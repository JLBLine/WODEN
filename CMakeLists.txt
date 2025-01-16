cmake_minimum_required(VERSION 3.21)

##Need to find CUDA compiler
include(CheckLanguage)
# check_language(CUDA)

option(USE_HIP "Enable the GPU acceleration for AMD GPUs." OFF)
##default to using CUDA to be somewhat backwards compatible
option(USE_CUDA "Enable the GPU acceleration for NVIDIA GPUs." ON)

if(USE_HIP)
  set(USE_CUDA OFF)
endif()

if(USE_CUDA)
    message("USING CUDA : USE_CUDA=ON")
    project(woden VERSION 2.3 LANGUAGES CXX C CUDA)

    ##Sets the -std=c++11 nvcc compiler flag I think
    set(CMAKE_CUDA_STANDARD 11)
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 ")
    ##Internal flag to let code know we're using nvcc
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__NVCC__ ")
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    if (NOT DEFINED ENV{CUDAARCHS})
        message("WARNING : CUDAARCHS not set, so CMAKE_CUDA_ARCHITECTURES has defaulted to whatever CMake has set, which is ${CMAKE_CUDA_ARCHITECTURES}")
    else()
        message("INFO : CUDAARCHS was set to $ENV{CUDAARCHS}, so now CMAKE_CUDA_ARCHITECTURES is set to ${CMAKE_CUDA_ARCHITECTURES}")
    endif()
    

elseif(USE_HIP)
    message("DEBUG : USE_HIP=ON")

    project(woden VERSION 2.3 LANGUAGES CXX C HIP)
    
    set(CMAKE_CXX_COMPILER hipcc)
    add_definitions("-D__HIP_PLATFORM_AMD__ -D__HIPCC__")

    ##TODO, remove this if -DGPU_TARGETS="gfx1032;gfx1035" works
    if(HIP_ARCH)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__HIP_PLATFORM_AMD__ -D__HIPCC__ --offload-arch=${HIP_ARCH}")
        message("-- Using user specified HIP ARCH of --offload-arch=${HIP_ARCH}")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__HIP_PLATFORM_AMD__ -D__HIPCC__")
    endif()

endif()

##Try and find the hyperbeam library
set(HINTS_HBEAM_INC "/usr/local/include/;${HBEAM_ROOT}/;")
set(HINTS_HBEAM_LIB "/usr/local/lib/;${HBEAM_ROOT}")

find_path(HBEAM_INC mwa_hyperbeam.h HINTS ${HINTS_HBEAM_INC})
find_library(HBEAM_LIB mwa_hyperbeam HINTS ${HINTS_HBEAM_LIB} REQUIRED)

##Find the source codes
file(GLOB CSOURCES "src/*.c")
file(GLOB GPUSOURCES "src/*.cpp")

if(USE_CUDA)
    set_source_files_properties( ${GPUSOURCES} PROPERTIES LANGUAGE CUDA)
endif()


include_directories("src/")

if(HBEAM_LIB)
  message(STATUS "HBEAM include path: ${HBEAM_INC}")
  message(STATUS "HBEAM library path: ${HBEAM_LIB}")
endif()

##Set the linkers for the C code
set(CC_LINKLIBS "-lm -pthread ${HDF5_LIB} ${HBEAM_LIB}")

##All the WODEN header dirs are in include dir
include_directories("${CMAKE_SOURCE_DIR}/include" ${HBEAM_INC})

##Choose whether we are created the main executable, or making the unit tests
set(TARGET_GROUP production CACHE STRING "Group to build")

##If making production code
if(TARGET_GROUP STREQUAL production)
    ##Do it for the float code first--------------------------------------------
    ##Generate library for the CUDA code
    add_library(wodenGPU_float SHARED ${GPUSOURCES})
    ##Compile CUDA code with all warnings
    if(USE_CUDA)
       target_compile_options(wodenGPU_float PRIVATE --compiler-options -Wall)
    elseif(USE_HIP)
       target_compile_options(wodenGPU_float PRIVATE -Wall)
    endif()

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_float SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_float PUBLIC wodenGPU_float ${CC_LINKLIBS})
    target_compile_options(woden_float PRIVATE -O3 -Wall)

    ##Copy things into the wodenpy dir, which means everything is found when
    ##you do a pip install
    add_custom_command(TARGET woden_float POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libwoden_float.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libwoden_float.so")

    ##Now do it for the double code by adding the -DDOUBLE_PRECISION compiler
    ##flag----------------------------------------------------------------------

    ##Generate library for the CUDA code
    add_library(wodenGPU_double SHARED ${GPUSOURCES})
    ##Compile CUDA code with all warnings
    if(USE_CUDA)
       target_compile_options(wodenGPU_double PRIVATE -DDOUBLE_PRECISION --compiler-options -Wall)
    elseif(USE_HIP)
        target_compile_options(wodenGPU_double PRIVATE -DDOUBLE_PRECISION -Wall)
    endif()

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_double SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_double PUBLIC wodenGPU_double ${CC_LINKLIBS})
    target_compile_options(woden_double PRIVATE -O3 -Wall -DDOUBLE_PRECISION)

    ##Copy things into the wodenpy dir, which means everything is found when
    ##you do a pip install
    add_custom_command(TARGET woden_double POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libwoden_double.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libwoden_double.so")


##If debugging build
elseif(TARGET_GROUP STREQUAL debug)
    ##Do it for the float code first--------------------------------------------

    ##Generate library for the CUDA code
    add_library(wodenGPU_float SHARED ${GPUSOURCES})
    ##Compile CUDA code with all warnings
    if(USE_CUDA)
        target_compile_options(wodenGPU_float PRIVATE -G --compiler-options -g)
    elseif(USE_HIP)
        target_compile_options(wodenGPU_float PRIVATE -G -g)
    endif()
    
    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_float SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_float PUBLIC wodenGPU_float ${CC_LINKLIBS})
    target_compile_options(woden_float PRIVATE -O0 -Wall)

    ##Copy things into the wodenpy dir, which means everything is found when
    ##you do a pip install
    add_custom_command(TARGET woden_float POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libwoden_float.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libwoden_float.so")

    ##Do it for the double code now---------------------------------------------
    ##Generate library for the CUDA code
    add_library(wodenGPU_double SHARED ${GPUSOURCES})
    ##Compile CUDA code with all warnings
    target_compile_options(wodenGPU_double PRIVATE -DDOUBLE_PRECISION
                                                       -G --compiler-options -g)
    if(USE_CUDA)
        target_compile_options(wodenGPU_double PRIVATE -DDOUBLE_PRECISION
                               -G --compiler-options -g)
    elseif(USE_HIP)
        target_compile_options(wodenGPU_double PRIVATE -DDOUBLE_PRECISION -G -g)
    endif()

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_double SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_double PUBLIC wodenGPU_double ${CC_LINKLIBS})
    target_compile_options(woden_double PRIVATE -O0 -g -Wall -DDOUBLE_PRECISION)

    ##Copy things into the wodenpy dir, which means everything is found when
    ##you do a pip install
    add_custom_command(TARGET woden_double POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libwoden_double.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libwoden_double.so")

##If making test code
elseif(TARGET_GROUP STREQUAL test)

      include(CTest)

      ##Find out where unity is hiding
      set(HINTS_UNITY "${UNITY_ROOT};${UNITY_ROOT}/src/;/usr/local/Unity/src/")
      find_path(UNITY_DIR unity.c HINTS ${HINTS_UNITY})
      if(UNITY_DIR)
        message(STATUS "UNITY src path: ${UNITY_DIR}")
      endif()

      ##Compile Unity
      add_library(Unity STATIC
          "${UNITY_DIR}/unity.c"
      )
      target_include_directories(Unity PUBLIC
          "${UNITY_DIR}"
      )
      ##Enable DOUBLE precision tests (whhhhhy is this not default CHRIST)
      target_compile_options(Unity PUBLIC
      -DUNITY_INCLUDE_DOUBLE -DUNITY_DOUBLE_PRECISION=1e-12)

      ##Inlcude all the testing that lives in ./cmake_testing
      add_subdirectory(cmake_testing)

      # coverage_evaluate()

else()
    message(FATAL_ERROR "Given TARGET_GROUP unknown")
endif()
