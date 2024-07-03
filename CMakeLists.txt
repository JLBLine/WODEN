cmake_minimum_required(VERSION 3.10)

##Need to find CUDA compiler
include(CheckLanguage)
check_language(CUDA)

option(USE_HIP "Enable the GPU acceleration for AMD GPUs." OFF)
option(USE_CUDA "Enable the GPU acceleration for NVIDIA GPUs." OFF)

if(USE_CUDA)
    message("DEBUG : USE_CUDA=ON")
    project(woden VERSION 2.2 LANGUAGES CXX C CUDA) 

    ##Gotta find CUDA and get it to find libraries
    find_package(CUDA)
    if(CUDA_FOUND)
       message(STATUS "Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")
       message(STATUS "CUDA library found; building GPU-enabled executables.")
       ##Sets the -std=c++11 nvcc compiler flag I think
       set(CMAKE_CUDA_STANDARD 11)
       #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 ")
       ##Internal flag to let code know we're using nvcc
       set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__NVCC__ ")
       set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif()

    ##Pull out arch version if specified by user, and pass to nvcc
    if (CUDA_ARCH EQUAL "1.3")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_13")
    elseif (CUDA_ARCH EQUAL "2.0")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_20")
    elseif (CUDA_ARCH EQUAL "2.1")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_21")
    elseif (CUDA_ARCH EQUAL "3.0")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_30")
    elseif (CUDA_ARCH EQUAL "3.2")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_32")
    elseif (CUDA_ARCH EQUAL "3.5")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_35")
    elseif (CUDA_ARCH EQUAL "3.7")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_37")
    elseif (CUDA_ARCH EQUAL "5.0")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_50")
    elseif (CUDA_ARCH EQUAL "5.2")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_52")
    elseif (CUDA_ARCH EQUAL "6.0")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_60")
    elseif (CUDA_ARCH EQUAL "6.1")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_61")
    elseif (CUDA_ARCH EQUAL "6.2")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_62")
    elseif (CUDA_ARCH EQUAL "7.0")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_70")
    elseif (CUDA_ARCH EQUAL "7.5")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch=sm_75")
    endif()

    set(HINTS_HBEAM_INC "/usr/local/include/;${HBEAM_ROOT}/;")
    set(HINTS_HBEAM_LIB "/usr/local/lib/;${HBEAM_ROOT}")
    
    find_path(HBEAM_INC mwa_hyperbeam.h HINTS ${HINTS_HBEAM_INC})
    find_library(HBEAM_LIB mwa_hyperbeam HINTS ${HINTS_HBEAM_LIB} REQUIRED)
elseif(USE_HIP)
    message("DEBUG : USE_HIP=ON")

    project(woden VERSION 2.2 LANGUAGES CXX C)
    set(CMAKE_CXX_COMPILER hipcc)
    add_definitions("-D__HIP_PLATFORM_AMD__ -D__HIPCC__")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__HIP_PLATFORM_AMD__ -D__HIPCC__ ")
    # TODO should try and get --offload-arch=gfx90a in there for HIP
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__HIP_PLATFORM_AMD__ -D__HIPCC__")

    ##Try and find the hyperbeam library
    set(HINTS_HBEAM_INC "/usr/local/include/;${HBEAM_ROOT}/;")
    set(HINTS_HBEAM_LIB "/usr/local/lib/;${HBEAM_ROOT}")

    find_path(HBEAM_INC mwa_hyperbeam.h HINTS ${HINTS_HBEAM_INC})
    find_library(HBEAM_LIB mwa_hyperbeam HINTS ${HINTS_HBEAM_LIB} REQUIRED)
endif()

##Find the source codes
file(GLOB CSOURCES "src/*.c")
file(GLOB CUDASOURCES "src/*.cpp")

if(USE_CUDA)
    set_source_files_properties( ${CUDASOURCES} PROPERTIES LANGUAGE CUDA)
endif()


include_directories("src/")

if(HBEAM_LIB)
  message(STATUS "HBEAM include path: ${HBEAM_INC}")
  message(STATUS "HBEAM library path: ${HBEAM_LIB}")
endif()

##Set the linkers for the C code
set(CC_LINKLIBS "-lm -pthread ${HDF5_LIB} ${HBEAM_LIB}") #

##All the WODEN header dirs are in include dir
include_directories("${CMAKE_SOURCE_DIR}/include" ${HBEAM_INC})

##Choose whether we are created the main executable, or making the unit tests
set(TARGET_GROUP production CACHE STRING "Group to build")

##If making production code
if(TARGET_GROUP STREQUAL production)
    ##Do it for the float code first--------------------------------------------
    ##Generate library for the CUDA code
    add_library(wodenCUDA_float SHARED ${CUDASOURCES})
    ##Compile CUDA code with all warnings
    if(USE_CUDA)
       target_compile_options(wodenCUDA_float PRIVATE --compiler-options -Wall)
    endif()
    ## For CMake versions > 3.18, the arch flag settings are handled by
    ## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
    ## replace the CUDA_FLAG lines above with something like
    set_property(TARGET wodenCUDA_float PROPERTY CUDA_ARCHITECTURES OFF)

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_float SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_float PUBLIC wodenCUDA_float ${CC_LINKLIBS})
    target_compile_options(woden_float PRIVATE -O3 -Wall)

    ##Copy things into the wodenpy dir, which means everything is found when
    ##you do a pip install
    add_custom_command(TARGET woden_float POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libwoden_float.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libwoden_float.so")

    ##Now do it for the double code by adding the -DDOUBLE_PRECISION compiler
    ##flag----------------------------------------------------------------------

    ##Generate library for the CUDA code
    add_library(wodenCUDA_double SHARED ${CUDASOURCES})
    ##Compile CUDA code with all warnings
    if(USE_CUDA)
       target_compile_options(wodenCUDA_double PRIVATE -DDOUBLE_PRECISION --compiler-options -Wall)
    endif()
    ## For CMake versions > 3.18, the arch flag settings are handled by
    ## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
    ## replace the CUDA_FLAG lines above with something like
    set_property(TARGET wodenCUDA_double PROPERTY CUDA_ARCHITECTURES OFF)

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_double SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_double PUBLIC wodenCUDA_double ${CC_LINKLIBS})
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
    add_library(wodenCUDA_float SHARED ${CUDASOURCES})
    ##Compile CUDA code with all warnings
    target_compile_options(wodenCUDA_float PRIVATE -G --compiler-options -g)
    ## For CMake versions > 3.18, the arch flag settings are handled by
    ## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
    ## replace the CUDA_FLAG lines above with something like
    set_property(TARGET wodenCUDA_float PROPERTY CUDA_ARCHITECTURES OFF)

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_float SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_float PUBLIC wodenCUDA_float ${CC_LINKLIBS})
    target_compile_options(woden_float PRIVATE -O3 -Wall)

    ##Do it for the double code now---------------------------------------------
    ##Generate library for the CUDA code
    add_library(wodenCUDA_double SHARED ${CUDASOURCES})
    ##Compile CUDA code with all warnings
    target_compile_options(wodenCUDA_double PRIVATE -DDOUBLE_PRECISION
                                                       -G --compiler-options -g)
    ## For CMake versions > 3.18, the arch flag settings are handled by
    ## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
    ## replace the CUDA_FLAG lines above with something like
    set_property(TARGET wodenCUDA_double PROPERTY CUDA_ARCHITECTURES OFF)

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_double SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_double PUBLIC wodenCUDA_double ${CC_LINKLIBS})
    target_compile_options(woden_double PRIVATE -O3 -g -Wall -DDOUBLE_PRECISION)

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
