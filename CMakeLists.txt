cmake_minimum_required(VERSION 3.21)

set(WODEN_VERSION 2.7)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckLanguage)

option(USE_HIP "Enable the GPU acceleration for AMD GPUs." OFF)
##default to using CUDA to be somewhat backwards compatible
option(USE_CUDA "Enable the GPU acceleration for NVIDIA GPUs." ON)

if(USE_HIP)
  set(USE_CUDA OFF)
endif()

if(USE_CUDA)
    message("USING CUDA : USE_CUDA=ON")

    if (NOT DEFINED ENV{CUDAARCHS} AND NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        message(STATUS "CUDAARCHS/CMAKE_CUDA_ARCHITECTURES so using CMAKE default architectures for CUDA")
        ##If you want to build a bunch of architectures by default, do this
        #set(CMAKE_CUDA_ARCHITECTURES "60;61;70;75;80;86;89;90")
    elseif(DEFINED ENV{CUDAARCHS})
        set(CMAKE_CUDA_ARCHITECTURES $ENV{CUDAARCHS})
        message(STATUS "CUDAARCHS was set to $ENV{CUDAARCHS}, so now CMAKE_CUDA_ARCHITECTURES is set to ${CMAKE_CUDA_ARCHITECTURES}")
    else()
        message(STATUS "CMAKE_CUDA_ARCHITECTURES was set to ${CMAKE_CUDA_ARCHITECTURES}, so using that")
    endif()

    project(woden VERSION ${WODEN_VERSION} LANGUAGES CXX C CUDA)

    ##Sets the -std=c++11 nvcc compiler flag I think
    set(CMAKE_CUDA_STANDARD 17)
    ##Internal flag to let code know we're using nvcc
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D__NVCC__ ")
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

elseif(USE_HIP)
    message("DEBUG : USE_HIP=ON")

    project(woden VERSION ${WODEN_VERSION} LANGUAGES CXX C HIP)
    
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

##Find the source codes
file(GLOB CSOURCES "src/*.c")
file(GLOB GPUSOURCES "src/*.cpp")

##because of HIP, we use .cpp for GPU code. Here, use .cc for c++
##Specify the language to the compiler
file(GLOB CXXSOURCES "src/*.cc")
set_source_files_properties( ${CXXSOURCES} PROPERTIES LANGUAGE CXX)

if(USE_CUDA)
    set_source_files_properties( ${GPUSOURCES} PROPERTIES LANGUAGE CUDA)
endif()

include_directories("src/")

##Try and find the hyperbeam library--------------------------------------------
set(HINTS_HBEAM_INC "/usr/local/include/;${HBEAM_ROOT}/;")
set(HINTS_HBEAM_LIB "/usr/local/lib/;${HBEAM_ROOT}")

find_path(HBEAM_INC mwa_hyperbeam.h HINTS ${HINTS_HBEAM_INC})
find_library(HBEAM_LIB mwa_hyperbeam HINTS ${HINTS_HBEAM_LIB} REQUIRED)

if(HBEAM_LIB)
  message(STATUS "HBEAM include path: ${HBEAM_INC}")
  message(STATUS "HBEAM library path: ${HBEAM_LIB}")
endif()

##Look for EveryBeam------------------------------------------------------------
##Finding of EveryBeam and casacore shamelessly stolen from WSClean CMake
##if it's not there; set a flag so we can not include it during compilation

set(HINTS_EVERYBEAM "/usr/local/lib/;${EBEAM_INSTALL}")

find_package(EveryBeam NO_MODULE QUIET HINTS ${HINTS_EVERYBEAM})

set(USE_EVERYBEAM FALSE)

if(${EVERYBEAM_FOUND})
    message(STATUS "EveryBeam library found.")

    ## Only look for casacore if we have EveryBeam
    set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_casacore)
    set(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)
    find_package(Casacore REQUIRED COMPONENTS casa ms tables measures fits)

    if(${CASACORE_FOUND})
        # include_directories(${CASACORE_INCLUDE_DIRS})
        message(STATUS "Casacore library found.")

        ## Now we have to find some aocommon headers because of course we do
        set(HINTS_AOCOM "${EBEAM_ROOT}/external/aocommon/include/aocommon/")
        message(STATUS "Looking for aocommon in ${HINTS_AOCOM}")
        find_path(AOCOM_INC matrix2x2.h HINTS ${HINTS_AOCOM})

        if(AOCOM_INC)
            message(STATUS "AOCOMMON include path found: ${AOCOM_INC}")
            get_filename_component(AOCOM_BASE ${AOCOM_INC} DIRECTORY)
            message(STATUS "AOCOM_BASE: ${AOCOM_BASE}")
            set(USE_EVERYBEAM TRUE)
        else(AOCOM_INC)
            message(
              STATUS
                "Path to aocommon include not found. WODEN needs it to use the EveryBeam library. Building without EveryBeam."
            )
        endif(AOCOM_INC)

        ##We want to directly call the EveryBeam MWA function with az,za
        ##as an input. This means we need to include the mwabeam headers
        ##which aren't installed by default
        set(HINTS_EBEAM_MWA "${EBEAM_ROOT}/cpp/mwabeam/")
        message(STATUS "Looking for beam2016implementation.h in ${HINTS_EBEAM_MWA}")
        find_path(EBEAM_MWA_DIR beam2016implementation.h HINTS ${HINTS_EBEAM_MWA})

        if(EBEAM_MWA_DIR)
            message(STATUS "EveryBeam MWA source code dir found: ${EBEAM_MWA_DIR}")
            # get_filename_component(EBEAM_MWA_DIR ${EBEAM_MWA_DIR} DIRECTORY)
            message(STATUS "EBEAM_MWA_DIR: ${EBEAM_MWA_DIR}")
            set(USE_EVERYBEAM TRUE)
        else(EBEAM_MWA_DIR)
            message(
              STATUS
                "EveryBeam MWA source code dir not found. WODEN needs it to use the EveryBeam library. Building without EveryBeam."
            )
            set(USE_EVERYBEAM FALSE)
        endif(EBEAM_MWA_DIR)


    else(${CASACORE_FOUND})
      message(
        STATUS
          "Casacore not found. WODEN needs to link to it to use the EveryBeam library. Building without EveryBeam."
      )
      set(CASACORE_LIBRARIES "")
    endif(${CASACORE_FOUND})


else(${EVERYBEAM_FOUND})
  message(
    STATUS
      "EveryBeam library not found: WODEN will not be compiled against the EveryBeam library. So no LOFAR for you."
  )
  set(EVERYBEAM_LIB "")
endif(${EVERYBEAM_FOUND})



##Set the linkers for the C code
set(CC_LINKLIBS "-lm -pthread ${HDF5_LIB} ${HBEAM_LIB}")

##All the WODEN header dirs are in include dir
include_directories("${CMAKE_SOURCE_DIR}/include" ${HBEAM_INC})

##Choose whether we are created the main executable, or making the unit tests
set(TARGET_GROUP production CACHE STRING "Group to build")

##If making production code
if(TARGET_GROUP STREQUAL production OR TARGET_GROUP STREQUAL debug)

    if(USE_EVERYBEAM)
        add_library(use_everybeam SHARED
            ${CXXSOURCES}
            ${CMAKE_SOURCE_DIR}/src/logger.c
            ${EBEAM_MWA_DIR}/beam2016implementation.cc
            )
        target_link_libraries(use_everybeam ${EVERYBEAM_LIB} ${CASACORE_LIBRARIES})

        target_include_directories(use_everybeam PRIVATE
                        ${EVERYBEAM_INCLUDE_DIRS}
                        ${AOCOM_BASE}
                        ${CASACORE_INCLUDE_DIRS}
                        ${EBEAM_MWA_DIR})

        ##Need the -march=native to tell aocommon to use the correct 
        ##matrix class, an internal data type to EveryBeam
        if(TARGET_GROUP STREQUAL production)
            target_compile_options(use_everybeam PRIVATE -DHAVE_EVERYBEAM -O3 -Wall
                                                -march=native)
        else()
            target_compile_options(use_everybeam PRIVATE -DHAVE_EVERYBEAM -O0 -g -Wall
                                                -march=native)
        endif()
        

    endif(USE_EVERYBEAM)

    ##Do it for the float code first--------------------------------------------
    ##Generate library for the CUDA code
    add_library(wodenGPU_float SHARED ${GPUSOURCES})
    ##Compile CUDA code with all warnings
    if(USE_CUDA)
        if(TARGET_GROUP STREQUAL production)
            target_compile_options(wodenGPU_float PRIVATE --compiler-options -Wall)
        else()
            target_compile_options(wodenGPU_float PRIVATE -G --compiler-options -g -Wall)
        endif()
    elseif(USE_HIP)
        if(TARGET_GROUP STREQUAL production)
            target_compile_options(wodenGPU_float PRIVATE -Wall)
        else()
            target_compile_options(wodenGPU_float PRIVATE -G -g -Wall)
        endif()
    endif()

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_float SHARED src/woden.c ${CSOURCES})
    target_link_libraries(woden_float PUBLIC wodenGPU_float ${CC_LINKLIBS})
    if(TARGET_GROUP STREQUAL production)
        target_compile_options(woden_float PRIVATE -O3 -Wall)
    else()
        target_compile_options(woden_float PRIVATE -O0 -g -Wall)
    endif()

    if(USE_EVERYBEAM)
        message(INFO " Linking EveryBeam to woden_float")
        target_link_libraries(woden_float PUBLIC use_everybeam)
        target_compile_options(woden_float PRIVATE -DHAVE_EVERYBEAM)
    endif(USE_EVERYBEAM)

    ##Copy things into the wodenpy dir, which means everything is found when
    ##you do a pip install
    add_custom_command(TARGET woden_float POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libwoden_float.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libwoden_float.so")

    ##Now do it for the double code by adding the -DDOUBLE_PRECISION compiler
    ##flag----------------------------------------------------------------------

    add_library(wodenGPU_double SHARED ${GPUSOURCES})
    target_compile_options(wodenGPU_double PRIVATE -DDOUBLE_PRECISION)
    ##Compile CUDA code with all warnings
    if(USE_CUDA)
        if(TARGET_GROUP STREQUAL production)
            target_compile_options(wodenGPU_double PRIVATE --compiler-options -Wall)
        else()
            target_compile_options(wodenGPU_double PRIVATE -G --compiler-options -g -Wall)
        endif()
    elseif(USE_HIP)
        if(TARGET_GROUP STREQUAL production)
            target_compile_options(wodenGPU_double PRIVATE -Wall)
        else()
            target_compile_options(wodenGPU_double PRIVATE -G -g -Wall)
        endif()
    endif()

    ##Add the total WODEN library, linking in the CUDA code
    add_library(woden_double SHARED src/woden.c ${CSOURCES})
    target_compile_options(woden_double PRIVATE -DDOUBLE_PRECISION)

    target_link_libraries(woden_double PUBLIC wodenGPU_double ${CC_LINKLIBS})
    if(TARGET_GROUP STREQUAL production)
        target_compile_options(woden_double PRIVATE -O3 -Wall)
    else()
        target_compile_options(woden_double PRIVATE -O0 -g -Wall)
    endif()

    if(USE_EVERYBEAM)
        message(INFO " Linking EveryBeam to woden_double")
        target_link_libraries(woden_double PUBLIC use_everybeam)
        target_compile_options(woden_double PRIVATE -DHAVE_EVERYBEAM)
    endif(USE_EVERYBEAM)

    ##Copy things into the wodenpy dir, which means everything is found when
    ##you do a pip install
    add_custom_command(TARGET woden_double POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libwoden_double.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libwoden_double.so")

    if(USE_EVERYBEAM)
        message(STATUS "Copying EveryBeam library to wodenpy directory")
    else(USE_EVERYBEAM)
        message(STATUS "EveryBeam not found, so creating a dummy library for wodenpy. This is just a text file so setup.py doesn't complain")
        file(WRITE "${CMAKE_BINARY_DIR}/libuse_everybeam.so" "False\n")
    endif(USE_EVERYBEAM)

    add_custom_command(TARGET woden_double POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy "libuse_everybeam.so"
                  "${CMAKE_CURRENT_LIST_DIR}/wodenpy/libuse_everybeam.so")

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
