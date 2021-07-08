cmake_minimum_required(VERSION 3.10)

##Need to find CUDA compiler
include(CheckLanguage)
check_language(CUDA)

project(woden VERSION 0.1 LANGUAGES CUDA CXX C)

##Gotta find CUDA and get it to find libraries
find_package(CUDA)
if(CUDA_FOUND)
  message(STATUS "Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")
  message(STATUS "CUDA library found; building GPU-enabled executables.")
  ##Sets the -std=c++11 nvcc compiler flag I think
  set(CMAKE_CUDA_STANDARD 11)
  #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 ")
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

##Set the -dc flag lets you pull device functions into separate files
##However, doing so actually causes a huge slowdown
##Leaving here in case things change in the future and we need this flag again
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -rdc=true")
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas -O3 -Xcompiler -O3")

##Find the source codes
file(GLOB CSOURCES "src/*.c")
file(GLOB CUDASOURCES "src/*.cu")

##Find the following libraries


set(HINTS_JSONC_LIB "${JSONC_ROOT}/lib;${JSONC_ROOT}/lib64;${JSONC_ROOT};/usr/local/lib/")
set(HINTS_JSONC_INC "${JSONC_ROOT}/include;${JSONC_ROOT}/include/json-c;${JSONC_ROOT};/usr/include/json-c/;/usr/local/include/json-c/")

find_path(JSONC_INC json.h HINTS ${HINTS_JSONC_INC})
find_library(JSONC_LIB json-c HINTS ${HINTS_JSONC_LIB} REQUIRED)
if(JSONC_LIB)
  message(STATUS "JSONC includes path: ${JSONC_INC}")
  message(STATUS "JSONC library path: ${JSONC_LIB}")
endif()

set(HINTS_ERFA_INC "/usr/local/include/;${ERFA_ROOT}/include/;${ERFA_ROOT}")
set(HINTS_ERFA_LIB "/usr/local/lib/;${ERFA_ROOT}/src/.libs/;${ERFA_ROOT}.libs/;${ERFA_ROOT};${ERFA_ROOT}/lib/")


find_path(ERFA_INC erfa.h HINTS ${HINTS_ERFA_INC})
find_library(ERFA_LIB erfa HINTS ${HINTS_ERFA_LIB} REQUIRED)
if(ERFA_LIB)
  message(STATUS "ERFA includes path: ${ERFA_INC}")
  message(STATUS "ERFA library path: ${ERFA_LIB}")
endif()


##Try and find the hdf5 library
set(HINTS_HDF5_INC "/usr/include/hdf5/serial/;${HDF5_ROOT}/include/")
set(HINTS_HDF5_LIB "/usr/lib/x86_64-linux-gnu/hdf5/serial/;${HDF5_ROOT}/lib/")

find_path(HDF5_INC hdf5.h HINTS ${HINTS_HDF5_INC})
find_library(HDF5_LIB hdf5 HINTS ${HINTS_HDF5_LIB} REQUIRED)
if(HDF5_LIB)
  message(STATUS "HDF5 include path: ${HDF5_INC}")
  message(STATUS "HDF5 library path: ${HDF5_LIB}")
endif()

##Try and find the pal library
set(HINTS_PAL_INC "${PAL_ROOT}/include/;/usr/local/include/star/")
set(HINTS_PAL_LIB "/usr/local/lib/;${PAL_ROOT}/lib/")

find_path(PAL_INC pal.h HINTS ${HINTS_PAL_INC})
find_library(PAL_LIB pal HINTS ${HINTS_PAL_LIB} REQUIRED)
if(PAL_LIB)
  message(STATUS "PAL include path: ${PAL_INC}")
  message(STATUS "PAL library path: ${PAL_LIB}")
endif()

##Set the linkers for the C code
set(CC_LINKLIBS "-lm -pthread ${JSONC_LIB} ${ERFA_LIB} ${HDF5_LIB} ${PAL_LIB}")

##All the WODEN header dirs are in include dir
include_directories("${CMAKE_SOURCE_DIR}/include" ${JSONC_INC}
                     ${ERFA_INC} ${HDF5_INC} ${PAL_INC})

##Choose whether we are created the main executable, or making the unit tests
set(TARGET_GROUP production CACHE STRING "Group to build")

# SET(GCC_COVERAGE_COMPILE_FLAGS "-Wall")
# add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

##If making production code
if(TARGET_GROUP STREQUAL production)
    ##Generate library from the C
    add_library(wodenC SHARED ${CSOURCES})
    ##Compile C code with all warnings and optimise
    target_compile_options(wodenC PRIVATE -O3 -Wall)

    ##Generate library for the CUDA code
    add_library(wodenCUDA SHARED ${CUDASOURCES})
    ##Compile CUDA code with all warnings
    target_compile_options(wodenCUDA PRIVATE --compiler-options -Wall)
    ## For CMake versions > 3.18, the arch flag settings are handled by
    ## CUDA_ARCHITECTURES. For now, just switch this warning off - in future,
    ## replace the CUDA_FLAG lines above with something like
    set_property(TARGET wodenCUDA PROPERTY CUDA_ARCHITECTURES OFF)

    ##make the final WODEN executable, linking the libraries
    add_executable(woden src/woden.c)
    target_link_libraries(woden PUBLIC wodenC wodenCUDA ${CC_LINKLIBS})
    target_compile_options(woden PRIVATE -O3 -Wall)


    ##Make a symlink to run_woden.py so it's easy to find
    set(link_woden "run_woden.py")
    find_path(WODENPY run_woden.py HINTS "${CMAKE_CURRENT_LIST_DIR}/src/")
    if(WODENPY)
      message(STATUS "WODENPY path: ${WODENPY}")
    endif()
    set(target_woden "${WODENPY}/run_woden.py")
    ADD_CUSTOM_TARGET(target_woden ALL
                      COMMAND ${CMAKE_COMMAND} -E create_symlink ${target_woden} ${link_woden})

    ##Make a symlink to convert_WSClean_list_to_WODEN.py so it's easy to find
    set(link_convert "convert_WSClean_list_to_WODEN.py")
    find_path(CONVERTWSCLEAN convert_WSClean_list_to_WODEN.py HINTS "${CMAKE_CURRENT_LIST_DIR}/src/")
    if(CONVERTWSCLEAN)
      message(STATUS "CONVERTWSCLEAN path: ${CONVERTWSCLEAN}")
    endif()
    set(target_convert "${CONVERTWSCLEAN}/convert_WSClean_list_to_WODEN.py")
    ADD_CUSTOM_TARGET(target_convert ALL
                      COMMAND ${CMAKE_COMMAND} -E create_symlink ${target_convert} ${link_convert})

    ##Make a symlink to uv2ms.py so it's easy to find
    set(link_uv2ms "uv2ms.py")
    find_path(UV2MS uv2ms.py HINTS "${CMAKE_CURRENT_LIST_DIR}/src/")
    if(UV2MS)
      message(STATUS "UV2MS path: ${UV2MS}")
    endif()
    set(target_uv2ms "${UV2MS}/uv2ms.py")
    ADD_CUSTOM_TARGET(target_uv2ms ALL
                      COMMAND ${CMAKE_COMMAND} -E create_symlink ${target_uv2ms} ${link_uv2ms})

    ##Copy the variable init script so this file can be put in bashrc
    find_path(INIT_WODEN init_WODEN.sh HINTS "${CMAKE_CURRENT_LIST_DIR}/templates/")
    if(INIT_WODEN)
      message(STATUS "INIT_WODEN path: ${INIT_WODEN}")
    endif()
    configure_file("${INIT_WODEN}/init_WODEN.sh" init_WODEN.sh COPYONLY)

##If making test code
elseif(TARGET_GROUP STREQUAL test)

      include(CTest)

      set(HINTS_UNITY "${UNITY_ROOT};${UNITY_ROOT}/src/;/usr/local/Unity/src/")

      find_path(UNITY_DIR unity.c HINTS ${HINTS_UNITY})

      if(UNITY_DIR)
        message(STATUS "UNITY src path: ${UNITY_DIR}")
      endif()

      add_library(Unity STATIC
          "${UNITY_DIR}/unity.c"
      )

      target_include_directories(Unity PUBLIC
          "${UNITY_DIR}"
      )
      add_subdirectory(cmake_testing)

else()
    message(FATAL_ERROR "Given TARGET_GROUP unknown")
endif()
