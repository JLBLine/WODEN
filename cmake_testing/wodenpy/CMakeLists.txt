##Set an environment variable for the python testing
##Means we don't have to pip install everytime we change python code
# set("WODENPY_TESTING" ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

add_subdirectory(array_layout)
add_subdirectory(observational)
add_subdirectory(phase_rotate)
add_subdirectory(skymodel)
add_subdirectory(uvfits)
add_subdirectory(use_libwoden)
add_subdirectory(wodenpy_setup)
