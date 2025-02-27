cp ../../build/cmake_testing/GPU_or_C_code/calculate_visibilities/profile_lofar_everybeam .

valgrind --tool=callgrind profile_lofar_everybeam

kcachegrind callgrind.out.*