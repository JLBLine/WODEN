# cp ../../build/cmake_testing/GPU_or_C_code/calculate_visibilities/profile_lofar_everybeam .

rm callgrind.out.*

start_dir=$(pwd)

cd ../../build/cmake_testing/GPU_or_C_code/calculate_visibilities

valgrind --tool=callgrind ./profile_lofar_everybeam

mv callgrind.out.* $start_dir

cd $start_dir

kcachegrind callgrind.out.*