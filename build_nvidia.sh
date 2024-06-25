mkdir build
cd build
cmake .. -DUSE_CUDA=ON -DHBEAM_INC=/home/msok/github/mwa_hyperbeam_gpu_compilation/include/ -DHBEAM_LIB=/home/msok/github/mwa_hyperbeam_gpu_compilation/target/release/libmwa_hyperbeam.so

make VERBOSE=1



