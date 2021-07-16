
gcc -fPIC -c run_Gaussian_beam.c \
  -I$WODEN_DIR/../include/ \
  -I/usr/include/hdf5/serial


g++ -o run_Gaussian_beam run_Gaussian_beam.o \
  -L$WODEN_DIR -lwodenCUDA -lwodenC \
  -I$WODEN_DIR/../include/ \
  -lerfa -ljson-c -lpal \
  -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WODEN_DIR
