WODEN_DIR=/home/jline/software/WODEN_new_freq_interp/build/

gcc -O3 -Wall -fPIC -c -DDOUBLE_PRECISION run_RTS_analy_beam.c \
  -I$WODEN_DIR/../include/ \
  -I/usr/include/hdf5/serial


g++ -o run_RTS_analy_beam run_RTS_analy_beam.o \
  -L$WODEN_DIR -lwodenCUDA_double -lwodenC_double \
  -I$WODEN_DIR/../include/ \
  -lerfa -ljson-c -lpal \
  -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5

./run_RTS_analy_beam

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WODEN_DIR
