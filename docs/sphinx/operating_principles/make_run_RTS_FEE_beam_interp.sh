WODEN_DIR=/home/jline/software/WODEN_new_freq_interp/build/


gcc -O3 -Wall -fPIC -c -DDOUBLE_PRECISION run_RTS_FEE_beam_interp.c \
  -I${WODEN_DIR}/../include/ \
  -I/usr/include/hdf5/serial


g++ -o run_RTS_FEE_beam_interp run_RTS_FEE_beam_interp.o \
  -L${WODEN_DIR} -lwodenCUDA_double -lwodenC_double \
  -I${WODEN_DIR}/../include/ \
  -lerfa -ljson-c -lpal -lm \
  -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5

./run_RTS_FEE_beam_interp

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WODEN_DIR
