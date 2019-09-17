# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = all woden clean

CFITSIO_INC := 
CFITSIO_LIB := 

##Compiling
CC = gcc-6                                                     #cc compiler to use
CC_INCLUDES = -I/usr/local/cuda/include -I/usr/include/json-c \
              -I/usr/local/include  #cc comp includes

NVCC = nvcc                       #nvcc compiler to use
NVCC_COMP_FLAGS = -std=c++11 -ccbin=/usr/bin/gcc-6 -arch=sm_61 #nvcc comp flags
NVCC_INCLUDES = -I/usr/local/cuda/include                      #nvcc comp includes 


##Linking
CC_LINK_FLAGS = -g -Wall                                            ##Linking flags
CC_LINKLIBS = -L/usr/local/cuda/lib64 -lm -lcuda -lcudart -lcufft -ljson-c \
              -L/usr/local/lib -lcfitsio ##Linking libs

##CC sources and objs
CC_SRCS := woden.c read_and_write.c shapelet_basis.c
CC_OBJS := $(CC_SRCS:%.c=%.o)

##CUDA sources and objs
CUDA_SRCS := woden_lib.cu
CUDA_OBJS := $(CUDA_SRCS:%.cu=%.o)

all: woden

woden: 
	@echo "Creating CUDA objects"
	${NVCC} -c ${NVCC_COMP_FLAGS} $(CUDA_SRCS)
	@echo "Creating C objects"
	${CC} -c $(CC_SRCS) ${CC_INCLUDES}
	@echo "Linking dem all"
	${CC} ${CC_LINK_FLAGS} ${CC_OBJS} ${CUDA_OBJS} -o woden ${CC_LINKLIBS}

clean:
	@echo "Cleaning up..."
	rm -rvf *.o woden
