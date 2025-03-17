#include <stdio.h>

#define exit mock_exit

void mock_exit(int status) {
    printf("Mock exit called with status %d\n", status);
}


#include "gpu_macros.h"

extern "C" void try_to_memcpy_to_null_pointer_gpu(double *cpu_array, int num_elements) {

  double *gpu_array = NULL;

  // gpuMalloc( (void**)&d_components->ras, num_coords*sizeof(double) );
  gpuMemcpy(gpu_array, cpu_array, num_elements*sizeof(double), gpuMemcpyHostToDevice );

  

}

