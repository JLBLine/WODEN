#ifndef CUBE_H
#define CUBE_H

// maximum number of kernels CUBE can track
#define CUBE_KERNEL_MAX 100
#define ulonglong unsigned long long int

#if defined(__cplusplus)
extern "C" {
#endif

  extern int CUBE_nKernel;
  extern ulonglong CUBE_Flops[CUBE_KERNEL_MAX];
  extern ulonglong CUBE_Bytes[CUBE_KERNEL_MAX];
  extern ulonglong CUBE_Hits[CUBE_KERNEL_MAX];
  extern ulonglong CUBE_Misses[CUBE_KERNEL_MAX];
  extern float CUBE_Times[CUBE_KERNEL_MAX];
  extern char *CUBE_Names[CUBE_KERNEL_MAX];

  void CUBE_Init(void);
  void CUBE_Write_Flops(void);
  void CUBE_Write_Benchmark(void);
  int CUBE_Get_Index(const char *);

#if defined(__cplusplus)
}
#endif

// All the C-preprocessor "magic" is below this

#define CUBE_FLOPS 1000
#define CUBE_TIME  2000

#if defined(CUBE_COUNT_MODE)
#define CUBE_MODE CUBE_FLOPS // enable for intial flops counting
#elif defined(CUBE_TIME_MODE)
#define CUBE_MODE CUBE_TIME  // enable for final benchmark
#endif

// Counting macro definitions

#define CUBE_KERNEL_CALL_COUNT(kernel, grid, threads, shared, ...)	\
  {									\
    ulonglong *d_CUBE_Flops = 0;						\
    ulonglong *d_CUBE_Bytes = 0;						\
    ulonglong *d_CUBE_Hits = 0;						\
    ulonglong *d_CUBE_Misses = 0;						\
    CUDA_SAFE_CALL( cudaMalloc( (void**)&d_CUBE_Flops, sizeof(ulonglong)) ); \
    CUDA_SAFE_CALL( cudaMemset( d_CUBE_Flops, 0, sizeof(ulonglong)) );	\
    CUDA_SAFE_CALL( cudaMalloc( (void**)&d_CUBE_Bytes, sizeof(ulonglong)) ); \
    CUDA_SAFE_CALL( cudaMemset( d_CUBE_Bytes, 0, sizeof(ulonglong)) );	\
    CUDA_SAFE_CALL( cudaMalloc( (void**)&d_CUBE_Hits, sizeof(ulonglong)) ); \
    CUDA_SAFE_CALL( cudaMemset( d_CUBE_Hits, 0, sizeof(ulonglong)) );	\
    CUDA_SAFE_CALL( cudaMalloc( (void**)&d_CUBE_Misses, sizeof(ulonglong)) ); \
    CUDA_SAFE_CALL( cudaMemset( d_CUBE_Misses, 0, sizeof(ulonglong)) );	\
									\
    kernel <<< grid , threads, shared >>> ( __VA_ARGS__ , d_CUBE_Flops, d_CUBE_Bytes, \
					    d_CUBE_Hits, d_CUBE_Misses);	\
									\
    ulonglong h_CUBE_Flops, h_CUBE_Bytes, h_CUBE_Hits, h_CUBE_Misses;	\
    CUDA_SAFE_CALL( cudaMemcpy(&h_CUBE_Flops, d_CUBE_Flops, sizeof(ulonglong), cudaMemcpyDeviceToHost) ); \
    CUDA_SAFE_CALL( cudaMemcpy(&h_CUBE_Bytes, d_CUBE_Bytes, sizeof(ulonglong), cudaMemcpyDeviceToHost) ); \
    CUDA_SAFE_CALL( cudaMemcpy(&h_CUBE_Hits, d_CUBE_Hits, sizeof(ulonglong), cudaMemcpyDeviceToHost) ); \
    CUDA_SAFE_CALL( cudaMemcpy(&h_CUBE_Misses, d_CUBE_Misses, sizeof(ulonglong), cudaMemcpyDeviceToHost) ); \
    CUDA_SAFE_CALL( cudaFree( d_CUBE_Flops ) );				\
    CUDA_SAFE_CALL( cudaFree( d_CUBE_Bytes ) );				\
    CUDA_SAFE_CALL( cudaFree( d_CUBE_Hits ) );				\
    CUDA_SAFE_CALL( cudaFree( d_CUBE_Misses ) );				\
    int index = CUBE_Get_Index(#kernel);				\
    CUBE_Flops[index] += h_CUBE_Flops;					\
    CUBE_Bytes[index] += h_CUBE_Bytes;					\
    CUBE_Hits[index] += h_CUBE_Hits;					\
    CUBE_Misses[index] += h_CUBE_Misses;					\
  }

#define CUBE_KERNEL_COUNT(kernel, ...)					\
  __global__ void kernel(__VA_ARGS__, ulonglong *CUBE_total_flops, ulonglong *CUBE_total_bytes,	\
			 ulonglong *CUBE_total_hits, ulonglong *CUBE_total_misses)

#define CUBE_DEVICE_CALL_COUNT(function, ...)	\
  function(__VA_ARGS__, CUBE_flops, CUBE_bytes, CUBE_hits, CUBE_misses)

#define CUBE_DEVICE_COUNT(rtn_type, function, ...)			\
  __device__ rtn_type function(__VA_ARGS__, ulonglong &CUBE_flops, ulonglong &CUBE_bytes,\
			       ulonglong &CUBE_hits, ulonglong &CUBE_misses)


// Timing macro definitions

#define CUBE_KERNEL_CALL_TIME(kernel, grid, threads, shared, ...)	\
  {									\
    cudaEvent_t start, end;						\
    cudaEventCreate(&start);						\
    cudaEventCreate(&end);						\
    cudaEventSynchronize(start);					\
    cudaEventRecord(start, 0);						\
    kernel <<< grid , threads, shared >>> ( __VA_ARGS__ );		\
    cudaEventRecord(end, 0);						\
    cudaEventSynchronize(end);						\
    float runTime;							\
    cudaEventElapsedTime(&runTime, start, end);				\
    cudaEventDestroy(start);						\
    cudaEventDestroy(end);						\
    int index = CUBE_Get_Index(#kernel);				\
    CUBE_Times[index] += runTime;					\
  }
  



// Default macro definitions

#define CUBE_KERNEL_CALL_DEFAULT(kernel, grid, threads, shared, ...)	\
  kernel <<< grid , threads, shared >>>(__VA_ARGS__);

#define CUBE_KERNEL_DEFAULT(kernel, ...)	\
  __global__ void kernel(__VA_ARGS__)

#define CUBE_DEVICE_CALL_DEFAULT(function, ...)	\
  function(__VA_ARGS__)

#define CUBE_DEVICE_DEFAULT(rtn_type, function, ...)	\
  __device__ rtn_type function(__VA_ARGS__)


// set macros according to mode

#if (CUBE_MODE == CUBE_FLOPS) // Count flops and bytes

#define CUBE_KERNEL_CALL CUBE_KERNEL_CALL_COUNT
#define CUBE_KERNEL CUBE_KERNEL_COUNT
#define CUBE_DEVICE_CALL CUBE_DEVICE_CALL_COUNT
#define CUBE_DEVICE CUBE_DEVICE_COUNT

// create flop and byte counters

#define CUBE_START				\
  ulonglong CUBE_flops = 0;			\
  ulonglong CUBE_bytes = 0;			\
  ulonglong CUBE_hits = 0;			\
  ulonglong CUBE_misses = 0;

#define CUBE_ADD_FLOPS(x) CUBE_flops += x;
#define CUBE_ADD_BYTES(x) CUBE_bytes += x;
#define CUBE_ADD_HIT CUBE_hits++;
#define CUBE_ADD_MISS CUBE_misses++;

#define CUBE_END				\
  atomicAdd(CUBE_total_flops, CUBE_flops);	\
  atomicAdd(CUBE_total_bytes, CUBE_bytes);	\
  atomicAdd(CUBE_total_hits, CUBE_hits);	\
  atomicAdd(CUBE_total_misses, CUBE_misses);

#define CUBE_INIT() CUBE_Init()
#define CUBE_WRITE() CUBE_Write_Flops()

#elif (CUBE_MODE == CUBE_TIME) // Measure time and calculate Gflop/s and GiByte/s

#define CUBE_KERNEL_CALL CUBE_KERNEL_CALL_TIME
#define CUBE_KERNEL CUBE_KERNEL_DEFAULT
#define CUBE_DEVICE_CALL CUBE_DEVICE_CALL_DEFAULT
#define CUBE_DEVICE CUBE_DEVICE_DEFAULT

#define CUBE_START
#define CUBE_ADD_FLOPS(x)
#define CUBE_ADD_BYTES(x)
#define CUBE_ADD_HIT
#define CUBE_ADD_MISS
#define CUBE_END

#define CUBE_INIT() CUBE_Init()
#define CUBE_WRITE() CUBE_Write_Benchmark()

#else // Default behaviour (do nothing)

#define CUBE_KERNEL_CALL CUBE_KERNEL_CALL_DEFAULT
#define CUBE_KERNEL CUBE_KERNEL_DEFAULT
#define CUBE_DEVICE_CALL CUBE_DEVICE_CALL_DEFAULT
#define CUBE_DEVICE CUBE_DEVICE_DEFAULT

#define CUBE_START
#define CUBE_ADD_FLOPS(x)
#define CUBE_ADD_BYTES(x)
#define CUBE_ADD_HIT
#define CUBE_ADD_MISS
#define CUBE_END

#define CUBE_INIT() 
#define CUBE_WRITE() 

#endif // end CUBE_MODE

#endif // CUBE
