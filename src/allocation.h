#ifndef SWIFT_ALLOCATION_H
#define SWIFT_ALLOCATION_H

#ifdef WITH_CUDA
  #ifndef __NVCC__
    #include <cuda.h>
    #include <cuda_runtime.h>
  #endif
#endif

#include <stdlib.h>

int swift_alloc(void ** memptr, size_t alignment, size_t size){

#ifdef WITH_CUDA

  int result = cudaMallocHost( memptr, size );
  if(result != cudaSuccess)
    result = 1;
  else
    result = 0;

#else

  int result = posix_memalign(memptr, alignment, size);

#endif
return result;
}

#endif
