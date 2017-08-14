#ifndef SWIFT_ALLOCATION_H
#define SWIFT_ALLOCATION_H

#ifdef WITH_CUDA
  #ifndef __NVCC__
    #include <cuda.h>
    #include <cuda_runtime.h>
  #endif
#endif

#include <stdlib.h>

int swift_alloc(void ** memptr, size_t alignment, size_t size);
void swift_free( void* memptr);
#endif
