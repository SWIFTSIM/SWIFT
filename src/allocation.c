#include "allocation.h"

int swift_alloc(void ** memptr, size_t alignment, size_t size){

#ifdef WITH_CUDA

  int result = cudaMallocHost( memptr, size );
  if(result != cudaSuccess)
    return 1;
  else
    return 0;

#else

  int result = posix_memalign(memptr, alignment, size);

#endif
return result;
}

void swift_free(void *memptr){
#ifdef WITH_CUDA
cudaFreeHost(memptr);
#else
free(memptr);
#endif
}

