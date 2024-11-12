#ifdef WITH_CUDA
#ifndef static
#define static
#endif
#ifndef restrict
#define restrict __restrict__
#endif
#endif

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "cuda_headers.h"
#ifdef __cplusplus
}
#endif

extern "C" {
void print_something_cu() { printf("In Here\n"); }
}

#ifdef __cplusplus
extern "C" {
#endif
void Initialise_GPU() {
  int devId = 0;
  // find and print device name
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  cudaSetDevice(devId);
}
#ifdef __cplusplus
}
#endif
