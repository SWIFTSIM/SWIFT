/*******************************************************************************
 * This file contains functions used to setup and execute GPU tasks from within
 *runner_main.c. Consider this a translator allowing .cu based functions to be
 *called from within runner_main.c
 ******************************************************************************/
#ifdef WITH_CUDA
#ifndef static
#define static
#endif
//#ifndef restrict
//#define restrict __restrict__
//#endif
#endif

/* Required header files */
#include <stdio.h>
/*ifdef __cplusplus prevents name mangling. C code sees exact names
 of functions rather than mangled template names produced by C++*/
#ifdef __cplusplus
extern "C" {
#endif
#include "cell_gpu.h"
#include "cuda_headers.h"
#ifdef __cplusplus
}
#endif

/* function to initialise and printout GPU name*/
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
  // cuda
}
#ifdef __cplusplus
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void CPU_runner_doself1_branch_gradient(struct cell_gpu *restrict ci_gpu) {
  int id = ci_gpu->hydro.parts[0].id;
  printf("id of first part %d\n", id);
  // Do stuff here for interactions on CPU but using the temporary GPU arrays
  //	const int count_i = ci_gpu->hydro.count;
  //  	const int count_j = cj_gpu->hydro.count;
  //	system("pause");
  /* Anything to do here? */
  //  	if (!cell_is_active_hydro(ci_gpu, e) && !cell_is_active_hydro(cj_gpu,
  //  e)) return;
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void GPU_runner_doself1_branch_gradient(struct cell_gpu *restrict ci_gpu) {
  int count = ci_gpu->hydro.count;
  int numBlocks = (count + BLOCK_SIZE - 1) / BLOCK_SIZE;

  struct cell_gpu *d_ci_gpu;
  cudaMalloc((void **)&d_ci_gpu, sizeof(cell_gpu));

  cudaMemcpy(d_ci_gpu, ci_gpu, sizeof(cell_gpu), cudaMemcpyHostToDevice);
  SPH_Sum_Self<<<numBlocks, BLOCK_SIZE>>>(d_ci_gpu);
  cudaMemcpy(ci_gpu, d_ci_gpu, sizeof(cell_gpu), cudaMemcpyDeviceToHost);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void SPH_Sum_Self(cell_gpu *d_ci_gpu) {
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int i = index;
  float sumLoc, xi, yi, zi;
  struct part_gpu *restrict parts = d_ci_gpu->hydro.parts;
  xi = parts[i].x[0];
  yi = parts[i].x[1];
  zi = parts[i].x[2];
  sumLoc = 0.f;
  float h = parts[i].h, mass = parts[i].mass, rho = parts[i].rho;
  const int count = d_ci_gpu->hydro.count;
  //__shared__ float sh_x[BLOCK_SIZE], sh_y[BLOCK_SIZE];
  // copy neighbour particles data to shared memory
  // for (unsigned int j1=0; j1<n; j1+=BLOCK_SIZE){
  // float xj=
  //}
  for (int j = 0; j < count; j++) {
    float xj = parts[j].x[0], yj = parts[j].x[1], zj = parts[j].x[2];
    float rad = sqrt((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi));
    float q = rad / h;
    float q4 = 1.f - 0.5f * q;
    q4 = q4 * q4 * q4 * q4;
    float w = q4 * (2.0f * q + 1.0f);
    float v = mass / rho;
    if (q < 2.0f)
      sumLoc = sumLoc + w * v * 7.0 * 7.0 / (4.0 * 22.0 * h * h);
  }
  // d_Particles[i].ker_sum=sumLoc;
}
#ifdef WITH_CUDA
}
#endif
