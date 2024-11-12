#include "hip/hip_runtime.h"
/*******************************************************************************
 * This file contains functions used to setup and execute GPU tasks from within
 *runner_main.c. Consider this a translator allowing .cu based functions to be
 *called from within runner_main.c
 ******************************************************************************/

/* Hacky method to make c++ compilers not die. */
#ifdef WITH_HIP
#ifndef static
#define static
#endif
#ifndef restrict
#define restrict __restrict__
#endif
#endif

/* Required header files */
#include <stdio.h>
/*ifdef WITH_HIP prevents name mangling. C code sees exact names
 of functions rather than mangled template names produced by C++*/

#ifdef __cplusplus
extern "C" {
#endif

#include "../../config.h"
#include "BLOCK_SIZE.h"
#include "HIP_runner_functions.h"
#include "hip/device_functions.h"
#include "part_gpu.h"

void Initialise_GPU() {
  int devId = 0;
  // find and print device name
  hipDeviceProp_t prop;
  hipGetDeviceProperties(&prop, devId);
  printf("Device : %s\n", prop.name);
  hipSetDevice(devId);
  // cuda
}
#ifdef __cplusplus
}
#endif

__global__ void runner_do_self_density_GPU(
    struct part_soa parts_soa, int *d_task_first_part, int *d_task_last_part,
    int *d_bundle_first_part, int *d_bundle_last_part, float d_a, float d_H,
    int bid, int tid, int count_tasks, int tasksperbundle, int nBlocks_per_task,
    int bundle_first_task, int max_parts, int time_bin_inhibited) {
  extern __shared__ float vars[];
  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
  __shared__ int first_part_in_task_blocks, last_part_in_task_blocks;
  first_part_in_task_blocks = d_task_first_part[task_id],
  last_part_in_task_blocks = d_task_last_part[task_id];
  __syncthreads();
  const int b_first_part = d_bundle_first_part[bid];
  const int pid = threadid + first_part_in_task_blocks;
  const int b_last_part = d_bundle_last_part[bid];

  int ttid = 0;
  int first_part = 0;
  int count = 0;
  int last_part = 0;
  float cellx = 0.0, celly = 0.0, cellz = 0.0;
  float hi = 0.0, hig2 = hi * hi * kernel_gamma2;
  float mi = 0.0;
  float uxi = 0.0;
  float uyi = 0.0;
  float uzi = 0.0;
  float pix = 0.0;
  float piy = 0.0;
  float piz = 0.0;
  float rhoi = 0.0;
  float rho_dhi = 0.0;
  float wcounti = 0.0;
  float wcount_dhi = 0.0;
  float div_vi = 0.0;
  float rot_uxi = 0.0;
  float rot_uyi = 0.0;
  float rot_uzi = 0.0;
  //	if(pid<b_last_part&&pid<last_part_in_task_blocks){
  if (pid < last_part_in_task_blocks) {
    ttid = parts_soa.tid_p[pid];
    first_part = d_task_first_part[ttid];
    last_part = d_task_last_part[ttid];
    count = last_part - first_part;
    cellx = parts_soa.locx[pid], celly = parts_soa.locy[pid],
    cellz = parts_soa.locz[pid];
    hi = parts_soa.h[pid], hig2 = hi * hi * kernel_gamma2;
    mi = parts_soa.mass[pid];
    uxi = parts_soa.ux[pid];
    uyi = parts_soa.uy[pid];
    uzi = parts_soa.uz[pid];
    pix = parts_soa.x_p[pid] - cellx;
    piy = parts_soa.y_p[pid] - celly;
    piz = parts_soa.z_p[pid] - cellz;
  }
  if (threadIdx.x == 0) {
    first_part_tid_0 = first_part;
    last_part_tid_0 = last_part;
  }
  __syncthreads();
  int n_neighbours = 0;
  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float *x_p_tmp = (float *)&vars[0];
  float *y_p_tmp = (float *)&vars[BLOCK_SIZE];
  float *z_p_tmp = (float *)&vars[BLOCK_SIZE * 2];
  float *h_tmp = (float *)&vars[BLOCK_SIZE * 3];
  float *mass_tmp = (float *)&vars[BLOCK_SIZE * 4];
  float *ux_tmp = (float *)&vars[BLOCK_SIZE * 5];
  float *uy_tmp = (float *)&vars[BLOCK_SIZE * 6];
  float *uz_tmp = (float *)&vars[BLOCK_SIZE * 7];
  timebin_t *timebin = (timebin_t *)&vars[BLOCK_SIZE * 8];
  /*Particles copied in blocks to shared memory*/
  for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
       b += BLOCK_SIZE) {
    int j = b + threadIdx.x;
    x_p_tmp[threadIdx.x] = parts_soa.x_p[j];
    y_p_tmp[threadIdx.x] = parts_soa.y_p[j];
    z_p_tmp[threadIdx.x] = parts_soa.z_p[j];
    h_tmp[threadIdx.x] = parts_soa.h[j];
    mass_tmp[threadIdx.x] = parts_soa.mass[j];
    ux_tmp[threadIdx.x] = parts_soa.ux[j];
    uy_tmp[threadIdx.x] = parts_soa.uy[j];
    uz_tmp[threadIdx.x] = parts_soa.uz[j];
    timebin[threadIdx.x] = parts_soa.time_bin[j];
    __syncthreads();
    for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
      j = j_block + b;
//      if ((j != pid) && (j < last_part_in_task_blocks) &&
//          timebin[j_block] != time_bin_inhibited) {
//      if ((j < last_part_in_task_blocks) &&
//    	  timebin[j_block] != time_bin_inhibited) {
      if (j < last_part_in_task_blocks) {
        /* Compute the pairwise distance. */
        const float pjx = x_p_tmp[j_block] - cellx;
        const float pjy = y_p_tmp[j_block] - celly;
        const float pjz = z_p_tmp[j_block] - cellz;
        const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
        const float r2 = xij * xij + yij * yij + zij * zij;
        const float hj = h_tmp[j_block], hjg2 = hj * hj * kernel_gamma2;
        //				if((hi < 0.0001f || hj < 0.0001f || r2 <
        //0.0000001f) && pid < last_part_in_task_blocks){ 					printf("very small
        //value for hi %f or hj %f or r2 %f\n", hi, hj, r2);
        //				}
        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
          const float r = sqrt(r2);
          /* Recover some data */
          const float mj = mass_tmp[j_block];
          /* Get the kernel for hi. */
          if(hi<1.f/128.f)printf("h < dx\n");
          const float h_inv = 1.f / hi;
          const float ui = r * h_inv;
          float wi, wi_dx;

          d_kernel_deval(ui, &wi, &wi_dx);

          rhoi += mj * wi;
          rho_dhi -= mj * (hydro_dimension * wi + ui * wi_dx);

          wcounti += wi;
          wcount_dhi -= (hydro_dimension * wi + ui * wi_dx);

          const float r_inv = 1.f / r;
          const float faci = mj * wi_dx * r_inv;

          /* Compute dv dot r */
          float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
                dvz = uzi - uz_tmp[j_block];
          const float dvdr = dvx * xij + dvy * yij + dvz * zij;

          div_vi -= faci * dvdr;

          /* Compute dv cross r */
          float curlvrx = dvy * zij - dvz * yij;
          float curlvry = dvz * xij - dvx * zij;
          float curlvrz = dvx * yij - dvy * xij;

          rot_uxi += faci * curlvrx;
          rot_uyi += faci * curlvry;
          rot_uzi += faci * curlvrz;
        }
      }
    }
    __syncthreads();
  }
  if (pid < last_part_in_task_blocks) {
	float wi, wi_dx;
	d_kernel_deval(0.f, &wi, &wi_dx);
//	printf("mass i %e, self rho %e sum rho %e\n", mi, mi*wi, rhoi);
    parts_soa.rho[pid] = rhoi, parts_soa.rho_dh[pid] = rho_dhi;
    parts_soa.wcount[pid] = wcounti, parts_soa.wcount_dh[pid] = wcount_dhi;
    parts_soa.div_v[pid] = div_vi;
    parts_soa.rot_ux[pid] = rot_uxi, parts_soa.rot_uy[pid] = rot_uyi,
    parts_soa.rot_uz[pid] = rot_uzi;
  }
}
#ifdef __cplusplus
extern "C" {
#endif
void launch_density_kernel(struct part_soa parts_soa, int *d_task_first_part,
                           int *d_task_last_part, int *d_bundle_first_part,
                           int *d_bundle_last_part, float d_a, float d_H,
                           const char *loop_type, hipStream_t stream, int bid,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y, int tid,
                           int offset, int bundle_first_task, int max_parts,
                           int time_bin_inhibited) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  runner_do_self_density_GPU<<<gridShape, BLOCK_SIZE,
                               8 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(timebin_t),
                               stream>>>(
      parts_soa, d_task_first_part, d_task_last_part, d_bundle_first_part,
      d_bundle_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef __cplusplus
}
#endif
