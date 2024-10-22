/*******************************************************************************
 * This file contains functions used to setup and execute GPU tasks from within
 *runner_main.c. Consider this a translator allowing .cu based functions to be
 *called from within runner_main.c
 ******************************************************************************/

/* Hacky method to make c++ compilers not die. */
#ifdef WITH_CUDA
#ifndef static
#define static
#endif
#ifndef restrict
#define restrict __restrict__
#endif
#endif

/* Required header files */
#include <stdio.h>
/*ifdef WITH_CUDA prevents name mangling. C code sees exact names
 of functions rather than mangled template names produced by C++*/
#ifdef WITH_CUDA
extern "C" {
#endif

#include "../../config.h"

#ifndef BLOCK_SIZE_H
#include "BLOCK_SIZE.h"
#endif

#include "GPU_runner_functions.h"
#include "device_functions.h"
#include "part_gpu.h"
#include <cuda_profiler_api.h>

#ifdef WITH_CUDA
}
#endif

/* function to initialise GPU and printout GPU name*/
#ifdef WITH_CUDA
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
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void tester(
    struct part_soa parts_soa, int *d_task_first_part, int *d_task_last_part, float d_a, float d_H,
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
  const int pid = threadid + first_part_in_task_blocks;

  if (pid < last_part_in_task_blocks) {
    parts_soa.tid_p[pid] = 1;
  }
//  if(parts_soa.tid_p[pid] == 1 && pid < last_part_in_task_blocks)
//	  printf("tid %i last_part_in_blocks %i\n", parts_soa.tid_p[pid], last_part_in_task_blocks);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_self_density_GPU(
    struct part_soa parts_soa, int *d_task_first_part, int *d_task_last_part, float d_a, float d_H,
    int count_tasks, int tasksperbundle, int nBlocks_per_task,
    int bundle_first_task, int max_parts) {
  extern __shared__ float vars[];
  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  __shared__ int first_part_in_task_blocks, last_part_in_task_blocks;
  int first_part_in_task_blocks, last_part_in_task_blocks;
  first_part_in_task_blocks = d_task_first_part[task_id],
  last_part_in_task_blocks = d_task_last_part[task_id];
//  __syncthreads();
  const int pid = threadid + first_part_in_task_blocks;

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
  int Found_neighbours = 0;
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
//  if (threadIdx.x == 0) {
//    first_part_tid_0 = first_part;
//    last_part_tid_0 = last_part;
//  }
//  __syncthreads();
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
        //				if((hi < 0.0001f || hj < 0.0001f || r2 <
        //0.0000001f) && pid < last_part_in_task_blocks){ 					printf("very small
        //value for hi %f or hj %f or r2 %f\n", hi, hj, r2);
        //				}
        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
//        if (r2 < hig2 && r2 > (0.01f/256.f)*(0.01f/256.f)) {
          Found_neighbours=1;
          const float r = sqrt(r2);
          /* Recover some data */
          const float mj = mass_tmp[j_block];
          /* Get the kernel for hi. */
          if(hi<1.f/256.f)printf("h < dx\n");
//          if(hi<1.f/256.f)printf("h < dx\n");
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
//	float wi, wi_dx;
//	d_kernel_deval(0.f, &wi, &wi_dx);
//	printf("mass i %e, self rho %e sum rho %e\n", mi, mi*wi, rhoi);
//    if(Found_neighbours == 0) printf("Not sure what's going on but no neighbours found in GPU loop\n");
    parts_soa.rho[pid] = rhoi, parts_soa.rho_dh[pid] = rho_dhi;
    parts_soa.wcount[pid] = wcounti, parts_soa.wcount_dh[pid] = wcount_dhi;
    parts_soa.div_v[pid] = div_vi;
    parts_soa.rot_ux[pid] = rot_uxi, parts_soa.rot_uy[pid] = rot_uyi,
    parts_soa.rot_uz[pid] = rot_uzi;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void DOSELF_GPU_AOS(
    struct part_aos *parts_aos, int *d_task_first_part, int *d_task_last_part, float d_a, float d_H,
    int count_tasks, int tasksperbundle, int nBlocks_per_task,
    int bundle_first_task, int max_parts, double * d_cell_x,
	double * d_cell_y, double * d_cell_z) {
  extern __shared__ float vars[];
  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  __shared__ int first_part_in_task_blocks, last_part_in_task_blocks;
  int first_part_in_task_blocks, last_part_in_task_blocks;
  first_part_in_task_blocks = d_task_first_part[task_id],
  last_part_in_task_blocks = d_task_last_part[task_id];
//  __syncthreads();
  const int pid = threadid + first_part_in_task_blocks;

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
  int Found_neighbours = 0;
  struct part_aos ipart = parts_aos[pid];
  //	if(pid<b_last_part&&pid<last_part_in_task_blocks){
  if (pid < last_part_in_task_blocks) {
    ttid = task_id;
    first_part = d_task_first_part[ttid];
    last_part = d_task_last_part[ttid];
    count = last_part - first_part;
    cellx = d_cell_x[ttid], celly = d_cell_y[ttid],
    cellz = d_cell_z[ttid];
    hi = ipart.h, hig2 = hi * hi * kernel_gamma2;
    mi = ipart.mass;
    uxi = ipart.ux;
    uyi = ipart.uy;
    uzi = ipart.uz;
    pix = ipart.x_p - cellx;
    piy = ipart.y_p - celly;
    piz = ipart.z_p - cellz;
  }
//  if (threadIdx.x == 0) {
//    first_part_tid_0 = first_part;
//    last_part_tid_0 = last_part;
//  }
//  __syncthreads();
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
  int *timebin = (int *)&vars[BLOCK_SIZE * 8];
  /*Particles copied in blocks to shared memory*/
//  struct parts_aos jparts[count];
  for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
       b += BLOCK_SIZE) {
    int j = b + threadIdx.x;
    struct part_aos jpart = parts_aos[j];
    x_p_tmp[threadIdx.x] = jpart.x_p;
    y_p_tmp[threadIdx.x] = jpart.y_p;
    z_p_tmp[threadIdx.x] = jpart.z_p;
    h_tmp[threadIdx.x] = jpart.h;
    mass_tmp[threadIdx.x] = jpart.mass;
    ux_tmp[threadIdx.x] = jpart.ux;
    uy_tmp[threadIdx.x] = jpart.uy;
    uz_tmp[threadIdx.x] = jpart.uz;
    timebin[threadIdx.x] = jpart.time_bin;
    __syncthreads();
    for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
      j = j_block + b;
      if (j < last_part_in_task_blocks) {
        /* Compute the pairwise distance. */
        const float pjx = x_p_tmp[j_block] - cellx;
        const float pjy = y_p_tmp[j_block] - celly;
        const float pjz = z_p_tmp[j_block] - cellz;
        const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
        const float r2 = xij * xij + yij * yij + zij * zij;
        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
          Found_neighbours=1;
          const float r = sqrt(r2);
          /* Recover some data */
          const float mj = mass_tmp[j_block];
          /* Get the kernel for hi. */
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
//	float wi, wi_dx;
//	d_kernel_deval(0.f, &wi, &wi_dx);
//	printf("mass i %e, self rho %e sum rho %e\n", mi, mi*wi, rhoi);
//    if(Found_neighbours == 0) printf("Not sure what's going on but no neighbours found in GPU loop\n");
    parts_aos[pid].rho = rhoi, parts_aos[pid].rho_dh = rho_dhi;
    parts_aos[pid].wcount = wcounti, parts_aos[pid].wcount_dh = wcount_dhi;
    parts_aos[pid].div_v = div_vi;
    parts_aos[pid].rot_ux = rot_uxi, parts_aos[pid].rot_uy = rot_uyi,
    parts_aos[pid].rot_uz = rot_uzi;
  }
}
#ifdef WITH_CUDA
}
#endif

//template <typename T>

#ifdef WITH_CUDA
extern "C" {
#endif
//#include <cuda/barrier>
__global__ void DOSELF_GPU_AOS_F4(
		struct part_aos_f4_send * __restrict__ parts_send, struct part_aos_f4_recv * __restrict__ parts_recv,
		const float d_a, const float d_H,
	const int bundle_first_task, const int2 * __restrict__ d_task_first_part_f4) {

  extern __shared__ float4 vars_f4[];

//  auto group = cooperative_groups::this_thread_block();
  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;
//  cuda::barrier<cuda::thread_scope_system> bar;

  int first_part_in_task_blocks, last_part_in_task_blocks;
  int2 first_last_parts = d_task_first_part_f4[task_id];
  first_part_in_task_blocks = first_last_parts.x;
  last_part_in_task_blocks = first_last_parts.y;

  const int pid = threadid + first_part_in_task_blocks;

  float4 res_rho = {0.0, 0.0, 0.0, 0.0};
  float4 res_rot = {0.0, 0.0, 0.0, 0.0};
  const part_aos_f4_send pi = parts_send[pid];
  const float4 x_pi = pi.x_p_h;
  const float4 ux_pi = pi.ux_m;
  const float hi = x_pi.w, hig2 = hi * hi * kernel_gamma2;
  int n_neighbours = 0;
  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float4 * __restrict__ x_p_h_tmp = (float4 *)&vars_f4[0];
  float4 * __restrict__ ux_m_tmp = (float4 *)&vars_f4[BLOCK_SIZE];
  /*Particles copied in blocks to shared memory*/
  for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
       b += BLOCK_SIZE) {
    int j = b + threadIdx.x;
    struct part_aos_f4_send pj = parts_send[j];
    x_p_h_tmp[threadIdx.x] = pj.x_p_h;
    ux_m_tmp[threadIdx.x] = pj.ux_m;
    __syncthreads();
    for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
      j = j_block + b;
      if (j < last_part_in_task_blocks) {
        /* Compute the pairwise distance. */
        const float4 x_p_h_j = x_p_h_tmp[j_block];
        const float4 ux_m_j = ux_m_tmp[j_block];
        const float xij = x_pi.x - x_p_h_j.x,
        		    yij = x_pi.y - x_p_h_j.y,
        		    zij = x_pi.z - x_p_h_j.z;
        const float r2 = xij * xij + yij * yij + zij * zij;
        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
          const float r = sqrtf(r2);
          /* Recover some data */
          const float mj = ux_m_j.w;
          /* Get the kernel for hi. */
          const float h_inv = 1.f / hi;
          const float ui = r * h_inv;
          float wi, wi_dx;

          d_kernel_deval(ui, &wi, &wi_dx);
          /*Add to sums of rho, rho_dh, wcount and wcount_dh*/
          res_rho.x += mj * wi;
          res_rho.y -= mj * (hydro_dimension * wi + ui * wi_dx);
          res_rho.z += wi;
          res_rho.w -= (hydro_dimension * wi + ui * wi_dx);

          const float r_inv = 1.f / r;
          const float faci = mj * wi_dx * r_inv;

          /* Compute dv dot r */
          const float dvx = ux_pi.x - ux_m_j.x, dvy = ux_pi.y - ux_m_j.y,
                dvz = ux_pi.z - ux_m_j.z;
          const float dvdr = dvx * xij + dvy * yij + dvz * zij;

          /* Compute dv cross r */
          const float curlvrx = dvy * zij - dvz * yij;
          const float curlvry = dvz * xij - dvx * zij;
          const float curlvrz = dvx * yij - dvy * xij;
          /*Add to sums of rot_u and div_v*/
          res_rot.x += faci * curlvrx;
          res_rot.y += faci * curlvry;
          res_rot.z += faci * curlvrz;
          res_rot.w -= faci * dvdr;
        }
      }
    }
    __syncthreads();
  }
  if (pid < last_part_in_task_blocks) {
    parts_recv[pid].rho_dh_wcount = res_rho;
    parts_recv[pid].rot_ux_div_v = res_rot;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_density_aos(struct part_aos *parts_aos, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts,
                           double * d_cell_x,
						   double * d_cell_y, double * d_cell_z) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  DOSELF_GPU_AOS<<<gridShape, BLOCK_SIZE,
                               8 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(int),
                               stream>>>(
      parts_aos, d_task_first_part, d_task_last_part, d_a, d_H, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, max_parts, d_cell_x,
	  d_cell_y, d_cell_z);
//  runner_do_self_density_GPU_naive<<<gridShape, BLOCK_SIZE, 0, stream>>>(
//        parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//        nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
struct first_part{
	int list[32];
};
void launch_density_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
		                   float d_a, float d_H, cudaStream_t stream, int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int2 *d_task_first_part_f4) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  DOSELF_GPU_AOS_F4<<<gridShape, BLOCK_SIZE,
                               2 * BLOCK_SIZE * sizeof(float4), stream>>>(
      parts_send, parts_recv, d_a, d_H, bundle_first_task, d_task_first_part_f4);
//  runner_do_self_density_GPU_naive<<<gridShape, BLOCK_SIZE, 0, stream>>>(
//        parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//        nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void DOSELF_GPU_AOS_G(
    struct part_aos_g *parts_aos, int *d_task_first_part, int *d_task_last_part, float d_a, float d_H,
    int count_tasks, int tasksperbundle, int nBlocks_per_task,
    int bundle_first_task, int max_parts, double * d_cell_x,
	double * d_cell_y, double * d_cell_z) {
  extern __shared__ float varsg[];
  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  __shared__ int first_part_in_task_blocks, last_part_in_task_blocks;
  int first_part_in_task_blocks, last_part_in_task_blocks;
  first_part_in_task_blocks = d_task_first_part[task_id],
  last_part_in_task_blocks = d_task_last_part[task_id];
//  __syncthreads();
  const int pid = threadid + first_part_in_task_blocks;

  int ttid = 0;
  int first_part = 0;
  int count = 0;
  int last_part = 0;
  float cellx = 0.0, celly = 0.0, cellz = 0.0;
  float ci = 0.0, cj = 0.0;
  float hi = 0.0, hig2 = hi * hi * kernel_gamma2;
  float mi = 0.0;
  float uxi = 0.0;
  float uyi = 0.0;
  float uzi = 0.0;
  float pix = 0.0;
  float piy = 0.0;
  float piz = 0.0;
  float rhoi = 0.0;
  float div_vi = 0.0;
  int Found_neighbours = 0;
  float v_sig;
  float u = 0.f;
  float laplace_u = 0.0;
  float alpha_visc_max_ngb = 0.0;
  if (pid < last_part_in_task_blocks) {
    ttid = task_id;
    first_part = d_task_first_part[ttid];
    last_part = d_task_last_part[ttid];
    count = last_part - first_part;
    cellx = d_cell_x[ttid], celly = d_cell_y[ttid],
    cellz = d_cell_z[ttid];
    hi = parts_aos[pid].h, hig2 = hi * hi * kernel_gamma2;
    mi = parts_aos[pid].mass;
    uxi = parts_aos[pid].ux;
    uyi = parts_aos[pid].uy;
    uzi = parts_aos[pid].uz;
    pix = parts_aos[pid].x_p - cellx;
    piy = parts_aos[pid].y_p - celly;
    piz = parts_aos[pid].z_p - cellz;
    ci = parts_aos[pid].soundspeed;
    v_sig = parts_aos[pid].v_sig;
    u = parts_aos[pid].u;
    laplace_u = parts_aos[pid].laplace_u;
    alpha_visc_max_ngb = parts_aos[pid].alpha_visc_max_ngb;
  }
//  if (threadIdx.x == 0) {
//    first_part_tid_0 = first_part;
//    last_part_tid_0 = last_part;
//  }
//  __syncthreads();
  int n_neighbours = 0;
  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float *x_p_tmp = (float *)&varsg[0];
  float *y_p_tmp = (float *)&varsg[BLOCK_SIZE];
  float *z_p_tmp = (float *)&varsg[BLOCK_SIZE * 2];
  float *h_tmp = (float *)&varsg[BLOCK_SIZE * 3];
  float *mass_tmp = (float *)&varsg[BLOCK_SIZE * 4];
  float *ux_tmp = (float *)&varsg[BLOCK_SIZE * 5];
  float *uy_tmp = (float *)&varsg[BLOCK_SIZE * 6];
  float *uz_tmp = (float *)&varsg[BLOCK_SIZE * 7];
  float *cj_tmp = (float *)&varsg[BLOCK_SIZE * 8];
  float *alpha_tmp = (float *)&varsg[BLOCK_SIZE * 9];
  float *u_tmp = (float *)&varsg[BLOCK_SIZE * 10];
  float *rho_tmp = (float *)&varsg[BLOCK_SIZE * 11];
  int *timebin = (int *)&varsg[BLOCK_SIZE * 12];
  /*Particles copied in blocks to shared memory*/
  for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
       b += BLOCK_SIZE) {
    int j = b + threadIdx.x;
    x_p_tmp[threadIdx.x] = parts_aos[j].x_p;
    y_p_tmp[threadIdx.x] = parts_aos[j].y_p;
    z_p_tmp[threadIdx.x] = parts_aos[j].z_p;
    h_tmp[threadIdx.x] = parts_aos[j].h;
    mass_tmp[threadIdx.x] = parts_aos[j].mass;
    ux_tmp[threadIdx.x] = parts_aos[j].ux;
    uy_tmp[threadIdx.x] = parts_aos[j].uy;
    uz_tmp[threadIdx.x] = parts_aos[j].uz;
    timebin[threadIdx.x] = parts_aos[j].time_bin;
    cj_tmp[threadIdx.x] = parts_aos[j].soundspeed;
    alpha_tmp[threadIdx.x] = parts_aos[j].visc_alpha;
    u_tmp[threadIdx.x] = parts_aos[j].u;
    rho_tmp[threadIdx.x] = parts_aos[j].rho;
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
        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
          Found_neighbours=1;
          const float r = sqrt(r2);
          const float r_inv = 1.f / r;
          /* Recover some data */
          const float mj = mass_tmp[j_block];
          /* Get the kernel for hi. */
          const float h_inv = 1.f / hi;
          float wi, wi_dx;
          /* Cosmology terms for the signal velocity */
          const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
          const float a2_Hubble = d_a * d_a * d_H;
          /* Compute dv dot r */
          float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
                dvz = uzi - uz_tmp[j_block];
          const float dvdr = dvx * xij + dvy * yij + dvz * zij;
          /* Add Hubble flow */
          const float dvdr_Hubble = dvdr + a2_Hubble * r2;
          /* Are the particles moving towards each others ? */
          const float omega_ij = min(dvdr_Hubble, 0.f);
          const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

          /* Signal velocity */
          const float new_v_sig = ci + cj_tmp[j_block] - const_viscosity_beta * mu_ij;
          /* Update if we need to */
          v_sig = max(v_sig, new_v_sig);
          /* Calculate Del^2 u for the thermal diffusion coefficient. */
          /* Need to get some kernel values F_ij = wi_dx */
          const float ui = r * h_inv;
          d_kernel_deval(ui, &wi, &wi_dx);

          const float delta_u_factor = (u - u_tmp[j_block]) * r_inv;
          laplace_u += mj * delta_u_factor * wi_dx / rho_tmp[j_block];

          /* Set the maximal alpha from the previous step over the neighbours
           * (this is used to limit the diffusion in hydro_prepare_force) */
          const float alpha_j = alpha_tmp[j_block];
          alpha_visc_max_ngb = max(alpha_visc_max_ngb, alpha_j);

        }
      }
    }
    __syncthreads();
  }
  if (pid < last_part_in_task_blocks) {
    parts_aos[pid].v_sig = v_sig, parts_aos[pid].laplace_u = laplace_u;
    parts_aos[pid].alpha_visc_max_ngb = alpha_visc_max_ngb;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void DOSELF_GPU_AOS_F4_G(
		struct part_aos_f4_g_send * __restrict__ parts_send, struct part_aos_f4_g_recv * __restrict__ parts_recv,
		const float d_a, const float d_H,
	    const int bundle_first_task, const int2 * __restrict__ d_task_first_part_f4) {

  extern __shared__ float4 varsf4_g[];

  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

//  __shared__ int first_part_in_task_blocks, last_part_in_task_blocks;
  int2 first_last_parts = d_task_first_part_f4[task_id];
  int first_part_in_task_blocks = first_last_parts.x;
  int last_part_in_task_blocks = first_last_parts.y;
//  __syncthreads();
  const int pid = threadid + first_part_in_task_blocks;

  /*Keep this*/
  float v_sig = 0.f;
  float alpha_visc_max_ngb = 0.f;
  /////////////

  struct part_aos_f4_g_send pi = parts_send[pid];
  float4 x_h_i = pi.x_h;
  float4 ux_m_i = pi.ux_m;
  float4 rho_avisc_u_c_i = pi.rho_avisc_u_c;
  float3 vsig_lapu_aviscmax_i = {0.f, 0.f, 0.f};

  const float hi = x_h_i.w, hig2 = hi * hi * kernel_gamma2;

  int n_neighbours = 0;
  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float4 * __restrict__ x_h_tmp = (float4 *)&varsf4_g[0];
  float4 * __restrict__ ux_m_tmp = (float4 *)&varsf4_g[BLOCK_SIZE];
  float4 * __restrict__ rho_avisc_u_c_tmp = (float4 *)&varsf4_g[BLOCK_SIZE * 2];

  /*Particles copied in blocks to shared memory*/
  for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
       b += BLOCK_SIZE) {

    int j = b + threadIdx.x;

    struct part_aos_f4_g_send pj = parts_send[j];
    x_h_tmp[threadIdx.x] = pj.x_h;
    ux_m_tmp[threadIdx.x] = pj.ux_m;
    rho_avisc_u_c_tmp[threadIdx.x] = pj.rho_avisc_u_c;

    __syncthreads();
    for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
      j = j_block + b;
      if (j < last_part_in_task_blocks) {
    	float4 x_h_j = x_h_tmp[j_block];
    	float4 ux_m_j = ux_m_tmp[j_block];
    	float4 rho_avisc_u_c_j = rho_avisc_u_c_tmp[j_block];
        /* Compute the pairwise distance. */
        const float xij = x_h_i.x - x_h_j.x, yij = x_h_i.y - x_h_j.y, zij = x_h_i.z - x_h_j.z;
        const float r2 = xij * xij + yij * yij + zij * zij;

        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
          const float r = sqrt(r2);
          const float r_inv = 1.f / r;
          /* Recover some data */
          const float mj = ux_m_j.w;
          /* Get the kernel for hi. */
          const float h_inv = 1.f / hi;
          float wi, wi_dx;
          /* Cosmology terms for the signal velocity */
          const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
          const float a2_Hubble = d_a * d_a * d_H;
          /* Compute dv dot r */
          float dvx = ux_m_i.x - ux_m_j.x, dvy = ux_m_i.y - ux_m_j.y,
                dvz = ux_m_i.z - ux_m_j.z;
          const float dvdr = dvx * xij + dvy * yij + dvz * zij;
          /* Add Hubble flow */
          const float dvdr_Hubble = dvdr + a2_Hubble * r2;
          /* Are the particles moving towards each others ? */
          const float omega_ij = min(dvdr_Hubble, 0.f);
          const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

          /* Signal velocity */
          const float new_v_sig = rho_avisc_u_c_i.w + rho_avisc_u_c_j.w - const_viscosity_beta * mu_ij;
          /* Update if we need to */
          vsig_lapu_aviscmax_i.x = fmaxf(vsig_lapu_aviscmax_i.x, new_v_sig);
          /* Calculate Del^2 u for the thermal diffusion coefficient. */
          /* Need to get some kernel values F_ij = wi_dx */
          const float ui = r * h_inv;
          d_kernel_deval(ui, &wi, &wi_dx);

          const float delta_u_factor = (rho_avisc_u_c_i.z - rho_avisc_u_c_j.z) * r_inv;
          vsig_lapu_aviscmax_i.y += mj * delta_u_factor * wi_dx / rho_avisc_u_c_j.x;

          /* Set the maximal alpha from the previous step over the neighbours
           * (this is used to limit the diffusion in hydro_prepare_force) */
          const float alpha_j = rho_avisc_u_c_j.y;
          vsig_lapu_aviscmax_i.z = fmaxf(vsig_lapu_aviscmax_i.z, alpha_j);

        }
      }
    }
    __syncthreads();
  }
  if (pid < last_part_in_task_blocks) {
//	  printf("v %f lap %f maxvisc %f\n", vsig_lapu_aviscmax_empty_i.x, vsig_lapu_aviscmax_empty_i.y, vsig_lapu_aviscmax_empty_i.z);
    parts_recv[pid].vsig_lapu_aviscmax = vsig_lapu_aviscmax_i;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void DOSELF_GPU_AOS_F(
    struct part_aos_f *parts_aos, int *d_task_first_part, int *d_task_last_part, float d_a, float d_H,
    int count_tasks, int tasksperbundle, int nBlocks_per_task,
    int bundle_first_task, int max_parts, double * d_cell_x,
	double * d_cell_y, double * d_cell_z) {
  extern __shared__ float varsf[];
  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  int first_part_in_task_blocks, last_part_in_task_blocks;
  first_part_in_task_blocks = d_task_first_part[task_id],
  last_part_in_task_blocks = d_task_last_part[task_id];

  const int pid = threadid + first_part_in_task_blocks;

  int ttid = 0;
  int first_part = 0;
  int count = 0;
  int last_part = 0;
  float cellx = 0.0, celly = 0.0, cellz = 0.0;
  float ci = 0.0, cj = 0.0;
  float hi = 0.0, hig2 = 0.0;
  float mi = 0.0;
  float uxi = 0.0;
  float uyi = 0.0;
  float uzi = 0.0;
  float pix = 0.0;
  float piy = 0.0;
  float piz = 0.0;
  float rhoi = 0.0;
  float div_vi = 0.0;
  int Found_neighbours = 0;
  float v_sigi;
  float ui = 0.f;
  float u_dti = 0.f;
  float laplace_ui = 0.0;
  float alpha_visc_max_ngb = 0.0;
  float pressurei = 0.0;
  float alphavisci = 0.0;
  float alphadiffi = 0.0;
  float fi = 0.0;
  float balsarai = 0.0;
  float ahydroxi = 0.0;
  float ahydroyi = 0.0;
  float ahydrozi = 0.0;
  float h_dti = 0.0;
  int min_ngb_time_bin = 0;
  if (pid < last_part_in_task_blocks) {
    ttid = task_id;
    first_part = d_task_first_part[ttid];
    last_part = d_task_last_part[ttid];
    count = last_part - first_part;
    cellx = d_cell_x[ttid], celly = d_cell_y[ttid],
    cellz = d_cell_z[ttid];
    hi = parts_aos[pid].h, hig2 = hi * hi * kernel_gamma2;
    mi = parts_aos[pid].mass;
    uxi = parts_aos[pid].ux;
    uyi = parts_aos[pid].uy;
    uzi = parts_aos[pid].uz;
    pix = parts_aos[pid].x_p - cellx;
    piy = parts_aos[pid].y_p - celly;
    piz = parts_aos[pid].z_p - cellz;
    ci = parts_aos[pid].soundspeed;
    fi = parts_aos[pid].f;
    v_sigi = parts_aos[pid].v_sig;
    ui = parts_aos[pid].u;
    rhoi = parts_aos[pid].rho;
    pressurei = parts_aos[pid].pressure;
    balsarai = parts_aos[pid].balsara;
    alphavisci = parts_aos[pid].alpha_visc;
    alphadiffi = parts_aos[pid].alpha_diff;
    min_ngb_time_bin = parts_aos[pid].min_ngb_time_bin;
//    laplace_u = parts_aos[pid].laplace_u;
//    alpha_visc_max_ngb = parts_aos[pid].alpha_visc_max_ngb;
  }
//  if (threadIdx.x == 0) {
//    first_part_tid_0 = first_part;
//    last_part_tid_0 = last_part;
//  }
//  __syncthreads();
  int n_neighbours = 0;
  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float *x_p_tmp = (float *)&varsf[0];
  float *y_p_tmp = (float *)&varsf[BLOCK_SIZE];
  float *z_p_tmp = (float *)&varsf[BLOCK_SIZE * 2];
  float *h_tmp = (float *)&varsf[BLOCK_SIZE * 3];
  float *mass_tmp = (float *)&varsf[BLOCK_SIZE * 4];
  float *ux_tmp = (float *)&varsf[BLOCK_SIZE * 5];
  float *uy_tmp = (float *)&varsf[BLOCK_SIZE * 6];
  float *uz_tmp = (float *)&varsf[BLOCK_SIZE * 7];
  float *cj_tmp = (float *)&varsf[BLOCK_SIZE * 8];
  float *alphavisc_tmp = (float *)&varsf[BLOCK_SIZE * 9];
  float *alphadiff_tmp = (float *)&varsf[BLOCK_SIZE * 10];
  float *u_tmp = (float *)&varsf[BLOCK_SIZE * 11];
  float *rho_tmp = (float *)&varsf[BLOCK_SIZE * 12];
  float *pressure_tmp = (float *)&varsf[BLOCK_SIZE * 13];
  float *f_tmp = (float *)&varsf[BLOCK_SIZE * 14];
  float *balsara_tmp = (float *)&varsf[BLOCK_SIZE * 15];
  int *timebin = (int *)&varsf[BLOCK_SIZE * 16];
  /*Particles copied in blocks to shared memory*/
  for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
       b += BLOCK_SIZE) {
    int j = b + threadIdx.x;
    x_p_tmp[threadIdx.x] = parts_aos[j].x_p;
    y_p_tmp[threadIdx.x] = parts_aos[j].y_p;
    z_p_tmp[threadIdx.x] = parts_aos[j].z_p;
    h_tmp[threadIdx.x] = parts_aos[j].h;
    mass_tmp[threadIdx.x] = parts_aos[j].mass;
    ux_tmp[threadIdx.x] = parts_aos[j].ux;
    uy_tmp[threadIdx.x] = parts_aos[j].uy;
    uz_tmp[threadIdx.x] = parts_aos[j].uz;
    timebin[threadIdx.x] = parts_aos[j].time_bin;
    cj_tmp[threadIdx.x] = parts_aos[j].soundspeed;
//    alpha_tmp[threadIdx.x] = parts_aos[j].visc_alpha;
    u_tmp[threadIdx.x] = parts_aos[j].u;
    rho_tmp[threadIdx.x] = parts_aos[j].rho;
    alphavisc_tmp[threadIdx.x] = parts_aos[j].alpha_visc;
    alphadiff_tmp[threadIdx.x] = parts_aos[j].alpha_diff;
    pressure_tmp[threadIdx.x] = parts_aos[j].pressure;
    f_tmp[threadIdx.x] = parts_aos[j].f;
    balsara_tmp[threadIdx.x] = parts_aos[j].balsara;
    __syncthreads();
    for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
      j = j_block + b;
      if (j < last_part_in_task_blocks) {
        /* Compute the pairwise distance. */
        const float pjx = x_p_tmp[j_block] - cellx;
        const float pjy = y_p_tmp[j_block] - celly;
        const float pjz = z_p_tmp[j_block] - cellz;
        const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
        const float r2 = xij * xij + yij * yij + zij * zij;
        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {

          //          /* Cosmology terms for the signal velocity */
          const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
          const float a2_Hubble = d_a * d_a * d_H;
          const float r = sqrt(r2);
          const float r_inv = 1.f / r;
//          /* Recover some data */
          const float mj = mass_tmp[j_block];
//          /* Get the kernel for hi. */
          const float hi_inv = 1.f / hi;
          const float hid_inv = d_pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
          const float xi = r * hi_inv;
          float wi, wi_dx;
          d_kernel_deval(xi, &wi, &wi_dx);
          const float wi_dr = hid_inv * wi_dx;
          /* Get the kernel for hj. */
          const float hj = h_tmp[j_block];
          const float hj_inv = 1.0f / hj;
          const float hjd_inv = d_pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
          const float xj = r * hj_inv;
          float wj, wj_dx;
          d_kernel_deval(xj, &wj, &wj_dx);
          const float wj_dr = hjd_inv * wj_dx;
//          /* Compute dv dot r */
          float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
                dvz = uzi - uz_tmp[j_block];
          const float dvdr = dvx * xij + dvy * yij + dvz * zij;
//          /* Add Hubble flow */
          const float dvdr_Hubble = dvdr + a2_Hubble * r2;
//          /* Are the particles moving towards each others ? */
          const float omega_ij = min(dvdr_Hubble, 0.f);
          const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */
//
//          /* Signal velocity */
          const float v_sig = ci + cj_tmp[j_block] - const_viscosity_beta * mu_ij;

          /* Variable smoothing length term */
          const float f_ij = 1.f - fi / mj;
          const float f_ji = 1.f - f_tmp[j_block] / mi;

          /* Balsara term */
          const float balsaraj = balsara_tmp[j_block];
          /* Construct the full viscosity term */
          const float rhoj = rho_tmp[j_block];
          const float pressurej = pressure_tmp[j_block];
          const float rho_ij = rhoi + rhoj;
          const float alpha = alphavisci + alphavisc_tmp[j_block];
          const float visc =
              -0.25f * alpha * v_sig * mu_ij * (balsarai + balsaraj) / rho_ij;
          /* Convolve with the kernel */
          const float visc_acc_term =
              0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;
          /* Compute gradient terms */
          const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
          const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

          /* SPH acceleration term */
          const float sph_acc_term =
              (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

          /* Assemble the acceleration */
          const float acc = sph_acc_term + visc_acc_term;
          /* Use the force Luke ! */
          ahydroxi -= mj * acc * xij;
          ahydroyi -= mj * acc * yij;
          ahydrozi -= mj * acc * zij;
//          if(rhoi == 0 || rhoj == 0 || pressurei == 0 || pressurej == 0)printf("ri %f rj %f pi %f pj %f\n", rhoi, rhoj, pressurei, pressurej);
          /* Get the time derivative for u. */
          const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

          /* Viscosity term */
          const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;
          const float press_sum = pressurei + pressurej;
          /* Diffusion term */
          /* Combine the alpha_diff into a pressure-based switch -- this allows the
           * alpha from the highest pressure particle to dominate, so that the
           * diffusion limited particles always take precedence - another trick to
           * allow the scheme to work with thermal feedback. */
          float alpha_diff =
              (pressurei * alphadiffi + pressurej * alphadiff_tmp[j_block]) /
              (press_sum);
          if (fabsf(press_sum) < 1e-10) alpha_diff = 0.f;
          const float v_diff = alpha_diff * 0.5f *
                               (sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij) +
                                fabsf(fac_mu * r_inv * dvdr_Hubble));
          /* wi_dx + wj_dx / 2 is F_ij */
          const float diff_du_term =
              v_diff * (ui - u_tmp[j_block]) * (f_ij * wi_dr / rhoi + f_ji * wj_dr / rhoj);

          /* Assemble the energy equation term */
          const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

          /* Internal energy time derivative */
          u_dti += du_dt_i * mj;
          if(mj == 0.f)printf("zero mass mj %f\n", mj);

          /* Get the time derivative for h. */
          h_dti -= mj * dvdr * r_inv / rhoj * wi_dr;

          /* Update if we need to; this should be guaranteed by the gradient loop but
           * due to some possible synchronisation problems this is here as a _quick
           * fix_. Added: 14th August 2019. To be removed by 1st Jan 2020. (JB) */
          v_sigi = max(v_sigi, v_sig);
          int time_bin_j = timebin[j_block];
          if(time_bin_j > 0)min_ngb_time_bin = min(min_ngb_time_bin, time_bin_j);
//          printf("Got in\n");
        }
      }
    }
    __syncthreads();
  }
  if (pid < last_part_in_task_blocks) {
    parts_aos[pid].v_sig = v_sigi;
    parts_aos[pid].h_dt = h_dti;
    parts_aos[pid].u_dt = u_dti;
    parts_aos[pid].a_hydrox = ahydroxi;
    parts_aos[pid].a_hydroy = ahydroyi;
    parts_aos[pid].a_hydroz = ahydrozi;
    parts_aos[pid].min_ngb_time_bin = min_ngb_time_bin;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void DOSELF_GPU_AOS_F4_F(
	struct part_aos_f4_f_send * __restrict__ parts_send, struct part_aos_f4_f_recv * __restrict__ parts_recv,
	const float d_a, const float d_H,
    const int bundle_first_task, const int2 * __restrict__ d_task_first_part_f4) {

  extern __shared__ float4 varsf4_f[];

  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  int first_part_in_task_blocks, last_part_in_task_blocks;
//  first_part_in_task_blocks = d_task_first_part[task_id],
//  last_part_in_task_blocks = d_task_last_part[task_id];
  int2 first_last_parts = d_task_first_part_f4[task_id];
  first_part_in_task_blocks = first_last_parts.x;
  last_part_in_task_blocks = first_last_parts.y;

  const int pid = threadid + first_part_in_task_blocks;

  int ttid = 0;
  int first_part = 0;
  int count = 0;
  int last_part = 0;
  const part_aos_f4_f_send pi = parts_send[pid];
  float4 x_h_i = pi.x_h;
  float4 ux_m_i = pi.ux_m;
  float4 f_b_t_mintbinngb_i = pi.f_bals_timebin_mintimebin_ngb;
  float4 rho_p_c_vsig_i = pi.rho_p_c_vsigi;
  float3 u_avisc_adiff_i = pi.u_alphavisc_alphadiff;

  const float mi = ux_m_i.w;
  int Found_neighbours = 0;
  float pressurei = rho_p_c_vsig_i.y;
  const float ci = rho_p_c_vsig_i.z;
  float3 ahydro = {0.0, 0.0, 0.0};
  float4 udt_hdt_vsig_mintbinngb = {0.0, 0.0, 0.0, 0.0};
  udt_hdt_vsig_mintbinngb.z = rho_p_c_vsig_i.w;
  udt_hdt_vsig_mintbinngb.w = f_b_t_mintbinngb_i.w;

  float hi = x_h_i.w;
  float hig2 = hi * hi * kernel_gamma2;

  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float4 * __restrict__ x_h_tmp = (float4 *)&varsf4_f[0];
  float4 * __restrict__ ux_m_tmp = (float4 *)&varsf4_f[BLOCK_SIZE];
  float4 * __restrict__ f_b_t_mintbinngb_tmp = (float4 *)&varsf4_f[BLOCK_SIZE * 2];
  float4 * __restrict__ rho_p_c_vsig_tmp = (float4 *)&varsf4_f[BLOCK_SIZE * 3];
  float3 * __restrict__ u_avisc_adiff_tmp = (float3 *)&varsf4_f[BLOCK_SIZE * 4];
  /*Particles copied in blocks to shared memory*/
  for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
       b += BLOCK_SIZE) {
    int j = b + threadIdx.x;
    struct part_aos_f4_f_send pj = parts_send[j];
    x_h_tmp[threadIdx.x] = pj.x_h;
    ux_m_tmp[threadIdx.x] = pj.ux_m;
    f_b_t_mintbinngb_tmp[threadIdx.x] = pj.f_bals_timebin_mintimebin_ngb;
    rho_p_c_vsig_tmp[threadIdx.x] = pj.rho_p_c_vsigi;
//    alpha_tmp[threadIdx.x] = parts_aos[j].visc_alpha;
    u_avisc_adiff_tmp[threadIdx.x] = pj.u_alphavisc_alphadiff;
    __syncthreads();
    for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
      j = j_block + b;
      if (j < last_part_in_task_blocks) {
        /* Compute the pairwise distance. */
    	float4 x_h_j = x_h_tmp[j_block];
    	float4 ux_m_j = ux_m_tmp[j_block];
    	float4 f_b_t_mintbinngb_j = f_b_t_mintbinngb_tmp[j_block];
    	float4 rho_p_c_vsig_j = rho_p_c_vsig_tmp[j_block];
    	float3 u_avisc_adiff_j = u_avisc_adiff_tmp[j_block];
        const float xij = x_h_i.x - x_h_j.x,
        		    yij = x_h_i.y - x_h_j.y,
        		    zij = x_h_i.z - x_h_j.z;
        const float r2 = xij * xij + yij * yij + zij * zij;
        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
          //          /* Cosmology terms for the signal velocity */
          const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
          const float a2_Hubble = d_a * d_a * d_H;
          const float r = sqrt(r2);
          const float r_inv = 1.f / r;
//          /* Recover some data */
          const float mj = ux_m_j.w;
//          /* Get the kernel for hi. */
          const float hi_inv = 1.f / hi;
          const float hid_inv = d_pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
          const float xi = r * hi_inv;
          float wi, wi_dx;
          d_kernel_deval(xi, &wi, &wi_dx);
          const float wi_dr = hid_inv * wi_dx;
          /* Get the kernel for hj. */
          const float hj = x_h_j.w;
          const float hj_inv = 1.0f / hj;
          const float hjd_inv = d_pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
          const float xj = r * hj_inv;
          float wj, wj_dx;
          d_kernel_deval(xj, &wj, &wj_dx);
          const float wj_dr = hjd_inv * wj_dx;
//          /* Compute dv dot r */
          float dvx = ux_m_i.x - ux_m_j.x, dvy = ux_m_i.y - ux_m_j.y,
                dvz = ux_m_i.z - ux_m_j.z;
          const float dvdr = dvx * xij + dvy * yij + dvz * zij;
//          /* Add Hubble flow */
          const float dvdr_Hubble = dvdr + a2_Hubble * r2;
//          /* Are the particles moving towards each others ? */
          const float omega_ij = min(dvdr_Hubble, 0.f);
          const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */
//
//          /* Signal velocity */
          const float cj = rho_p_c_vsig_j.z;
          const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

          /* Variable smoothing length term */
          const float f_ij = 1.f - f_b_t_mintbinngb_i.x / mj;
          const float f_ji = 1.f - f_b_t_mintbinngb_j.x / mi;

          /* Construct the full viscosity term */
          const float pressurej = rho_p_c_vsig_j.y;
          const float rho_ij = rho_p_c_vsig_i.x + rho_p_c_vsig_j.x;
          const float alpha = u_avisc_adiff_i.y + u_avisc_adiff_j.y;
          const float visc =
              -0.25f * alpha * v_sig * mu_ij * (f_b_t_mintbinngb_i.y + f_b_t_mintbinngb_j.y) / rho_ij;
          /* Convolve with the kernel */
          const float visc_acc_term =
              0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;
          /* Compute gradient terms */
          const float rhoi2 = rho_p_c_vsig_i.x * rho_p_c_vsig_i.x;
          const float rhoj2 = rho_p_c_vsig_j.x * rho_p_c_vsig_j.x;
          const float P_over_rho2_i = pressurei / (rhoi2) * f_ij;
          const float P_over_rho2_j = pressurej / (rhoj2) * f_ji;

          /* SPH acceleration term */
          const float sph_acc_term =
              (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

          /* Assemble the acceleration */
          const float acc = sph_acc_term + visc_acc_term;
          /* Use the force Luke ! */
          ahydro.x -= mj * acc * xij;
          ahydro.y -= mj * acc * yij;
          ahydro.z -= mj * acc * zij;
          /* Get the time derivative for u. */
          const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

          /* Viscosity term */
          const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;
          /* Diffusion term */
          /* Combine the alpha_diff into a pressure-based switch -- this allows the
           * alpha from the highest pressure particle to dominate, so that the
           * diffusion limited particles always take precedence - another trick to
           * allow the scheme to work with thermal feedback. */
          float alpha_diff =
              (pressurei * u_avisc_adiff_i.z + pressurej * u_avisc_adiff_j.z) /
              (pressurei + pressurej);
          if (fabsf(pressurei + pressurej) < 1e-10) alpha_diff = 0.f;
          const float v_diff = alpha_diff * 0.5f *
                               (sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij) +
                                fabsf(fac_mu * r_inv * dvdr_Hubble));
          /* wi_dx + wj_dx / 2 is F_ij */
          const float diff_du_term =
              v_diff * (u_avisc_adiff_i.x - u_avisc_adiff_j.x) *
			  (f_ij * wi_dr / rho_p_c_vsig_i.x + f_ji * wj_dr / rho_p_c_vsig_j.x);

          /* Assemble the energy equation term */
          const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

          /* Internal energy time derivative */
          udt_hdt_vsig_mintbinngb.x += du_dt_i * mj;

          /* Get the time derivative for h. */
          udt_hdt_vsig_mintbinngb.y -= mj * dvdr * r_inv / rho_p_c_vsig_j.x * wi_dr;

          /* Update if we need to; this should be guaranteed by the gradient loop but
           * due to some possible synchronisation problems this is here as a _quick
           * fix_. Added: 14th August 2019. To be removed by 1st Jan 2020. (JB) */
          udt_hdt_vsig_mintbinngb.z = fmaxf(udt_hdt_vsig_mintbinngb.z, v_sig);
          unsigned int time_bin_j = (f_b_t_mintbinngb_j.z + 0.5f);
          unsigned int min_tb_i = (f_b_t_mintbinngb_i.w + 0.5f);
          if(time_bin_j > 0)f_b_t_mintbinngb_i.w =
        		  min(min_tb_i, time_bin_j);
//          printf("Got in\n");
        }
      }
    }
    __syncthreads();
  }
  if (pid < last_part_in_task_blocks) {
	udt_hdt_vsig_mintbinngb.w = f_b_t_mintbinngb_i.w;
    parts_recv[pid].udt_hdt_vsig_mintimebin_ngb = udt_hdt_vsig_mintbinngb;
    parts_recv[pid].a_hydro = ahydro;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_density_GPU_naive(
	struct part_soa parts_soa_ci, struct part_soa parts_soa_cj, int *d_task_first_part_ci, int *d_task_first_part_cj,
	int *d_task_last_part_ci, int *d_task_last_part_cj, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, int time_bin_inhibited) {

  extern __shared__ float vars[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
  __shared__ int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
  __shared__ int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;

  first_part_in_task_blocks_ci = d_task_first_part_ci[task_id];
  last_part_in_task_blocks_ci = d_task_last_part_ci[task_id];
  first_part_in_task_blocks_cj = d_task_first_part_cj[task_id];
  last_part_in_task_blocks_cj = d_task_last_part_cj[task_id];

  __syncthreads();
  // Now we start calculations for particles in cell i
  const int pid = threadid + first_part_in_task_blocks_ci;

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves
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
  int Found_neighbours = 0;
  //	if(pid<b_last_part&&pid<last_part_in_task_blocks){
  if (pid < last_part_in_task_blocks_ci) {
    ttid = parts_soa_ci.tid_p[pid];
    first_part = d_task_first_part_ci[ttid];
    last_part = d_task_last_part_ci[ttid];
    count = last_part - first_part;
    cellx = parts_soa_ci.locx[pid], celly = parts_soa_ci.locy[pid],
    cellz = parts_soa_ci.locz[pid];
    hi = parts_soa_ci.h[pid], hig2 = hi * hi * kernel_gamma2;
    mi = parts_soa_ci.mass[pid];
    uxi = parts_soa_ci.ux[pid];
    uyi = parts_soa_ci.uy[pid];
    uzi = parts_soa_ci.uz[pid];
    pix = parts_soa_ci.x_p[pid] - cellx;
    piy = parts_soa_ci.y_p[pid] - celly;
    piz = parts_soa_ci.z_p[pid] - cellz;
  }
//  if (threadIdx.x == 0) {
//    first_part_tid_0 = first_part;
//    last_part_tid_0 = last_part;
//  }
//  __syncthreads();
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
  for (int b = first_part_in_task_blocks_cj; b < last_part_in_task_blocks_cj;
       b += BLOCK_SIZE) {
    int j = b + threadIdx.x;
    x_p_tmp[threadIdx.x] = parts_soa_cj.x_p[j];
    y_p_tmp[threadIdx.x] = parts_soa_cj.y_p[j];
    z_p_tmp[threadIdx.x] = parts_soa_cj.z_p[j];
    h_tmp[threadIdx.x] = parts_soa_cj.h[j];
    mass_tmp[threadIdx.x] = parts_soa_cj.mass[j];
    ux_tmp[threadIdx.x] = parts_soa_cj.ux[j];
    uy_tmp[threadIdx.x] = parts_soa_cj.uy[j];
    uz_tmp[threadIdx.x] = parts_soa_cj.uz[j];
    timebin[threadIdx.x] = parts_soa_cj.time_bin[j];
    __syncthreads();
    for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
      j = j_block + b;
//      if ((j != pid) && (j < last_part_in_task_blocks) &&
//          timebin[j_block] != time_bin_inhibited) {
//      if ((j < last_part_in_task_blocks) &&
//    	  timebin[j_block] != time_bin_inhibited) {
      if (j < last_part_in_task_blocks_cj) {
        /* Compute the pairwise distance. */
        const float pjx = x_p_tmp[j_block] - cellx;
        const float pjy = y_p_tmp[j_block] - celly;
        const float pjz = z_p_tmp[j_block] - cellz;
        const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
        const float r2 = xij * xij + yij * yij + zij * zij;
//        const float hj = h_tmp[j_block], hjg2 = hj * hj * kernel_gamma2;
        //				if((hi < 0.0001f || hj < 0.0001f || r2 <
        //0.0000001f) && pid < last_part_in_task_blocks){ 					printf("very small
        //value for hi %f or hj %f or r2 %f\n", hi, hj, r2);
        //				}
//        if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
        if (r2 < hig2 && r2 > (0.01f/dx)*(0.01f/dx)) {
          Found_neighbours=1;
          const float r = sqrt(r2);
          /* Recover some data */
          const float mj = mass_tmp[j_block];
          /* Get the kernel for hi. */
          if(hi<1.f/dx)printf("h < dx\n");
//          if(hi<1.f/256.f)printf("h < dx\n");
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
  if (pid < last_part_in_task_blocks_ci) {
    parts_soa_ci.rho[pid] = rhoi, parts_soa_ci.rho_dh[pid] = rho_dhi;
    parts_soa_ci.wcount[pid] = wcounti, parts_soa_ci.wcount_dh[pid] = wcount_dhi;
    parts_soa_ci.div_v[pid] = div_vi;
    parts_soa_ci.rot_ux[pid] = rot_uxi, parts_soa_ci.rot_uy[pid] = rot_uyi;
    parts_soa_ci.rot_uz[pid] = rot_uzi;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_density_pair_two_kernels(struct part_soa parts_soa_ci, struct part_soa parts_soa_cj, int *d_task_first_part_ci,
	          int *d_task_first_part_cj, int *d_task_last_part_ci, int *d_task_last_part_cj, float d_a, float d_H,
              const char *loop_type, cudaStream_t stream, int bid, int block_size, int count_tasks, int tasksperbundle,
              int max_parts_i, int max_parts_j, int numBlocks_y, int tid, int offset, int bundle_first_task, int time_bin_inhibited) {


  int max_parts = max(max_parts_j, max_parts_i);
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;

  /*Do ci*/
  runner_do_pair_density_GPU_naive<<<gridShape, BLOCK_SIZE,
                               8 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(timebin_t),
                               stream>>>(
      parts_soa_ci, parts_soa_cj, d_task_first_part_ci, d_task_first_part_cj, d_task_last_part_ci,
	  d_task_last_part_cj, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, time_bin_inhibited);

//  numBlocks_x = (max_parts_i + BLOCK_SIZE - 1) / BLOCK_SIZE;
//  gridShape = dim3(numBlocks_x, numBlocks_y);
//  nBlocks_per_task = numBlocks_x;
  /*Now do cj*/
  runner_do_pair_density_GPU_naive<<<gridShape, BLOCK_SIZE,
                               8 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(timebin_t),
                               stream>>>(
      parts_soa_cj, parts_soa_ci, d_task_first_part_cj, d_task_first_part_ci, d_task_last_part_cj,
	  d_task_last_part_ci, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIRGPU(
	struct part_soa parts_soa, int pid,
	int last_part_in_task_blocks_ci, int first_part_in_task_blocks_cj,
	int last_part_in_task_blocks_cj, float d_a, float d_H,
	int time_bin_inhibited, float *vars) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

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
  int Found_neighbours = 0;

  if (pid < last_part_in_task_blocks_ci) {
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
  timebin_t *timebin = (timebin_t *)&uz_tmp[BLOCK_SIZE];
  /*Particles copied in blocks to shared memory*/
  for (int b = first_part_in_task_blocks_cj; b < last_part_in_task_blocks_cj;
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
	  if (j < last_part_in_task_blocks_cj) {
		/* Compute the pairwise distance. */
		const float pjx = x_p_tmp[j_block] - cellx;
		const float pjy = y_p_tmp[j_block] - celly;
		const float pjz = z_p_tmp[j_block] - cellz;
		const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
		const float r2 = xij * xij + yij * yij + zij * zij;

		if (r2 < hig2 && r2 > (0.01f/dx)*(0.01f/dx)) {
		  Found_neighbours=1;
		  const float r = sqrt(r2);
		  /* Recover some data */
		  const float mj = mass_tmp[j_block];
		  /* Get the kernel for hi. */
		  if(hi<1.f/dx)printf("h < dx\n");
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
  if (pid < last_part_in_task_blocks_ci) {
	parts_soa.rho[pid] = rhoi, parts_soa.rho_dh[pid] = rho_dhi;
	parts_soa.wcount[pid] = wcounti, parts_soa.wcount_dh[pid] = wcount_dhi;
	parts_soa.div_v[pid] = div_vi;
	parts_soa.rot_ux[pid] = rot_uxi, parts_soa.rot_uy[pid] = rot_uyi;
	parts_soa.rot_uz[pid] = rot_uzi;
  }

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NONSYMGPU(
	struct part_soa parts_soa, int pid, const int ci_start,
	const int ci_end, const int cj_start,
	const int cj_end, float d_a, float d_H,
	float *vars_pair, double *d_shift_x, double *d_shift_y, double *d_shift_z, const int task_id_tmp, int flip_order) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

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
  int Found_neighbours = 0;
  int count_i = cj_start;
//  printf("first_part_in_task_blocks_cj %i last_part_in_task_blocks_cj %i last_part_in_task_blocks_ci %i\n",
//		  first_part_in_task_blocks_cj, last_part_in_task_blocks_cj, last_part_in_task_blocks_ci);
  if (pid < ci_end) {
	hi = parts_soa.h[pid], hig2 = hi * hi * kernel_gamma2;
	mi = parts_soa.mass[pid];
	uxi = parts_soa.ux[pid];
	uyi = parts_soa.uy[pid];
	uzi = parts_soa.uz[pid];
	pix = parts_soa.x_p[pid] - d_shift_x[task_id_tmp];
	piy = parts_soa.y_p[pid] - d_shift_y[task_id_tmp];
	piz = parts_soa.z_p[pid] - d_shift_z[task_id_tmp];
  }

  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float *x_p_tmp = (float *)&vars_pair[0];
  float *y_p_tmp = (float *)&x_p_tmp[BLOCK_SIZE];
  float *z_p_tmp = (float *)&y_p_tmp[BLOCK_SIZE];
  float *h_tmp = (float *)&z_p_tmp[BLOCK_SIZE];
  float *mass_tmp = (float *)&h_tmp[BLOCK_SIZE];
  float *ux_tmp = (float *)&mass_tmp[BLOCK_SIZE];
  float *uy_tmp = (float *)&ux_tmp[BLOCK_SIZE];
  float *uz_tmp = (float *)&uy_tmp[BLOCK_SIZE];
  timebin_t *timebin = (timebin_t *)&uz_tmp[BLOCK_SIZE];

  /*Particles copied in blocks to shared memory*/
  for (int b = cj_start; b < cj_end;
	   b += BLOCK_SIZE) {
	const int tid_x = threadIdx.x;
	int j = b + tid_x;
	x_p_tmp[tid_x] = parts_soa.x_p[j];
	y_p_tmp[tid_x] = parts_soa.y_p[j];
	z_p_tmp[tid_x] = parts_soa.z_p[j];
//	h_tmp[tid_x] = parts_soa.h[j];
	mass_tmp[tid_x] = parts_soa.mass[j];
	ux_tmp[tid_x] = parts_soa.ux[j];
	uy_tmp[tid_x] = parts_soa.uy[j];
	uz_tmp[tid_x] = parts_soa.uz[j];
	timebin[tid_x] = parts_soa.time_bin[j];

	__syncthreads();
	const float shift_x_j = d_shift_x[task_id_tmp + flip_order];
	const float shift_y_j = d_shift_y[task_id_tmp + flip_order];
	const float shift_z_j = d_shift_z[task_id_tmp + flip_order];
	/*j_block is the particle's index in the block. Loop through particles in shared memory one by one*/
	for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
	  int jj = b + j_block;
	  if (jj < cj_end && pid < ci_end && pid >= ci_start) {

		const float pjx = x_p_tmp[j_block] - shift_x_j;
		const float pjy = y_p_tmp[j_block] - shift_y_j;
		const float pjz = z_p_tmp[j_block] - shift_z_j;

		const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
//		const float xij = (pix - pjx) * flip_order, yij = (piy - pjy) * flip_order, zij = (piz - pjz) * flip_order;
		const float r2 = xij * xij + yij * yij + zij * zij;
		if (r2 < hig2) {
		  /* Recover some data */
		  const float mj = mass_tmp[j_block];
		  const float r = sqrt(r2);
		  /* Get the kernel for hi. */
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
		  const float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
		  dvz = uzi - uz_tmp[j_block];
		  const float dvdr = dvx * xij + dvy * yij + dvz * zij;
		  /* Compute dv cross r */
		  const float curlvrx = dvy * zij - dvz * yij;
		  const float curlvry = dvz * xij - dvx * zij;
		  const float curlvrz = dvx * yij - dvy * xij;

		  div_vi -= faci * dvdr;

		  rot_uxi += faci * curlvrx;
		  rot_uyi += faci * curlvry;
		  rot_uzi += faci * curlvrz;
		}
	  } /*if (jj < cj_end && pid < ci_end && pid >= ci_start)*/
	} /*End of looping through particles in shared memory---Shared arrays zero'ed for next step in outer loop*/
	__syncthreads();
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
  if (pid >= ci_start && pid < ci_end) {
	parts_soa.rho[pid] = rhoi, parts_soa.rho_dh[pid] = rho_dhi;
	parts_soa.wcount[pid] = wcounti, parts_soa.wcount_dh[pid] = wcount_dhi;
	parts_soa.div_v[pid] = div_vi;
	parts_soa.rot_ux[pid] = rot_uxi, parts_soa.rot_uy[pid] = rot_uyi;
	parts_soa.rot_uz[pid] = rot_uzi;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NONSYMGPUAOS(
	struct part_aos *parts_aos, int pid, const int ci_start,
	const int ci_end, const int cj_start,
	const int cj_end, float d_a, float d_H,
	float *vars_pair_aos, double *d_shift_x, double *d_shift_y, double *d_shift_z, const int task_id_tmp, int flip_order) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

  float hi = 0.0, hig2 = 0.0;

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
  int Found_neighbours = 0;
  int count_i = cj_start;
//  printf("first_part_in_task_blocks_cj %i last_part_in_task_blocks_cj %i last_part_in_task_blocks_ci %i\n",
//		  first_part_in_task_blocks_cj, last_part_in_task_blocks_cj, last_part_in_task_blocks_ci);
  if (pid < ci_end) {
	hi = parts_aos[pid].h, hig2 = hi * hi * kernel_gamma2;
	mi = parts_aos[pid].mass;
	uxi = parts_aos[pid].ux;
	uyi = parts_aos[pid].uy;
	uzi = parts_aos[pid].uz;
	pix = parts_aos[pid].x_p;// - d_shift_x[task_id_tmp];
	piy = parts_aos[pid].y_p;// - d_shift_y[task_id_tmp];
	piz = parts_aos[pid].z_p;// - d_shift_z[task_id_tmp];
  }

  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float *x_p_tmp = (float *)&vars_pair_aos[0];
  float *y_p_tmp = (float *)&x_p_tmp[BLOCK_SIZE];
  float *z_p_tmp = (float *)&y_p_tmp[BLOCK_SIZE];
  float *h_tmp = (float *)&z_p_tmp[BLOCK_SIZE];
  float *mass_tmp = (float *)&h_tmp[BLOCK_SIZE];
  float *ux_tmp = (float *)&mass_tmp[BLOCK_SIZE];
  float *uy_tmp = (float *)&ux_tmp[BLOCK_SIZE];
  float *uz_tmp = (float *)&uy_tmp[BLOCK_SIZE];
  int *timebin = (int *)&uz_tmp[BLOCK_SIZE];

  /*Particles copied in blocks to shared memory*/
  for (int b = cj_start; b < cj_end;
	   b += BLOCK_SIZE) {
	const int tid_x = threadIdx.x;
	int j = b + tid_x;
	x_p_tmp[tid_x] = parts_aos[j].x_p;
	y_p_tmp[tid_x] = parts_aos[j].y_p;
	z_p_tmp[tid_x] = parts_aos[j].z_p;
//	h_tmp[tid_x] = parts_aos[j].h;
	mass_tmp[tid_x] = parts_aos[j].mass;
	ux_tmp[tid_x] = parts_aos[j].ux;
	uy_tmp[tid_x] = parts_aos[j].uy;
	uz_tmp[tid_x] = parts_aos[j].uz;
	timebin[tid_x] = parts_aos[j].time_bin;
//	const float shift_x_j = d_shift_x[task_id_tmp + flip_order];
//	const float shift_y_j = d_shift_y[task_id_tmp + flip_order];
//	const float shift_z_j = d_shift_z[task_id_tmp + flip_order];
	__syncthreads();
	/*j_block is the particle's index in the block. Loop through particles in shared memory one by one*/
	for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
	  int jj = b + j_block;
	  if (jj < cj_end && pid < ci_end && pid >= ci_start) {

		const float pjx = x_p_tmp[j_block];// - shift_x_j;
		const float pjy = y_p_tmp[j_block];// - shift_y_j;
		const float pjz = z_p_tmp[j_block];// - shift_z_j;

		const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
//		const float xij = (pix - pjx) * flip_order, yij = (piy - pjy) * flip_order, zij = (piz - pjz) * flip_order;
		const float r2 = xij * xij + yij * yij + zij * zij;
		if (r2 < hig2) {
		  /* Recover some data */
		  const float mj = mass_tmp[j_block];
		  const float r = sqrt(r2);
		  /* Get the kernel for hi. */
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
		  const float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
		  dvz = uzi - uz_tmp[j_block];
		  const float dvdr = dvx * xij + dvy * yij + dvz * zij;
		  /* Compute dv cross r */
		  const float curlvrx = dvy * zij - dvz * yij;
		  const float curlvry = dvz * xij - dvx * zij;
		  const float curlvrz = dvx * yij - dvy * xij;

		  div_vi -= faci * dvdr;

		  rot_uxi += faci * curlvrx;
		  rot_uyi += faci * curlvry;
		  rot_uzi += faci * curlvrz;
//		  if(timebin[j_block] != 1000 && timebin[j_block] != 20)printf("incorrect timebin %i\n", timebin[j_block]);
		}
	  } /*if (jj < cj_end && pid < ci_end && pid >= ci_start)*/
	} /*End of looping through particles in shared memory---Shared arrays zero'ed for next step in outer loop*/
	__syncthreads();
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
  if (pid >= ci_start && pid < ci_end) {
//	printf("timebin %i\n", parts_aos[pid].time_bin);
	parts_aos[pid].rho = rhoi, parts_aos[pid].rho_dh = rho_dhi;
	parts_aos[pid].wcount = wcounti, parts_aos[pid].wcount_dh = wcount_dhi;
	parts_aos[pid].div_v = div_vi;
	parts_aos[pid].rot_ux = rot_uxi, parts_aos[pid].rot_uy = rot_uyi;
	parts_aos[pid].rot_uz = rot_uzi;
	parts_aos[pid].time_bin = 20;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NONSYMGPUAOSF4(
	struct part_aos_f4_send * __restrict__ parts_send, struct part_aos_f4_recv * __restrict__ parts_recv, int pid,
	const int ci_start, const int ci_end, const int cj_start, const int cj_end, float d_a, float d_H, float4 *vars_pair_aos_f4) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

  float hi = 0.0, hig2 = 0.0;

  int Found_neighbours = 0;
  int count_i = cj_start;

  float4 res_rho = {0.0, 0.0, 0.0, 0.0};
  float4 res_rot = {0.0, 0.0, 0.0, 0.0};
  const part_aos_f4_send pi = parts_send[pid];
  const float4 x_pi = pi.x_p_h;
  const float4 ux_pi = pi.ux_m;
//  printf("first_part_in_task_blocks_cj %i last_part_in_task_blocks_cj %i last_part_in_task_blocks_ci %i\n",
//		  first_part_in_task_blocks_cj, last_part_in_task_blocks_cj, last_part_in_task_blocks_ci);
//  if (pid < ci_end) {
	hi = x_pi.w, hig2 = hi * hi * kernel_gamma2;
//  }

  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float4 * __restrict__ x_p_h_tmp = (float4 *)&vars_pair_aos_f4[0];
  float4 * __restrict__ ux_m_tmp = (float4 *)&vars_pair_aos_f4[BLOCK_SIZE];

  /*Particles copied in blocks to shared memory*/
  for (int b = cj_start; b < cj_end;
	   b += BLOCK_SIZE) {
	const int tid_x = threadIdx.x;
	int j = b + tid_x;
    struct part_aos_f4_send pj = parts_send[j];
	x_p_h_tmp[tid_x] = pj.x_p_h;
	ux_m_tmp[tid_x] = pj.ux_m;
	__syncthreads();
	/*j_block is the particle's index in the block. Loop through particles in shared memory one by one*/
	for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
	  int jj = b + j_block;
	  if (jj < cj_end && pid < ci_end && pid >= ci_start) {

		const float4 x_p_h_j = x_p_h_tmp[j_block];
        const float4 ux_m_j = ux_m_tmp[j_block];

		const float xij = x_pi.x - x_p_h_j.x, yij = x_pi.y - x_p_h_j.y,
		zij = x_pi.z - x_p_h_j.z;
		const float r2 = xij * xij + yij * yij + zij * zij;
		if (r2 < hig2) {
		  /* Recover some data */
		  const float mj = ux_m_j.w;
		  const float r = sqrt(r2);
		  /* Get the kernel for hi. */
		  const float h_inv = 1.f / hi;
		  const float ui = r * h_inv;
		  float wi, wi_dx;

		  d_kernel_deval(ui, &wi, &wi_dx);
		  /*Add to sums of rho, rho_dh, wcount and wcount_dh*/
		  res_rho.x += mj * wi;
		  res_rho.y -= mj * (hydro_dimension * wi + ui * wi_dx);
		  res_rho.z += wi;
		  res_rho.w -= (hydro_dimension * wi + ui * wi_dx);

		  const float r_inv = 1.f / r;
		  const float faci = mj * wi_dx * r_inv;
		  /* Compute dv dot r */
		  const float dvx = ux_pi.x - ux_m_j.x, dvy = ux_pi.y - ux_m_j.y,
		  dvz = ux_pi.z - ux_m_j.z;
		  const float dvdr = dvx * xij + dvy * yij + dvz * zij;
		  /* Compute dv cross r */
		  const float curlvrx = dvy * zij - dvz * yij;
		  const float curlvry = dvz * xij - dvx * zij;
		  const float curlvrz = dvx * yij - dvy * xij;

		  res_rot.x += faci * curlvrx;
		  res_rot.y += faci * curlvry;
		  res_rot.z += faci * curlvrz;
		  res_rot.w -= faci * dvdr;
		}
	  } /*if (jj < cj_end && pid < ci_end && pid >= ci_start)*/
	} /*End of looping through particles in shared memory---Shared arrays zero'ed for next step in outer loop*/
	__syncthreads();
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
  if (pid >= ci_start && pid < ci_end) {
    parts_recv[pid].rho_dh_wcount = res_rho;
	parts_recv[pid].rot_ux_div_v = res_rot;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NAIVEGPUAOSF4(const struct part_aos_f4_send pi,
	struct part_aos_f4_send * __restrict__ parts_send, struct part_aos_f4_recv * __restrict__ parts_recv, int pid,
	const int cj_start, const int cj_end, float d_a, float d_H) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

  float hi = 0.0, hig2 = 0.0;

  int Found_neighbours = 0;
  int count_i = cj_start;

  float4 res_rho = {0.0, 0.0, 0.0, 0.0};
  float4 res_rot = {0.0, 0.0, 0.0, 0.0};
//  const part_aos_f4_send pi = parts_send[pid];
  const float4 x_pi = pi.x_p_h;
  const float4 ux_pi = pi.ux_m;
//  printf("first_part_in_task_blocks_cj %i last_part_in_task_blocks_cj %i last_part_in_task_blocks_ci %i\n",
//		  first_part_in_task_blocks_cj, last_part_in_task_blocks_cj, last_part_in_task_blocks_ci);
//  if (pid < ci_end) {
	hi = x_pi.w, hig2 = hi * hi * kernel_gamma2;
//  }

//  printf("js %i je %i\n", cj_start, cj_end);
  /*Particles copied in blocks to shared memory*/
  for (int j = cj_start; j < cj_end; j ++) {
    struct part_aos_f4_send pj = parts_send[j];

	const float4 x_p_h_j = pj.x_p_h;
	const float4 ux_m_j = pj.ux_m;

	const float xij = x_pi.x - x_p_h_j.x, yij = x_pi.y - x_p_h_j.y,
	zij = x_pi.z - x_p_h_j.z;
	const float r2 = xij * xij + yij * yij + zij * zij;
//	printf("r2 %f \n", r2);
	if (r2 < hig2) {
	  /* Recover some data */
	  const float mj = ux_m_j.w;
	  const float r = sqrt(r2);
	  /* Get the kernel for hi. */
	  const float h_inv = 1.f / hi;
	  const float ui = r * h_inv;
	  float wi, wi_dx;

	  d_kernel_deval(ui, &wi, &wi_dx);
	  /*Add to sums of rho, rho_dh, wcount and wcount_dh*/
	  res_rho.x += mj * wi;
	  res_rho.y -= mj * (hydro_dimension * wi + ui * wi_dx);
	  res_rho.z += wi;
	  res_rho.w -= (hydro_dimension * wi + ui * wi_dx);

	  const float r_inv = 1.f / r;
	  const float faci = mj * wi_dx * r_inv;
	  /* Compute dv dot r */
	  const float dvx = ux_pi.x - ux_m_j.x, dvy = ux_pi.y - ux_m_j.y,
	  dvz = ux_pi.z - ux_m_j.z;
	  const float dvdr = dvx * xij + dvy * yij + dvz * zij;
	  /* Compute dv cross r */
	  const float curlvrx = dvy * zij - dvz * yij;
	  const float curlvry = dvz * xij - dvx * zij;
	  const float curlvrz = dvx * yij - dvy * xij;

	  res_rot.x += faci * curlvrx;
	  res_rot.y += faci * curlvry;
	  res_rot.z += faci * curlvrz;
	  res_rot.w -= faci * dvdr;
	}
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
//  if (pid >= ci_start && pid < ci_end) {
    parts_recv[pid].rho_dh_wcount = res_rho;
	parts_recv[pid].rot_ux_div_v = res_rot;
//  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NONSYMGPUAOSG(
	struct part_aos_g *parts_aos, int pid, const int ci_start,
	const int ci_end, const int cj_start,
	const int cj_end, float d_a, float d_H,
	float *vars_pair_aosg, double *d_shift_x, double *d_shift_y, double *d_shift_z, const int task_id_tmp, int flip_order) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

  float hi = 0.0, hig2 = 0.0;

  float mi = 0.0;
  float uxi = 0.0;
  float uyi = 0.0;
  float uzi = 0.0;
  float pix = 0.0;
  float piy = 0.0;
  float piz = 0.0;
  float rhoi = 0.0;
  float div_vi = 0.0;
  int Found_neighbours = 0;
  float v_sig;
  float u = 0.f;
  float laplace_u = 0.0;
  float alpha_visc_max_ngb = 0.0;
  float ci = 0.0;

  int count_i = cj_start;
  if (pid < ci_end) {
	hi = parts_aos[pid].h, hig2 = hi * hi * kernel_gamma2;
	mi = parts_aos[pid].mass;
	uxi = parts_aos[pid].ux;
	uyi = parts_aos[pid].uy;
	uzi = parts_aos[pid].uz;
	ci = parts_aos[pid].soundspeed;
	v_sig = parts_aos[pid].v_sig;
	u = parts_aos[pid].u;
	laplace_u = parts_aos[pid].laplace_u;
	alpha_visc_max_ngb = parts_aos[pid].alpha_visc_max_ngb;

	pix = parts_aos[pid].x_p - d_shift_x[task_id_tmp];
	piy = parts_aos[pid].y_p - d_shift_y[task_id_tmp];
	piz = parts_aos[pid].z_p - d_shift_z[task_id_tmp];
  }

  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float *x_p_tmp = (float *)&vars_pair_aosg[0];
  float *y_p_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE];
  float *z_p_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 2];
  float *h_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 3];
  float *mass_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 4];
  float *ux_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 5];
  float *uy_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 6];
  float *uz_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 7];
  float *cj_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 8];
  float *alpha_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 9];
  float *u_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 10];
  float *rho_tmp = (float *)&vars_pair_aosg[BLOCK_SIZE * 11];
  int *timebin = (int *)&vars_pair_aosg[BLOCK_SIZE * 12];

  /*Particles copied in blocks to shared memory*/
  for (int b = cj_start; b < cj_end;
	   b += BLOCK_SIZE) {
	const int tid_x = threadIdx.x;
	int j = b + tid_x;
    x_p_tmp[threadIdx.x] = parts_aos[j].x_p;
    y_p_tmp[threadIdx.x] = parts_aos[j].y_p;
    z_p_tmp[threadIdx.x] = parts_aos[j].z_p;
    h_tmp[threadIdx.x] = parts_aos[j].h;
    mass_tmp[threadIdx.x] = parts_aos[j].mass;
    ux_tmp[threadIdx.x] = parts_aos[j].ux;
    uy_tmp[threadIdx.x] = parts_aos[j].uy;
    uz_tmp[threadIdx.x] = parts_aos[j].uz;
    timebin[threadIdx.x] = parts_aos[j].time_bin;
    cj_tmp[threadIdx.x] = parts_aos[j].soundspeed;
    alpha_tmp[threadIdx.x] = parts_aos[j].visc_alpha;
    u_tmp[threadIdx.x] = parts_aos[j].u;
    rho_tmp[threadIdx.x] = parts_aos[j].rho;
	const float shift_x_j = d_shift_x[task_id_tmp + flip_order];
	const float shift_y_j = d_shift_y[task_id_tmp + flip_order];
	const float shift_z_j = d_shift_z[task_id_tmp + flip_order];
    __syncthreads();
	/*j_block is the particle's index in the block. Loop through particles in shared memory one by one*/
	for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
	  int jj = b + j_block;
	  if (jj < cj_end && pid < ci_end && pid >= ci_start) {

		const float pjx = x_p_tmp[j_block] - shift_x_j;
		const float pjy = y_p_tmp[j_block] - shift_y_j;
		const float pjz = z_p_tmp[j_block] - shift_z_j;
		const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
		const float r2 = xij * xij + yij * yij + zij * zij;
		if (r2 < hig2) {
		  /* Recover some data */
		  const float mj = mass_tmp[j_block];
		  const float r = sqrt(r2);
          const float r_inv = 1.f / r;
		  /* Get the kernel for hi. */
		  const float h_inv = 1.f / hi;
		  float wi, wi_dx;
          /* Cosmology terms for the signal velocity */
          const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
          const float a2_Hubble = d_a * d_a * d_H;
		  /* Compute dv dot r */
		  const float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
		              dvz = uzi - uz_tmp[j_block];
		  const float dvdr = dvx * xij + dvy * yij + dvz * zij;
          /* Add Hubble flow */
          const float dvdr_Hubble = dvdr + a2_Hubble * r2;
          /* Are the particles moving towards each others ? */
          const float omega_ij = min(dvdr_Hubble, 0.f);
          const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

          /* Signal velocity */
          const float new_v_sig = ci + cj_tmp[j_block] - const_viscosity_beta * mu_ij;
          /* Update if we need to */
          v_sig = max(v_sig, new_v_sig);
          /* Calculate Del^2 u for the thermal diffusion coefficient. */
          /* Need to get some kernel values F_ij = wi_dx */
		  const float ui = r * h_inv;
		  d_kernel_deval(ui, &wi, &wi_dx);

          const float delta_u_factor = (u - u_tmp[j_block]) * r_inv;
          laplace_u += mj * delta_u_factor * wi_dx / rho_tmp[j_block];

          /* Set the maximal alpha from the previous step over the neighbours
           * (this is used to limit the diffusion in hydro_prepare_force) */
          const float alpha_j = alpha_tmp[j_block];
          alpha_visc_max_ngb = max(alpha_visc_max_ngb, alpha_j);
		}
	  } /*if (jj < cj_end && pid < ci_end && pid >= ci_start)*/
	} /*End of looping through particles in shared memory---Shared arrays zero'ed for next step in outer loop*/
	__syncthreads();
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
  if (pid >= ci_start && pid < ci_end) {
    parts_aos[pid].v_sig = v_sig, parts_aos[pid].laplace_u = laplace_u;
    parts_aos[pid].alpha_visc_max_ngb = alpha_visc_max_ngb;
  }
}
#ifdef WITH_CUDA
}
#endif


#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NAIVEGPUAOSF4G(const struct part_aos_f4_g_send pi,
	struct part_aos_f4_g_send * __restrict__ parts_send, struct part_aos_f4_g_recv * __restrict__ parts_recv, int pid,
	const int cj_start, const int cj_end, float d_a, float d_H) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

  float hi = 0.0, hig2 = 0.0;

  int Found_neighbours = 0;
  int count_i = cj_start;

  float4 res_rho = {0.0, 0.0, 0.0, 0.0};
  float4 res_rot = {0.0, 0.0, 0.0, 0.0};
//  const part_aos_f4_send pi = parts_send[pid];
  const float4 x_h_i = pi.x_h;
  const float4 ux_m_i = pi.ux_m;
  const float4 rho_avisc_u_c_i = pi.rho_avisc_u_c;
  float3 vsig_lapu_aviscmax_i = {0.f, 0.f, 0.f};

//  printf("first_part_in_task_blocks_cj %i last_part_in_task_blocks_cj %i last_part_in_task_blocks_ci %i\n",
//		  first_part_in_task_blocks_cj, last_part_in_task_blocks_cj, last_part_in_task_blocks_ci);
//  if (pid < ci_end) {
	hi = x_h_i.w, hig2 = hi * hi * kernel_gamma2;
//  }

//  printf("js %i je %i\n", cj_start, cj_end);
  /*Particles copied in blocks to shared memory*/
  for (int j = cj_start; j < cj_end; j ++) {
    struct part_aos_f4_g_send pj = parts_send[j];

	const float4 x_h_j = pj.x_h;
	const float4 ux_m_j = pj.ux_m;
    const float4 rho_avisc_u_c_j = pj.rho_avisc_u_c;
	const float xij = x_h_i.x - x_h_j.x, yij = x_h_i.y - x_h_j.y, zij = x_h_i.z - x_h_j.z;
	const float r2 = xij * xij + yij * yij + zij * zij;
//	printf("r2 %f \n", r2);
	if (r2 < hig2) {
        const float r = sqrt(r2);
        const float r_inv = 1.f / r;
        /* Recover some data */
        const float mj = ux_m_j.w;
        /* Get the kernel for hi. */
        const float h_inv = 1.f / hi;
        float wi, wi_dx;
        /* Cosmology terms for the signal velocity */
        const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
        const float a2_Hubble = d_a * d_a * d_H;
        /* Compute dv dot r */
        float dvx = ux_m_i.x - ux_m_j.x, dvy = ux_m_i.y - ux_m_j.y,
              dvz = ux_m_i.z - ux_m_j.z;
        const float dvdr = dvx * xij + dvy * yij + dvz * zij;
        /* Add Hubble flow */
        const float dvdr_Hubble = dvdr + a2_Hubble * r2;
        /* Are the particles moving towards each others ? */
        const float omega_ij = min(dvdr_Hubble, 0.f);
        const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */

        /* Signal velocity */
        const float new_v_sig = rho_avisc_u_c_i.w + rho_avisc_u_c_j.w - const_viscosity_beta * mu_ij;
        /* Update if we need to */
        vsig_lapu_aviscmax_i.x = fmaxf(vsig_lapu_aviscmax_i.x, new_v_sig);
        /* Calculate Del^2 u for the thermal diffusion coefficient. */
        /* Need to get some kernel values F_ij = wi_dx */
        const float ui = r * h_inv;
        d_kernel_deval(ui, &wi, &wi_dx);

        const float delta_u_factor = (rho_avisc_u_c_i.z - rho_avisc_u_c_j.z) * r_inv;
        vsig_lapu_aviscmax_i.y += mj * delta_u_factor * wi_dx / rho_avisc_u_c_j.x;

        /* Set the maximal alpha from the previous step over the neighbours
         * (this is used to limit the diffusion in hydro_prepare_force) */
        const float alpha_j = rho_avisc_u_c_j.y;
        vsig_lapu_aviscmax_i.z = fmaxf(vsig_lapu_aviscmax_i.z, alpha_j);
	}
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
//  if (pid >= ci_start && pid < ci_end) {
  parts_recv[pid].vsig_lapu_aviscmax = vsig_lapu_aviscmax_i;
//  }
}
#ifdef WITH_CUDA
}
#endif


#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NONSYMGPUAOSF(
	struct part_aos_f *parts_aos, int pid, const int ci_start,
	const int ci_end, const int cj_start,
	const int cj_end, float d_a, float d_H,
	float *vars_pair_aosf, double *d_shift_x, double *d_shift_y, double *d_shift_z, const int task_id_tmp, int flip_order) {

  float ci = 0.0, cj = 0.0;
  float hi = 0.0, hig2 = 0.0;
  float mi = 0.0;
  float uxi = 0.0;
  float uyi = 0.0;
  float uzi = 0.0;
  float pix = 0.0;
  float piy = 0.0;
  float piz = 0.0;
  float rhoi = 0.0;
  float div_vi = 0.0;
  int Found_neighbours = 0;
  float v_sigi;
  float ui = 0.f;
  float u_dti = 0.f;
  float laplace_ui = 0.0;
  float alpha_visc_max_ngb = 0.0;
  float pressurei = 0.0;
  float alphavisci = 0.0;
  float alphadiffi = 0.0;
  float fi = 0.0;
  float balsarai = 0.0;
  float ahydroxi = 0.0;
  float ahydroyi = 0.0;
  float ahydrozi = 0.0;
  float h_dti = 0.0;
  int min_ngb_time_bin = 0;
  if (pid < ci_end) {
	hi = parts_aos[pid].h, hig2 = hi * hi * kernel_gamma2;
	mi = parts_aos[pid].mass;
	uxi = parts_aos[pid].ux;
	uyi = parts_aos[pid].uy;
	uzi = parts_aos[pid].uz;
	ci = parts_aos[pid].soundspeed;
	fi = parts_aos[pid].f;
	v_sigi = parts_aos[pid].v_sig;
	ui = parts_aos[pid].u;
	rhoi = parts_aos[pid].rho;
	pressurei = parts_aos[pid].pressure;
	balsarai = parts_aos[pid].balsara;
	alphavisci = parts_aos[pid].alpha_visc;
	alphadiffi = parts_aos[pid].alpha_diff;
	min_ngb_time_bin = parts_aos[pid].min_ngb_time_bin;
	pix = parts_aos[pid].x_p - d_shift_x[task_id_tmp];
	piy = parts_aos[pid].y_p - d_shift_y[task_id_tmp];
	piz = parts_aos[pid].z_p - d_shift_z[task_id_tmp];
  }
//  if (threadIdx.x == 0) {
//    first_part_tid_0 = first_part;
//    last_part_tid_0 = last_part;
//  }
//  __syncthreads();
  int n_neighbours = 0;
  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  float *x_p_tmp = (float *)&vars_pair_aosf[0];
  float *y_p_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE];
  float *z_p_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 2];
  float *h_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 3];
  float *mass_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 4];
  float *ux_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 5];
  float *uy_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 6];
  float *uz_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 7];
  float *cj_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 8];
  float *alphavisc_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 9];
  float *alphadiff_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 10];
  float *u_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 11];
  float *rho_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 12];
  float *pressure_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 13];
  float *f_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 14];
  float *balsara_tmp = (float *)&vars_pair_aosf[BLOCK_SIZE * 15];
  int *timebin = (int *)&vars_pair_aosf[BLOCK_SIZE * 16];
  /*Particles copied in blocks to shared memory*/
  for (int b = cj_start; b < cj_end;
	   b += BLOCK_SIZE) {
	int j = b + threadIdx.x;
	x_p_tmp[threadIdx.x] = parts_aos[j].x_p;
	y_p_tmp[threadIdx.x] = parts_aos[j].y_p;
	z_p_tmp[threadIdx.x] = parts_aos[j].z_p;
	h_tmp[threadIdx.x] = parts_aos[j].h;
	mass_tmp[threadIdx.x] = parts_aos[j].mass;
	ux_tmp[threadIdx.x] = parts_aos[j].ux;
	uy_tmp[threadIdx.x] = parts_aos[j].uy;
	uz_tmp[threadIdx.x] = parts_aos[j].uz;
	timebin[threadIdx.x] = parts_aos[j].time_bin;
	cj_tmp[threadIdx.x] = parts_aos[j].soundspeed;
//    alpha_tmp[threadIdx.x] = parts_aos[j].visc_alpha;
	u_tmp[threadIdx.x] = parts_aos[j].u;
	rho_tmp[threadIdx.x] = parts_aos[j].rho;
	alphavisc_tmp[threadIdx.x] = parts_aos[j].alpha_visc;
	alphadiff_tmp[threadIdx.x] = parts_aos[j].alpha_diff;
	pressure_tmp[threadIdx.x] = parts_aos[j].pressure;
	f_tmp[threadIdx.x] = parts_aos[j].f;
	balsara_tmp[threadIdx.x] = parts_aos[j].balsara;
	const float shift_x_j = d_shift_x[task_id_tmp + flip_order];
	const float shift_y_j = d_shift_y[task_id_tmp + flip_order];
	const float shift_z_j = d_shift_z[task_id_tmp + flip_order];
	__syncthreads();
	/*j_block is the particle's index in the block. Loop through particles in shared memory one by one*/
	for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
	  int jj = b + j_block;
	  if (jj < cj_end && pid < ci_end && pid >= ci_start) {
		/* Compute the pairwise distance. */
		const float pjx = x_p_tmp[j_block] - shift_x_j;
		const float pjy = y_p_tmp[j_block] - shift_y_j;
		const float pjz = z_p_tmp[j_block] - shift_z_j;
		const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
		const float r2 = xij * xij + yij * yij + zij * zij;
		if (r2 < hig2) {

		  //          /* Cosmology terms for the signal velocity */
		  const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
		  const float a2_Hubble = d_a * d_a * d_H;
		  const float r = sqrt(r2);
		  const float r_inv = 1.f / r;
//          /* Recover some data */
		  const float mj = mass_tmp[j_block];
//          /* Get the kernel for hi. */
		  const float hi_inv = 1.f / hi;
		  const float hid_inv = d_pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
		  const float xi = r * hi_inv;
		  float wi, wi_dx;
		  d_kernel_deval(xi, &wi, &wi_dx);
		  const float wi_dr = hid_inv * wi_dx;
		  /* Get the kernel for hj. */
		  const float hj = h_tmp[j_block];
		  const float hj_inv = 1.0f / hj;
		  const float hjd_inv = d_pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
		  const float xj = r * hj_inv;
		  float wj, wj_dx;
		  d_kernel_deval(xj, &wj, &wj_dx);
		  const float wj_dr = hjd_inv * wj_dx;
//          /* Compute dv dot r */
		  float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
				dvz = uzi - uz_tmp[j_block];
		  const float dvdr = dvx * xij + dvy * yij + dvz * zij;
//          /* Add Hubble flow */
		  const float dvdr_Hubble = dvdr + a2_Hubble * r2;
//          /* Are the particles moving towards each others ? */
		  const float omega_ij = min(dvdr_Hubble, 0.f);
		  const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */
//
//          /* Signal velocity */
		  const float v_sig = ci + cj_tmp[j_block] - const_viscosity_beta * mu_ij;

		  /* Variable smoothing length term */
		  const float f_ij = 1.f - fi / mj;
		  const float f_ji = 1.f - f_tmp[j_block] / mi;

		  /* Balsara term */
		  const float balsaraj = balsara_tmp[j_block];
		  /* Construct the full viscosity term */
		  const float rhoj = rho_tmp[j_block];
		  const float pressurej = pressure_tmp[j_block];
		  const float rho_ij = rhoi + rhoj;
		  const float alpha = alphavisci + alphavisc_tmp[j_block];
		  const float visc =
			  -0.25f * alpha * v_sig * mu_ij * (balsarai + balsaraj) / rho_ij;
		  /* Convolve with the kernel */
		  const float visc_acc_term =
			  0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;
		  /* Compute gradient terms */
		  const float P_over_rho2_i = pressurei / (rhoi * rhoi) * f_ij;
		  const float P_over_rho2_j = pressurej / (rhoj * rhoj) * f_ji;

		  /* SPH acceleration term */
		  const float sph_acc_term =
			  (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

		  /* Assemble the acceleration */
		  const float acc = sph_acc_term + visc_acc_term;
		  /* Use the force Luke ! */
		  ahydroxi -= mj * acc * xij;
		  ahydroyi -= mj * acc * yij;
		  ahydrozi -= mj * acc * zij;
//          if(rhoi == 0 || rhoj == 0 || pressurei == 0 || pressurej == 0)printf("ri %f rj %f pi %f pj %f\n", rhoi, rhoj, pressurei, pressurej);
		  /* Get the time derivative for u. */
		  const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

		  /* Viscosity term */
		  const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;
		  const float press_sum = pressurei + pressurej;
		  /* Diffusion term */
		  /* Combine the alpha_diff into a pressure-based switch -- this allows the
		   * alpha from the highest pressure particle to dominate, so that the
		   * diffusion limited particles always take precedence - another trick to
		   * allow the scheme to work with thermal feedback. */
		  float alpha_diff =
			  (pressurei * alphadiffi + pressurej * alphadiff_tmp[j_block]) /
			  (press_sum);
		  if (fabsf(press_sum) < 1e-10) alpha_diff = 0.f;
		  const float v_diff = alpha_diff * 0.5f *
							   (sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij) +
								fabsf(fac_mu * r_inv * dvdr_Hubble));
		  /* wi_dx + wj_dx / 2 is F_ij */
		  const float diff_du_term =
			  v_diff * (ui - u_tmp[j_block]) * (f_ij * wi_dr / rhoi + f_ji * wj_dr / rhoj);

		  /* Assemble the energy equation term */
		  const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

		  /* Internal energy time derivative */
		  u_dti += du_dt_i * mj;
		  if(mj == 0.f)printf("zero mass mj %f\n", mj);

		  /* Get the time derivative for h. */
		  h_dti -= mj * dvdr * r_inv / rhoj * wi_dr;

		  /* Update if we need to; this should be guaranteed by the gradient loop but
		   * due to some possible synchronisation problems this is here as a _quick
		   * fix_. Added: 14th August 2019. To be removed by 1st Jan 2020. (JB) */
		  v_sigi = max(v_sigi, v_sig);
		  int time_bin_j = timebin[j_block];
		  if(time_bin_j > 0)min_ngb_time_bin = min(min_ngb_time_bin, time_bin_j);
//          printf("Got in\n");
		}
	  }
	}
	__syncthreads();
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
  if (pid >= ci_start && pid < ci_end) {
	parts_aos[pid].v_sig = v_sigi;
	parts_aos[pid].h_dt = h_dti;
	parts_aos[pid].u_dt = u_dti;
	parts_aos[pid].a_hydrox = ahydroxi;
	parts_aos[pid].a_hydroy = ahydroyi;
	parts_aos[pid].a_hydroz = ahydrozi;
	parts_aos[pid].min_ngb_time_bin = min_ngb_time_bin;
//    printf("%f %f %f %f %f %f\n", v_sigi, h_dti, u_dti, ahydroxi, ahydroyi, ahydrozi);
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2NAIVEGPUAOSF4F(const struct part_aos_f4_f_send pi,
	struct part_aos_f4_f_send * __restrict__ parts_send, struct part_aos_f4_f_recv * __restrict__ parts_recv, int pid,
	const int cj_start, const int cj_end, float d_a, float d_H) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

  int Found_neighbours = 0;

//  const part_aos_f4_send pi = parts_send[pid];
  const float4 x_h_i = pi.x_h;
  const float4 ux_m_i = pi.ux_m;

  float4 f_b_t_mintbinngb_i = pi.f_bals_timebin_mintimebin_ngb;
  const float4 rho_p_c_vsig_i = pi.rho_p_c_vsigi;
  const float3 u_avisc_adiff_i = pi.u_alphavisc_alphadiff;

  const float mi = ux_m_i.w;
  const float pressurei = rho_p_c_vsig_i.y;
  const float ci = rho_p_c_vsig_i.z;
  float3 ahydro = {0.0, 0.0, 0.0};
  float4 udt_hdt_vsig_mintbinngb = {0.0, 0.0, 0.0, 0.0};
  udt_hdt_vsig_mintbinngb.z = rho_p_c_vsig_i.w;
  udt_hdt_vsig_mintbinngb.w = f_b_t_mintbinngb_i.w;

  const float hi = x_h_i.w;
  const float hig2 = hi * hi * kernel_gamma2;

//  printf("js %i je %i\n", cj_start, cj_end);
  /*Particles copied in blocks to shared memory*/
  for (int j = cj_start; j < cj_end; j ++) {
    struct part_aos_f4_f_send pj = parts_send[j];
	const float4 x_h_j = pj.x_h;
	const float4 ux_m_j = pj.ux_m;
    const float4 f_b_t_mintbinngb_j = pj.f_bals_timebin_mintimebin_ngb;
    const float4 rho_p_c_vsig_j = pj.rho_p_c_vsigi;
//    alpha_tmp[threadIdx.x] = parts_aos[j].visc_alpha;
    const float3 u_avisc_adiff_j = pj.u_alphavisc_alphadiff;
	const float xij = x_h_i.x - x_h_j.x,
			    yij = x_h_i.y - x_h_j.y,
			    zij = x_h_i.z - x_h_j.z;
	const float r2 = xij * xij + yij * yij + zij * zij;
//	printf("r2 %f \n", r2);
	if (r2 < hig2) {
        //          /* Cosmology terms for the signal velocity */
        const float fac_mu = d_pow_three_gamma_minus_five_over_two(d_a);
        const float a2_Hubble = d_a * d_a * d_H;
        const float r = sqrt(r2);
        const float r_inv = 1.f / r;
//          /* Recover some data */
        const float mj = ux_m_j.w;
//          /* Get the kernel for hi. */
        const float hi_inv = 1.f / hi;
        const float hid_inv = d_pow_dimension_plus_one(hi_inv); /* 1/h^(d+1) */
        const float xi = r * hi_inv;
        float wi, wi_dx;
        d_kernel_deval(xi, &wi, &wi_dx);
        const float wi_dr = hid_inv * wi_dx;
        /* Get the kernel for hj. */
        const float hj = x_h_j.w;
        const float hj_inv = 1.0f / hj;
        const float hjd_inv = d_pow_dimension_plus_one(hj_inv); /* 1/h^(d+1) */
        const float xj = r * hj_inv;
        float wj, wj_dx;
        d_kernel_deval(xj, &wj, &wj_dx);
        const float wj_dr = hjd_inv * wj_dx;
//          /* Compute dv dot r */
        float dvx = ux_m_i.x - ux_m_j.x, dvy = ux_m_i.y - ux_m_j.y,
              dvz = ux_m_i.z - ux_m_j.z;
        const float dvdr = dvx * xij + dvy * yij + dvz * zij;
//          /* Add Hubble flow */
        const float dvdr_Hubble = dvdr + a2_Hubble * r2;
//          /* Are the particles moving towards each others ? */
        const float omega_ij = min(dvdr_Hubble, 0.f);
        const float mu_ij = fac_mu * r_inv * omega_ij; /* This is 0 or negative */
//
//          /* Signal velocity */
        const float cj = rho_p_c_vsig_j.z;
        const float v_sig = ci + cj - const_viscosity_beta * mu_ij;

        /* Variable smoothing length term */
        const float f_ij = 1.f - f_b_t_mintbinngb_i.x / mj;
        const float f_ji = 1.f - f_b_t_mintbinngb_j.x / mi;

        /* Construct the full viscosity term */
        const float pressurej = rho_p_c_vsig_j.y;
        const float rho_ij = rho_p_c_vsig_i.x + rho_p_c_vsig_j.x;
        const float alpha = u_avisc_adiff_i.y + u_avisc_adiff_j.y;
        const float visc =
            -0.25f * alpha * v_sig * mu_ij * (f_b_t_mintbinngb_i.y + f_b_t_mintbinngb_j.y) / rho_ij;
        /* Convolve with the kernel */
        const float visc_acc_term =
            0.5f * visc * (wi_dr * f_ij + wj_dr * f_ji) * r_inv;
        /* Compute gradient terms */
        const float rhoi2 = rho_p_c_vsig_i.x * rho_p_c_vsig_i.x;
        const float rhoj2 = rho_p_c_vsig_j.x * rho_p_c_vsig_j.x;
        const float P_over_rho2_i = pressurei / (rhoi2) * f_ij;
        const float P_over_rho2_j = pressurej / (rhoj2) * f_ji;

        /* SPH acceleration term */
        const float sph_acc_term =
            (P_over_rho2_i * wi_dr + P_over_rho2_j * wj_dr) * r_inv;

        /* Assemble the acceleration */
        const float acc = sph_acc_term + visc_acc_term;
        /* Use the force Luke ! */
        ahydro.x -= mj * acc * xij;
        ahydro.y -= mj * acc * yij;
        ahydro.z -= mj * acc * zij;
        /* Get the time derivative for u. */
        const float sph_du_term_i = P_over_rho2_i * dvdr * r_inv * wi_dr;

        /* Viscosity term */
        const float visc_du_term = 0.5f * visc_acc_term * dvdr_Hubble;
        /* Diffusion term */
        /* Combine the alpha_diff into a pressure-based switch -- this allows the
         * alpha from the highest pressure particle to dominate, so that the
         * diffusion limited particles always take precedence - another trick to
         * allow the scheme to work with thermal feedback. */
        float alpha_diff =
            (pressurei * u_avisc_adiff_i.z + pressurej * u_avisc_adiff_j.z) /
            (pressurei + pressurej);
        if (fabsf(pressurei + pressurej) < 1e-10) alpha_diff = 0.f;
        const float v_diff = alpha_diff * 0.5f *
                             (sqrtf(2.f * fabsf(pressurei - pressurej) / rho_ij) +
                              fabsf(fac_mu * r_inv * dvdr_Hubble));
        /* wi_dx + wj_dx / 2 is F_ij */
        const float diff_du_term =
            v_diff * (u_avisc_adiff_i.x - u_avisc_adiff_j.x) *
			  (f_ij * wi_dr / rho_p_c_vsig_i.x + f_ji * wj_dr / rho_p_c_vsig_j.x);

        /* Assemble the energy equation term */
        const float du_dt_i = sph_du_term_i + visc_du_term + diff_du_term;

        /* Internal energy time derivative */
        udt_hdt_vsig_mintbinngb.x += du_dt_i * mj;

        /* Get the time derivative for h. */
        udt_hdt_vsig_mintbinngb.y -= mj * dvdr * r_inv / rho_p_c_vsig_j.x * wi_dr;

        /* Update if we need to; this should be guaranteed by the gradient loop but
         * due to some possible synchronisation problems this is here as a _quick
         * fix_. Added: 14th August 2019. To be removed by 1st Jan 2020. (JB) */
        udt_hdt_vsig_mintbinngb.z = fmaxf(udt_hdt_vsig_mintbinngb.z, v_sig);
        unsigned int time_bin_j = (f_b_t_mintbinngb_j.z + 0.5f);
        unsigned int min_tb_i = (f_b_t_mintbinngb_i.w + 0.5f);
        if(time_bin_j > 0)f_b_t_mintbinngb_i.w =
      		  min(min_tb_i, time_bin_j);
//          printf("Got in\n");
	}
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
//  if (pid >= ci_start && pid < ci_end) {
  udt_hdt_vsig_mintbinngb.w = f_b_t_mintbinngb_i.w;
  parts_recv[pid].udt_hdt_vsig_mintimebin_ngb = udt_hdt_vsig_mintbinngb;
  parts_recv[pid].a_hydro = ahydro;
//  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__device__ void DOPAIR2GPU(
	struct part_soa parts_soa, int pid, const int ci_start,
	const int ci_end, const int cj_start,
	const int cj_end, float d_a, float d_H,
	int time_bin_inhibited, float *vars_pair, double *d_shift_x, double *d_shift_y, double *d_shift_z, const int task_id_tmp) {

  float dx = 1.f/64.f; //Value used to avoid interacting parts with themselves

  float cellx = 0.0, celly = 0.0, cellz = 0.0;
  float cellxj = 0.0, cellyj = 0.0, cellzj = 0.0;
  float hi = 0.0, hig2 = hi * hi * kernel_gamma2;

  float mi = 0.0;
  float uxi = 0.0;
  float uyi = 0.0;
  float uzi = 0.0;
  double pix = 0.0;
  double piy = 0.0;
  double piz = 0.0;
  float rhoi = 0.0;
  float rho_dhi = 0.0;
  float wcounti = 0.0;
  float wcount_dhi = 0.0;
  float div_vi = 0.0;
  float rot_uxi = 0.0;
  float rot_uyi = 0.0;
  float rot_uzi = 0.0;
  int Found_neighbours = 0;
  int count_i = cj_start;
//  printf("first_part_in_task_blocks_cj %i last_part_in_task_blocks_cj %i last_part_in_task_blocks_ci %i\n",
//		  first_part_in_task_blocks_cj, last_part_in_task_blocks_cj, last_part_in_task_blocks_ci);

  if (pid < ci_end) {
	cellx = parts_soa.locx[pid];
	celly = parts_soa.locy[pid];
	cellz = parts_soa.locz[pid];
	const int j = cj_start;
	cellxj = parts_soa.locx[j];
	cellyj = parts_soa.locy[j];
	cellzj = parts_soa.locz[j];
	hi = parts_soa.h[pid], hig2 = hi * hi * kernel_gamma2;
	mi = parts_soa.mass[pid];
	uxi = parts_soa.ux[pid];
	uyi = parts_soa.uy[pid];
	uzi = parts_soa.uz[pid];
	pix = parts_soa.x_p[pid] - d_shift_x[task_id_tmp];
	piy = parts_soa.y_p[pid] - d_shift_y[task_id_tmp];
	piz = parts_soa.z_p[pid] - d_shift_z[task_id_tmp];
  }

  int n_neighbours = 0;
  float av_dist = 0.f;
  float av_distx = 0.f;
  float av_disty = 0.f;
  float av_distz = 0.f;
  float distby2h = 0.f;
  /*Here we use different pointers "x_p_tmp", etc. to point to different regions
   * of the single shared memory space "vars" which we allocate in kernel
   * invocation*/
  double *x_p_tmp = (double *)&vars_pair[0];
  double *y_p_tmp = (double *)&x_p_tmp[BLOCK_SIZE];
  double *z_p_tmp = (double *)&y_p_tmp[BLOCK_SIZE];
  float *h_tmp = (float *)&z_p_tmp[BLOCK_SIZE];
  float *mass_tmp = (float *)&h_tmp[BLOCK_SIZE];
  float *ux_tmp = (float *)&mass_tmp[BLOCK_SIZE];
  float *uy_tmp = (float *)&ux_tmp[BLOCK_SIZE];
  float *uz_tmp = (float *)&uy_tmp[BLOCK_SIZE];
  timebin_t *timebin = (timebin_t *)&uz_tmp[BLOCK_SIZE];
  float *rho_tmp = (float *)&timebin[BLOCK_SIZE];
  float *rho_dh_tmp = (float *)&rho_tmp[BLOCK_SIZE];
  float *wcount_tmp = (float *)&rho_dh_tmp[BLOCK_SIZE];
  float *wcount_dh_tmp = (float *)&wcount_tmp[BLOCK_SIZE];
  float *div_v_tmp = (float *)&wcount_dh_tmp[BLOCK_SIZE];
  float *rot_ux_tmp = (float *)&div_v_tmp[BLOCK_SIZE];
  float *rot_uy_tmp = (float *)&rot_ux_tmp[BLOCK_SIZE];
  float *rot_uz_tmp = (float *)&rot_uy_tmp[BLOCK_SIZE];

  /*Particles copied in blocks to shared memory*/
  for (int b = cj_start; b < cj_end;
	   b += BLOCK_SIZE) {
	const int tid_x = threadIdx.x;
	int j = b + tid_x;
	x_p_tmp[tid_x] = parts_soa.x_p[j];
	y_p_tmp[tid_x] = parts_soa.y_p[j];
	z_p_tmp[tid_x] = parts_soa.z_p[j];
	h_tmp[tid_x] = parts_soa.h[j];
	mass_tmp[tid_x] = parts_soa.mass[j];
	ux_tmp[tid_x] = parts_soa.ux[j];
	uy_tmp[tid_x] = parts_soa.uy[j];
	uz_tmp[tid_x] = parts_soa.uz[j];
	timebin[tid_x] = parts_soa.time_bin[j];
	rho_tmp[tid_x] = 0.f;
	rho_dh_tmp[tid_x] = 0.f;
	wcount_tmp[tid_x] = 0.f;
	wcount_dh_tmp[tid_x] = 0.f;
	div_v_tmp[tid_x] = 0.f;
	rot_ux_tmp[tid_x] = 0.f;
	rot_uy_tmp[tid_x] = 0.f;
	rot_uz_tmp[tid_x] = 0.f;
	__syncthreads();
	const double shift_x_j = d_shift_x[task_id_tmp + 1];
	const double shift_y_j = d_shift_y[task_id_tmp + 1];
	const double shift_z_j = d_shift_z[task_id_tmp + 1];
	/*j_block is the particle's index in the block. Loop through particles in shared memory one by one*/
	for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
	  int jj = b + j_block;
	  if (jj < cj_end && pid < ci_end && pid >= ci_start) {

		const double pjx = x_p_tmp[j_block] - shift_x_j;
		const double pjy = y_p_tmp[j_block] - shift_y_j;
		const double pjz = z_p_tmp[j_block] - shift_z_j;

		const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
//		const float xij = pjx - pix, yij = pjy - piy, zij = pjz - piz;
		const float r2 = xij * xij + yij * yij + zij * zij;
		const float hj =  h_tmp[j_block];
		const float hjg2 = hj * hj * kernel_gamma2;
//		if(r2 > 32.f * hig2 && hig2 != 0.f) printf("x %f y %f z %f r %f hig2 %f\n", xij/dx, yij/dx, zij/dx, sqrt(r2)/dx);
		/* Compute dv dot r */
		const float dvx = uxi - ux_tmp[j_block], dvy = uyi - uy_tmp[j_block],
		dvz = uzi - uz_tmp[j_block];
		const float dvdr = dvx * xij + dvy * yij + dvz * zij;
		/* Compute dv cross r */
		const float curlvrx = dvy * zij - dvz * yij;
		const float curlvry = dvz * xij - dvx * zij;
		const float curlvrz = dvx * yij - dvy * xij;

		const float r = sqrt(r2);
		if (r2 < hig2) {
		  /* Recover some data */
		  const float mj = mass_tmp[j_block];
		  /* Get the kernel for hi. */
//		  if(hi<1.f/dx)printf("h < dx\n");
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

		  div_vi -= faci * dvdr;

		  rot_uxi += faci * curlvrx;
		  rot_uyi += faci * curlvry;
		  rot_uzi += faci * curlvrz;
//
		}
		if (r2 < hjg2) {
		  /* Recover some data */
		  /* Get the kernel for hi. */
		  const float hj_inv = 1.f / hj;
		  const float uj = r * hj_inv;
		  float wj, wj_dx;

		  d_kernel_deval(uj, &wj, &wj_dx);

//		  atomicAdd(&rho_tmp[j_block], mi * wj);
		  atomicAdd(&parts_soa.rho[j], mi * wj);
//		  atomicAdd(&rho_dh_tmp[j_block], -mi * (hydro_dimension * wj + uj * wj_dx));
		  atomicAdd(&parts_soa.rho_dh[j], -mi * (hydro_dimension * wj + uj * wj_dx));

//		  atomicAdd(&wcount_tmp[j_block], wj);
		  atomicAdd(&parts_soa.wcount[j], wj);
//		  atomicAdd(&wcount_dh_tmp[j_block], -(hydro_dimension * wj + uj * wj_dx));
		  atomicAdd(&parts_soa.wcount_dh[j], -(hydro_dimension * wj + uj * wj_dx));

		  const float r_inv = 1.f / r;
		  const float facj = mi * wj_dx * r_inv;

//		  atomicAdd(&div_v_tmp[j_block], -facj * dvdr);
		  atomicAdd(&parts_soa.div_v[j], -facj * dvdr);

//		  atomicAdd(&rot_ux_tmp[j_block], facj * curlvrx);
//		  atomicAdd(&rot_uy_tmp[j_block], facj * curlvry);
//		  atomicAdd(&rot_uz_tmp[j_block], facj * curlvrz);
		  atomicAdd(&parts_soa.rot_ux[j], facj * curlvrx);
		  atomicAdd(&parts_soa.rot_uy[j], facj * curlvry);
		  atomicAdd(&parts_soa.rot_uz[j], facj * curlvrz);
//		  printf("rho %f rho_dh %f wcount %f wcount_dh %f div_v %f rotux %f rotuy %f rotuz %f\n"
//				 ,rhoi, rho_dhi, wcounti, wcount_dhi, div_vi, rot_uxi, rot_uyi, rot_uzi);
		} /*if r2<hjg2 */
	  } /*if (jj < cj_end && pid < ci_end && pid >= ci_start)*/
	} /*End of looping through particles in shared memory---Shared arrays zero'ed for next step in outer loop*/
	__syncthreads();
//	if(j < cj_end){
//	  atomicAdd(&parts_soa.rho[j], rho_tmp[threadIdx.x]);
//	  atomicAdd(&parts_soa.rho_dh[j], rho_dh_tmp[threadIdx.x]);
//	  atomicAdd(&parts_soa.wcount[j], wcount_tmp[threadIdx.x]);
//	  atomicAdd(&parts_soa.wcount_dh[j], wcount_dh_tmp[threadIdx.x]);
//	  atomicAdd(&parts_soa.div_v[j], div_v_tmp[threadIdx.x]);
//	  atomicAdd(&parts_soa.rot_ux[j], rot_ux_tmp[threadIdx.x]);
//	  atomicAdd(&parts_soa.rot_uy[j], rot_uy_tmp[threadIdx.x]);
//	  atomicAdd(&parts_soa.rot_uz[j], rot_uz_tmp[threadIdx.x]);
//	}
//	__syncthreads();
//	parts_soa.rho[j] += rho_tmp[threadIdx.x];
//	parts_soa.rho_dh[j] += rho_dh_tmp[threadIdx.x];
//	parts_soa.wcount[j] += wcount_tmp[threadIdx.x];
//	parts_soa.wcount_dh[j] =+ wcount_dh_tmp[threadIdx.x];
//	parts_soa.div_v[j] += div_v_tmp[threadIdx.x];
//	parts_soa.rot_ux[j] += rot_ux_tmp[threadIdx.x];
//	parts_soa.rot_uy[j] =+ rot_uy_tmp[threadIdx.x];
//	parts_soa.rot_uz[j] += rot_uz_tmp[threadIdx.x];
  } /*Loop through parts in cell j one BLOCK_SIZE at a time*/
  if (pid >= ci_start && pid < ci_end) {
//	if(n_neighbours > 0){
//		distby2h = distby2h/n_neighbours;
//		av_dist = av_dist/(n_neighbours*dx);
//	}
//    av_distx = av_distx/(n_neighbours*dx);
//    av_disty = av_disty/(n_neighbours*dx);
//    av_distz = av_distz/(n_neighbours*dx);
	parts_soa.rho[pid] = rhoi, parts_soa.rho_dh[pid] = rho_dhi;
	parts_soa.wcount[pid] = wcounti, parts_soa.wcount_dh[pid] = wcount_dhi;
	parts_soa.div_v[pid] = div_vi;
	parts_soa.rot_ux[pid] = rot_uxi, parts_soa.rot_uy[pid] = rot_uyi;
	parts_soa.rot_uz[pid] = rot_uzi;
//	if(rhoi != 0.f)printf("rho i %f, rho_dh i %f\n", rhoi, rho_dhi);
  }

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_density_GPU(
	struct part_soa parts_soa, int *d_task_first_part_ci, int *d_task_first_part_cj,
	int *d_task_last_part_ci, int *d_task_last_part_cj, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, int time_bin_inhibited) {

  extern __shared__ float vars[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;

  first_part_in_task_blocks_ci = d_task_first_part_ci[task_id];
  last_part_in_task_blocks_ci = d_task_last_part_ci[task_id];
  first_part_in_task_blocks_cj = d_task_first_part_cj[task_id];
  last_part_in_task_blocks_cj = d_task_last_part_cj[task_id];

  // Now we start calculations for particles in cell i
  const int pid = threadid + first_part_in_task_blocks_ci;

  /*Don't ever put me in an if statement. I've got __syncthreads inside*/
  DOPAIRGPU(
  		parts_soa, pid, last_part_in_task_blocks_ci,
		first_part_in_task_blocks_cj, last_part_in_task_blocks_cj, d_a, d_H,
		time_bin_inhibited, vars);
//  __syncthreads();
  // Now we start calculations for particles in cell i
  const int pjd = threadid + last_part_in_task_blocks_ci;
  /*Don't ever put me in an if statement. I've got __syncthreads inside*/
  DOPAIRGPU(
  		parts_soa, pjd, last_part_in_task_blocks_cj,
  		first_part_in_task_blocks_ci, last_part_in_task_blocks_ci, d_a, d_H,
  		time_bin_inhibited, vars);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_sym_density_GPU(
	struct part_soa parts_soa, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, int time_bin_inhibited, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  // Now we start calculations for particles in cell i
  const int pid = threadid + ci_start;

  /*Don't ever put me in an if statement. I've got __syncthreads inside*/
  DOPAIR2GPU(
  		parts_soa, pid, ci_start, ci_end,
		cj_start, cj_end, d_a, d_H,
		time_bin_inhibited, vars_pair, d_shift_x, d_shift_y, d_shift_z, task_id_tmp);
//  __syncthreads();

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_nonsym_density_GPU(
	struct part_soa parts_soa, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, int time_bin_inhibited, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /* Start calculations for particles in cell i
  * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pid = threadid + ci_start;
  const int flip_i = 1;
  DOPAIR2NONSYMGPU(
  		parts_soa, pid, ci_start, ci_end,
		cj_start, cj_end, d_a, d_H,
		vars_pair, d_shift_x, d_shift_y, d_shift_z, task_id_tmp, flip_i);

  /*Necessary evil to stop parts from j and i co-existing on shared memory for sums*/
  __syncthreads();

  /*Now do cj
   * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pjd = threadid + cj_start;
  const int flip_j = -1;
  DOPAIR2NONSYMGPU(
  		parts_soa, pjd, cj_start, cj_end,
		ci_start, ci_end, d_a, d_H,
		vars_pair, d_shift_x, d_shift_y, d_shift_z, task_id_tmp + 1, flip_j);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_ci_density_GPU(
	struct part_soa parts_soa, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /* Start calculations for particles in cell i
  * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pid = threadid + ci_start;
  const int flip_i = 1;
  DOPAIR2NONSYMGPU(
  		parts_soa, pid, ci_start, ci_end,
		cj_start, cj_end, d_a, d_H,
		vars_pair, d_shift_x, d_shift_y, d_shift_z, task_id_tmp, flip_i);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_cj_density_GPU(
	struct part_soa parts_soa, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /*Now do cj
   * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pjd = threadid + cj_start;
  const int flip_j = -1;
  DOPAIR2NONSYMGPU(
  		parts_soa, pjd, cj_start, cj_end,
		ci_start, ci_end, d_a, d_H,
		vars_pair, d_shift_x, d_shift_y, d_shift_z, task_id_tmp + 1, flip_j);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_ci_density_GPU_aos(
	struct part_aos *parts_aos, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair_aos[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /* Start calculations for particles in cell i
  * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pid = threadid + ci_start;
  const int flip_i = 1;
  DOPAIR2NONSYMGPUAOS(
  		parts_aos, pid, ci_start, ci_end,
		cj_start, cj_end, d_a, d_H,
		vars_pair_aos, d_shift_x, d_shift_y, d_shift_z, task_id_tmp, flip_i);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_cj_density_GPU_aos(
	struct part_aos *parts_aos, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair_aos[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /*Now do cj
   * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pjd = threadid + cj_start;
  const int flip_j = -1;
  DOPAIR2NONSYMGPUAOS(
  		parts_aos, pjd, cj_start, cj_end,
		ci_start, ci_end, d_a, d_H,
		vars_pair_aos, d_shift_x, d_shift_y, d_shift_z, task_id_tmp + 1, flip_j);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_ci_density_GPU_aos_f4(
	struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv, int4 *fparti_fpartj_lparti_lpartj_dens,
	float d_a, float d_H, int bundle_first_task) {

  extern __shared__ float4 vars_pair_i_f4[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int ci_start = fparti_fpartj_lparti_lpartj_dens[task_id].x;
  const int cj_start = fparti_fpartj_lparti_lpartj_dens[task_id].y;
  const int ci_end = fparti_fpartj_lparti_lpartj_dens[task_id].z;
  const int cj_end = fparti_fpartj_lparti_lpartj_dens[task_id].w;

  /* Start calculations for particles in cell i
  * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pid = threadid + ci_start;

  DOPAIR2NONSYMGPUAOSF4(parts_send, parts_recv, pid, ci_start, ci_end, cj_start, cj_end, d_a, d_H, vars_pair_i_f4);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_cj_density_GPU_aos_f4(
		struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv, int4 *fparti_fpartj_lparti_lpartj_dens,
			float d_a, float d_H, int bundle_first_task) {

  extern __shared__ float4 vars_pair_j_f4[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  const int ci_start = fparti_fpartj_lparti_lpartj_dens[task_id].x;
  const int cj_start = fparti_fpartj_lparti_lpartj_dens[task_id].y;
  const int ci_end = fparti_fpartj_lparti_lpartj_dens[task_id].z;
  const int cj_end =fparti_fpartj_lparti_lpartj_dens[task_id].w;

  /*Now do cj
   * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pjd = threadid + cj_start;
  DOPAIR2NONSYMGPUAOSF4(parts_send, parts_recv, pjd, cj_start, cj_end, ci_start, ci_end, d_a, d_H, vars_pair_j_f4);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_density_GPU_aos_f4(
	struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
	float d_a, float d_H, int bundle_first_part, int bundle_n_parts) {

//  extern __shared__ float4 vars_pair_i_f4[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int pid = bundle_first_part + threadid;
//  const int task_id = bundle_first_part + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  if(pid < bundle_first_part + bundle_n_parts){
    const struct part_aos_f4_send pi = parts_send[pid];
    const int cj_start = pi.cjs_cje.x;
    const int cj_end = pi.cjs_cje.y;

  /* Start calculations for particles in cell i*/
    DOPAIR2NAIVEGPUAOSF4(pi, parts_send, parts_recv, pid, cj_start, cj_end, d_a, d_H);
  }

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_ci_density_GPU_aos_g(
	struct part_aos_g *parts_aos, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair_aosg[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /* Start calculations for particles in cell i
  * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pid = threadid + ci_start;
  const int flip_i = 1;
  DOPAIR2NONSYMGPUAOSG(
  		parts_aos, pid, ci_start, ci_end,
		cj_start, cj_end, d_a, d_H,
		vars_pair_aosg, d_shift_x, d_shift_y, d_shift_z, task_id_tmp, flip_i);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_cj_density_GPU_aos_g(
	struct part_aos_g *parts_aos, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair_aosg[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /*Now do cj
   * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pjd = threadid + cj_start;
  const int flip_j = -1;
  DOPAIR2NONSYMGPUAOSG(
  		parts_aos, pjd, cj_start, cj_end,
		ci_start, ci_end, d_a, d_H,
		vars_pair_aosg, d_shift_x, d_shift_y, d_shift_z, task_id_tmp + 1, flip_j);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_gradient_GPU_aos_f4(
	struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv,
	float d_a, float d_H, int bundle_first_part, int bundle_n_parts) {

//  extern __shared__ float4 vars_pair_i_f4[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int pid = bundle_first_part + threadid;
//  const int task_id = bundle_first_part + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
  if(pid < bundle_first_part + bundle_n_parts){
    const struct part_aos_f4_g_send pi = parts_send[pid];
    const int cj_start = pi.cjs_cje.x;
    const int cj_end = pi.cjs_cje.y;
  /* Start calculations for particles in cell i*/
    DOPAIR2NAIVEGPUAOSF4G(pi, parts_send, parts_recv, pid, cj_start, cj_end, d_a, d_H);
  }

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_ci_density_GPU_aos_f(
	struct part_aos_f *parts_aos, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair_aosf[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /* Start calculations for particles in cell i
  * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pid = threadid + ci_start;
  const int flip_i = 1;
  DOPAIR2NONSYMGPUAOSF(
  		parts_aos, pid, ci_start, ci_end,
		cj_start, cj_end, d_a, d_H,
		vars_pair_aosf, d_shift_x, d_shift_y, d_shift_z, task_id_tmp, flip_i);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_cj_density_GPU_aos_f(
	struct part_aos_f *parts_aos, int *d_task_first_parts_pair,
	int *d_task_last_parts_pair, float d_a, float d_H, int bid, int tid, int count_tasks,
	int tasksperbundle, int nBlocks_per_task, int bundle_first_task, double *d_shift_x
	, double *d_shift_y, double *d_shift_z) {

  extern __shared__ float vars_pair_aosf[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
//  int first_part_in_task_blocks_ci, last_part_in_task_blocks_ci;
//  int first_part_in_task_blocks_cj, last_part_in_task_blocks_cj;
  const int task_id_tmp = 2 * task_id;
  const int ci_start = d_task_first_parts_pair[task_id_tmp];
  const int ci_end = d_task_last_parts_pair[task_id_tmp];
  const int cj_start = d_task_first_parts_pair[task_id_tmp + 1];
  const int cj_end = d_task_last_parts_pair[task_id_tmp + 1];

  /*Now do cj
   * Don't ever put me in an if statement. I've got __syncthreads inside*/
  const int pjd = threadid + cj_start;
  const int flip_j = -1;
  DOPAIR2NONSYMGPUAOSF(
  		parts_aos, pjd, cj_start, cj_end,
		ci_start, ci_end, d_a, d_H,
		vars_pair_aosf, d_shift_x, d_shift_y, d_shift_z, task_id_tmp + 1, flip_j);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
__global__ void runner_do_pair_force_GPU_aos_f4(
	struct part_aos_f4_f_send *parts_send, struct part_aos_f4_f_recv *parts_recv,
	float d_a, float d_H, int bundle_first_part, int bundle_n_parts) {

//  extern __shared__ float4 vars_pair_i_f4[];
//  __shared__ int first_part_tid_0, last_part_tid_0;
  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int pid = bundle_first_part + threadid;
//  const int task_id = bundle_first_part + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
  if(pid < bundle_first_part + bundle_n_parts){
    const struct part_aos_f4_f_send pi = parts_send[pid];
    const int cj_start = pi.cjs_cje.x;
    const int cj_end = pi.cjs_cje.y;
  /* Start calculations for particles in cell i */
    DOPAIR2NAIVEGPUAOSF4F(pi, parts_send, parts_recv, pid, cj_start, cj_end, d_a, d_H);
  }

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopair1_branch_density_gpu(struct part_soa parts_soa, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, int time_bin_inhibited, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max(max_parts_j, max_parts_i);
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
//  fprintf(stderr,"max_parts %i, max_partsi %i, max_partsj %i\n, "
//		  "numBlocks_x %i, numBlocks_y %i, BLOCK_SIZE %i\n", max_parts,
//		  max_parts_i, max_parts_j, numBlocks_x, numBlocks_y, BLOCK_SIZE);

  /*Do ci & cj*/
//  fprintf(stderr, "BLOCK_SIZE %i max parts %i num idle threads %i\n", BLOCK_SIZE, max_parts, numBlocks_x * BLOCK_SIZE - max_parts);

//  runner_do_pair_sym_density_GPU<<<gridShape, BLOCK_SIZE,
//		  13 * BLOCK_SIZE * sizeof(float) +
//		  3 * BLOCK_SIZE * sizeof(double) +
//              BLOCK_SIZE * sizeof(timebin_t),
//          stream>>>(
//      parts_soa, d_task_first_parts_pair, d_task_last_parts_pair,
//      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//      nBlocks_per_task, bundle_first_task, time_bin_inhibited, d_shift_x, d_shift_y, d_shift_z);

  runner_do_pair_nonsym_density_GPU<<<gridShape, BLOCK_SIZE,
		  5 * BLOCK_SIZE * sizeof(float) +
		  3 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(timebin_t),
          stream>>>(
      parts_soa, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, time_bin_inhibited, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopairci_branch_density_gpu(struct part_soa parts_soa, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_i;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_ci_density_GPU<<<gridShape, BLOCK_SIZE,
		  5 * BLOCK_SIZE * sizeof(float) +
		  3 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(timebin_t),
          stream>>>(
      parts_soa, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopaircj_branch_density_gpu(struct part_soa parts_soa, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_j;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_cj_density_GPU<<<gridShape, BLOCK_SIZE,
		  5 * BLOCK_SIZE * sizeof(float) +
		  3 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(timebin_t),
          stream>>>(
      parts_soa, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopairci_branch_density_gpu_aos(struct part_aos *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_i;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_ci_density_GPU_aos<<<gridShape, BLOCK_SIZE,
		  5 * BLOCK_SIZE * sizeof(float) +
		  3 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(int),
          stream>>>(
      parts_aos, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopaircj_branch_density_gpu_aos(struct part_aos *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_j;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_cj_density_GPU_aos<<<gridShape, BLOCK_SIZE,
		  5 * BLOCK_SIZE * sizeof(float) +
		  3 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(int),
          stream>>>(
      parts_aos, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif


#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopairci_branch_density_gpu_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_task, int4 *fparti_fpartj_lparti_lpartj_dens){

	  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
	  int nBlocks_per_task = numBlocks_x;


	  runner_do_pair_ci_density_GPU_aos_f4<<<gridShape, BLOCK_SIZE, 2 * BLOCK_SIZE * sizeof(float4), stream>>>(
			  parts_send, parts_recv, fparti_fpartj_lparti_lpartj_dens, d_a, d_H, bundle_first_task);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopaircj_branch_density_gpu_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
	      float d_a, float d_H, cudaStream_t stream,
		  int numBlocks_x, int numBlocks_y, int bundle_first_task, int4 *fparti_fpartj_lparti_lpartj_dens) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;

  runner_do_pair_cj_density_GPU_aos_f4<<<gridShape, BLOCK_SIZE, 2 * BLOCK_SIZE * sizeof(float4), stream>>>(
		  parts_send, parts_recv, fparti_fpartj_lparti_lpartj_dens, d_a, d_H, bundle_first_task);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopair_branch_density_gpu_aos_f4(struct part_aos_f4_send *parts_send, struct part_aos_f4_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_part, int bundle_n_parts){

	  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
	  int nBlocks_per_task = numBlocks_x;

//	  fprintf(stderr, "nblocks %i\n", numBlocks_x);
	  runner_do_pair_density_GPU_aos_f4<<<numBlocks_x, BLOCK_SIZE, 0, stream>>>(
			  parts_send, parts_recv, d_a, d_H, bundle_first_part, bundle_n_parts);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopairci_branch_density_gpu_aos_g(struct part_aos_g *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_i;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_ci_density_GPU_aos_g<<<gridShape, BLOCK_SIZE,
		  12 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(int),
          stream>>>(
      parts_aos, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopaircj_branch_density_gpu_aos_g(struct part_aos_g *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_j;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_cj_density_GPU_aos_g<<<gridShape, BLOCK_SIZE,
		  12 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(int),
          stream>>>(
      parts_aos, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopair_branch_gradient_gpu_aos_f4(struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_part, int bundle_n_parts){

	  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
	  int nBlocks_per_task = numBlocks_x;

//	  fprintf(stderr, "nblocks %i\n", numBlocks_x);
	  runner_do_pair_gradient_GPU_aos_f4<<<numBlocks_x, BLOCK_SIZE, 0, stream>>>(
			  parts_send, parts_recv, d_a, d_H, bundle_first_part, bundle_n_parts);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopairci_branch_density_gpu_aos_f(struct part_aos_f *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_i;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_ci_density_GPU_aos_f<<<gridShape, BLOCK_SIZE,
		  17 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(int),
          stream>>>(
      parts_aos, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopaircj_branch_density_gpu_aos_f(struct part_aos_f *parts_aos, int *d_task_first_parts_pair,
	          int *d_task_last_parts_pair, float d_a, float d_H, const char *loop_type, cudaStream_t stream,
			  int bid, int block_size, int count_tasks, int tasksperbundle,int max_parts_i, int max_parts_j,
			  int numBlocks_y, int tid, int offset, int bundle_first_task, double * d_shift_x
			  , double * d_shift_y, double * d_shift_z) {


  int max_parts = max_parts_j;
  int numBlocks_x = (max_parts + BLOCK_SIZE - 1) / BLOCK_SIZE;
  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;


  runner_do_pair_cj_density_GPU_aos_f<<<gridShape, BLOCK_SIZE,
		  17 * BLOCK_SIZE * sizeof(float) +
              BLOCK_SIZE * sizeof(int),
          stream>>>(
      parts_aos, d_task_first_parts_pair, d_task_last_parts_pair,
      d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, d_shift_x, d_shift_y, d_shift_z);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void runner_dopair_branch_force_gpu_aos_f4(struct part_aos_f4_f_send *parts_send, struct part_aos_f4_f_recv *parts_recv,
		      float d_a, float d_H, cudaStream_t stream,
			  int numBlocks_x, int numBlocks_y, int bundle_first_part, int bundle_n_parts){

	  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
	  int nBlocks_per_task = numBlocks_x;

//	  fprintf(stderr, "nblocks %i\n", numBlocks_x);
	  runner_do_pair_force_GPU_aos_f4<<<numBlocks_x, BLOCK_SIZE, 0, stream>>>(
			  parts_send, parts_recv, d_a, d_H, bundle_first_part, bundle_n_parts);

}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif

__global__ void runner_do_self_density_GPU_naive(
	    struct part_soa parts_soa, int *d_task_first_part, int *d_task_last_part, float d_a, float d_H,
	    int bid, int tid, int count_tasks, int tasksperbundle, int nBlocks_per_task,
	    int bundle_first_task, int max_parts, int time_bin_inhibited) {

  const int threadid = blockDim.x * blockIdx.x + threadIdx.x;
  const int task_id = bundle_first_task + blockIdx.y;

  //	printf("task_id is %i, count_tasks is %i\n", task_id, count_tasks);
  __shared__ int first_part_in_task_blocks, last_part_in_task_blocks;
  first_part_in_task_blocks = d_task_first_part[task_id];
  last_part_in_task_blocks = d_task_last_part[task_id];

  const int pid = threadid + first_part_in_task_blocks;

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
  int Found_neighbours = 0;

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

    int n_neighbours = 0;

    /*Naive loop over neighbours*/
    for (int b = first_part_in_task_blocks; b < last_part_in_task_blocks;
    	       b += BLOCK_SIZE) {
      for (int j_block = 0; j_block < BLOCK_SIZE; j_block++) {
    	int j = j_block + b;
    	if (j < last_part_in_task_blocks) {
    	  const float x_p_tmp = parts_soa.x_p[j];
    	  const float y_p_tmp = parts_soa.y_p[j];
    	  const float z_p_tmp = parts_soa.z_p[j];
    	  const float h_tmp = parts_soa.h[j];
    	  const float mass_tmp = parts_soa.mass[j];
    	  const float ux_tmp = parts_soa.ux[j];
    	  const float uy_tmp = parts_soa.uy[j];
    	  const float uz_tmp = parts_soa.uz[j];
    	  const timebin_t timebin = parts_soa.time_bin[j];

		  /* Compute the pairwise distance. */
		  const float pjx = x_p_tmp - cellx;
		  const float pjy = y_p_tmp - celly;
		  const float pjz = z_p_tmp - cellz;
		  const float xij = pix - pjx, yij = piy - pjy, zij = piz - pjz;
		  const float r2 = xij * xij + yij * yij + zij * zij;
		  const float hj = h_tmp, hjg2 = hj * hj * kernel_gamma2;
		  if (r2 < hig2 && r2 > (0.01f/128.f)*(0.01f/128.f)) {
			Found_neighbours=1;
			const float r = sqrt(r2);
			/* Recover some data */
			const float mj = mass_tmp;
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
			float dvx = uxi - ux_tmp, dvy = uyi - uy_tmp,
				dvz = uzi - uz_tmp;
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
    }
//    float wi, wi_dx;
//    d_kernel_deval(0.f, &wi, &wi_dx);
    if(Found_neighbours == 0) printf("Not sure what's going on but no neighbours found in GPU loop\n");
    parts_soa.rho[pid] = rhoi, parts_soa.rho_dh[pid] = rho_dhi;
    parts_soa.wcount[pid] = wcounti, parts_soa.wcount_dh[pid] = wcount_dhi;
    parts_soa.div_v[pid] = div_vi;
    parts_soa.rot_ux[pid] = rot_uxi, parts_soa.rot_uy[pid] = rot_uyi,
    parts_soa.rot_uz[pid] = rot_uzi;
  }
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_tester_kernel(struct part_soa parts_soa, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream, int bid,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y, int tid,
                           int offset, int bundle_first_task, int max_parts,
                           int time_bin_inhibited) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  tester<<<gridShape, BLOCK_SIZE,
                               8 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(timebin_t),
                               stream>>>(
      parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_density_kernel(struct part_soa parts_soa, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  runner_do_self_density_GPU<<<gridShape, BLOCK_SIZE,
                               8 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(timebin_t),
                               stream>>>(
      parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, max_parts);
//  runner_do_self_density_GPU_naive<<<gridShape, BLOCK_SIZE, 0, stream>>>(
//        parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//        nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_gradient_aos(struct part_aos_g *parts_aos, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts,
                           double * d_cell_x,
						   double * d_cell_y, double * d_cell_z) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  DOSELF_GPU_AOS_G<<<gridShape, BLOCK_SIZE,
                               12 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(int),
                               stream>>>(
      parts_aos, d_task_first_part, d_task_last_part, d_a, d_H, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, max_parts, d_cell_x,
	  d_cell_y, d_cell_z);
//  runner_do_self_density_GPU_naive<<<gridShape, BLOCK_SIZE, 0, stream>>>(
//        parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//        nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_gradient_aos_f4(struct part_aos_f4_g_send *parts_send, struct part_aos_f4_g_recv *parts_recv, float d_a, float d_H,
                           cudaStream_t stream, int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int2 * d_task_first_part_f4) {

	dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
	int nBlocks_per_task = numBlocks_x;
	DOSELF_GPU_AOS_F4_G<<<gridShape, BLOCK_SIZE,
	                        3 * BLOCK_SIZE * sizeof(float4), stream>>>(
	      parts_send, parts_recv, d_a, d_H, bundle_first_task, d_task_first_part_f4);
//  runner_do_self_density_GPU_naive<<<gridShape, BLOCK_SIZE, 0, stream>>>(
//        parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//        nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_force_aos(struct part_aos_f *parts_aos, int *d_task_first_part,
                           int *d_task_last_part, float d_a, float d_H,
                           const char *loop_type, cudaStream_t stream,
                           int block_size, int count_tasks, int tasksperbundle,
                           int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int max_parts,
                           double * d_cell_x,
						   double * d_cell_y, double * d_cell_z) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  DOSELF_GPU_AOS_F<<<gridShape, BLOCK_SIZE,
                               16 * BLOCK_SIZE * sizeof(float) +
                                   BLOCK_SIZE * sizeof(int),
                               stream>>>(
      parts_aos, d_task_first_part, d_task_last_part, d_a, d_H, count_tasks, tasksperbundle,
      nBlocks_per_task, bundle_first_task, max_parts, d_cell_x,
	  d_cell_y, d_cell_z);
//  runner_do_self_density_GPU_naive<<<gridShape, BLOCK_SIZE, 0, stream>>>(
//        parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//        nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif

#ifdef WITH_CUDA
extern "C" {
#endif
void launch_force_aos_f4(struct part_aos_f4_f_send *d_parts_send, struct part_aos_f4_f_recv *d_parts_recv, float d_a, float d_H,
                           cudaStream_t stream, int numBlocks_x, int numBlocks_y,
                           int bundle_first_task, int2 * d_task_first_part_f4) {

  dim3 gridShape = dim3(numBlocks_x, numBlocks_y);
  int nBlocks_per_task = numBlocks_x;
  DOSELF_GPU_AOS_F4_F<<<gridShape, BLOCK_SIZE,
                                   4 * BLOCK_SIZE * sizeof(float4) +
                                   BLOCK_SIZE * sizeof(float3),
                                   stream>>>(
      d_parts_send, d_parts_recv, d_a, d_H, bundle_first_task, d_task_first_part_f4);
//  runner_do_self_density_GPU_naive<<<gridShape, BLOCK_SIZE, 0, stream>>>(
//        parts_soa, d_task_first_part, d_task_last_part, d_a, d_H, bid, tid, count_tasks, tasksperbundle,
//        nBlocks_per_task, bundle_first_task, max_parts, time_bin_inhibited);
}
#ifdef WITH_CUDA
}
#endif
