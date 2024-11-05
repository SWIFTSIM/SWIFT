#include "cuda/part_gpu.h"

#include <cuda.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <stdio.h>

#ifdef WITH_CUDA
extern "C" {
#endif

#include "arrays_malloc.h"

void allocate_host(struct part_soa *parts_soa, int count_max_parts_tmp) {
  ///////////Malloc Host arrays
  cudaMallocHost((void **)&parts_soa->tid_p, count_max_parts_tmp * sizeof(int));
  cudaMallocHost((void **)&parts_soa->id,
                 count_max_parts_tmp * sizeof(long long));
  cudaMallocHost((void **)&parts_soa->mass,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->h, count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->u, count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->u_dt,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->rho, count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->SPH_sum,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->x_p,
                 count_max_parts_tmp * sizeof(double));
  cudaMallocHost((void **)&parts_soa->y_p,
                 count_max_parts_tmp * sizeof(double));
  cudaMallocHost((void **)&parts_soa->z_p,
                 count_max_parts_tmp * sizeof(double));
  cudaMallocHost((void **)&parts_soa->ux, count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->uy, count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->uz, count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->a_hydrox,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->a_hydroy,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->a_hydroz,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->locx,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->locy,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->locz,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->widthx,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->widthy,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->widthz,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->h_max,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->count_p,
                 count_max_parts_tmp * sizeof(int));
  cudaMallocHost((void **)&parts_soa->wcount,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->wcount_dh,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->rho_dh,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->rot_ux,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->rot_uy,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->rot_uz,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->div_v,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->div_v_previous_step,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->alpha_visc,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->v_sig,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->laplace_u,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->alpha_diff,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->f, count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->soundspeed,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->h_dt,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->balsara,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->pressure,
                 count_max_parts_tmp * sizeof(float));
  cudaMallocHost((void **)&parts_soa->alpha_visc_max_ngb,
                 count_max_parts_tmp * sizeof(float));
  /* timestep stuff */
  cudaMallocHost((void **)&parts_soa->time_bin,
                 count_max_parts_tmp * sizeof(timebin_t));
  cudaMallocHost((void **)&parts_soa->wakeup,
                 count_max_parts_tmp * sizeof(timebin_t));
  cudaMallocHost((void **)&parts_soa->min_ngb_time_bin,
                 count_max_parts_tmp * sizeof(timebin_t));
  cudaMallocHost((void **)&parts_soa->to_be_synchronized,
                 count_max_parts_tmp * sizeof(char));
}

void allocate_device(struct part_soa d_parts_soa, int count_max_parts_tmp) {
  ////////now malloc variables for particle data on the GPU. Sheesh
  fprintf(stderr, "before malloc\n");
  cudaMalloc((void **)&(d_parts_soa.tid_p), sizeof(int) * count_max_parts_tmp);
  fprintf(stderr, "after malloc\n");
  cudaMalloc((void **)&(d_parts_soa.id),
             sizeof(long long) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.x_p), sizeof(double) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.y_p), sizeof(double) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.z_p), sizeof(double) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.ux), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.uy), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.uz), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.a_hydrox),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.a_hydroy),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.a_hydroz),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.mass), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.h), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.u), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.u_dt), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.rho), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.SPH_sum),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.locx), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.locy), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.locz), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.widthx),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.widthy),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.widthz),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.h_max),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.count_p),
             sizeof(int) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.wcount),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.wcount_dh),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.rho_dh),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.rot_ux),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.rot_uy),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.rot_uz),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.div_v),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.div_v_previous_step),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.alpha_visc),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.v_sig),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.laplace_u),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.alpha_diff),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.f), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.soundspeed),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.h_dt), sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.balsara),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.pressure),
             sizeof(float) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.alpha_visc_max_ngb),
             sizeof(float) * count_max_parts_tmp);
  /* timestep stuff */
  cudaMalloc((void **)&(d_parts_soa.time_bin),
             sizeof(timebin_t) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.wakeup),
             sizeof(timebin_t) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.min_ngb_time_bin),
             sizeof(timebin_t) * count_max_parts_tmp);
  cudaMalloc((void **)&(d_parts_soa.to_be_synchronized),
             sizeof(char) * count_max_parts_tmp);
}

cudaError_t cudaAllocInt(int **d_var, int elements) {
  return cudaMalloc((void **)d_var, elements * sizeof(int));
}
cudaError_t cudaAllocFloat(float **d_var, int elements) {
  return cudaMalloc((void **)d_var, elements * sizeof(float));
}
cudaError_t cudaAllocDouble(double **d_var, int elements) {
  return cudaMalloc((void **)d_var, elements * sizeof(double));
}
cudaError_t cudaAllocLonglong(long long **d_var, int elements) {
  return cudaMalloc((void **)d_var, elements * sizeof(long long));
}
cudaError_t cudaAllocChar(char **d_var, int elements) {
  return cudaMalloc((void **)d_var, elements * sizeof(char));
}
cudaError_t cudaAllocTimebin(timebin_t **d_var, int elements) {
  return cudaMalloc((void **)d_var, elements * sizeof(timebin_t));
}

void allocate_device_dirty(
    int **d_tid_p, long long **d_id, double **d_x_p, double **d_y_p,
    double **d_z_p, float **d_ux, float **d_uy, float **d_uz,
    float **d_a_hydrox, float **d_a_hydroy, float **d_a_hydroz, float **d_mass,
    float **d_h, float **d_u, float **d_u_dt, float **d_rho, float **d_locx,
    float **d_locy, float **d_locz, float **d_widthx, float **d_widthy,
    float **d_widthz, float **d_h_max, int **d_count_p, float **d_wcount,
    float **d_wcount_dh, float **d_rho_dh, float **d_rot_ux, float **d_rot_uy,
    float **d_rot_uz, float **d_div_v, float **d_div_v_previous_step,
    float **d_alpha_visc, float **d_v_sig, float **d_laplace_u,
    float **d_alpha_diff, float **d_f, float **d_soundspeed, float **d_h_dt,
    float **d_balsara, float **d_pressure, float **d_alpha_visc_max_ngb,
    timebin_t **d_time_bin, timebin_t **d_wakeup,
    timebin_t **d_min_ngb_time_bin, char **d_to_be_synchronized,
    int count_max_parts_tmp) {
  ////////Malloc variables for particle data on the GPU. Sheesh, that's a lot

  size_t free_byte;
  size_t total_byte;

  cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);
  double free = (double)free_byte;
  double available = (double)total_byte;
  double used = (available - free);
  //          message("free %lf used %lf", free/10.E8, used/10.E8);

  cudaError_t cu_error = cudaAllocInt(d_tid_p, count_max_parts_tmp);
  cu_error = cudaAllocLonglong(d_id, count_max_parts_tmp);
  cu_error = cudaAllocDouble(d_x_p, count_max_parts_tmp);
  cu_error = cudaAllocDouble(d_y_p, count_max_parts_tmp);
  cu_error = cudaAllocDouble(d_z_p, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_ux, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_uy, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_uz, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_a_hydrox, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_a_hydroy, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_a_hydroz, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_mass, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_h, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_u, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_u_dt, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_rho, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_locx, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_locy, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_locz, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_widthx, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_widthy, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_widthz, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_h_max, count_max_parts_tmp);
  cu_error = cudaAllocInt(d_count_p, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_wcount, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_wcount_dh, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_rho_dh, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_rot_ux, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_rot_uy, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_rot_uz, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_div_v, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_div_v_previous_step, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_alpha_visc, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_v_sig, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_laplace_u, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_alpha_diff, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_f, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_soundspeed, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_h_dt, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_balsara, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_pressure, count_max_parts_tmp);
  cu_error = cudaAllocFloat(d_alpha_visc_max_ngb, count_max_parts_tmp);
  /* timestep stuff */
  cu_error = cudaAllocTimebin(d_time_bin, count_max_parts_tmp);
  cu_error = cudaAllocTimebin(d_wakeup, count_max_parts_tmp);
  cu_error = cudaAllocTimebin(d_min_ngb_time_bin, count_max_parts_tmp);
  cu_error = cudaAllocChar(d_to_be_synchronized, count_max_parts_tmp);
//		  cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
//		  double free_end = (double)free_byte;
//		  available = (double)total_byte;
//		  double used_end = (available - free_end);
//          message("cuda malloc self free %lf GB used %lf GB used to allocate
//          self"
//        		  " data %lf MB", free_end/10.E8, used_end/10.E8,
//        (used_end - used)/10.E5);
//          message("at end of malloc dirty: %s",
//		  	       cudaGetErrorString(cu_error));
#ifdef CUDA_DEBUG
  if (cu_error != cudaSuccess) {
    fprintf(stderr, "CUDA error at end of malloc dirty: %s\n",
            cudaGetErrorString(cu_error));
    exit(0);
  }
#endif
}

void allocate_device_test(int **tid_test, int count_max_parts_tmp) {
  ////////now malloc variables for particle data on the GPU. Sheesh

  cudaMalloc((void **)tid_test, sizeof(int) * count_max_parts_tmp);

  cudaError_t cu_error = cudaPeekAtLastError();  // Get error code
  fprintf(stderr, "malloc tid: %s\n", cudaGetErrorString(cu_error));

  if (cu_error != cudaSuccess) {
    fprintf(stderr, "CUDA error with malloc tid: %s\n",
            cudaGetErrorString(cu_error));
    exit(0);
  }
}
/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_malloc(struct part_soa *parts_soa, int alloc_type,
                 int count_max_parts_tmp) {
  allocate_host(parts_soa, count_max_parts_tmp);
}
/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_malloc(struct part_soa d_parts_soa, int alloc_type,
                   int count_max_parts_tmp) {
  allocate_device(d_parts_soa, count_max_parts_tmp);
}

void device_malloc_dirty(
    int **d_tid_p, long long **d_id, double **d_x_p, double **d_y_p,
    double **d_z_p, float **d_ux, float **d_uy, float **d_uz,
    float **d_a_hydrox, float **d_a_hydroy, float **d_a_hydroz, float **d_mass,
    float **d_h, float **d_u, float **d_u_dt, float **d_rho, float **d_locx,
    float **d_locy, float **d_locz, float **d_widthx, float **d_widthy,
    float **d_widthz, float **d_h_max, int **d_count_p, float **d_wcount,
    float **d_wcount_dh, float **d_rho_dh, float **d_rot_ux, float **d_rot_uy,
    float **d_rot_uz, float **d_div_v, float **d_div_v_previous_step,
    float **d_alpha_visc, float **d_v_sig, float **d_laplace_u,
    float **d_alpha_diff, float **d_f, float **d_soundspeed, float **d_h_dt,
    float **d_balsara, float **d_pressure, float **d_alpha_visc_max_ngb,
    timebin_t **d_time_bin, timebin_t **d_wakeup,
    timebin_t **d_min_ngb_time_bin, char **d_to_be_synchronized,
    int count_max_parts_tmp) {

  allocate_device_dirty(
      d_tid_p, d_id, d_x_p, d_y_p, d_z_p, d_ux, d_uy, d_uz, d_a_hydrox,
      d_a_hydroy, d_a_hydroz, d_mass, d_h, d_u, d_u_dt, d_rho, d_locx, d_locy,
      d_locz, d_widthx, d_widthy, d_widthz, d_h_max, d_count_p, d_wcount,
      d_wcount_dh, d_rho_dh, d_rot_ux, d_rot_uy, d_rot_uz, d_div_v,
      d_div_v_previous_step, d_alpha_visc, d_v_sig, d_laplace_u, d_alpha_diff,
      d_f, d_soundspeed, d_h_dt, d_balsara, d_pressure, d_alpha_visc_max_ngb,
      d_time_bin, d_wakeup, d_min_ngb_time_bin, d_to_be_synchronized,
      count_max_parts_tmp);
}

void device_malloc_test(int **tid_test, int count_max_parts_tmp) {

  allocate_device_test(tid_test, count_max_parts_tmp);
}

#ifdef WITH_CUDA
}
#endif
