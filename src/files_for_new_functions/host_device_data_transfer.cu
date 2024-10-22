#include "cuda/part_gpu.h"
#include <cuda.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <stdio.h>

#ifdef WITH_CUDA
extern "C" {
#endif

void host2device_test(int *d_tid_p, int *tid_h, int count_max_parts_tmp){
//	int * tid_h;
//	cudaMallocHost((void **)&tid_h,
//			count_max_parts_tmp * sizeof(int));
	for (int i =0; i< count_max_parts_tmp; i++){
		tid_h[i] = 100;
//		fprintf(stderr,"tid_h %i\n", tid_h[i]);
	}

	cudaMemcpy(d_tid_p, tid_h, count_max_parts_tmp * sizeof(int), cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
//	cudaFree(tid_h);
}

void device2host_test(struct part_soa parts_soa, int *tid_h, int count_max_parts_tmp){
	int *tid_p = parts_soa.tid_p;
	cudaMemcpy(tid_h, tid_p, count_max_parts_tmp * sizeof(int), cudaMemcpyDeviceToHost);
	for (int i =0; i< count_max_parts_tmp; i++){
		fprintf(stderr,"tid is %i\n", tid_h[i]);
	}
}

void device2device_test(int *tid_p, struct part_soa parts_soa, int count_max_parts_tmp){
	cudaMemcpy(tid_p, parts_soa.tid_p, sizeof(int *), cudaMemcpyHostToDevice);
}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_device_test(int *d_tid_p, int *tid_h, int count_max_parts_tmp){

	host2device_test(d_tid_p, tid_h, count_max_parts_tmp);

}

void device_host_test(struct part_soa parts_soa, int *tid_h, int count_max_parts_tmp){

	device2host_test(parts_soa, tid_h, count_max_parts_tmp);

}

void device_device_test(int *tid_p, struct part_soa parts_soa, int count_max_parts_tmp){

	device2device_test(tid_p, parts_soa, count_max_parts_tmp);

}

void device2host_density(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int count_max_parts_tmp){
  cudaMemcpy(parts_soa_buffer.tid_p, tid_p, count_max_parts_tmp * sizeof(int), cudaMemcpyDeviceToHost);
}
void device_host_cpy(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int count_max_parts_tmp){

	device2host_density(parts_soa_buffer, tid_p, id, x_p, y_p, z_p,
			 ux, uy, uz, a_hydrox, a_hydroy, a_hydroz,
			 mass, h , u, u_dt, rho, locx, locy,
			 locz, widthx, widthy, widthz, h_max, count_p,
			 wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
			 div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
			 alpha_diff, f, soundspeed, h_dt, balsara, pressure,
			 alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
	         to_be_synchronized, count_max_parts_tmp);

}

void device2device_density(struct part_soa *parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int count_max_parts_tmp, cudaStream_t stream){

  cudaMemcpyAsync(&(parts_soa_buffer->tid_p), &tid_p,
				  sizeof(int *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->locx), &locx,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->locy), &locy,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->locz), &locz,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->h), &h,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->mass), &mass,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->x_p), &x_p,
				  sizeof(double *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->y_p), &y_p,
				  sizeof(double *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->z_p), &z_p,
				  sizeof(double *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->ux), &ux,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->uy), &uy,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->uz), &uz,
				  sizeof(float *), cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&(parts_soa_buffer->time_bin), &time_bin,
				  sizeof(timebin_t *), cudaMemcpyHostToDevice, stream);

}


void host2device_density(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int count_max_parts_tmp){
  cudaError_t cu_error;
  cudaMemcpy(&tid_p, &(parts_soa_buffer.tid_p),
		  count_max_parts_tmp * sizeof(int), cudaMemcpyHostToDevice);
}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_device_cpy(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int count_max_parts_tmp){

	host2device_density(parts_soa_buffer, tid_p, id, x_p, y_p, z_p,
			 ux, uy, uz, a_hydrox, a_hydroy, a_hydroz,
			 mass, h , u, u_dt, rho, locx, locy,
			 locz, widthx, widthy, widthz, h_max, count_p,
			 wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
			 div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
			 alpha_diff, f, soundspeed, h_dt, balsara, pressure,
			 alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
	         to_be_synchronized, count_max_parts_tmp);

}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_device_bind(struct part_soa *parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int count_max_parts_tmp, cudaStream_t stream){

	device2device_density(parts_soa_buffer, tid_p, id, x_p, y_p, z_p,
				 ux, uy, uz, a_hydrox, a_hydroy, a_hydroz,
				 mass, h , u, u_dt, rho, locx, locy,
				 locz, widthx, widthy, widthz, h_max, count_p,
				 wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
				 div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
				 alpha_diff, f, soundspeed, h_dt, balsara, pressure,
				 alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
		         to_be_synchronized, count_max_parts_tmp, stream);

}

void host2device_async_density(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp, int bundle_n_parts, cudaStream_t stream){
  cudaError_t cu_error;
  cudaMemcpyAsync(&tid_p[first_part_tmp], &(parts_soa_buffer.tid_p[first_part_tmp]),
				  bundle_n_parts * sizeof(int), cudaMemcpyHostToDevice,
				  stream);
  cudaMemcpyAsync(&locx[first_part_tmp], &(parts_soa_buffer.locx[first_part_tmp]),
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&locy[first_part_tmp], &(parts_soa_buffer.locy[first_part_tmp]),
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&locz[first_part_tmp], &parts_soa_buffer.locz[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&h[first_part_tmp], &parts_soa_buffer.h[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&mass[first_part_tmp], &parts_soa_buffer.mass[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&x_p[first_part_tmp], &parts_soa_buffer.x_p[first_part_tmp],
				  bundle_n_parts * sizeof(double),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&y_p[first_part_tmp], &parts_soa_buffer.y_p[first_part_tmp],
				  bundle_n_parts * sizeof(double),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&z_p[first_part_tmp], &parts_soa_buffer.z_p[first_part_tmp],
				  bundle_n_parts * sizeof(double),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&ux[first_part_tmp], &parts_soa_buffer.ux[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&uy[first_part_tmp], &parts_soa_buffer.uy[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&uz[first_part_tmp], &parts_soa_buffer.uz[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&time_bin[first_part_tmp],
				  &parts_soa_buffer.time_bin[first_part_tmp],
				  bundle_n_parts * sizeof(timebin_t),
				  cudaMemcpyHostToDevice, stream);
}

void host2device_async_density_pair(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp, int bundle_n_parts, cudaStream_t stream){

//  int bundle_n_parts = bundle_n_parts_i + bundle_n_parts_j;
  cudaError_t cu_error;
//  cudaMemcpyAsync(&tid_p[first_part_tmp], &(parts_soa_buffer.tid_p[first_part_tmp]),
//				  bundle_n_parts * sizeof(int), cudaMemcpyHostToDevice,
//				  stream);
//  cudaMemcpyAsync(&locx[first_part_tmp], &(parts_soa_buffer.locx[first_part_tmp]),
//				  bundle_n_parts * sizeof(float),
//				  cudaMemcpyHostToDevice, stream);
//  cudaMemcpyAsync(&locy[first_part_tmp], &(parts_soa_buffer.locy[first_part_tmp]),
//				  bundle_n_parts * sizeof(float),
//				  cudaMemcpyHostToDevice, stream);
//  cudaMemcpyAsync(&locz[first_part_tmp], &parts_soa_buffer.locz[first_part_tmp],
//				  bundle_n_parts * sizeof(float),
//				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&h[first_part_tmp], &parts_soa_buffer.h[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&mass[first_part_tmp], &parts_soa_buffer.mass[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&x_p[first_part_tmp], &parts_soa_buffer.x_p[first_part_tmp],
				  bundle_n_parts * sizeof(double),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&y_p[first_part_tmp], &parts_soa_buffer.y_p[first_part_tmp],
				  bundle_n_parts * sizeof(double),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&z_p[first_part_tmp], &parts_soa_buffer.z_p[first_part_tmp],
				  bundle_n_parts * sizeof(double),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&ux[first_part_tmp], &parts_soa_buffer.ux[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&uy[first_part_tmp], &parts_soa_buffer.uy[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&uz[first_part_tmp], &parts_soa_buffer.uz[first_part_tmp],
				  bundle_n_parts * sizeof(float),
				  cudaMemcpyHostToDevice, stream);
  cudaMemcpyAsync(&time_bin[first_part_tmp],
				  &parts_soa_buffer.time_bin[first_part_tmp],
				  bundle_n_parts * sizeof(timebin_t),
				  cudaMemcpyHostToDevice, stream);
}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_device_async_cpy(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp, int bundle_n_parts, cudaStream_t stream){

	host2device_async_density(parts_soa_buffer, tid_p, id, x_p, y_p, z_p,
			 ux, uy, uz, a_hydrox, a_hydroy, a_hydroz,
			 mass, h , u, u_dt, rho, locx, locy,
			 locz, widthx, widthy, widthz, h_max, count_p,
			 wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
			 div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
			 alpha_diff, f, soundspeed, h_dt, balsara, pressure,
			 alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
	         to_be_synchronized, first_part_tmp, bundle_n_parts, stream);

}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_device_async_cpy_pair(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp_i, int bundle_n_parts, cudaStream_t stream){

	host2device_async_density_pair(parts_soa_buffer, tid_p, id, x_p, y_p, z_p,
			 ux, uy, uz, a_hydrox, a_hydroy, a_hydroz,
			 mass, h , u, u_dt, rho, locx, locy,
			 locz, widthx, widthy, widthz, h_max, count_p,
			 wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
			 div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
			 alpha_diff, f, soundspeed, h_dt, balsara, pressure,
			 alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
	         to_be_synchronized, first_part_tmp_i, bundle_n_parts,
			 stream);

}

void device2host_async_density(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp, int bundle_n_parts, cudaStream_t stream){
  cudaError_t cu_error;

  cudaMemcpyAsync(&parts_soa_buffer.rho[first_part_tmp], &rho[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rho_dh[first_part_tmp],
				&rho_dh[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.wcount[first_part_tmp],
				&wcount[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.wcount_dh[first_part_tmp],
				&wcount_dh[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.div_v[first_part_tmp], &div_v[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rot_ux[first_part_tmp],
				&rot_ux[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rot_uy[first_part_tmp],
				&rot_uy[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rot_uz[first_part_tmp],
				&rot_uz[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
}

void device2host_async_density_pair(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp, int bundle_n_parts, cudaStream_t stream){
  cudaError_t cu_error;
//  fprintf(stderr, "parts i %i parts j %i\n", bundle_n_parts_i, bundle_n_parts_j);
//  int bundle_n_parts = bundle_n_parts_i + bundle_n_parts_j;

  cudaMemcpyAsync(&parts_soa_buffer.rho[first_part_tmp], &rho[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rho_dh[first_part_tmp],
				&rho_dh[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.wcount[first_part_tmp],
				&wcount[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.wcount_dh[first_part_tmp],
				&wcount_dh[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.div_v[first_part_tmp], &div_v[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rot_ux[first_part_tmp],
				&rot_ux[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rot_uy[first_part_tmp],
				&rot_uy[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
  cudaMemcpyAsync(&parts_soa_buffer.rot_uz[first_part_tmp],
				&rot_uz[first_part_tmp],
				bundle_n_parts * sizeof(float),
				cudaMemcpyDeviceToHost, stream);
}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_host_async_cpy(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp, int bundle_n_parts, cudaStream_t stream){

	device2host_async_density(parts_soa_buffer, tid_p, id, x_p, y_p, z_p,
			 ux, uy, uz, a_hydrox, a_hydroy, a_hydroz,
			 mass, h , u, u_dt, rho, locx, locy,
			 locz, widthx, widthy, widthz, h_max, count_p,
			 wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
			 div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
			 alpha_diff, f, soundspeed, h_dt, balsara, pressure,
			 alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
	         to_be_synchronized, first_part_tmp, bundle_n_parts, stream);

}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_host_async_cpy_pair(struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized, int first_part_tmp, int bundle_n_parts, cudaStream_t stream){

	device2host_async_density_pair(parts_soa_buffer, tid_p, id, x_p, y_p, z_p,
			 ux, uy, uz, a_hydrox, a_hydroy, a_hydroz,
			 mass, h , u, u_dt, rho, locx, locy,
			 locz, widthx, widthy, widthz, h_max, count_p,
			 wcount, wcount_dh, rho_dh, rot_ux, rot_uy, rot_uz,
			 div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
			 alpha_diff, f, soundspeed, h_dt, balsara, pressure,
			 alpha_visc_max_ngb, time_bin, wakeup, min_ngb_time_bin,
	         to_be_synchronized, first_part_tmp, bundle_n_parts, stream);

}

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_device_async_bind(struct part_soa *parts_soa, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
		float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
		float *mass, float *h ,float *u, float *u_dt, float *rho, float *locx, float *locy,
		float *locz, float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
		float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux, float *rot_uy, float *rot_uz,
		float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig, float *laplace_u,
		float *alpha_diff, float *f, float *soundspeed, float *h_dt, float *balsara,float *pressure,
		float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
        char *to_be_synchronized){

    parts_soa->tid_p = tid_p;
    parts_soa->locx = locx;
    parts_soa->locy = locy;
    parts_soa->locz = locz;
    parts_soa->h = h;
    parts_soa->mass = mass;
    parts_soa->x_p = x_p;
    parts_soa->y_p = y_p;
    parts_soa->z_p = z_p;
    parts_soa->rho = rho;
    parts_soa->rho_dh = rho_dh;
    parts_soa->wcount = wcount;
    parts_soa->wcount_dh = wcount_dh;
    parts_soa->ux = ux;
    parts_soa->uy = uy;
    parts_soa->uz = uz;
    parts_soa->div_v = div_v;
    parts_soa->rot_ux = rot_ux;
    parts_soa->rot_uy = rot_uy;
    parts_soa->rot_uz = rot_uz;
    parts_soa->time_bin = time_bin;

}

#ifdef WITH_CUDA
}
#endif
