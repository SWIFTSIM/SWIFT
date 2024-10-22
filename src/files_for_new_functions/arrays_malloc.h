#include "cuda/part_gpu.h"
#include <cuda.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <error.h>

cudaError_t cudaAllocInt(int ** d_var, int elements);
cudaError_t cudaAllocFloat(float ** d_var, int elements);
cudaError_t cudaAllocDouble(double ** d_var, int elements);
cudaError_t cudaAllocLonglong(long long ** d_var, int elements);
cudaError_t cudaAllocChar(char ** d_var, int elements);
cudaError_t cudaAllocTimebin(timebin_t ** d_var, int elements);

void allocate_host(struct part_soa *parts_soa, int count_max_parts_tmp);

void allocate_device(struct part_soa d_parts_soa, int count_max_parts_tmp);

void allocate_device_dirty(int **d_tid_p, long long **d_id, double **d_x_p, double **d_y_p, double **d_z_p,
		float **d_ux, float **d_uy, float **d_uz, float **d_a_hydrox, float **d_a_hydroy, float **d_a_hydroz,
		float **d_mass, float **d_h ,float **d_u, float **d_u_dt, float **d_rho, float **d_locx, float **d_locy,
		float **d_locz, float **d_widthx, float **d_widthy, float **d_widthz, float **d_h_max, int **d_count_p,
		float **d_wcount, float **d_wcount_dh, float **d_rho_dh, float **d_rot_ux, float **d_rot_uy, float **d_rot_uz,
		float **d_div_v, float **d_div_v_previous_step, float **d_alpha_visc, float **d_v_sig, float **d_laplace_u,
		float **d_alpha_diff, float **d_f, float **d_soundspeed, float **d_h_dt, float **d_balsara,float **d_pressure,
		float **d_alpha_visc_max_ngb, timebin_t **d_time_bin, timebin_t **d_wakeup, timebin_t **d_min_ngb_time_bin,
        char **d_to_be_synchronized, int count_max_parts_tmp);

void allocate_device_test(int **tid_test, int count_max_parts_tmp);
/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_malloc(struct part_soa *parts_soa, int alloc_type, int count_max_parts_tmp);

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_malloc(struct part_soa d_parts_soa, int alloc_type, int count_max_parts_tmp);

void device_malloc_dirty(int **d_tid_p, long long **d_id, double **d_x_p, double **d_y_p, double **d_z_p,
		float **d_ux, float **d_uy, float **d_uz, float **d_a_hydrox, float **d_a_hydroy, float **d_a_hydroz,
		float **d_mass, float **d_h ,float **d_u, float **d_u_dt, float **d_rho, float **d_locx, float **d_locy,
		float **d_locz, float **d_widthx, float **d_widthy, float **d_widthz, float **d_h_max, int **d_count_p,
		float **d_wcount, float **d_wcount_dh, float **d_rho_dh, float **d_rot_ux, float **d_rot_uy, float **d_rot_uz,
		float **d_div_v, float **d_div_v_previous_step, float **d_alpha_visc, float **d_v_sig, float **d_laplace_u,
		float **d_alpha_diff, float **d_f, float **d_soundspeed, float **d_h_dt, float **d_balsara,float **d_pressure,
		float **d_alpha_visc_max_ngb, timebin_t **d_time_bin, timebin_t **d_wakeup, timebin_t **d_min_ngb_time_bin,
        char **d_to_be_synchronized, int count_max_parts_tmp);

void device_malloc_test(int **tid_test, int count_max_parts_tmp);

