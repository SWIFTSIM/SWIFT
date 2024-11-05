#include "cuda/part_gpu.h"

#include <cuda.h>
#include <cuda_profiler_api.h>
#include <cuda_runtime.h>
#include <stdio.h>

void host2device_test(int *d_tid_p, int *tid_h, int count_max_parts_tmp);

void device2host_test(struct part_soa parts_soa, int *tid_h,
                      int count_max_parts_tmp);

void device2device_test(int *tid_p, struct part_soa parts_soa,
                        int count_max_parts_tmp);

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_device_test(int *d_tid_p, int *tid_h, int count_max_parts_tmp);

void device_host_test(struct part_soa parts_soa, int *tid_h,
                      int count_max_parts_tmp);

void device_device_test(int *tid_p, struct part_soa parts_soa,
                        int count_max_parts_tmp);

void device2host_density(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized,
    int count_max_parts_tmp);

void device_host_cpy(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized,
    int count_max_parts_tmp);

void device2device_density(
    struct part_soa *parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized,
    int count_max_parts_tmp, cudaStream_t stream);

void host2device_density(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized,
    int count_max_parts_tmp);

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_device_cpy(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized,
    int count_max_parts_tmp);

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_device_bind(
    struct part_soa *parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized,
    int count_max_parts_tmp, cudaStream_t stream);

void host2device_async_density(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized, int first_part_tmp,
    int bundle_n_parts, cudaStream_t stream);

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void host_device_async_cpy(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized, int first_part_tmp,
    int bundle_n_parts, cudaStream_t stream);

void device2host_async_density(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized, int first_part_tmp,
    int bundle_n_parts, cudaStream_t stream);
/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_host_async_cpy(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized, int first_part_tmp,
    int bundle_n_parts, cudaStream_t stream);

/*Function to be overloaded using different part_soa structs
 * and allocate their internal arrays
 * alloc_type 0 for density, 1 for force, 2 for gradient*/
void device_device_async_bind(
    struct part_soa *parts_soa, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized);

void host_device_async_cpy_pair(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized, int first_part_tmp,
    int bundle_n_parts, cudaStream_t stream);

void device_host_async_cpy_pair(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized, int first_part_tmp,
    int bundle_n_parts, cudaStream_t stream);

void device2host_async_density_pair(
    struct part_soa parts_soa_buffer, int *tid_p, long long *id, double *x_p,
    double *y_p, double *z_p, float *ux, float *uy, float *uz, float *a_hydrox,
    float *a_hydroy, float *a_hydroz, float *mass, float *h, float *u,
    float *u_dt, float *rho, float *locx, float *locy, float *locz,
    float *widthx, float *widthy, float *widthz, float *h_max, int *count_p,
    float *wcount, float *wcount_dh, float *rho_dh, float *rot_ux,
    float *rot_uy, float *rot_uz, float *div_v, float *div_v_previous_step,
    float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
    float *f, float *soundspeed, float *h_dt, float *balsara, float *pressure,
    float *alpha_visc_max_ngb, timebin_t *time_bin, timebin_t *wakeup,
    timebin_t *min_ngb_time_bin, char *to_be_synchronized, int first_part_tmp,
    int bundle_n_parts_i, int bundle_n_parts_j, cudaStream_t stream);
