#ifndef PART_GPU_H
#define PART_GPU_H
/* Config parameters. */
#include "../../config.h"
#include "../align.h"
typedef int8_t timebin_t;

#ifdef __WITH_CUDA
extern "C" {
#endif

#include "/usr/local/cuda-12.3/targets/x86_64-linux/include/vector_types.h"

typedef struct part_soa {
  /*Task ID*/
  int *tid_p;
  /*bundle ID*/
  int *bid_p;
  /*! Particle unique ID. */
  long long *id;
  /*! Pointer to corresponding gravity part. */
  //	struct gpu_gpart* gpart;
  /*! Particle position. */
  double *x_p;
  double *y_p;
  double *z_p;
  /*! Particle predicted velocity. */
  float *ux;
  float *uy;
  float *uz;
  /*! Particle acceleration. */
  float *a_hydrox;
  float *a_hydroy;
  float *a_hydroz;
  /*! Particle mass. */
  float *mass;
  /*! Particle smoothing length. */
  float *h;
  /*! Particle internal energy. */
  float *u;
  /*! Time derivative of the internal energy. */
  float *u_dt;
  /*! Particle density. */
  float *rho;
  /*! Kernel summation (For testing/debugging). */
  float *SPH_sum;

  /* Cell information */
  /*! The cell location on the grid (corner nearest to the origin). */
  float *locx;
  float *locy;
  float *locz;
  /*! The cell dimensions. */
  float *widthx;
  float *widthy;
  float *widthz;
  float *h_max;
  int *count_p;
  int *count_test;
  /* Density information */

  /*! Neighbour number count. */
  float *wcount;

  /*! Derivative of the neighbour number with respect to h. */
  float *wcount_dh;

  /*! Derivative of density with respect to h */
  float *rho_dh;

  /*! Particle velocity curl. */
  float *rot_ux;
  float *rot_uy;
  float *rot_uz;

  /* viscosity information */

  /*! Particle velocity divergence */
  float *div_v;

  /*! Particle velocity divergence from previous step */
  float *div_v_previous_step;

  /*! Artificial viscosity parameter */
  float *alpha_visc;

  /*! Signal velocity */
  float *v_sig;

  /* thermal diffusion information  */

  /*! del^2 u, a smoothed quantity */
  float *laplace_u;

  /*! Thermal diffusion coefficient */
  float *alpha_diff;

  /* force information  */

  /*! "Grad h" term -- only partial in P-U */
  float *f;

  /*! Particle soundspeed. */
  float *soundspeed;

  /*! Time derivative of smoothing length  */
  float *h_dt;

  /*! Balsara switch */
  float *balsara;

  /*! Particle pressure. */
  float *pressure;
  /*! Maximal alpha (viscosity) over neighbours */
  float *alpha_visc_max_ngb;

  /* timestep stuff */

  /*! Time-step length */
  timebin_t *time_bin;

  /*all part of struct timestep_limiter_data, we had to destruct it
   as GPUs don't like pointer chasing especially when memcpying*/
  /* Need waking-up ? */
  timebin_t *wakeup;

  /*! Minimal time-bin across all neighbours */
  timebin_t *min_ngb_time_bin;

  /* Do we want this particle to be synched back on the time-line? */
  char *to_be_synchronized;
} part_soa;
/*Container for particle data requierd for density calcs*/
typedef struct part_aos {

  /*! Particle position. */
  double x_p;
  double y_p;
  double z_p;

  /*! Particle position. */
  double locx;
  double locy;
  double locz;

  /*! Particle predicted velocity. */
  float ux;
  float uy;
  float uz;
  /*! Particle mass. */
  float mass;
  /*! Particle smoothing length. */
  float h;
  /*! Particle density. */
  float rho;

  /* Density information */
  /*! Neighbour number count. */
  float wcount;
  /*! Derivative of the neighbour number with respect to h. */
  float wcount_dh;
  /*! Derivative of density with respect to h */
  float rho_dh;
  /*! Particle velocity curl. */
  float rot_ux;
  float rot_uy;
  float rot_uz;

  /* viscosity information */
  /*! Particle velocity divergence */
  float div_v;

  /* timestep stuff */
  /*! Time-step length */
  int time_bin;
} part_aos;

/*Container for particle data requierd for density calcs*/
typedef struct part_aos_f4_send {
  /*! Particle position and h -> x, y, z, h */
  float4 x_p_h;

  /*! Particle predicted velocity and mass -> ux, uy, uz, m */
  float4 ux_m;
  /*Markers for where neighbour cell j starts and stops in array indices for pair tasks*/
  int2 cjs_cje;
}part_aos_f4_send __attribute__((aligned(SWIFT_STRUCT_ALIGNMENT)));

typedef struct part_aos_f4_recv{
  /* Density information; rho */
  /*! Derivative of density with respect to h; rho_dh,
  * Neighbour number count; w_count
  * * Derivative of the neighbour number with respect to h; w_count_dh */
  float4 rho_dh_wcount;
  /*! Particle velocity curl; rot_ux and
  * velocity divergence; div_v */
  float4 rot_ux_div_v;
} part_aos_f4_recv;

/*Container for particle data required for density calcs*/
typedef struct part_aos_f4 {
  /*! Particle position and h -> x, y, z, h */
  float4 x_p_h;

  /*! Particle predicted velocity and mass -> ux, uy, uz, m */
  float4 ux_m;
  /* Density information; rho */
  /*! Derivative of density with respect to h; rho_dh,
  * Neighbour number count; w_count
  * * Derivative of the neighbour number with respect to h; w_count_dh */
  float4 rho_dh_wcount;

  /*! Particle velocity curl; rot_ux and
  * velocity divergence; div_v */
  float4 rot_ux_div_v;

} part_aos_f4;

/*Container for particle data required for force calcs*/
typedef struct part_aos_f {

  /*! Particle position. */
  double x_p;
  double y_p;
  double z_p;

  /*! Particle predicted velocity. */
  float ux;
  float uy;
  float uz;
  /*! Particle mass. */
  float mass;
  /*! Particle smoothing length. */
  float h;
  /*! Particle density. */
  float rho;
  /*! Particle pressure. */
  float pressure;

  /* Density information */
  /*! Speed of sound. */
  float soundspeed;
  /*! Variable smoothing length term */
  float f;
  /*! Derivative of density with respect to h */
  float balsara;
  /*! Particle velocity curl. */
  float alpha_visc;
  float a_hydrox;
  float a_hydroy;
  float a_hydroz;
  float alpha_diff;

  /* viscosity information */
  /*! Internal energy */
  float u;
  float u_dt;
  /*! h time derivative */
  float h_dt;
  float v_sig;

  /* timestep stuff */
  /*! Time-step length */
  int time_bin;
  int min_ngb_time_bin;
} part_aos_f;

/*Container for particle data requierd for force calcs*/
typedef struct part_aos_f4_f {

  /*Data required for the calculation:
  Values read to local GPU memory*/
  /*! Particle position smoothing length */
  float4 x_h;
  /*! Particle predicted velocity and mass */
  float4 ux_m;
  /*! Variable smoothing length term f, balsara, timebin
   * and initial value of min neighbour timebin */
  float4 f_bals_timebin_mintimebin_ngb;
  /*! Particle density, pressure, speed of sound & v_sig to read*/
  float4 rho_p_c_vsigi;
  /*! Particle Internal energy u, alpha constants for visc and diff */
  float3 u_alphavisc_alphadiff;

  /*Result: Values output to global GPU memory*/
  /* change of u and h with dt, v_sig and returned value of
   * minimum neighbour timebin */
  float4 udt_hdt_vsig_mintimebin_ngb;
  /*Particle acceleration vector*/
  float3 a_hydro;

} part_aos_f4_f;

/*Container for particle data requierd for force calcs*/
typedef struct part_aos_f4_f_send {

  /*Data required for the calculation:
  Values read to local GPU memory*/
  /*! Particle position smoothing length */
  float4 x_h;
  /*! Particle predicted velocity and mass */
  float4 ux_m;
  /*! Variable smoothing length term f, balsara, timebin
   * and initial value of min neighbour timebin */
  float4 f_bals_timebin_mintimebin_ngb;
  /*! Particle density, pressure, speed of sound & v_sig to read*/
  float4 rho_p_c_vsigi;
  /*! Particle Internal energy u, alpha constants for visc and diff */
  float3 u_alphavisc_alphadiff;

  int2 cjs_cje;

} part_aos_f4_f_send;

/*Container for particle data requierd for force calcs*/
typedef struct part_aos_f4_f_recv {

  /*Result: Values output to global GPU memory*/
  /* change of u and h with dt, v_sig and returned value of
   * minimum neighbour timebin */
  float4 udt_hdt_vsig_mintimebin_ngb;
  /*Particle acceleration vector*/
  float3 a_hydro;

} part_aos_f4_f_recv;

/*Container for particle data requierd for gradient calcs*/
typedef struct part_aos_g {

  /*! Particle position. */
  double x_p;
  double y_p;
  double z_p;

  /*! Particle velocity. */
  float ux;
  float uy;
  float uz;
  /*! Particle mass. */
  float mass;
  /*! Particle smoothing length. */
  float h;
  /*! Particle density. */
  float rho;

  /* viscosity information */
  float visc_alpha;
  float laplace_u;
  float alpha_visc_max_ngb;
  float v_sig;

  float u;

  float soundspeed;

  /* timestep stuff */
  /*! Time-step length */
  int time_bin;
} part_aos_g;

/*Container for particle data requierd for gradient calcs*/
typedef struct part_aos_f4_g {

  /*! Particle position & smoothing length */
  float4 x_h;

  /*! Particle velocity and mass */
  float4 ux_m;

  /*! Particle density alpha visc internal energy u and speed of sound c */
  float4 rho_avisc_u_c;

  /* viscosity information results */
  float3 vsig_lapu_aviscmax_empty;

} part_aos_f4_g;

/*Container for particle data requierd for gradient calcs*/
typedef struct part_aos_f4_g_send {

  /*! Particle position & smoothing length */
  float4 x_h;

  /*! Particle velocity and mass */
  float4 ux_m;

  /*! Particle density alpha visc internal energy u and speed of sound c */
  float4 rho_avisc_u_c;

  /* viscosity information results */
  float3 vsig_lapu_aviscmax;

  /*Data for cell start and end*/
  int2 cjs_cje;

} part_aos_f4_g_send;

/*Container for particle data requierd for gradient calcs*/
typedef struct part_aos_f4_g_recv {

  /* viscosity information results */
  float3 vsig_lapu_aviscmax;

} part_aos_f4_g_recv;


#ifdef __WITH_CUDA
}
#endif

#endif // PART_GPU_H
