#ifndef CELL_GPU_H
#define CELL_GPU_H
/* Config parameters. */
#include "../config.h"
typedef int8_t timebin_t;
struct xpart_gpu {
  /*! Offset between current position and position at last tree rebuild. */
  float x_diff[3];
  /*! Offset between the current position and position at the last sort. */
  float x_diff_sort[3];
  /*! Velocity at the last full step. */
  float v_full[3];
  /*! Internal energy at the last full step. */
  float u_full;
};
struct part_gpu {
  /*Task ID*/
  int tid;
  /*! Particle unique ID. */
  long long id;
  /*! Pointer to corresponding gravity part. */
  //	struct gpu_gpart* gpart;
  /*! Particle position. */
  float x[3];
  /*! Particle predicted velocity. */
  float v[3];
  /*! Particle acceleration. */
  float a_hydro[3];
  /*! Particle mass. */
  float mass;
  /*! Particle smoothing length. */
  float h;
  /*! Particle internal energy. */
  float u;
  /*! Time derivative of the internal energy. */
  float u_dt;
  /*! Particle density. */
  float rho;
  /*! Kernel summation (For testing/debugging). */
  float SPH_sum;

  /* Cell information */
  /*! The cell location on the grid (corner nearest to the origin). */
  float loc[3];
  /*! The cell dimensions. */
  float width[3];
  float h_max;
  int count;
  /* Density information */

  /*! Neighbour number count. */
  float wcount;

  /*! Derivative of the neighbour number with respect to h. */
  float wcount_dh;

  /*! Derivative of density with respect to h */
  float rho_dh;

  /*! Particle velocity curl. */
  float rot_v[3];

  /* viscosity information */

  /*! Particle velocity divergence */
  float div_v;

  /*! Particle velocity divergence from previous step */
  float div_v_previous_step;

  /*! Artificial viscosity parameter */
  float alpha_visc;

  /*! Signal velocity */
  float v_sig;

  /* thermal diffusion information  */

  /*! del^2 u, a smoothed quantity */
  float laplace_u;

  /*! Thermal diffusion coefficient */
  float alpha_diff;

  /* force information  */

  /*! "Grad h" term -- only partial in P-U */
  float f;

  /*! Particle soundspeed. */
  float soundspeed;

  /*! Time derivative of smoothing length  */
  float h_dt;

  /*! Balsara switch */
  float balsara;

  /*! Particle pressure. */
  float pressure;
  /*! Maximal alpha (viscosity) over neighbours */
  float alpha_visc_max_ngb;

  /* timestep stuff */

  /*! Time-step length */
  timebin_t time_bin;

  /*all part of struct timestep_limiter_data, we had to destruct it
   as GPUs don't like pointer chasing especially when memcpying*/
  /* Need waking-up ? */
  timebin_t wakeup;

  /*! Minimal time-bin across all neighbours */
  timebin_t min_ngb_time_bin;

  /* Do we want this particle to be synched back on the time-line? */
  char to_be_synchronized;
};

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

struct task_cell {
  struct part_gpu *parts;
};
// struct parts_gpu_SoA{
//	struct task_cell *tasks;
// };

struct cell_hydro_gpu {
  //	struct part_gpu *parts;
  //	struct xpart_gpu *xparts;
  float h_max;
  int count;
};
struct cell_gpu {
  /*! The cell location on the grid (corner nearest to the origin). */
  float loc[3];
  /*! The cell dimensions. */
  float width[3];
  /*Details of contents (particles) and properties*/
  struct cell_hydro_gpu hydro;
};
struct cell_gpu_flat {
  /*! The cell location on the grid (corner nearest to the origin). */
  float loc[3];
  /*! The cell dimensions. */
  float width[3];
  float h_max;
  int count;
};

struct cells_gpu_flat {
  float *locx;
  float *locy;
  float *locz;
  /*! The cell dimensions. */
  float *widthx;
  float *widthy;
  float *widthz;
  /*! The cell location on the grid (corner nearest to the origin). */
  /*	float *loc[3];*/
  /*! The cell dimensions. */
  /*	float *width[3];*/
  float *h_max;
  int *count;
};

struct cells_gpu_flat_test {
  float *locx;
};

#endif // CELL_GPU_H
