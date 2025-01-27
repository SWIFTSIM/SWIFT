/*
 * Data_and_GPU_prep_functions.cu
 *
 *  Created on: 17 Apr 2022
 *      Author: abouzied
 */

/*ifdef WITH_CUDA prevents name mangling. C code sees exact names
 of functions rather than mangled template names produced by C++*/
// #ifdef WITH_CUDA
//	extern "C"{
// #endif

// #include "cuda/cuda_headers.h"
// #include "device_functions.h"
// #include "cuda/cell_gpu.h"
#include <cuda_profiler_api.h>
#include <vector.h>
// #include "../config.h"

void populate_parts_list(struct cell *ci, struct part_gpu *parts) {
  ////////////////////////////////////////////
  ///*****Copy variables for cell i (self interaction)*****/
  int count = ci->hydro.count;

  //	   fprintf(stderr,"Tester 111\n");
  for (int p = 0; p < count; p++) {

    parts[p].id = ci->hydro.parts[p].id;

    //		   fprintf(stderr,"Tester 222\n");
    parts[p].count = count;
    parts[p].h_max = ci->hydro.h_max;

    for (int d = 0; d < 3; d++) {
      parts[p].x[d] = ci->hydro.parts[p].x[d];
      parts[p].v[d] = ci->hydro.parts[p].v[d];
      parts[p].a_hydro[d] = ci->hydro.parts[p].a_hydro[d];
      parts[p].loc[d] = ci->loc[d];
    }
    parts[p].mass = ci->hydro.parts[p].mass;
    parts[p].h = ci->hydro.parts[p].h;
    parts[p].u = ci->hydro.parts[p].u;
    parts[p].u_dt = ci->hydro.parts[p].u_dt;
    parts[p].rho = ci->hydro.parts[p].rho;
    parts[p].div_v = ci->hydro.parts[p].viscosity.div_v;
    parts[p].div_v_previous_step =
        ci->hydro.parts[p].viscosity.div_v_previous_step;
    parts[p].alpha_visc = ci->hydro.parts[p].viscosity.alpha;
    parts[p].v_sig = ci->hydro.parts[p].viscosity.v_sig;
    parts[p].laplace_u = ci->hydro.parts[p].diffusion.laplace_u;
    parts[p].alpha_diff = ci->hydro.parts[p].diffusion.alpha;
    parts[p].f = ci->hydro.parts[p].force.f;
    parts[p].soundspeed = ci->hydro.parts[p].force.soundspeed;
    parts[p].h_dt = ci->hydro.parts[p].force.h_dt;
    parts[p].balsara = ci->hydro.parts[p].force.balsara;
    parts[p].pressure = ci->hydro.parts[p].force.pressure;
    parts[p].time_bin = ci->hydro.parts[p].time_bin;
    parts[p].wakeup = ci->hydro.parts[p].limiter_data.wakeup;
    parts[p].min_ngb_time_bin =
        ci->hydro.parts[p].limiter_data.min_ngb_time_bin;
    parts[p].to_be_synchronized =
        ci->hydro.parts[p].limiter_data.to_be_synchronized;
    parts[p].wcount = ci->hydro.parts[p].density.wcount;
    parts[p].wcount_dh = ci->hydro.parts[p].density.wcount_dh;
    parts[p].rho_dh = ci->hydro.parts[p].density.rho_dh;
    parts[p].div_v = ci->hydro.parts[p].viscosity.div_v;
    parts[p].rot_v[0] = ci->hydro.parts[p].density.rot_v[0];
    parts[p].rot_v[1] = ci->hydro.parts[p].density.rot_v[1];
    parts[p].rot_v[2] = ci->hydro.parts[p].density.rot_v[2];
    parts[p].SPH_sum = 0.f;
  }
}

void populate_parts_list_soa(
    int count_all_parts, struct cell *ci, int first_part_tmp, int count,
    int tid, int *tid_p, long long *id, double *x_p, double *y_p, double *z_p,
    float *ux, float *uy, float *uz, float *a_hydrox, float *a_hydroy,
    float *a_hydroz, float *mass, float *h, float *u, float *u_dt, float *rho,
    float *SPH_sum, float *locx, float *locy, float *locz, float *widthx,
    float *widthy, float *widthz, float *h_max, int *count_p, float *wcount,
    float *wcount_dh, float *rho_dh, float *rot_u, float *rot_v, float *rot_w,
    float *div_v, float *div_v_previous_step, float *alpha_visc, float *v_sig,
    float *laplace_u, float *alpha_diff, float *f, float *soundspeed,
    float *h_dt, float *balsara, float *pressure, float *alpha_visc_max_ngb,
    timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
    char *to_be_synchronized) {
  ////////////////////////////////////////////
  struct part *ptmps;
  ptmps = ci->hydro.parts;
  //	   fprintf(stderr,"Tester 111\n");
#pragma unroll
  for (int p = 0; p < count; p++) {
    int p_gid = p + first_part_tmp;
    //		    if(p_gid>=count_all_parts){
    //		    	fprintf(stderr,"p>all parts");
    //		    	exit(0);
    //		    }
    id[p_gid] = ptmps[p].id;
    count_p[p_gid] = count;
    tid_p[p_gid] = tid;
    h_max[p_gid] = ci->hydro.h_max;
    x_p[p_gid] = ptmps[p].x[0];
    y_p[p_gid] = ptmps[p].x[1];
    z_p[p_gid] = ptmps[p].x[2];
    ux[p_gid] = ptmps[p].v[0];
    uy[p_gid] = ptmps[p].v[1];
    uz[p_gid] = ptmps[p].v[2];
    a_hydrox[p_gid] = ptmps[p].a_hydro[0];
    a_hydroy[p_gid] = ptmps[p].a_hydro[1];
    a_hydroz[p_gid] = ptmps[p].a_hydro[2];
    locx[p_gid] = ci->loc[0];
    locy[p_gid] = ci->loc[1];
    locz[p_gid] = ci->loc[2];

    mass[p_gid] = ptmps[p].mass;
    h[p_gid] = ptmps[p].h;
    u[p_gid] = ptmps[p].u;
    u_dt[p_gid] = ptmps[p].u_dt;
    rho[p_gid] = ptmps[p].rho;
    div_v[p_gid] = ptmps[p].viscosity.div_v;
    div_v_previous_step[p_gid] = ptmps[p].viscosity.div_v_previous_step;
    alpha_visc[p_gid] = ptmps[p].viscosity.alpha;
    v_sig[p_gid] = ptmps[p].viscosity.v_sig;
    laplace_u[p_gid] = ptmps[p].diffusion.laplace_u;
    alpha_diff[p_gid] = ptmps[p].diffusion.alpha;
    f[p_gid] = ptmps[p].force.f;
    soundspeed[p_gid] = ptmps[p].force.soundspeed;
    h_dt[p_gid] = ptmps[p].force.h_dt;
    balsara[p_gid] = ptmps[p].force.balsara;
    pressure[p_gid] = ptmps[p].force.pressure;
    time_bin[p_gid] = ptmps[p].time_bin;
    wakeup[p_gid] = ptmps[p].limiter_data.wakeup;
    min_ngb_time_bin[p_gid] = ptmps[p].limiter_data.min_ngb_time_bin;
    to_be_synchronized[p_gid] = ptmps[p].limiter_data.to_be_synchronized;
    wcount[p_gid] = ptmps[p].density.wcount;
    wcount_dh[p_gid] = ptmps[p].density.wcount_dh;
    rho_dh[p_gid] = ptmps[p].density.rho_dh;
    div_v[p_gid] = ptmps[p].viscosity.div_v;
    rot_u[p_gid] = ptmps[p].density.rot_v[0];
    rot_v[p_gid] = ptmps[p].density.rot_v[1];
    rot_w[p_gid] = ptmps[p].density.rot_v[2];
    SPH_sum[p_gid] = 0.f;
    //			fprintf(stderr,"tid is %i\n",tid_p[p]);
    //			fprintf(stderr,"Tester 222, count=%i, p=%i\n", count,
    // id[p_gid]);
  }
}

void pack_data_soa(int count_all_parts, struct cell *ci, int first_part_tmp,
                   int count, int tid, int *tid_p, long long *id, double *x_p,
                   double *y_p, double *z_p, float *ux, float *uy, float *uz,
                   float *a_hydrox, float *a_hydroy, float *a_hydroz,
                   float *mass, float *h, float *u, float *u_dt, float *rho,
                   float *SPH_sum, float *locx, float *locy, float *locz,
                   float *widthx, float *widthy, float *widthz, float *h_max,
                   int *count_p, float *wcount, float *wcount_dh, float *rho_dh,
                   float *rot_u, float *rot_v, float *rot_w, float *div_v,
                   float *div_v_previous_step, float *alpha_visc, float *v_sig,
                   float *laplace_u, float *alpha_diff, float *f,
                   float *soundspeed, float *h_dt, float *balsara,
                   float *pressure, float *alpha_visc_max_ngb,
                   timebin_t *time_bin, timebin_t *wakeup,
                   timebin_t *min_ngb_time_bin, char *to_be_synchronized) {
  ////////////////////////////////////////////
  struct part *ptmps;
  ptmps = ci->hydro.parts;
  //	   fprintf(stderr,"Tester 111\n");
#pragma unroll
  for (int p = 0; p < count; p++) {
    int p_gid = p + first_part_tmp;
    //		    if(p_gid>=count_all_parts){
    //		    	fprintf(stderr,"p>all parts");
    //		    	exit(0);
    //		    }
    id[p_gid] = ptmps[p].id;
    count_p[p_gid] = count;
    tid_p[p_gid] = tid;
    h_max[p_gid] = ci->hydro.h_max;
    x_p[p_gid] = ptmps[p].x[0];
    y_p[p_gid] = ptmps[p].x[1];
    z_p[p_gid] = ptmps[p].x[2];
    ux[p_gid] = ptmps[p].v[0];
    uy[p_gid] = ptmps[p].v[1];
    uz[p_gid] = ptmps[p].v[2];
    a_hydrox[p_gid] = ptmps[p].a_hydro[0];
    a_hydroy[p_gid] = ptmps[p].a_hydro[1];
    a_hydroz[p_gid] = ptmps[p].a_hydro[2];
    locx[p_gid] = ci->loc[0];
    locy[p_gid] = ci->loc[1];
    locz[p_gid] = ci->loc[2];

    mass[p_gid] = ptmps[p].mass;
    h[p_gid] = ptmps[p].h;
    u[p_gid] = ptmps[p].u;
    u_dt[p_gid] = ptmps[p].u_dt;
    rho[p_gid] = ptmps[p].rho;
    div_v[p_gid] = ptmps[p].viscosity.div_v;
    div_v_previous_step[p_gid] = ptmps[p].viscosity.div_v_previous_step;
    alpha_visc[p_gid] = ptmps[p].viscosity.alpha;
    v_sig[p_gid] = ptmps[p].viscosity.v_sig;
    laplace_u[p_gid] = ptmps[p].diffusion.laplace_u;
    alpha_diff[p_gid] = ptmps[p].diffusion.alpha;
    f[p_gid] = ptmps[p].force.f;
    soundspeed[p_gid] = ptmps[p].force.soundspeed;
    h_dt[p_gid] = ptmps[p].force.h_dt;
    balsara[p_gid] = ptmps[p].force.balsara;
    pressure[p_gid] = ptmps[p].force.pressure;
    time_bin[p_gid] = ptmps[p].time_bin;
    wakeup[p_gid] = ptmps[p].limiter_data.wakeup;
    min_ngb_time_bin[p_gid] = ptmps[p].limiter_data.min_ngb_time_bin;
    to_be_synchronized[p_gid] = ptmps[p].limiter_data.to_be_synchronized;
    wcount[p_gid] = ptmps[p].density.wcount;
    wcount_dh[p_gid] = ptmps[p].density.wcount_dh;
    rho_dh[p_gid] = ptmps[p].density.rho_dh;
    div_v[p_gid] = ptmps[p].viscosity.div_v;
    rot_u[p_gid] = ptmps[p].density.rot_v[0];
    rot_v[p_gid] = ptmps[p].density.rot_v[1];
    rot_w[p_gid] = ptmps[p].density.rot_v[2];
    SPH_sum[p_gid] = 0.f;
    //			fprintf(stderr,"tid is %i\n",tid_p[p]);
    //			fprintf(stderr,"Tester 222, count=%i, p=%i\n", count,
    // id[p_gid]);
  }
}

// #ifdef WITH_CUDA
//	}
// #endif
