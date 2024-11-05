// #include "active.h"
// #include <cuda_runtime.h>
// #include <vector>
// #include "cuda/cell_gpu.h"
// #include "runner_gpu_functions.cuh"
/* This object's header. */
#include "runner.h"
/* Local headers. */
#include "active.h"
#include "engine.h"
#include "runner_gpu_pack_functions.h"
#include "scheduler.h"
#include "space_getsid.h"
#include "timers.h"

// #ifdef WITHCUDA
// extern "C" {
// #endif

void runner_doself1_gpu_pack_neat(struct runner *r, struct cell *c,
                                  struct part_soa parts_soa_buffer, int timer,
                                  int *pack_length, int tid,
                                  int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat(c, parts_soa_buffer, tid, local_pack_position, count);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_doself1_gpu_pack_neat_aos(struct runner *r, struct cell *c,
                                      struct part_aos *parts_aos_buffer,
                                      int timer, int *pack_length, int tid,
                                      int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    error("0");
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos(c, parts_aos_buffer, tid, local_pack_position, count);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_doself1_gpu_pack_neat_aos_f4(
    struct runner *r, struct cell *__restrict__ c,
    struct part_aos_f4_send *__restrict__ parts_aos_buffer, int timer,
    int *pack_length, int tid, int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    error("0");
  }
#endif
  int2 frst_lst_prts = {local_pack_position, local_pack_position + count};
  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_f4(c, parts_aos_buffer, tid, local_pack_position, count,
                   frst_lst_prts);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_doself1_gpu_pack_neat_aos_g(struct runner *r, struct cell *c,
                                        struct part_aos_g *parts_aos_buffer,
                                        int timer, int *pack_length, int tid,
                                        int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_g(c, parts_aos_buffer, tid, local_pack_position, count);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_doself1_gpu_pack_neat_aos_f4_g(
    struct runner *r, struct cell *c,
    struct part_aos_f4_g_send *parts_aos_buffer, int timer, int *pack_length,
    int tid, int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_f4_g(c, parts_aos_buffer, tid, local_pack_position, count);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_doself1_gpu_pack_neat_aos_f(struct runner *r, struct cell *c,
                                        struct part_aos_f *parts_aos_buffer,
                                        int timer, int *pack_length, int tid,
                                        int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_f(c, parts_aos_buffer, tid, local_pack_position, count);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_doself1_gpu_pack_neat_aos_f4_f(
    struct runner *r, struct cell *restrict c,
    struct part_aos_f4_f_send *restrict parts_aos_buffer, int timer,
    int *pack_length, int tid, int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_f4_f(c, parts_aos_buffer, tid, local_pack_position, count);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void pack_neat(struct cell *c, struct part_soa parts_soa_buffer, int tid,
               int local_pack_position, int count) {

  const struct part *ptmps;
  ptmps = c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    const struct part p = ptmps[i];
    /*Data to be copied to GPU*/
    parts_soa_buffer.x_p[id_in_pack] = p.x[0];
    parts_soa_buffer.y_p[id_in_pack] = p.x[1];
    parts_soa_buffer.z_p[id_in_pack] = p.x[2];
    parts_soa_buffer.tid_p[id_in_pack] = tid;
    parts_soa_buffer.ux[id_in_pack] = p.v[0];
    parts_soa_buffer.uy[id_in_pack] = p.v[1];
    parts_soa_buffer.uz[id_in_pack] = p.v[2];
    parts_soa_buffer.locx[id_in_pack] = c->loc[0];
    parts_soa_buffer.locy[id_in_pack] = c->loc[1];
    parts_soa_buffer.locz[id_in_pack] = c->loc[2];
    parts_soa_buffer.mass[id_in_pack] = p.mass;
    parts_soa_buffer.h[id_in_pack] = p.h;
    //    parts_soa_buffer.time_bin[id_in_pack] = p.time_bin;
    /*Initialise sums to zero before CPU/GPU copy*/
    parts_soa_buffer.rho[id_in_pack] = 0.f;        // p.rho;
    parts_soa_buffer.rho_dh[id_in_pack] = 0.f;     // p.density.rho_dh;
    parts_soa_buffer.wcount[id_in_pack] = 0.f;     // p.density.wcount;
    parts_soa_buffer.wcount_dh[id_in_pack] = 0.f;  // p.density.wcount_dh;
    parts_soa_buffer.div_v[id_in_pack] = 0.f;      // p.viscosity.div_v;
    parts_soa_buffer.rot_ux[id_in_pack] = 0.f;     // p.density.rot_v[0];
    parts_soa_buffer.rot_uy[id_in_pack] = 0.f;     // p.density.rot_v[1];
    parts_soa_buffer.rot_uz[id_in_pack] = 0.f;     // p.density.rot_v[2];
  }
}

void pack_neat_aos(struct cell *c, struct part_aos *parts_aos_buffer, int tid,
                   int local_pack_position, int count) {

  const struct part *ptmps;
  ptmps = c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    const struct part p = ptmps[i];
    /*Data to be copied to GPU*/
    parts_aos_buffer[id_in_pack].x_p = p.x[0];
    parts_aos_buffer[id_in_pack].y_p = p.x[1];
    parts_aos_buffer[id_in_pack].z_p = p.x[2];
    parts_aos_buffer[id_in_pack].ux = p.v[0];
    parts_aos_buffer[id_in_pack].uy = p.v[1];
    parts_aos_buffer[id_in_pack].uz = p.v[2];
    parts_aos_buffer[id_in_pack].mass = p.mass;
    parts_aos_buffer[id_in_pack].h = p.h;
    parts_aos_buffer[id_in_pack].time_bin = 1000;  // p.time_bin;
    /*Initialise sums to zero before CPU/GPU copy*/
    parts_aos_buffer[id_in_pack].rho = 0.f;        // p.rho;
    parts_aos_buffer[id_in_pack].rho_dh = 0.f;     // p.density.rho_dh;
    parts_aos_buffer[id_in_pack].wcount = 0.f;     // p.density.wcount;
    parts_aos_buffer[id_in_pack].wcount_dh = 0.f;  // p.density.wcount_dh;
    parts_aos_buffer[id_in_pack].div_v = 0.f;      // p.viscosity.div_v;
    parts_aos_buffer[id_in_pack].rot_ux = 0.f;     // p.density.rot_v[0];
    parts_aos_buffer[id_in_pack].rot_uy = 0.f;     // p.density.rot_v[1];
    parts_aos_buffer[id_in_pack].rot_uz = 0.f;     // p.density.rot_v[2];
  }
}

void pack_neat_pair_aos(struct cell *c, struct part_aos *parts_aos_buffer,
                        int tid, int local_pack_position, int count,
                        float3 shift) {

  const struct part *ptmps;
  ptmps = c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    const struct part p = ptmps[i];
    /*Data to be copied to GPU*/
    parts_aos_buffer[id_in_pack].x_p = p.x[0] - shift.x;
    parts_aos_buffer[id_in_pack].y_p = p.x[1] - shift.y;
    parts_aos_buffer[id_in_pack].z_p = p.x[2] - shift.z;
    parts_aos_buffer[id_in_pack].ux = p.v[0];
    parts_aos_buffer[id_in_pack].uy = p.v[1];
    parts_aos_buffer[id_in_pack].uz = p.v[2];
    parts_aos_buffer[id_in_pack].mass = p.mass;
    parts_aos_buffer[id_in_pack].h = p.h;
    parts_aos_buffer[id_in_pack].time_bin = 1000;  // p.time_bin;
    /*Initialise sums to zero before CPU/GPU copy*/
    parts_aos_buffer[id_in_pack].rho = 0.f;        // p.rho;
    parts_aos_buffer[id_in_pack].rho_dh = 0.f;     // p.density.rho_dh;
    parts_aos_buffer[id_in_pack].wcount = 0.f;     // p.density.wcount;
    parts_aos_buffer[id_in_pack].wcount_dh = 0.f;  // p.density.wcount_dh;
    parts_aos_buffer[id_in_pack].div_v = 0.f;      // p.viscosity.div_v;
    parts_aos_buffer[id_in_pack].rot_ux = 0.f;     // p.density.rot_v[0];
    parts_aos_buffer[id_in_pack].rot_uy = 0.f;     // p.density.rot_v[1];
    parts_aos_buffer[id_in_pack].rot_uz = 0.f;     // p.density.rot_v[2];
  }
}

extern inline void pack_neat_pair_aos_f4(
    struct cell *__restrict c,
    struct part_aos_f4_send *__restrict parts_aos_buffer, int tid,
    const int local_pack_position, const int count, const float3 shift,
    const int2 cstarts) {
  /*Data to be copied to GPU*/
  for (int i = 0; i < count; i++) {
    const int id_in_pack = i + local_pack_position;
    parts_aos_buffer[id_in_pack].x_p_h.x = c->hydro.parts[i].x[0] - shift.x;
    parts_aos_buffer[id_in_pack].x_p_h.y = c->hydro.parts[i].x[1] - shift.y;
    parts_aos_buffer[id_in_pack].x_p_h.z = c->hydro.parts[i].x[2] - shift.z;
    parts_aos_buffer[id_in_pack].x_p_h.w = c->hydro.parts[i].h;
    parts_aos_buffer[id_in_pack].ux_m.x = c->hydro.parts[i].v[0];
    parts_aos_buffer[id_in_pack].ux_m.y = c->hydro.parts[i].v[1];
    parts_aos_buffer[id_in_pack].ux_m.z = c->hydro.parts[i].v[2];
    parts_aos_buffer[id_in_pack].ux_m.w = c->hydro.parts[i].mass;
    parts_aos_buffer[id_in_pack].cjs_cje.x = cstarts.x;
    parts_aos_buffer[id_in_pack].cjs_cje.y = cstarts.y;
  }
}

void pack_neat_aos_f4(struct cell *__restrict__ c,
                      struct part_aos_f4_send *__restrict__ parts_aos_buffer,
                      int tid, int local_pack_position, int count,
                      int2 frst_lst_prts) {

  struct part ptmps[count];
  memcpy(ptmps, (c->hydro.parts), count * sizeof(struct part));
  //  ptmps = c->hydro.parts;
  const float cellx = c->loc[0], celly = c->loc[1], cellz = c->loc[2];
  for (int i = 0; i < count; i++) {
    const int id_in_pack = i + local_pack_position;
    //    const struct part p = ptmps[i];
    /*Data to be copied to GPU*/
    parts_aos_buffer[id_in_pack].x_p_h.x = ptmps[i].x[0] - cellx;
    parts_aos_buffer[id_in_pack].x_p_h.y = ptmps[i].x[1] - celly;
    parts_aos_buffer[id_in_pack].x_p_h.z = ptmps[i].x[2] - cellz;
    parts_aos_buffer[id_in_pack].x_p_h.w = ptmps[i].h;
    parts_aos_buffer[id_in_pack].ux_m.x = ptmps[i].v[0];
    parts_aos_buffer[id_in_pack].ux_m.y = ptmps[i].v[1];
    parts_aos_buffer[id_in_pack].ux_m.z = ptmps[i].v[2];
    parts_aos_buffer[id_in_pack].ux_m.w = ptmps[i].mass;
    //    /*Initialise sums to zero before CPU/GPU copy*/
    //    const float4 zeroes = {0.0, 0.0, 0.0, 0.0};
    //    parts_aos_buffer[id_in_pack].rho_dh_wcount = zeroes;
    //    parts_aos_buffer[id_in_pack].rot_ux_div_v = zeroes;
  }
}

void pack_neat_aos_g(struct cell *c, struct part_aos_g *parts_aos_buffer,
                     int tid, int local_pack_position, int count) {

  const struct part *ptmps;
  ptmps = c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    const struct part p = ptmps[i];
    /*Data to be copied to GPU*/
    parts_aos_buffer[id_in_pack].x_p = p.x[0];
    parts_aos_buffer[id_in_pack].y_p = p.x[1];
    parts_aos_buffer[id_in_pack].z_p = p.x[2];
    parts_aos_buffer[id_in_pack].ux = p.v[0];
    parts_aos_buffer[id_in_pack].uy = p.v[1];
    parts_aos_buffer[id_in_pack].uz = p.v[2];
    parts_aos_buffer[id_in_pack].mass = p.mass;
    parts_aos_buffer[id_in_pack].h = p.h;
    parts_aos_buffer[id_in_pack].time_bin = 1000;
    parts_aos_buffer[id_in_pack].rho = p.rho;
    parts_aos_buffer[id_in_pack].visc_alpha = p.viscosity.alpha;
    parts_aos_buffer[id_in_pack].alpha_visc_max_ngb =
        p.force.alpha_visc_max_ngb;  // p.density.wcount_dh;
    parts_aos_buffer[id_in_pack].v_sig =
        p.viscosity.v_sig;  // p.viscosity.div_v;
    parts_aos_buffer[id_in_pack].soundspeed =
        p.force.soundspeed;                // p.density.rot_v[0];
    parts_aos_buffer[id_in_pack].u = p.u;  // p.density.rot_v[0];
    /*Initialise sums to zero before CPU/GPU copy*/
    parts_aos_buffer[id_in_pack].laplace_u = 0.f;  // p.density.wcount;
  }
}

void pack_neat_aos_f4_g(struct cell *c,
                        struct part_aos_f4_g_send *parts_aos_buffer, int tid,
                        int local_pack_position, int count) {

  const struct part *ptmps;
  ptmps = c->hydro.parts;
  const float cellx = c->loc[0], celly = c->loc[1], cellz = c->loc[2];
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    const struct part p = ptmps[i];
    /*Data to be copied to GPU*/
    parts_aos_buffer[id_in_pack].x_h.x = p.x[0] - cellx;
    parts_aos_buffer[id_in_pack].x_h.y = p.x[1] - celly;
    parts_aos_buffer[id_in_pack].x_h.z = p.x[2] - cellz;
    parts_aos_buffer[id_in_pack].x_h.w = p.h;
    parts_aos_buffer[id_in_pack].ux_m.x = p.v[0];
    parts_aos_buffer[id_in_pack].ux_m.y = p.v[1];
    parts_aos_buffer[id_in_pack].ux_m.z = p.v[2];
    parts_aos_buffer[id_in_pack].ux_m.w = p.mass;
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.x = p.rho;
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.y = p.viscosity.alpha;
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.z = p.u;  // p.density.rot_v[0];
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.w =
        p.force.soundspeed;  // p.density.rot_v[0];
  }
}

extern inline void pack_neat_pair_aos_f4_g(
    struct cell *__restrict c,
    struct part_aos_f4_g_send *__restrict parts_aos_buffer, int tid,
    const int local_pack_position, const int count, const float3 shift,
    const int2 cstarts) {
  /*Data to be copied to GPU*/
  for (int i = 0; i < count; i++) {
    const int id_in_pack = i + local_pack_position;
    parts_aos_buffer[id_in_pack].x_h.x = c->hydro.parts[i].x[0] - shift.x;
    parts_aos_buffer[id_in_pack].x_h.y = c->hydro.parts[i].x[1] - shift.y;
    parts_aos_buffer[id_in_pack].x_h.z = c->hydro.parts[i].x[2] - shift.z;
    parts_aos_buffer[id_in_pack].x_h.w = c->hydro.parts[i].h;
    parts_aos_buffer[id_in_pack].ux_m.x = c->hydro.parts[i].v[0];
    parts_aos_buffer[id_in_pack].ux_m.y = c->hydro.parts[i].v[1];
    parts_aos_buffer[id_in_pack].ux_m.z = c->hydro.parts[i].v[2];
    parts_aos_buffer[id_in_pack].ux_m.w = c->hydro.parts[i].mass;
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.x = c->hydro.parts[i].rho;
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.y =
        c->hydro.parts[i].viscosity.alpha;
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.z =
        c->hydro.parts[i].u;  // p.density.rot_v[0];
    parts_aos_buffer[id_in_pack].rho_avisc_u_c.w =
        c->hydro.parts[i].force.soundspeed;  // p.density.rot_v[0];
    parts_aos_buffer[id_in_pack].cjs_cje.x = cstarts.x;
    parts_aos_buffer[id_in_pack].cjs_cje.y = cstarts.y;
  }
}

void pack_neat_aos_f(struct cell *c, struct part_aos_f *parts_aos, int tid,
                     int local_pack_position, int count) {

  const struct part *ptmps;
  ptmps = c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    const struct part p = ptmps[i];
    /*Data to be copied to GPU*/
    parts_aos[id_in_pack].x_p = p.x[0];
    parts_aos[id_in_pack].y_p = p.x[1];
    parts_aos[id_in_pack].z_p = p.x[2];
    parts_aos[id_in_pack].ux = p.v[0];
    parts_aos[id_in_pack].uy = p.v[1];
    parts_aos[id_in_pack].uz = p.v[2];
    parts_aos[id_in_pack].mass = p.mass;
    parts_aos[id_in_pack].h = p.h;
    parts_aos[id_in_pack].time_bin = p.time_bin;
    parts_aos[id_in_pack].min_ngb_time_bin = p.limiter_data.min_ngb_time_bin;
    parts_aos[id_in_pack].rho = p.rho;
    parts_aos[id_in_pack].pressure = p.force.pressure;
    parts_aos[id_in_pack].soundspeed = p.force.soundspeed;
    parts_aos[id_in_pack].f = p.force.f;
    parts_aos[id_in_pack].balsara = p.force.balsara;
    parts_aos[id_in_pack].alpha_visc = p.viscosity.alpha;
    parts_aos[id_in_pack].a_hydrox = 0.0;
    parts_aos[id_in_pack].a_hydroy = 0.0;
    parts_aos[id_in_pack].a_hydroz = 0.0;
    parts_aos[id_in_pack].alpha_diff = p.diffusion.alpha;
    parts_aos[id_in_pack].u = p.u;
    parts_aos[id_in_pack].u_dt = 0.0;
    parts_aos[id_in_pack].h_dt = 0.0;
    /*Initialise sums to zero before CPU/GPU copy*/
    parts_aos[id_in_pack].v_sig = p.viscosity.v_sig;
  }
}

void pack_neat_aos_f4_f(const struct cell *restrict c,
                        struct part_aos_f4_f_send *restrict parts_aos, int tid,
                        int local_pack_position, int count) {

  //  const struct part *restrict ptmps;
  //  ptmps = c->hydro.parts;
  const int pp = local_pack_position;
  const float cellx = c->loc[0];
  const float celly = c->loc[1];
  const float cellz = c->loc[2];
  /*Data to be copied to GPU local memory*/
  for (int i = 0; i < count; i++) {
    parts_aos[i + pp].x_h.x = c->hydro.parts[i].x[0] - cellx;
    parts_aos[i + pp].x_h.y = c->hydro.parts[i].x[1] - celly;
    parts_aos[i + pp].x_h.z = c->hydro.parts[i].x[2] - cellz;
    parts_aos[i + pp].x_h.w = c->hydro.parts[i].h;
  }
  for (int i = 0; i < count; i++) {
    parts_aos[i + pp].ux_m.x = c->hydro.parts[i].v[0];
    parts_aos[i + pp].ux_m.y = c->hydro.parts[i].v[1];
    parts_aos[i + pp].ux_m.z = c->hydro.parts[i].v[2];
    parts_aos[i + pp].ux_m.w = c->hydro.parts[i].mass;
  }
  for (int i = 0; i < count; i++) {
    parts_aos[i + pp].f_bals_timebin_mintimebin_ngb.x =
        c->hydro.parts[i].force.f;
    parts_aos[i + pp].f_bals_timebin_mintimebin_ngb.y =
        c->hydro.parts[i].force.balsara;
    parts_aos[i + pp].f_bals_timebin_mintimebin_ngb.z =
        c->hydro.parts[i].time_bin;
    parts_aos[i + pp].f_bals_timebin_mintimebin_ngb.w =
        c->hydro.parts[i].limiter_data.min_ngb_time_bin;
  }
  for (int i = 0; i < count; i++) {
    parts_aos[i + pp].rho_p_c_vsigi.x = c->hydro.parts[i].rho;
    parts_aos[i + pp].rho_p_c_vsigi.y = c->hydro.parts[i].force.pressure;
    parts_aos[i + pp].rho_p_c_vsigi.z = c->hydro.parts[i].force.soundspeed;
    parts_aos[i + pp].rho_p_c_vsigi.w = c->hydro.parts[i].viscosity.v_sig;
  }
  for (int i = 0; i < count; i++) {
    parts_aos[i + pp].u_alphavisc_alphadiff.x = c->hydro.parts[i].u;
    parts_aos[i + pp].u_alphavisc_alphadiff.y =
        c->hydro.parts[i].viscosity.alpha;
    parts_aos[i + pp].u_alphavisc_alphadiff.z =
        c->hydro.parts[i].diffusion.alpha;
  }
}

extern inline void pack_neat_pair_aos_f4_f(
    struct cell *__restrict c, struct part_aos_f4_f_send *__restrict parts_aos,
    int tid, const int local_pack_position, const int count, const float3 shift,
    const int2 cstarts) {
  //  const struct part *restrict ptmps;
  //  ptmps = c->hydro.parts;
  const int pp = local_pack_position;
  /*Data to be copied to GPU local memory*/
  for (int i = 0; i < count; i++) {
    const int id = i + pp;
    parts_aos[id].x_h.x = c->hydro.parts[i].x[0] - shift.x;
    parts_aos[id].x_h.y = c->hydro.parts[i].x[1] - shift.y;
    parts_aos[id].x_h.z = c->hydro.parts[i].x[2] - shift.z;
    parts_aos[id].x_h.w = c->hydro.parts[i].h;
    parts_aos[id].ux_m.x = c->hydro.parts[i].v[0];
    parts_aos[id].ux_m.y = c->hydro.parts[i].v[1];
    parts_aos[id].ux_m.z = c->hydro.parts[i].v[2];
    parts_aos[id].ux_m.w = c->hydro.parts[i].mass;
    parts_aos[id].f_bals_timebin_mintimebin_ngb.x = c->hydro.parts[i].force.f;
    parts_aos[id].f_bals_timebin_mintimebin_ngb.y =
        c->hydro.parts[i].force.balsara;
    parts_aos[id].f_bals_timebin_mintimebin_ngb.z = c->hydro.parts[i].time_bin;
    parts_aos[id].f_bals_timebin_mintimebin_ngb.w =
        c->hydro.parts[i].limiter_data.min_ngb_time_bin;
    parts_aos[id].rho_p_c_vsigi.x = c->hydro.parts[i].rho;
    parts_aos[id].rho_p_c_vsigi.y = c->hydro.parts[i].force.pressure;
    parts_aos[id].rho_p_c_vsigi.z = c->hydro.parts[i].force.soundspeed;
    parts_aos[id].rho_p_c_vsigi.w = c->hydro.parts[i].viscosity.v_sig;
    parts_aos[id].u_alphavisc_alphadiff.x = c->hydro.parts[i].u;
    parts_aos[id].u_alphavisc_alphadiff.y = c->hydro.parts[i].viscosity.alpha;
    parts_aos[id].u_alphavisc_alphadiff.z = c->hydro.parts[i].diffusion.alpha;
    parts_aos[id].cjs_cje.x = cstarts.x;
    parts_aos[id].cjs_cje.y = cstarts.y;
  }
}

void runner_doself1_gpu_unpack_neat(struct runner *r, struct cell *c,
                                    struct part_soa parts_soa_buffer, int timer,
                                    int *pack_length, int tid,
                                    int count_max_parts_tmp, struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) {
    message("Inactive cell\n");
    return;
  }
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count, e);
  }
#endif

  /* Copy particle data from CPU buffers to cells */
  unpack_neat(c, parts_soa_buffer, tid, local_pack_position, count, e);
  // Increment pack length accordingly
  (*pack_length) += count;
}

void runner_doself1_gpu_unpack_neat_aos(struct runner *r, struct cell *c,
                                        struct part_aos *parts_aos_buffer,
                                        int timer, int *pack_length, int tid,
                                        int count_max_parts_tmp,
                                        struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) {
    message("Inactive cell\n");
    return;
  }
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count, e);
  }
#endif

  /* Copy particle data from CPU buffers to cells */
  unpack_neat_aos(c, parts_aos_buffer, tid, local_pack_position, count, e);
  // Increment pack length accordingly
  (*pack_length) += count;
}

void runner_doself1_gpu_unpack_neat_aos_f4(
    struct runner *r, struct cell *c, struct part_aos_f4_recv *parts_aos_buffer,
    int timer, int *pack_length, int tid, int count_max_parts_tmp,
    struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) {
    message("Inactive cell\n");
    return;
  }
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count, e);
  }
#endif

  /* Copy particle data from CPU buffers to cells */
  unpack_neat_aos_f4(c, parts_aos_buffer, tid, local_pack_position, count, e);
  // Increment pack length accordingly
  (*pack_length) += count;
}

void runner_doself1_gpu_unpack_neat_aos_g(struct runner *r, struct cell *c,
                                          struct part_aos_g *parts_aos_buffer,
                                          int timer, int *pack_length, int tid,
                                          int count_max_parts_tmp,
                                          struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) {
    message("Inactive cell\n");
    return;
  }
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count, e);
  }
#endif

  /* Copy particle data from CPU buffers to cells */
  unpack_neat_aos_g(c, parts_aos_buffer, tid, local_pack_position, count, e);
  // Increment pack length accordingly
  (*pack_length) += count;
}

void runner_doself1_gpu_unpack_neat_aos_f4_g(
    struct runner *r, struct cell *c,
    struct part_aos_f4_g_recv *parts_aos_buffer, int timer, int *pack_length,
    int tid, int count_max_parts_tmp, struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) {
    message("Inactive cell\n");
    return;
  }
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count, e);
  }
#endif

  /* Copy particle data from CPU buffers to cells */
  unpack_neat_aos_f4_g(c, parts_aos_buffer, tid, local_pack_position, count, e);
  // Increment pack length accordingly
  (*pack_length) += count;
}

void runner_doself1_gpu_unpack_neat_aos_f(struct runner *r, struct cell *c,
                                          struct part_aos_f *parts_aos_buffer,
                                          int timer, int *pack_length, int tid,
                                          int count_max_parts_tmp,
                                          struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) {
    message("Inactive cell\n");
    return;
  }
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count, e);
  }
#endif

  /* Copy particle data from CPU buffers to cells */
  unpack_neat_aos_f(c, parts_aos_buffer, tid, local_pack_position, count, e);
  // Increment pack length accordingly
  (*pack_length) += count;
}

void runner_doself1_gpu_unpack_neat_aos_f4_f(
    struct runner *r, struct cell *c,
    struct part_aos_f4_f_recv *parts_aos_buffer, int timer, int *pack_length,
    int tid, int count_max_parts_tmp, struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) {
    message("Inactive cell\n");
    return;
  }
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count, e);
  }
#endif

  /* Copy particle data from CPU buffers to cells */
  unpack_neat_aos_f4_f(c, parts_aos_buffer, tid, local_pack_position, count, e);
  // Increment pack length accordingly
  (*pack_length) += count;
}

void unpack_neat(struct cell *c, struct part_soa parts_soa_buffer, int tid,
                 int local_pack_position, int count, struct engine *e) {

  //  struct part *ptmps;
  //  ptmps=c->hydro.parts;

  //  memcpy(&rho[0], &parts_soa_buffer.rho[local_pack_position], count *
  //  sizeof(float));
  //	fprintf(stderr, "count %i\n", count);
  //  memcpy(rho, &parts_soa_buffer.rho[local_pack_position], count *
  //  sizeof(float)); memcpy(rho_dh,
  //  &parts_soa_buffer.rho_dh[local_pack_position], count * sizeof(float));
  //  memcpy(wcount, &parts_soa_buffer.wcount[local_pack_position], count *
  //  sizeof(float)); memcpy(wcount_dh,
  //  &parts_soa_buffer.wcount_dh[local_pack_position], count * sizeof(float));
  //  memcpy(div_v, &parts_soa_buffer.div_v[local_pack_position], count *
  //  sizeof(float)); memcpy(rot_ux,
  //  &parts_soa_buffer.rot_ux[local_pack_position], count * sizeof(float));
  //  memcpy(rot_uy, &parts_soa_buffer.rot_uy[local_pack_position], count *
  //  sizeof(float)); memcpy(rot_uz,
  //  &parts_soa_buffer.rot_uz[local_pack_position], count * sizeof(float));
  float *rho =
      &parts_soa_buffer
           .rho[local_pack_position];  // = calloc(count, sizeof(float));//
  float *rho_dh =
      &parts_soa_buffer
           .rho_dh[local_pack_position];  // = calloc(count, sizeof(float));//
  float *wcount =
      &parts_soa_buffer
           .wcount[local_pack_position];  // = calloc(count, sizeof(float));//
  float *wcount_dh =
      &parts_soa_buffer.wcount_dh[local_pack_position];  // = calloc(count,
                                                         // sizeof(float));//
  float *div_v =
      &parts_soa_buffer
           .div_v[local_pack_position];  // = calloc(count, sizeof(float));//
  float *rot_ux =
      &parts_soa_buffer
           .rot_ux[local_pack_position];  // = calloc(count, sizeof(float));//
  float *rot_uy =
      &parts_soa_buffer
           .rot_uy[local_pack_position];  // = calloc(count, sizeof(float));//
  float *rot_uz =
      &parts_soa_buffer
           .rot_uz[local_pack_position];  // = calloc(count, sizeof(float));//

  //  fprintf(stderr, "rho %f rho %f\n", rho[1],
  //  parts_soa_buffer.rho[local_pack_position+1]);
  for (int i = 0; i < count; i++) {
    //    int id_in_pack = i + local_pack_position;
    //    struct part *part_cpu = &c->hydro.parts[i];
    struct part *pi = &c->hydro.parts[i];
    //    if (part_is_inhibited(pi, e)) {
    //      fprintf(stderr, "inhibited part\n");
    //      continue;
    //    }
    //    const int pi_active = part_is_active(pi, e);
    //    if (pi_active) {
    pi->rho += rho[i];
    //      c->hydro.parts[i].rho += parts_soa_buffer.rho[id_in_pack];
    pi->density.rho_dh += rho_dh[i];
    pi->density.wcount += wcount[i];
    pi->density.wcount_dh += wcount_dh[i];
    pi->viscosity.div_v += div_v[i];
    pi->density.rot_v[0] += rot_ux[i];
    pi->density.rot_v[1] += rot_uy[i];
    pi->density.rot_v[2] += rot_uz[i];

    //      c->hydro.parts[i].rho += rho[i];
    //      c->hydro.parts[i].density.rho_dh += rho_dh[i];
    //      c->hydro.parts[i].density.wcount += wcount[i];
    //      c->hydro.parts[i].density.wcount_dh += wcount_dh[i];
    //      c->hydro.parts[i].viscosity.div_v += div_v[i];
    //      c->hydro.parts[i].density.rot_v[0] += rot_ux[i];
    //      c->hydro.parts[i].density.rot_v[1] += rot_uy[i];
    //      c->hydro.parts[i].density.rot_v[2] += rot_uz[i];

    //      c->hydro.parts[i].rho += parts_tmp->rho[i];
    //      c->hydro.parts[i].density.rho_dh += parts_tmp->rho_dh[i];
    //      c->hydro.parts[i].density.wcount += parts_tmp->wcount[i];
    //      c->hydro.parts[i].density.wcount_dh += parts_tmp->wcount_dh[i];
    //      c->hydro.parts[i].viscosity.div_v += parts_tmp->div_v[i];
    //      c->hydro.parts[i].density.rot_v[0] += parts_tmp->rot_ux[i];
    //      c->hydro.parts[i].density.rot_v[1] += parts_tmp->rot_uy[i];
    //      c->hydro.parts[i].density.rot_v[2] += parts_tmp->rot_uz[i];

    //      part_cpu[i].rho += parts_soa_buffer.rho[i];
    //      part_cpu[i].density.rho_dh += parts_soa_buffer.rho_dh[i];
    //      part_cpu[i].density.wcount += parts_soa_buffer.wcount[i];
    //      part_cpu[i].density.wcount_dh += parts_soa_buffer.wcount_dh[i];
    //      part_cpu[i].viscosity.div_v += parts_soa_buffer.div_v[i];
    //      part_cpu[i].density.rot_v[0] += parts_soa_buffer.rot_ux[i];
    //      part_cpu[i].density.rot_v[1] += parts_soa_buffer.rot_uy[i];
    //      part_cpu[i].density.rot_v[2] += parts_soa_buffer.rot_uz[i];
    //    }
    //    else fprintf(stderr,"a part is not active\n");
  }
  //  c->hydro.parts=ptmps;
}
void unpack_neat_aos(struct cell *c, struct part_aos *parts_aos_buffer, int tid,
                     int local_pack_position, int count, struct engine *e) {

  //  float *rho = &parts_aos_buffer[local_pack_position].rho;// = calloc(count,
  //  sizeof(float));// float *rho_dh  =
  //  &parts_aos_buffer[local_pack_position].rho_dh;// = calloc(count,
  //  sizeof(float));// float *wcount  =
  //  &parts_aos_buffer[local_pack_position].wcount;// = calloc(count,
  //  sizeof(float));// float *wcount_dh =
  //  &parts_aos_buffer[local_pack_position].wcount_dh;// = calloc(count,
  //  sizeof(float));// float *div_v  =
  //  &parts_aos_buffer[local_pack_position].div_v;// = calloc(count,
  //  sizeof(float));// float *rot_ux  =
  //  &parts_aos_buffer[local_pack_position].rot_ux;// = calloc(count,
  //  sizeof(float));// float *rot_uy  =
  //  &parts_aos_buffer[local_pack_position].rot_uy;// = calloc(count,
  //  sizeof(float));// float *rot_uz  =
  //  &parts_aos_buffer[local_pack_position].rot_uz;// = calloc(count,
  //  sizeof(float));//
  struct part_aos *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {

    struct part_aos p_tmp = parts_tmp[i];
    struct part *p = &c->hydro.parts[i];
    p->rho += p_tmp.rho;
    p->density.rho_dh += p_tmp.rho_dh;
    p->density.wcount += p_tmp.wcount;
    p->density.wcount_dh += p_tmp.wcount_dh;
    p->viscosity.div_v += p_tmp.div_v;
    p->density.rot_v[0] += p_tmp.rot_ux;
    p->density.rot_v[1] += p_tmp.rot_uy;
    p->density.rot_v[2] += p_tmp.rot_uz;
  }
}
#include <stdatomic.h>
void unpack_neat_aos_f4(struct cell *c,
                        struct part_aos_f4_recv *parts_aos_buffer, int tid,
                        int local_pack_position, int count, struct engine *e) {

  //  float *rho = &parts_aos_buffer[local_pack_position].rho;// = calloc(count,
  //  sizeof(float));// float *rho_dh  =
  //  &parts_aos_buffer[local_pack_position].rho_dh;// = calloc(count,
  //  sizeof(float));// float *wcount  =
  //  &parts_aos_buffer[local_pack_position].wcount;// = calloc(count,
  //  sizeof(float));// float *wcount_dh =
  //  &parts_aos_buffer[local_pack_position].wcount_dh;// = calloc(count,
  //  sizeof(float));// float *div_v  =
  //  &parts_aos_buffer[local_pack_position].div_v;// = calloc(count,
  //  sizeof(float));// float *rot_ux  =
  //  &parts_aos_buffer[local_pack_position].rot_ux;// = calloc(count,
  //  sizeof(float));// float *rot_uy  =
  //  &parts_aos_buffer[local_pack_position].rot_uy;// = calloc(count,
  //  sizeof(float));// float *rot_uz  =
  //  &parts_aos_buffer[local_pack_position].rot_uz;// = calloc(count,
  //  sizeof(float));//
  struct part_aos_f4_recv *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {

    struct part_aos_f4_recv p_tmp = parts_tmp[i];
    float4 rho_dh_wcount = p_tmp.rho_dh_wcount;
    float4 rot_ux_div_v = p_tmp.rot_ux_div_v;
    struct part *p = &c->hydro.parts[i];

    p->rho += rho_dh_wcount.x;
    p->density.rho_dh += rho_dh_wcount.y;
    p->density.wcount += rho_dh_wcount.z;
    p->density.wcount_dh += rho_dh_wcount.w;
    p->density.rot_v[0] += rot_ux_div_v.x;
    p->density.rot_v[1] += rot_ux_div_v.y;
    p->density.rot_v[2] += rot_ux_div_v.z;
    p->viscosity.div_v += rot_ux_div_v.w;
    //	      fprintf(stderr, "rho %f div_v %f\n", p_tmp.rho_dh_wcount.x,
    //p_tmp.rot_ux_div_v.w);
  }
}

void unpack_neat_aos_g(struct cell *c, struct part_aos_g *parts_aos_buffer,
                       int tid, int local_pack_position, int count,
                       struct engine *e) {

  struct part_aos_g *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {
    struct part_aos_g p_tmp = parts_tmp[i];
    struct part *p = &c->hydro.parts[i];
    const float v_sig = p->viscosity.v_sig;
    p->viscosity.v_sig = max(p_tmp.v_sig, v_sig);
    p->diffusion.laplace_u += p_tmp.laplace_u;
    const float max_ngb = p->force.alpha_visc_max_ngb;
    p->force.alpha_visc_max_ngb = max(p_tmp.alpha_visc_max_ngb, max_ngb);
  }
}

void unpack_neat_aos_f4_g(struct cell *c,
                          struct part_aos_f4_g_recv *parts_aos_buffer, int tid,
                          int local_pack_position, int count,
                          struct engine *e) {

  struct part_aos_f4_g_recv *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {
    struct part_aos_f4_g_recv p_tmp = parts_tmp[i];
    struct part *p = &c->hydro.parts[i];
    const float v_sig = p->viscosity.v_sig;
    p->viscosity.v_sig = fmaxf(p_tmp.vsig_lapu_aviscmax.x, v_sig);
    p->diffusion.laplace_u += p_tmp.vsig_lapu_aviscmax.y;
    const float max_ngb = p->force.alpha_visc_max_ngb;
    p->force.alpha_visc_max_ngb = fmaxf(p_tmp.vsig_lapu_aviscmax.z, max_ngb);
  }
}

void unpack_neat_aos_f(struct cell *c, struct part_aos_f *parts_aos_buffer,
                       int tid, int local_pack_position, int count,
                       struct engine *e) {

  struct part_aos_f *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {
    struct part_aos_f p_tmp = parts_tmp[i];
    struct part *p = &c->hydro.parts[i];
    p->a_hydro[0] += p_tmp.a_hydrox;
    p->a_hydro[1] += p_tmp.a_hydroy;
    p->a_hydro[2] += p_tmp.a_hydroz;
    p->u_dt += p_tmp.u_dt;
    p->force.h_dt += p_tmp.h_dt;
    //	      p->limiter_data.min_ngb_time_bin = min(p_tmp.min_ngb_time_bin,
    //p->limiter_data.min_ngb_time_bin);
    p->limiter_data.min_ngb_time_bin = p_tmp.min_ngb_time_bin;
    const float v_sig = p->viscosity.v_sig;
    p->viscosity.v_sig = max(p_tmp.v_sig, v_sig);
    //	      p->viscosity.v_sig = p_tmp.v_sig;

    //          fprintf(stderr, "ax %f ay %f az %f\n", p_tmp.a_hydrox,
    //          p_tmp.a_hydroy, p_tmp.a_hydroz);
  }
}

void unpack_neat_aos_f4_f(struct cell *restrict c,
                          struct part_aos_f4_f_recv *restrict parts_aos_buffer,
                          int tid, int local_pack_position, int count,
                          struct engine *e) {

  //	  struct part_aos_f4_f_recv *restrict parts_tmp =
  //&parts_aos_buffer[local_pack_position];
  int pp = local_pack_position;
  for (int i = 0; i < count; i++) {
    //	      struct part_aos_f4_f_recv p_tmp = parts_tmp[i];
    //	      struct part *restrict p = &c->hydro.parts[i];
    c->hydro.parts[i].a_hydro[0] += parts_aos_buffer[i + pp].a_hydro.x;
    c->hydro.parts[i].a_hydro[1] += parts_aos_buffer[i + pp].a_hydro.y;
    c->hydro.parts[i].a_hydro[2] += parts_aos_buffer[i + pp].a_hydro.z;
  }
  for (int i = 0; i < count; i++) {
    c->hydro.parts[i].viscosity.v_sig =
        fmaxf(parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.z,
              c->hydro.parts[i].viscosity.v_sig);
    c->hydro.parts[i].limiter_data.min_ngb_time_bin =
        (int)(parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.w + 0.5f);
  }
  for (int i = 0; i < count; i++) {
    c->hydro.parts[i].u_dt +=
        parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.x;
    c->hydro.parts[i].force.h_dt +=
        parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.y;
  }
}

void unpack_neat_pair(struct runner *r, struct cell *c,
                      struct part_soa parts_soa_buffer, int tid,
                      int local_pack_position, int count, struct engine *e) {

  //  struct part *ptmps;
  //  ptmps=c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    //    struct part *pi = &c->hydro.parts[i];
    //    if (part_is_inhibited(pi, e)) {
    //      fprintf(stderr, "inhibited part\n");
    //      continue;
    //    }
    //    const int pi_active = part_is_active(pi, e);
    //    if (pi_active) {
    c->hydro.parts[i].rho += parts_soa_buffer.rho[id_in_pack];
    c->hydro.parts[i].density.rho_dh += parts_soa_buffer.rho_dh[id_in_pack];
    c->hydro.parts[i].density.wcount += parts_soa_buffer.wcount[id_in_pack];
    c->hydro.parts[i].density.wcount_dh +=
        parts_soa_buffer.wcount_dh[id_in_pack];
    c->hydro.parts[i].viscosity.div_v += parts_soa_buffer.div_v[id_in_pack];
    c->hydro.parts[i].density.rot_v[0] += parts_soa_buffer.rot_ux[id_in_pack];
    c->hydro.parts[i].density.rot_v[1] += parts_soa_buffer.rot_uy[id_in_pack];
    c->hydro.parts[i].density.rot_v[2] += parts_soa_buffer.rot_uz[id_in_pack];
    //      if(r->cpuid == 0)fprintf(stderr, "i %i rho %lf\n", i,
    //      parts_soa_buffer.rho[id_in_pack]);
    //    }
    //    else fprintf(stderr,"a part is not active\n");
  }
  //  c->hydro.parts=ptmps;
}

void unpack_neat_pair_aos(struct runner *r, struct cell *c,
                          struct part_aos *parts_aos_buffer, int tid,
                          int local_pack_position, int count,
                          struct engine *e) {

  //  float *rho = &parts_aos_buffer[local_pack_position].rho;// = calloc(count,
  //  sizeof(float));// float *rho_dh  =
  //  &parts_aos_buffer[local_pack_position].rho_dh;// = calloc(count,
  //  sizeof(float));// float *wcount  =
  //  &parts_aos_buffer[local_pack_position].wcount;// = calloc(count,
  //  sizeof(float));// float *wcount_dh =
  //  &parts_aos_buffer[local_pack_position].wcount_dh;// = calloc(count,
  //  sizeof(float));// float *div_v  =
  //  &parts_aos_buffer[local_pack_position].div_v;// = calloc(count,
  //  sizeof(float));// float *rot_ux  =
  //  &parts_aos_buffer[local_pack_position].rot_ux;// = calloc(count,
  //  sizeof(float));// float *rot_uy  =
  //  &parts_aos_buffer[local_pack_position].rot_uy;// = calloc(count,
  //  sizeof(float));// float *rot_uz  =
  //  &parts_aos_buffer[local_pack_position].rot_uz;// = calloc(count,
  //  sizeof(float));//
  struct part_aos *parts_tmp = &parts_aos_buffer[local_pack_position];
  //  struct part *ptmps;
  //  ptmps=c->hydro.parts;
  //  struct part *part_cpu = c->hydro.parts;
  for (int i = 0; i < count; i++) {
    //    int id_in_pack = i + local_pack_position;
    //      struct part_aos part_gpu = parts_aos_buffer[id_in_pack];
    //    struct part *pi = &c->hydro.parts[i];
    //    if (part_is_inhibited(pi, e)) {
    //      fprintf(stderr, "inhibited part\n");
    //      continue;
    //    }
    //    const int pi_active = part_is_active(pi, e);
    //    if (pi_active) {
    //      if(parts_aos_buffer[id_in_pack].time_bin == 1000)(*count1000)++
    //      ;//fprintf(stderr, "timebin %i\n",
    //      parts_aos_buffer[id_in_pack].time_bin); else
    //      if(parts_aos_buffer[id_in_pack].time_bin == 20)(*count20)++
    //      ;//fprintf(stderr, "timebin %i\n",
    //      parts_aos_buffer[id_in_pack].time_bin); else fprintf(stderr, "not 20
    //      or 1000\n");
    //
    struct part_aos p_tmp = parts_tmp[i];
    struct part *p = &c->hydro.parts[i];
    p->rho += p_tmp.rho;
    p->density.rho_dh += p_tmp.rho_dh;
    p->density.wcount += p_tmp.wcount;
    p->density.wcount_dh += p_tmp.wcount_dh;
    p->viscosity.div_v += p_tmp.div_v;
    p->density.rot_v[0] += p_tmp.rot_ux;
    p->density.rot_v[1] += p_tmp.rot_uy;
    p->density.rot_v[2] += p_tmp.rot_uz;

    //      c->hydro.parts[i].rho += parts_aos_buffer[id_in_pack].rho;
    //      c->hydro.parts[i].density.rho_dh +=
    //      parts_aos_buffer[id_in_pack].rho_dh;
    //      c->hydro.parts[i].density.wcount +=
    //      parts_aos_buffer[id_in_pack].wcount;
    //      c->hydro.parts[i].density.wcount_dh +=
    //      parts_aos_buffer[id_in_pack].wcount_dh;
    //      c->hydro.parts[i].viscosity.div_v +=
    //      parts_aos_buffer[id_in_pack].div_v;
    //      c->hydro.parts[i].density.rot_v[0] +=
    //      parts_aos_buffer[id_in_pack].rot_ux;
    //      c->hydro.parts[i].density.rot_v[1] +=
    //      parts_aos_buffer[id_in_pack].rot_uy;
    //      c->hydro.parts[i].density.rot_v[2] +=
    //      parts_aos_buffer[id_in_pack].rot_uz;

    //      part_cpu[i].rho += part_gpu.rho;
    //      part_cpu[i].density.rho_dh += part_gpu.rho_dh;
    //      part_cpu[i].density.wcount += part_gpu.wcount;
    //      part_cpu[i].density.wcount_dh += part_gpu.wcount_dh;
    //      part_cpu[i].viscosity.div_v += part_gpu.div_v;
    //      part_cpu[i].density.rot_v[0] += part_gpu.rot_ux;
    //      part_cpu[i].density.rot_v[1] += part_gpu.rot_uy;
    //      part_cpu[i].density.rot_v[2] += part_gpu.rot_uz;
    //      if(r->cpuid == 0)fprintf(stderr, "i %i rho %lf\n", i,
    //      parts_soa_buffer.rho[id_in_pack]);
    //    }
    //    else fprintf(stderr,"a part is not active\n");
  }
  //  c->hydro.parts=ptmps;
}

void unpack_neat_pair_aos_f4(struct runner *r, struct cell *restrict c,
                             struct part_aos_f4_recv *restrict parts_aos_buffer,
                             int tid, int local_pack_position, int count,
                             struct engine *e) {

  //  struct part_aos_f4_recv * restrict parts_tmp =
  //  &parts_aos_buffer[local_pack_position];
  if (cell_is_active_hydro(c, e)) {
    int pp = local_pack_position;
    for (int i = 0; i < count; i++) {
      int j = i + pp;
      c->hydro.parts[i].rho += parts_aos_buffer[j].rho_dh_wcount.x;
      c->hydro.parts[i].density.rho_dh += parts_aos_buffer[j].rho_dh_wcount.y;
      c->hydro.parts[i].density.wcount += parts_aos_buffer[j].rho_dh_wcount.z;
      c->hydro.parts[i].density.wcount_dh +=
          parts_aos_buffer[j].rho_dh_wcount.w;
      c->hydro.parts[i].density.rot_v[0] += parts_aos_buffer[j].rot_ux_div_v.x;
      c->hydro.parts[i].density.rot_v[1] += parts_aos_buffer[j].rot_ux_div_v.y;
      c->hydro.parts[i].density.rot_v[2] += parts_aos_buffer[j].rot_ux_div_v.z;
      c->hydro.parts[i].viscosity.div_v += parts_aos_buffer[j].rot_ux_div_v.w;
    }
  }
}

void unpack_neat_pair_aos_g(struct runner *r, struct cell *c,
                            struct part_aos_g *parts_aos_buffer, int tid,
                            int local_pack_position, int count,
                            struct engine *e) {
  struct part_aos_g *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {
    struct part_aos_g p_tmp = parts_tmp[i];
    struct part *p = &c->hydro.parts[i];
    p->viscosity.v_sig = p_tmp.v_sig;
    p->diffusion.laplace_u += p_tmp.laplace_u;
    p->force.alpha_visc_max_ngb = p_tmp.alpha_visc_max_ngb;
  }
}

void unpack_neat_pair_aos_f4_g(
    struct runner *r, struct cell *restrict c,
    struct part_aos_f4_g_recv *restrict parts_aos_buffer, int tid,
    int local_pack_position, int count, struct engine *e) {
  //  struct part_aos_f4_recv * restrict parts_tmp =
  //  &parts_aos_buffer[local_pack_position]; int pp = local_pack_position; for
  //  (int i = 0; i < count; i++) {
  //	  int j = i + pp;
  //	  c->hydro.parts[i].viscosity.v_sig =
  //parts_aos_buffer[j].vsig_lapu_aviscmax.x;
  //	  c->hydro.parts[i].diffusion.laplace_u +=
  //parts_aos_buffer[j].vsig_lapu_aviscmax.y;
  //	  c->hydro.parts[i].force.alpha_visc_max_ngb =
  //parts_aos_buffer[j].vsig_lapu_aviscmax.z;
  //  }
  if (cell_is_active_hydro(c, e)) {

    struct part_aos_f4_g_recv *parts_tmp =
        &parts_aos_buffer[local_pack_position];
    for (int i = 0; i < count; i++) {
      struct part_aos_f4_g_recv p_tmp = parts_tmp[i];
      struct part *p = &c->hydro.parts[i];
      const float v_sig = p->viscosity.v_sig;
      p->viscosity.v_sig = fmaxf(p_tmp.vsig_lapu_aviscmax.x, v_sig);
      p->diffusion.laplace_u += p_tmp.vsig_lapu_aviscmax.y;
      const float max_ngb = p->force.alpha_visc_max_ngb;
      p->force.alpha_visc_max_ngb = fmaxf(p_tmp.vsig_lapu_aviscmax.z, max_ngb);
    }
  }
}

void unpack_neat_pair_aos_f(struct runner *r, struct cell *c,
                            struct part_aos_f *parts_aos_buffer, int tid,
                            int local_pack_position, int count,
                            struct engine *e) {
  struct part_aos_f *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {
    struct part_aos_f p_tmp = parts_tmp[i];
    struct part *p = &c->hydro.parts[i];
    p->a_hydro[0] += p_tmp.a_hydrox;
    p->a_hydro[1] += p_tmp.a_hydroy;
    p->a_hydro[2] += p_tmp.a_hydroz;
    p->u_dt += p_tmp.u_dt;
    p->force.h_dt += p_tmp.h_dt;
    const float v_sig = p->viscosity.v_sig;
    p->viscosity.v_sig = max(p_tmp.v_sig, v_sig);
    p->limiter_data.min_ngb_time_bin = p_tmp.min_ngb_time_bin;
    //	      p->viscosity.v_sig = p_tmp.v_sig;
  }
}

void unpack_neat_pair_aos_f4_f(
    struct runner *r, struct cell *restrict c,
    struct part_aos_f4_f_recv *restrict parts_aos_buffer, int tid,
    int local_pack_position, int count, struct engine *e) {
  //	  struct part_aos_f4_f_recv *restrict parts_tmp =
  //&parts_aos_buffer[local_pack_position];
  if (cell_is_active_hydro(c, e)) {
    int pp = local_pack_position;
    for (int i = 0; i < count; i++) {
      //	      struct part_aos_f4_f_recv p_tmp = parts_tmp[i];
      //	      struct part *restrict p = &c->hydro.parts[i];
      int j = i + pp;
      c->hydro.parts[i].a_hydro[0] += parts_aos_buffer[j].a_hydro.x;
      c->hydro.parts[i].a_hydro[1] += parts_aos_buffer[j].a_hydro.y;
      c->hydro.parts[i].a_hydro[2] += parts_aos_buffer[j].a_hydro.z;
      c->hydro.parts[i].viscosity.v_sig =
          fmaxf(parts_aos_buffer[j].udt_hdt_vsig_mintimebin_ngb.z,
                c->hydro.parts[i].viscosity.v_sig);
      c->hydro.parts[i].limiter_data.min_ngb_time_bin =
          (int)(parts_aos_buffer[j].udt_hdt_vsig_mintimebin_ngb.w + 0.5f);
      c->hydro.parts[i].u_dt +=
          parts_aos_buffer[j].udt_hdt_vsig_mintimebin_ngb.x;
      c->hydro.parts[i].force.h_dt +=
          parts_aos_buffer[j].udt_hdt_vsig_mintimebin_ngb.y;
    }
  }
}

void runner_do_ci_cj_gpu_unpack_neat(struct runner *r, struct cell *ci,
                                     struct cell *cj,
                                     struct part_soa parts_soa_buffer,
                                     int timer, int *pack_length, int tid,
                                     int count_max_parts_tmp,
                                     struct engine *e) {
  TIMER_TIC;

  /* Anything to do here? */
  //  if (c->hydro.count == 0)
  //    return;
  if (!cell_is_active_hydro(ci, e)) {
    message("Inactive cell\n");
    return;
  }
  int count_ci = ci->hydro.count;
  int count_cj = cj->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count_ci, e);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair(r, ci, parts_soa_buffer, tid, local_pack_position, count_ci,
                   e);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair(r, cj, parts_soa_buffer, tid, local_pack_position, count_cj,
                   e);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
  //  if(r->cpuid == 0)exit(0);
}

void runner_do_ci_cj_gpu_unpack_neat_aos(struct runner *r, struct cell *ci,
                                         struct cell *cj,
                                         struct part_aos *parts_aos_buffer,
                                         int timer, int *pack_length, int tid,
                                         int count_max_parts_tmp,
                                         struct engine *e) {

  /* Anything to do here? */
  //  if (c->hydro.count == 0)
  //    return;
  if (!cell_is_active_hydro(ci, e)) {
    message("Inactive cell\n");
    return;
  }
  int count_ci = ci->hydro.count;
  int count_cj = cj->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count_ci, e);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos(r, ci, parts_aos_buffer, tid, local_pack_position,
                       count_ci, e);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos(r, cj, parts_aos_buffer, tid, local_pack_position,
                       count_cj, e);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;
  //  if(r->cpuid == 0)exit(0);
}

void runner_do_ci_cj_gpu_unpack_neat_aos_f4(
    struct runner *r, struct cell *ci, struct cell *cj,
    struct part_aos_f4_recv *parts_aos_buffer, int timer, int *pack_length,
    int tid, int count_max_parts_tmp, struct engine *e) {

  /* Anything to do here? */
  //  if (c->hydro.count == 0)
  //    return;
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) {
    message("Inactive cell\n");
    return;
  }
  int count_ci = ci->hydro.count;
  int count_cj = cj->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count_ci, e);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f4(r, ci, parts_aos_buffer, tid, local_pack_position,
                          count_ci, e);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f4(r, cj, parts_aos_buffer, tid, local_pack_position,
                          count_cj, e);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;
  //  if(r->cpuid == 0)exit(0);
}

void runner_do_ci_cj_gpu_unpack_neat_aos_g(struct runner *r, struct cell *ci,
                                           struct cell *cj,
                                           struct part_aos_g *parts_aos_buffer,
                                           int timer, int *pack_length, int tid,
                                           int count_max_parts_tmp,
                                           struct engine *e) {

  /* Anything to do here? */
  //  if (c->hydro.count == 0)
  //    return;
  if (!cell_is_active_hydro(ci, e)) {
    message("Inactive cell\n");
    return;
  }
  int count_ci = ci->hydro.count;
  int count_cj = cj->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count_ci, e);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_g(r, ci, parts_aos_buffer, tid, local_pack_position,
                         count_ci, e);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_g(r, cj, parts_aos_buffer, tid, local_pack_position,
                         count_cj, e);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;
  //  if(r->cpuid == 0)exit(0);
}

void runner_do_ci_cj_gpu_unpack_neat_aos_f4_g(
    struct runner *r, struct cell *ci, struct cell *cj,
    struct part_aos_f4_g_recv *parts_aos_buffer, int timer, int *pack_length,
    int tid, int count_max_parts_tmp, struct engine *e) {

  /* Anything to do here? */
  //  if (c->hydro.count == 0)
  //    return;
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) {
    message("Inactive cell\n");
    return;
  }
  int count_ci = ci->hydro.count;
  int count_cj = cj->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count_ci, e);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f4_g(r, ci, parts_aos_buffer, tid, local_pack_position,
                            count_ci, e);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f4_g(r, cj, parts_aos_buffer, tid, local_pack_position,
                            count_cj, e);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;
  //  if(r->cpuid == 0)exit(0);
}

void runner_do_ci_cj_gpu_unpack_neat_aos_f(struct runner *r, struct cell *ci,
                                           struct cell *cj,
                                           struct part_aos_f *parts_aos_buffer,
                                           int timer, int *pack_length, int tid,
                                           int count_max_parts_tmp,
                                           struct engine *e) {

  /* Anything to do here? */
  //  if (c->hydro.count == 0)
  //    return;
  if (!cell_is_active_hydro(ci, e)) {
    message("Inactive cell\n");
    return;
  }
  int count_ci = ci->hydro.count;
  int count_cj = cj->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count_ci, e);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f(r, ci, parts_aos_buffer, tid, local_pack_position,
                         count_ci, e);
  local_pack_position += count_ci;
  //  for (int i = 0; i < count_ci; i++){
  //    struct part *p = &ci->hydro.parts[i];
  //    fprintf(stderr, "ax %f, ay %f, az %f, u_dt %f, h_dt %f\n",
  //    p->a_hydro[0], p->a_hydro[1], p->a_hydro[2], p->u_dt, p->force.h_dt);
  //  }
  //	      p->viscosity.v_sig = p_tmp.v_sig;
  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f(r, cj, parts_aos_buffer, tid, local_pack_position,
                         count_cj, e);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;
  //  if(r->cpuid == 0)exit(0);
}

void runner_do_ci_cj_gpu_unpack_neat_aos_f4_f(
    struct runner *r, struct cell *ci, struct cell *cj,
    struct part_aos_f4_f_recv *parts_aos_buffer, int timer, int *pack_length,
    int tid, int count_max_parts_tmp, struct engine *e) {

  /* Anything to do here? */
  //  if (c->hydro.count == 0)
  //    return;
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) {
    message("Inactive cell\n");
    return;
  }
  int count_ci = ci->hydro.count;
  int count_cj = cj->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i pointer to pack_length is %i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), pack_length, local_pack_position, count_ci, e);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f4_f(r, ci, parts_aos_buffer, tid, local_pack_position,
                            count_ci, e);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  //  if(r->cpuid == 0)fprintf(stderr, "unpacking ci l_pos %i count_i %i count_j
  //  %i\n", local_pack_position, count_ci, count_cj);
  unpack_neat_pair_aos_f4_f(r, cj, parts_aos_buffer, tid, local_pack_position,
                            count_cj, e);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;
  //  if(r->cpuid == 0)exit(0);
}

void runner_dopair_gpu_pack_neat(struct runner *r, struct cell *c,
                                 struct part_soa parts_soa_buffer, int timer,
                                 int *pack_length, int tid,
                                 int count_max_parts_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;

  int count = c->hydro.count;
  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat(c, parts_soa_buffer, tid, local_pack_position, count);
  /* Increment pack length accordingly */
  (*pack_length) += count;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_do_ci_cj_gpu_pack_neat(struct runner *r, struct cell *ci,
                                   struct cell *cj,
                                   struct part_soa parts_soa_buffer, int timer,
                                   int *pack_length, int tid,
                                   int count_max_parts_tmp, int count_ci,
                                   int count_cj) {

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;

  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= 2 * count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat(ci, parts_soa_buffer, tid, local_pack_position, count_ci);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  pack_neat(cj, parts_soa_buffer, tid, local_pack_position, count_cj);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_do_ci_cj_gpu_pack_neat_aos(struct runner *r, struct cell *ci,
                                       struct cell *cj,
                                       struct part_aos *parts_aos_buffer,
                                       int timer, int *pack_length, int tid,
                                       int count_max_parts_tmp, int count_ci,
                                       int count_cj, float3 shift_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;

  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= 2 * count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! Pack pos %i"
            "ci %i cj %i count_max %i\n",
            local_pack_position, count_ci, count_cj, count_max_parts_tmp);
    error();
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  float3 shift_i = {shift_tmp.x + cj->loc[0], shift_tmp.y + cj->loc[1],
                    shift_tmp.z + cj->loc[2]};
  pack_neat_pair_aos(ci, parts_aos_buffer, tid, local_pack_position, count_ci,
                     shift_i);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  float3 shift_j = {cj->loc[0], cj->loc[1], cj->loc[2]};
  pack_neat_pair_aos(cj, parts_aos_buffer, tid, local_pack_position, count_cj,
                     shift_j);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_do_ci_cj_gpu_pack_neat_aos_f4(
    struct runner *r, struct cell *restrict ci, struct cell *restrict cj,
    struct part_aos_f4_send *restrict parts_aos_buffer, int timer,
    int *pack_length, int tid, int count_max_parts_tmp, const int count_ci,
    const int count_cj, float3 shift_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;

  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= 2 * count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! Pack pos %i"
            "ci %i cj %i count_max %i\n",
            local_pack_position, count_ci, count_cj, count_max_parts_tmp);
    error();
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  const float3 shift_i = {shift_tmp.x + cj->loc[0], shift_tmp.y + cj->loc[1],
                          shift_tmp.z + cj->loc[2]};
  const int lpp1 = local_pack_position;

  const int2 cis_cie = {local_pack_position, local_pack_position + count_ci};

  const int2 cjs_cje = {local_pack_position + count_ci,
                        local_pack_position + count_ci + count_cj};

  pack_neat_pair_aos_f4(ci, parts_aos_buffer, tid, lpp1, count_ci, shift_i,
                        cjs_cje);

  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  const float3 shift_j = {cj->loc[0], cj->loc[1], cj->loc[2]};
  const int lpp2 = local_pack_position;

  pack_neat_pair_aos_f4(cj, parts_aos_buffer, tid, lpp2, count_cj, shift_j,
                        cis_cie);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_do_ci_cj_gpu_pack_neat_aos_g(struct runner *r, struct cell *ci,
                                         struct cell *cj,
                                         struct part_aos_g *parts_aos_buffer,
                                         int timer, int *pack_length, int tid,
                                         int count_max_parts_tmp, int count_ci,
                                         int count_cj) {

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;

  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= 2 * count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! Pack pos %i"
            "ci %i cj %i count_max %i\n",
            local_pack_position, count_ci, count_cj, count_max_parts_tmp);
    error();
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_g(ci, parts_aos_buffer, tid, local_pack_position, count_ci);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_g(cj, parts_aos_buffer, tid, local_pack_position, count_cj);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_do_ci_cj_gpu_pack_neat_aos_f4_g(
    struct runner *r, struct cell *restrict ci, struct cell *restrict cj,
    struct part_aos_f4_g_send *restrict parts_aos_buffer, int timer,
    int *pack_length, int tid, int count_max_parts_tmp, const int count_ci,
    const int count_cj, float3 shift_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;

  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= 2 * count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! Pack pos %i"
            "ci %i cj %i count_max %i\n",
            local_pack_position, count_ci, count_cj, count_max_parts_tmp);
    error();
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  const float3 shift_i = {shift_tmp.x + cj->loc[0], shift_tmp.y + cj->loc[1],
                          shift_tmp.z + cj->loc[2]};
  const int lpp1 = local_pack_position;

  const int2 cis_cie = {local_pack_position, local_pack_position + count_ci};

  const int2 cjs_cje = {local_pack_position + count_ci,
                        local_pack_position + count_ci + count_cj};

  pack_neat_pair_aos_f4_g(ci, parts_aos_buffer, tid, lpp1, count_ci, shift_i,
                          cjs_cje);

  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  const float3 shift_j = {cj->loc[0], cj->loc[1], cj->loc[2]};
  const int lpp2 = local_pack_position;

  pack_neat_pair_aos_f4_g(cj, parts_aos_buffer, tid, lpp2, count_cj, shift_j,
                          cis_cie);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_do_ci_cj_gpu_pack_neat_aos_f(struct runner *r, struct cell *ci,
                                         struct cell *cj,
                                         struct part_aos_f *parts_aos_buffer,
                                         int timer, int *pack_length, int tid,
                                         int count_max_parts_tmp, int count_ci,
                                         int count_cj) {

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;

  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= 2 * count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! Pack pos %i"
            "ci %i cj %i count_max %i\n",
            local_pack_position, count_ci, count_cj, count_max_parts_tmp);
    error();
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_f(ci, parts_aos_buffer, tid, local_pack_position, count_ci);
  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  pack_neat_aos_f(cj, parts_aos_buffer, tid, local_pack_position, count_cj);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_do_ci_cj_gpu_pack_neat_aos_f4_f(
    struct runner *r, struct cell *restrict ci, struct cell *restrict cj,
    struct part_aos_f4_f_send *restrict parts_aos_buffer, int timer,
    int *pack_length, int tid, int count_max_parts_tmp, const int count_ci,
    const int count_cj, float3 shift_tmp) {

  TIMER_TIC;

  /* Anything to do here? */
  if (ci->hydro.count == 0) return;

  int local_pack_position = (*pack_length);

#ifdef SWIFT_DEBUG_CHECKS
  if (local_pack_position + count_ci + count_cj >= 2 * count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! Pack pos %i"
            "ci %i cj %i count_max %i\n",
            local_pack_position, count_ci, count_cj, count_max_parts_tmp);
    error();
  }
#endif

  /* Pack the particle data into CPU-side buffers*/
  const float3 shift_i = {shift_tmp.x + cj->loc[0], shift_tmp.y + cj->loc[1],
                          shift_tmp.z + cj->loc[2]};
  const int lpp1 = local_pack_position;

  const int2 cis_cie = {local_pack_position, local_pack_position + count_ci};

  const int2 cjs_cje = {local_pack_position + count_ci,
                        local_pack_position + count_ci + count_cj};

  pack_neat_pair_aos_f4_f(ci, parts_aos_buffer, tid, lpp1, count_ci, shift_i,
                          cjs_cje);

  local_pack_position += count_ci;
  /* Pack the particle data into CPU-side buffers*/
  const float3 shift_j = {cj->loc[0], cj->loc[1], cj->loc[2]};
  const int lpp2 = local_pack_position;

  pack_neat_pair_aos_f4_f(cj, parts_aos_buffer, tid, lpp2, count_cj, shift_j,
                          cis_cie);
  /* Increment pack length accordingly */
  (*pack_length) += count_ci + count_cj;

  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void runner_doself1_gpu_pack(
    struct runner *r, struct cell *c, int timer, int *pack_length, double *x_p,
    double *y_p, double *z_p, int tid, int *tid_p, long long *id, float *ux,
    float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
    float *mass, float *h, float *u, float *u_dt, float *rho, float *SPH_sum,
    float *locx, float *locy, float *locz, float *widthx, float *widthy,
    float *widthz, float *h_max, int *count_p, float *wcount, float *wcount_dh,
    float *rho_dh, float *rot_u, float *rot_v, float *rot_w, float *div_v,
    float *div_v_previous_step, float *alpha_visc, float *v_sig,
    float *laplace_u, float *alpha_diff, float *f, float *soundspeed,
    float *h_dt, float *balsara, float *pressure, float *alpha_visc_max_ngb,
    timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
    char *to_be_synchronized, int count_max_parts_tmp) {

  TIMER_TIC;
  //  fprintf(stderr,"Entered outer packing code!\n");

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  /* Recurse? */
  //  if (c->split) {
  ////	fprintf(stderr,"Entered recursive packing code!\n");
  //    for (int k = 0; k < 8; k++){
  //      if (c->progeny[k] != NULL){
  //    	  runner_doself1_gpu_pack(r, c, timer, pack_length,
  //    	  x_p, y_p, z_p, tid, tid_p, id, ux, uy, uz, a_hydrox,
  //    	  a_hydroy, a_hydroz, mass, h, u, u_dt, rho, SPH_sum, locx,
  //    locy, locz, 	  widthx, widthy, widthz, h_max, count_p, wcount,
  //    wcount_dh, rho_dh, rot_u, rot_v, 	  rot_w, div_v,
  //    div_v_previous_step, alpha_visc, v_sig, laplace_u, alpha_diff, f,
  //    soundspeed, 	  h_dt, balsara, pressure, alpha_visc_max_ngb, time_bin,
  //    wakeup, min_ngb_time_bin, 	  to_be_synchronized, count_max_parts_tmp,
  //    fgpuin); 	  fprintf(stderr,"working on a split cell\n");
  //      }
  //    }
  //  }
  //  else {
  //	    fprintf(stderr,"Entered inner packing code!\n");
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr, "Exceeded count_max_parts_tmp. Make arrays bigger!\n");
    exit(0);
  }
  // Pack the particle data
  pack(c, x_p, y_p, z_p, tid, tid_p, id, ux, uy, uz, a_hydrox, a_hydroy,
       a_hydroz, mass, h, u, u_dt, rho, SPH_sum, locx, locy, locz, widthx,
       widthy, widthz, h_max, count_p, wcount, wcount_dh, rho_dh, rot_u, rot_v,
       rot_w, div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
       alpha_diff, f, soundspeed, h_dt, balsara, pressure, alpha_visc_max_ngb,
       time_bin, wakeup, min_ngb_time_bin, to_be_synchronized,
       local_pack_position, count);
  // Increment pack length accordingly
  (*pack_length) += count;
  //  }
  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void pack(struct cell *c, double *x_p, double *y_p, double *z_p, int tid,
          int *tid_p, long long *id, float *ux, float *uy, float *uz,
          float *a_hydrox, float *a_hydroy, float *a_hydroz, float *mass,
          float *h, float *u, float *u_dt, float *rho, float *SPH_sum,
          float *locx, float *locy, float *locz, float *widthx, float *widthy,
          float *widthz, float *h_max, int *count_p, float *wcount,
          float *wcount_dh, float *rho_dh, float *rot_u, float *rot_v,
          float *rot_w, float *div_v, float *div_v_previous_step,
          float *alpha_visc, float *v_sig, float *laplace_u, float *alpha_diff,
          float *f, float *soundspeed, float *h_dt, float *balsara,
          float *pressure, float *alpha_visc_max_ngb, timebin_t *time_bin,
          timebin_t *wakeup, timebin_t *min_ngb_time_bin,
          char *to_be_synchronized, int local_pack_position, int count) {

  const struct part *ptmps;
  ptmps = c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    x_p[id_in_pack] = ptmps[i].x[0];
    y_p[id_in_pack] = ptmps[i].x[1];
    z_p[id_in_pack] = ptmps[i].x[2];
    //    id[id_in_pack]=ptmps[i].id;
    //    count_p[id_in_pack]=count;
    tid_p[id_in_pack] = tid;
    //    h_max[id_in_pack]=c->hydro.h_max;
    ux[id_in_pack] = ptmps[i].v[0];
    uy[id_in_pack] = ptmps[i].v[1];
    uz[id_in_pack] = ptmps[i].v[2];
    //	a_hydrox[id_in_pack]=ptmps[i].a_hydro[0];
    //	a_hydroy[id_in_pack]=ptmps[i].a_hydro[1];
    //	a_hydroz[id_in_pack]=ptmps[i].a_hydro[2];
    locx[id_in_pack] = c->loc[0];
    locy[id_in_pack] = c->loc[1];
    locz[id_in_pack] = c->loc[2];
    mass[id_in_pack] = ptmps[i].mass;
    h[id_in_pack] = ptmps[i].h;
    //	u[id_in_pack]=ptmps[i].u;
    //	u_dt[id_in_pack]=ptmps[i].u_dt;
    //////////////////////////////////////////////////////
    rho[id_in_pack] = 0.f;  // ptmps[i].rho;
    /////////////////////////////////////////////////////
    //	div_v_previous_step[id_in_pack]=ptmps[i].viscosity.div_v_previous_step;
    //	alpha_visc[id_in_pack]=ptmps[i].viscosity.alpha;
    //	v_sig[id_in_pack]=ptmps[i].viscosity.v_sig;
    //	laplace_u[id_in_pack]=ptmps[i].diffusion.laplace_u;
    //	alpha_diff[id_in_pack]=ptmps[i].diffusion.alpha;
    //	f[id_in_pack]=ptmps[i].force.f;
    //	soundspeed[id_in_pack]=ptmps[i].force.soundspeed;
    //	h_dt[id_in_pack]=ptmps[i].force.h_dt;
    //	balsara[id_in_pack]=ptmps[i].force.balsara;
    //	pressure[id_in_pack]=ptmps[i].force.pressure;
    //    time_bin[id_in_pack] = ptmps[i].time_bin;
    //	wakeup[id_in_pack]=ptmps[i].limiter_data.wakeup;
    //	min_ngb_time_bin[id_in_pack]=ptmps[i].limiter_data.min_ngb_time_bin;
    //	to_be_synchronized[id_in_pack]=ptmps[i].limiter_data.to_be_synchronized;
    ///////////////////////////////////////////////////////////////////
    wcount[id_in_pack] = 0.f;     // ptmps[i].density.wcount;
    wcount_dh[id_in_pack] = 0.f;  // ptmps[i].density.wcount_dh;
    rho_dh[id_in_pack] = 0.f;     // ptmps[i].density.rho_dh;
    div_v[id_in_pack] = 0.f;      // ptmps[i].viscosity.div_v;
    rot_u[id_in_pack] = 0.f;      // ptmps[i].density.rot_v[0];
    rot_v[id_in_pack] = 0.f;      // ptmps[i].density.rot_v[1];
    rot_w[id_in_pack] = 0.f;      // ptmps[i].density.rot_v[2];
    ///////////////////////////////////////////////////////////////////
  }
}

void runner_doself1_gpu_unpack(
    struct runner *r, struct cell *c, int timer, int *pack_length, double *x_p,
    double *y_p, double *z_p, int tid, int *tid_p, long long *id, float *ux,
    float *uy, float *uz, float *a_hydrox, float *a_hydroy, float *a_hydroz,
    float *mass, float *h, float *u, float *u_dt, float *rho, float *SPH_sum,
    float *locx, float *locy, float *locz, float *widthx, float *widthy,
    float *widthz, float *h_max, int *count_p, float *wcount, float *wcount_dh,
    float *rho_dh, float *rot_u, float *rot_v, float *rot_w, float *div_v,
    float *div_v_previous_step, float *alpha_visc, float *v_sig,
    float *laplace_u, float *alpha_diff, float *f, float *soundspeed,
    float *h_dt, float *balsara, float *pressure, float *alpha_visc_max_ngb,
    timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
    char *to_be_synchronized, int count_max_parts_tmp, struct engine *e) {
  TIMER_TIC;
  //  fprintf(stderr, "got into pack function\n");
  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) return;
  /* Anything to do here? */
  /* Recurse? */
  //  if (c->split) {
  //	  fprintf(stderr,"working on a split cell\n");
  //    for (int k = 0; k < 8; k++){
  //      if (c->progeny[k] != NULL){
  //    	  runner_doself1_gpu_unpack(r, c, timer, pack_length,
  //    	  x_p, y_p, z_p, tid, tid_p, id, ux, uy, uz, a_hydrox,
  //    	  a_hydroy, a_hydroz, mass, h, u, u_dt, rho, SPH_sum, locx,
  //    locy, locz, 	  widthx, widthy, widthz, h_max, count_p, wcount,
  //    wcount_dh, rho_dh, rot_u, rot_v, 	  rot_w, div_v,
  //    div_v_previous_step, alpha_visc, v_sig, laplace_u, alpha_diff, f,
  //    soundspeed, 	  h_dt, balsara, pressure, alpha_visc_max_ngb, time_bin,
  //    wakeup, min_ngb_time_bin, 	  to_be_synchronized, count_max_parts_tmp,
  //    fgpuin); 	  fprintf(stderr,"working on a split cell\n");
  //      }
  //    }
  //  } else {
  int count = c->hydro.count;
  int local_pack_position = (*pack_length);
  if (local_pack_position + count >= count_max_parts_tmp) {
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! pack_length is "
            "%i, local_pack_position is % i, "
            "count is %i\n",
            (*pack_length), local_pack_position, count);
    //	      exit(0);
  }
  // Pack the particle data
  unpack(c, x_p, y_p, z_p, tid, tid_p, id, ux, uy, uz, a_hydrox, a_hydroy,
         a_hydroz, mass, h, u, u_dt, rho, SPH_sum, locx, locy, locz, widthx,
         widthy, widthz, h_max, count_p, wcount, wcount_dh, rho_dh, rot_u,
         rot_v, rot_w, div_v, div_v_previous_step, alpha_visc, v_sig, laplace_u,
         alpha_diff, f, soundspeed, h_dt, balsara, pressure, alpha_visc_max_ngb,
         time_bin, wakeup, min_ngb_time_bin, to_be_synchronized,
         local_pack_position, count, e);
  //  for (int i = *pack_length; i < count+*pack_length; i++) {
  //  for (int i = 0; i < count; i++) {
  //	  message("wcount is %f", c->hydro.parts[i].density.wcount);
  //  }
  // Increment pack length accordingly
  (*pack_length) += count;
  //  }
  if (timer) TIMER_TOC(timer_doself_gpu_pack);
}

void unpack(struct cell *c, double *x_p, double *y_p, double *z_p, int tid,
            int *tid_p, long long *id, float *ux, float *uy, float *uz,
            float *a_hydrox, float *a_hydroy, float *a_hydroz, float *mass,
            float *h, float *u, float *u_dt, float *rho, float *SPH_sum,
            float *locx, float *locy, float *locz, float *widthx, float *widthy,
            float *widthz, float *h_max, int *count_p, float *wcount,
            float *wcount_dh, float *rho_dh, float *rot_u, float *rot_v,
            float *rot_w, float *div_v, float *div_v_previous_step,
            float *alpha_visc, float *v_sig, float *laplace_u,
            float *alpha_diff, float *f, float *soundspeed, float *h_dt,
            float *balsara, float *pressure, float *alpha_visc_max_ngb,
            timebin_t *time_bin, timebin_t *wakeup, timebin_t *min_ngb_time_bin,
            char *to_be_synchronized, int local_pack_position, int count,
            struct engine *e) {

  //  struct part *ptmps;
  //  ptmps=c->hydro.parts;
  for (int i = 0; i < count; i++) {
    int id_in_pack = i + local_pack_position;
    struct part *pi = &c->hydro.parts[i];
    if (part_is_inhibited(pi, e)) {
      fprintf(stderr, "inhibited part\n");
      continue;
    }
    const int pi_active = part_is_active(pi, e);
    if (!pi_active)
      fprintf(stderr, "Inactive part\n");
    else if (pi_active) {
      //    c->hydro.parts[i].rho = rho[id_in_pack];
      //    c->hydro.parts[i].viscosity.div_v = div_v[id_in_pack];
      //    c->hydro.parts[i].density.rho_dh = rho_dh[id_in_pack];
      //    c->hydro.parts[i].density.wcount = wcount[id_in_pack];
      //    c->hydro.parts[i].density.wcount_dh = wcount_dh[id_in_pack];
      //    c->hydro.parts[i].density.rot_v[0] = rot_u[id_in_pack];
      //    c->hydro.parts[i].density.rot_v[1] = rot_v[id_in_pack];
      //    c->hydro.parts[i].density.rot_v[2] = rot_w[id_in_pack];
      pi->rho += rho[id_in_pack];
      pi->viscosity.div_v += div_v[id_in_pack];
      pi->density.rho_dh += rho_dh[id_in_pack];
      pi->density.wcount += wcount[id_in_pack];
      pi->density.wcount_dh += wcount_dh[id_in_pack];
      pi->density.rot_v[0] += rot_u[id_in_pack];
      pi->density.rot_v[1] += rot_v[id_in_pack];
      pi->density.rot_v[2] += rot_w[id_in_pack];
    }
    //    else fprintf(stderr,"a part is not active\n");
  }
  //  c->hydro.parts=ptmps;
}
// #ifdef WITHCUDA
// }
// #endif
