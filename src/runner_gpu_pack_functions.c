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
#include "runner_doiact_hydro.h"

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
    fprintf(stderr,
            "Exceeded count_max_parts_tmp. Make arrays bigger! count_max %i "
            "count %i\n",
            count_max_parts_tmp, local_pack_position + count);
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

#include <stdatomic.h>
void unpack_neat_aos_f4(struct cell *c,
                        struct part_aos_f4_recv *parts_aos_buffer, int tid,
                        int local_pack_position, int count, struct engine *e) {

  struct part_aos_f4_recv *parts_tmp = &parts_aos_buffer[local_pack_position];
  for (int i = 0; i < count; i++) {

    struct part_aos_f4_recv p_tmp = parts_tmp[i];
    float4 rho_dh_wcount = p_tmp.rho_dh_wcount;
    float4 rot_ux_div_v = p_tmp.rot_ux_div_v;
    struct part *p = &c->hydro.parts[i];
    if(!PART_IS_ACTIVE(p, e))continue;
    p->rho += rho_dh_wcount.x;
    p->density.rho_dh += rho_dh_wcount.y;
    p->density.wcount += rho_dh_wcount.z;
    p->density.wcount_dh += rho_dh_wcount.w;
    p->density.rot_v[0] += rot_ux_div_v.x;
    p->density.rot_v[1] += rot_ux_div_v.y;
    p->density.rot_v[2] += rot_ux_div_v.z;
    p->viscosity.div_v += rot_ux_div_v.w;
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
    if(!PART_IS_ACTIVE(p, e))continue;
    const float v_sig = p->viscosity.v_sig;
    p->viscosity.v_sig = fmaxf(p_tmp.vsig_lapu_aviscmax.x, v_sig);
    p->diffusion.laplace_u += p_tmp.vsig_lapu_aviscmax.y;
    const float max_ngb = p->force.alpha_visc_max_ngb;
    p->force.alpha_visc_max_ngb = fmaxf(p_tmp.vsig_lapu_aviscmax.z, max_ngb);
  }
}

void unpack_neat_aos_f4_f(struct cell *restrict c,
                          struct part_aos_f4_f_recv *restrict parts_aos_buffer,
                          int tid, int local_pack_position, int count,
                          struct engine *e) {
  int pp = local_pack_position;
  for (int i = 0; i < count; i++) {
	if(!PART_IS_ACTIVE(&c->hydro.parts[i], e))continue;
    c->hydro.parts[i].a_hydro[0] += parts_aos_buffer[i + pp].a_hydro.x;
    c->hydro.parts[i].a_hydro[1] += parts_aos_buffer[i + pp].a_hydro.y;
    c->hydro.parts[i].a_hydro[2] += parts_aos_buffer[i + pp].a_hydro.z;
  }
  for (int i = 0; i < count; i++) {
	if(!PART_IS_ACTIVE(&c->hydro.parts[i], e))continue;
    c->hydro.parts[i].viscosity.v_sig =
        fmaxf(parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.z,
              c->hydro.parts[i].viscosity.v_sig);
    c->hydro.parts[i].limiter_data.min_ngb_time_bin =
        (int)(parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.w + 0.5f);
  }
  for (int i = 0; i < count; i++) {
    if(!PART_IS_ACTIVE(&c->hydro.parts[i], e))continue;
    c->hydro.parts[i].u_dt +=
        parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.x;
    c->hydro.parts[i].force.h_dt +=
        parts_aos_buffer[i + pp].udt_hdt_vsig_mintimebin_ngb.y;
  }
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

void unpack_neat_pair_aos_f4_g(
    struct runner *r, struct cell *restrict c,
    struct part_aos_f4_g_recv *restrict parts_aos_buffer, int tid,
    int local_pack_position, int count, struct engine *e) {
  //  struct part_aos_f4_recv * restrict parts_tmp =
  //  &parts_aos_buffer[local_pack_position]; int pp = local_pack_position; for
  //  (int i = 0; i < count; i++) {
  //	  int j = i + pp;
  //	  c->hydro.parts[i].viscosity.v_sig =
  // parts_aos_buffer[j].vsig_lapu_aviscmax.x;
  //	  c->hydro.parts[i].diffusion.laplace_u +=
  // parts_aos_buffer[j].vsig_lapu_aviscmax.y;
  //	  c->hydro.parts[i].force.alpha_visc_max_ngb =
  // parts_aos_buffer[j].vsig_lapu_aviscmax.z;
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

void runner_do_ci_cj_gpu_unpack_neat_aos_f4(
    struct runner *r, struct cell *ci, struct cell *cj,
    struct part_aos_f4_recv *parts_aos_buffer, int timer, int *pack_length,
    int tid, int count_max_parts_tmp, struct engine *e) {

  /* Anything to do here? */
//    if (ci->hydro.count == 0 || cj->hydro.count == 0)
//      return;
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
// #ifdef WITHCUDA
// }
// #endif
