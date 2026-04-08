/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_ACTIVE_H
#define SWIFT_ACTIVE_H

/* Config parameters. */
#include <config.h>

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "feedback_new_stars.h"
#include "part.h"
#include "proxy.h"
#include "task.h"
#include "timeline.h"

#ifdef SWIFT_DEBUG_CHECKS
__attribute__((always_inline)) INLINE static void
cell_active_gravity_debug_error(const struct cell *c, const struct engine *e) {

#ifdef WITH_MPI
  int proxy_in = 0, proxy_out = 0;
  int proxy_in_grav = 0, proxy_out_grav = 0;
  int top_proxy_in = 0, top_proxy_out = 0;
  int top_proxy_in_grav = 0, top_proxy_out_grav = 0;
  int first_proxy_in_node = -1, first_proxy_in_type = 0;
  int first_proxy_out_node = -1, first_proxy_out_type = 0;
#endif
  int task_ref = 0, task_ref_active = 0;
  int pair_ref = 0, pair_ref_active = 0;
  int mm_ref = 0, mm_ref_active = 0;
  int tend_ref = 0, tend_ref_active = 0;
  int top_pair_ref = 0, top_mm_ref = 0;
  int void_mm_ref = 0;
  int first_task_type = task_type_none;
  int first_task_subtype = task_subtype_none;
  int first_task_skip = -1;
  long long first_task_flags = -1;

#ifdef WITH_MPI
  for (int pid = 0; pid < e->nr_proxies; pid++) {
    const struct proxy *p = &e->proxies[pid];
    for (int k = 0; k < p->nr_cells_in; k++) {
      if (p->cells_in[k] == c) {
        proxy_in += 1;
        if (p->cells_in_type[k] & proxy_cell_type_gravity) proxy_in_grav += 1;
        if (first_proxy_in_node < 0) {
          first_proxy_in_node = p->nodeID;
          first_proxy_in_type = p->cells_in_type[k];
        }
      }
      if (c->top != NULL && p->cells_in[k] == c->top) {
        top_proxy_in += 1;
        if (p->cells_in_type[k] & proxy_cell_type_gravity)
          top_proxy_in_grav += 1;
      }
    }
    for (int k = 0; k < p->nr_cells_out; k++) {
      if (p->cells_out[k] == c) {
        proxy_out += 1;
        if (p->cells_out_type[k] & proxy_cell_type_gravity) proxy_out_grav += 1;
        if (first_proxy_out_node < 0) {
          first_proxy_out_node = p->nodeID;
          first_proxy_out_type = p->cells_out_type[k];
        }
      }
      if (c->top != NULL && p->cells_out[k] == c->top) {
        top_proxy_out += 1;
        if (p->cells_out_type[k] & proxy_cell_type_gravity)
          top_proxy_out_grav += 1;
      }
    }
  }
#endif

  for (int i = 0; i < e->sched.nr_tasks; i++) {
    const struct task *t = &e->sched.tasks[i];
    const int touches_c = (t->ci == c || t->cj == c);
    const int touches_top =
        (c->top != NULL && (t->ci == c->top || t->cj == c->top));
    const int touches_void_parent =
        (c->void_parent != NULL &&
         (t->ci == c->void_parent || t->cj == c->void_parent));

    if (touches_c) {
      task_ref += 1;
      if (!t->skip) task_ref_active += 1;
      if (first_task_skip < 0) {
        first_task_type = t->type;
        first_task_subtype = t->subtype;
        first_task_skip = t->skip;
        first_task_flags = t->flags;
      }
      if (t->type == task_type_pair && t->subtype == task_subtype_grav) {
        pair_ref += 1;
        if (!t->skip) pair_ref_active += 1;
      } else if (t->type == task_type_grav_mm) {
        mm_ref += 1;
        if (!t->skip) mm_ref_active += 1;
      } else if ((t->type == task_type_send || t->type == task_type_recv) &&
                 t->subtype == task_subtype_tend) {
        tend_ref += 1;
        if (!t->skip) tend_ref_active += 1;
      }
    }

    if (touches_top) {
      if (t->type == task_type_pair && t->subtype == task_subtype_grav)
        top_pair_ref += 1;
      else if (t->type == task_type_grav_mm)
        top_mm_ref += 1;
    }

    if (touches_void_parent && t->type == task_type_grav_mm) void_mm_ref += 1;
  }

  error(
      "cell in an impossible time-zone! c->ti_end_min=%lld (t=%e) and "
      "e->ti_current=%lld (t=%e, a=%e) cell: %s/%s depth=%d nodeID=%d "
#ifdef WITH_MPI
      "recv=%p send=%p tag=%d proxy_in=%d/%d proxy_out=%d/%d "
      "top_proxy_in=%d/%d "
      "top_proxy_out=%d/%d first_in=(node=%d,type=%d) "
      "first_out=(node=%d,type=%d) "
#endif
      "top=%p super=%p grav.super=%p void_parent=%p contains_zoom=%d "
      "below_diff=%d "
      "task_ref=%d/%d pair_ref=%d/%d mm_ref=%d/%d tend_ref=%d/%d top_pair=%d "
      "top_mm=%d "
      "void_parent_mm=%d first_task=%s/%s skip=%d flags=%lld timestep=%p "
      "collect=%p "
      "grav.drift=%p grav.grav=%p grav.mm=%p grav.init=%p grav.long_range=%p "
      "grav.down=%p "
      "grav.count=%d ti_old_part=%lld ti_old_multipole=%lld m_pole=%g",
      c->grav.ti_end_min, c->grav.ti_end_min * e->time_base, e->ti_current,
      e->ti_current * e->time_base, e->cosmology->a, cellID_names[c->type],
      subcellID_names[c->subtype], c->depth, c->nodeID,
#ifdef WITH_MPI
      c->mpi.recv, c->mpi.send, c->mpi.tag, proxy_in, proxy_in_grav, proxy_out,
      proxy_out_grav, top_proxy_in, top_proxy_in_grav, top_proxy_out,
      top_proxy_out_grav, first_proxy_in_node, first_proxy_in_type,
      first_proxy_out_node, first_proxy_out_type,
#endif
      c->top, c->super, c->grav.super, c->void_parent, c->contains_zoom_cells,
      c->grav.tasks_below_diff_grav_depth, task_ref, task_ref_active, pair_ref,
      pair_ref_active, mm_ref, mm_ref_active, tend_ref, tend_ref_active,
      top_pair_ref, top_mm_ref, void_mm_ref, taskID_names[first_task_type],
      subtaskID_names[first_task_subtype], first_task_skip, first_task_flags,
      c->timestep, c->timestep_collect, c->grav.drift, c->grav.grav, c->grav.mm,
      c->grav.init, c->grav.long_range, c->grav.down, c->grav.count,
      c->grav.ti_old_part, c->grav.ti_old_multipole,
      c->grav.multipole != NULL ? c->grav.multipole->m_pole.M_000 : -1.);
}
#endif

/**
 * @brief Check that the #part in a #cell have been drifted to the current time.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell has been drifted to the current time, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_are_part_drifted(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.ti_old_part > e->ti_current)
    error(
        "Cell has been drifted too far forward in time! c->ti_old_part=%lld "
        "(t=%e) "
        "and e->ti_current=%lld (t=%e, a=%e)",
        c->hydro.ti_old_part, c->hydro.ti_old_part * e->time_base,
        e->ti_current, e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->hydro.ti_old_part == e->ti_current);
}

/**
 * @brief Check that the #gpart in a #cell have been drifted to the current
 * time.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell has been drifted to the current time, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_are_gpart_drifted(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grav.ti_old_part > e->ti_current)
    error(
        "Cell has been drifted too far forward in time! c->ti_old=%lld (t=%e) "
        "and e->ti_current=%lld (t=%e)",
        c->grav.ti_old_part, c->grav.ti_old_part * e->time_base, e->ti_current,
        e->ti_current * e->time_base);
#endif

  return (c->grav.ti_old_part == e->ti_current);
}

/**
 * @brief Check that the #spart in a #cell have been drifted to the current
 * time.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell has been drifted to the current time, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_are_spart_drifted(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->stars.ti_old_part > e->ti_current)
    error(
        "Cell has been drifted too far forward in time! c->ti_old=%lld (t=%e) "
        "and e->ti_current=%lld (t=%e)",
        c->stars.ti_old_part, c->stars.ti_old_part * e->time_base,
        e->ti_current, e->ti_current * e->time_base);
#endif

  return (c->stars.ti_old_part == e->ti_current);
}

/**
 * @brief Check that the #sink in a #cell have been drifted to the current
 * time.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell has been drifted to the current time, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_are_sink_drifted(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->sinks.ti_old_part > e->ti_current)
    error(
        "Cell has been drifted too far forward in time! c->ti_old=%lld (t=%e) "
        "and e->ti_current=%lld (t=%e)",
        c->sinks.ti_old_part, c->sinks.ti_old_part * e->time_base,
        e->ti_current, e->ti_current * e->time_base);
#endif

  return (c->sinks.ti_old_part == e->ti_current);
}

/**
 * @brief Check that the #bpart in a #cell have been drifted to the current
 * time.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell has been drifted to the current time, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_are_bpart_drifted(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->black_holes.ti_old_part > e->ti_current)
    error(
        "Cell has been drifted too far forward in time! c->ti_old=%lld (t=%e) "
        "and e->ti_current=%lld (t=%e)",
        c->black_holes.ti_old_part, c->black_holes.ti_old_part * e->time_base,
        e->ti_current, e->ti_current * e->time_base);
#endif

  return (c->black_holes.ti_old_part == e->ti_current);
}

/**
 * @brief Check that the #part in a #cell have been drifted to the current time.
 * This is just a prototype function to keep the iact functions clean. As we
 * don't care about the drifts during the RT sub-cycling, this always just
 * returns true.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell has been drifted to the current time, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int
cell_are_part_drifted_rt_sub_cycle(const struct cell *c,
                                   const struct engine *e) {

  /* Note: we can't just use "cell_are_part_drifted" in the hydro_iact
   * functions, because an RT sub-cycle may be called during a main
   * step for a cell that is hydro inactive and thus may be not drifted. */
  return 1;
}

/* Are cells / particles active for regular tasks ? */

/**
 * @brief Does a cell contain any particle finishing their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_active_hydro(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.ti_end_min < e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_end_min=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e) c->nodeID=%d",
        c->hydro.ti_end_min, c->hydro.ti_end_min * e->time_base, e->ti_current,
        e->ti_current * e->time_base, e->cosmology->a, c->nodeID);
#endif

  return (c->hydro.ti_end_min == e->ti_current);
}

/**
 * @brief Does a cell contain any particle finishing their RT time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_rt_active(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->rt.ti_rt_end_min < e->ti_current_subcycle) {
    error(
        "cell %lld in an impossible time-zone! c->ti_rt_end_min=%lld (t=%e) "
        "and e->ti_current=%lld (t=%e, a=%e) c->nodeID=%d ACT=%d count=%d",
        c->cellID, c->rt.ti_rt_end_min, c->rt.ti_rt_end_min * e->time_base,
        e->ti_current_subcycle, e->ti_current_subcycle * e->time_base,
        e->cosmology->a, c->nodeID, c->rt.rt_advance_cell_time != NULL,
        c->hydro.count);
  }
#endif

  /* If there are no sub-cycles, e->ti_current_subcycle = e->ti_current.
   * This is also the case if we're currently doing a normal SWIFT step. */
  return (c->rt.ti_rt_end_min == e->ti_current_subcycle);
}

/**
 * @brief Does a cell contain any g-particle finishing their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_active_gravity(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grav.ti_end_min < e->ti_current) {
    cell_active_gravity_debug_error(c, e);
  }
#endif

  return (c->grav.ti_end_min == e->ti_current);
}

/**
 * @brief Does a cell contain any multipole requiring calculation ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_active_gravity_mm(
    const struct cell *c, const struct engine *e) {

  return (c->grav.ti_end_min == e->ti_current);
}

/**
 * @brief Does a cell contain any s-particle finishing their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_active_stars(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->stars.ti_end_min < e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_end_min=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->stars.ti_end_min, c->stars.ti_end_min * e->time_base, e->ti_current,
        e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->stars.ti_end_min == e->ti_current);
}

/**
 * @brief Does a cell contain any sink-particle finishing their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_active_sinks(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->sinks.ti_end_min < e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_end_min=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->sinks.ti_end_min, c->sinks.ti_end_min * e->time_base, e->ti_current,
        e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->sinks.ti_end_min == e->ti_current);
}

/**
 * @brief Does a cell contain any s-particle finishing their time-step now ?
 *
 * This also considers additional physics modules interacting with stars.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @param with_star_formation Are we running with standard star formation?
 * @param with_star_formation_sink Are we running with star formation from
 * sinks?
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_need_activating_stars(
    const struct cell *c, const struct engine *e, const int with_star_formation,
    const int with_star_formation_sink) {

  return cell_is_active_stars(c, e) ||
         (feedback_use_newborn_stars && with_star_formation &&
          cell_is_active_hydro(c, e)) ||
         (with_star_formation_sink &&
          (cell_is_active_sinks(c, e) || cell_is_active_hydro(c, e)));
}

/**
 * @brief Does a cell contain any b-particle finishing their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_active_black_holes(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->black_holes.ti_end_min < e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_end_min=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->black_holes.ti_end_min, c->black_holes.ti_end_min * e->time_base,
        e->ti_current, e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->black_holes.ti_end_min == e->ti_current);
}

/**
 * @brief Is this particle finishing its time-step now ?
 *
 * @param p The #part.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #part is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int part_is_active(
    const struct part *p, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t part_bin = p->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, p->time_bin);
  if (ti_end < ti_current)
    error(
        "particle in an impossible time-zone! p->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (part_bin <= max_active_bin);
}

__attribute__((always_inline)) INLINE static int part_is_active_no_debug(
    const struct part *p, const timebin_t max_active_bin) {

  const timebin_t part_bin = p->time_bin;

  return (part_bin <= max_active_bin);
}

/**
 * @brief Is this particle finishing its RT time-step now ?
 *
 * @param p The #part.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #part is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int part_is_rt_active(
    const struct part *p, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin_subcycle;
  const timebin_t part_bin = p->rt_time_data.time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current_subcycle = e->ti_current_subcycle;
  const integertime_t ti_end =
      get_integer_time_end(ti_current_subcycle, p->rt_time_data.time_bin);
  if (ti_end < ti_current_subcycle)
    error(
        "particle in an impossible time-zone! p->ti_end_subcycle=%lld "
        "e->ti_current_subcycle=%lld",
        ti_end, ti_current_subcycle);
#endif

  return (part_bin <= max_active_bin);
}

/**
 * @brief Is this g-particle finishing its time-step now ?
 *
 * @param gp The #gpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #gpart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int gpart_is_active(
    const struct gpart *gp, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t gpart_bin = gp->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, gp->time_bin);

  if (ti_end < ti_current)
    error(
        "g-particle in an impossible time-zone! gp->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (gpart_bin <= max_active_bin);
}

/**
 * @brief Is this s-particle finishing its time-step now ?
 *
 * @param sp The #spart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #spart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int spart_is_active(
    const struct spart *sp, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t spart_bin = sp->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, sp->time_bin);

  if (ti_end < ti_current)
    error(
        "s-particle in an impossible time-zone! sp->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (spart_bin <= max_active_bin);
}

/**
 * @brief Is this sink-particle finishing its time-step now ?
 *
 * @param sink The #sink.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #bpart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int sink_is_active(
    const struct sink *sink, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t sink_bin = sink->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, sink->time_bin);

  if (ti_end < ti_current)
    error(
        "sink-particle in an impossible time-zone! bp->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (sink_bin <= max_active_bin);
}

/**
 * @brief Is this b-particle finishing its time-step now ?
 *
 * @param bp The #bpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #bpart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int bpart_is_active(
    const struct bpart *bp, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t bpart_bin = bp->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, bp->time_bin);

  if (ti_end < ti_current)
    error(
        "b-particle in an impossible time-zone! bp->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (bpart_bin <= max_active_bin);
}

/**
 * @brief Has this particle been inhibited?
 *
 * @param p The #part.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #part is inhibited, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int part_is_inhibited(
    const struct part *p, const struct engine *e) {
  return p->time_bin == time_bin_inhibited;
}

/**
 * @brief Has this gravity particle been inhibited?
 *
 * @param gp The #gpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #gpart is inhibited, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int gpart_is_inhibited(
    const struct gpart *gp, const struct engine *e) {
  return gp->time_bin == time_bin_inhibited;
}

/**
 * @brief Has this foreign gravity particle been inhibited?
 *
 * @param gp The #gpart_part.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #gpart_foreign is inhibited, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int gpart_foreign_is_inhibited(
    const struct gpart_foreign *gp, const struct engine *e) {
  return gp->time_bin == time_bin_inhibited;
}

/**
 * @brief Has this star particle been inhibited?
 *
 * @param sp The #spart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #spart is inhibited, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int spart_is_inhibited(
    const struct spart *sp, const struct engine *e) {
  return sp->time_bin == time_bin_inhibited;
}

/**
 * @brief Has this sink particle been inhibited?
 *
 * @param sink The #sink.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #sink is inhibited, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int sink_is_inhibited(
    const struct sink *sink, const struct engine *e) {
  return sink->time_bin == time_bin_inhibited;
}

/**
 * @brief Has this black hole particle been inhibited?
 *
 * @param bp The #bpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #bpart is inhibited, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int bpart_is_inhibited(
    const struct bpart *bp, const struct engine *e) {
  return bp->time_bin == time_bin_inhibited;
}

/* Are cells / particles active for kick1 tasks ? */

/**
 * @brief Does a cell contain any particle starting their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_starting_hydro(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.ti_beg_max > e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_beg_max=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->hydro.ti_beg_max, c->hydro.ti_beg_max * e->time_base, e->ti_current,
        e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->hydro.ti_beg_max == e->ti_current);
}

/**
 * @brief Does a cell contain any g-particle starting their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_starting_gravity(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->grav.ti_beg_max > e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_beg_max=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->grav.ti_beg_max, c->grav.ti_beg_max * e->time_base, e->ti_current,
        e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->grav.ti_beg_max == e->ti_current);
}

/**
 * @brief Does a cell contain any s-particle starting their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_starting_stars(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->stars.ti_beg_max > e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_beg_max=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->stars.ti_beg_max, c->stars.ti_beg_max * e->time_base, e->ti_current,
        e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->stars.ti_beg_max == e->ti_current);
}

/**
 * @brief Does a cell contain any sink-particle starting their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_starting_sinks(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->sinks.ti_beg_max > e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_beg_max=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->sinks.ti_beg_max, c->sinks.ti_beg_max * e->time_base, e->ti_current,
        e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->sinks.ti_beg_max == e->ti_current);
}

/**
 * @brief Does a cell contain any b-particle starting their time-step now ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_starting_black_holes(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->black_holes.ti_beg_max > e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_beg_max=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e, a=%e)",
        c->black_holes.ti_beg_max, c->black_holes.ti_beg_max * e->time_base,
        e->ti_current, e->ti_current * e->time_base, e->cosmology->a);
#endif

  return (c->black_holes.ti_beg_max == e->ti_current);
}

/**
 * @brief Is this particle starting its time-step now ?
 *
 * @param p The #part.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #part is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int part_is_starting(
    const struct part *p, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t part_bin = p->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_beg =
      get_integer_time_begin(ti_current + 1, p->time_bin);

  if (ti_beg > ti_current)
    error(
        "particle in an impossible time-zone! p->ti_beg=%lld "
        "e->ti_current=%lld",
        ti_beg, ti_current);
#endif

  return (part_bin <= max_active_bin);
}

/**
 * @brief Is this g-particle starting its time-step now ?
 *
 * @param gp The #gpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #gpart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int gpart_is_starting(
    const struct gpart *gp, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t gpart_bin = gp->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_beg =
      get_integer_time_begin(ti_current + 1, gp->time_bin);

  if (ti_beg > ti_current)
    error(
        "g-particle in an impossible time-zone! gp->ti_beg=%lld "
        "e->ti_current=%lld",
        ti_beg, ti_current);
#endif

  return (gpart_bin <= max_active_bin);
}

/**
 * @brief Is this s-particle starting its time-step now ?
 *
 * @param sp The #spart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #spart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int spart_is_starting(
    const struct spart *sp, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t spart_bin = sp->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_beg =
      get_integer_time_begin(ti_current + 1, sp->time_bin);

  if (ti_beg > ti_current)
    error(
        "s-particle in an impossible time-zone! sp->ti_beg=%lld "
        "e->ti_current=%lld",
        ti_beg, ti_current);
#endif

  return (spart_bin <= max_active_bin);
}

/**
 * @brief Is this b-particle starting its time-step now ?
 *
 * @param bp The #bpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #bpart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int bpart_is_starting(
    const struct bpart *bp, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t bpart_bin = bp->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_beg =
      get_integer_time_begin(ti_current + 1, bp->time_bin);

  if (ti_beg > ti_current)
    error(
        "s-particle in an impossible time-zone! bp->ti_beg=%lld "
        "e->ti_current=%lld",
        ti_beg, ti_current);
#endif

  return (bpart_bin <= max_active_bin);
}

/**
 * @brief Is this sink-particle starting its time-step now ?
 *
 * @param sink The #sink.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #sink is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int sink_is_starting(
    const struct sink *sink, const struct engine *e) {

  const timebin_t max_active_bin = e->max_active_bin;
  const timebin_t sink_bin = sink->time_bin;

#ifdef SWIFT_DEBUG_CHECKS
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_beg =
      get_integer_time_begin(ti_current + 1, sink->time_bin);

  if (ti_beg > ti_current)
    error(
        "sink-particle in an impossible time-zone! sink->ti_beg=%lld "
        "e->ti_current=%lld",
        ti_beg, ti_current);
#endif

  return (sink_bin <= max_active_bin);
}

#endif /* SWIFT_ACTIVE_H */
