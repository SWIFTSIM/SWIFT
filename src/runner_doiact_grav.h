/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_RUNNER_DOIACT_GRAV_H
#define SWIFT_RUNNER_DOIACT_GRAV_H

/* Config */
#include <config.h>

/* Local includes */
#include "active.h"
#include "atomic.h"
#include "cell.h"
#include "inline.h"
#include "multipole.h"
#include "timers.h"

/* Avoid cyclic inclusions. */
struct runner;
struct cell;

#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_GRAVITY_FORCE_CHECKS)
enum runner_debug_source_class {
  runner_debug_source_class_zoom = 0,
  runner_debug_source_class_bkg_void = 1,
  runner_debug_source_class_bkg_neigh = 2,
  runner_debug_source_class_other = 3,
};

__attribute__((always_inline)) INLINE static enum runner_debug_source_class
runner_debug_get_source_class(const struct cell *source) {

  if (source->type == cell_type_zoom) return runner_debug_source_class_zoom;

  if (source->type == cell_type_bkg && source->subtype == cell_subtype_void)
    return runner_debug_source_class_bkg_void;

  if (source->type == cell_type_bkg &&
      source->subtype == cell_subtype_neighbour)
    return runner_debug_source_class_bkg_neigh;

  return runner_debug_source_class_other;
}

__attribute__((always_inline)) INLINE static void
runner_debug_add_gpart_interactions_by_type(long long counts[4],
                                            const struct cell *source,
                                            const long long delta) {
  counts[runner_debug_get_source_class(source)] += delta;
}

__attribute__((always_inline)) INLINE static void
runner_debug_inc_gpart_interactions_by_type(long long counts[4],
                                            const struct cell *source) {
  counts[runner_debug_get_source_class(source)] += 1;
}
#endif

#ifdef SWIFT_DEBUG_CHECKS
__attribute__((always_inline)) INLINE static void
runner_debug_add_tensor_interactions_by_type(long long counts[4],
                                             const struct cell *source,
                                             const long long delta) {
  counts[runner_debug_get_source_class(source)] += delta;
}

enum runner_debug_coverage_kind {
  runner_debug_coverage_kind_p2p = 0,
  runner_debug_coverage_kind_mm = 1,
  runner_debug_coverage_kind_pm = 2,
};

__attribute__((always_inline)) INLINE static double
runner_debug_get_cell_fractional_weight(const struct cell *source) {

  double weight = 1.;
  for (int depth = source->depth; depth > 0; --depth) weight *= 0.125;
  return weight;
}

__attribute__((always_inline)) INLINE static void runner_debug_add_cell_coverage(
    struct cell *recipient, const struct cell *source,
    const enum runner_debug_coverage_kind kind) {

  const int source_class = runner_debug_get_source_class(source);
  const double weight = runner_debug_get_cell_fractional_weight(source);

  atomic_add_d(&recipient->num_interacted_cells_total, weight);
  atomic_add_d(&recipient->num_interacted_cells_total_by_type[source_class],
               weight);

  if (kind == runner_debug_coverage_kind_p2p) {
    atomic_add_d(&recipient->num_interacted_cells_p2p, weight);
    atomic_add_d(&recipient->num_interacted_cells_p2p_by_type[source_class],
                 weight);
  } else if (kind == runner_debug_coverage_kind_mm) {
    atomic_add_d(&recipient->num_interacted_cells_mm, weight);
    atomic_add_d(&recipient->num_interacted_cells_mm_by_type[source_class],
                 weight);
  } else if (kind == runner_debug_coverage_kind_pm) {
    atomic_add_d(&recipient->num_interacted_cells_pm, weight);
    atomic_add_d(&recipient->num_interacted_cells_pm_by_type[source_class],
                 weight);
  }
}

__attribute__((always_inline)) INLINE static void
runner_debug_inherit_cell_interactions(struct cell *child,
                                       const struct cell *parent) {

  child->num_interacted_cells_total += parent->num_interacted_cells_mm;
  child->num_interacted_cells_total += parent->num_interacted_cells_pm;
  child->num_interacted_cells_mm += parent->num_interacted_cells_mm;
  child->num_interacted_cells_pm += parent->num_interacted_cells_pm;

  for (int i = 0; i < 4; ++i) {
    child->num_interacted_cells_total_by_type[i] +=
        parent->num_interacted_cells_mm_by_type[i];
    child->num_interacted_cells_total_by_type[i] +=
        parent->num_interacted_cells_pm_by_type[i];
    child->num_interacted_cells_mm_by_type[i] +=
        parent->num_interacted_cells_mm_by_type[i];
    child->num_interacted_cells_pm_by_type[i] +=
        parent->num_interacted_cells_pm_by_type[i];
  }
}

void runner_debug_get_top_level_methods_by_type(const struct engine *e,
                                                const struct cell *ci,
                                                long long mesh_counts[4],
                                                long long mm_counts[4],
                                                long long p2p_counts[4],
                                                int mesh_nr_cells[4],
                                                int mm_nr_cells[4],
                                                int p2p_nr_cells[4]);

int runner_debug_dump_recursive_mesh_budget_overlaps(const struct cell *ci);
int runner_debug_dump_zoom_pair_recursive_revisits(const struct cell *ci);
int runner_debug_dump_zoom_pair_handoff_overlaps(const struct cell *ci);
int runner_debug_dump_path_pair_recursive_sources(const struct cell *ci);
#endif

void runner_do_grav_down(struct runner *r, struct cell *c, int timer);

void runner_dopair_grav_pp(struct runner *r, struct cell *ci, struct cell *cj,
                           const int symmetric, const int allow_mpole);

void runner_doself_recursive_grav(struct runner *r, struct cell *c,
                                  int gettimer);

void runner_dopair_recursive_grav(struct runner *r, struct cell *ci,
                                  struct cell *cj, int gettimer);
void runner_dopair_grav_mm_progenies(struct runner *r, const long long flags,
                                     struct cell *restrict ci,
                                     struct cell *restrict cj);

void runner_do_grav_long_range(struct runner *r, struct cell *ci, int timer);

/* Internal functions (for unit tests and debugging) */

void runner_doself_grav_pp(struct runner *r, struct cell *c);

void runner_dopair_grav_pp(struct runner *r, struct cell *ci, struct cell *cj,
                           const int symmetric, const int allow_mpole);

/**
 * @brief Computes the interaction of the field tensor in a cell with the
 * multipole of another cell.
 *
 * Defined here to enable inlining while this function is called from both
 * runner_doiact_grav.c and runner_doiact_long_range_grav.c.
 *
 * @param r The #runner.
 * @param ci The #cell with field tensor to interact.
 * @param cj The #cell with the multipole.
 */
static INLINE void runner_dopair_grav_mm_nonsym(struct runner *r,
                                                struct cell *restrict ci,
                                                struct cell *restrict cj) {

  /* Some constants */
  const struct engine *e = r->e;
  const struct gravity_props *props = e->gravity_properties;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};
  const float r_s_inv = e->mesh->r_s_inv;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity_mm(ci, e) || ci->nodeID != engine_rank) return;

  /* Short-cut to the multipole */
  const struct multipole *multi_j = &cj->grav.multipole->m_pole;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Interacting a cell with itself using M2L");

  if (multi_j->num_gpart == 0)
    error("Multipole does not seem to have been set.");

  if (ci->grav.multipole->pot.ti_init != e->ti_current)
    error("ci->grav tensor not initialised.");

  if (cj->grav.ti_old_multipole != e->ti_current)
    error(
        "Undrifted multipole cj->grav.ti_old_multipole=%lld cj->nodeID=%d "
        "ci->nodeID=%d e->ti_current=%lld",
        cj->grav.ti_old_multipole, cj->nodeID, ci->nodeID, e->ti_current);
#endif

#ifndef SWIFT_TASKS_WITHOUT_ATOMICS
  /* Lock the multipoles
   * Note we impose a hierarchy to solve the dining philosopher problem */
  if (ci < cj) {
    lock_lock(&ci->grav.mlock);
    lock_lock(&cj->grav.mlock);
  } else {
    lock_lock(&cj->grav.mlock);
    lock_lock(&ci->grav.mlock);
  }
#endif

  /* Let's interact at this level */
  gravity_M2L_nonsym(&ci->grav.multipole->pot, multi_j, ci->grav.multipole->CoM,
                     cj->grav.multipole->CoM, props, periodic, dim, r_s_inv);

#ifdef SWIFT_DEBUG_CHECKS
  runner_debug_add_tensor_interactions_by_type(
      ci->grav.multipole->pot.num_interacted_tree_by_type, cj,
      multi_j->num_gpart);
  runner_debug_add_cell_coverage(ci, cj, runner_debug_coverage_kind_mm);
#endif

#ifndef SWIFT_TASKS_WITHOUT_ATOMICS
  /* Unlock the multipoles */
  if (lock_unlock(&ci->grav.mlock) != 0) error("Failed to unlock multipole");
  if (lock_unlock(&cj->grav.mlock) != 0) error("Failed to unlock multipole");
#endif

  TIMER_TOC(timer_dopair_grav_mm);
}

/**
 * @brief Computes the interaction of the field tensor and multipole
 * of two cells symmetrically.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
static INLINE void runner_dopair_grav_mm_symmetric(struct runner *r,
                                                   struct cell *restrict ci,
                                                   struct cell *restrict cj) {

  /* Some constants */
  const struct engine *e = r->e;
  const struct gravity_props *props = e->gravity_properties;
  const int periodic = e->mesh->periodic;
  const double dim[3] = {e->mesh->dim[0], e->mesh->dim[1], e->mesh->dim[2]};
  const float r_s_inv = e->mesh->r_s_inv;

  TIMER_TIC;

  /* Anything to do here? */
  if ((!cell_is_active_gravity_mm(ci, e) || ci->nodeID != engine_rank) ||
      (!cell_is_active_gravity_mm(cj, e) || cj->nodeID != engine_rank))
    error("Invalid state in symmetric M-M calculation!");

  /* Short-cut to the multipole */
  const struct multipole *multi_i = &ci->grav.multipole->m_pole;
  const struct multipole *multi_j = &cj->grav.multipole->m_pole;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci == cj) error("Interacting a cell with itself using M2L");

  if (multi_i->num_gpart == 0)
    error("Multipole i does not seem to have been set.");

  if (multi_j->num_gpart == 0)
    error("Multipole j does not seem to have been set.");

  if (ci->grav.multipole->pot.ti_init != e->ti_current)
    error("ci->grav tensor not initialised.");

  if (ci->grav.multipole->pot.ti_init != e->ti_current)
    error("cj->grav tensor not initialised.");

  if (ci->grav.ti_old_multipole != e->ti_current)
    error(
        "Undrifted multipole ci->grav.ti_old_multipole=%lld ci->nodeID=%d "
        "cj->nodeID=%d e->ti_current=%lld",
        ci->grav.ti_old_multipole, ci->nodeID, cj->nodeID, e->ti_current);

  if (cj->grav.ti_old_multipole != e->ti_current)
    error(
        "Undrifted multipole cj->grav.ti_old_multipole=%lld cj->nodeID=%d "
        "ci->nodeID=%d e->ti_current=%lld",
        cj->grav.ti_old_multipole, cj->nodeID, ci->nodeID, e->ti_current);
#endif

#ifndef SWIFT_TASKS_WITHOUT_ATOMICS
  /* Lock the multipoles
   * Note we impose a hierarchy to solve the dining philosopher problem */
  if (ci < cj) {
    lock_lock(&ci->grav.mlock);
    lock_lock(&cj->grav.mlock);
  } else {
    lock_lock(&cj->grav.mlock);
    lock_lock(&ci->grav.mlock);
  }
#endif

  /* Let's interact at this level */
  gravity_M2L_symmetric(&ci->grav.multipole->pot, &cj->grav.multipole->pot,
                        multi_i, multi_j, ci->grav.multipole->CoM,
                        cj->grav.multipole->CoM, props, periodic, dim, r_s_inv);

#ifdef SWIFT_DEBUG_CHECKS
  runner_debug_add_tensor_interactions_by_type(
      ci->grav.multipole->pot.num_interacted_tree_by_type, cj,
      multi_j->num_gpart);
  runner_debug_add_tensor_interactions_by_type(
      cj->grav.multipole->pot.num_interacted_tree_by_type, ci,
      multi_i->num_gpart);
  runner_debug_add_cell_coverage(ci, cj, runner_debug_coverage_kind_mm);
  runner_debug_add_cell_coverage(cj, ci, runner_debug_coverage_kind_mm);
#endif

#ifndef SWIFT_TASKS_WITHOUT_ATOMICS
  /* Unlock the multipoles */
  if (lock_unlock(&ci->grav.mlock) != 0) error("Failed to unlock multipole");
  if (lock_unlock(&cj->grav.mlock) != 0) error("Failed to unlock multipole");
#endif

  TIMER_TOC(timer_dopair_grav_mm);
}

/**
 * @brief Call the M-M calculation on two cells if active.
 *
 * @param r The #runner object.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
static INLINE void runner_dopair_grav_mm(struct runner *r,
                                         struct cell *restrict ci,
                                         struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* What do we need to do? */
  const int do_i =
      cell_is_active_gravity_mm(ci, e) && (ci->nodeID == e->nodeID);
  const int do_j =
      cell_is_active_gravity_mm(cj, e) && (cj->nodeID == e->nodeID);

  /* Do we need drifting first? */
  if (ci->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(ci, e);
  if (cj->grav.ti_old_multipole < e->ti_current) cell_drift_multipole(cj, e);

  /* Interact! */
  if (do_i && do_j)
    runner_dopair_grav_mm_symmetric(r, ci, cj);
  else if (do_i)
    runner_dopair_grav_mm_nonsym(r, ci, cj);
  else if (do_j)
    runner_dopair_grav_mm_nonsym(r, cj, ci);
}

#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
