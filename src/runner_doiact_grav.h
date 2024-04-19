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
#include "inline.h"

/* Avoid cyclic inclusions. */
struct runner;
struct cell;

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
 * Defined here to enable inlining while this function is called from
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

#ifndef SWIFT_TASKS_WITHOUT_ATOMICS
  /* Unlock the multipoles */
  if (lock_unlock(&ci->grav.mlock) != 0) error("Failed to unlock multipole");
  if (lock_unlock(&cj->grav.mlock) != 0) error("Failed to unlock multipole");
#endif

  TIMER_TOC(timer_dopair_grav_mm);
}

#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
