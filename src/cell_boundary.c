/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 SWIFT developers.
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
#include <config.h>

#include "cell_boundary.h"

#include "active.h"
#include "cell.h"
#include "csds.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "lock.h"

/** Apply open box boundaries to gas particles in a cell hierarchy. */
void cell_boundary_part(struct cell *c, const struct engine *e) {
  if (c->hydro.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_part(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  for (int k = 0; k < c->hydro.count; k++) {
    struct part *p = &c->hydro.parts[k];
    struct xpart *xp = &c->hydro.xparts[k];
    if (part_is_inhibited(p, e) || !cell_boundary_outside(p->x, dim)) continue;

    lock_lock(&e->s->lock);
    if (!part_is_inhibited(p, e)) {
#ifdef WITH_CSDS
      if (e->policy & engine_policy_csds)
        csds_log_part(e->csds, p, xp, e, /*log_all=*/1, csds_flag_delete,
                      /*data=*/0);
#endif
      hydro_remove_part(p, xp, e->time);
      cell_remove_part(e, c, p, xp);
    }
    if (lock_unlock(&e->s->lock) != 0) error("Failed to unlock the space!");
  }
}

/** Apply open box boundaries to gravity particles in a cell hierarchy. */
void cell_boundary_gpart(struct cell *c, const struct engine *e) {
  if (c->grav.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_gpart(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  for (int k = 0; k < c->grav.count; k++) {
    struct gpart *gp = &c->grav.parts[k];
    if (gpart_is_inhibited(gp, e) || !cell_boundary_outside(gp->x, dim))
      continue;
    /* Baryonic gparts are removed with their counterparts. */
    if (gp->type != swift_type_dark_matter) continue;

    lock_lock(&e->s->lock);
    if (!gpart_is_inhibited(gp, e)) {
#ifdef WITH_CSDS
      if (e->policy & engine_policy_csds)
        csds_log_gpart(e->csds, gp, e, /*log_all=*/1, csds_flag_delete,
                       /*data=*/0);
#endif
      cell_remove_gpart(e, c, gp);
    }
    if (lock_unlock(&e->s->lock) != 0) error("Failed to unlock the space!");
  }
}

/** Apply open box boundaries to star particles in a cell hierarchy. */
void cell_boundary_spart(struct cell *c, const struct engine *e) {
  if (c->stars.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_spart(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  for (int k = 0; k < c->stars.count; k++) {
    struct spart *sp = &c->stars.parts[k];
    if (spart_is_inhibited(sp, e) || !cell_boundary_outside(sp->x, dim))
      continue;

    lock_lock(&e->s->lock);
    if (!spart_is_inhibited(sp, e)) {
#ifdef WITH_CSDS
      if (e->policy & engine_policy_csds)
        csds_log_spart(e->csds, sp, e, /*log_all=*/1, csds_flag_delete,
                       /*data=*/0);
#endif
      cell_remove_spart(e, c, sp);
    }
    if (lock_unlock(&e->s->lock) != 0) error("Failed to unlock the space!");
  }
}

/** Apply open box boundaries to black hole particles in a cell hierarchy. */
void cell_boundary_bpart(struct cell *c, const struct engine *e) {
  if (c->black_holes.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_bpart(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  for (int k = 0; k < c->black_holes.count; k++) {
    struct bpart *bp = &c->black_holes.parts[k];
    if (bpart_is_inhibited(bp, e) || !cell_boundary_outside(bp->x, dim))
      continue;

    lock_lock(&e->s->lock);
    if (!bpart_is_inhibited(bp, e)) {
#ifdef WITH_CSDS
      if (e->policy & engine_policy_csds)
        error("Logging of black hole particles is not yet implemented.");
#endif
      cell_remove_bpart(e, c, bp);
    }
    if (lock_unlock(&e->s->lock) != 0) error("Failed to unlock the space!");
  }
}

/** Apply open box boundaries to sink particles in a cell hierarchy. */
void cell_boundary_sink(struct cell *c, const struct engine *e) {
  if (c->sinks.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_sink(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  for (int k = 0; k < c->sinks.count; k++) {
    struct sink *sink = &c->sinks.parts[k];
    if (sink_is_inhibited(sink, e) || !cell_boundary_outside(sink->x, dim))
      continue;

    lock_lock(&e->s->lock);
    if (!sink_is_inhibited(sink, e)) {
#ifdef WITH_CSDS
      if (e->policy & engine_policy_csds)
        error("Logging of sink particles is not yet implemented.");
#endif
      cell_remove_sink(e, c, sink);
    }
    if (lock_unlock(&e->s->lock) != 0) error("Failed to unlock the space!");
  }
}
