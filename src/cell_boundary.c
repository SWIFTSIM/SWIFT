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
#include "cell_boundary.h"

#include "active.h"
#include "cell.h"
#include "csds.h"
#include "engine.h"
#include "error.h"
#include "hydro.h"
#include "lock.h"
#include "potential.h"

#include <config.h>

/** Apply open and reflective box boundaries to gas particles in a cell
 * hierarchy. */
void cell_boundary_part(struct cell *c, const struct engine *e) {
  if (c->hydro.count == 0 || c->subtype == cell_subtype_mirror) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_part(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  const unsigned int walls = e->s->boundary_conditions.walls;
  double shift[3];
  for (int k = 0; k < c->hydro.count; k++) {
    struct part *p = &c->hydro.parts[k];
    struct xpart *xp = &c->hydro.xparts[k];
    if (part_is_inhibited(p, e)) continue;

    /* Mirror particles that crossed a wall, unless they also left the box
     * through an open face. */
    if (!cell_boundary_outside(p->x, dim, walls)) {
      if (cell_boundary_reflect(p->x, dim, walls, shift)) {
        for (int j = 0; j < 3; j++) {
          if (shift[j] == 0.) continue;
          p->v[j] = -p->v[j];
          xp->v_full[j] = -xp->v_full[j];
          xp->x_diff[j] -= shift[j];
          xp->x_diff_sort[j] -= shift[j];
        }
      }
      continue;
    }

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

/** Apply open and reflective box boundaries to gravity particles in a cell
 * hierarchy. */
void cell_boundary_gpart(struct cell *c, const struct engine *e) {
  if (c->grav.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_gpart(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  const unsigned int walls = e->s->boundary_conditions.walls;
  double shift[3];
  for (int k = 0; k < c->grav.count; k++) {
    struct gpart *gp = &c->grav.parts[k];
    if (gpart_is_inhibited(gp, e)) continue;

    /* Mirror particles that crossed a wall, unless they also left the box
     * through an open face. */
    if (!cell_boundary_outside(gp->x, dim, walls)) {
      if (cell_boundary_reflect(gp->x, dim, walls, shift)) {
        for (int j = 0; j < 3; j++) {
          if (shift[j] == 0.) continue;
          gp->v_full[j] = -gp->v_full[j];
        }
      }
      continue;
    }
    /* Baryonic gparts are removed with their counterparts (but are mirrored
     * here, as their drift matches their counterparts'). */
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

/** Apply open and reflective box boundaries to star particles in a cell
 * hierarchy. */
void cell_boundary_spart(struct cell *c, const struct engine *e) {
  if (c->stars.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_spart(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  const unsigned int walls = e->s->boundary_conditions.walls;
  double shift[3];
  for (int k = 0; k < c->stars.count; k++) {
    struct spart *sp = &c->stars.parts[k];
    if (spart_is_inhibited(sp, e)) continue;

    /* Mirror particles that crossed a wall, unless they also left the box
     * through an open face. */
    if (!cell_boundary_outside(sp->x, dim, walls)) {
      if (cell_boundary_reflect(sp->x, dim, walls, shift)) {
        for (int j = 0; j < 3; j++) {
          if (shift[j] == 0.) continue;
          sp->v[j] = -sp->v[j];
          sp->x_diff[j] -= shift[j];
          sp->x_diff_sort[j] -= shift[j];
        }
      }
      continue;
    }

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

/** Apply open and reflective box boundaries to black hole particles in a cell
 * hierarchy. */
void cell_boundary_bpart(struct cell *c, const struct engine *e) {
  if (c->black_holes.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_bpart(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  const unsigned int walls = e->s->boundary_conditions.walls;
  double shift[3];
  for (int k = 0; k < c->black_holes.count; k++) {
    struct bpart *bp = &c->black_holes.parts[k];
    if (bpart_is_inhibited(bp, e)) continue;

    /* Mirror particles that crossed a wall, unless they also left the box
     * through an open face. */
    if (!cell_boundary_outside(bp->x, dim, walls)) {
      if (cell_boundary_reflect(bp->x, dim, walls, shift)) {
        for (int j = 0; j < 3; j++) {
          if (shift[j] == 0.) continue;
          bp->v[j] = -bp->v[j];
          bp->x_diff[j] -= shift[j];
        }
      }
      continue;
    }

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

/** Apply open and reflective box boundaries to sink particles in a cell
 * hierarchy. */
void cell_boundary_sink(struct cell *c, const struct engine *e) {
  if (c->sinks.count == 0) return;
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_boundary_sink(c->progeny[k], e);
    return;
  }

  const double *dim = e->s->dim;
  const unsigned int walls = e->s->boundary_conditions.walls;
  double shift[3];
  for (int k = 0; k < c->sinks.count; k++) {
    struct sink *sink = &c->sinks.parts[k];
    if (sink_is_inhibited(sink, e)) continue;

    /* Mirror particles that crossed a wall, unless they also left the box
     * through an open face. */
    if (!cell_boundary_outside(sink->x, dim, walls)) {
      if (cell_boundary_reflect(sink->x, dim, walls, shift)) {
        for (int j = 0; j < 3; j++) {
          if (shift[j] == 0.) continue;
          sink->v[j] = -sink->v[j];
          sink->x_diff[j] -= shift[j];
        }
      }
      continue;
    }

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

/**
 * @brief Refresh the mirror particles of a mirror cell hierarchy.
 *
 * Each mirror particle becomes a copy of the particle with the same index in
 * the leaf cell this one reflects, with its position reflected in the walls
 * and its vectors normal to them reversed (free slip). The cells take the
 * smoothing lengths, displacements and drift time of the cells they reflect.
 *
 * @param c The mirror #cell.
 * Under external gravity the pressure of an image continues the hydrostatic
 * profile of the fluid through the wall, P_image = P + rho g . (x_image - x),
 * or a fluid at rest would sink into the walls it lies on.
 *
 * @param s The #space.
 * @param moved Have the particles moved since the last refresh? If so the
 * sorts of the mirror cells are stale.
 */
void cell_mirror_refresh(struct cell *c, const struct space *s,
                         const int moved) {

  /* ponytail: mirror cells are re-sorted every step they are used; track
   * their displacement since the last sort if that ever shows in a profile. */
  if (moved) c->hydro.sorted = 0;

  if (c->split) {
    float h_max = 0.f, h_max_active = 0.f, dx_max_part = 0.f;
    integertime_t ti_old_part = max_nr_timesteps;
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      if (cp == NULL) continue;
      cell_mirror_refresh(cp, s, moved);
      h_max = max(h_max, cp->hydro.h_max);
      h_max_active = max(h_max_active, cp->hydro.h_max_active);
      dx_max_part = max(dx_max_part, cp->hydro.dx_max_part);
      ti_old_part = min(ti_old_part, cp->hydro.ti_old_part);
    }
    c->hydro.h_max = h_max;
    c->hydro.h_max_active = h_max_active;
    c->hydro.dx_max_part = dx_max_part;
    c->hydro.ti_old_part = ti_old_part;
    return;
  }

  const struct cell *real = c->mirror_of;
  c->hydro.h_max = real->hydro.h_max;
  c->hydro.h_max_active = real->hydro.h_max_active;
  c->hydro.dx_max_part = real->hydro.dx_max_part;
  c->hydro.ti_old_part = real->hydro.ti_old_part;

  const struct part *restrict parts = real->hydro.parts;
  struct part *restrict mirrors = c->hydro.parts;
  const int count = c->hydro.count;
  const unsigned int walls = c->mirror_walls;
  const double *dim = s->dim;
  const struct engine *e = s->e;
  const int with_gravity =
      (e != NULL) && (e->policy & engine_policy_external_gravity);
  const float G_newton =
      with_gravity ? e->physical_constants->const_newton_G : 0.f;

  for (int i = 0; i < count; i++) {

    /* The external acceleration at the particle (stored without G) */
    struct gpart gp;
    gp.a_grav[0] = gp.a_grav[1] = gp.a_grav[2] = 0.f;
    if (with_gravity) {
      gp.x[0] = parts[i].x[0];
      gp.x[1] = parts[i].x[1];
      gp.x[2] = parts[i].x[2];
      external_gravity_acceleration(e->time, e->external_potential,
                                    e->physical_constants, &gp);
    }

    /* Reflect a copy and store it in one go: the ghosts of neighbouring
     * cells may be reading this mirror particle's (unchanged) position. */
    struct part p = parts[i];
    p.gpart = NULL;
    for (int k = 0; k < 3; k++) {
      if (walls & (1u << (2 * k))) {
        p.x[k] = -p.x[k];
      } else if (walls & (1u << (2 * k + 1))) {
        p.x[k] = 2. * dim[k] - p.x[k];
      } else {
        continue;
      }
      p.v[k] = -p.v[k];
      hydro_mirror_part(
          &p, k, G_newton * gp.a_grav[k] * (float)(p.x[k] - parts[i].x[k]));
    }
    mirrors[i] = p;
  }
}
