/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Config parameters. */
#include <config.h>

/* This object's header. */
#include "cell.h"

/* Local headers */
#include "active.h"
#include "part_init.h"

/**
 * @brief Recursively initialise all gas particles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine.
 */
void cell_init_part(struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID) error("Calling function on foreign cell!");
#endif

  if (c->split) {
    /* Loop over the progeny and recurse. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_init_part(c->progeny[k], e);
      }
    }
  } else {

    /* Loop over all the gas particles in the cell */
    const size_t nr_parts = c->hydro.count;
    for (size_t k = 0; k < nr_parts; k++) {
      /* Get a handle on the part. */
      struct part *const p = &c->hydro.parts[k];
      struct xpart *const xp = &c->hydro.xparts[k];

      /* Ignore inhibited particles */
      if (part_is_inhibited(p, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (p->ti_drift != e->ti_current) {
        error(
            "Trying to initialise a particle not drifted to the current time!");
      }
#endif

      if (part_is_active(p, e)) {
        part_init(p, xp, e);
      }
    }
  }
}

/**
 * @brief Recursively initialise all gravity particles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine.
 */
void cell_init_gpart(struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID) error("Calling function on foreign cell!");
#endif

  if (c->split) {
    /* Loop over the progeny and recurse. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_init_gpart(c->progeny[k], e);
      }
    }
  } else {

    /* Loop over all the gas particles in the cell */
    const size_t nr_gparts = c->grav.count;
    for (size_t k = 0; k < nr_gparts; k++) {
      /* Get a handle on the part. */
      struct gpart *const gp = &c->grav.parts[k];

      /* Ignore inhibited particles */
      if (gpart_is_inhibited(gp, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (gp->ti_drift != e->ti_current) {
        error(
            "Trying to initialise a particle not drifted to the current time!");
      }
#endif

      if (gpart_is_active(gp, e)) {
        gpart_init(gp, e);
      }
    }
  }
}

/**
 * @brief Recursively initialise all star particles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine.
 */
void cell_init_spart(struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID) error("Calling function on foreign cell!");
#endif

  if (c->split) {
    /* Loop over the progeny and recurse. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_init_spart(c->progeny[k], e);
      }
    }
  } else {

    /* Loop over all the gas particles in the cell */
    const size_t nr_sparts = c->stars.count;
    for (size_t k = 0; k < nr_sparts; k++) {
      /* Get a handle on the part. */
      struct spart *const sp = &c->stars.parts[k];

      /* Ignore inhibited particles */
      if (spart_is_inhibited(sp, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (sp->ti_drift != e->ti_current) {
        error(
            "Trying to initialise a particle not drifted to the current time!");
      }
#endif

      if (spart_is_active(sp, e)) {
        spart_init(sp, e);
      }
    }
  }
}

/**
 * @brief Recursively initialise all black hole particles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine.
 */
void cell_init_bpart(struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID) error("Calling function on foreign cell!");
#endif

  if (c->split) {
    /* Loop over the progeny and recurse. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_init_bpart(c->progeny[k], e);
      }
    }
  } else {

    /* Loop over all the gas particles in the cell */
    const size_t nr_bparts = c->black_holes.count;
    for (size_t k = 0; k < nr_bparts; k++) {
      /* Get a handle on the part. */
      struct bpart *const bp = &c->black_holes.parts[k];

      /* Ignore inhibited particles */
      if (bpart_is_inhibited(bp, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (bp->ti_drift != e->ti_current) {
        error(
            "Trying to initialise a particle not drifted to the current time!");
      }
#endif

      if (bpart_is_active(bp, e)) {
        bpart_init(bp, e);
      }
    }
  }
}

/**
 * @brief Recursively initialise all sink particles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine.
 */
void cell_init_sink(struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID) error("Calling function on foreign cell!");
#endif

  if (c->split) {
    /* Loop over the progeny and recurse. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        cell_init_sink(c->progeny[k], e);
      }
    }
  } else {

    /* Loop over all the gas particles in the cell */
    const size_t nr_sinks = c->sinks.count;
    for (size_t k = 0; k < nr_sinks; k++) {
      /* Get a handle on the part. */
      struct sink *const si = &c->sinks.parts[k];

      /* Ignore inhibited particles */
      if (sink_is_inhibited(si, e)) continue;

#ifdef SWIFT_DEBUG_CHECKS
      if (si->ti_drift != e->ti_current) {
        error(
            "Trying to initialise a particle not drifted to the current time!");
      }
#endif

      if (sink_is_active(si, e)) {
        sink_init(si, e);
      }
    }
  }
}
