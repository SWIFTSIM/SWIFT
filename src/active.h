/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* Local includes. */
#include "cell.h"
#include "engine.h"
#include "part.h"
#include "timeline.h"

/**
 * @brief Check that a cell been drifted to the current time.
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell has been drifted to the current time, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_drifted(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->ti_old > e->ti_current)
    error(
        "Cell has been drifted too far forward in time! c->ti_old=%lld (t=%e) "
        "and e->ti_current=%lld (t=%e)",
        c->ti_old, c->ti_old * e->timeBase, e->ti_current,
        e->ti_current * e->timeBase);
#endif

  return (c->ti_old == e->ti_current);
}

/**
 * @brief Does a cell contain any active particle ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #cell contains at least an active particle, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_active(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->ti_end_min < e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_end_min=%lld (t=%e) and "
        "e->ti_current=%lld (t=%e)",
        c->ti_end_min, c->ti_end_min * e->timeBase, e->ti_current,
        e->ti_current * e->timeBase);
#endif

  return (c->ti_end_min == e->ti_current);
}

/**
 * @brief Are *all* particles in a cell active ?
 *
 * @param c The #cell.
 * @param e The #engine containing information about the current time.
 * @return 1 if all particles in a #cell are active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int cell_is_all_active(
    const struct cell *c, const struct engine *e) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->ti_end_max < e->ti_current)
    error(
        "cell in an impossible time-zone! c->ti_end_max=%lld "
        "e->ti_current=%lld",
        c->ti_end_max, e->ti_current);
#endif

  return (c->ti_end_max == e->ti_current);
}

/**
 * @brief Is this particle active ?
 *
 * @param p The #part.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #part is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int part_is_active(
    const struct part *p, const struct engine *e) {

  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, p->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_current)
    error(
        "particle in an impossible time-zone! p->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (ti_end == ti_current);
}

/**
 * @brief Is this g-particle active ?
 *
 * @param gp The #gpart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #gpart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int gpart_is_active(
    const struct gpart *gp, const struct engine *e) {

  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, gp->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_current)
    error(
        "g-particle in an impossible time-zone! gp->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (ti_end == ti_current);
}

/**
 * @brief Is this s-particle active ?
 *
 * @param sp The #spart.
 * @param e The #engine containing information about the current time.
 * @return 1 if the #spart is active, 0 otherwise.
 */
__attribute__((always_inline)) INLINE static int spart_is_active(
    const struct spart *sp, const struct engine *e) {

  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = get_integer_time_end(ti_current, sp->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_end < ti_current)
    error(
        "s-particle in an impossible time-zone! gp->ti_end=%lld "
        "e->ti_current=%lld",
        ti_end, ti_current);
#endif

  return (ti_end == ti_current);
}

#endif /* SWIFT_ACTIVE_H */
