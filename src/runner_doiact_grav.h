/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

/* Includes. */
#include "cell.h"
#include "clocks.h"
#include "part.h"

/**
 * @brief Compute the recursive upward sweep, i.e. construct the
 *        multipoles in a cell hierarchy.
 *
 * @param r The #runner.
 * @param c The top-level #cell.
 */
void runner_dograv_up(struct runner *r, struct cell *c) {

  if (c->split) {/* Regular node */

    /* Recurse. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_dograv_up(r, c->progeny[k]);

    /* Collect the multipoles from the progeny. */
    multipole_reset(&c->multipole);
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL)
        multipole_add(&c->multipole, &c->progeny[k]->multipole);
    }

  } else {/* Leaf node. */

    /* Just construct the multipole from the gparts. */
    multipole_init(&c->multipole, c->gparts, c->gcount);
  }
}

/**
 * @brief Checks whether the cells are direct neighbours ot not. Both cells have
 * to be of the same size
 * 
 */
__attribute__((always_inline)) INLINE static
int are_neighbours(const struct cell *restrict ci, const struct cell *restrict cj) {

#ifdef SANITY_CHECKS
  if (ci->h[0] != cj->h[0])
    error(" Cells of different size in distance calculation.");
#endif

  /* Maximum allowed distance */
  const float min_dist = ci->h[0];

  /* (Manhattan) Distance between the cells */
  for (int k = 0; k < 3; k++) {
    const float center_i = ci->loc[k];
    const float center_j = cj->loc[k];
    if (fabsf(center_i - center_j) > min_dist) return 0;
  }

  return 1;
}



/**
 * @brief Computes the interaction of all the particles in a cell with the
 * multipole of another cell.
 *
 * @param r The #runner.
 * @param ci The #cell with particles to interct.
 * @param cj The #cell with the multipole.
 */
__attribute__((always_inline)) INLINE static void runner_dograv_pair_pm(
	  const struct runner *r, const struct cell *restrict ci, const struct cell *restrict cj) {

  const struct engine *e = r->e;
  const int gcount = ci->gcount;
  struct gpart *restrict gparts = ci->gparts;
  const struct multipole multi = cj->multipole;
  const int ti_current = e->ti_current;

  TIMER_TIC;

#ifdef SANITY_CHECK
  if (gcount == 0) error("Empty cell!");  // MATTHIEU sanity check

  if (multi.mass == 0.0)  // MATTHIEU sanity check
    error("Multipole does not seem to have been set.");
#endif

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current) return;

  /* Loop over every particle in leaf. */
  for (int pid = 0; pid < gcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gp = &gparts[pid];

    if (gp->ti_end > ti_current) continue;

    /* Compute the pairwise distance. */
    const float dx[3] = {multi.CoM[0] - gp->x[0],   // x
                         multi.CoM[1] - gp->x[1],   // y
                         multi.CoM[2] - gp->x[2]};  // z
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    /* Apply the gravitational acceleration. */
    const float ir = 1.0f / sqrtf(r2);
    const float mrinv3 = multi.mass * ir * ir * ir;

#if multipole_order < 2

    /* 0th and 1st order terms */
    gp->a_grav[0] += mrinv3 * dx[0];
    gp->a_grav[1] += mrinv3 * dx[1];
    gp->a_grav[2] += mrinv3 * dx[2];

#elif multipole_order == 2

    /* Terms up to 2nd order (quadrupole) */

    /* Follows the notation in Bonsai */
    const float mrinv5 = mrinv3 * ir * ir;
    const float mrinv7 = mrinv5 * ir * ir;

    const float D1 = -mrinv3;
    const float D2 = 3.f * mrinv5;
    const float D3 = -15.f * mrinv7;

    const float q = multi.I_xx + multi.I_yy + multi.I_zz;
    const float qRx =
        multi.I_xx * dx[0] + multi.I_xy * dx[1] + multi.I_xz * dx[2];
    const float qRy =
        multi.I_xy * dx[0] + multi.I_yy * dx[1] + multi.I_yz * dx[2];
    const float qRz =
        multi.I_xz * dx[0] + multi.I_yz * dx[1] + multi.I_zz * dx[2];
    const float qRR = qRx * dx[0] + qRy * dx[1] + qRz * dx[2];
    const float C = D1 + 0.5f * D2 * q + 0.5f * D3 * qRR;

    gp->a_grav[0] -= C * dx[0] + D2 * qRx;
    gp->a_grav[1] -= C * dx[1] + D2 * qRy;
    gp->a_grav[2] -= C * dx[2] + D2 * qRz;

#else
#error "Multipoles of order >2 not yet implemented."
#endif
  }

  TIMER_TOC(TIMER_DOPAIR);  // MATTHIEU
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 *
 * @todo Use a local cache for the particles.
 */
__attribute__((always_inline)) INLINE static void runner_dograv_pair_pp(
    struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *e = r->e;
  const int gcount_i = ci->gcount;
  const int gcount_j = cj->gcount;
  struct gpart *restrict gparts_i = ci->gparts;
  struct gpart *restrict gparts_j = cj->gparts;
  const int ti_current = e->ti_current;

  TIMER_TIC;

#ifdef SANITY_CHECK
  if (ci->h[0] != cj->h[0])  // MATTHIEU sanity check
    error("Non matching cell sizes !! h_i=%f h_j=%f\n", ci->h[0], cj->h[0]);
#endif

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount_i; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gpi = &gparts_i[pid];
    const float mi = gpi->mass;

    /* Loop over every particle in the other cell. */
    for (int pjd = 0; pjd < gcount_j; pjd++) {

      /* Get a hold of the jth part in cj. */
      struct gpart *restrict gpj = &gparts_j[pjd];
      const float mj = gpj->mass;

      /* Compute the pairwise distance. */
      const float dx[3] = {gpi->x[0] - gpj->x[0],   // x
                           gpi->x[1] - gpj->x[1],   // y
                           gpi->x[2] - gpj->x[2]};  // z
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Apply the gravitational acceleration. */
      const float ir = 1.0f / sqrtf(r2);
      const float w = ir * ir * ir;
      const float wdx[3] = {w * dx[0], w * dx[1], w * dx[2]};

      if (gpi->ti_end <= ti_current) {
        gpi->a_grav[0] -= wdx[0] * mj;
        gpi->a_grav[1] -= wdx[1] * mj;
        gpi->a_grav[2] -= wdx[2] * mj;
      }
      if (gpj->ti_end <= ti_current) {
        gpj->a_grav[0] += wdx[0] * mi;
        gpj->a_grav[1] += wdx[1] * mi;
        gpj->a_grav[2] += wdx[2] * mi;
      }
    }
  }

  TIMER_TOC(TIMER_DOPAIR);  // MATTHIEU
}

/**
 * @brief Computes the interaction of all the particles in a cell directly
 *
 * @param r The #runner.
 * @param c The #cell.
 *
 * @todo Use a local cache for the particles.
 */
__attribute__((always_inline))
    INLINE static void runner_dograv_self_pp(struct runner *r, struct cell *c) {

  const struct engine *e = r->e;
  const int gcount = c->gcount;
  struct gpart *restrict gparts = c->gparts;
  const int ti_current = e->ti_current;

  TIMER_TIC;

#ifdef SANITY_CHECK
  if (c->gcount == 0)  // MATTHIEU sanity check
    error("Empty cell !");
#endif

  /* Anything to do here? */
  if (c->ti_end_min > ti_current) return;

  /* Loop over all particles in ci... */
  for (int pid = 0; pid < gcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct gpart *restrict gpi = &gparts[pid];
    const float mi = gpi->mass;

    /* Loop over every particle in the other cell. */
    for (int pjd = pid + 1; pjd < gcount; pjd++) {

      /* Get a hold of the jth part in ci. */
      struct gpart *restrict gpj = &gparts[pjd];
      const float mj = gpj->mass;

      /* Compute the pairwise distance. */
      const float dx[3] = {gpi->x[0] - gpj->x[0],   // x
                           gpi->x[1] - gpj->x[1],   // y
                           gpi->x[2] - gpj->x[2]};  // z
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      /* Apply the gravitational acceleration. */
      const float ir = 1.0f / sqrtf(r2);
      const float w = ir * ir * ir;
      const float wdx[3] = {w * dx[0], w * dx[1], w * dx[2]};

      if (gpi->ti_end <= ti_current) {
        gpi->a_grav[0] -= wdx[0] * mj;
        gpi->a_grav[1] -= wdx[1] * mj;
        gpi->a_grav[2] -= wdx[2] * mj;
      }
      if (gpj->ti_end <= ti_current) {
        gpj->a_grav[0] += wdx[0] * mi;
        gpj->a_grav[1] += wdx[1] * mi;
        gpj->a_grav[2] += wdx[2] * mi;
      }
    }
  }

  TIMER_TOC(TIMER_DOSELF);  // MATTHIEU
}

#endif /* SWIFT_RUNNER_DOIACT_GRAV_H */
