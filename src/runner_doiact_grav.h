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
 * @brief Compute the sorted gravity interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */

void runner_dopair_grav_new(struct runner *r, struct cell *ci,
                            struct cell *cj) {

  struct engine *restrict e = r->e;
  int pid, pjd, k, sid;
  double rshift, shift[3] = {0.0, 0.0, 0.0}, nshift[3];
  struct entry *restrict sort_i, *restrict sort_j;
  struct gpart *restrict pi, *restrict pj, *restrict parts_i, *restrict parts_j;
  double pix[3];
  float dx[3], r2, h_max, di, dj;
  int count_i, count_j, cnj, cnj_new;
  const int ti_current = e->ti_current;
  struct multipole m;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct gpart *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

  /* Get the sort ID. */
  sid = space_getsid(e->s, &ci, &cj, shift);

  /* Make sure the cells are sorted. */
  runner_dogsort(r, ci, (1 << sid), 0);
  runner_dogsort(r, cj, (1 << sid), 0);

  /* Have the cells been sorted? */
  if (!(ci->gsorted & (1 << sid)) || !(cj->gsorted & (1 << sid)))
    error("Trying to interact unsorted cells.");

  /* Get the cutoff shift. */
  for (rshift = 0.0, k = 0; k < 3; k++)
    rshift += shift[k] * runner_shift[3 * sid + k];

  /* Pick-out the sorted lists. */
  sort_i = &ci->gsort[sid * (ci->count + 1)];
  sort_j = &cj->gsort[sid * (cj->count + 1)];

  /* Get some other useful values. */
  h_max =
      sqrtf(ci->h[0] * ci->h[0] + ci->h[1] * ci->h[1] + ci->h[2] * ci->h[2]) *
      const_theta_max;
  count_i = ci->gcount;
  count_j = cj->gcount;
  parts_i = ci->gparts;
  parts_j = cj->gparts;
  cnj = count_j;
  multipole_reset(&m);
  nshift[0] = -shift[0];
  nshift[1] = -shift[1];
  nshift[2] = -shift[2];

  /* Loop over the parts in ci. */
  for (pid = count_i - 1; pid >= 0; pid--) {

    /* Get a hold of the ith part in ci. */
    pi = &parts_i[sort_i[pid].i];
    if (pi->ti_end > ti_current) continue;
    di = sort_i[pid].d + h_max - rshift;

    for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];

    /* Loop over the parts in cj. */
    for (pjd = 0; pjd < cnj && sort_j[pjd].d < di; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts_j[sort_j[pjd].i];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

#ifndef VECTORIZE

      // if ( pi->part->id == 3473472412525 || pj->part->id == 3473472412525 )
      //     message( "interacting particles pi=%lli and pj=%lli with r=%.3e in
      // cells %lli/%lli." , pi->part->id , pj->part->id , sqrtf(r2) , ((long
      // long int)ci) / sizeof(struct cell) , ((long long int)cj) /
      // sizeof(struct cell) );

      runner_iact_grav(r2, dx, pi, pj);

#else

      /* Add this interaction to the queue. */
      r2q[icount] = r2;
      dxq[3 * icount + 0] = dx[0];
      dxq[3 * icount + 1] = dx[1];
      dxq[3 * icount + 2] = dx[2];
      piq[icount] = pi;
      pjq[icount] = pj;
      icount += 1;

      /* Flush? */
      if (icount == VEC_SIZE) {
        runner_iact_vec_grav(r2q, dxq, piq, pjq);
        icount = 0;
      }

#endif

    } /* loop over the parts in cj. */

    /* Set the new limit. */
    cnj_new = pjd;

    /* Add trailing parts to the multipole. */
    for (pjd = cnj_new; pjd < cnj; pjd++) {

      /* Add the part to the multipole. */
      multipole_addpart(&m, &parts_j[sort_j[pjd].i]);

    } /* add trailing parts to the multipole. */

    /* Set the new cnj. */
    cnj = cnj_new;

    /* Interact the ith particle with the multipole. */
    multipole_iact_mp(&m, pi, nshift);

  } /* loop over the parts in ci. */

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      runner_iact_grav(r2q[k], &dxq[3 * k], piq[k], pjq[k]);
#endif

  /* Re-set the multipole. */
  multipole_reset(&m);

  /* Loop over the parts in cj and interact with the multipole in ci. */
  for (pid = count_i - 1, pjd = 0; pjd < count_j; pjd++) {

    /* Get the position of pj along the axis. */
    dj = sort_j[pjd].d - h_max + rshift;

    /* Add any left-over parts in cell_i to the multipole. */
    while (pid >= 0 && sort_i[pid].d < dj) {

      /* Add this particle to the multipole. */
      multipole_addpart(&m, &parts_i[sort_i[pid].i]);

      /* Decrease pid. */
      pid -= 1;
    }

    /* Interact pj with the multipole. */
    multipole_iact_mp(&m, &parts_j[sort_j[pjd].i], shift);

  } /* loop over the parts in cj and interact with the multipole. */

  TIMER_TOC(TIMER_DOPAIR);
}

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
 * @brief Computes the interaction of all the particles in a cell with the
 * multipole of another cell.
 *
 * @param r The #runner.
 * @param ci The #cell with particles to interct.
 * @param cj The #cell with the multipole.
 */
__attribute__((always_inline)) INLINE static void runner_dograv_pair_pm(
    struct runner *r, struct cell *ci, struct cell *cj) {

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
 * @param ci The first #cell.
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
