/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "runner_doiact_dark_matter_limiter.h"

/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR1_NAIVE(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_dark_matter(ci, e) && !cell_is_starting_dark_matter(cj, e)) return;

  const int count_i = ci->dark_matter.count;
  const int count_j = cj->dark_matter.count;
  struct dmpart *restrict dmparts_i = ci->dark_matter.parts;
  struct dmpart *restrict dmparts_j = cj->dark_matter.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count_i; pid++) {

    /* Get a hold of the ith part in ci. */
    struct dmpart *restrict pi = &dmparts_i[pid];

    /* Skip inhibited particles. */
    if (dmpart_is_inhibited(pi, e)) continue;

    const int pi_active = dmpart_is_starting(pi, e);
    const float hi = pi->h;
    const float hig2 = hi * hi * dm_kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
                          (float)(pi->x[1] - (cj->loc[1] + shift[1])),
                          (float)(pi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct dmpart *restrict pj = &dmparts_j[pjd];

      /* Skip inhibited particles. */
      if (dmpart_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * dm_kernel_gamma2;
      const int pj_active = dmpart_is_starting(pj, e);

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2 && pi_active) {

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
      }
      if (r2 < hjg2 && pj_active) {

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
      }

    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1_NAIVE(struct runner *r, struct cell *restrict c) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_dark_matter(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count = c->dark_matter.count;
  struct dmpart *restrict dmparts = c->dark_matter.parts;

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct dmpart *restrict pi = &dmparts[pid];

    /* Skip inhibited particles. */
    if (dmpart_is_inhibited(pi, e)) continue;

    const int pi_active = dmpart_is_starting(pi, e);
    const float hi = pi->h;
    const float hig2 = hi * hi * dm_kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
                          (float)(pi->x[1] - c->loc[1]),
                          (float)(pi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct dmpart *restrict pj = &dmparts[pjd];

      /* Skip inhibited particles. */
      if (dmpart_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * dm_kernel_gamma2;
      const int pj_active = dmpart_is_starting(pj, e);

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int doi = pi_active && (r2 < hig2);
      const int doj = pj_active && (r2 < hjg2);

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (doi && doj) {

        IACT(r2, dx, hi, hj, pi, pj, a, H);
      } else if (doi) {

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
      } else if (doj) {

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions between a cell pair (non-symmetric).
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1(struct runner *r, struct cell *ci, struct cell *cj, const int sid,
             const double *shift) {}

/**
 * @brief Determine which version of DOPAIR1 needs to be called depending on the
 * orientation of the cells or whether DOPAIR1 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (ci->dark_matter.count == 0 || cj->dark_matter.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_starting_dark_matter(ci, e) && !cell_is_starting_dark_matter(cj, e)) return;

  /* Check that cells are drifted. */
  if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
    error("Interacting undrifted cells.");

  DOPAIR1_NAIVE(r, ci, cj);

/*#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOPAIR1_NAIVE(r, ci, cj);
#else
  DOPAIR1(r, ci, cj, sid, shift);
#endif*/
}

/**
 * @brief Compute the cell self-interaction (non-symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1(struct runner *r, struct cell *restrict c) {}

/**
 * @brief Determine which version of DOSELF1 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->dark_matter.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_starting_dark_matter(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->dark_matter.h_max_old * dm_kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  /* Check that cells are drifted. */
  if (!cell_are_dmpart_drifted(c, e)) error("Interacting undrifted cell.");

  DOSELF1_NAIVE(r, c);
/*#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF1_NAIVE(r, c);
#else
  DOSELF1(r, c);
#endif*/
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj,
                 int gettimer) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  TIMER_TIC;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  /* Should we even bother? */
  const int do_i = cell_get_flag(ci, cell_flag_do_dark_matter_limiter);
  const int do_j = cell_get_flag(cj, cell_flag_do_dark_matter_limiter);
  const int do_sub_i = cell_get_flag(ci, cell_flag_do_dark_matter_sub_limiter);
  const int do_sub_j = cell_get_flag(cj, cell_flag_do_dark_matter_sub_limiter);

  if (!do_i && !do_j && !do_sub_i && !do_sub_j) return;
  if (!cell_is_starting_dark_matter(ci, e) && !cell_is_starting_dark_matter(cj, e)) return;
  if (ci->dark_matter.count == 0 || cj->dark_matter.count == 0) return;

  /* Recurse? */
  if (cell_can_recurse_in_pair_dark_matter_task(ci) &&
      cell_can_recurse_in_pair_dark_matter_task(cj)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else if ((cell_is_starting_dark_matter(ci, e) && (do_i || do_sub_i)) ||
           (cell_is_starting_dark_matter(cj, e) && (do_j || do_sub_j))) {

    /* Make sure both cells are drifted to the current timestep. */
    if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
      error("Interacting undrifted cells.");

    /* Compute the interactions. */
    DOPAIR1_BRANCH(r, ci, cj);
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1(struct runner *r, struct cell *ci, int gettimer) {

  TIMER_TIC;

  /* Should we even bother? */
  const int do_i = cell_get_flag(ci, cell_flag_do_dark_matter_limiter);
  const int do_sub_i = cell_get_flag(ci, cell_flag_do_dark_matter_sub_limiter);

  if (!do_i && !do_sub_i) return;
  if (!cell_is_starting_dark_matter(ci, r->e)) return;
  if (ci->dark_matter.count == 0) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_dark_matter_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1(r, ci->progeny[k], ci->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Drift the cell to the current timestep if needed. */
    if (!cell_are_dmpart_drifted(ci, r->e)) error("Interacting undrifted cell.");

    DOSELF1_BRANCH(r, ci);
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}
