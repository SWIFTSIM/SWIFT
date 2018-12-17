/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOPAIR1_BRANCH(f) PASTE(runner_dopair1_branch, f)
#define DOPAIR1_BRANCH _DOPAIR1_BRANCH(FUNCTION)

#define _DOPAIR1(f) PASTE(runner_dopair1, f)
#define DOPAIR1 _DOPAIR1(FUNCTION)

#define _DOPAIR2_BRANCH(f) PASTE(runner_dopair2_branch, f)
#define DOPAIR2_BRANCH _DOPAIR2_BRANCH(FUNCTION)

#define _DOPAIR2(f) PASTE(runner_dopair2, f)
#define DOPAIR2 _DOPAIR2(FUNCTION)

#define _DOPAIR_SUBSET(f) PASTE(runner_dopair_subset, f)
#define DOPAIR_SUBSET _DOPAIR_SUBSET(FUNCTION)

#define _DOPAIR_SUBSET_BRANCH(f) PASTE(runner_dopair_subset_branch, f)
#define DOPAIR_SUBSET_BRANCH _DOPAIR_SUBSET_BRANCH(FUNCTION)

#define _DOPAIR_SUBSET_NOSORT(f) PASTE(runner_dopair_subset_nosort, f)
#define DOPAIR_SUBSET_NOSORT _DOPAIR_SUBSET_NOSORT(FUNCTION)

#define _DOPAIR_SUBSET_NAIVE(f) PASTE(runner_dopair_subset_naive, f)
#define DOPAIR_SUBSET_NAIVE _DOPAIR_SUBSET_NAIVE(FUNCTION)

#define _DOPAIR1_NAIVE(f) PASTE(runner_dopair1_naive, f)
#define DOPAIR1_NAIVE _DOPAIR1_NAIVE(FUNCTION)

#define _DOPAIR2_NAIVE(f) PASTE(runner_dopair2_naive, f)
#define DOPAIR2_NAIVE _DOPAIR2_NAIVE(FUNCTION)

#define _DOSELF1_NAIVE(f) PASTE(runner_doself1_naive, f)
#define DOSELF1_NAIVE _DOSELF1_NAIVE(FUNCTION)

#define _DOSELF2_NAIVE(f) PASTE(runner_doself2_naive, f)
#define DOSELF2_NAIVE _DOSELF2_NAIVE(FUNCTION)

#define _DOSELF1_BRANCH(f) PASTE(runner_doself1_branch, f)
#define DOSELF1_BRANCH _DOSELF1_BRANCH(FUNCTION)

#define _DOSELF1(f) PASTE(runner_doself1, f)
#define DOSELF1 _DOSELF1(FUNCTION)

#define _DOSELF2_BRANCH(f) PASTE(runner_doself2_branch, f)
#define DOSELF2_BRANCH _DOSELF2_BRANCH(FUNCTION)

#define _DOSELF2(f) PASTE(runner_doself2, f)
#define DOSELF2 _DOSELF2(FUNCTION)

#define _DOSELF_SUBSET(f) PASTE(runner_doself_subset, f)
#define DOSELF_SUBSET _DOSELF_SUBSET(FUNCTION)

#define _DOSELF_SUBSET_BRANCH(f) PASTE(runner_doself_subset_branch, f)
#define DOSELF_SUBSET_BRANCH _DOSELF_SUBSET_BRANCH(FUNCTION)

#define _DOSUB_SELF1(f) PASTE(runner_dosub_self1, f)
#define DOSUB_SELF1 _DOSUB_SELF1(FUNCTION)

#define _DOSUB_PAIR1(f) PASTE(runner_dosub_pair1, f)
#define DOSUB_PAIR1 _DOSUB_PAIR1(FUNCTION)

#define _DOSUB_SELF2(f) PASTE(runner_dosub_self2, f)
#define DOSUB_SELF2 _DOSUB_SELF2(FUNCTION)

#define _DOSUB_PAIR2(f) PASTE(runner_dosub_pair2, f)
#define DOSUB_PAIR2 _DOSUB_PAIR2(FUNCTION)

#define _DOSUB_SUBSET(f) PASTE(runner_dosub_subset, f)
#define DOSUB_SUBSET _DOSUB_SUBSET(FUNCTION)

#define _IACT_NONSYM(f) PASTE(runner_iact_nonsym, f)
#define IACT_NONSYM _IACT_NONSYM(FUNCTION)

#define _IACT(f) PASTE(runner_iact, f)
#define IACT _IACT(FUNCTION)

#define _IACT_NONSYM_VEC(f) PASTE(runner_iact_nonsym_vec, f)
#define IACT_NONSYM_VEC _IACT_NONSYM_VEC(FUNCTION)

#define _IACT_VEC(f) PASTE(runner_iact_vec, f)
#define IACT_VEC _IACT_VEC(FUNCTION)

#define _TIMER_DOSELF(f) PASTE(timer_doself, f)
#define TIMER_DOSELF _TIMER_DOSELF(FUNCTION)

#define _TIMER_DOPAIR(f) PASTE(timer_dopair, f)
#define TIMER_DOPAIR _TIMER_DOPAIR(FUNCTION)

#define _TIMER_DOSUB_SELF(f) PASTE(timer_dosub_self, f)
#define TIMER_DOSUB_SELF _TIMER_DOSUB_SELF(FUNCTION)

#define _TIMER_DOSUB_PAIR(f) PASTE(timer_dosub_pair, f)
#define TIMER_DOSUB_PAIR _TIMER_DOSUB_PAIR(FUNCTION)

#define _TIMER_DOSELF_SUBSET(f) PASTE(timer_doself_subset, f)
#define TIMER_DOSELF_SUBSET _TIMER_DOSELF_SUBSET(FUNCTION)

#define _TIMER_DOPAIR_SUBSET(f) PASTE(timer_dopair_subset, f)
#define TIMER_DOPAIR_SUBSET _TIMER_DOPAIR_SUBSET(FUNCTION)

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
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

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
    struct part *restrict pi = &parts_i[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = part_is_active(pi, e);
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
                          (float)(pi->x[1] - (cj->loc[1] + shift[1])),
                          (float)(pi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;
      const int pj_active = part_is_active(pj, e);

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
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
      }
      if (r2 < hjg2 && pj_active) {

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
      }

    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions between a cell pair (symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR2_NAIVE(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

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
    struct part *restrict pi = &parts_i[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = part_is_active(pi, e);
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
                          (float)(pi->x[1] - (cj->loc[1] + shift[1])),
                          (float)(pi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const int pj_active = part_is_active(pj, e);
      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

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
      if (r2 < hig2 || r2 < hjg2) {

        if (pi_active && pj_active) {

          IACT(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
        } else if (pi_active) {

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
        } else if (pj_active) {

          dx[0] = -dx[0];
          dx[1] = -dx[1];
          dx[2] = -dx[2];

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
        }
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
  if (!cell_is_active_hydro(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = part_is_active(pi, e);
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
                          (float)(pi->x[1] - c->loc[1]),
                          (float)(pi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;
      const int pj_active = part_is_active(pj, e);

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
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
      } else if (doi) {

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
      } else if (doj) {

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions within a cell (symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF2_NAIVE(struct runner *r, struct cell *restrict c) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const int pi_active = part_is_active(pi, e);
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
                          (float)(pi->x[1] - c->loc[1]),
                          (float)(pi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;
      const int pj_active = part_is_active(pj, e);

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      const int doi = pi_active && ((r2 < hig2) || (r2 < hjg2));
      const int doj = pj_active && ((r2 < hig2) || (r2 < hjg2));

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
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
      } else if (doi) {

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
      } else if (doj) {

        dx[0] = -dx[0];
        dx[1] = -dx[1];
        dx[2] = -dx[2];

        IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR_SUBSET_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct part *restrict parts_i, int *restrict ind,
                         int count, struct cell *restrict cj,
                         const double *shift) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over the parts_i. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[ind[pid]];
    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(pi, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pi->ti_drift != e->ti_current)
        error("Particle pi not drifted to current time");
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {

        IACT_NONSYM(r2, dx, hi, pj->h, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, pj->h, pi, pj, a, H);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(timer_dopair_subset_naive);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param flipped Flag to check whether the cells have been flipped or not.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR_SUBSET(struct runner *r, struct cell *restrict ci,
                   struct part *restrict parts_i, int *restrict ind, int count,
                   struct cell *restrict cj, const int sid, const int flipped,
                   const double *shift) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_j = cj->hydro.sort[sid];
  const float dxj = cj->hydro.dx_max_sort;

  /* Parts are on the left? */
  if (!flipped) {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - (shift[0]);
      const double piy = pi->x[1] - (shift[1]);
      const double piz = pi->x[2] - (shift[2]);
      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = hi * kernel_gamma + dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }

  /* Parts are on the right. */
  else {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      const double pix = pi->x[0] - (shift[0]);
      const double piy = pi->x[1] - (shift[1]);
      const double piz = pi->x[2] - (shift[2]);
      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const double di = -hi * kernel_gamma - dxj + pix * runner_shift[sid][0] +
                        piy * runner_shift[sid][1] + piz * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const double pjx = pj->x[0];
        const double pjy = pj->x[1];
        const double pjz = pj->x[2];

        /* Compute the pairwise distance. */
        float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                       (float)(piz - pjz)};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }

  TIMER_TOC(timer_dopair_subset);
}

/**
 * @brief Determine which version of DOPAIR_SUBSET needs to be called depending
 * on the
 * orientation of the cells or whether DOPAIR_SUBSET needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR_SUBSET_BRANCH(struct runner *r, struct cell *restrict ci,
                          struct part *restrict parts_i, int *restrict ind,
                          int count, struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

#if !defined(SWIFT_USE_NAIVE_INTERACTIONS)
  /* Get the sorting index. */
  int sid = 0;
  for (int k = 0; k < 3; k++)
    sid = 3 * sid + ((cj->loc[k] - ci->loc[k] + shift[k] < 0)
                         ? 0
                         : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2 : 1);

  /* Switch the cells around? */
  const int flipped = runner_flip[sid];
  sid = sortlistID[sid];

  /* Has the cell cj been sorted? */
  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");
#endif

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOPAIR_SUBSET_NAIVE(r, ci, parts_i, ind, count, cj, shift);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH)
  if (sort_is_face(sid))
    runner_dopair_subset_density_vec(r, ci, parts_i, ind, count, cj, sid,
                                     flipped, shift);
  else
    DOPAIR_SUBSET(r, ci, parts_i, ind, count, cj, sid, flipped, shift);
#else
  DOPAIR_SUBSET(r, ci, parts_i, ind, count, cj, sid, flipped, shift);
#endif
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 */
void DOSELF_SUBSET(struct runner *r, struct cell *restrict ci,
                   struct part *restrict parts, int *restrict ind, int count) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *pi = &parts[ind[pid]];
    const float pix[3] = {(float)(pi->x[0] - ci->loc[0]),
                          (float)(pi->x[1] - ci->loc[1]),
                          (float)(pi->x[2] - ci->loc[2])};
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!part_is_active(pi, e)) error("Inactive particle in subset function!");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_i; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - ci->loc[0]),
                            (float)(pj->x[1] - ci->loc[1]),
                            (float)(pj->x[2] - ci->loc[2])};
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
      if (r2 > 0.f && r2 < hig2) {

        IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
        runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */

  TIMER_TOC(timer_doself_subset);
}

/**
 * @brief Determine which version of DOSELF_SUBSET needs to be called depending
 * on the optimisation level.

 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 */
void DOSELF_SUBSET_BRANCH(struct runner *r, struct cell *restrict ci,
                          struct part *restrict parts, int *restrict ind,
                          int count) {

#if defined(WITH_VECTORIZATION) && defined(GADGET2_SPH)
  runner_doself_subset_density_vec(r, ci, parts, ind, count);
#else
  DOSELF_SUBSET(r, ci, parts, ind, count);
#endif
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
             const double *shift) {

  const struct engine *restrict e = r->e;
  const struct cosmology *restrict cosmo = e->cosmology;

  TIMER_TIC;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_i = ci->hydro.sort[sid];
  const struct entry *restrict sort_j = cj->hydro.sort[sid];

#ifdef SWIFT_DEBUG_CHECKS
  /* Some constants used to checks that the parts are in the right frame */
  const float shift_threshold_x =
      2. * ci->width[0] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_y =
      2. * ci->width[1] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_z =
      2. * ci->width[2] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get some other useful values. */
  const double hi_max = ci->hydro.h_max * kernel_gamma - rshift;
  const double hj_max = cj->hydro.h_max * kernel_gamma;
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;
  const float dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  if (cell_is_active_hydro(ci, e)) {

    /* Loop over the parts in ci. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[sort_i[pid].i];
      const float hi = pi->h;

      /* Skip inactive particles */
      if (!part_is_active(pi, e)) continue;

      /* Is there anything we need to interact with ? */
      const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
      if (di < dj_min) continue;

      /* Get some additional information about pi */
      const float hig2 = hi * hi * kernel_gamma2;
      const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
      const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
      const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;
        const float pjx = pj->x[0] - cj->loc[0];
        const float pjy = pj->x[1] - cj->loc[1];
        const float pjz = pj->x[2] - cj->loc[2];

        /* Compute the pairwise distance. */
        float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2) {

          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
        }
      } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
  }     /* Cell ci is active */

  if (cell_is_active_hydro(cj, e)) {

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth part in cj. */
      struct part *pj = &parts_j[sort_j[pjd].i];
      const float hj = pj->h;

      /* Skip inactive particles */
      if (!part_is_active(pj, e)) continue;

      /* Is there anything we need to interact with ? */
      const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max + rshift;
      if (dj - rshift > di_max) continue;

      /* Get some additional information about pj */
      const float hjg2 = hj * hj * kernel_gamma2;
      const float pjx = pj->x[0] - cj->loc[0];
      const float pjy = pj->x[1] - cj->loc[1];
      const float pjz = pj->x[2] - cj->loc[2];

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
        const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
        const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

        /* Compute the pairwise distance. */
        float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hjg2) {

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
        }
      } /* loop over the parts in ci. */
    }   /* loop over the parts in cj. */
  }     /* Cell cj is active */

  TIMER_TOC(TIMER_DOPAIR);
}

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
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells.");

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells.");
  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

#ifdef SWIFT_DEBUG_CHECKS
  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_i = ci->hydro.sort[sid];
  const struct entry *restrict sort_j = cj->hydro.sort[sid];

  /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
  for (int pid = 0; pid < ci->hydro.count; pid++) {
    const struct part *p = &ci->hydro.parts[sort_i[pid].i];
    if (part_is_inhibited(p, e)) continue;

    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_i[pid].d) - ci->hydro.dx_max_sort >
            1.0e-4 * max(fabsf(d), ci->hydro.dx_max_sort_old) &&
        fabsf(d - sort_i[pid].d) - ci->hydro.dx_max_sort >
            ci->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
          "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->hydro.dx_max_sort=%e "
          "ci->hydro.dx_max_sort_old=%e",
          ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->hydro.dx_max_sort,
          ci->hydro.dx_max_sort_old);
  }
  for (int pjd = 0; pjd < cj->hydro.count; pjd++) {
    const struct part *p = &cj->hydro.parts[sort_j[pjd].i];
    if (part_is_inhibited(p, e)) continue;

    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if ((fabsf(d - sort_j[pjd].d) - cj->hydro.dx_max_sort) >
            1.0e-4 * max(fabsf(d), cj->hydro.dx_max_sort_old) &&
        (fabsf(d - sort_j[pjd].d) - cj->hydro.dx_max_sort) >
            cj->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
          "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->hydro.dx_max_sort=%e "
          "cj->hydro.dx_max_sort_old=%e",
          cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->hydro.dx_max_sort,
          cj->hydro.dx_max_sort_old);
  }
#endif /* SWIFT_DEBUG_CHECKS */

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOPAIR1_NAIVE(r, ci, cj);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  if (!sort_is_corner(sid))
    runner_dopair1_density_vec(r, ci, cj, sid, shift);
  else
    DOPAIR1(r, ci, cj, sid, shift);
#else
  DOPAIR1(r, ci, cj, sid, shift);
#endif
}

/**
 * @brief Compute the interactions between a cell pair (symmetric)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR2(struct runner *r, struct cell *ci, struct cell *cj, const int sid,
             const double *shift) {

  const struct engine *restrict e = r->e;
  const struct cosmology *restrict cosmo = e->cosmology;

  TIMER_TIC;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  struct entry *restrict sort_i = ci->hydro.sort[sid];
  struct entry *restrict sort_j = cj->hydro.sort[sid];

#ifdef SWIFT_DEBUG_CHECKS
  /* Some constants used to checks that the parts are in the right frame */
  const float shift_threshold_x =
      2. * ci->width[0] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_y =
      2. * ci->width[1] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
  const float shift_threshold_z =
      2. * ci->width[2] +
      2. * max(ci->hydro.dx_max_part, cj->hydro.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get some other useful values. */
  const double hi_max = ci->hydro.h_max;
  const double hj_max = cj->hydro.h_max;
  const int count_i = ci->hydro.count;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_i = ci->hydro.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Maximal displacement since last rebuild */
  const double dx_max = (ci->hydro.dx_max_sort + cj->hydro.dx_max_sort);

  /* Position on the axis of the particles closest to the interface */
  const double di_max = sort_i[count_i - 1].d;
  const double dj_min = sort_j[0].d;

  /* Shifts to apply to the particles to be in a good frame */
  const double shift_i[3] = {cj->loc[0] + shift[0], cj->loc[1] + shift[1],
                             cj->loc[2] + shift[2]};
  const double shift_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

  int count_active_i = 0, count_active_j = 0;
  struct entry *restrict sort_active_i = NULL, *restrict sort_active_j = NULL;

  if (cell_is_all_active_hydro(ci, e)) {
    /* If everybody is active don't bother copying */
    sort_active_i = sort_i;
    count_active_i = count_i;
  } else if (cell_is_active_hydro(ci, e)) {
    if (posix_memalign((void **)&sort_active_i, SWIFT_CACHE_ALIGNMENT,
                       sizeof(struct entry) * count_i) != 0)
      error("Failed to allocate active sortlists.");

    /* Collect the active particles in ci */
    for (int k = 0; k < count_i; k++) {
      if (part_is_active(&parts_i[sort_i[k].i], e)) {
        sort_active_i[count_active_i] = sort_i[k];
        count_active_i++;
      }
    }
  }

  if (cell_is_all_active_hydro(cj, e)) {
    /* If everybody is active don't bother copying */
    sort_active_j = sort_j;
    count_active_j = count_j;
  } else if (cell_is_active_hydro(cj, e)) {
    if (posix_memalign((void **)&sort_active_j, SWIFT_CACHE_ALIGNMENT,
                       sizeof(struct entry) * count_j) != 0)
      error("Failed to allocate active sortlists.");

    /* Collect the active particles in cj */
    for (int k = 0; k < count_j; k++) {
      if (part_is_active(&parts_j[sort_j[k].i], e)) {
        sort_active_j[count_active_j] = sort_j[k];
        count_active_j++;
      }
    }
  }

  /* Loop over *all* the parts in ci starting from the centre until
     we are out of range of anything in cj (using the maximal hi). */
  for (int pid = count_i - 1;
       pid >= 0 &&
       sort_i[pid].d + hi_max * kernel_gamma + dx_max - rshift > dj_min;
       pid--) {

    /* Get a hold of the ith part in ci. */
    struct part *pi = &parts_i[sort_i[pid].i];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    const float hi = pi->h;

    /* Is there anything we need to interact with (for this specific hi) ? */
    const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
    if (di < dj_min) continue;

    /* Get some additional information about pi */
    const float hig2 = hi * hi * kernel_gamma2;
    const float pix = pi->x[0] - shift_i[0];
    const float piy = pi->x[1] - shift_i[1];
    const float piz = pi->x[2] - shift_i[2];

    /* Do we need to only check active parts in cj
       (i.e. pi does not need updating) ? */
    if (!part_is_active(pi, e)) {

      /* Loop over the *active* parts in cj within range of pi */
      for (int pjd = 0; pjd < count_active_j && sort_active_j[pjd].d < di;
           pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_active_j[pjd].i];
        const float hj = pj->h;

        /* Get the position of pj in the right frame */
        const float pjx = pj->x[0] - shift_j[0];
        const float pjy = pj->x[1] - shift_j[1];
        const float pjz = pj->x[2] - shift_j[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss?
           (note that we will do the other condition in the reverse loop) */
        if (r2 < hig2) {
          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
        }
      } /* loop over the active parts in cj. */
    }

    else { /* pi is active, we may need to update pi and pj */

      /* Loop over *all* the parts in cj in range of pi. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Recover pj */
        struct part *pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

        /* Get the position of pj in the right frame */
        const float pjx = pj->x[0] - shift_j[0];
        const float pjy = pj->x[1] - shift_j[1];
        const float pjz = pj->x[2] - shift_j[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif
        /* Hit or miss?
           (note that we will do the other condition in the reverse loop) */
        if (r2 < hig2) {

          /* Does pj need to be updated too? */
          if (part_is_active(pj, e)) {
            IACT(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
          } else {
            IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
          }
        }
      } /* loop over the parts in cj. */
    }   /* Is pi active? */
  }     /* Loop over all ci */

  /* Loop over *all* the parts in cj starting from the centre until
     we are out of range of anything in ci (using the maximal hj). */
  for (int pjd = 0;
       pjd < count_j &&
       sort_j[pjd].d - hj_max * kernel_gamma - dx_max < di_max - rshift;
       pjd++) {

    /* Get a hold of the jth part in cj. */
    struct part *pj = &parts_j[sort_j[pjd].i];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pj, e)) continue;

    const float hj = pj->h;

    /* Is there anything we need to interact with (for this specific hj) ? */
    const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max;
    if (dj > di_max - rshift) continue;

    /* Get some additional information about pj */
    const float hjg2 = hj * hj * kernel_gamma2;
    const float pjx = pj->x[0] - shift_j[0];
    const float pjy = pj->x[1] - shift_j[1];
    const float pjz = pj->x[2] - shift_j[2];

    /* Do we need to only check active parts in ci
       (i.e. pj does not need updating) ? */
    if (!part_is_active(pj, e)) {

      /* Loop over the *active* parts in ci. */
      for (int pid = count_active_i - 1;
           pid >= 0 && sort_active_i[pid].d - rshift > dj; pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_active_i[pid].i];
        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;

        /* Get the position of pi in the right frame */
        const float pix = pi->x[0] - shift_i[0];
        const float piy = pi->x[1] - shift_i[1];
        const float piz = pi->x[2] - shift_i[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss?
           (note that we must avoid the r2 < hig2 cases we already processed) */
        if (r2 < hjg2 && r2 >= hig2) {
          IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
        }
      } /* loop over the active parts in ci. */
    }

    else { /* pj is active, we may need to update pj and pi */

      /* Loop over *all* the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d - rshift > dj;
           pid--) {

        /* Recover pi */
        struct part *pi = &parts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pi, e)) continue;

        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;

        /* Get the position of pi in the right frame */
        const float pix = pi->x[0] - shift_i[0];
        const float piy = pi->x[1] - shift_i[1];
        const float piz = pi->x[2] - shift_i[2];

        /* Compute the pairwise distance. */
        const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles are in the correct frame after the shifts */
        if (pix > shift_threshold_x || pix < -shift_threshold_x)
          error(
              "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
              pix, ci->width[0]);
        if (piy > shift_threshold_y || piy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
              piy, ci->width[1]);
        if (piz > shift_threshold_z || piz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
              piz, ci->width[2]);
        if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
          error(
              "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
              pjx, ci->width[0]);
        if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
          error(
              "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
              pjy, ci->width[1]);
        if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
          error(
              "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
              pjz, ci->width[2]);

        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss?
           (note that we must avoid the r2 < hig2 cases we already processed) */
        if (r2 < hjg2 && r2 >= hig2) {

          /* Does pi need to be updated too? */
          if (part_is_active(pi, e)) {
            IACT(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
          } else {
            IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
          }
        }
      } /* loop over the parts in ci. */
    }   /* Is pj active? */
  }     /* Loop over all cj */

  /* Clean-up if necessary */
  if (cell_is_active_hydro(ci, e) && !cell_is_all_active_hydro(ci, e))
    free(sort_active_i);
  if (cell_is_active_hydro(cj, e) && !cell_is_all_active_hydro(cj, e))
    free(sort_active_j);

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Determine which version of DOPAIR2 needs to be called depending on the
 * orientation of the cells or whether DOPAIR2 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR2_BRANCH(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells.");

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->hydro.sorted & (1 << sid)) ||
      ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells.");
  if (!(cj->hydro.sorted & (1 << sid)) ||
      cj->hydro.dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

#ifdef SWIFT_DEBUG_CHECKS
  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_i = ci->hydro.sort[sid];
  const struct entry *restrict sort_j = cj->hydro.sort[sid];

  /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
  for (int pid = 0; pid < ci->hydro.count; pid++) {
    const struct part *p = &ci->hydro.parts[sort_i[pid].i];
    if (part_is_inhibited(p, e)) continue;

    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_i[pid].d) - ci->hydro.dx_max_sort >
            1.0e-4 * max(fabsf(d), ci->hydro.dx_max_sort_old) &&
        fabsf(d - sort_i[pid].d) - ci->hydro.dx_max_sort >
            ci->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
          "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->hydro.dx_max_sort=%e "
          "ci->hydro.dx_max_sort_old=%e",
          ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->hydro.dx_max_sort,
          ci->hydro.dx_max_sort_old);
  }
  for (int pjd = 0; pjd < cj->hydro.count; pjd++) {
    const struct part *p = &cj->hydro.parts[sort_j[pjd].i];
    if (part_is_inhibited(p, e)) continue;

    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_j[pjd].d) - cj->hydro.dx_max_sort >
            1.0e-4 * max(fabsf(d), cj->hydro.dx_max_sort_old) &&
        fabsf(d - sort_j[pjd].d) - cj->hydro.dx_max_sort >
            cj->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
          "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->hydro.dx_max_sort=%e "
          "cj->hydro.dx_max_sort_old=%e",
          cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->hydro.dx_max_sort,
          cj->hydro.dx_max_sort_old);
  }
#endif /* SWIFT_DEBUG_CHECKS */

#ifdef SWIFT_USE_NAIVE_INTERACTIONS
  DOPAIR2_NAIVE(r, ci, cj);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  if (!sort_is_corner(sid))
    runner_dopair2_force_vec(r, ci, cj, sid, shift);
  else
    DOPAIR2(r, ci, cj, sid, shift);
#else
  DOPAIR2(r, ci, cj, sid, shift);
#endif
}

/**
 * @brief Compute the cell self-interaction (non-symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1(struct runner *r, struct cell *restrict c) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  struct part *restrict parts = c->hydro.parts;
  const int count = c->hydro.count;

  /* Set up indt. */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < count; k++)
    if (part_is_active(&parts[k], e)) {
      indt[countdt] = k;
      countdt += 1;
    }

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    /* Get the particle position and radius. */
    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k];
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle inactive? */
    if (!part_is_active(pi, e)) {

      /* Loop over the other particles .*/
      for (int pjd = firstdt; pjd < countdt; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[indt[pjd]];
        const float hj = pj->h;

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pj->x[k] - pix[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hj * hj * kernel_gamma2) {

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
        }
      } /* loop over all other particles. */
    }

    /* Otherwise, interact with all candidates. */
    else {

      /* We caught a live one! */
      firstdt += 1;

      /* Loop over the other particles .*/
      for (int pjd = pid + 1; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[pjd];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }
        const int doj =
            (part_is_active(pj, e)) && (r2 < hj * hj * kernel_gamma2);

        const int doi = (r2 < hig2);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (doi || doj) {

          /* Which parts need to be updated? */
          if (doi && doj) {

            IACT(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
          } else if (doi) {

            IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
          } else if (doj) {

            dx[0] = -dx[0];
            dx[1] = -dx[1];
            dx[2] = -dx[2];
            IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
          }
        }
      } /* loop over all other particles. */
    }
  } /* loop over all particles. */

  free(indt);

  TIMER_TOC(TIMER_DOSELF);
}

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
  if (!cell_is_active_hydro(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->hydro.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF1_NAIVE(r, c);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  runner_doself1_density_vec(r, c);
#else
  DOSELF1(r, c);
#endif
}

/**
 * @brief Compute the cell self-interaction (symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF2(struct runner *r, struct cell *restrict c) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  struct part *restrict parts = c->hydro.parts;
  const int count = c->hydro.count;

  /* Set up indt. */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < count; k++)
    if (part_is_active(&parts[k], e)) {
      indt[countdt] = k;
      countdt += 1;
    }

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts[pid];

    /* Skip inhibited particles. */
    if (part_is_inhibited(pi, e)) continue;

    /* Get the particle position and radius. */
    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k];
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle not active? */
    if (!part_is_active(pi, e)) {

      /* Loop over the other particles .*/
      for (int pjd = firstdt; pjd < countdt; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[indt[pjd]];
        const float hj = pj->h;

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pj->x[k] - pix[k];
          r2 += dx[k] * dx[k];
        }

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2 || r2 < hj * hj * kernel_gamma2) {

          IACT_NONSYM(r2, dx, hj, hi, pj, pi, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
          runner_iact_nonsym_chemistry(r2, dx, hj, hi, pj, pi, a, H);
#endif
        }
      } /* loop over all other particles. */
    }

    /* Otherwise, interact with all candidates. */
    else {

      /* We caught a live one! */
      firstdt += 1;

      /* Loop over the other particles .*/
      for (int pjd = pid + 1; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[pjd];

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2 || r2 < hj * hj * kernel_gamma2) {

          /* Does pj need to be updated too? */
          if (part_is_active(pj, e)) {
            IACT(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
          } else {
            IACT_NONSYM(r2, dx, hi, hj, pi, pj, a, H);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
            runner_iact_nonsym_chemistry(r2, dx, hi, hj, pi, pj, a, H);
#endif
          }
        }
      } /* loop over all other particles. */
    }
  } /* loop over all particles. */

  free(indt);

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Determine which version of DOSELF2 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF2_BRANCH(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->hydro.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
  DOSELF2_NAIVE(r, c);
#elif defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (FUNCTION_TASK_LOOP == TASK_LOOP_FORCE)
  runner_doself2_force_vec(r, c);
#else
  DOSELF2(r, c);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction linking the cells
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj, int sid,
                 int gettimer) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  TIMER_TIC;

  /* Should we even bother? */
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Get the type of pair if not specified explicitly. */
  double shift[3];
  sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_hydro_task(ci) &&
      cell_can_recurse_in_pair_hydro_task(cj)) {

    /* Different types of flags. */
    switch (sid) {

      /* Regular sub-cell interactions of a single cell. */
      case 0: /* (  1 ,  1 ,  1 ) */
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[0], -1, 0);
        break;

      case 1: /* (  1 ,  1 ,  0 ) */
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[1], -1, 0);
        break;

      case 2: /* (  1 ,  1 , -1 ) */
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[1], -1, 0);
        break;

      case 3: /* (  1 ,  0 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[2], -1, 0);
        break;

      case 4: /* (  1 ,  0 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[0], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[1], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[2], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[1], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[3], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[2], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[3], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[3], -1, 0);
        break;

      case 5: /* (  1 ,  0 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[1], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[3], -1, 0);
        break;

      case 6: /* (  1 , -1 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[2], -1, 0);
        break;

      case 7: /* (  1 , -1 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[2], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[3], -1, 0);
        break;

      case 8: /* (  1 , -1 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR1(r, ci->progeny[4], cj->progeny[3], -1, 0);
        break;

      case 9: /* (  0 ,  1 ,  1 ) */
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[4], -1, 0);
        break;

      case 10: /* (  0 ,  1 ,  0 ) */
        if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[2], cj->progeny[0], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[2], cj->progeny[1], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[2], cj->progeny[4], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1(r, ci->progeny[2], cj->progeny[5], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[1], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[5], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[4], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[5], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[5], -1, 0);
        break;

      case 11: /* (  0 ,  1 , -1 ) */
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[2], cj->progeny[1], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1(r, ci->progeny[2], cj->progeny[5], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR1(r, ci->progeny[6], cj->progeny[5], -1, 0);
        break;

      case 12: /* (  0 ,  0 ,  1 ) */
        if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[1], cj->progeny[0], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[1], cj->progeny[2], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[1], cj->progeny[4], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1(r, ci->progeny[1], cj->progeny[6], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[2], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1(r, ci->progeny[3], cj->progeny[6], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[4], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1(r, ci->progeny[5], cj->progeny[6], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR1(r, ci->progeny[7], cj->progeny[6], -1, 0);
        break;
    }

  }

  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e)) {

    /* Make sure both cells are drifted to the current timestep. */
    if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
      error("Interacting undrifted cells.");

    /* Do any of the cells need to be sorted first? */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx)
      error(
          "Interacting unsorted cell. ci->hydro.dx_max_sort_old=%e ci->dmin=%e "
          "ci->sorted=%d sid=%d",
          ci->hydro.dx_max_sort_old, ci->dmin, ci->hydro.sorted, sid);
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx)
      error(
          "Interacting unsorted cell. cj->hydro.dx_max_sort_old=%e cj->dmin=%e "
          "cj->sorted=%d sid=%d",
          cj->hydro.dx_max_sort_old, cj->dmin, cj->hydro.sorted, sid);

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
  if (ci->hydro.count == 0 || !cell_is_active_hydro(ci, r->e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_hydro_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1(r, ci->progeny[k], ci->progeny[j], -1, 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Drift the cell to the current timestep if needed. */
    if (!cell_are_part_drifted(ci, r->e)) error("Interacting undrifted cell.");

    DOSELF1_BRANCH(r, ci);
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}

/**
 * @brief Compute grouped sub-cell interactions for pairs (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction linking the cells
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR2(struct runner *r, struct cell *ci, struct cell *cj, int sid,
                 int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Should we even bother? */
  if (!cell_is_active_hydro(ci, e) && !cell_is_active_hydro(cj, e)) return;
  if (ci->hydro.count == 0 || cj->hydro.count == 0) return;

  /* Get the type of pair if not specified explicitly. */
  double shift[3];
  sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_hydro_task(ci) &&
      cell_can_recurse_in_pair_hydro_task(cj)) {

    /* Different types of flags. */
    switch (sid) {

      /* Regular sub-cell interactions of a single cell. */
      case 0: /* (  1 ,  1 ,  1 ) */
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[0], -1, 0);
        break;

      case 1: /* (  1 ,  1 ,  0 ) */
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[1], -1, 0);
        break;

      case 2: /* (  1 ,  1 , -1 ) */
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[1], -1, 0);
        break;

      case 3: /* (  1 ,  0 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[2], -1, 0);
        break;

      case 4: /* (  1 ,  0 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[0], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[1], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[2], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[1], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[3], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[2], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[3], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[3], -1, 0);
        break;

      case 5: /* (  1 ,  0 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[1], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[3], -1, 0);
        break;

      case 6: /* (  1 , -1 ,  1 ) */
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[2], -1, 0);
        break;

      case 7: /* (  1 , -1 ,  0 ) */
        if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[2], -1, 0);
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[3], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[3], -1, 0);
        break;

      case 8: /* (  1 , -1 , -1 ) */
        if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
          DOSUB_PAIR2(r, ci->progeny[4], cj->progeny[3], -1, 0);
        break;

      case 9: /* (  0 ,  1 ,  1 ) */
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[4], -1, 0);
        break;

      case 10: /* (  0 ,  1 ,  0 ) */
        if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[2], cj->progeny[0], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[2], cj->progeny[1], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[2], cj->progeny[4], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR2(r, ci->progeny[2], cj->progeny[5], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[1], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[5], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[0], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[4], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[5], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[1], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[5], -1, 0);
        break;

      case 11: /* (  0 ,  1 , -1 ) */
        if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[2], cj->progeny[1], -1, 0);
        if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR2(r, ci->progeny[2], cj->progeny[5], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[1], -1, 0);
        if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
          DOSUB_PAIR2(r, ci->progeny[6], cj->progeny[5], -1, 0);
        break;

      case 12: /* (  0 ,  0 ,  1 ) */
        if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[1], cj->progeny[0], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[1], cj->progeny[2], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[1], cj->progeny[4], -1, 0);
        if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR2(r, ci->progeny[1], cj->progeny[6], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[0], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[2], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[4], -1, 0);
        if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR2(r, ci->progeny[3], cj->progeny[6], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[0], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[2], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[4], -1, 0);
        if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR2(r, ci->progeny[5], cj->progeny[6], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[0], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[2], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[4], -1, 0);
        if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
          DOSUB_PAIR2(r, ci->progeny[7], cj->progeny[6], -1, 0);
        break;
    }

  }

  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e)) {

    /* Make sure both cells are drifted to the current timestep. */
    if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
      error("Interacting undrifted cells.");

    /* Do any of the cells need to be sorted first? */
    if (!(ci->hydro.sorted & (1 << sid)) ||
        ci->hydro.dx_max_sort_old > ci->dmin * space_maxreldx)
      error(
          "Interacting unsorted cell. ci->hydro.dx_max_sort_old=%e ci->dmin=%e "
          "ci->sorted=%d sid=%d",
          ci->hydro.dx_max_sort_old, ci->dmin, ci->hydro.sorted, sid);
    if (!(cj->hydro.sorted & (1 << sid)) ||
        cj->hydro.dx_max_sort_old > cj->dmin * space_maxreldx)
      error(
          "Interacting unsorted cell. cj->hydro.dx_max_sort_old=%e cj->dmin=%e "
          "cj->sorted=%d sid=%d",
          cj->hydro.dx_max_sort_old, cj->dmin, cj->hydro.sorted, sid);

    /* Compute the interactions. */
    DOPAIR2_BRANCH(r, ci, cj);
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF2(struct runner *r, struct cell *ci, int gettimer) {

  TIMER_TIC;

  /* Should we even bother? */
  if (ci->hydro.count == 0 || !cell_is_active_hydro(ci, r->e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_hydro_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF2(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR2(r, ci->progeny[k], ci->progeny[j], -1, 0);
      }

  }

  /* Otherwise, compute self-interaction. */
  else {
    DOSELF2_BRANCH(r, ci);
  }
  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}

void DOSUB_SUBSET(struct runner *r, struct cell *ci, struct part *parts,
                  int *ind, int count, struct cell *cj, int sid, int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Should we even bother? */
  if (!cell_is_active_hydro(ci, e) &&
      (cj == NULL || !cell_is_active_hydro(cj, e)))
    return;
  if (ci->hydro.count == 0 || (cj != NULL && cj->hydro.count == 0)) return;

  /* Find out in which sub-cell of ci the parts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&parts[ind[0]] >= &ci->progeny[k]->hydro.parts[0] &&
            &parts[ind[0]] <
                &ci->progeny[k]->hydro.parts[ci->progeny[k]->hydro.count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_hydro_task(ci)) {

      /* Loop over all progeny. */
      DOSUB_SUBSET(r, sub, parts, ind, count, NULL, -1, 0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          DOSUB_SUBSET(r, sub, parts, ind, count, ci->progeny[j], -1, 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF_SUBSET_BRANCH(r, ci, parts, ind, count);
  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_hydro_task(ci) &&
        cell_can_recurse_in_pair_hydro_task(cj)) {

      /* Get the type of pair if not specified explicitly. */
      double shift[3] = {0.0, 0.0, 0.0};
      sid = space_getsid(s, &ci, &cj, shift);

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[7],
                         -1, 0);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[7],
                         -1, 0);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[6],
                         -1, 0);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[7],
                         -1, 0);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[7],
                         -1, 0);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[6],
                         -1, 0);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[5],
                         -1, 0);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[4],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[5],
                         -1, 0);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] == sub && cj->progeny[3] != NULL)
            DOSUB_SUBSET(r, ci->progeny[4], parts, ind, count, cj->progeny[3],
                         -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] == sub)
            DOSUB_SUBSET(r, cj->progeny[3], parts, ind, count, ci->progeny[4],
                         -1, 0);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[7],
                         -1, 0);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[2], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[2],
                         -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[2], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[2],
                         -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[2], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[2],
                         -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET(r, ci->progeny[2], parts, ind, count, cj->progeny[5],
                         -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET(r, cj->progeny[5], parts, ind, count, ci->progeny[2],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[5],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET(r, cj->progeny[5], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[5],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET(r, cj->progeny[5], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[5],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET(r, cj->progeny[5], parts, ind, count, ci->progeny[7],
                         -1, 0);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[2], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[2],
                         -1, 0);
          if (ci->progeny[2] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET(r, ci->progeny[2], parts, ind, count, cj->progeny[5],
                         -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET(r, cj->progeny[5], parts, ind, count, ci->progeny[2],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[1] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[1],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] == sub)
            DOSUB_SUBSET(r, cj->progeny[1], parts, ind, count, ci->progeny[6],
                         -1, 0);
          if (ci->progeny[6] == sub && cj->progeny[5] != NULL)
            DOSUB_SUBSET(r, ci->progeny[6], parts, ind, count, cj->progeny[5],
                         -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] == sub)
            DOSUB_SUBSET(r, cj->progeny[5], parts, ind, count, ci->progeny[6],
                         -1, 0);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[1], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[1],
                         -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[1], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[1],
                         -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[1], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[1],
                         -1, 0);
          if (ci->progeny[1] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET(r, ci->progeny[1], parts, ind, count, cj->progeny[6],
                         -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET(r, cj->progeny[6], parts, ind, count, ci->progeny[1],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[3] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET(r, ci->progeny[3], parts, ind, count, cj->progeny[6],
                         -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET(r, cj->progeny[6], parts, ind, count, ci->progeny[3],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[5] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET(r, ci->progeny[5], parts, ind, count, cj->progeny[6],
                         -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET(r, cj->progeny[6], parts, ind, count, ci->progeny[5],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, cj->progeny[0], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[2] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[2],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] == sub)
            DOSUB_SUBSET(r, cj->progeny[2], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[4] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[4],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] == sub)
            DOSUB_SUBSET(r, cj->progeny[4], parts, ind, count, ci->progeny[7],
                         -1, 0);
          if (ci->progeny[7] == sub && cj->progeny[6] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[6],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[6] == sub)
            DOSUB_SUBSET(r, cj->progeny[6], parts, ind, count, ci->progeny[7],
                         -1, 0);
          break;
      }

    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_hydro(ci, e) || cell_is_active_hydro(cj, e)) {

      /* Do any of the cells need to be drifted first? */
      if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");

      DOPAIR_SUBSET_BRANCH(r, ci, parts, ind, count, cj);
    }

  } /* otherwise, pair interaction. */

  if (gettimer) TIMER_TOC(timer_dosub_subset);
}
