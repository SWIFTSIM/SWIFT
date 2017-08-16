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

#define _DOPAIR2(f) PASTE(runner_dopair2, f)
#define DOPAIR2 _DOPAIR2(FUNCTION)

#define _DOPAIR_SUBSET(f) PASTE(runner_dopair_subset, f)
#define DOPAIR_SUBSET _DOPAIR_SUBSET(FUNCTION)

#define _DOPAIR_SUBSET_NOSORT(f) PASTE(runner_dopair_subset_nosort, f)
#define DOPAIR_SUBSET_NOSORT _DOPAIR_SUBSET_NOSORT(FUNCTION)

#define _DOPAIR_SUBSET_NAIVE(f) PASTE(runner_dopair_subset_naive, f)
#define DOPAIR_SUBSET_NAIVE _DOPAIR_SUBSET_NAIVE(FUNCTION)

#define _DOPAIR1_NAIVE(f) PASTE(runner_dopair1_naive, f)
#define DOPAIR1_NAIVE _DOPAIR1_NAIVE(FUNCTION)

#define _DOPAIR2_NAIVE(f) PASTE(runner_dopair2_naive, f)
#define DOPAIR2_NAIVE _DOPAIR2_NAIVE(FUNCTION)

#define _DOSELF2_NAIVE(f) PASTE(runner_doself2_naive, f)
#define DOSELF2_NAIVE _DOSELF2_NAIVE(FUNCTION)

#define _DOSELF1(f) PASTE(runner_doself1, f)
#define DOSELF1 _DOSELF1(FUNCTION)

#define _DOSELF2(f) PASTE(runner_doself2, f)
#define DOSELF2 _DOSELF2(FUNCTION)

#define _DOSELF_SUBSET(f) PASTE(runner_doself_subset, f)
#define DOSELF_SUBSET _DOSELF_SUBSET(FUNCTION)

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

#ifndef SWIFT_DEBUG_CHECKS
  error("Don't use in actual runs ! Slow code !");
#endif

#ifdef WITH_VECTORIZATION
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  const int count_i = ci->count;
  const int count_j = cj->count;
  struct part *restrict parts_i = ci->parts;
  struct part *restrict parts_j = cj->parts;

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
    const float hi = pi->h;

    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    const float hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {

#ifndef WITH_VECTORIZATION

        IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);

#else

        /* Add this interaction to the queue. */
        r2q[icount] = r2;
        dxq[3 * icount + 0] = dx[0];
        dxq[3 * icount + 1] = dx[1];
        dxq[3 * icount + 2] = dx[2];
        hiq[icount] = hi;
        hjq[icount] = pj->h;
        piq[icount] = pi;
        pjq[icount] = pj;
        icount += 1;

        /* Flush? */
        if (icount == VEC_SIZE) {
          IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
          icount = 0;
        }

#endif
      }
      if (r2 < pj->h * pj->h * kernel_gamma2) {

#ifndef WITH_VECTORIZATION

        for (int k = 0; k < 3; k++) dx[k] = -dx[k];
        IACT_NONSYM(r2, dx, pj->h, hi, pj, pi);

#else

        /* Add this interaction to the queue. */
        r2q[icount] = r2;
        dxq[3 * icount + 0] = -dx[0];
        dxq[3 * icount + 1] = -dx[1];
        dxq[3 * icount + 2] = -dx[2];
        hiq[icount] = pj->h;
        hjq[icount] = hi;
        piq[icount] = pj;
        pjq[icount] = pi;
        icount += 1;

        /* Flush? */
        if (icount == VEC_SIZE) {
          IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
          icount = 0;
        }

#endif
      }

    } /* loop over the parts in cj. */

  } /* loop over the parts in ci. */

#ifdef WITH_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount > 0)
    for (int k = 0; k < icount; k++)
      IACT_NONSYM(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

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

#ifndef SWIFT_DEBUG_CHECKS
  error("Don't use in actual runs ! Slow code !");
#endif

#ifdef WITH_OLD_VECTORIZATION
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  const int count_i = ci->count;
  const int count_j = cj->count;
  struct part *restrict parts_i = ci->parts;
  struct part *restrict parts_j = cj->parts;

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
    const float hi = pi->h;

    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    const float hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < pj->h * pj->h * kernel_gamma2) {

#ifndef WITH_OLD_VECTORIZATION

        IACT(r2, dx, hi, pj->h, pi, pj);

#else

        /* Add this interaction to the queue. */
        r2q[icount] = r2;
        dxq[3 * icount + 0] = dx[0];
        dxq[3 * icount + 1] = dx[1];
        dxq[3 * icount + 2] = dx[2];
        hiq[icount] = hi;
        hjq[icount] = pj->h;
        piq[icount] = pi;
        pjq[icount] = pj;
        icount += 1;

        /* Flush? */
        if (icount == VEC_SIZE) {
          IACT_VEC(r2q, dxq, hiq, hjq, piq, pjq);
          icount = 0;
        }

#endif
      }

    } /* loop over the parts in cj. */

  } /* loop over the parts in ci. */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount > 0)
    for (int k = 0; k < icount; k++)
      IACT(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

  TIMER_TOC(TIMER_DOPAIR);
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

#ifndef SWIFT_DEBUG_CHECKS
  error("Don't use in actual runs ! Slow code !");
#endif

#ifdef WITH_OLD_VECTORIZATION
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(c, e)) return;

  const int count = c->count;
  struct part *restrict parts = c->parts;

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts[pid];
    const double pix[3] = {pi->x[0], pi->x[1], pi->x[2]};
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (int pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts[pjd];

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < pj->h * pj->h * kernel_gamma2) {

#ifndef WITH_OLD_VECTORIZATION

        IACT(r2, dx, hi, pj->h, pi, pj);

#else

        /* Add this interaction to the queue. */
        r2q[icount] = r2;
        dxq[3 * icount + 0] = dx[0];
        dxq[3 * icount + 1] = dx[1];
        dxq[3 * icount + 2] = dx[2];
        hiq[icount] = hi;
        hjq[icount] = pj->h;
        piq[icount] = pi;
        pjq[icount] = pj;
        icount += 1;

        /* Flush? */
        if (icount == VEC_SIZE) {
          IACT_VEC(r2q, dxq, hiq, hjq, piq, pjq);
          icount = 0;
        }

#endif
      }

    } /* loop over the parts in cj. */

  } /* loop over the parts in ci. */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount > 0)
    for (int k = 0; k < icount; k++)
      IACT(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

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
 */
void DOPAIR_SUBSET_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct part *restrict parts_i, int *restrict ind,
                         int count, struct cell *restrict cj) {

  const struct engine *e = r->e;

#ifdef WITH_OLD_VECTORIZATION
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif

  TIMER_TIC;

  const int count_j = cj->count;
  struct part *restrict parts_j = cj->parts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

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

#ifndef WITH_OLD_VECTORIZATION

        IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);

#else

        /* Add this interaction to the queue. */
        r2q[icount] = r2;
        dxq[3 * icount + 0] = dx[0];
        dxq[3 * icount + 1] = dx[1];
        dxq[3 * icount + 2] = dx[2];
        hiq[icount] = hi;
        hjq[icount] = pj->h;
        piq[icount] = pi;
        pjq[icount] = pj;
        icount += 1;

        /* Flush? */
        if (icount == VEC_SIZE) {
          IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
          icount = 0;
        }

#endif
      }

    } /* loop over the parts in cj. */

  } /* loop over the parts in ci. */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount > 0)
    for (int k = 0; k < icount; k++)
      IACT_NONSYM(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

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
 */
void DOPAIR_SUBSET(struct runner *r, struct cell *restrict ci,
                   struct part *restrict parts_i, int *restrict ind, int count,
                   struct cell *restrict cj) {

  struct engine *e = r->e;

#ifdef WITH_OLD_VECTORIZATION
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif

  TIMER_TIC;

  const int count_j = cj->count;
  struct part *restrict parts_j = cj->parts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

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
  if (!(cj->sorted & (1 << sid)) ||
      cj->dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_j = cj->sort[sid];
  const float dxj = cj->dx_max_sort;

  /* Parts are on the left? */
  if (!flipped) {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      double pix[3];
      for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];

      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float di = hi * kernel_gamma + dxj + pix[0] * runner_shift[sid][0] +
                       pix[1] * runner_shift[sid][1] +
                       pix[2] * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2) {

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);

#else

          /* Add this interaction to the queue. */
          r2q[icount] = r2;
          dxq[3 * icount + 0] = dx[0];
          dxq[3 * icount + 1] = dx[1];
          dxq[3 * icount + 2] = dx[2];
          hiq[icount] = hi;
          hjq[icount] = pj->h;
          piq[icount] = pi;
          pjq[icount] = pj;
          icount += 1;

          /* Flush? */
          if (icount == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
            icount = 0;
          }

#endif
        }

      } /* loop over the parts in cj. */

    } /* loop over the parts in ci. */

  }

  /* Parts are on the right. */
  else {

    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[ind[pid]];
      double pix[3];
      for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
      const float hi = pi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float di =
          -hi * kernel_gamma - dxj + pix[0] * runner_shift[sid][0] +
          pix[1] * runner_shift[sid][1] + pix[2] * runner_shift[sid][2];

      /* Loop over the parts in cj. */
      for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2) {

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);

#else

          /* Add this interaction to the queue. */
          r2q[icount] = r2;
          dxq[3 * icount + 0] = dx[0];
          dxq[3 * icount + 1] = dx[1];
          dxq[3 * icount + 2] = dx[2];
          hiq[icount] = hi;
          hjq[icount] = pj->h;
          piq[icount] = pi;
          pjq[icount] = pj;
          icount += 1;

          /* Flush? */
          if (icount == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
            icount = 0;
          }

#endif
        }

      } /* loop over the parts in cj. */

    } /* loop over the parts in ci. */
  }

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount > 0)
    for (int k = 0; k < icount; k++)
      IACT_NONSYM(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

  TIMER_TOC(timer_dopair_subset);
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

#ifdef WITH_OLD_VECTORIZATION
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif

  TIMER_TIC;

  const int count_i = ci->count;
  struct part *restrict parts_j = ci->parts;

  /* Loop over the parts in ci. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts[ind[pid]];
    const double pix[3] = {pi->x[0], pi->x[1], pi->x[2]};
    const float hi = pi->h;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_i; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      float r2 = 0.0f;
      float dx[3];
      for (int k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 > 0.0f && r2 < hig2) {

#ifndef WITH_OLD_VECTORIZATION

        IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);

#else

        /* Add this interaction to the queue. */
        r2q[icount] = r2;
        dxq[3 * icount + 0] = dx[0];
        dxq[3 * icount + 1] = dx[1];
        dxq[3 * icount + 2] = dx[2];
        hiq[icount] = hi;
        hjq[icount] = pj->h;
        piq[icount] = pi;
        pjq[icount] = pj;
        icount += 1;

        /* Flush? */
        if (icount == VEC_SIZE) {
          IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
          icount = 0;
        }

#endif
      }

    } /* loop over the parts in cj. */

  } /* loop over the parts in ci. */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount > 0)
    for (int k = 0; k < icount; k++)
      IACT_NONSYM(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

  TIMER_TOC(timer_doself_subset);
}

/**
 * @brief Compute the interactions between a cell pair (non-symmetric).
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1(struct runner *r, struct cell *ci, struct cell *cj, const int sid,
             const double *shift) {

  const struct engine *restrict e = r->e;

#ifdef WITH_OLD_VECTORIZATION
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif

  TIMER_TIC;

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  const struct entry *restrict sort_i = ci->sort[sid];
  const struct entry *restrict sort_j = cj->sort[sid];

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
  for (int pid = 0; pid < ci->count; pid++) {
    const struct part *p = &ci->parts[sort_i[pid].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_i[pid].d) - ci->dx_max_sort >
        1.0e-4 * max(fabsf(d), ci->dx_max_sort_old))
      error(
          "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
          "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->dx_max_sort=%e "
          "ci->dx_max_sort_old=%e",
          ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->dx_max_sort,
          ci->dx_max_sort_old);
  }
  for (int pjd = 0; pjd < cj->count; pjd++) {
    const struct part *p = &cj->parts[sort_j[pjd].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_j[pjd].d) - cj->dx_max_sort >
        1.0e-4 * max(fabsf(d), cj->dx_max_sort_old))
      error(
          "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
          "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->dx_max_sort=%e "
          "cj->dx_max_sort_old=%e",
          cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->dx_max_sort,
          cj->dx_max_sort_old);
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get some other useful values. */
  const double hi_max = ci->h_max * kernel_gamma - rshift;
  const double hj_max = cj->h_max * kernel_gamma;
  const int count_i = ci->count;
  const int count_j = cj->count;
  struct part *restrict parts_i = ci->parts;
  struct part *restrict parts_j = cj->parts;
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;
  const float dx_max = (ci->dx_max_sort + cj->dx_max_sort);

  if (cell_is_active(ci, e)) {

    /* Loop over the parts in ci. */
    for (int pid = count_i - 1;
         pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

      /* Get a hold of the ith part in ci. */
      struct part *restrict pi = &parts_i[sort_i[pid].i];
      if (!part_is_active(pi, e)) continue;
      const float hi = pi->h;
      const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
      if (di < dj_min) continue;

      double pix[3];
      for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
      const float hig2 = hi * hi * kernel_gamma2;

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];

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

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hi, pj->h, pi, pj);

#else

          /* Add this interaction to the queue. */
          r2q[icount] = r2;
          dxq[3 * icount + 0] = dx[0];
          dxq[3 * icount + 1] = dx[1];
          dxq[3 * icount + 2] = dx[2];
          hiq[icount] = hi;
          hjq[icount] = pj->h;
          piq[icount] = pi;
          pjq[icount] = pj;
          icount += 1;

          /* Flush? */
          if (icount == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
            icount = 0;
          }

#endif
        }

      } /* loop over the parts in cj. */

    } /* loop over the parts in ci. */

  } /* Cell ci is active */

  if (cell_is_active(cj, e)) {

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
         pjd++) {

      /* Get a hold of the jth part in cj. */
      struct part *restrict pj = &parts_j[sort_j[pjd].i];
      if (!part_is_active(pj, e)) continue;
      const float hj = pj->h;
      const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max + rshift;
      if (dj - rshift > di_max) continue;
      
      double pjx[3];
      for (int k = 0; k < 3; k++) pjx[k] = pj->x[k] + shift[k];
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pi = &parts_i[sort_i[pid].i];

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pjx[k] - pi->x[k];
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
        if (r2 < hjg2) {

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hj, pi->h, pj, pi);

#else

          /* Add this interaction to the queue. */
          r2q[icount] = r2;
          dxq[3 * icount + 0] = dx[0];
          dxq[3 * icount + 1] = dx[1];
          dxq[3 * icount + 2] = dx[2];
          hiq[icount] = hj;
          hjq[icount] = pi->h;
          piq[icount] = pj;
          pjq[icount] = pi;
          icount += 1;

          /* Flush? */
          if (icount == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q, dxq, hiq, hjq, piq, pjq);
            icount = 0;
          }

#endif
        }

      } /* loop over the parts in ci. */

    } /* loop over the parts in cj. */

  } /* Cell cj is active */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount > 0)
    for (int k = 0; k < icount; k++)
      IACT_NONSYM(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

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
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells.");

  /* Get the sort ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->sorted & (1 << sid)) ||
      ci->dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells.");
  if (!(cj->sorted & (1 << sid)) ||
      cj->dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

#if defined(WITH_VECTORIZATION) && defined(GADGET2_SPH) && \
    (DOPAIR1_BRANCH == runner_dopair1_density_branch)
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
 */
void DOPAIR2(struct runner *r, struct cell *ci, struct cell *cj) {

  struct engine *restrict e = r->e;

#ifdef WITH_OLD_VECTORIZATION
  int icount1 = 0;
  float r2q1[VEC_SIZE] __attribute__((aligned(16)));
  float hiq1[VEC_SIZE] __attribute__((aligned(16)));
  float hjq1[VEC_SIZE] __attribute__((aligned(16)));
  float dxq1[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq1[VEC_SIZE], *pjq1[VEC_SIZE];
  int icount2 = 0;
  float r2q2[VEC_SIZE] __attribute__((aligned(16)));
  float hiq2[VEC_SIZE] __attribute__((aligned(16)));
  float hjq2[VEC_SIZE] __attribute__((aligned(16)));
  float dxq2[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq2[VEC_SIZE], *pjq2[VEC_SIZE];
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;

  if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
    error("Interacting undrifted cells.");

  /* Get the shift ID. */
  double shift[3] = {0.0, 0.0, 0.0};
  const int sid = space_getsid(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->sorted & (1 << sid)) ||
      ci->dx_max_sort_old > space_maxreldx * ci->dmin)
    error("Interacting unsorted cells.");
  if (!(cj->sorted & (1 << sid)) ||
      cj->dx_max_sort_old > space_maxreldx * cj->dmin)
    error("Interacting unsorted cells.");

  /* Get the cutoff shift. */
  double rshift = 0.0;
  for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

  /* Pick-out the sorted lists. */
  struct entry *restrict sort_i = ci->sort[sid];
  struct entry *restrict sort_j = cj->sort[sid];

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
  for (int pid = 0; pid < ci->count; pid++) {
    const struct part *p = &ci->parts[sort_i[pid].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_i[pid].d) - ci->dx_max_sort >
        1.0e-4 * max(fabsf(d), ci->dx_max_sort_old))
      error(
          "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
          "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->dx_max_sort=%e "
          "ci->dx_max_sort_old=%e",
          ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->dx_max_sort,
          ci->dx_max_sort_old);
  }
  for (int pjd = 0; pjd < cj->count; pjd++) {
    const struct part *p = &cj->parts[sort_j[pjd].i];
    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_j[pjd].d) - cj->dx_max_sort >
        1.0e-4 * max(fabsf(d), cj->dx_max_sort_old))
      error(
          "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
          "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->dx_max_sort=%e "
          "cj->dx_max_sort_old=%e",
          cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->dx_max_sort,
          cj->dx_max_sort_old);
  }
#endif /* SWIFT_DEBUG_CHECKS */

  /* Get some other useful values. */
  const double hi_max = ci->h_max * kernel_gamma - rshift;
  const double hj_max = cj->h_max * kernel_gamma;
  const int count_i = ci->count;
  const int count_j = cj->count;
  struct part *restrict parts_i = ci->parts;
  struct part *restrict parts_j = cj->parts;
  const double di_max = sort_i[count_i - 1].d - rshift;
  const double dj_min = sort_j[0].d;
  const double dx_max = (ci->dx_max_sort + cj->dx_max_sort);

  /* Collect the number of parts left and right below dt. */
  int countdt_i = 0, countdt_j = 0;
  struct entry *restrict sortdt_i = NULL, *restrict sortdt_j = NULL;
  if (cell_is_all_active(ci, e)) {
    sortdt_i = sort_i;
    countdt_i = count_i;
  } else if (cell_is_active(ci, e)) {
    if (posix_memalign((void *)&sortdt_i, VEC_SIZE * sizeof(float),
                       sizeof(struct entry) * count_i) != 0)
      error("Failed to allocate dt sortlists.");
    for (int k = 0; k < count_i; k++)
      if (part_is_active(&parts_i[sort_i[k].i], e)) {
        sortdt_i[countdt_i] = sort_i[k];
        countdt_i += 1;
      }
  }
  if (cell_is_all_active(cj, e)) {
    sortdt_j = sort_j;
    countdt_j = count_j;
  } else if (cell_is_active(cj, e)) {
    if (posix_memalign((void *)&sortdt_j, VEC_SIZE * sizeof(float),
                       sizeof(struct entry) * count_j) != 0)
      error("Failed to allocate dt sortlists.");
    for (int k = 0; k < count_j; k++)
      if (part_is_active(&parts_j[sort_j[k].i], e)) {
        sortdt_j[countdt_j] = sort_j[k];
        countdt_j += 1;
      }
  }

  /* Loop over the parts in ci. */
  for (int pid = count_i - 1;
       pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

    /* Get a hold of the ith part in ci. */
    struct part *restrict pi = &parts_i[sort_i[pid].i];
    const float hi = pi->h;
    const double di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
    if (di < dj_min) continue;

    double pix[3];
    for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    const float hig2 = hi * hi * kernel_gamma2;

    /* Look at valid dt parts only? */
    if (!part_is_active(pi, e)) {

      /* Loop over the parts in cj within dt. */
      for (int pjd = 0; pjd < countdt_j && sortdt_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sortdt_j[pjd].i];
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
        if (r2 < hig2) {

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hj, hi, pj, pi);

#else

          /* Add this interaction to the queue. */
          r2q1[icount1] = r2;
          dxq1[3 * icount1 + 0] = dx[0];
          dxq1[3 * icount1 + 1] = dx[1];
          dxq1[3 * icount1 + 2] = dx[2];
          hiq1[icount1] = hj;
          hjq1[icount1] = hi;
          piq1[icount1] = pj;
          pjq1[icount1] = pi;
          icount1 += 1;

          /* Flush? */
          if (icount1 == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
            icount1 = 0;
          }

#endif
        }

      } /* loop over the parts in cj. */

    }

    /* Otherwise, look at all parts. */
    else {

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[sort_j[pjd].i];
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
        if (r2 < hig2) {

#ifndef WITH_OLD_VECTORIZATION

          /* Does pj need to be updated too? */
          if (part_is_active(pj, e))
            IACT(r2, dx, hi, hj, pi, pj);
          else
            IACT_NONSYM(r2, dx, hi, hj, pi, pj);

#else

          /* Does pj need to be updated too? */
          if (part_is_active(pj, e)) {

            /* Add this interaction to the symmetric queue. */
            r2q2[icount2] = r2;
            dxq2[3 * icount2 + 0] = dx[0];
            dxq2[3 * icount2 + 1] = dx[1];
            dxq2[3 * icount2 + 2] = dx[2];
            hiq2[icount2] = hi;
            hjq2[icount2] = hj;
            piq2[icount2] = pi;
            pjq2[icount2] = pj;
            icount2 += 1;

            /* Flush? */
            if (icount2 == VEC_SIZE) {
              IACT_VEC(r2q2, dxq2, hiq2, hjq2, piq2, pjq2);
              icount2 = 0;
            }

          } else {

            /* Add this interaction to the non-symmetric queue. */
            r2q1[icount1] = r2;
            dxq1[3 * icount1 + 0] = dx[0];
            dxq1[3 * icount1 + 1] = dx[1];
            dxq1[3 * icount1 + 2] = dx[2];
            hiq1[icount1] = hi;
            hjq1[icount1] = hj;
            piq1[icount1] = pi;
            pjq1[icount1] = pj;
            icount1 += 1;

            /* Flush? */
            if (icount1 == VEC_SIZE) {
              IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
              icount1 = 0;
            }
          }

#endif
        }

      } /* loop over the parts in cj. */
    }

  } /* loop over the parts in ci. */

  /* Loop over the parts in cj. */
  for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
       pjd++) {

    /* Get a hold of the jth part in cj. */
    struct part *restrict pj = &parts_j[sort_j[pjd].i];
    const float hj = pj->h;
    const double dj = sort_j[pjd].d - hj * kernel_gamma - dx_max + rshift;
    if (dj - rshift > di_max) continue;
    
    double pjx[3];
    for (int k = 0; k < 3; k++) pjx[k] = pj->x[k] + shift[k];
    const float hjg2 = hj * hj * kernel_gamma2;

    /* Is this particle outside the dt? */
    if (!part_is_active(pj, e)) {

      /* Loop over the parts in ci. */
      for (int pid = countdt_i - 1; pid >= 0 && sortdt_i[pid].d > dj; pid--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pi = &parts_i[sortdt_i[pid].i];
        const float hi = pi->h;

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pi->x[k] - pjx[k];
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
        if (r2 < hjg2 && r2 > hi * hi * kernel_gamma2) {

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hi, hj, pi, pj);

#else

          /* Add this interaction to the queue. */
          r2q1[icount1] = r2;
          dxq1[3 * icount1 + 0] = dx[0];
          dxq1[3 * icount1 + 1] = dx[1];
          dxq1[3 * icount1 + 2] = dx[2];
          hiq1[icount1] = hi;
          hjq1[icount1] = hj;
          piq1[icount1] = pi;
          pjq1[icount1] = pj;
          icount1 += 1;

          /* Flush? */
          if (icount1 == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
            icount1 = 0;
          }

#endif
        }

      } /* loop over the parts in cj. */
    }

    /* Otherwise, interact with all particles in cj. */
    else {

      /* Loop over the parts in ci. */
      for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pi = &parts_i[sort_i[pid].i];
        const float hi = pi->h;

        /* Compute the pairwise distance. */
        float r2 = 0.0f;
        float dx[3];
        for (int k = 0; k < 3; k++) {
          dx[k] = pjx[k] - pi->x[k];
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
        if (r2 < hjg2 && r2 > hi * hi * kernel_gamma2) {

#ifndef WITH_OLD_VECTORIZATION

          /* Does pi need to be updated too? */
          if (part_is_active(pi, e))
            IACT(r2, dx, hj, hi, pj, pi);
          else
            IACT_NONSYM(r2, dx, hj, hi, pj, pi);

#else

          /* Does pi need to be updated too? */
          if (part_is_active(pi, e)) {

            /* Add this interaction to the symmetric queue. */
            r2q2[icount2] = r2;
            dxq2[3 * icount2 + 0] = dx[0];
            dxq2[3 * icount2 + 1] = dx[1];
            dxq2[3 * icount2 + 2] = dx[2];
            hiq2[icount2] = hj;
            hjq2[icount2] = hi;
            piq2[icount2] = pj;
            pjq2[icount2] = pi;
            icount2 += 1;

            /* Flush? */
            if (icount2 == VEC_SIZE) {
              IACT_VEC(r2q2, dxq2, hiq2, hjq2, piq2, pjq2);
              icount2 = 0;
            }

          } else {

            /* Add this interaction to the non-symmetric queue. */
            r2q1[icount1] = r2;
            dxq1[3 * icount1 + 0] = dx[0];
            dxq1[3 * icount1 + 1] = dx[1];
            dxq1[3 * icount1 + 2] = dx[2];
            hiq1[icount1] = hj;
            hjq1[icount1] = hi;
            piq1[icount1] = pj;
            pjq1[icount1] = pi;
            icount1 += 1;

            /* Flush? */
            if (icount1 == VEC_SIZE) {
              IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
              icount1 = 0;
            }
          }

#endif
        }

      } /* loop over the parts in cj. */
    }

  } /* loop over the parts in ci. */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount1 > 0)
    for (int k = 0; k < icount1; k++)
      IACT_NONSYM(r2q1[k], &dxq1[3 * k], hiq1[k], hjq1[k], piq1[k], pjq1[k]);
  if (icount2 > 0)
    for (int k = 0; k < icount2; k++)
      IACT(r2q2[k], &dxq2[3 * k], hiq2[k], hjq2[k], piq2[k], pjq2[k]);
#endif

  /* Clean-up if necessary */
  if (cell_is_active(ci, e) && !cell_is_all_active(ci, e)) free(sortdt_i);
  if (cell_is_active(cj, e) && !cell_is_all_active(cj, e)) free(sortdt_j);

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the cell self-interaction (non-symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1(struct runner *r, struct cell *restrict c) {

  const struct engine *e = r->e;

#ifdef WITH_OLD_VECTORIZATION
  int icount1 = 0;
  float r2q1[VEC_SIZE] __attribute__((aligned(16)));
  float hiq1[VEC_SIZE] __attribute__((aligned(16)));
  float hjq1[VEC_SIZE] __attribute__((aligned(16)));
  float dxq1[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq1[VEC_SIZE], *pjq1[VEC_SIZE];
  int icount2 = 0;
  float r2q2[VEC_SIZE] __attribute__((aligned(16)));
  float hiq2[VEC_SIZE] __attribute__((aligned(16)));
  float hjq2[VEC_SIZE] __attribute__((aligned(16)));
  float dxq2[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq2[VEC_SIZE], *pjq2[VEC_SIZE];
#endif

  TIMER_TIC;

  if (!cell_is_active(c, e)) return;

  if (!cell_are_part_drifted(c, e)) error("Interacting undrifted cell.");

  struct part *restrict parts = c->parts;
  const int count = c->count;

  /* Set up indt. */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void *)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < count; k++)
    if (part_is_active(&parts[k], e)) {
      indt[countdt] = k;
      countdt += 1;
    }

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts[pid];

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

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hj, hi, pj, pi);

#else

          /* Add this interaction to the queue. */
          r2q1[icount1] = r2;
          dxq1[3 * icount1 + 0] = dx[0];
          dxq1[3 * icount1 + 1] = dx[1];
          dxq1[3 * icount1 + 2] = dx[2];
          hiq1[icount1] = hj;
          hjq1[icount1] = hi;
          piq1[icount1] = pj;
          pjq1[icount1] = pi;
          icount1 += 1;

          /* Flush? */
          if (icount1 == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
            icount1 = 0;
          }

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

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pi->ti_drift != e->ti_current)
          error("Particle pi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        /* Hit or miss? */
        if (r2 < hig2 || doj) {

#ifndef WITH_OLD_VECTORIZATION

          /* Which parts need to be updated? */
          if (r2 < hig2 && doj)
            IACT(r2, dx, hi, hj, pi, pj);
          else if (!doj)
            IACT_NONSYM(r2, dx, hi, hj, pi, pj);
          else {
            dx[0] = -dx[0];
            dx[1] = -dx[1];
            dx[2] = -dx[2];
            IACT_NONSYM(r2, dx, hj, hi, pj, pi);
          }

#else

          /* Does pj need to be updated too? */
          if (r2 < hig2 && doj) {

            /* Add this interaction to the symmetric queue. */
            r2q2[icount2] = r2;
            dxq2[3 * icount2 + 0] = dx[0];
            dxq2[3 * icount2 + 1] = dx[1];
            dxq2[3 * icount2 + 2] = dx[2];
            hiq2[icount2] = hi;
            hjq2[icount2] = hj;
            piq2[icount2] = pi;
            pjq2[icount2] = pj;
            icount2 += 1;

            /* Flush? */
            if (icount2 == VEC_SIZE) {
              IACT_VEC(r2q2, dxq2, hiq2, hjq2, piq2, pjq2);
              icount2 = 0;
            }

          } else if (!doj) {

            /* Add this interaction to the non-symmetric queue. */
            r2q1[icount1] = r2;
            dxq1[3 * icount1 + 0] = dx[0];
            dxq1[3 * icount1 + 1] = dx[1];
            dxq1[3 * icount1 + 2] = dx[2];
            hiq1[icount1] = hi;
            hjq1[icount1] = hj;
            piq1[icount1] = pi;
            pjq1[icount1] = pj;
            icount1 += 1;

            /* Flush? */
            if (icount1 == VEC_SIZE) {
              IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
              icount1 = 0;
            }

          } else {

            /* Add this interaction to the non-symmetric queue. */
            r2q1[icount1] = r2;
            dxq1[3 * icount1 + 0] = -dx[0];
            dxq1[3 * icount1 + 1] = -dx[1];
            dxq1[3 * icount1 + 2] = -dx[2];
            hiq1[icount1] = hj;
            hjq1[icount1] = hi;
            piq1[icount1] = pj;
            pjq1[icount1] = pi;
            icount1 += 1;

            /* Flush? */
            if (icount1 == VEC_SIZE) {
              IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
              icount1 = 0;
            }
          }

#endif
        }

      } /* loop over all other particles. */
    }

  } /* loop over all particles. */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount1 > 0)
    for (int k = 0; k < icount1; k++)
      IACT_NONSYM(r2q1[k], &dxq1[3 * k], hiq1[k], hjq1[k], piq1[k], pjq1[k]);
  if (icount2 > 0)
    for (int k = 0; k < icount2; k++)
      IACT(r2q2[k], &dxq2[3 * k], hiq2[k], hjq2[k], piq2[k], pjq2[k]);
#endif

  free(indt);

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the cell self-interaction (symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF2(struct runner *r, struct cell *restrict c) {

  const struct engine *e = r->e;

#ifdef WITH_OLD_VECTORIZATION
  int icount1 = 0;
  float r2q1[VEC_SIZE] __attribute__((aligned(16)));
  float hiq1[VEC_SIZE] __attribute__((aligned(16)));
  float hjq1[VEC_SIZE] __attribute__((aligned(16)));
  float dxq1[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq1[VEC_SIZE], *pjq1[VEC_SIZE];
  int icount2 = 0;
  float r2q2[VEC_SIZE] __attribute__((aligned(16)));
  float hiq2[VEC_SIZE] __attribute__((aligned(16)));
  float hjq2[VEC_SIZE] __attribute__((aligned(16)));
  float dxq2[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq2[VEC_SIZE], *pjq2[VEC_SIZE];
#endif

  TIMER_TIC;

  if (!cell_is_active(c, e)) return;

  if (!cell_are_part_drifted(c, e)) error("Cell is not drifted");

  struct part *restrict parts = c->parts;
  const int count = c->count;

  /* Set up indt. */
  int *indt = NULL;
  int countdt = 0, firstdt = 0;
  if (posix_memalign((void *)&indt, VEC_SIZE * sizeof(int),
                     count * sizeof(int)) != 0)
    error("Failed to allocate indt.");
  for (int k = 0; k < count; k++)
    if (part_is_active(&parts[k], e)) {
      indt[countdt] = k;
      countdt += 1;
    }

  /* Loop over the particles in the cell. */
  for (int pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    struct part *restrict pi = &parts[pid];

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

#ifndef WITH_OLD_VECTORIZATION

          IACT_NONSYM(r2, dx, hj, hi, pj, pi);

#else

          /* Add this interaction to the queue. */
          r2q1[icount1] = r2;
          dxq1[3 * icount1 + 0] = dx[0];
          dxq1[3 * icount1 + 1] = dx[1];
          dxq1[3 * icount1 + 2] = dx[2];
          hiq1[icount1] = hj;
          hjq1[icount1] = hi;
          piq1[icount1] = pj;
          pjq1[icount1] = pi;
          icount1 += 1;

          /* Flush? */
          if (icount1 == VEC_SIZE) {
            IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
            icount1 = 0;
          }

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

#ifndef WITH_OLD_VECTORIZATION

          /* Does pj need to be updated too? */
          if (part_is_active(pj, e))
            IACT(r2, dx, hi, hj, pi, pj);
          else
            IACT_NONSYM(r2, dx, hi, hj, pi, pj);

#else

          /* Does pj need to be updated too? */
          if (part_is_active(pj, e)) {

            /* Add this interaction to the symmetric queue. */
            r2q2[icount2] = r2;
            dxq2[3 * icount2 + 0] = dx[0];
            dxq2[3 * icount2 + 1] = dx[1];
            dxq2[3 * icount2 + 2] = dx[2];
            hiq2[icount2] = hi;
            hjq2[icount2] = hj;
            piq2[icount2] = pi;
            pjq2[icount2] = pj;
            icount2 += 1;

            /* Flush? */
            if (icount2 == VEC_SIZE) {
              IACT_VEC(r2q2, dxq2, hiq2, hjq2, piq2, pjq2);
              icount2 = 0;
            }

          } else {

            /* Add this interaction to the non-symmetric queue. */
            r2q1[icount1] = r2;
            dxq1[3 * icount1 + 0] = dx[0];
            dxq1[3 * icount1 + 1] = dx[1];
            dxq1[3 * icount1 + 2] = dx[2];
            hiq1[icount1] = hi;
            hjq1[icount1] = hj;
            piq1[icount1] = pi;
            pjq1[icount1] = pj;
            icount1 += 1;

            /* Flush? */
            if (icount1 == VEC_SIZE) {
              IACT_NONSYM_VEC(r2q1, dxq1, hiq1, hjq1, piq1, pjq1);
              icount1 = 0;
            }
          }

#endif
        }

      } /* loop over all other particles. */
    }

  } /* loop over all particles. */

#ifdef WITH_OLD_VECTORIZATION
  /* Pick up any leftovers. */
  if (icount1 > 0)
    for (int k = 0; k < icount1; k++)
      IACT_NONSYM(r2q1[k], &dxq1[3 * k], hiq1[k], hjq1[k], piq1[k], pjq1[k]);
  if (icount2 > 0)
    for (int k = 0; k < icount2; k++)
      IACT(r2q2[k], &dxq2[3 * k], hiq2[k], hjq2[k], piq2[k], pjq2[k]);
#endif

  free(indt);

  TIMER_TOC(TIMER_DOSELF);
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
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;
  if (ci->count == 0 || cj->count == 0) return;

  /* Get the type of pair if not specified explicitly. */
  double shift[3];
  sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_task(ci) && cell_can_recurse_in_pair_task(cj)) {

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
  else if (cell_is_active(ci, e) || cell_is_active(cj, e)) {

    /* Make sure both cells are drifted to the current timestep. */
    if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
      error("Interacting undrifted cells.");

    /* Do any of the cells need to be sorted first? */
    if (!(ci->sorted & (1 << sid)) ||
        ci->dx_max_sort_old > ci->dmin * space_maxreldx)
      error("Interacting unsorted cell.");
    if (!(cj->sorted & (1 << sid)) ||
        cj->dx_max_sort_old > cj->dmin * space_maxreldx)
      error("Interacting unsorted cell.");

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
  if (ci->count == 0 || !cell_is_active(ci, r->e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_task(ci)) {

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

#if (DOSELF1 == runner_doself1_density) && defined(WITH_VECTORIZATION) && \
    defined(GADGET2_SPH)
    runner_doself1_density_vec(r, ci);
#else
    DOSELF1(r, ci);
#endif
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
  if (!cell_is_active(ci, e) && !cell_is_active(cj, e)) return;
  if (ci->count == 0 || cj->count == 0) return;

  /* Get the type of pair if not specified explicitly. */
  double shift[3];
  sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_task(ci) && cell_can_recurse_in_pair_task(cj)) {

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
  else if (cell_is_active(ci, e) || cell_is_active(cj, e)) {

    /* Make sure both cells are drifted to the current timestep. */
    if (!cell_are_part_drifted(ci, e) || !cell_are_part_drifted(cj, e))
      error("Interacting undrifted cells.");

    /* Do any of the cells need to be sorted first? */
    if (!(ci->sorted & (1 << sid)) ||
        ci->dx_max_sort_old > ci->dmin * space_maxreldx)
      error("Interacting unsorted cells.");
    if (!(cj->sorted & (1 << sid)) ||
        cj->dx_max_sort_old > cj->dmin * space_maxreldx)
      error("Interacting unsorted cells.");

    /* Compute the interactions. */
    DOPAIR2(r, ci, cj);
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
  if (ci->count == 0 || !cell_is_active(ci, r->e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_task(ci)) {

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
  else
    DOSELF2(r, ci);

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}

void DOSUB_SUBSET(struct runner *r, struct cell *ci, struct part *parts,
                  int *ind, int count, struct cell *cj, int sid, int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  TIMER_TIC;

  /* Should we even bother? */
  if (!cell_is_active(ci, e) && (cj == NULL || !cell_is_active(cj, e))) return;
  if (ci->count == 0 || (cj != NULL && cj->count == 0)) return;

  /* Find out in which sub-cell of ci the parts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&parts[ind[0]] >= &ci->progeny[k]->parts[0] &&
            &parts[ind[0]] < &ci->progeny[k]->parts[ci->progeny[k]->count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_task(ci)) {

      /* Loop over all progeny. */
      DOSUB_SUBSET(r, sub, parts, ind, count, NULL, -1, 0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          DOSUB_SUBSET(r, sub, parts, ind, count, ci->progeny[j], -1, 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF_SUBSET(r, ci, parts, ind, count);

  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_task(ci) &&
        cell_can_recurse_in_pair_task(cj)) {

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
    else if (cell_is_active(ci, e) || cell_is_active(cj, e)) {

      /* Do any of the cells need to be drifted first? */
      if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");

      DOPAIR_SUBSET(r, ci, parts, ind, count, cj);
    }

  } /* otherwise, pair interaction. */

  if (gettimer) TIMER_TOC(timer_dosub_subset);
}
