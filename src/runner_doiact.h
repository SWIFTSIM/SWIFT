
/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Includes. */
#include "cell.h"
#include "part.h"

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOPAIR1(f) PASTE(runner_dopair1, f)
#define DOPAIR1 _DOPAIR1(FUNCTION)

#define _DOPAIR2(f) PASTE(runner_dopair2, f)
#define DOPAIR2 _DOPAIR2(FUNCTION)

#define _DOPAIR_SUBSET(f) PASTE(runner_dopair_subset, f)
#define DOPAIR_SUBSET _DOPAIR_SUBSET(FUNCTION)

#define _DOPAIR_SUBSET_NAIVE(f) PASTE(runner_dopair_subset_naive, f)
#define DOPAIR_SUBSET_NAIVE _DOPAIR_SUBSET_NAIVE(FUNCTION)

#define _DOPAIR_NAIVE(f) PASTE(runner_dopair_naive, f)
#define DOPAIR_NAIVE _DOPAIR_NAIVE(FUNCTION)

#define _DOSELF_NAIVE(f) PASTE(runner_doself_naive, f)
#define DOSELF_NAIVE _DOSELF_NAIVE(FUNCTION)

#define _DOSELF1(f) PASTE(runner_doself1, f)
#define DOSELF1 _DOSELF1(FUNCTION)

#define _DOSELF2(f) PASTE(runner_doself2, f)
#define DOSELF2 _DOSELF2(FUNCTION)

#define _DOSELF_SUBSET(f) PASTE(runner_doself_subset, f)
#define DOSELF_SUBSET _DOSELF_SUBSET(FUNCTION)

#define _DOSUB1(f) PASTE(runner_dosub1, f)
#define DOSUB1 _DOSUB1(FUNCTION)

#define _DOSUB2(f) PASTE(runner_dosub2, f)
#define DOSUB2 _DOSUB2(FUNCTION)

#define _DOSUB_SUBSET(f) PASTE(runner_dosub_subset, f)
#define DOSUB_SUBSET _DOSUB_SUBSET(FUNCTION)

#define _IACT_NONSYM(f) PASTE(runner_iact_nonsym, f)
#define IACT_NONSYM _IACT_NONSYM(FUNCTION)

#define _IACT(f) PASTE(runner_iact, f)
#define IACT _IACT(FUNCTION)

#define _TIMER_DOSELF(f) PASTE(timer_doself, f)
#define TIMER_DOSELF _TIMER_DOSELF(FUNCTION)

#define _TIMER_DOPAIR(f) PASTE(timer_dopair, f)
#define TIMER_DOPAIR _TIMER_DOPAIR(FUNCTION)

#define _TIMER_DOSUB(f) PASTE(timer_dosub, f)
#define TIMER_DOSUB _TIMER_DOSUB(FUNCTION)

#define _TIMER_DOSELF_SUBSET(f) PASTE(timer_doself_subset, f)
#define TIMER_DOSELF_SUBSET _TIMER_DOSELF_SUBSET(FUNCTION)

#define _TIMER_DOPAIR_SUBSET(f) PASTE(timer_dopair_subset, f)
#define TIMER_DOPAIR_SUBSET _TIMER_DOPAIR_SUBSET(FUNCTION)

#define _IACT_NONSYM_VEC(f) PASTE(runner_iact_nonsym_vec, f)
#define IACT_NONSYM_VEC _IACT_NONSYM_VEC(FUNCTION)

#define _IACT_VEC(f) PASTE(runner_iact_vec, f)
#define IACT_VEC _IACT_VEC(FUNCTION)

/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */

void DOPAIR_NAIVE(struct runner *r, struct cell *restrict ci,
                  struct cell *restrict cj) {

  struct engine *e = r->e;
  int pid, pjd, k, count_i = ci->count, count_j = cj->count;
  double shift[3] = {0.0, 0.0, 0.0};
  struct part *restrict parts_i = ci->parts, *restrict parts_j = cj->parts;
  struct part *restrict pi, *restrict pj;
  double pix[3];
  float dx[3], hi, hig2, r2;
  const int ti_current = e->ti_current;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

  /* Get the relative distance between the pairs, wrapping. */
  for (k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* printf( "runner_dopair_naive: doing pair [ %g %g %g ]/[ %g %g %g ] with
  %i/%i parts and shift = [ %g %g %g ].\n" ,
      ci->loc[0] , ci->loc[1] , ci->loc[2] , cj->loc[0] , cj->loc[1] ,
  cj->loc[2] ,
      ci->count , cj->count , shift[0] , shift[1] , shift[2] ); fflush(stdout);
  tic = getticks(); */

  /* Loop over the parts in ci. */
  for (pid = 0; pid < count_i; pid++) {

    /* Get a hold of the ith part in ci. */
    pi = &parts_i[pid];
    for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < pj->h * pj->h * kernel_gamma2) {

#ifndef VECTORIZE

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      IACT(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

  TIMER_TOC(TIMER_DOPAIR);
}

void DOSELF_NAIVE(struct runner *r, struct cell *restrict c) {

  int pid, pjd, k, count = c->count;
  struct part *restrict parts = c->parts;
  struct part *restrict pi, *restrict pj;
  double pix[3] = {0.0, 0.0, 0.0};
  float dx[3], hi, hig2, r2;
  const int ti_current = r->e->ti_current;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Anything to do here? */
  if (c->ti_end_min > ti_current) return;

  /* printf( "runner_dopair_naive: doing pair [ %g %g %g ]/[ %g %g %g ] with
  %i/%i parts and shift = [ %g %g %g ].\n" ,
      ci->loc[0] , ci->loc[1] , ci->loc[2] , cj->loc[0] , cj->loc[1] ,
  cj->loc[2] ,
      ci->count , cj->count , shift[0] , shift[1] , shift[2] ); fflush(stdout);
  tic = getticks(); */

  /* Loop over the parts in ci. */
  for (pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    pi = &parts[pid];
    pix[0] = pi->x[0];
    pix[1] = pi->x[1];
    pix[2] = pi->x[2];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (pjd = pid + 1; pjd < count; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts[pjd];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2 || r2 < pj->h * pj->h * kernel_gamma2) {

#ifndef VECTORIZE

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      IACT(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

  TIMER_TOC(TIMER_DOSELF);
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
  int pid, pjd, sid, k, count_j = cj->count, flipped;
  double shift[3] = {0.0, 0.0, 0.0};
  struct part *restrict pi, *restrict pj, *restrict parts_j = cj->parts;
  double pix[3];
  float dx[3], hi, hig2, r2, di, dxj;
  struct entry *sort_j;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Get the relative distance between the pairs, wrapping. */
  for (k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Get the sorting index. */
  for (sid = 0, k = 0; k < 3; k++)
    sid = 3 * sid + ((cj->loc[k] - ci->loc[k] + shift[k] < 0)
                         ? 0
                         : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2 : 1);

  /* Switch the cells around? */
  flipped = runner_flip[sid];
  sid = sortlistID[sid];

  /* Have the cells been sorted? */
  if (!(cj->sorted & (1 << sid))) error("Trying to interact unsorted cells.");

  /* printf( "runner_dopair_naive: doing pair [ %g %g %g ]/[ %g %g %g ] with
  %i/%i parts and shift = [ %g %g %g ].\n" ,
      ci->loc[0] , ci->loc[1] , ci->loc[2] , cj->loc[0] , cj->loc[1] ,
  cj->loc[2] ,
      ci->count , cj->count , shift[0] , shift[1] , shift[2] ); fflush(stdout);
  tic = getticks(); */

  /* Pick-out the sorted lists. */
  sort_j = &cj->sort[sid * (cj->count + 1)];
  dxj = cj->dx_max;

  /* Parts are on the left? */
  if (!flipped) {

    /* Loop over the parts_i. */
    for (pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      pi = &parts_i[ind[pid]];
      for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
      hi = pi->h;
      hig2 = hi * hi * kernel_gamma2;
      di = hi * kernel_gamma + dxj + pix[0] * runner_shift[3 * sid + 0] +
           pix[1] * runner_shift[3 * sid + 1] +
           pix[2] * runner_shift[3 * sid + 2];

      /* Loop over the parts in cj. */
      for (pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        pj = &parts_j[sort_j[pjd].i];

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2) {

#ifndef VECTORIZE

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
    for (pid = 0; pid < count; pid++) {

      /* Get a hold of the ith part in ci. */
      pi = &parts_i[ind[pid]];
      for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
      hi = pi->h;
      hig2 = hi * hi * kernel_gamma2;
      di = -hi * kernel_gamma - dxj + pix[0] * runner_shift[3 * sid + 0] +
           pix[1] * runner_shift[3 * sid + 1] +
           pix[2] * runner_shift[3 * sid + 2];

      /* Loop over the parts in cj. */
      for (pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {

        /* Get a pointer to the jth particle. */
        pj = &parts_j[sort_j[pjd].i];

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2) {

#ifndef VECTORIZE

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
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
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 */

void DOPAIR_SUBSET_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct part *restrict parts_i, int *restrict ind,
                         int count, struct cell *restrict cj) {

  struct engine *e = r->e;
  int pid, pjd, k, count_j = cj->count;
  double shift[3] = {0.0, 0.0, 0.0};
  struct part *restrict pi, *restrict pj, *restrict parts_j = cj->parts;
  double pix[3];
  float dx[3], hi, hig2, r2;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Get the relative distance between the pairs, wrapping. */
  for (k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* printf( "runner_dopair_naive: doing pair [ %g %g %g ]/[ %g %g %g ] with
  %i/%i parts and shift = [ %g %g %g ].\n" ,
      ci->loc[0] , ci->loc[1] , ci->loc[2] , cj->loc[0] , cj->loc[1] ,
  cj->loc[2] ,
      ci->count , cj->count , shift[0] , shift[1] , shift[2] ); fflush(stdout);
  tic = getticks(); */

  /* Loop over the parts_i. */
  for (pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    pi = &parts_i[ind[pid]];
    for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {

#ifndef VECTORIZE

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
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

  int pid, pjd, k, count_i = ci->count;
  struct part *restrict parts_j = ci->parts;
  struct part *restrict pi, *restrict pj;
  double pix[3] = {0.0, 0.0, 0.0};
  float dx[3], hi, hig2, r2;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* printf( "runner_dopair_naive: doing pair [ %g %g %g ]/[ %g %g %g ] with
  %i/%i parts and shift = [ %g %g %g ].\n" ,
      ci->loc[0] , ci->loc[1] , ci->loc[2] , cj->loc[0] , cj->loc[1] ,
  cj->loc[2] ,
      ci->count , cj->count , shift[0] , shift[1] , shift[2] ); fflush(stdout);
  tic = getticks(); */

  /* Loop over the parts in ci. */
  for (pid = 0; pid < count; pid++) {

    /* Get a hold of the ith part in ci. */
    pi = &parts[ind[pid]];
    pix[0] = pi->x[0];
    pix[1] = pi->x[1];
    pix[2] = pi->x[2];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Loop over the parts in cj. */
    for (pjd = 0; pjd < count_i; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts_j[pjd];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 > 0.0f && r2 < hig2) {

#ifndef VECTORIZE

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      IACT_NONSYM(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

  TIMER_TOC(timer_dopair_subset);
}

/**
 * @brief Compute the interactions between a cell pair.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */

void DOPAIR1(struct runner *r, struct cell *ci, struct cell *cj) {

  struct engine *restrict e = r->e;
  int pid, pjd, k, sid;
  double rshift, shift[3] = {0.0, 0.0, 0.0};
  struct entry *restrict sort_i, *restrict sort_j;
  struct part *restrict pi, *restrict pj, *restrict parts_i, *restrict parts_j;
  double pix[3], pjx[3], di, dj;
  float dx[3], hi, hig2, hj, hjg2, r2, dx_max;
  double hi_max, hj_max;
  double di_max, dj_min;
  int count_i, count_j;
  const int ti_current = e->ti_current;
#ifdef VECTORIZE
  int icount = 0;
  float r2q[VEC_SIZE] __attribute__((aligned(16)));
  float hiq[VEC_SIZE] __attribute__((aligned(16)));
  float hjq[VEC_SIZE] __attribute__((aligned(16)));
  float dxq[3 * VEC_SIZE] __attribute__((aligned(16)));
  struct part *piq[VEC_SIZE], *pjq[VEC_SIZE];
#endif
  TIMER_TIC

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

  /* Get the sort ID. */
  sid = space_getsid(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->sorted & (1 << sid)) || !(cj->sorted & (1 << sid)))
    error("Trying to interact unsorted cells.");

  /* Get the cutoff shift. */
  for (rshift = 0.0, k = 0; k < 3; k++)
    rshift += shift[k] * runner_shift[3 * sid + k];

  /* Pick-out the sorted lists. */
  sort_i = &ci->sort[sid * (ci->count + 1)];
  sort_j = &cj->sort[sid * (cj->count + 1)];

  /* Get some other useful values. */
  hi_max = ci->h_max * kernel_gamma - rshift;
  hj_max = cj->h_max * kernel_gamma;
  count_i = ci->count;
  count_j = cj->count;
  parts_i = ci->parts;
  parts_j = cj->parts;
  di_max = sort_i[count_i - 1].d - rshift;
  dj_min = sort_j[0].d;
  dx_max = (ci->dx_max + cj->dx_max);

  /* Loop over the parts in ci. */
  for (pid = count_i - 1; pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min;
       pid--) {

    /* Get a hold of the ith part in ci. */
    pi = &parts_i[sort_i[pid].i];
    if (pi->ti_end > ti_current) continue;
    hi = pi->h;
    di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
    if (di < dj_min) continue;

    hig2 = hi * hi * kernel_gamma2;
    for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];

    /* Loop over the parts in cj. */
    for (pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

      /* Get a pointer to the jth particle. */
      pj = &parts_j[sort_j[pjd].i];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pix[k] - pj->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hig2) {

#ifndef VECTORIZE

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

  /* printf( "runner_dopair: first half took %.3f %s...\n" ,
  clocks_from_ticks(getticks() - tic), clocks_getunit());
  tic = getticks(); */

  /* Loop over the parts in cj. */
  for (pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
       pjd++) {

    /* Get a hold of the jth part in cj. */
    pj = &parts_j[sort_j[pjd].i];
    if (pj->ti_end > ti_current) continue;
    hj = pj->h;
    dj = sort_j[pjd].d - hj * kernel_gamma - dx_max - rshift;
    if (dj > di_max) continue;

    for (k = 0; k < 3; k++) pjx[k] = pj->x[k] + shift[k];
    hjg2 = hj * hj * kernel_gamma2;

    /* Loop over the parts in ci. */
    for (pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

      /* Get a pointer to the jth particle. */
      pi = &parts_i[sort_i[pid].i];

      /* Compute the pairwise distance. */
      r2 = 0.0f;
      for (k = 0; k < 3; k++) {
        dx[k] = pjx[k] - pi->x[k];
        r2 += dx[k] * dx[k];
      }

      /* Hit or miss? */
      if (r2 < hjg2) {

#ifndef VECTORIZE

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

    } /* loop over the parts in cj. */

  } /* loop over the parts in ci. */

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount > 0)
    for (k = 0; k < icount; k++)
      IACT_NONSYM(r2q[k], &dxq[3 * k], hiq[k], hjq[k], piq[k], pjq[k]);
#endif

  TIMER_TOC(TIMER_DOPAIR);
}

void DOPAIR2(struct runner *r, struct cell *ci, struct cell *cj) {

  struct engine *restrict e = r->e;
  int pid, pjd, k, sid;
  double rshift, shift[3] = {0.0, 0.0, 0.0};
  struct entry *sort_i, *sort_j;
  struct entry *sortdt_i = NULL, *sortdt_j = NULL;
  int countdt_i = 0, countdt_j = 0;
  struct part *restrict pi, *restrict pj, *restrict parts_i, *restrict parts_j;
  double pix[3], pjx[3], di, dj;
  float dx[3], hi, hig2, hj, hjg2, r2, dx_max;
  double hi_max, hj_max;
  double di_max, dj_min;
  int count_i, count_j;
  const int ti_current = e->ti_current;
#ifdef VECTORIZE
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
  TIMER_TIC

  /* Anything to do here? */
  if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

  /* Get the shift ID. */
  sid = space_getsid(e->s, &ci, &cj, shift);

  /* Have the cells been sorted? */
  if (!(ci->sorted & (1 << sid)) || !(cj->sorted & (1 << sid)))
    error("Trying to interact unsorted cells.");

  /* Get the cutoff shift. */
  for (rshift = 0.0, k = 0; k < 3; k++)
    rshift += shift[k] * runner_shift[3 * sid + k];

  /* Pick-out the sorted lists. */
  sort_i = &ci->sort[sid * (ci->count + 1)];
  sort_j = &cj->sort[sid * (cj->count + 1)];

  /* Get some other useful values. */
  hi_max = ci->h_max * kernel_gamma - rshift;
  hj_max = cj->h_max * kernel_gamma;
  count_i = ci->count;
  count_j = cj->count;
  parts_i = ci->parts;
  parts_j = cj->parts;
  di_max = sort_i[count_i - 1].d - rshift;
  dj_min = sort_j[0].d;
  dx_max = (ci->dx_max + cj->dx_max);

  /* Collect the number of parts left and right below dt. */
  if (ci->ti_end_max <= ti_current) {
    sortdt_i = sort_i;
    countdt_i = count_i;
  } else if (ci->ti_end_min <= ti_current) {
    if ((sortdt_i = (struct entry *)alloca(sizeof(struct entry) * count_i)) ==
        NULL)
      error("Failed to allocate dt sortlists.");
    for (k = 0; k < count_i; k++)
      if (parts_i[sort_i[k].i].ti_end <= ti_current) {
        sortdt_i[countdt_i] = sort_i[k];
        countdt_i += 1;
      }
  }
  if (cj->ti_end_max <= ti_current) {
    sortdt_j = sort_j;
    countdt_j = count_j;
  } else if (cj->ti_end_min <= ti_current) {
    if ((sortdt_j = (struct entry *)alloca(sizeof(struct entry) * count_j)) ==
        NULL)
      error("Failed to allocate dt sortlists.");
    for (k = 0; k < count_j; k++)
      if (parts_j[sort_j[k].i].ti_end <= ti_current) {
        sortdt_j[countdt_j] = sort_j[k];
        countdt_j += 1;
      }
  }

  /* Loop over the parts in ci. */
  for (pid = count_i - 1; pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min;
       pid--) {

    /* Get a hold of the ith part in ci. */
    pi = &parts_i[sort_i[pid].i];
    hi = pi->h;
    di = sort_i[pid].d + hi * kernel_gamma + dx_max - rshift;
    if (di < dj_min) continue;

    hig2 = hi * hi * kernel_gamma2;
    for (k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];

    /* Look at valid dt parts only? */
    if (pi->ti_end > ti_current) {

      /* Loop over the parts in cj within dt. */
      for (pjd = 0; pjd < countdt_j && sortdt_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        pj = &parts_j[sortdt_j[pjd].i];
        hj = pj->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pj->x[k] - pix[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2) {

#ifndef VECTORIZE

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
      for (pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

        /* Get a pointer to the jth particle. */
        pj = &parts_j[sort_j[pjd].i];
        hj = pj->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2) {

#ifndef VECTORIZE

          /* Does pj need to be updated too? */
          if (pj->ti_end <= ti_current)
            IACT(r2, dx, hi, hj, pi, pj);
          else
            IACT_NONSYM(r2, dx, hi, hj, pi, pj);

#else

          /* Does pj need to be updated too? */
          if (pj->ti_end <= ti_current) {

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

  /* printf( "runner_dopair: first half took %.3f %s...\n" ,
  clocks_from_ticks(getticks() - tic), clocks_getunit());
  tic = getticks(); */

  /* Loop over the parts in cj. */
  for (pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
       pjd++) {

    /* Get a hold of the jth part in cj. */
    pj = &parts_j[sort_j[pjd].i];
    hj = pj->h;
    dj = sort_j[pjd].d - hj * kernel_gamma - dx_max - rshift;
    if (dj > di_max) continue;

    for (k = 0; k < 3; k++) pjx[k] = pj->x[k] + shift[k];
    hjg2 = hj * hj * kernel_gamma2;

    /* Is this particle outside the dt? */
    if (pj->ti_end > ti_current) {

      /* Loop over the parts in ci. */
      for (pid = countdt_i - 1; pid >= 0 && sortdt_i[pid].d > dj; pid--) {

        /* Get a pointer to the jth particle. */
        pi = &parts_i[sortdt_i[pid].i];
        hi = pi->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pi->x[k] - pjx[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hjg2 && r2 > hi * hi * kernel_gamma2) {

#ifndef VECTORIZE

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
      for (pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

        /* Get a pointer to the jth particle. */
        pi = &parts_i[sort_i[pid].i];
        hi = pi->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pjx[k] - pi->x[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hjg2 && r2 > hi * hi * kernel_gamma2) {

#ifndef VECTORIZE

          /* Does pi need to be updated too? */
          if (pi->ti_end <= ti_current)
            IACT(r2, dx, hj, hi, pj, pi);
          else
            IACT_NONSYM(r2, dx, hj, hi, pj, pi);

#else

          /* Does pi need to be updated too? */
          if (pi->ti_end <= ti_current) {

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount1 > 0)
    for (k = 0; k < icount1; k++)
      IACT_NONSYM(r2q1[k], &dxq1[3 * k], hiq1[k], hjq1[k], piq1[k], pjq1[k]);
  if (icount2 > 0)
    for (k = 0; k < icount2; k++)
      IACT(r2q2[k], &dxq2[3 * k], hiq2[k], hjq2[k], piq2[k], pjq2[k]);
#endif

  TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the cell self-interaction.
 *
 * @param r The #runner.
 * @param c The #cell.
 */

void DOSELF1(struct runner *r, struct cell *restrict c) {

  int k, pid, pjd, count = c->count;
  double pix[3];
  float dx[3], hi, hj, hig2, r2;
  struct part *restrict parts = c->parts, *restrict pi, *restrict pj;
  const int ti_current = r->e->ti_current;
  int firstdt = 0, countdt = 0, *indt = NULL, doj;
#ifdef VECTORIZE
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
  TIMER_TIC

  /* Set up indt if needed. */
  if (c->ti_end_min > ti_current)
    return;
  else if (c->ti_end_max > ti_current) {
    if ((indt = (int *)alloca(sizeof(int) * count)) == NULL)
      error("Failed to allocate indt.");
    for (k = 0; k < count; k++)
      if (parts[k].ti_end <= ti_current) {
        indt[countdt] = k;
        countdt += 1;
      }
  }

  /* Loop over the particles in the cell. */
  for (pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    pi = &parts[pid];

    /* Get the particle position and radius. */
    for (k = 0; k < 3; k++) pix[k] = pi->x[k];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle inactive? */
    if (pi->ti_end > ti_current) {

      /* Loop over the other particles .*/
      for (pjd = firstdt; pjd < countdt; pjd++) {

        /* Get a pointer to the jth particle. */
        pj = &parts[indt[pjd]];
        hj = pj->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pj->x[k] - pix[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hj * hj * kernel_gamma2) {

#ifndef VECTORIZE

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
      for (pjd = pid + 1; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        pj = &parts[pjd];
        hj = pj->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }
        doj = (pj->ti_end <= ti_current) && (r2 < hj * hj * kernel_gamma2);

        /* Hit or miss? */
        if (r2 < hig2 || doj) {

#ifndef VECTORIZE

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount1 > 0)
    for (k = 0; k < icount1; k++)
      IACT_NONSYM(r2q1[k], &dxq1[3 * k], hiq1[k], hjq1[k], piq1[k], pjq1[k]);
  if (icount2 > 0)
    for (k = 0; k < icount2; k++)
      IACT(r2q2[k], &dxq2[3 * k], hiq2[k], hjq2[k], piq2[k], pjq2[k]);
#endif

  TIMER_TOC(TIMER_DOSELF);
}

void DOSELF2(struct runner *r, struct cell *restrict c) {

  int k, pid, pjd, count = c->count;
  double pix[3];
  float dx[3], hi, hj, hig2, r2;
  struct part *restrict parts = c->parts, *restrict pi, *restrict pj;
  const int ti_current = r->e->ti_current;
  int firstdt = 0, countdt = 0, *indt = NULL;
#ifdef VECTORIZE
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
  TIMER_TIC

  /* Set up indt if needed. */
  if (c->ti_end_min > ti_current)
    return;
  else if (c->ti_end_max > ti_current) {
    if ((indt = (int *)alloca(sizeof(int) * count)) == NULL)
      error("Failed to allocate indt.");
    for (k = 0; k < count; k++)
      if (parts[k].ti_end <= ti_current) {
        indt[countdt] = k;
        countdt += 1;
      }
  }

  /* Loop over the particles in the cell. */
  for (pid = 0; pid < count; pid++) {

    /* Get a pointer to the ith particle. */
    pi = &parts[pid];

    /* Get the particle position and radius. */
    for (k = 0; k < 3; k++) pix[k] = pi->x[k];
    hi = pi->h;
    hig2 = hi * hi * kernel_gamma2;

    /* Is the ith particle not active? */
    if (pi->ti_end > ti_current) {

      /* Loop over the other particles .*/
      for (pjd = firstdt; pjd < countdt; pjd++) {

        /* Get a pointer to the jth particle. */
        pj = &parts[indt[pjd]];
        hj = pj->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pj->x[k] - pix[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2 || r2 < hj * hj * kernel_gamma2) {

#ifndef VECTORIZE

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
      for (pjd = pid + 1; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        pj = &parts[pjd];
        hj = pj->h;

        /* Compute the pairwise distance. */
        r2 = 0.0f;
        for (k = 0; k < 3; k++) {
          dx[k] = pix[k] - pj->x[k];
          r2 += dx[k] * dx[k];
        }

        /* Hit or miss? */
        if (r2 < hig2 || r2 < hj * hj * kernel_gamma2) {

#ifndef VECTORIZE

          /* Does pj need to be updated too? */
          if (pj->ti_end <= ti_current)
            IACT(r2, dx, hi, hj, pi, pj);
          else
            IACT_NONSYM(r2, dx, hi, hj, pi, pj);

#else

          /* Does pj need to be updated too? */
          if (pj->ti_end <= ti_current) {

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

#ifdef VECTORIZE
  /* Pick up any leftovers. */
  if (icount1 > 0)
    for (k = 0; k < icount1; k++)
      IACT_NONSYM(r2q1[k], &dxq1[3 * k], hiq1[k], hjq1[k], piq1[k], pjq1[k]);
  if (icount2 > 0)
    for (k = 0; k < icount2; k++)
      IACT(r2q2[k], &dxq2[3 * k], hiq2[k], hjq2[k], piq2[k], pjq2[k]);
#endif

  TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute grouped sub-cell interactions
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

void DOSUB1(struct runner *r, struct cell *ci, struct cell *cj, int sid,
            int gettimer) {

  int j = 0, k;
  double shift[3];
  float h;
  struct space *s = r->e->s;
  const int ti_current = r->e->ti_current;

  TIMER_TIC

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Should we even bother? */
    if (ci->ti_end_min > ti_current) return;

    /* Recurse? */
    if (ci->split) {

      /* Loop over all progeny. */
      for (k = 0; k < 8; k++)
        if (ci->progeny[k] != NULL) {
          DOSUB1(r, ci->progeny[k], NULL, -1, 0);
          for (j = k + 1; j < 8; j++)
            if (ci->progeny[j] != NULL)
              DOSUB1(r, ci->progeny[k], ci->progeny[j], -1, 0);
        }

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF1(r, ci);

  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Should we even bother? */
    if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

    /* Get the cell dimensions. */
    h = fmin(ci->h[0], fmin(ci->h[1], ci->h[2]));

    /* Get the type of pair if not specified explicitly. */
    // if ( sid < 0 )
    sid = space_getsid(s, &ci, &cj, shift);

    /* Recurse? */
    if (ci->split && cj->split &&
        fmaxf(ci->h_max, cj->h_max) * kernel_gamma + ci->dx_max + cj->dx_max <
            h / 2) {

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[0], -1, 0);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[1], -1, 0);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[1], -1, 0);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[2], -1, 0);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[0], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[1], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[2], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[3], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[3], -1, 0);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[3], -1, 0);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[2], -1, 0);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[3], -1, 0);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB1(r, ci->progeny[4], cj->progeny[3], -1, 0);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[4], -1, 0);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[2], cj->progeny[0], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[2], cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[2], cj->progeny[4], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            DOSUB1(r, ci->progeny[2], cj->progeny[5], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[1], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[4], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[5], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[5], -1, 0);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[2], cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            DOSUB1(r, ci->progeny[2], cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            DOSUB1(r, ci->progeny[6], cj->progeny[5], -1, 0);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[1], cj->progeny[0], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[1], cj->progeny[2], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[1], cj->progeny[4], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
            DOSUB1(r, ci->progeny[1], cj->progeny[6], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[2], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
            DOSUB1(r, ci->progeny[3], cj->progeny[6], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[4], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
            DOSUB1(r, ci->progeny[5], cj->progeny[6], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
            DOSUB1(r, ci->progeny[7], cj->progeny[6], -1, 0);
          break;
      }

    }

    /* Otherwise, compute the pair directly. */
    else if (ci->ti_end_min <= ti_current || cj->ti_end_min <= ti_current) {

      /* Do any of the cells need to be sorted first? */
      if (!(ci->sorted & (1 << sid))) runner_dosort(r, ci, (1 << sid), 1);
      if (!(cj->sorted & (1 << sid))) runner_dosort(r, cj, (1 << sid), 1);

      /* Compute the interactions. */
      DOPAIR1(r, ci, cj);
    }

  } /* otherwise, pair interaction. */

  if (gettimer) TIMER_TOC(TIMER_DOSUB);
}

void DOSUB2(struct runner *r, struct cell *ci, struct cell *cj, int sid,
            int gettimer) {

  int j, k;
  double shift[3];
  float h;
  struct space *s = r->e->s;
  const int ti_current = r->e->ti_current;

  TIMER_TIC

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Should we even bother? */
    if (ci->ti_end_min > ti_current) return;

    /* Recurse? */
    if (ci->split) {

      /* Loop over all progeny. */
      for (k = 0; k < 8; k++)
        if (ci->progeny[k] != NULL) {
          DOSUB2(r, ci->progeny[k], NULL, -1, 0);
          for (j = k + 1; j < 8; j++)
            if (ci->progeny[j] != NULL)
              DOSUB2(r, ci->progeny[k], ci->progeny[j], -1, 0);
        }

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF2(r, ci);

  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Should we even bother? */
    if (ci->ti_end_min > ti_current && cj->ti_end_min > ti_current) return;

    /* Get the cell dimensions. */
    h = fmin(ci->h[0], fmin(ci->h[1], ci->h[2]));

    /* Get the type of pair if not specified explicitly. */
    // if ( sid < 0 )
    sid = space_getsid(s, &ci, &cj, shift);

    /* Recurse? */
    if (ci->split && cj->split &&
        fmaxf(ci->h_max, cj->h_max) * kernel_gamma + ci->dx_max + cj->dx_max <
            h / 2) {

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[0], -1, 0);
          break;

        case 1: /* (  1 ,  1 ,  0 ) */
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[1], -1, 0);
          break;

        case 2: /* (  1 ,  1 , -1 ) */
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[1], -1, 0);
          break;

        case 3: /* (  1 ,  0 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[2], -1, 0);
          break;

        case 4: /* (  1 ,  0 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[0], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[1], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[2], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[3], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[3], -1, 0);
          break;

        case 5: /* (  1 ,  0 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[1], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[3], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[3], -1, 0);
          break;

        case 6: /* (  1 , -1 ,  1 ) */
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[2], -1, 0);
          break;

        case 7: /* (  1 , -1 ,  0 ) */
          if (ci->progeny[4] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[2], -1, 0);
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[3], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[3], -1, 0);
          break;

        case 8: /* (  1 , -1 , -1 ) */
          if (ci->progeny[4] != NULL && cj->progeny[3] != NULL)
            DOSUB2(r, ci->progeny[4], cj->progeny[3], -1, 0);
          break;

        case 9: /* (  0 ,  1 ,  1 ) */
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[4], -1, 0);
          break;

        case 10: /* (  0 ,  1 ,  0 ) */
          if (ci->progeny[2] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[2], cj->progeny[0], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[2], cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[2], cj->progeny[4], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            DOSUB2(r, ci->progeny[2], cj->progeny[5], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[1], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[5] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[0], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[4], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[5], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[1], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[5] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[5], -1, 0);
          break;

        case 11: /* (  0 ,  1 , -1 ) */
          if (ci->progeny[2] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[2], cj->progeny[1], -1, 0);
          if (ci->progeny[2] != NULL && cj->progeny[5] != NULL)
            DOSUB2(r, ci->progeny[2], cj->progeny[5], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[1] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[1], -1, 0);
          if (ci->progeny[6] != NULL && cj->progeny[5] != NULL)
            DOSUB2(r, ci->progeny[6], cj->progeny[5], -1, 0);
          break;

        case 12: /* (  0 ,  0 ,  1 ) */
          if (ci->progeny[1] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[1], cj->progeny[0], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[1], cj->progeny[2], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[1], cj->progeny[4], -1, 0);
          if (ci->progeny[1] != NULL && cj->progeny[6] != NULL)
            DOSUB2(r, ci->progeny[1], cj->progeny[6], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[0], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[2], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[4], -1, 0);
          if (ci->progeny[3] != NULL && cj->progeny[6] != NULL)
            DOSUB2(r, ci->progeny[3], cj->progeny[6], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[0], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[2], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[4], -1, 0);
          if (ci->progeny[5] != NULL && cj->progeny[6] != NULL)
            DOSUB2(r, ci->progeny[5], cj->progeny[6], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[0], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[2] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[2], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[4] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[4], -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[6] != NULL)
            DOSUB2(r, ci->progeny[7], cj->progeny[6], -1, 0);
          break;
      }

    }

    /* Otherwise, compute the pair directly. */
    else if (ci->ti_end_min <= ti_current || cj->ti_end_min <= ti_current) {

      /* Do any of the cells need to be sorted first? */
      if (!(ci->sorted & (1 << sid))) runner_dosort(r, ci, (1 << sid), 1);
      if (!(cj->sorted & (1 << sid))) runner_dosort(r, cj, (1 << sid), 1);

      /* Compute the interactions. */
      DOPAIR2(r, ci, cj);
    }

  } /* otherwise, pair interaction. */

  if (gettimer) TIMER_TOC(TIMER_DOSUB);
}

void DOSUB_SUBSET(struct runner *r, struct cell *ci, struct part *parts,
                  int *ind, int count, struct cell *cj, int sid, int gettimer) {

  int j, k;
  double shift[3];
  float h;
  struct space *s = r->e->s;
  struct cell *sub = NULL;
  const int ti_current = r->e->ti_current;

  TIMER_TIC

  /* Find out in which sub-cell of ci the parts are. */
  for (k = 0; k < 8; k++)
    if (ci->progeny[k] != NULL) {
      // if ( parts[ ind[ 0 ] ].x[0] >= ci->progeny[k]->loc[0] &&
      //      parts[ ind[ 0 ] ].x[0] <= ci->progeny[k]->loc[0] +
      // ci->progeny[k]->h[0] &&
      //      parts[ ind[ 0 ] ].x[1] >= ci->progeny[k]->loc[1] &&
      //      parts[ ind[ 0 ] ].x[1] <= ci->progeny[k]->loc[1] +
      // ci->progeny[k]->h[1] &&
      //      parts[ ind[ 0 ] ].x[2] >= ci->progeny[k]->loc[2] &&
      //      parts[ ind[ 0 ] ].x[2] <= ci->progeny[k]->loc[2] +
      // ci->progeny[k]->h[2] ) {
      if (&parts[ind[0]] >= &ci->progeny[k]->parts[0] &&
          &parts[ind[0]] < &ci->progeny[k]->parts[ci->progeny[k]->count]) {
        sub = ci->progeny[k];
        break;
      }
    }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (ci->split) {

      /* Loop over all progeny. */
      DOSUB_SUBSET(r, sub, parts, ind, count, NULL, -1, 0);
      for (j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          DOSUB_SUBSET(r, sub, parts, ind, count, ci->progeny[j], -1, 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF_SUBSET(r, ci, parts, ind, count);

  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Get the cell dimensions. */
    h = fmin(ci->h[0], fmin(ci->h[1], ci->h[2]));

    /* Recurse? */
    if (ci->split && cj->split &&
        fmaxf(ci->h_max, cj->h_max) * kernel_gamma + ci->dx_max + cj->dx_max <
            h / 2) {

      /* Get the type of pair if not specified explicitly. */
      sid = space_getsid(s, &ci, &cj, shift);

      /* Different types of flags. */
      switch (sid) {

        /* Regular sub-cell interactions of a single cell. */
        case 0: /* (  1 ,  1 ,  1 ) */
          if (ci->progeny[7] == sub && cj->progeny[0] != NULL)
            DOSUB_SUBSET(r, ci->progeny[7], parts, ind, count, cj->progeny[0],
                         -1, 0);
          if (ci->progeny[7] != NULL && cj->progeny[0] == sub)
            DOSUB_SUBSET(r, ci->progeny[0], parts, ind, count, cj->progeny[7],
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
    else if (ci->ti_end_min <= ti_current || cj->ti_end_min <= ti_current) {

      /* Get the relative distance between the pairs, wrapping. */
      for (k = 0; k < 3; k++) {
        if (cj->loc[k] - ci->loc[k] < -s->dim[k] / 2)
          shift[k] = s->dim[k];
        else if (cj->loc[k] - ci->loc[k] > s->dim[k] / 2)
          shift[k] = -s->dim[k];
      }

      /* Get the sorting index. */
      for (sid = 0, k = 0; k < 3; k++)
        sid =
            3 * sid + ((cj->loc[k] - ci->loc[k] + shift[k] < 0)
                           ? 0
                           : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2 : 1);
      sid = sortlistID[sid];

      /* Do any of the cells need to be sorted first? */
      if (!(cj->sorted & (1 << sid))) runner_dosort(r, cj, (1 << sid), 1);

      /* Compute the interactions. */
      DOPAIR_SUBSET(r, ci, parts, ind, count, cj);
    }

  } /* otherwise, pair interaction. */

  if (gettimer) TIMER_TOC(TIMER_DOSUB);
}
