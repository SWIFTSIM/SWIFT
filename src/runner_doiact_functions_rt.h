/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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

#include "runner_doiact_rt.h"

/**
 * @brief Function for self-type interaction between stars and hydro particles
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_RT(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0 || c->stars.count == 0) return;

  const struct engine *e = r->e;
  /* TODO: will be looked into in next MR */
  /* if (!cell_are_part_drifted(c, e) || !cell_are_spart_drifted(c, e)) */
  /*   error("Cell should be drifted!"); */

  struct spart *restrict sparts = c->stars.parts;
  struct part *restrict parts = c->hydro.parts;

  const int scount = c->stars.count;
  const int count = c->hydro.count;

  /* Loop over the sparts in cell */
  for (int sid = 0; sid < scount; sid++) {

    struct spart *restrict si = &sparts[sid];

    /* Skip inhibited particles. */
    if (spart_is_inhibited(si, e)) continue;

    const float hi = si->h;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the parts in cell */
    for (int pid = 0; pid < count; pid++) {
      struct part *restrict pj = &parts[pid];

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                            (float)(pj->x[1] - c->loc[1]),
                            (float)(pj->x[2] - c->loc[2])};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (r2 < hjg2) IACT_RT(r2, dx, hi, hj, si, pj);
    }
  }

  if (timer) TIMER_TOC(TIMER_DOSELF_RT);
}

/**
 * @brief Function for non-symmetric pair-type interaction between stars
 *        and hydro particles. Will interact star particles of cell i
 *        with hydro particles of cell j.
 *
 *
 * @param r runner task
 * @param ci the first cell, where we take star particles from
 * @param cj the second cell, where we take hydro particles from
 */
void DOPAIR1_NONSYM_RT(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *e = r->e;

  const int scount_i = ci->stars.count;
  const int count_j = cj->hydro.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount_i; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts_i[sid];

    /* Skip inhibited particles. */
    if (spart_is_inhibited(si, e)) continue;

    const float hi = si->h;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      const float hj = pj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Skip inhibited particles. */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                            (float)(pj->x[1] - cj->loc[1]),
                            (float)(pj->x[2] - cj->loc[2])};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

      if (r2 < hjg2) IACT_RT(r2, dx, hi, hj, si, pj);

    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Function for pair-type interaction between stars
 *        and hydro particles. Will interact hydro particles of cell i
 *        with star particles of cell j.
 *
 * @param r runner task
 * @param ci the first cell
 * @param cj the second cell
 * @param timer 1 if the time is to be recorded.
 */
void DOPAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj, int timer) {

  TIMER_TIC;
  /* const struct engine *restrict e = r->e; */

  const int do_stars_in_ci = (cj->nodeID == r->e->nodeID) &&
                             (ci->stars.count != 0) && (cj->hydro.count != 0);
  if (do_stars_in_ci) {
    /* TODO: will be looked into in next MR */
    /* if (!cell_are_spart_drifted(ci, e) || !cell_are_part_drifted(cj, e)) */
    /*   error("Cell should be drifted!"); */
    DOPAIR1_NONSYM_RT(r, ci, cj);
  }

  const int do_stars_in_cj = (ci->nodeID == r->e->nodeID) &&
                             (cj->stars.count != 0) && (ci->hydro.count != 0);
  if (do_stars_in_cj) {
    /* TODO: will be looked into in next MR */
    /* if (!cell_are_spart_drifted(cj, e) || !cell_are_part_drifted(ci, e)) */
    /*   error("Cell should be drifted!"); */
    DOPAIR1_NONSYM_RT(r, cj, ci);
  }

  if (timer) TIMER_TOC(TIMER_DOPAIR_RT);
}

/**
 * @brief Determine which version of DOSELF1_RT needs to be called
 *
 * @param r #runner
 * @param c #cell c
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_BRANCH_RT(struct runner *r, struct cell *c, int timer) {
  DOSELF1_RT(r, c, timer);
}

/**
 * @brief Determine which version of DOPAIR1_RT needs to be called
 *
 * @param r #runner
 * @param c #cell c
 * @param timer 1 if the time is to be recorded.
 */
void DOPAIR1_BRANCH_RT(struct runner *r, struct cell *ci, struct cell *cj,
                       int timer) {
  DOPAIR1_RT(r, ci, cj, timer);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_RT(struct runner *r, struct cell *c, int timer) {
  DOSELF1_RT(r, c, timer);
}

/**
 * @brief Compute grouped sub-cell interactions for pair tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_PAIR1_RT(struct runner *r, struct cell *ci, struct cell *cj,
                    int timer) {
  DOPAIR1_RT(r, ci, cj, timer);
}
