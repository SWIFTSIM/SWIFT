/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

#include "runner_doiact_black_holes.h"

/**
 * @brief Calculate the number density of #part around the #bpart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_BH(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int bi_is_local = 1; /* SELF tasks are always local */

  /* Anything to do here? */
  if (c->black_holes.count == 0) return;
  if (!cell_is_active_black_holes(c, e)) return;

  const int bcount = c->black_holes.count;
  const int count = c->hydro.count;
  struct bpart *restrict bparts = c->black_holes.parts;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  /* Do we actually have any gas neighbours? */
  if (c->hydro.count != 0) {

    /* Loop over the bparts in ci. */
    for (int bid = 0; bid < bcount; bid++) {

      /* Get a hold of the ith bpart in ci. */
      struct bpart *restrict bi = &bparts[bid];

      /* Skip inactive particles */
      if (!bpart_is_active(bi, e)) continue;

      const float hi = bi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float bix[3] = {(float)(bi->x[0] - c->loc[0]),
                            (float)(bi->x[1] - c->loc[1]),
                            (float)(bi->x[2] - c->loc[2])};

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[pjd];
        struct xpart *restrict xpj = &xparts[pjd];
        const float hj = pj->h;

        /* Early abort? */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                              (float)(pj->x[1] - c->loc[1]),
                              (float)(pj->x[2] - c->loc[2])};
        const float dx[3] = {bix[0] - pjx[0], bix[1] - pjx[1], bix[2] - pjx[2]};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (bi->ti_drift != e->ti_current)
          error("Particle bi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        if (r2 < hig2) {
          IACT_BH_GAS(r2, dx, hi, hj, bi, pj, xpj, with_cosmology, cosmo,
                      e->gravity_properties, e->black_holes_properties,
                      e->entropy_floor, ti_current, e->time);

          if (bi_is_local) {
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
            runner_iact_nonsym_bh_gas_repos(
                r2, dx, hi, pj->h, bi, pj, xpj, with_cosmology, cosmo,
                e->gravity_properties, e->black_holes_properties,
                e->entropy_floor, ti_current, e->time);
#endif
          }
        }
      } /* loop over the parts in ci. */
    }   /* loop over the bparts in ci. */
  }     /* Do we have gas particles in the cell? */

  /* When doing BH swallowing, we need a quick loop also over the BH
   * neighbours */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)

  /* Loop over the bparts in ci. */
  for (int bid = 0; bid < bcount; bid++) {

    /* Get a hold of the ith bpart in ci. */
    struct bpart *restrict bi = &bparts[bid];

    /* Skip inactive particles */
    if (!bpart_is_active(bi, e)) continue;

    const float hi = bi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float bix[3] = {(float)(bi->x[0] - c->loc[0]),
                          (float)(bi->x[1] - c->loc[1]),
                          (float)(bi->x[2] - c->loc[2])};

    /* Loop over the parts in cj. */
    for (int bjd = 0; bjd < bcount; bjd++) {

      /* Skip self interaction */
      if (bid == bjd) continue;

      /* Get a pointer to the jth particle. */
      struct bpart *restrict bj = &bparts[bjd];
      const float hj = bj->h;

      /* Early abort? */
      if (bpart_is_inhibited(bj, e)) continue;

      /* Compute the pairwise distance. */
      const float bjx[3] = {(float)(bj->x[0] - c->loc[0]),
                            (float)(bj->x[1] - c->loc[1]),
                            (float)(bj->x[2] - c->loc[2])};
      const float dx[3] = {bix[0] - bjx[0], bix[1] - bjx[1], bix[2] - bjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (bi->ti_drift != e->ti_current)
        error("Particle bi not drifted to current time");
      if (bj->ti_drift != e->ti_current)
        error("Particle bj not drifted to current time");
#endif

      if (r2 < hig2) {
        IACT_BH_BH(r2, dx, hi, hj, bi, bj, cosmo, e->gravity_properties,
                   e->black_holes_properties, ti_current);

        if (bi_is_local) {
          runner_iact_nonsym_bh_bh_repos(r2, dx, hi, hj, bi, bj, cosmo,
                                         e->gravity_properties,
                                         e->black_holes_properties, ti_current);
        }
      }
    } /* loop over the bparts in ci. */
  }   /* loop over the bparts in ci. */

#endif /* (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW) */

  TIMER_TOC(TIMER_DOSELF_BH);
}

/**
 * @brief Calculate the number density of cj #part around the ci #bpart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_BH_NAIVE(struct runner *r, struct cell *restrict ci,
                              struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
  if (cj->nodeID != engine_rank) error("Should be run on a different node");
#endif
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int bi_is_local = ci->nodeID == e->nodeID;

  /* Anything to do here? */
  if (ci->black_holes.count == 0) return;
  if (!cell_is_active_black_holes(ci, e)) return;

  const int bcount_i = ci->black_holes.count;
  const int count_j = cj->hydro.count;
  struct bpart *restrict bparts_i = ci->black_holes.parts;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Do we actually have any gas neighbours? */
  if (cj->hydro.count != 0) {

    /* Loop over the bparts in ci. */
    for (int bid = 0; bid < bcount_i; bid++) {

      /* Get a hold of the ith bpart in ci. */
      struct bpart *restrict bi = &bparts_i[bid];

      /* Skip inactive particles */
      if (!bpart_is_active(bi, e)) continue;

      const float hi = bi->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float bix[3] = {(float)(bi->x[0] - (cj->loc[0] + shift[0])),
                            (float)(bi->x[1] - (cj->loc[1] + shift[1])),
                            (float)(bi->x[2] - (cj->loc[2] + shift[2]))};

      /* Loop over the parts in cj. */
      for (int pjd = 0; pjd < count_j; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[pjd];
        struct xpart *restrict xpj = &xparts_j[pjd];
        const float hj = pj->h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                              (float)(pj->x[1] - cj->loc[1]),
                              (float)(pj->x[2] - cj->loc[2])};
        const float dx[3] = {bix[0] - pjx[0], bix[1] - pjx[1], bix[2] - pjx[2]};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (bi->ti_drift != e->ti_current)
          error("Particle bi not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        if (r2 < hig2) {
          IACT_BH_GAS(r2, dx, hi, hj, bi, pj, xpj, with_cosmology, cosmo,
                      e->gravity_properties, e->black_holes_properties,
                      e->entropy_floor, ti_current, e->time);

          if (bi_is_local) {
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
            runner_iact_nonsym_bh_gas_repos(
                r2, dx, hi, hj, bi, pj, xpj, with_cosmology, cosmo,
                e->gravity_properties, e->black_holes_properties,
                e->entropy_floor, ti_current, e->time);
#endif
          }
        }
      } /* loop over the parts in cj. */
    }   /* loop over the bparts in ci. */
  }     /* Do we have gas particles in the cell? */

  /* When doing BH swallowing, we need a quick loop also over the BH
   * neighbours */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)

  const int bcount_j = cj->black_holes.count;
  struct bpart *restrict bparts_j = cj->black_holes.parts;

  /* Loop over the bparts in ci. */
  for (int bid = 0; bid < bcount_i; bid++) {

    /* Get a hold of the ith bpart in ci. */
    struct bpart *restrict bi = &bparts_i[bid];

    /* Skip inactive particles */
    if (!bpart_is_active(bi, e)) continue;

    const float hi = bi->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float bix[3] = {(float)(bi->x[0] - (cj->loc[0] + shift[0])),
                          (float)(bi->x[1] - (cj->loc[1] + shift[1])),
                          (float)(bi->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the bparts in cj. */
    for (int bjd = 0; bjd < bcount_j; bjd++) {

      /* Get a pointer to the jth particle. */
      struct bpart *restrict bj = &bparts_j[bjd];
      const float hj = bj->h;

      /* Skip inhibited particles. */
      if (bpart_is_inhibited(bj, e)) continue;

      /* Compute the pairwise distance. */
      const float bjx[3] = {(float)(bj->x[0] - cj->loc[0]),
                            (float)(bj->x[1] - cj->loc[1]),
                            (float)(bj->x[2] - cj->loc[2])};
      const float dx[3] = {bix[0] - bjx[0], bix[1] - bjx[1], bix[2] - bjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (bi->ti_drift != e->ti_current)
        error("Particle bi not drifted to current time");
      if (bj->ti_drift != e->ti_current)
        error("Particle bj not drifted to current time");
#endif

      if (r2 < hig2) {
        IACT_BH_BH(r2, dx, hi, hj, bi, bj, cosmo, e->gravity_properties,
                   e->black_holes_properties, ti_current);

        if (bi_is_local) {
          runner_iact_nonsym_bh_bh_repos(r2, dx, hi, hj, bi, bj, cosmo,
                                         e->gravity_properties,
                                         e->black_holes_properties, ti_current);
        }
      }
    } /* loop over the bparts in cj. */
  }   /* loop over the bparts in ci. */

#endif /* (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW) */
}

void DOPAIR1_BH_NAIVE(struct runner *r, struct cell *restrict ci,
                      struct cell *restrict cj, int timer) {

  TIMER_TIC;

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_bh = ci->nodeID == r->e->nodeID;
  const int do_cj_bh = cj->nodeID == r->e->nodeID;
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
  /* here we are updating the hydro -> switch ci, cj */
  const int do_ci_bh = cj->nodeID == r->e->nodeID;
  const int do_cj_bh = ci->nodeID == r->e->nodeID;
#else
  /* The swallow task is executed on both sides */
  const int do_ci_bh = 1;
  const int do_cj_bh = 1;
#endif

  if (do_ci_bh) DO_NONSYM_PAIR1_BH_NAIVE(r, ci, cj);
  if (do_cj_bh) DO_NONSYM_PAIR1_BH_NAIVE(r, cj, ci);

  TIMER_TOC(TIMER_DOPAIR_BH);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts_i The #bpart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1_SUBSET_BH_NAIVE(struct runner *r, struct cell *restrict ci,
                             struct bpart *restrict bparts_i, int *restrict ind,
                             const int bcount, struct cell *restrict cj,
                             const double *shift) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int bi_is_local = ci->nodeID == e->nodeID;

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Loop over the parts_i. */
  for (int bid = 0; bid < bcount; bid++) {

    /* Get a hold of the ith part in ci. */
    struct bpart *restrict bi = &bparts_i[ind[bid]];

    const double bix = bi->x[0] - (shift[0]);
    const double biy = bi->x[1] - (shift[1]);
    const double biz = bi->x[2] - (shift[2]);
    const float hi = bi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!bpart_is_active(bi, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      struct xpart *restrict xpj = &xparts_j[pjd];

      /* Skip inhibited particles */
      if (part_is_inhibited(pj, e)) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float dx[3] = {(float)(bix - pjx), (float)(biy - pjy),
                           (float)(biz - pjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_BH_GAS(r2, dx, hi, hj, bi, pj, xpj, with_cosmology, cosmo,
                    e->gravity_properties, e->black_holes_properties,
                    e->entropy_floor, ti_current, e->time);
        if (bi_is_local) {
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
          runner_iact_nonsym_bh_gas_repos(
              r2, dx, hi, hj, bi, pj, xpj, with_cosmology, cosmo,
              e->gravity_properties, e->black_holes_properties,
              e->entropy_floor, ti_current, e->time);
#endif
        }
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts The #bpart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_BH(struct runner *r, struct cell *restrict ci,
                       struct bpart *restrict bparts, int *restrict ind,
                       const int bcount) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int bi_is_local = 1; /* SELF tasks are always local */

  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;
  struct xpart *restrict xparts_j = ci->hydro.xparts;

  /* Early abort? */
  if (count_i == 0) return;

  /* Loop over the parts in ci. */
  for (int bid = 0; bid < bcount; bid++) {

    /* Get a hold of the ith part in ci. */
    struct bpart *bi = &bparts[ind[bid]];
    const float bix[3] = {(float)(bi->x[0] - ci->loc[0]),
                          (float)(bi->x[1] - ci->loc[1]),
                          (float)(bi->x[2] - ci->loc[2])};
    const float hi = bi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!bpart_is_active(bi, e)) error("Inactive particle in subset function!");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_i; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];
      struct xpart *restrict xpj = &xparts_j[pjd];

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - ci->loc[0]),
                            (float)(pj->x[1] - ci->loc[1]),
                            (float)(pj->x[2] - ci->loc[2])};
      const float dx[3] = {bix[0] - pjx[0], bix[1] - pjx[1], bix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_BH_GAS(r2, dx, hi, pj->h, bi, pj, xpj, with_cosmology, cosmo,
                    e->gravity_properties, e->black_holes_properties,
                    e->entropy_floor, ti_current, e->time);

        if (bi_is_local) {
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
          runner_iact_nonsym_bh_gas_repos(
              r2, dx, hi, pj->h, bi, pj, xpj, with_cosmology, cosmo,
              e->gravity_properties, e->black_holes_properties,
              e->entropy_floor, ti_current, e->time);
#endif
        }
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

/**
 * @brief Determine which version of DOSELF1_SUBSET_BH needs to be called
 * depending on the optimisation level.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts The #bpart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_BRANCH_BH(struct runner *r, struct cell *restrict ci,
                              struct bpart *restrict bparts, int *restrict ind,
                              const int bcount) {

  DOSELF1_SUBSET_BH(r, ci, bparts, ind, bcount);
}

/**
 * @brief Determine which version of DOPAIR1_SUBSET_BH needs to be called
 * depending on the orientation of the cells or whether DOPAIR1_SUBSET_BH
 * needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param bparts_i The #bpart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param bcount The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR1_SUBSET_BRANCH_BH(struct runner *r, struct cell *restrict ci,
                              struct bpart *restrict bparts_i,
                              int *restrict ind, int const bcount,
                              struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (cj->hydro.count == 0) return;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  DOPAIR1_SUBSET_BH_NAIVE(r, ci, bparts_i, ind, bcount, cj, shift);
}

void DOSUB_SUBSET_BH(struct runner *r, struct cell *ci, struct bpart *bparts,
                     int *ind, const int bcount, struct cell *cj,
                     int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  /* Should we even bother? */
  if (!cell_is_active_black_holes(ci, e) &&
      (cj == NULL || !cell_is_active_black_holes(cj, e)))
    return;

  /* Find out in which sub-cell of ci the parts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&bparts[ind[0]] >= &ci->progeny[k]->black_holes.parts[0] &&
            &bparts[ind[0]] <
                &ci->progeny[k]
                     ->black_holes.parts[ci->progeny[k]->black_holes.count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_black_holes_task(ci)) {

      /* Loop over all progeny. */
      DOSUB_SUBSET_BH(r, sub, bparts, ind, bcount, NULL, 0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          DOSUB_SUBSET_BH(r, sub, bparts, ind, bcount, ci->progeny[j], 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF1_SUBSET_BRANCH_BH(r, ci, bparts, ind, bcount);
  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_black_holes_task(ci, cj) &&
        cell_can_recurse_in_pair_black_holes_task(cj, ci)) {

      /* Get the type of pair and flip ci/cj if needed. */
      double shift[3] = {0.0, 0.0, 0.0};
      const int sid = space_getsid(s, &ci, &cj, shift);

      struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
          DOSUB_SUBSET_BH(r, ci->progeny[pid], bparts, ind, bcount,
                          cj->progeny[pjd], 0);
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
          DOSUB_SUBSET_BH(r, cj->progeny[pjd], bparts, ind, bcount,
                          ci->progeny[pid], 0);
      }
    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_black_holes(ci, e) && cj->hydro.count > 0) {

      /* Do any of the cells need to be drifted first? */
      if (cell_is_active_black_holes(ci, e)) {
        if (!cell_are_bpart_drifted(ci, e)) error("Cell should be drifted!");
        if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");
      }

      DOPAIR1_SUBSET_BRANCH_BH(r, ci, bparts, ind, bcount, cj);
    }

  } /* otherwise, pair interaction. */
}

/**
 * @brief Determine which version of DOSELF1_BH needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_BH(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->black_holes.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_black_holes(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->black_holes.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  DOSELF1_BH(r, c, 1);
}

/**
 * @brief Determine which version of DOPAIR1_BH needs to be called depending
 * on the orientation of the cells or whether DOPAIR1_BH needs to be called
 * at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH_BH(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *restrict e = r->e;

  const int ci_active = cell_is_active_black_holes(ci, e);
  const int cj_active = cell_is_active_black_holes(cj, e);
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_bh = ci->nodeID == e->nodeID;
  const int do_cj_bh = cj->nodeID == e->nodeID;
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
  /* here we are updating the hydro -> switch ci, cj */
  const int do_ci_bh = cj->nodeID == e->nodeID;
  const int do_cj_bh = ci->nodeID == e->nodeID;
#else
  /* The swallow task is executed on both sides */
  const int do_ci_bh = 1;
  const int do_cj_bh = 1;
#endif

  const int do_ci = (ci->black_holes.count != 0 && cj->hydro.count != 0 &&
                     ci_active && do_ci_bh);
  const int do_cj = (cj->black_holes.count != 0 && ci->hydro.count != 0 &&
                     cj_active && do_cj_bh);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci &&
      (!cell_are_bpart_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj &&
      (!cell_are_part_drifted(ci, e) || !cell_are_bpart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* No sorted intreactions here -> use the naive ones */
  DOPAIR1_BH_NAIVE(r, ci, cj, 1);
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
void DOSUB_PAIR1_BH(struct runner *r, struct cell *ci, struct cell *cj,
                    int gettimer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother?
   * In the swallow case we care about BH-BH and BH-gas
   * interactions.
   * In all other cases only BH-gas so we can abort if there is
   * is no gas in the cell */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
  const int should_do_ci =
      ci->black_holes.count != 0 && cell_is_active_black_holes(ci, e);
  const int should_do_cj =
      cj->black_holes.count != 0 && cell_is_active_black_holes(cj, e);
#else
  const int should_do_ci = ci->black_holes.count != 0 && cj->hydro.count != 0 &&
                           cell_is_active_black_holes(ci, e);
  const int should_do_cj = cj->black_holes.count != 0 && ci->hydro.count != 0 &&
                           cell_is_active_black_holes(cj, e);

#endif

  if (!should_do_ci && !should_do_cj) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_black_holes_task(ci, cj) &&
      cell_can_recurse_in_pair_black_holes_task(cj, ci)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_BH(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
    const int do_ci_bh = ci->nodeID == e->nodeID;
    const int do_cj_bh = cj->nodeID == e->nodeID;
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
    /* Here we are updating the hydro -> switch ci, cj */
    const int do_ci_bh = cj->nodeID == e->nodeID;
    const int do_cj_bh = ci->nodeID == e->nodeID;
#else
    /* Here we perform the task on both sides */
    const int do_ci_bh = 1;
    const int do_cj_bh = 1;
#endif

    const int do_ci = ci->black_holes.count != 0 &&
                      cell_is_active_black_holes(ci, e) && do_ci_bh;
    const int do_cj = cj->black_holes.count != 0 &&
                      cell_is_active_black_holes(cj, e) && do_cj_bh;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_bpart_drifted(ci, e))
        error("Interacting undrifted cells (bparts).");

      if (cj->hydro.count != 0 && !cell_are_part_drifted(cj, e))
        error("Interacting undrifted cells (parts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (ci->hydro.count != 0 && !cell_are_part_drifted(ci, e))
        error("Interacting undrifted cells (parts).");

      if (!cell_are_bpart_drifted(cj, e))
        error("Interacting undrifted cells (bparts).");
    }

    if (do_ci || do_cj) DOPAIR1_BRANCH_BH(r, ci, cj);
  }

  TIMER_TOC(TIMER_DOSUB_PAIR_BH);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_BH(struct runner *r, struct cell *ci, int gettimer) {

  TIMER_TIC;

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

    /* Should we even bother?
     * In the swallow case we care about BH-BH and BH-gas
     * interactions.
     * In all other cases only BH-gas so we can abort if there is
     * is no gas in the cell */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
  const int should_do_ci =
      ci->black_holes.count != 0 && cell_is_active_black_holes(ci, e);
#else
  const int should_do_ci = ci->black_holes.count != 0 && ci->hydro.count != 0 &&
                           cell_is_active_black_holes(ci, e);
#endif

  if (!should_do_ci) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_black_holes_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_BH(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1_BH(r, ci->progeny[k], ci->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Check we did drift to the current time */
    if (!cell_are_bpart_drifted(ci, e)) error("Interacting undrifted cell.");

    if (ci->hydro.count != 0 && !cell_are_part_drifted(ci, e))
      error("Interacting undrifted cells (bparts).");

    DOSELF1_BRANCH_BH(r, ci);
  }

  TIMER_TOC(TIMER_DOSUB_SELF_BH);
}
