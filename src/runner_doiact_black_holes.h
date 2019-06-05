/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2018 Loic Hausammann (loic.hausammann@epfl.ch)
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
   runner_dopair_FUNCTION, runner_doself_FUNCTION and runner_dosub_FUNCTION
   calling the pairwise interaction function runner_iact_FUNCTION. */

#define PASTE(x, y) x##_##y

#define _DOSELF1_BH(f) PASTE(runner_doself_bh, f)
#define DOSELF1_BH _DOSELF1_BH(FUNCTION)

#define _DO_SYM_PAIR1_BH(f) PASTE(runner_do_sym_pair_bh, f)
#define DO_SYM_PAIR1_BH _DO_SYM_PAIR1_BH(FUNCTION)

#define _DO_NONSYM_PAIR1_BH_NAIVE(f) PASTE(runner_do_nonsym_pair_bh_naive, f)
#define DO_NONSYM_PAIR1_BH_NAIVE _DO_NONSYM_PAIR1_BH_NAIVE(FUNCTION)

#define _DOPAIR1_BH_NAIVE(f) PASTE(runner_dopair_bh_naive, f)
#define DOPAIR1_BH_NAIVE _DOPAIR1_BH_NAIVE(FUNCTION)

#define _DOPAIR1_SUBSET_BH(f) PASTE(runner_dopair_subset_bh, f)
#define DOPAIR1_SUBSET_BH _DOPAIR1_SUBSET_BH(FUNCTION)

#define _DOPAIR1_SUBSET_BH_NAIVE(f) PASTE(runner_dopair_subset_bh_naive, f)
#define DOPAIR1_SUBSET_BH_NAIVE _DOPAIR1_SUBSET_BH_NAIVE(FUNCTION)

#define _DOSELF1_SUBSET_BH(f) PASTE(runner_doself_subset_bh, f)
#define DOSELF1_SUBSET_BH _DOSELF1_SUBSET_BH(FUNCTION)

#define _DOSELF1_SUBSET_BRANCH_BH(f) PASTE(runner_doself_subset_branch_bh, f)
#define DOSELF1_SUBSET_BRANCH_BH _DOSELF1_SUBSET_BRANCH_BH(FUNCTION)

#define _DOPAIR1_SUBSET_BRANCH_BH(f) PASTE(runner_dopair_subset_branch_bh, f)
#define DOPAIR1_SUBSET_BRANCH_BH _DOPAIR1_SUBSET_BRANCH_BH(FUNCTION)

#define _DOSUB_SUBSET_BH(f) PASTE(runner_dosub_subset_bh, f)
#define DOSUB_SUBSET_BH _DOSUB_SUBSET_BH(FUNCTION)

#define _DOSELF1_BRANCH_BH(f) PASTE(runner_doself_branch_bh, f)
#define DOSELF1_BRANCH_BH _DOSELF1_BRANCH_BH(FUNCTION)

#define _DOPAIR1_BRANCH_BH(f) PASTE(runner_dopair_branch_bh, f)
#define DOPAIR1_BRANCH_BH _DOPAIR1_BRANCH_BH(FUNCTION)

#define _DOSUB_PAIR1_BH(f) PASTE(runner_dosub_pair_bh, f)
#define DOSUB_PAIR1_BH _DOSUB_PAIR1_BH(FUNCTION)

#define _DOSUB_SELF1_BH(f) PASTE(runner_dosub_self_bh, f)
#define DOSUB_SELF1_BH _DOSUB_SELF1_BH(FUNCTION)

#define _TIMER_DOSELF_BH(f) PASTE(timer_doself_bh, f)
#define TIMER_DOSELF_BH _TIMER_DOSELF_BH(FUNCTION)

#define _TIMER_DOPAIR_BH(f) PASTE(timer_dopair_bh, f)
#define TIMER_DOPAIR_BH _TIMER_DOPAIR_BH(FUNCTION)

#define _TIMER_DOSUB_SELF_BH(f) PASTE(timer_dosub_self_bh, f)
#define TIMER_DOSUB_SELF_BH _TIMER_DOSUB_SELF_BH(FUNCTION)

#define _TIMER_DOSUB_PAIR_BH(f) PASTE(timer_dosub_pair_bh, f)
#define TIMER_DOSUB_PAIR_BH _TIMER_DOSUB_PAIR_BH(FUNCTION)

#define _IACT_BH(f) PASTE(runner_iact_nonsym_bh, f)
#define IACT_BH _IACT_BH(FUNCTION)

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

  /* Anything to do here? */
  if (c->hydro.count == 0 || c->black_holes.count == 0) return;
  if (!cell_is_active_black_holes(c, e)) return;

  const int bcount = c->black_holes.count;
  const int count = c->hydro.count;
  struct bpart *restrict bparts = c->black_holes.parts;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  /* Loop over the bparts in ci. */
  for (int sid = 0; sid < bcount; sid++) {

    /* Get a hold of the ith bpart in ci. */
    struct bpart *restrict si = &bparts[sid];

    /* Skip inactive particles */
    if (!bpart_is_active(si, e)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

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
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      if (r2 < hig2) {
        IACT_BH(r2, dx, hi, hj, si, pj, xpj, cosmo, ti_current);
      }
    } /* loop over the parts in ci. */
  }   /* loop over the bparts in ci. */

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
#else
  if (cj->nodeID != engine_rank) error("Should be run on a different node");
#endif
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (cj->hydro.count == 0 || ci->black_holes.count == 0) return;
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

  /* Loop over the bparts in ci. */
  for (int sid = 0; sid < bcount_i; sid++) {

    /* Get a hold of the ith bpart in ci. */
    struct bpart *restrict si = &bparts_i[sid];

    /* Skip inactive particles */
    if (!bpart_is_active(si, e)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

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
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      if (r2 < hig2) {
        IACT_BH(r2, dx, hi, hj, si, pj, xpj, cosmo, ti_current);
      }
    } /* loop over the parts in cj. */
  }   /* loop over the parts in ci. */
}

void DOPAIR1_BH_NAIVE(struct runner *r, struct cell *restrict ci,
                      struct cell *restrict cj, int timer) {

  TIMER_TIC;

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_bh = ci->nodeID == r->e->nodeID;
  const int do_cj_bh = cj->nodeID == r->e->nodeID;
#else
  /* here we are updating the hydro -> switch ci, cj */
  const int do_ci_bh = cj->nodeID == r->e->nodeID;
  const int do_cj_bh = ci->nodeID == r->e->nodeID;
#endif
  if (do_ci_bh && ci->black_holes.count != 0 && cj->hydro.count != 0)
    DO_NONSYM_PAIR1_BH_NAIVE(r, ci, cj);
  if (do_cj_bh && cj->black_holes.count != 0 && ci->hydro.count != 0)
    DO_NONSYM_PAIR1_BH_NAIVE(r, cj, ci);

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

  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Loop over the parts_i. */
  for (int pid = 0; pid < bcount; pid++) {

    /* Get a hold of the ith part in ci. */
    struct bpart *restrict bpi = &bparts_i[ind[pid]];

    const double pix = bpi->x[0] - (shift[0]);
    const double piy = bpi->x[1] - (shift[1]);
    const double piz = bpi->x[2] - (shift[2]);
    const float hi = bpi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!bpart_is_active(bpi, e))
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
      float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                     (float)(piz - pjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_BH(r2, dx, hi, hj, bpi, pj, xpj, cosmo, ti_current);
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

  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;
  struct xpart *restrict xparts_j = ci->hydro.xparts;

  /* Early abort? */
  if (count_i == 0) return;

  /* Loop over the parts in ci. */
  for (int bpid = 0; bpid < bcount; bpid++) {

    /* Get a hold of the ith part in ci. */
    struct bpart *bpi = &bparts[ind[bpid]];
    const float bpix[3] = {(float)(bpi->x[0] - ci->loc[0]),
                           (float)(bpi->x[1] - ci->loc[1]),
                           (float)(bpi->x[2] - ci->loc[2])};
    const float hi = bpi->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!bpart_is_active(bpi, e))
      error("Inactive particle in subset function!");
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
      float dx[3] = {bpix[0] - pjx[0], bpix[1] - pjx[1], bpix[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_BH(r2, dx, hi, pj->h, bpi, pj, xpj, cosmo, ti_current);
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
#else
  /* here we are updating the hydro -> switch ci, cj */
  const int do_ci_bh = cj->nodeID == e->nodeID;
  const int do_cj_bh = ci->nodeID == e->nodeID;
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

  /* Should we even bother? */
  const int should_do_ci = ci->black_holes.count != 0 && cj->hydro.count != 0 &&
                           cell_is_active_black_holes(ci, e);
  const int should_do_cj = cj->black_holes.count != 0 && ci->hydro.count != 0 &&
                           cell_is_active_black_holes(cj, e);
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
#else
    /* here we are updating the hydro -> switch ci, cj */
    const int do_ci_bh = cj->nodeID == e->nodeID;
    const int do_cj_bh = ci->nodeID == e->nodeID;
#endif
    const int do_ci = ci->black_holes.count != 0 && cj->hydro.count != 0 &&
                      cell_is_active_black_holes(ci, e) && do_ci_bh;
    const int do_cj = cj->black_holes.count != 0 && ci->hydro.count != 0 &&
                      cell_is_active_black_holes(cj, e) && do_cj_bh;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_bpart_drifted(ci, e))
        error("Interacting undrifted cells (bparts).");

      if (!cell_are_part_drifted(cj, e))
        error("Interacting undrifted cells (parts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_part_drifted(ci, e))
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

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (ci->hydro.count == 0 || ci->black_holes.count == 0 ||
      !cell_is_active_black_holes(ci, r->e))
    return;

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

    /* Drift the cell to the current timestep if needed. */
    if (!cell_are_bpart_drifted(ci, r->e)) error("Interacting undrifted cell.");

    DOSELF1_BRANCH_BH(r, ci);
  }

  TIMER_TOC(TIMER_DOSUB_SELF_BH);
}
