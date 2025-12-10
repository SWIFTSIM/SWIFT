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

#include "runner_doiact_stars_sidm.h"

/**
 * @brief Calculate the number density of #sipart around the #spart
 *
 * @param r runner task
 * @param c cell
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 * @param offset First particle in the cell to treat (for split tasks).
 * @param increment Interval between successive particles that are treated.
 */
void DOSELF1_STARS_SIDM(struct runner *r, const struct cell *c,
                   const int limit_min_h, const int limit_max_h) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (c->sidm.count == 0 || c->stars.count == 0) return;
  if (!cell_is_active_stars(c, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->stars.count;
  const int sicount = c->sidm.count;
  struct spart *restrict sparts = c->stars.parts;
  struct sipart *restrict siparts = c->sidm.parts;

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? c->depth : 0;
  const char max_depth = limit_min_h ? c->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? c->h_min_allowed : 0.;
  const float h_max = limit_max_h ? c->h_max_allowed : FLT_MAX;
#endif

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *si = &sparts[sid];
    const char depth_i = si->depth_h_sidm;
    const float hi = si->h_sidm;
    const float hig2 = hi * hi * kernel_gamma2;

    /* Skip inhibited particles */
    if (spart_is_inhibited(si, e)) continue;

    /* Skip inactive particles */
    if (!spart_is_active(si, e)) continue;

    /* Skip particles not in the range of h we care about */
    if (depth_i > max_depth) continue;
    if (depth_i < min_depth) continue;

    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the siparts in cj. */
    for (int sipjd = 0; sipjd < sicount; sipjd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts[sipjd];
      const float hj = sipj->h;

      /* Early abort? */
      if (sipart_is_inhibited(sipj, e)) continue;

      /* Compute the pairwise distance. */
      const float sipjx[3] = {(float)(sipj->x[0] - c->loc[0]),
                            (float)(sipj->x[1] - c->loc[1]),
                            (float)(sipj->x[2] - c->loc[2])};
      const float dx[3] = {six[0] - sipjx[0], six[1] - sipjx[1], six[2] - sipjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif

      if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_STARS_SIDM(r2, dx, hi, hj, si, sipj, a, H);

      }
    } /* loop over the siparts in ci. */
  } /* loop over the sparts in ci. */

  TIMER_TOC(TIMER_DOSELF_STARS_SIDM);
}

/**
 * @brief Calculate the number density of cj #sipart around the ci #spart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 * @param limit_min_h Only consider particles with h >= c->dmin/2.
 * @param limit_max_h Only consider particles with h < c->dmin.
 * @param offset First particle in the cell to treat (for split tasks).
 * @param increment Interval between successive particles that are treated.
 */
void DO_NONSYM_PAIR1_STARS_SIDM_NAIVE(struct runner *r,
                                 const struct cell *restrict ci,
                                 const struct cell *restrict cj,
                                 const int limit_min_h, const int limit_max_h) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (cj->sidm.count == 0 || ci->stars.count == 0) return;
  if (!cell_is_active_stars(ci, e)) return;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->stars.count;
  const int sicount_j = cj->sidm.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct sipart *restrict siparts_j = cj->sidm.parts;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->dmin != cj->dmin) error("Cells of different size!");
#endif

  /* Get the depth limits (if any) */
  const char min_depth = limit_max_h ? ci->depth : 0;
  const char max_depth = limit_min_h ? ci->depth : CHAR_MAX;

#ifdef SWIFT_DEBUG_CHECKS
  /* Get the limits in h (if any) */
  const float h_min = limit_min_h ? ci->h_min_allowed : 0.;
  const float h_max = limit_max_h ? ci->h_max_allowed : FLT_MAX;
#endif

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
    struct spart *si = &sparts_i[sid];
    const char depth_i = si->depth_h_sidm;

    /* Skip inhibited particles */
    if (spart_is_inhibited(si, e)) continue;

    /* Skip inactive particles */
    if (!spart_is_active(si, e)) continue;

    const float hi = si->h_sidm;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

#ifdef SWIFT_DEBUG_CHECKS
    if (hi > ci->stars.sidm.h_max_active)
      error("Particle has h larger than h_max_active");
#endif

    /* Skip particles not in the range of h we care about */
    if (depth_i > max_depth) continue;
    if (depth_i < min_depth) continue;

    /* Loop over the siparts in cj. */
    for (int sipjd = 0; sipjd < sicount_j; sipjd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts_j[sipjd];

      const float hj = sipj->h;

      /* Skip inhibited particles. */
      if (sipart_is_inhibited(sipj, e)) continue;

      /* Compute the pairwise distance. */
      const float sipjx[3] = {(float)(sipj->x[0] - cj->loc[0]),
                            (float)(sipj->x[1] - cj->loc[1]),
                            (float)(sipj->x[2] - cj->loc[2])};
      const float dx[3] = {six[0] - sipjx[0], six[1] - sipjx[1], six[2] - sipjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif

      if (r2 < hig2) {

#ifdef SWIFT_DEBUG_CHECKS
        if (hi < h_min || hi >= h_max) error("Inappropriate h for this level!");
#endif

        IACT_STARS_SIDM(r2, dx, hi, hj, si, sipj, a, H);

      }

    } /* loop over the siparts in cj. */
  } /* loop over the sparts in ci. */
}

void DOPAIR1_STARS_SIDM_NAIVE(struct runner *r, const struct cell *restrict ci,
                         const struct cell *restrict cj, const int limit_min_h,
                         const int limit_max_h) {

  TIMER_TIC;

  const int do_ci_stars = ci->nodeID == r->e->nodeID;
  const int do_cj_stars = cj->nodeID == r->e->nodeID;

  if (do_ci_stars && ci->stars.count != 0 && cj->sidm.count != 0)
    DO_NONSYM_PAIR1_STARS_SIDM_NAIVE(r, ci, cj, limit_min_h, limit_max_h);
  if (do_cj_stars && cj->stars.count != 0 && ci->sidm.count != 0)
    DO_NONSYM_PAIR1_STARS_SIDM_NAIVE(r, cj, ci, limit_min_h, limit_max_h);

  TIMER_TOC(TIMER_DOPAIR_STARS_SIDM);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts_i The #spart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1_SUBSET_STARS_SIDM_NAIVE(struct runner *r,
                                const struct cell *restrict ci,
                                struct spart *restrict sparts_i, const int *ind,
                                const int scount, struct cell *restrict cj,
                                const double shift[3]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int count_j = cj->sidm.count;
  struct sipart *restrict siparts_j = cj->sidm.parts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Loop over the sparts_i. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *restrict spi = &sparts_i[ind[sid]];

    const double pix = spi->x[0] - (shift[0]);
    const double piy = spi->x[1] - (shift[1]);
    const double piz = spi->x[2] - (shift[2]);
    const float hi = spi->h_sidm;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(spi, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the siparts in cj. */
    for (int sipjd = 0; sipjd < count_j; sipjd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts_j[sipjd];

      /* Skip inhibited particles */
      if (sipart_is_inhibited(sipj, e)) continue;

      const double sipjx = sipj->x[0];
      const double sipjy = sipj->x[1];
      const double sipjz = sipj->x[2];
      const float hj = sipj->h;

      /* Compute the pairwise distance. */
      float dx[3] = {(float)(pix - sipjx), (float)(piy - sipjy),
                     (float)(piz - sipjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sipj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS_SIDM(r2, dx, hi, hj, spi, sipj, a, H);
      }
    } /* loop over the siparts in cj. */
  } /* loop over the sparts in ci. */
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts The #spart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_STARS_SIDM(struct runner *r, const struct cell *ci,
                          struct spart *restrict sparts, const int *const ind,
                          const int scount) {
#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int sicount_i = ci->sidm.count;
  struct sipart *restrict siparts_j = ci->sidm.parts;

  /* Early abort? */
  if (sicount_i == 0) return;

  /* Loop over the sparts in ci. */
  for (int spid = 0; spid < scount; spid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *spi = &sparts[ind[spid]];
    const float spix[3] = {(float)(spi->x[0] - ci->loc[0]),
                           (float)(spi->x[1] - ci->loc[1]),
                           (float)(spi->x[2] - ci->loc[2])};
    const float hi = spi->h_sidm;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(spi, e))
      error("Inactive particle in subset function!");
#endif

    /* Loop over the siparts in cj. */
    for (int sipjd = 0; sipjd < sicount_i; sipjd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts_j[sipjd];

      /* Early abort? */
      if (sipart_is_inhibited(sipj, e)) continue;

      /* Compute the pairwise distance. */
      const float sipjx[3] = {(float)(sipj->x[0] - ci->loc[0]),
                            (float)(sipj->x[1] - ci->loc[1]),
                            (float)(sipj->x[2] - ci->loc[2])};
      float dx[3] = {spix[0] - sipjx[0], spix[1] - sipjx[1], spix[2] - sipjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sipj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS_SIDM(r2, dx, hi, sipj->h, spi, sipj, a, H);
      }
    } /* loop over the siparts in cj. */
  } /* loop over the sparts in ci. */
}

/**
 * @brief Determine which version of DOSELF1_SUBSET_STARS_SIDM needs to be called
 * depending on the optimisation level.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts The #spart to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_BRANCH_STARS_SIDM(struct runner *r, const struct cell *ci,
                                 struct spart *restrict sparts,
                                 const int *const ind, const int scount) {

  DOSELF1_SUBSET_STARS_SIDM(r, ci, sparts, ind, scount);
}

/**
 * @brief Determine which version of DOPAIR1_SUBSET_STARS_SIDM needs to be called
 * depending on the orientation of the cells or whether DOPAIR1_SUBSET_STARS_SIDM
 * needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sparts_i The #spart to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR1_SUBSET_BRANCH_STARS_SIDM(struct runner *r,
                                 const struct cell *restrict ci,
                                 struct spart *restrict sparts_i,
                                 const int *ind, const int scount,
                                 struct cell *restrict cj) {

  const struct engine *e = r->e;

  /* Anything to do here? */
  if (cj->sidm.count == 0) return;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  DOPAIR1_SUBSET_STARS_SIDM_NAIVE(r, ci, sparts_i, ind, scount, cj, shift);

}

/**
 * @brief Determine which version of DOSELF1_STARS_SIDM needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 * @param offset First particle in the cell to treat (for split tasks).
 * @param increment Interval between successive particles that are treated.
 */
void DOSELF1_BRANCH_STARS_SIDM(struct runner *r, const struct cell *c,
                          const int limit_min_h, const int limit_max_h) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->stars.count == 0) return;

  /* Anything to do here? */
  if (c->sidm.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e)) return;

#ifdef SWIFT_DEBUG_CHECKS

  /* Did we mess up the recursion? */
  if (c->stars.sidm.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  if (!limit_max_h && c->stars.sidm.h_max_active * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  /* Did we mess up the recursion? */
  if (limit_min_h && !limit_max_h)
    error("Fundamental error in the recursion logic");
#endif

  /* Check that cells are drifted. */
  if (!cell_are_part_drifted(c, e))
    error("Interacting undrifted cell (sidm).");
  if (!cell_are_spart_drifted(c, e))
    error("Interacting undrifted cell (stars).");

  DOSELF1_STARS_SIDM(r, c, limit_min_h, limit_max_h);
}

/**
 * @brief Determine which version of DOPAIR1_STARS_SIDM needs to be called depending
 * on the orientation of the cells or whether DOPAIR1_STARS_SIDM needs to be called
 * at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 * @param offset First particle in the cell to treat (for split tasks).
 * @param increment Interval between successive particles that are treated.
 */
void DOPAIR1_BRANCH_STARS_SIDM(struct runner *r, struct cell *ci, struct cell *cj,
                          const int limit_min_h, const int limit_max_h) {

  const struct engine *restrict e = r->e;

  /* Here we update the stars --> the star cell must be local */
  const int ci_local = (ci->nodeID == e->nodeID);
  const int cj_local = (cj->nodeID == e->nodeID);

  const int do_ci = ci->stars.count != 0 && cj->sidm.count != 0 &&
                    cell_is_active_stars(ci, e) && ci_local;
  const int do_cj = cj->stars.count != 0 && ci->sidm.count != 0 &&
                    cell_is_active_stars(cj, e) && cj_local;

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci &&
      (!cell_are_spart_drifted(ci, e) || !cell_are_sipart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj &&
      (!cell_are_sipart_drifted(ci, e) || !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  DOPAIR1_STARS_SIDM_NAIVE(r, ci, cj, limit_min_h, limit_max_h);

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
void DOSUB_PAIR1_STARS_SIDM(struct runner *r, struct cell *ci, struct cell *cj,
                       int recurse_below_h_max, const int gettimer) {
  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* Here we update the stars --> the star cell must be local */
  const int ci_local = (ci->nodeID == e->nodeID);
  const int cj_local = (cj->nodeID == e->nodeID);

  /* What kind of pair are we doing here? */
  const int do_ci = ci->stars.count != 0 && cj->sidm.count != 0 &&
                    cell_is_active_stars(ci, e) && ci_local;
  const int do_cj = cj->stars.count != 0 && ci->sidm.count != 0 &&
                    cell_is_active_stars(cj, e) && cj_local;

  /* Should we even bother? */
  if (!do_ci && !do_cj) return;

  /* We reached a leaf OR a cell small enough to be processed quickly */
  if (!ci->split || ci->stars.count < space_recurse_size_pair_stars ||
      !cj->split || cj->stars.count < space_recurse_size_pair_stars) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOPAIR1_BRANCH_STARS_SIDM(r, ci, cj, /*limit_h_min=*/0,
                         /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Both ci and cj are split */

    /* Should we change the recursion regime because we encountered a large
       particle? */
    if (!recurse_below_h_max && (!cell_can_recurse_in_subpair_stars_sidm_task(ci) ||
                                 !cell_can_recurse_in_subpair_stars_sidm_task(cj))) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* message("Multi-level PAIR! ci->count=%d cj->count=%d",
       * ci->sidm.count, cj->sidm.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOPAIR1_BRANCH_STARS_SIDM(r, ci, cj, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    const struct cell_split_pair *const csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL) {
        DOSUB_PAIR1_STARS_SIDM(r, ci->progeny[pid], cj->progeny[pjd],
                          recurse_below_h_max, /*gettimer=*/0);
      }
    }
  }

  TIMER_TOC(TIMER_DOSUB_PAIR_STARS_SIDM);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 * @param offset First particle in the cell to treat (for split tasks).
 * @param increment Interval between successive particles that are treated.
 */
void DOSUB_SELF1_STARS_SIDM(struct runner *r, struct cell *c,
                       int recurse_below_h_max, const int gettimer) {

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (c->sidm.count == 0 || c->stars.count == 0 ||
      !cell_is_active_stars(c, r->e))
    return;

  /* We reached a leaf OR a cell small enough to process quickly */
  if (!c->split || c->stars.count < space_recurse_size_self_stars) {

    /* We interact all particles in that cell:
       - No limit on the smallest h
       - Apply the max h limit if we are recursing below the level
       where h is smaller than the cell size */
    DOSELF1_BRANCH_STARS_SIDM(r, c, /*limit_h_min=*/0,
                         /*limit_h_max=*/recurse_below_h_max);

  } else {

    /* Should we change the recursion regime because we encountered a large
       particle at this level? */
    if (!recurse_below_h_max && !cell_can_recurse_in_subself_stars_sidm_task(c)) {
      recurse_below_h_max = 1;
    }

    /* If some particles are larger than the daughter cells, we must
       process them at this level before going deeper */
    if (recurse_below_h_max) {

      /* message("Multi-level SELF! c->count=%d", c->sidm.count); */

      /* Interact all *active* particles with h in the range [dmin/2, dmin)
         with all their neighbours */
      DOSELF1_BRANCH_STARS_SIDM(r, c, /*limit_h_min=*/1, /*limit_h_max=*/1);
    }

    /* Recurse to the lower levels. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        DOSUB_SELF1_STARS_SIDM(r, c->progeny[k], recurse_below_h_max,
                          /*gettimer=*/0);
        for (int j = k + 1; j < 8; j++) {
          if (c->progeny[j] != NULL) {
            DOSUB_PAIR1_STARS_SIDM(r, c->progeny[k], c->progeny[j],
                              recurse_below_h_max,
                              /*gettimer=*/0);
          }
        }
      }
    }
  }

  if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF_STARS_SIDM);
}

/**
 * @brief Find which sub-cell of a cell contains the subset of particles given
 * by the list of indices.
 *
 * Will throw an error if the sub-cell can't be found.
 *
 * @param c The #cell
 * @param sparts An array of #spart.
 * @param ind Index of the #spart's in the particle array to find in the subs.
 */
struct cell *FIND_SUB_STARS_SIDM(const struct cell *const c,
                            const struct spart *const sparts,
                            const int *const ind) {

#ifdef SWIFT_DEBUG_CHECKS
  if (!c->split) error("Can't search for subs in a non-split cell");
#endif

  /* Find out in which sub-cell of ci the sparts are.
   *
   * Note: We only need to check the first particle in the list */
  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL) {
      if (&sparts[ind[0]] >= &c->progeny[k]->stars.parts[0] &&
          &sparts[ind[0]] <
              &c->progeny[k]->stars.parts[c->progeny[k]->stars.count]) {
        return c->progeny[k];
        break;
      }
    }
  }
  error("Invalid sub!");
  return NULL;
}

void DOSUB_PAIR_SUBSET_STARS_SIDM(struct runner *r, struct cell *ci,
                             struct spart *sparts, const int *ind,
                             const int scount, struct cell *cj,
                             const int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  /* Should we even bother? */
  if (cj->sidm.count == 0) return;
  if (ci->stars.count == 0) return;
  if (!cell_is_active_stars(ci, e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_pair_stars_sidm_task(ci, cj) &&
      cell_can_recurse_in_pair_stars_sidm_task(cj, ci)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB_STARS_SIDM(ci, sparts, ind);

    /* Get the type of pair and flip ci/cj if needed. */
    double shift[3];
    const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
        DOSUB_PAIR_SUBSET_STARS_SIDM(r, ci->progeny[pid], sparts, ind, scount,
                                cj->progeny[pjd], /*gettimer=*/0);
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
        DOSUB_PAIR_SUBSET_STARS_SIDM(r, cj->progeny[pjd], sparts, ind, scount,
                                ci->progeny[pid], /*gettimer=*/0);
    }

  }
  /* Otherwise, compute the pair directly. */
  else if (cell_is_active_stars(ci, e)) {

    /* Do any of the cells need to be drifted first? */
    if (!cell_are_sipart_drifted(cj, e)) error("Cell should be drifted!");

    DOPAIR1_SUBSET_BRANCH_STARS_SIDM(r, ci, sparts, ind, scount, cj);
  }
}

void DOSUB_SELF_SUBSET_STARS_SIDM(struct runner *r, struct cell *ci,
                             struct spart *sparts, const int *ind,
                             const int scount, const int gettimer) {

  const struct engine *e = r->e;

  /* Should we even bother? */
  if (ci->sidm.count == 0) return;
  if (ci->stars.count == 0) return;
  if (!cell_is_active_stars(ci, e)) return;

  /* Recurse? */
  if (ci->split && cell_can_recurse_in_self_stars_sidm_task(ci)) {

    /* Find in which sub-cell of ci the particles are */
    struct cell *const sub = FIND_SUB_STARS_SIDM(ci, sparts, ind);

    /* Loop over all progeny. */
    DOSUB_SELF_SUBSET_STARS_SIDM(r, sub, sparts, ind, scount, /*gettimer=*/0);
    for (int j = 0; j < 8; j++)
      if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
        DOSUB_PAIR_SUBSET_STARS_SIDM(r, sub, sparts, ind, scount, ci->progeny[j],
                                /*gettimer=*/0);
  }

  /* Otherwise, compute self-interaction. */
  else
    DOSELF1_SUBSET_BRANCH_STARS_SIDM(r, ci, sparts, ind, scount);
}
