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
 * @brief Calculate the number density of #part around the #bpart
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_STARS_SIDM(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int si_is_local = 1; /* SELF tasks are always local */

  /* Anything to do here? */
  if (c->stars.count == 0) return;
  if (!cell_is_active_stars(c, e)) return;

  const int scount = c->stars.count;
  const int count = c->sidm.count;
  struct spart *restrict sparts = c->stars.parts;
  struct sipart *restrict siparts = c->sidm.parts;

  /* Do we actually have any sidm neighbours? */
  if (count != 0) {

    /* Loop over the sparts in ci. */
    for (int sid = 0; sid < scount; sid++) {

      /* Get a hold of the ith spart in ci. */
      struct spart *restrict si = &sparts[sid];

      /* Skip inactive particles */
      if (!spart_is_active(si, e)) continue;

      const float hi = si->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float six[3] = {(float)(si->x[0] - c->loc[0]),
                            (float)(si->x[1] - c->loc[1]),
                            (float)(si->x[2] - c->loc[2])};

      /* Loop over the siparts in cj. */
      for (int sipjd = 0; sipjd < count; sipjd++) {

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
        if (si->ti_drift != e->ti_current)
          error("Particle si not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        if (r2 < hig2) {
          IACT_STARS_SIDM(r2, dx, hi, hj, si, sipj, with_cosmology, cosmo,
                      e->gravity_properties, e->stars_properties, 
                      ti_current, e->time);
        }
      } /* loop over the siparts in ci. */
    } /* loop over the sparts in ci. */
  } /* Do we have star particles in the cell? */

  TIMER_TOC(TIMER_DOSELF_STARS_SIDM);
}

/**
 * @brief Calculate the number density of cj #part around the ci #bpart
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_STARS_SIDM_NAIVE(struct runner *r, struct cell *restrict ci,
                              struct cell *restrict cj) {

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int si_is_local = ci->nodeID == e->nodeID;

  /* Anything to do here? */
  if (ci->stars.count == 0) return;
  if (!cell_is_active_stars(ci, e)) return;

  const int scount_i = ci->stars.count;
  const int count_j = cj->sidm.count;
  struct spart *restrict sparts_i = ci->stars.parts;
  struct sipart *restrict siparts_j = cj->sidm.parts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Do we actually have any sidm neighbours? */
  if (cj->sidm.count != 0) {

    /* Loop over the sparts in ci. */
    for (int sid = 0; sid < scount_i; sid++) {

      /* Get a hold of the ith bpart in ci. */
      struct spart *restrict si = &sparts_i[sid];

      /* Skip inactive particles */
      if (!spart_is_active(si, e)) continue;

      const float hi = si->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                            (float)(si->x[1] - (cj->loc[1] + shift[1])),
                            (float)(si->x[2] - (cj->loc[2] + shift[2]))};

      /* Loop over the siparts in cj. */
      for (int sipjd = 0; sipjd < count_j; sipjd++) {

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
        if (si->ti_drift != e->ti_current)
          error("Particle si not drifted to current time");
        if (sipj->ti_drift != e->ti_current)
          error("Particle sipj not drifted to current time");
#endif

        if (r2 < hig2) {
          IACT_STARS_SIDM(r2, dx, hi, hj, si, sipj, with_cosmology, cosmo,
                      e->gravity_properties, e->stars_properties, ti_current, e->time);
        }
      } /* loop over the siparts in cj. */
    } /* loop over the sparts in ci. */
  } /* Do we have gas particles in the cell? */
}

void DOPAIR1_STARS_SIDM_NAIVE(struct runner *r, struct cell *restrict ci,
                      struct cell *restrict cj, int timer) {

  TIMER_TIC;

  const int do_ci_bh = ci->nodeID == r->e->nodeID;
  const int do_cj_bh = cj->nodeID == r->e->nodeID;

  if (do_ci_bh) DO_NONSYM_PAIR1_STARS_SIDM_NAIVE(r, ci, cj);
  if (do_cj_bh) DO_NONSYM_PAIR1_STARS_SIDM_NAIVE(r, cj, ci);

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
void DOPAIR1_SUBSET_STARS_SIDM_NAIVE(struct runner *r, struct cell *restrict ci,
                             struct spart *restrict sparts_i, int *restrict ind,
                             const int scount, struct cell *restrict cj,
                             const double *shift) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int si_is_local = ci->nodeID == e->nodeID;

  const int count_j = cj->sidm.count;
  struct sipart *restrict siparts_j = cj->sidm.parts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Loop over the siparts_i. */
  for (int sid = 0; sid < bcount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *restrict si = &sparts_i[ind[sid]];

    const double six = si->x[0] - (shift[0]);
    const double siy = si->x[1] - (shift[1]);
    const double siz = si->x[2] - (shift[2]);
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(si, e))
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
      const float dx[3] = {(float)(six - sipjx), (float)(siy - sipjy),
                           (float)(siz - sipjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS_SIDM(r2, dx, hi, hj, si, sipj, with_cosmology, cosmo,
                    e->gravity_properties, e->stars_properties, ti_current, e->time);
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
void DOSELF1_SUBSET_STARS_SIDM(struct runner *r, struct cell *restrict ci,
                       struct spart *restrict sparts, int *restrict ind,
                       const int scount) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int si_is_local = 1; /* SELF tasks are always local */

  const int count_i = ci->stars.count;
  struct sipart *restrict siparts_j = ci->sidm.parts;

  /* Early abort? */
  if (count_i == 0) return;

  /* Loop over the sparts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith part in ci. */
    struct spart *si = &sparts[ind[sid]];
    const float six[3] = {(float)(si->x[0] - ci->loc[0]),
                          (float)(si->x[1] - ci->loc[1]),
                          (float)(si->x[2] - ci->loc[2])};
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!spart_is_active(si, e)) error("Inactive particle in subset function!");
#endif

    /* Loop over the siparts in cj. */
    for (int sipjd = 0; sipjd < count_i; sipjd++) {

      /* Get a pointer to the jth particle. */
      struct sipart *restrict sipj = &siparts_j[sipjd];

      /* Early abort? */
      if (sipart_is_inhibited(sipj, e)) continue;

      /* Compute the pairwise distance. */
      const float sipjx[3] = {(float)(sipj->x[0] - ci->loc[0]),
                            (float)(sipj->x[1] - ci->loc[1]),
                            (float)(sipj->x[2] - ci->loc[2])};
      const float dx[3] = {six[0] - sipjx[0], six[1] - sipjx[1], six[2] - sipjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sipj->ti_drift != e->ti_current)
        error("Particle sipj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_STARS_SIDM(r2, dx, hi, pj->h, si, sipj, with_cosmology, cosmo,
                    e->gravity_properties, e->stars_properties, ti_current, e->time);
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
void DOSELF1_SUBSET_BRANCH_STARS_SIDM(struct runner *r, struct cell *restrict ci,
                              struct spart *restrict sparts, int *restrict ind,
                              const int scount) {

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
void DOPAIR1_SUBSET_BRANCH_STARS_SIDM(struct runner *r, struct cell *restrict ci,
                              struct spart *restrict sparts_i,
                              int *restrict ind, int const scount,
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

void DOSUB_SUBSET_STARS_SIDM(struct runner *r, struct cell *ci, struct spart *sparts,
                     int *ind, const int scount, struct cell *cj,
                     int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  /* Should we even bother? */
  if (!cell_is_active_stars(ci, e) &&
      (cj == NULL || !cell_is_active_stars(cj, e)))
    return;

  /* Find out in which sub-cell of ci the sparts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&sparts[ind[0]] >= &ci->progeny[k]->stars.parts[0] &&
            &sparts[ind[0]] <
                &ci->progeny[k]
                     ->stars.parts[ci->progeny[k]->stars.count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_stars_sidm_task(ci)) {

      /* Loop over all progeny. */
      DOSUB_SUBSET_STARS_SIDM(r, sub, sparts, ind, scount, NULL, 0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          DOSUB_SUBSET_STARS_SIDM(r, sub, sparts, ind, scount, ci->progeny[j], 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF1_SUBSET_BRANCH_STARS_SIDM(r, ci, sparts, ind, scount);
  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_stars_sidm_task(ci, cj) &&
        cell_can_recurse_in_pair_stars_sidm_task(cj, ci)) {

      /* Get the type of pair and flip ci/cj if needed. */
      double shift[3] = {0.0, 0.0, 0.0};
      const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

      struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
          DOSUB_SUBSET_STARS_SIDM(r, ci->progeny[pid], sparts, ind, scount,
                          cj->progeny[pjd], 0);
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
          DOSUB_SUBSET_STARS_SIDM(r, cj->progeny[pjd], sparts, ind, scount,
                          ci->progeny[pid], 0);
      }
    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_stars(ci, e) && cj->sidm.count > 0) {

      /* Do any of the cells need to be drifted first? */
      if (cell_is_active_stars(ci, e)) {
        if (!cell_are_spart_drifted(ci, e)) error("Cell should be drifted!");
        if (!cell_are_sipart_drifted(cj, e)) error("Cell should be drifted!");
      }

      DOPAIR1_SUBSET_BRANCH_STARS_SIDM(r, ci, sparts, ind, scount, cj);
    }

  } /* otherwise, pair interaction. */
}

/**
 * @brief Determine which version of DOSELF1_STARS_SIDM needs to be called depending
 * on the optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_STARS_SIDM(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->stars.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_stars(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->stars.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than smoothing length");

  DOSELF1_STARS_SIDM(r, c, 1);
}

/**
 * @brief Determine which version of DOPAIR1_STARS_SIDM needs to be called depending
 * on the orientation of the cells or whether DOPAIR1_STARS_SIDM needs to be called
 * at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH_STARS_SIDM(struct runner *r, struct cell *ci, struct cell *cj) {

  const struct engine *restrict e = r->e;

  const int ci_active = cell_is_active_stars(ci, e);
  const int cj_active = cell_is_active_stars(cj, e);
  
  const int do_ci_stars = ci->nodeID == e->nodeID;
  const int do_cj_stars = cj->nodeID == e->nodeID;

  const int do_ci = (ci->stars.count != 0 && cj->sidm.count != 0 &&
                     ci_active && do_ci_stars);
  const int do_cj = (cj->stars.count != 0 && ci->sidm.count != 0 &&
                     cj_active && do_cj_stars);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci &&
      (!cell_are_spart_drifted(ci, e) || !cell_are_sipart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj &&
      (!cell_are_sipart_drifted(ci, e) || !cell_are_spart_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* No sorted intreactions here -> use the naive ones */
  DOPAIR1_STARS_SIDM_NAIVE(r, ci, cj, 1);
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
                    int gettimer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  const int should_do_ci = ci->stars.count != 0 && cj->sidm.count != 0 &&
                           cell_is_active_stars(ci, e);
  const int should_do_cj = cj->stars.count != 0 && ci->sidm.count != 0 &&
                           cell_is_active_stars(cj, e);

  if (!should_do_ci && !should_do_cj) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_stars_sidm_task(ci, cj) &&
      cell_can_recurse_in_pair_stars_sidm_task(cj, ci)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_STARS_SIDM(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    const int do_ci_stars = ci->nodeID == e->nodeID;
    const int do_cj_stars = cj->nodeID == e->nodeID;

    const int do_ci = ci->stars.count != 0 &&
                      cell_is_active_stars(ci, e) && do_ci_stars;
    const int do_cj = cj->stars.count != 0 &&
                      cell_is_active_stars(cj, e) && do_cj_stars;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_spart_drifted(ci, e))
        error("Interacting undrifted cells (sparts).");

      if (cj->sidm.count != 0 && !cell_are_sipart_drifted(cj, e))
        error("Interacting undrifted cells (siparts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (ci->sidm.count != 0 && !cell_are_sipart_drifted(ci, e))
        error("Interacting undrifted cells (siparts).");

      if (!cell_are_spart_drifted(cj, e))
        error("Interacting undrifted cells (sparts).");
    }

    if (do_ci || do_cj) DOPAIR1_BRANCH_STARS_SIDM(r, ci, cj);
  }

  TIMER_TOC(TIMER_DOSUB_PAIR_STARS_SIDM);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_STARS_SIDM(struct runner *r, struct cell *ci, int gettimer) {

  TIMER_TIC;

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

    /* Should we even bother? */
  const int should_do_ci = ci->stars.count != 0 && ci->sidm.count != 0 &&
                           cell_is_active_stars(ci, e);

  if (!should_do_ci) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_stars_sidm_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_STARS_SIDM(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1_STARS_SIDM(r, ci->progeny[k], ci->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Check we did drift to the current time */
    if (!cell_are_spart_drifted(ci, e)) error("Interacting undrifted cell.");

    if (ci->sidm.count != 0 && !cell_are_sipart_drifted(ci, e))
      error("Interacting undrifted cells (siparts).");

    DOSELF1_BRANCH_STARS_SIDM(r, ci);
  }

  TIMER_TOC(TIMER_DOSUB_SELF_STARS_SIDM);
}
