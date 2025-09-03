/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Jonathan Davies (j.j.davies@ljmu.ac.uk)
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

#include "runner_doiact_sinks.h"

/**
 * @brief Calculate gas and sink interaction around #sink
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_SINKS(struct runner *r, struct cell *c, int timer) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  TIMER_TIC;

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;

  /* Anything to do here? */
  if (c->sinks.count == 0) return;
  if (!cell_is_active_sinks(c, e)) return;

  const int scount = c->sinks.count;
  const int count = c->hydro.count;
  struct sink *restrict sinks = c->sinks.parts;
  struct part *restrict parts = c->hydro.parts;

  /* Do we actually have any gas neighbours? */
  if (c->hydro.count != 0) {

    /* Loop over the sinks in ci. */
    for (int sid = 0; sid < scount; sid++) {

      /* Get a hold of the ith sinks in ci. */
      struct sink *restrict si = &sinks[sid];

      /* Skip inactive particles */
      if (!sink_is_active(si, e)) continue;

      const float hi = si->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float six[3] = {(float)(si->x[0] - c->loc[0]),
                            (float)(si->x[1] - c->loc[1]),
                            (float)(si->x[2] - c->loc[2])};

      /* Loop over the parts (gas) in cj. */
      for (int pjd = 0; pjd < count; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts[pjd];
        const float hj = pj->h;

        /* Early abort? */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                              (float)(pj->x[1] - c->loc[1]),
                              (float)(pj->x[2] - c->loc[2])};
        const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (si->ti_drift != e->ti_current)
          error("Particle si not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        if (r2 < hig2) {
          IACT_SINKS_GAS(r2, dx, hi, hj, si, pj, with_cosmology, cosmo,
                         e->gravity_properties, e->sink_properties,
                         e->ti_current, e->time);
        }
      } /* loop over the parts in ci. */
    } /* loop over the sinks in ci. */
  } /* Do we have gas particles in the cell? */

  /* When doing sink swallowing, we need a quick loop also over the sink
   * neighbours */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)

  /* Loop over the sinks in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith sink in ci. */
    struct sink *restrict si = &sinks[sid];

    /* Skip inactive particles */
    if (!sink_is_active(si, e)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the sinks in cj. */
    for (int sjd = 0; sjd < scount; sjd++) {

      /* Skip self interaction */
      if (sid == sjd) continue;

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks[sjd];
      const float hj = sj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Early abort? */
      if (sink_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - c->loc[0]),
                            (float)(sj->x[1] - c->loc[1]),
                            (float)(sj->x[2] - c->loc[2])};
      const float dx[3] = {six[0] - sjx[0], six[1] - sjx[1], six[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (si->ti_drift != e->ti_current)
        error("Particle si not drifted to current time");
      if (sj->ti_drift != e->ti_current)
        error("Particle bj not drifted to current time");
#endif

      if (r2 < hig2 || r2 < hjg2) {
        IACT_SINKS_SINK(r2, dx, hi, hj, si, sj, with_cosmology, cosmo,
                        e->gravity_properties, e->sink_properties,
                        e->ti_current, e->time);
      }
    } /* loop over the sinks in ci. */
  } /* loop over the sinks in ci. */

#endif /* (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW) */

  if (timer) TIMER_TOC(TIMER_DOSELF_SINKS);
}

/**
 * @brief Calculate gas and sink interaction around #sink
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_SINKS_NAIVE(struct runner *r, struct cell *restrict ci,
                                 struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;

  /* Anything to do here? */
  if (ci->sinks.count == 0) return;
  if (!cell_is_active_sinks(ci, e)) return;

  const int scount_i = ci->sinks.count;
  const int count_j = cj->hydro.count;
  struct sink *restrict sinks_i = ci->sinks.parts;
  struct part *restrict parts_j = cj->hydro.parts;

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

    /* Loop over the sinks in ci. */
    for (int sid = 0; sid < scount_i; sid++) {

      /* Get a hold of the ith sink in ci. */
      struct sink *restrict si = &sinks_i[sid];

      /* Skip inactive particles */
      if (!sink_is_active(si, e)) continue;

      const float hi = si->h;
      const float hig2 = hi * hi * kernel_gamma2;
      const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                            (float)(si->x[1] - (cj->loc[1] + shift[1])),
                            (float)(si->x[2] - (cj->loc[2] + shift[2]))};

      /* Loop over the parts (gas) in cj. */
      for (int pjd = 0; pjd < count_j; pjd++) {

        /* Get a pointer to the jth particle. */
        struct part *restrict pj = &parts_j[pjd];
        const float hj = pj->h;

        /* Skip inhibited particles. */
        if (part_is_inhibited(pj, e)) continue;

        /* Compute the pairwise distance. */
        const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                              (float)(pj->x[1] - cj->loc[1]),
                              (float)(pj->x[2] - cj->loc[2])};
        const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
        const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (si->ti_drift != e->ti_current)
          error("Particle si not drifted to current time");
        if (pj->ti_drift != e->ti_current)
          error("Particle pj not drifted to current time");
#endif

        if (r2 < hig2) {
          IACT_SINKS_GAS(r2, dx, hi, hj, si, pj, with_cosmology, cosmo,
                         e->gravity_properties, e->sink_properties,
                         e->ti_current, e->time);
        }
      } /* loop over the parts in cj. */
    } /* loop over the sinks in ci. */
  } /* Do we have gas particles in the cell? */

  /* When doing sink swallowing, we need a quick loop also over the sinks
   * neighbours */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)

  const int scount_j = cj->sinks.count;
  struct sink *restrict sinks_j = cj->sinks.parts;

  /* Loop over the sinks in ci. */
  for (int sid = 0; sid < scount_i; sid++) {

    /* Get a hold of the ith sink in ci. */
    struct sink *restrict si = &sinks_i[sid];

    /* Skip inactive particles */
    if (!sink_is_active(si, e)) continue;

    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the sinks in cj. */
    for (int sjd = 0; sjd < scount_j; sjd++) {

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks_j[sjd];
      const float hj = sj->h;
      const float hjg2 = hj * hj * kernel_gamma2;

      /* Skip inhibited particles. */
      if (sink_is_inhibited(sj, e)) continue;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - cj->loc[0]),
                            (float)(sj->x[1] - cj->loc[1]),
                            (float)(sj->x[2] - cj->loc[2])};
      const float dx[3] = {six[0] - sjx[0], six[1] - sjx[1], six[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (si->ti_drift != e->ti_current)
        error("Particle si not drifted to current time");
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
#endif

      if (r2 < hig2 || r2 < hjg2) {
        IACT_SINKS_SINK(r2, dx, hi, hj, si, sj, with_cosmology, cosmo,
                        e->gravity_properties, e->sink_properties,
                        e->ti_current, e->time);
      }
    } /* loop over the sinks in cj. */
  } /* loop over the sinks in ci. */

#endif /* (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW) */
}

/**
 * @brief Calculate swallow for ci #sink part around the cj #part and sinks and
 *                              cj #sink part around the ci #part and sinks
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DOPAIR1_SINKS_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct cell *restrict cj, int timer) {

  TIMER_TIC;

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_sink = ci->nodeID == r->e->nodeID;
  const int do_cj_sink = cj->nodeID == r->e->nodeID;
#else
  /* The swallow task is executed on both sides */
  const int do_ci_sink = 1;
  const int do_cj_sink = 1;
#endif

  if (do_ci_sink) DO_NONSYM_PAIR1_SINKS_NAIVE(r, ci, cj);
  if (do_cj_sink) DO_NONSYM_PAIR1_SINKS_NAIVE(r, cj, ci);

  if (timer) TIMER_TOC(TIMER_DOPAIR_SINKS);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sinks_i The #sink to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1_SUBSET_SINKS_NAIVE(struct runner *r, struct cell *restrict ci,
                                struct sink *restrict sinks_i,
                                int *restrict ind, const int scount,
                                struct cell *restrict cj, const double *shift) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int count_j = cj->hydro.count;
  struct part *restrict parts_j = cj->hydro.parts;

  /* Early abort? */
  if (count_j == 0) return;

  /* Loop over the parts_i. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith part in ci. */
    struct sink *restrict si = &sinks_i[ind[sid]];

    const double six = si->x[0] - (shift[0]);
    const double siy = si->x[1] - (shift[1]);
    const double siz = si->x[2] - (shift[2]);
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!sink_is_active(si, e))
      error("Trying to correct smoothing length of inactive particle !");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_j; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Skip inhibited particles */
      if (part_is_inhibited(pj, e)) continue;

      const double pjx = pj->x[0];
      const double pjy = pj->x[1];
      const double pjz = pj->x[2];
      const float hj = pj->h;

      /* Compute the pairwise distance. */
      const float dx[3] = {(float)(six - pjx), (float)(siy - pjy),
                           (float)(siz - pjz)};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif
      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_SINKS_GAS(r2, dx, hi, hj, si, pj, with_cosmology, cosmo,
                       e->gravity_properties, e->sink_properties, e->ti_current,
                       e->time);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sinks The #sink to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_SINKS(struct runner *r, struct cell *restrict ci,
                          struct sink *restrict sinks, int *restrict ind,
                          const int scount) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const int count_i = ci->hydro.count;
  struct part *restrict parts_j = ci->hydro.parts;

  /* Early abort? */
  if (count_i == 0) return;

  /* Loop over the parts in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith part in ci. */
    struct sink *si = &sinks[ind[sid]];
    const float six[3] = {(float)(si->x[0] - ci->loc[0]),
                          (float)(si->x[1] - ci->loc[1]),
                          (float)(si->x[2] - ci->loc[2])};
    const float hi = si->h;
    const float hig2 = hi * hi * kernel_gamma2;

#ifdef SWIFT_DEBUG_CHECKS
    if (!sink_is_active(si, e)) error("Inactive particle in subset function!");
#endif

    /* Loop over the parts in cj. */
    for (int pjd = 0; pjd < count_i; pjd++) {

      /* Get a pointer to the jth particle. */
      struct part *restrict pj = &parts_j[pjd];

      /* Early abort? */
      if (part_is_inhibited(pj, e)) continue;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(pj->x[0] - ci->loc[0]),
                            (float)(pj->x[1] - ci->loc[1]),
                            (float)(pj->x[2] - ci->loc[2])};
      const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (pj->ti_drift != e->ti_current)
        error("Particle pj not drifted to current time");
#endif

      /* Hit or miss? */
      if (r2 < hig2) {
        IACT_SINKS_GAS(r2, dx, hi, pj->h, si, pj, with_cosmology, cosmo,
                       e->gravity_properties, e->sink_properties, e->ti_current,
                       e->time);
      }
    } /* loop over the parts in cj. */
  } /* loop over the parts in ci. */
}

/**
 * @brief Determine which version of DOSELF1_SUBSET_SINKS needs to be called
 * depending on the optimisation level.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sinks The #sink to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 */
void DOSELF1_SUBSET_BRANCH_SINKS(struct runner *r, struct cell *restrict ci,
                                 struct sink *restrict sinks, int *restrict ind,
                                 const int scount) {

  DOSELF1_SUBSET_SINKS(r, ci, sinks, ind, scount);
}

/**
 * @brief Determine which version of DOPAIR1_SUBSET_SINKS needs to be called
 * depending on the orientation of the cells or whether DOPAIR1_SUBSET_SINKS
 * needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param sinks_i The #sink to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param scount The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR1_SUBSET_BRANCH_SINKS(struct runner *r, struct cell *restrict ci,
                                 struct sink *restrict sinks_i,
                                 int *restrict ind, int const scount,
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

  DOPAIR1_SUBSET_SINKS_NAIVE(r, ci, sinks_i, ind, scount, cj, shift);
}

void DOSUB_SUBSET_SINKS(struct runner *r, struct cell *ci, struct sink *sinks,
                        int *ind, const int scount, struct cell *cj,
                        int gettimer) {

  const struct engine *e = r->e;
  struct space *s = e->s;

  /* Should we even bother? */
  if (!cell_is_active_sinks(ci, e) &&
      (cj == NULL || !cell_is_active_sinks(cj, e)))
    return;

  /* Find out in which sub-cell of ci the parts are. */
  struct cell *sub = NULL;
  if (ci->split) {
    for (int k = 0; k < 8; k++) {
      if (ci->progeny[k] != NULL) {
        if (&sinks[ind[0]] >= &ci->progeny[k]->sinks.parts[0] &&
            &sinks[ind[0]] <
                &ci->progeny[k]->sinks.parts[ci->progeny[k]->sinks.count]) {
          sub = ci->progeny[k];
          break;
        }
      }
    }
  }

  /* Is this a single cell? */
  if (cj == NULL) {

    /* Recurse? */
    if (cell_can_recurse_in_self_sinks_task(ci)) {

      /* Loop over all progeny. */
      DOSUB_SUBSET_SINKS(r, sub, sinks, ind, scount, NULL, 0);
      for (int j = 0; j < 8; j++)
        if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
          DOSUB_SUBSET_SINKS(r, sub, sinks, ind, scount, ci->progeny[j], 0);

    }

    /* Otherwise, compute self-interaction. */
    else
      DOSELF1_SUBSET_BRANCH_SINKS(r, ci, sinks, ind, scount);
  } /* self-interaction. */

  /* Otherwise, it's a pair interaction. */
  else {

    /* Recurse? */
    if (cell_can_recurse_in_pair_sinks_task(ci, cj) &&
        cell_can_recurse_in_pair_sinks_task(cj, ci)) {

      /* Get the type of pair and flip ci/cj if needed. */
      double shift[3] = {0.0, 0.0, 0.0};
      const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

      struct cell_split_pair *csp = &cell_split_pairs[sid];
      for (int k = 0; k < csp->count; k++) {
        const int pid = csp->pairs[k].pid;
        const int pjd = csp->pairs[k].pjd;
        if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
          DOSUB_SUBSET_SINKS(r, ci->progeny[pid], sinks, ind, scount,
                             cj->progeny[pjd], 0);
        if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
          DOSUB_SUBSET_SINKS(r, cj->progeny[pjd], sinks, ind, scount,
                             ci->progeny[pid], 0);
      }
    }

    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_sinks(ci, e) && cj->hydro.count > 0) {

      /* Do any of the cells need to be drifted first? */
      if (cell_is_active_sinks(ci, e)) {
        if (!cell_are_sink_drifted(ci, e)) error("Cell should be drifted!");
        if (!cell_are_part_drifted(cj, e)) error("Cell should be drifted!");
      }

      DOPAIR1_SUBSET_BRANCH_SINKS(r, ci, sinks, ind, scount, cj);
    }

  } /* otherwise, pair interaction. */
}

/**
 * @brief Wrapper to runner_doself_sinks_swallow
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_SINKS(struct runner *r, struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->sinks.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_sinks(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->sinks.h_max_old * kernel_gamma > c->dmin)
    error("Cell smaller than the cut off radius or smoothing length");

  DOSELF1_SINKS(r, c, 1);
}

/**
 * @brief Wrapper for runner_dopair_sinks_naive_swallow.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH_SINKS(struct runner *r, struct cell *ci, struct cell *cj) {

#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

  const struct engine *restrict e = r->e;

  const int ci_active = cell_is_active_sinks(ci, e);
  const int cj_active = cell_is_active_sinks(cj, e);

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
  const int do_ci_sink = ci->nodeID == e->nodeID;
  const int do_cj_sink = cj->nodeID == e->nodeID;
#else
  /* The swallow task is executed on both sides */
  const int do_ci_sink = 1;
  const int do_cj_sink = 1;
#endif

  const int do_ci =
      (ci->sinks.count != 0 && cj->hydro.count != 0 && ci_active && do_ci_sink);
  const int do_cj =
      (cj->sinks.count != 0 && ci->hydro.count != 0 && cj_active && do_cj_sink);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci && (!cell_are_sink_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj && (!cell_are_part_drifted(ci, e) || !cell_are_sink_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* No sorted interactions here -> use the naive ones */
  DOPAIR1_SINKS_NAIVE(r, ci, cj, 1);
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
void DOSUB_PAIR1_SINKS(struct runner *r, struct cell *ci, struct cell *cj,
                       int timer) {
#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother? */
  const int should_do_ci = ci->sinks.count != 0 && cell_is_active_sinks(ci, e);
  const int should_do_cj = cj->sinks.count != 0 && cell_is_active_sinks(cj, e);

  if (!should_do_ci && !should_do_cj) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid_and_swap_cells(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_sinks_task(ci, cj) &&
      cell_can_recurse_in_pair_sinks_task(cj, ci)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_SINKS(r, ci->progeny[pid], cj->progeny[pjd], 0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
    const int do_ci_sink = ci->nodeID == e->nodeID;
    const int do_cj_sink = cj->nodeID == e->nodeID;
#else
    /* Here we perform the task on both sides */
    const int do_ci_sink = 1;
    const int do_cj_sink = 1;
#endif

    const int do_ci =
        ci->sinks.count != 0 && cell_is_active_sinks(ci, e) && do_ci_sink;
    const int do_cj =
        cj->sinks.count != 0 && cell_is_active_sinks(cj, e) && do_cj_sink;

    if (do_ci) {

      /* Make sure both cells are drifted to the current timestep. */
      if (!cell_are_sink_drifted(ci, e))
        error("Interacting undrifted cells (sinks).");

      if (cj->hydro.count != 0 && !cell_are_part_drifted(cj, e))
        error("Interacting undrifted cells (parts).");
    }

    if (do_cj) {

      /* Make sure both cells are drifted to the current timestep. */
      if (ci->hydro.count != 0 && !cell_are_part_drifted(ci, e))
        error("Interacting undrifted cells (parts).");

      if (!cell_are_sink_drifted(cj, e))
        error("Interacting undrifted cells (sinks).");
    }

    if (do_ci || do_cj) DOPAIR1_BRANCH_SINKS(r, ci, cj);
  }

  if (timer) TIMER_TOC(TIMER_DOSUB_PAIR_SINKS);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_SINKS(struct runner *r, struct cell *ci, int timer) {
#ifdef SWIFT_DEBUG_CHECKS_MPI_DOMAIN_DECOMPOSITION
  return;
#endif

  TIMER_TIC;

  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  const int should_do_ci = ci->sinks.count != 0 && cell_is_active_sinks(ci, e);

  if (!should_do_ci) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_sinks_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_SINKS(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1_SINKS(r, ci->progeny[k], ci->progeny[j], 0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Check we did drift to the current time */
    if (!cell_are_sink_drifted(ci, e)) error("Interacting undrifted cell.");

    if (ci->hydro.count != 0 && !cell_are_part_drifted(ci, e))
      error("Interacting undrifted cells (parts).");

    DOSELF1_BRANCH_SINKS(r, ci);
  }

  if (timer) TIMER_TOC(TIMER_DOSUB_SELF_SINKS);
}
