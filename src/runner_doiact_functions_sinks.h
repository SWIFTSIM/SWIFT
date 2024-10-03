/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
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
 * @brief Calculate gas and sink interaction around #sinks
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void DOSELF1_SINK(struct runner *r, struct cell *c, int timer) {

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

      const float ri = si->r_cut;
      const float ri2 = ri * ri;
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

        if (r2 < ri2) {
          IACT_SINK_GAS(
              r2, dx, ri, hj, si, pj, with_cosmology, cosmo,
              e->gravity_properties, e->sink_properties);
        }
      } /* loop over the parts in ci. */
    } /* loop over the bparts in ci. */
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

    const float ri = si->r_cut;
    const float ri2 = ri * ri;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the sinks in cj. */
    for (int sjd = 0; sjd < scount; sjd++) {

      /* Skip self interaction */
      if (sid == sjd) continue;

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks[sjd];
      const float rj = sj->r_cut;
      const float rj2 = rj * rj;

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
        error("Particle bi not drifted to current time");
      if (sj->ti_drift != e->ti_current)
        error("Particle bj not drifted to current time");
#endif

      if (r2 < ri2 || r2 < rj2) {
        IACT_SINK_SINK(r2, dx, ri, rj, si, sj,
                                              with_cosmology, cosmo,
                                              e->gravity_properties);
      }
    } /* loop over the sinks in ci. */
  } /* loop over the sinks in ci. */

#endif /* (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW) */

  if (timer) TIMER_TOC(TIMER_DOSELF_SINK);
}

/**
 * @brief Calculate gas and sink interaction around #sinks
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_NONSYM_PAIR1_SINK_NAIVE(struct runner *r,
                                struct cell *restrict ci,
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

      /* Get a hold of the ith bpart in ci. */
      struct sink *restrict si = &sinks_i[sid];

      /* Skip inactive particles */
      if (!sink_is_active(si, e)) continue;

      const float ri = si->r_cut;
      const float ri2 = ri * ri;
      const float six[3] = {(float)(si->x[0] - cj->loc[0]),
                            (float)(si->x[1] - cj->loc[1]),
                            (float)(si->x[2] - cj->loc[2])};

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

        if (r2 < ri2) {
          IACT_SINK_GAS(
              r2, dx, ri, hj, si, pj, with_cosmology, cosmo,
              e->gravity_properties, e->sink_properties);
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

    /* Get a hold of the ith bpart in ci. */
    struct sink *restrict si = &sinks_i[sid];

    /* Skip inactive particles */
    if (!sink_is_active(si, e)) continue;

    const float ri = si->r_cut;
    const float ri2 = ri * ri;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the sinks in cj. */
    for (int sjd = 0; sjd < scount_j; sjd++) {

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks_j[sjd];
      const float rj = sj->r_cut;
      const float rj2 = rj * rj;

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

      if (r2 < ri2 || r2 < rj2) {
        IACT_SINK_SINK(r2, dx, ri, rj, si, sj,
                                              with_cosmology, cosmo,
                                              e->gravity_properties);
      }
    } /* loop over the sinks in cj. */
  } /* loop over the sinks in ci. */

#endif /* (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW) */
}

/**
 * @brief Calculate swallow for ci #sinks part around the cj #gas and sinks and
 *                              cj #sinks part around the ci #gas and sinks
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DOPAIR1_SINK_NAIVE(struct runner *r,
                        struct cell *restrict ci,
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

  if (do_ci_sink) DO_NONSYM_PAIR1_SINK_NAIVE(r, ci, cj);
  if (do_cj_sink) DO_NONSYM_PAIR1_SINK_NAIVE(r, cj, ci);

  if (timer) TIMER_TOC(TIMER_DOPAIR_SINK);
}

/**
 * @brief Wrapper to runner_doself_sinks_swallow
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH_SINK(struct runner *r, struct cell *c) {

  const struct engine *restrict e = r->e;

  /* Anything to do here? */
  if (c->sinks.count == 0) return;

  /* Anything to do here? */
  if (!cell_is_active_sinks(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->sinks.r_cut_max_old > c->dmin)
    error("Cell smaller than the cut off radius");

  DOSELF1_SINK(r, c, 1);
}

/**
 * @brief Wrapper for runner_dopair_sinks_naive_swallow.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH_SINK(struct runner *r, struct cell *ci,
                          struct cell *cj) {

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

  const int do_ci = (ci->sinks.count != 0 && cj->hydro.count != 0 &&
                     ci_active && do_ci_sink);
  const int do_cj = (cj->sinks.count != 0 && ci->hydro.count != 0 &&
                     cj_active && do_cj_sink);

  /* Anything to do here? */
  if (!do_ci && !do_cj) return;

  /* Check that cells are drifted. */
  if (do_ci && (!cell_are_sink_drifted(ci, e) || !cell_are_part_drifted(cj, e)))
    error("Interacting undrifted cells.");

  if (do_cj && (!cell_are_part_drifted(ci, e) || !cell_are_sink_drifted(cj, e)))
    error("Interacting undrifted cells.");

  /* No sorted interactions here -> use the naive ones */
  DOPAIR1_SINK_NAIVE(r, ci, cj, 1);
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
void DOSUB_PAIR1_SINK(struct runner *r, struct cell *ci,
                                     struct cell *cj, int timer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const struct engine *e = r->e;

  /* Should we even bother?
   * In the swallow case we care about sink-sink and sink-gas
   * interactions.
   * In all other cases only sink-gas so we can abort if there is
   * is no gas in the cell
   * ^^ this is from BHs - do we care about sink-sink density interactions? */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
  const int should_do_ci =
      ci->sinks.count != 0 && cell_is_active_sinks(ci, e);
  const int should_do_cj =
      cj->sinks.count != 0 && cell_is_active_sinks(cj, e);
#else
  const int should_do_ci = ci->sinks.count != 0 && cj->hydro.count != 0 &&
                           cell_is_active_sinks(ci, e);
  const int should_do_cj = cj->sinks.count != 0 && ci->hydro.count != 0 &&
                           cell_is_active_sinks(cj, e);

#endif

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
        DOSUB_PAIR1_SINK(r, ci->progeny[pid], cj->progeny[pjd],
                                        0);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

#if (FUNCTION_TASK_LOOP == TASK_LOOP_DENSITY)
    const int do_ci_sink = ci->nodeID == e->nodeID;
    const int do_cj_sink = cj->nodeID == e->nodeID;
#elif (FUNCTION_TASK_LOOP == TASK_LOOP_FEEDBACK)
    /* Here we are updating the hydro -> switch ci, cj */
    const int do_ci_sink = cj->nodeID == e->nodeID;
    const int do_cj_sink = ci->nodeID == e->nodeID;
#else
    /* Here we perform the task on both sides */
    const int do_ci_sink = 1;
    const int do_cj_sink = 1;
#endif

    const int do_ci = ci->sinks.count != 0 &&
                      cell_is_active_sinks(ci, e) && do_ci_sink;
    const int do_cj = cj->sinks.count != 0 &&
                      cell_is_active_sinks(cj, e) && do_cj_sink;

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

    if (do_ci || do_cj) DOPAIR1_BRANCH_SINK(r, ci, cj);
  }

  if (timer) TIMER_TOC(TIMER_DOSUB_PAIR_SINK);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1_SINK(struct runner *r, struct cell *ci,
                                     int timer) {

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
     * is no gas in the cell 
     * NB - same concern as above */
#if (FUNCTION_TASK_LOOP == TASK_LOOP_SWALLOW)
  const int should_do_ci =
      ci->sinks.count != 0 && cell_is_active_sinks(ci, e);
#else
  const int should_do_ci = ci->sinks.count != 0 && ci->hydro.count != 0 &&
                           cell_is_active_sinks(ci, e);
#endif

  /* Should we even bother?
   * In the swallow case we care about sink-sink and sink-gas
   * interactions. */

  const int should_do_ci = ci->sinks.count != 0 && cell_is_active_sinks(ci, e);

  if (!should_do_ci) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_sinks_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_SINK(r, ci->progeny[k], 0);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_SELF1_SINK(r, ci->progeny[k], ci->progeny[j],
                                            0);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Check we did drift to the current time */
    if (!cell_are_sink_drifted(ci, e)) error("Interacting undrifted cell.");

    if (ci->hydro.count != 0 && !cell_are_part_drifted(ci, e))
      error("Interacting undrifted cells (parts).");

    DOSELF1_BRANCH_SINK(r, ci);
  }

  if (timer) TIMER_TOC(TIMER_DOSUB_SELF_SINK);
}
