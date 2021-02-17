/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
   name of the interaction function.
   This creates the required interaction functions. */

#include "active.h"
#include "runner_doiact_sinks_merger.h"

/**
 * @brief Merge the sink particles that are too close to each other.
 *
 * @param r runner task
 * @param c cell
 */
void DOSELF1_SINKS_MERGER(struct runner *r, struct cell *c) {

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != engine_rank) error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (c->sinks.count == 0) return;
  if (!cell_is_active_sinks(c, e)) return;

  /* Did we mess up the recursion? */
  if (c->sinks.r_cut_max_old > c->dmin)
    error("Cell smaller than the cut off radius");

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount = c->sinks.count;
  struct sink *restrict sinks = c->sinks.parts;

  /* Loop over the sinks in ci. */
  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith sink in ci. */
    struct sink *restrict si = &sinks[sid];

    /* Skip inhibited particles */
    if (sink_is_inhibited(si, e)) continue;

    const float ri = si->r_cut;
    const float ri2 = ri * ri;
    const float six[3] = {(float)(si->x[0] - c->loc[0]),
                          (float)(si->x[1] - c->loc[1]),
                          (float)(si->x[2] - c->loc[2])};

    /* Loop over the sinks in cj. */
    for (int sjd = sid + 1; sjd < scount; sjd++) {

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks[sjd];

      /* Early abort? */
      if (!sink_is_active(si, e) && !sink_is_active(sj, e)) continue;
      if (sink_is_inhibited(sj, e)) continue;

      /* Get the cutoff radius */
      const float rj = sj->r_cut;
      const float rj2 = rj * rj;

      /* Compute the pairwise distance. */
      const float sjx[3] = {(float)(sj->x[0] - c->loc[0]),
                            (float)(sj->x[1] - c->loc[1]),
                            (float)(sj->x[2] - c->loc[2])};
      float dx[3] = {six[0] - sjx[0], six[1] - sjx[1], six[2] - sjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
      if (si->ti_drift != e->ti_current)
        error("Particle si not drifted to current time");
#endif
      if (r2 < ri2 || r2 < rj2) {
        enum sink_merger_remove remove =
            IACT_SINK_MERGER(r2, dx, ri, rj, si, sj, a, H);

        /* Remove the particle. */
        switch (remove) {
          case sink_merger_remove_none:
            break;
          case sink_merger_remove_first:
            cell_remove_sink(e, c, si);
            break;
          case sink_merger_remove_second:
            cell_remove_sink(e, c, sj);
            break;
          default:
            error("Unknown value, please check your iact function.");
        }

        /* Can we continue the loop on the current particle? */
        if (remove == sink_merger_remove_first) {
          break;
        }
      }
    } /* loop over the sinks in ci. */
  }   /* loop over the sinks in ci. */
}

/**
 * @brief Merge the sink particles that are too close to each other.
 *
 * @param r runner task
 * @param ci The first #cell
 * @param cj The second #cell
 */
void DO_SYM_PAIR1_SINKS_MERGER(struct runner *r, struct cell *restrict ci,
                               struct cell *restrict cj) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank && cj->nodeID != engine_rank)
    error("Should be run on a different node");
#endif

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  /* Anything to do here? */
  if (cj->sinks.count == 0 || ci->sinks.count == 0) return;
  if (!cell_is_active_sinks(ci, e) && !cell_is_active_sinks(cj, e)) return;

  const int ci_local = ci->nodeID == e->nodeID;
  const int cj_local = cj->nodeID == e->nodeID;

  if (!ci_local || !cj_local) error("TODO");

  /* Cosmological terms */
  const float a = cosmo->a;
  const float H = cosmo->H;

  const int scount_i = ci->sinks.count;
  const int scount_j = cj->sinks.count;
  struct sink *restrict sinks_i = ci->sinks.parts;
  struct sink *restrict sinks_j = cj->sinks.parts;

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  /* Loop over the sinks in ci. */
  for (int sid = 0; sid < scount_i; sid++) {

    /* Get a hold of the ith sink in ci. */
    struct sink *restrict si = &sinks_i[sid];

    /* Skip inhibited particles */
    if (sink_is_inhibited(si, e)) continue;

    const float ri = si->r_cut;
    const float ri2 = ri * ri;
    const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                          (float)(si->x[1] - (cj->loc[1] + shift[1])),
                          (float)(si->x[2] - (cj->loc[2] + shift[2]))};

    /* Loop over the sinks in cj. */
    for (int sjd = 0; sjd < scount_j; sjd++) {

      /* Get a pointer to the jth particle. */
      struct sink *restrict sj = &sinks_j[sjd];

      /* Early abort? */
      if (!sink_is_active(si, e) && !sink_is_active(sj, e)) continue;
      if (sink_is_inhibited(sj, e)) continue;

      /* Get the cutoff radius */
      const float rj = sj->r_cut;
      const float rj2 = rj * rj;

      /* Compute the pairwise distance. */
      const float pjx[3] = {(float)(sj->x[0] - cj->loc[0]),
                            (float)(sj->x[1] - cj->loc[1]),
                            (float)(sj->x[2] - cj->loc[2])};
      float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
      const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
      /* Check that particles have been drifted to the current time */
      if (sj->ti_drift != e->ti_current)
        error("Particle sj not drifted to current time");
      if (si->ti_drift != e->ti_current)
        error("Particle si not drifted to current time");
#endif

      if (r2 < ri2 || r2 < rj2) {
        enum sink_merger_remove remove =
            IACT_SINK_MERGER(r2, dx, ri, rj, si, sj, a, H);

        /* Remove the particle. */
        switch (remove) {
          case sink_merger_remove_none:
            break;
          case sink_merger_remove_first:
            cell_remove_sink(e, ci, si);
            break;
          case sink_merger_remove_second:
            cell_remove_sink(e, cj, sj);
            break;
          default:
            error("Unknown value, please check your iact function.");
        }

        /* Can we continue the loop on the current particle? */
        if (remove == sink_merger_remove_first) {
          break;
        }
      }
    } /* loop over the sinks in cj. */
  }   /* loop over the sinks in ci. */
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1_SINKS_MERGER(struct runner *r, struct cell *ci,
                              struct cell *cj) {

  struct space *s = r->e->s;
  const struct engine *e = r->e;

#ifdef WITH_MPI
  if (ci->nodeID != e->nodeID && cj->nodeID != e->nodeID)
    error("Should not be running on this node");
#endif

  /* Should we even bother? */
  if (!cell_is_active_sinks(ci, e) && !cell_is_active_sinks(cj, e)) return;
  if (ci->sinks.count == 0 || cj->sinks.count == 0) return;

  /* Get the type of pair and flip ci/cj if needed. */
  double shift[3];
  const int sid = space_getsid(s, &ci, &cj, shift);

  /* Recurse? */
  if (cell_can_recurse_in_pair_sinks_task(ci, cj) &&
      cell_can_recurse_in_pair_sinks_task(cj, ci)) {
    struct cell_split_pair *csp = &cell_split_pairs[sid];
    for (int k = 0; k < csp->count; k++) {
      const int pid = csp->pairs[k].pid;
      const int pjd = csp->pairs[k].pjd;
      if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
        DOSUB_PAIR1_SINKS_MERGER(r, ci->progeny[pid], cj->progeny[pjd]);
    }
  }

  /* Otherwise, compute the pair directly. */
  else {

    /* Make sure both cells are drifted to the current timestep. */
    if (!cell_are_sink_drifted(ci, e))
      error("Interacting undrifted cells (ci).");

    if (!cell_are_sink_drifted(cj, e))
      error("Interacting undrifted cells (cj).");

    DO_SYM_PAIR1_SINKS_MERGER(r, ci, cj);
  }
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 */
void DOSUB_SELF1_SINKS_MERGER(struct runner *r, struct cell *ci) {

#ifdef SWIFT_DEBUG_CHECKS
  if (ci->nodeID != engine_rank)
    error("This function should not be called on foreign cells");
#endif

  /* Should we even bother? */
  if (ci->sinks.count == 0 || !cell_is_active_sinks(ci, r->e)) return;

  /* Recurse? */
  if (cell_can_recurse_in_self_sinks_task(ci)) {

    /* Loop over all progeny. */
    for (int k = 0; k < 8; k++)
      if (ci->progeny[k] != NULL) {
        DOSUB_SELF1_SINKS_MERGER(r, ci->progeny[k]);
        for (int j = k + 1; j < 8; j++)
          if (ci->progeny[j] != NULL)
            DOSUB_PAIR1_SINKS_MERGER(r, ci->progeny[k], ci->progeny[j]);
      }
  }

  /* Otherwise, compute self-interaction. */
  else {

    /* Drift the cell to the current timestep if needed. */
    if (!cell_are_sink_drifted(ci, r->e)) {
      error("Interacting undrifted cell.");
    }

    DOSELF1_SINKS_MERGER(r, ci);
  }
}
